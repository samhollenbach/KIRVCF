import csv
import vcf
import os
from joblib import Parallel, delayed
import math
import itertools


class VCFparser:

    def __init__(self):
        self.n_jobs = 4
        self.output_dir = "../results/"

    def set_njobs(self, n):
        self.n_jobs = n

    def set_output_dir(self, directory):
        self.output_dir = directory


    # Given a list of VCF records from a single file, returns all records with DP above *dp_lim*
    # If less than *record_lim* records fit the specifications, returns None
    def filter_records(self, record_list, dp_lim=21, record_lim=6):
        approved_recs = [r for r in record_list if r.INFO['DP'] >= dp_lim]
        return approved_recs if len(approved_recs) >= record_lim else None


    # Returns list of all filepaths ending with *file_ending* in *directory*
    def get_files_in_dir(self, directory, file_ending = ""):
        return [os.path.join(directory, file) for file in os.listdir(directory) if file.endswith(file_ending)]

    def get_dirs_in_dir(self, directory, exclude_dirs=None, make_dict=False, sub_folder=""):
        if exclude_dirs is None:
            exclude_dirs = ["Vcf"]
        if make_dict:
            return { dir : os.path.join(directory, dir, sub_folder) for dir in os.listdir(directory) if os.path.isdir(os.path.join(directory, dir)) and dir not in exclude_dirs and not dir.startswith(".")}
        else:
            return [dir for dir in os.listdir(directory) if os.path.isdir(os.path.join(directory, dir)) and dir not in exclude_dirs and not dir.startswith(".")]

    # Returns list of all records in specified vcf file
    # TODO: ADD PARAMS FOR FILTER IF NECESSARY
    def read_vcf(self, file, filter_recs=True):
        recs = [line for line in vcf.Reader(open(file), 'r')]
        if filter_recs:
            recs = self.filter_records(recs)
        return recs


    # Returns a dictionary of position:genotype for all records in *record_list*
    def process_genotypes(self, record_list, only_variants=False):
        if record_list is None:
            return {}
        genos = {}
        for record in record_list:
            ref = record.REF
            alts = record.ALT
            if None not in alts:
                alt = str(alts[0])
                unformatted_geno = [sample['GT'] for sample in record.samples][0]

                # Replace 0 with REF, 1 with ALT
                genotype = unformatted_geno.replace('0', ref).replace('1', alt)
            else:
                if only_variants:
                    continue
                genotype = '{}/{}'.format(ref, ref)
            genos[record.POS] = genotype
        return genos


    def apply_function_in_parallel(self, list, func, *args, n_jobs=None):
        if n_jobs is None:
            n_jobs = self.n_jobs
        with Parallel(n_jobs=n_jobs) as parallelizer:
            group_size = math.ceil(len(list) / n_jobs)

            # Splits the filenames list into groups of size group_size
            buffer_names = [iter(list)] * group_size
            split_names = itertools.zip_longest(fillvalue=None, *buffer_names)

            # Assings split_names groups to different
            task_iterator = (delayed(func)(file,*args) for file in split_names)
            results = parallelizer(task_iterator)
            return results

    def write_output_file(self, filename, col_names, data, delim = "\t", output_dir=None,):
        if output_dir is None:
            output_dir = self.output_dir
        print("\nWriting output to {}\n".format(output_dir+filename))
        with open(output_dir+filename, 'w') as w:
            writer = csv.writer(w, delimiter=delim)
            writer.writerow(col_names)
            for row in data:
                writer.writerow(row)


    def gather_private_alleles(self, sample_paths_dict):
        #### { Sample : { "pos:var" : count, "pos2:var2": count2 }, Sample2 : { "pos3:var3" : count3 } }
        private_alleles = { sn : {} for sn in sample_paths_dict.keys() }
        public_list = []
        for sample_name, sample_path in sample_paths_dict.items():

            #dir = os.path.join(path_to_samples, sample, sub_folder)
            files = self.get_files_in_dir(sample_path, file_ending=".vcf")
            for file in files:
                # Will return empty dict if recs = None
                genos = self.process_genotypes(self.read_vcf(file), only_variants=True)
                print(sample_name,genos)
                for pos, var in genos.items():
                    g = "{}:{}".format(pos, var)
                    if g in public_list:
                        continue
                    private = True

                    for s in sample_paths_dict.keys():
                        if s is sample_name:
                            continue
                        if g in private_alleles[s]:
                            private_alleles[s].pop(g, None)
                            public_list.append(g)
                            private = False

                    if private:
                        if g not in private_alleles[sample_name]:
                            private_alleles[sample_name][g] = 0
                        private_alleles[sample_name][g] += 1
        print("Private Alleles: ", private_alleles)
        return private_alleles



