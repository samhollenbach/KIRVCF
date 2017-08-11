import csv
import vcf
import os
from joblib import Parallel, delayed
import math
import itertools


class VCFparser:

    def __init__(self):
        self.n_jobs = 4
        self.output_dir = "results/"

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

    def get_dirs_in_dir(self, directory, exclude_dirs=None, full_path=False, sub_folder=""):
        if exclude_dirs is None:
            exclude_dirs = ["Vcf"]
        if full_path:
            return [os.path.join(directory, dir, sub_folder) for dir in os.listdir(directory) if dir not in exclude_dirs and not dir.startswith(".")]
        else:
            return [dir for dir in os.listdir(directory) if dir not in exclude_dirs and not dir.startswith(".")]

    # Returns list of all records in specified vcf file
    # TODO: ADD PARAMS FOR FILTER IF NECESSARY
    def read_vcf(self, file, filter_recs=True):
        recs = [line for line in vcf.Reader(open(file), 'r')]
        if filter_recs:
            recs = self.filter_records(recs)
        return recs


    # Returns a dictionary of position:genotype for all records in *record_list*
    def process_genotypes(self, record_list, only_variants=False):
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
        with open(output_dir+filename, 'w') as w:
            writer = csv.writer(w, delimiter=delim)
            writer.writerow(col_names)
            for row in data:
                writer.writerow(row)


    def gather_private_alleles(self, sample_names, path_to_samples = "", sub_folder=""):

        private_alleles = {}
        for sample in sample_names:
            dir = os.path.join(path_to_samples, sample, sub_folder)
            files = self.get_files_in_dir(dir, file_ending=".vcf")
            for file in files:
                recs = self.read_vcf(file)
                if recs is not None:
                    genos = self.process_genotypes(recs, only_variants=True)
                    print(sample,genos)




vcfp = VCFparser()
inputfolder = "../input_test"
sample_names = vcfp.get_dirs_in_dir(inputfolder, full_path=False)
vcfp.gather_private_alleles(sample_names, path_to_samples=inputfolder)