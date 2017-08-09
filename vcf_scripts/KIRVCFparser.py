import vcf
import os
import csv
import time
import math
from collections import defaultdict
from joblib import Parallel, delayed
import argparse
import itertools

# Directory containing the input .vcf files (including '/')
input_dir = 'input/3DL2'
# Output file
output_dir = '3DL2_output.csv'
# Valid VCF file ending
vcf_file_ending = '.vcf'
# Master dictionary containing all found positions with genotype frequencies for each parallel system
master = {}
# Number of jobs to use in multiprocessing
n_jobs = 3

parser = argparse.ArgumentParser(
    description="Compiles the genotype frequencies for each position from selected VCF files")
parser.add_argument("-i", "--input", metavar="INPUT_DIR", help="Directory containing VCF files")
parser.add_argument("-o", "--output", metavar="OUTPUT_DIR", help="Desired path of output file (with .csv ending)")
parser.add_argument("-j", "--jobs", metavar="NUM_JOBS", help="How many cores to use in multiprocessing", type=int)
parser.add_argument("-s", "--sample", metavar="SAMPLE_NAME", help="Name of the sample (Defaults to name of input dir)")
parser.add_argument("-f", "--vcf-format", metavar="FILE_ENDING", dest="format",
                    help="Only process VCF files with this ending (.vcf by default)")
args = parser.parse_args()

if args.input:
    input_dir = args.input
sample_name = input_dir.split("/", 1)[-1]
if args.sample:
    sample_name = args.sample
if args.output:
    output_dir = args.output
if args.jobs:
    n_jobs = args.jobs
if args.format:
    vcf_file_ending = args.format

input_dir = input_dir + "/"


# Creates a list of vcf readers along with the names of the files
def open_files():
    print('Reading VCF files...\n')
    # Stores all the files with specified ending in the input directory

    vcf_filenames = [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith(vcf_file_ending)]
    print(
        "Found {} Vcf files in {}, checking validity and processing records...\n".format(len(vcf_filenames), input_dir))

    # Stores number of valid files found for frequency calculation
    master['good_files'] = 0
    master['bad_files'] = 0

    print("Starting parallel processing loop with {} jobs...\n".format(n_jobs))

    with Parallel(n_jobs=n_jobs) as parallelizer:
        group_size = math.ceil(len(vcf_filenames) / n_jobs)

        # Splits the filenames list into groups of size group_size
        buffer_names = [iter(vcf_filenames)] * group_size
        split_names = itertools.zip_longest(fillvalue=None, *buffer_names)

        # Assings split_names groups to different
        task_iterator = (delayed(read_files)(files) for files in split_names)
        results = parallelizer(task_iterator)
        write_file(merge_results(results))


# This method is called seperately from each multiprocessing pool to read its set of files
def read_files(filenames):
    files_remaining = len(filenames)
    for filename in filenames:
        if filename is None:
            continue
        files_remaining -= 1
        vcf_values = [line for line in vcf.Reader(open(filename), 'r')]
        # Returns true if the record is valid
        if process_records(vcf_values):
            master['good_files'] += 1
            print("Vcf file \'{}\' finished processesing. {} valid files processed / "
                  "{} remaining for this pool".format(filename, master['good_files'], files_remaining))
        else:
            master['bad_files'] += 1
            print("Vcf file \'{}\' is invalid, skipping...".format(filename))
    return master


# Processes all files into list of records from the vcf readers
def process_records(vcf_vals):
    vcf_output = check_vcf(vcf_vals)
    if vcf_output is not None:
        for rec in vcf_output:
            process_genotype(rec)
        return True
    else:
        return False


# Checks a single VCF file for quality
def check_vcf(single_vcf_records):
    approved_records = [r for r in single_vcf_records if r.INFO['DP'] > 20]
    return approved_records if len(approved_records) > 5 else None


# Processes one line of a file and converts 0/0 format genotype to REF/REF, 0/1 to REF/ALT and 1/1 to ALT/ALT
def process_genotype(record):
    ref = record.REF
    alts = record.ALT
    if None not in alts:
        alt = str(alts[0])
        unformatted_geno = [sample['GT'] for sample in record.samples][0]

        # Replace 0 with REF, 1 with ALT
        genotype = unformatted_geno.replace('0', ref).replace('1', alt)
    else:
        genotype = '{}/{}'.format(ref, ref)
    add_to_master_table(record.POS, genotype)


# Adds one vcf files line to the master list
def add_to_master_table(pos, genotype):
    if pos not in master:
        master[pos] = defaultdict(int)
        master[pos]['appears'] = 0
    master[pos][genotype] += 1.
    master[pos]['appears'] += 1.


# Merges the parallelized results into a final merged master table to write
def merge_results(results):
    merged_master = {'good_files': 0, 'bad_files': 0}
    for res in results:
        for pos, geno_dict in res.items():
            if pos == 'good_files':
                merged_master['good_files'] += geno_dict
            elif pos == 'bad_files':
                merged_master['bad_files'] += geno_dict
            else:
                if pos not in merged_master:
                    merged_master[pos] = defaultdict(int)
                for geno, val in geno_dict.items():
                    merged_master[pos][geno] += val
    return merged_master


# Writes an output file with ordered positions showing the frequency of relevant genotypes
def write_file(master_dict):
    print("\nFinished compiling all data ({} files), writing output to {}...\n".format(master_dict['good_files'],
                                                                                       output_dir))
    with open(output_dir, "w") as f:
        w = csv.writer(f, delimiter=",")
        numfiles = master_dict.pop('good_files')
        num_bad_files = master_dict.pop('bad_files')
        print("{} out of {} files were invalid\n".format(num_bad_files, (numfiles + num_bad_files)))
        w.writerow(["pos", sample_name])
        for position in sorted(master_dict.keys()):
            genos_at_pos = master_dict[position]
            w.writerow(
                [position] + ["{}:{}".format(geno, (master_dict[position][geno] / genos_at_pos['appears'])) for geno in
                              sorted(genos_at_pos.keys()) if geno != 'appears'])


# Start of main program run
if __name__ == '__main__':
    start_time = time.time()
    open_files()
    elapsed = time.time() - start_time
    print("\nDone! ({:05.3f} seconds)".format(elapsed))
