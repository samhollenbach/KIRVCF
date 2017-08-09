import csv
import os
import argparse

output_files_dir = "."
output_files_format = "_out.csv"
merged_output_name = "merged_output.csv"
file_counter_file = "file_counter.csv"

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", metavar="OUTPUT_FILES_DIR", help="Directory containing _output.csv files to merge")
parser.add_argument("-o", "--output", metavar="MERGED_OUTPUT_NAME", help="Name of the merged output file")
parser.add_argument("-f", "--format", metavar="OUTPUT_FILES_FORMAT", help="The ending of the desired output files")

args = parser.parse_args()
if args.input:
    output_files_dir = args.input
if args.format:
    output_files_format = args.format
if args.output:
    merged_output_name = args.output

output_file_names = [os.path.join(output_files_dir, file) for file in os.listdir(output_files_dir) if file.endswith(output_files_format)]

sample_list = ['pos']
genotype_list = ['A/A', 'A/C', 'A/G', 'A/T', 'C/A', 'C/C', 'C/G', 'C/T', 'G/A', 'G/C', 'G/G', 'G/T', 'T/A', 'T/C', 'T/G', 'T/T']
master_merged = {}
files_dict = {}

for file in output_file_names:
    for row in csv.reader(open(file),delimiter=","):
        pos = row[0]
        if pos == 'pos':
            sample = row[1]
            sample_list.append(sample)
            files_dict[sample] = num_files
            continue
        if pos == 'files':
            num_files = row[1]
            continue
        if pos not in master_merged:
            master_merged[pos] = {sample : row[1:]}
        else:
            master_merged[pos][sample] = row[1:]

max_col_width = max([len(l) for l_samples in master_merged.values() for l in l_samples.values()])

with open(file_counter_file,"w") as ff:
    ww = csv.writer(ff,delimiter=",")
    ww.writerow(["Sample", "Number of Files"])
    for sample, file_num in files_dict.items():
        ww.writerow([sample, file_num])


with open(merged_output_name,"w") as f:
    w = csv.writer(f, delimiter=",")
    header1 = ['Genotypes']
    header2 = [sample_list.pop(0)]
    sample_list.sort()
    num_samples = len(sample_list)

    # Write first row of all genotypes
    for g in genotype_list:
        header1 += ([g] + ['']*(num_samples-1))
    w.writerow(header1)
    # Write second row with all sampels for each genotype
    header2 += sample_list*len(genotype_list)
    w.writerow(header2)

    for position in sorted(master_merged.keys()):
        pos_dict = master_merged[position]
        print(pos_dict)
        row_no_pos = ['']*(len(header1)-1)
        for sample_at_pos, g_list_at_pos in pos_dict.items():
            sample_index = sample_list.index(sample_at_pos)
            print(sample_at_pos, sample_index)
            for geno_at_pos_combo in g_list_at_pos:
                geno_at_pos = geno_at_pos_combo.split(':')
                geno_index = genotype_list.index(geno_at_pos[0])
                if geno_at_pos[1] != '1.0':
                    row_no_pos[geno_index*num_samples+sample_index] = geno_at_pos[1]
        if all(x is '' for x in row_no_pos):
            continue
        row = [position] + row_no_pos
        w.writerow(row)

