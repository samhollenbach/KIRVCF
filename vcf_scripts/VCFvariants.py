import vcf
import os
from joblib import Parallel, delayed
import math
import itertools
import multiprocessing
import csv

vcf_file_ending = ".vcf"


variant_cap = 5
use_variant_cap = False

include_filenames = False






def open_files(country, main_dir, res_dir):
    input_dir = "{}/{}/Vcf".format(main_dir, country)
    #input_dir = "{}/{}".format(main_dir, country)
    vars_dir = "{}/vars".format(res_dir)
    vars_outfile = "{}/{}_variants.txt".format(vars_dir, country)

    all_vars = {}
    if not os.path.exists(vars_dir):
        os.makedirs(vars_dir)

    print('Reading VCF files...\n')
    # Stores all the files with specified ending in the input directory

    vcf_filenames = [os.path.join(input_dir, file) for file in os.listdir(input_dir) if file.endswith(vcf_file_ending)]

    print(
        "Found {} Vcf files in {}, checking validity and processing records...\n".format(len(vcf_filenames), input_dir))
    total_files_num = len(vcf_filenames)
    n_jobs = 40

    with Parallel(n_jobs=n_jobs) as parallelizer:
        group_size = math.ceil(len(vcf_filenames) / n_jobs)

        # Splits the filenames list into groups of size group_size
        buffer_names = [iter(vcf_filenames)] * group_size
        split_names = itertools.zip_longest(fillvalue=None, *buffer_names)

        # Assings split_names groups to different
        task_iterator = (delayed(read_files)(files) for files in split_names)
        results = parallelizer(task_iterator)
        for r in results:
            for p, d in r.items():
                if p in all_vars:
                    all_vars[p].update(d)
                else:
                    all_vars[p] = d

        return process_variants(all_vars, total_files_num, vars_outfile)


def read_files(filenames):
    files_remaining = len(filenames)
    temp_vars = {}
    for filename in filenames:
        files_remaining -= 1

        print("{}: {} files remaining".format(multiprocessing.current_process().name, files_remaining))
        if filename is None:
            continue
        vcf_values = [line for line in vcf.Reader(open(filename), 'r')]
        new_vars = check_variants(filename.split("/")[-1].split(".vcf")[0], vcf_values)
        for p, d in new_vars.items():
            if p in temp_vars:
                temp_vars[p].update(d)
            else:
                temp_vars[p] = d
    return temp_vars


def check_variants(filename, vcf_vals):
    variants = {}

    for record in vcf_vals:
        ref = record.REF
        alts = record.ALT
        if None in alts:
            continue
        alt = str(alts[0])
        unformatted_geno = [sample['GT'] for sample in record.samples][0]
        genotype = unformatted_geno.replace('0', ref).replace('1', alt)
        if record.POS not in variants.keys():
            variants[record.POS] = {filename: genotype}
        else:
            variants[record.POS].update({filename: genotype})

    return variants


def process_variants(avar, total_files_num, outfile):
    count_ciel = variant_cap
    single_variant_files = {}
    multi_variant_files = {}

    for pos, variants in avar.items():
        num_variants = len(variants.values())
        if num_variants > count_ciel:
            continue
        elif num_variants == 1:
            single_variant_files[pos] = variants
        else:
            multi_variant_files[pos] = variants

    return print_variants(single_variant_files, multi_variant_files, total_files_num, outfile)


def print_variants(single, multiple, total_num, outfile):
    single_num = len(single.keys())
    multi_num = len(multiple.keys())
    print("Found {} single variant files".format(single_num))
    print("Found {} multi variant files".format(multi_num))
    print("Writing results to {}".format(outfile))

    with open(outfile, 'w') as w:
        w.write("Single Variant Files ({}/{})\n".format(single_num, total_num))

        for pos, vars in single.items():
            w.write("Position: {}".format(pos))
            for f, g in vars.items():
                w.write("\tFile: {}\tVariant: {}".format(f, g))
            w.write("\n\n")

        w.write("\n\nMulti Variant Files ({}/{})\n".format(multi_num, total_num))

        for pos, vars in multiple.items():
            w.write("Position: {}".format(pos))
            for f, g in vars.items():
                w.write("\tFile: {}\tVariant: {}".format(f, g))
            w.write("\n\n")

    return {"Single": single_num, "Multiple": multi_num, "Total": total_num}


def run_all_countries():
    exclude = ["Vcf"]
    main_dir = "/home/jurgen/Results/"

    #main_dir = "input"
    countries = [file for file in os.listdir(main_dir) if os.path.isdir(os.path.join(main_dir, file))]
    results = "results"
    if not os.path.exists(results):
        os.makedirs(results)
    with open("{}/variant_counts.txt".format(results), 'w') as w:
        writer = csv.writer(w, delimiter="\t")
        writer.writerow(["Country", "Single", "Multiple", "Total"])
        for c in countries:
            if c in exclude:
                continue
            print("\nPreparing variant counter for {}\n\n".format(c))
            counts_dict = open_files(c, main_dir, results)
            single_percent = counts_dict["Single"] / counts_dict["Total"]
            multi_percent = counts_dict["Multiple"] / counts_dict["Total"]
            line = [c, "{} ({})".format(counts_dict["Single"], single_percent),
                    "{} ({})".format(counts_dict["Multiple"], multi_percent), counts_dict["Total"]]
            writer.writerow(line)


run_all_countries()
