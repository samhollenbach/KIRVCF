import os
import csv
from VCFparser import VCFparser



vcfp = VCFparser()
inputfolder = "/home/jurgen/Results"
sample_paths = vcfp.get_dirs_in_dir(directory=inputfolder, make_dict=True, sub_folder="Vcf")
print(sample_paths)
pa = vcfp.gather_private_alleles(sample_paths)

header = ["Country", "Position", "Variant", "Count"]
rows = []
for sample, private_alleles in pa.items():
    for priv, count in private_alleles.items():
        p = priv.split(":")
        rows.append([sample, p[0], p[1], count])

vcfp.write_output_file("private_alleles.txt", header, rows)
