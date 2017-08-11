import os
import csv
from VCFparser import VCFparser



vcfp = VCFparser()
inputfolder = "../input_test"
sample_names = vcfp.get_dirs_in_dir(directory=inputfolder, full_path=False)
vcfp.gather_private_alleles(sample_names, path_to_samples=inputfolder)

header = ["Country", "Position", "Variant", "Count"]
