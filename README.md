# KIRVCF

A collection of scripts to parse VCF files and run analysis. KIRVCFparser compiles the genotype frequencies for each position from selected VCF files. VCFparser is a more general parser for unpacking .vcf files. The remaining scripts run specialized analysis. 


## Future Work

In the future this collection of scripts will be more succinctly wrapped into a Python package, but for now there are command line options for many of the scripts. 

## Tech Specs

joblib is the only external python package that these scripts depend on. Otherwise all scripts should run in Python 3. They have not been tested yet in Python 2.
