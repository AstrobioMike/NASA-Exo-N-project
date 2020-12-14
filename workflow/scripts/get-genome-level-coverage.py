#!/usr/bin/env python

"""
This is an ad hoc script for the corresponding workflow. It takes output from pileup.sh and returns genome-level coverage and detection info tables. 
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
from statistics import mean
import numpy

parser = argparse.ArgumentParser(description="This is an ad hoc script for the corresponding workflow. It takes output from pileup.sh and returns genome-level coverage and detection info tables.")

required = parser.add_argument_group('required arguments')

required.add_argument("-i", "--input-table", help="Output table from pileup.sh", action="store", dest="input_table")
required.add_argument("-g", "--unique-genome-IDs-file", help="Single-column file with unique genome identifiers", action="store", dest="unique_genomes")

parser.add_argument("-d", "--detection-threshold", help="Detection threshold. Lower than this, and the genome coverage will be set to 0 in the 'filtered-coverage' output table. Should be within 0 and 1 (default: 0.5)", type=float, action="store", dest="detection_threshold", default=0.5)
parser.add_argument("-o", "--output-prefix", help='Output prefix (default: "Sample")', action="store", default="Sample", dest="output_prefix")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

################################################################################


def main():

    # reading target genomes into list
    target_genomes = [line.strip() for line in open(args.unique_genomes)]

    # initializing coverage and detection dataframes
    cov_df = pd.DataFrame({"Genome_ID" : target_genomes})
    det_df = pd.DataFrame({"Genome_ID" : target_genomes})
    filtered_cov_df = pd.DataFrame({"Genome_ID" : target_genomes})

    
    # reading input table
    tab = pd.read_csv(args.input_table, sep="\t", usecols = [0,1,2,4,5])
    tab.rename(columns = {"#ID":"ID"}, inplace = True)

    # initializing coverage and detection lists
    cov_list = []
    det_list = []

    # looping through target genomes and summarizing their values from their contig-level info
    for genome in target_genomes:
        curr_sub = tab[tab['ID'].str.startswith(genome + "_")]
        cov_list.append(mean(curr_sub.Avg_fold))
        det_list.append(sum(curr_sub.Covered_bases) / sum(curr_sub.Length))

    # adding summary values to coverage and detection dataframes
    cov_df['coverage'] = cov_list
    det_df['detection'] = det_list

    # creating detection filtered table
    # getting which from the coverage table have less than the specified
    boolean_vec_of_those_with_too_low_detection = list(map(lambda x: x < args.detection_threshold, det_list))
    indices_to_change = [i for i, x in enumerate(boolean_vec_of_those_with_too_low_detection) if x]

    # making array so easier to manipulate values based on indices
    cov_array = numpy.array(cov_list)
    cov_array[indices_to_change] = 0

    filtered_cov_df['coverage'] = cov_array

    # writing out tables
    cov_df.to_csv(args.output_prefix + "-genome-level-coverages.tsv", index=False, sep="\t", na_rep = "NA")
    det_df.to_csv(args.output_prefix + "-genome-level-detections.tsv", index=False, sep="\t", na_rep = "NA")
    filtered_cov_df.to_csv(args.output_prefix + "-genome-level-detection-filtered-coverages.tsv", index=False, sep="\t", na_rep = "NA")


if __name__ == "__main__":
    main()
