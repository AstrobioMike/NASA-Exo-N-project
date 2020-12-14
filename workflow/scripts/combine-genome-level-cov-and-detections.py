#!/usr/bin/env python

"""
This is an ad hoc script for the corresponding workflow. It combines the individual sample genome-level coverage, detections, and filtered-coverage files into 3 individual tables.
It also produces CPM-normalized coverage tables for detection-filtered and the base coverage table.

Modified from my `bit-GL-combine-KO-and-tax-tables`: https://github.com/AstrobioMike/bioinf_tools
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
from math import isnan
from numpy import NaN

parser = argparse.ArgumentParser(description="This is an ad hoc script for the corresponding workflow. It combines the individual sample genome-level coverage, detection, and filtered-coverage files into 3 individual tables. \
                                              It also produces CPM-normalized coverage tables for detection-filtered and the base coverage table.")

required = parser.add_argument_group('required arguments')

required.add_argument("--input-cov-tables", metavar="genome-level-coverages.tsv", type=str, nargs="+", help="Input genome-level coverage tables (as written, expected to end with extension '.tsv'.")
required.add_argument("--input-det-tables", metavar="genome-level-detections.tsv", type=str, nargs="+", help="Input genome-level detection tables (as written, expected to end with extension '.tsv'.")
required.add_argument("--input-filtered-cov-tables", metavar="genome-level-detection-filtered-coverages.tsv", type=str, nargs="+", help="Input genome-level, detection-filtered coverage tables (as written, expected to end with extension '.tsv'.")

parser.add_argument("-o", "--output-prefix", help='Desired output prefix (default: "Combined")', action="store", default="Combined", dest="output_prefix")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

################################################################################


def main():

    check_all_inputs_exist(args.input_cov_tables, args.input_det_tables, args.input_filtered_cov_tables)

    input_cov_files, input_det_files, input_filt_cov_files, sample_names = setup_input_lists(args.input_cov_tables, args.input_det_tables, args.input_filtered_cov_tables)

    cov_combined_tab, detection_combined_tab, filtered_cov_combined_tab, norm_cov_combined_tab, norm_filt_cov_combined_tab = process_each_table(input_cov_files, input_det_files, input_filt_cov_files, sample_names)

    # writing out
    cov_combined_tab.to_csv(args.output_prefix + "-genome-coverages.tsv", index = False, sep = "\t", na_rep = "NA")
    detection_combined_tab.to_csv(args.output_prefix + "-genome-detections.tsv", index = False, sep = "\t", na_rep = "NA")
    filtered_cov_combined_tab.to_csv(args.output_prefix + "-genome-detection-filtered-coverages.tsv", index = False, sep = "\t", na_rep = "NA")
    norm_cov_combined_tab.to_csv(args.output_prefix + "-CPM-normalized-genome-coverages.tsv", index = False, sep = "\t", na_rep = "NA")
    norm_filt_cov_combined_tab.to_csv(args.output_prefix + "-CPM-normalized-genome-detection-filtered-coverages.tsv", index = False, sep = "\t", na_rep = "NA")


################################################################################


# setting some colors
tty_colors = {
    'green' : '\033[0;32m%s\033[0m',
    'yellow' : '\033[0;33m%s\033[0m',
    'red' : '\033[0;31m%s\033[0m'
}


### functions ###
def color_text(text, color='green'):
    if sys.stdout.isatty():
        return tty_colors[color] % text
    else:
        return text


def wprint(text):
    """ print wrapper """

    print(textwrap.fill(text, width=80, initial_indent="  ",
          subsequent_indent="  ", break_on_hyphens=False))


def check_all_inputs_exist(input_cov_tables, input_det_tables, input_filtered_cov_tables):

    for file in input_cov_tables + input_det_tables + input_filtered_cov_tables:
        if not os.path.exists(file):
            print("")
            wprint(color_text("It seems the specified input file '" + str(file) + "' can't be found.", "yellow"))
            print("\nExiting for now.\n")
            sys.exit(1)


def setup_input_lists(input_cov_tables, input_det_tables, input_filtered_cov_tables):
    """ setting up input lists for file locations and sample names """

    # checking for all three to be sure they all appear in all 3 groups of input files
    input_cov_files = []
    cov_sample_names = []

    for sample in input_cov_tables:
        input_cov_files.append(sample)
        cov_sample_names.append(os.path.splitext(os.path.basename(sample))[0].replace("-genome-level-coverages", ""))


    input_det_files = []
    det_sample_names = []

    for sample in input_det_tables:
        input_det_files.append(sample)
        det_sample_names.append(os.path.splitext(os.path.basename(sample))[0].replace("-genome-level-detections", ""))


    input_filt_cov_files = []
    filt_cov_sample_names = []

    for sample in input_filtered_cov_tables:
        input_filt_cov_files.append(sample)
        filt_cov_sample_names.append(os.path.splitext(os.path.basename(sample))[0].replace("-genome-level-detection-filtered-coverages", ""))



    # making sure all 3 are equal and in same order
    if not cov_sample_names == det_sample_names or not cov_sample_names == filt_cov_sample_names:
            print("")
            wprint(color_text("It seems the 3 groups of specified input files either don't all hold the same files, or aren't in the same relative order.", "yellow"))
            print("\nExiting for now.\n")
            sys.exit(1)

    return(input_cov_files, input_det_files, input_filt_cov_files, filt_cov_sample_names)


def process_each_table(input_cov_files, input_det_files, input_filt_cov_files, sample_names):
    """ For each group of input files, reads in each table and creates combined table, also makes CPM normalzied of coverage tabs """

    cov_tabs = []
    det_tabs = []
    filt_cov_tabs = []
    cov_norm_tabs = []
    filt_cov_norm_tabs = []
    

    # iterator to access the same input file and sample name
    for i in range(len(sample_names)):

        cov_tab = pd.read_csv(input_cov_files[i], sep="\t")
        det_tab = pd.read_csv(input_det_files[i], sep="\t")
        filt_cov_tab = pd.read_csv(input_filt_cov_files[i], sep="\t")

        # making CPM normalized versions of coverage tables
        cov_norm_tab = cov_tab.copy()
        cov_norm_tab.coverage = cov_norm_tab.coverage / cov_norm_tab.coverage.sum() * 1000000

        filt_cov_norm_tab = filt_cov_tab.copy()
        filt_cov_norm_tab.coverage = filt_cov_norm_tab.coverage / filt_cov_norm_tab.coverage.sum() * 1000000


        # changing column names to be sample name
        cov_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)
        det_tab.rename(columns = {"detection":sample_names[i]}, inplace = True)
        filt_cov_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)
        cov_norm_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)
        filt_cov_norm_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)


        # adding to lists
        cov_tabs.append(cov_tab)
        det_tabs.append(det_tab)
        filt_cov_tabs.append(filt_cov_tab)
        cov_norm_tabs.append(cov_norm_tab)
        filt_cov_norm_tabs.append(filt_cov_norm_tab)

    # combining tables
    cov_combined_tab = pd.concat(cov_tabs, axis = 1).T.drop_duplicates().T
    det_combined_tab = pd.concat(det_tabs, axis = 1).T.drop_duplicates().T
    filt_cov_combined_tab = pd.concat(filt_cov_tabs, axis = 1).T.drop_duplicates().T
    norm_cov_combined_tab = pd.concat(cov_norm_tabs, axis = 1).T.drop_duplicates().T
    norm_filt_cov_combined_tab = pd.concat(filt_cov_norm_tabs, axis = 1).T.drop_duplicates().T


    return(cov_combined_tab, det_combined_tab, filt_cov_combined_tab, norm_cov_combined_tab, norm_filt_cov_combined_tab)


if __name__ == "__main__":
    main()
