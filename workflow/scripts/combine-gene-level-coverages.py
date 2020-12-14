#!/usr/bin/env python

"""
This is an ad hoc script for the corresponding workflow. It combines the individual sample gene-level coverage files and KO annotations into one table.
It produces 2 output files, one normalized to coverage-per-million, one not normalized.

Modified from my `bit-GL-combine-KO-and-tax-tables`: https://github.com/AstrobioMike/bioinf_tools
"""

import os
import sys
import argparse
import textwrap
import pandas as pd
from math import isnan
from numpy import NaN

parser = argparse.ArgumentParser(description="This is an ad hoc script for the corresponding workflow. It combines the individual sample gene-level coverage files and KO annotations into one table. \
                                              It produces 2 output files, one normalized to coverage-per-million, one not normalized.")

required = parser.add_argument_group('required arguments')

required.add_argument("input_tables", metavar="input-tables", type=str, nargs="+", help="Input gene-level coverage tables (as written, expected to end with extension '.tsv'.")
required.add_argument("-a", "--KO-annotations-file", help="Input KO annotation table", action="store", dest="KO_tab")
parser.add_argument("-o", "--output-prefix", help='Desired output prefix (default: "Combined")', action="store", default="Combined", dest="output_prefix")

if len(sys.argv)==1:
    parser.print_help(sys.stderr)
    sys.exit(0)

args = parser.parse_args()

################################################################################


def main():

    check_all_inputs_exist(args.input_tables, args.KO_tab)

    input_files, sample_names = setup_input_lists(args.input_tables)

    unnormd_combined_tab, normd_combined_tab = process_each_table(input_files, sample_names)

    unnormd_combined_with_annots_tab, normd_combined_with_annots_tab = add_KO_annots(unnormd_combined_tab, normd_combined_tab, args.KO_tab)

    unnormd_combined_with_annots_tab.to_csv(args.output_prefix + "-gene-coverages-and-KO-annots.tsv", index=False, sep="\t", na_rep = "NA")

    normd_combined_with_annots_tab.to_csv(args.output_prefix + "-CPM-normalized-gene-coverages-and-KO-annots.tsv", index=False, sep="\t", na_rep = "NA")


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


def check_all_inputs_exist(input_tables, KO_tab):

    for file in input_tables:
        if not os.path.exists(file):
            print("")
            wprint(color_text("It seems the specified input file '" + str(file) + "' can't be found.", "yellow"))
            print("\nExiting for now.\n")
            sys.exit(1)

    if not os.path.exists(KO_tab):
        print("")
        wprint(color_text("It seems the specified input file '" + str(KO_tab) + "' can't be found.", "yellow"))
        print("\nExiting for now.\n")
        sys.exit(1)


def setup_input_lists(input_tables):
    """ setting up input lists for file locations and sample names """

    input_files = []
    sample_names = []

    for sample in input_tables:
        input_files.append(sample)
        sample_names.append(os.path.splitext(os.path.basename(sample))[0].replace("-gene-coverages", ""))

    return(input_files, sample_names)


def process_each_table(input_files, sample_names):
    """ reads in each table, creates combined tables, one normalized to coverage-per-million, one not normalized  """

    normd_tabs = []
    unnormd_tabs = []

    # iterator to access the same input file and sample name
    for i in range(len(input_files)):

        unnormd_tab = pd.read_csv(input_files[i], sep="\t")

        # generating a normalized version
        normd_tab = unnormd_tab.copy()
        normd_tab.coverage = normd_tab.coverage / normd_tab.coverage.sum() * 1000000

        # changing coverage column headers to be sample name
        unnormd_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)
        normd_tab.rename(columns = {"coverage":sample_names[i]}, inplace = True)

        # adding to lists
        unnormd_tabs.append(unnormd_tab)
        normd_tabs.append(normd_tab)


    # combining tables
    unnormd_combined_tab = pd.concat(unnormd_tabs, axis = 1).T.drop_duplicates().T
    normd_combined_tab = pd.concat(normd_tabs, axis = 1).T.drop_duplicates().T

    return(unnormd_combined_tab, normd_combined_tab)

def add_KO_annots(unnormd_combined_tab, normd_combined_tab, KO_tab):
    """ adds KO annotations to tables """
    
    KO_tab = pd.read_csv(KO_tab, sep = "\t")

    unnormd_combined_with_annots_tab = KO_tab.merge(unnormd_combined_tab, on = "gene_ID")
    normd_combined_with_annots_tab = KO_tab.merge(normd_combined_tab, on = "gene_ID")

    return(unnormd_combined_with_annots_tab, normd_combined_with_annots_tab)

if __name__ == "__main__":
    main()
