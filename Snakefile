############################################################################################
## Snakefile for Nitrogen Exobio metagenomics                                             ##
## Developed by Michael D. Lee (Mike.Lee@nasa.gov)                                        ##
############################################################################################

import os
import pandas as pd

########################################
############# General Info #############
########################################

"""
See the corresponding 'config.yaml' file for general use information. 
Variables that may need to be adjusted should be changed there, not here.
"""

###############################################################################
######## Pulling variables from config file for easier calling in here ########
###############################################################################
configfile: "config.yaml"

sample_IDs_file = config["sample_IDs_file"]
group_IDs_file = config["group_IDs_file"]
groups_map_file = config["groups_map_file"]

raw_R1_suffix = config["raw_R1_suffix"]
raw_R2_suffix = config["raw_R2_suffix"]
trimmed_R1_suffix = config["trimmed_R1_suffix"]
trimmed_R2_suffix = config["trimmed_R2_suffix"]

raw_reads_dir = config["raw_reads_dir"]
fastqc_out_dir = config["fastqc_out_dir"]
filtered_reads_dir = config["filtered_reads_dir"]
assemblies_dir = config["assemblies_dir"]
genes_dir = config["genes_dir"]
annotations_and_tax_dir = config["annotations_and_tax_dir"]
mapping_dir = config["mapping_dir"]
combined_output_dir = config["combined_output_dir"]
bins_dir = config["bins_dir"]
MAGs_dir = config["MAGs_dir"]

dirs_to_create = [fastqc_out_dir, filtered_reads_dir, assemblies_dir, genes_dir, 
                 annotations_and_tax_dir, mapping_dir, combined_output_dir, bins_dir, MAGs_dir]

num_threads = config["num_threads"]
num_cpus = config["num_cpus"]
gtdb_tk_num_cpus = config["gtdb_tk_num_cpus"]
max_mem = config["max_mem"]

###########################################
########### Ref database stuff ############
###########################################

REF_DB_ROOT_DIR = config["REF_DB_ROOT_DIR"]

# database triggers (for downloading and setting up if they aren't already present)
kofamscan_db_trigger = REF_DB_ROOT_DIR + config["KOFAMSCAN_DIR"] + "/KO_DB_SETUP"
cat_db_trigger = REF_DB_ROOT_DIR + config["CAT_DIR"] + "/CAT_DB_SETUP"
checkm_db_trigger = REF_DB_ROOT_DIR + config["CHECKM_DIR"] + "/CHECKM_DB_SETUP"
gtdbtk_db_trigger = REF_DB_ROOT_DIR + config["GTDB_DATA_PATH"] + "/GTDBTK_DB_SETUP"



#############################################
###### Reading info files into objects ######
#############################################

sample_ID_list = [line.strip() for line in open(sample_IDs_file)]

group_ID_list = [line.strip() for line in open(group_IDs_file)]

group_map_tab = pd.read_csv(groups_map_file, sep="\t", names = ["sample", "group"])


################################################################################
#### Helper functions and setting up variables to handle co-assembly groups ####
################################################################################

def group_R1_reads(target_group, group_map_tab = group_map_tab):
    target_samples_list = list(group_map_tab[group_map_tab.group.eq(target_group)]['sample'].values)
    return(expand(filtered_reads_dir + "{ID}" + trimmed_R1_suffix, ID = target_samples_list))

def group_R2_reads(target_group, group_map_tab = group_map_tab):
    target_samples_list = list(group_map_tab[group_map_tab.group.eq(target_group)]['sample'].values)
    return(expand(filtered_reads_dir + "{ID}" + trimmed_R2_suffix, ID = target_samples_list))

R1_group_dict, R2_group_dict = {}, {}
for group in group_ID_list:
    R1_group_dict[group] = group_R1_reads(group)
    R2_group_dict[group] = group_R2_reads(group)


group_sample_bams_list, group_bam_dict = [], {}
def group_sample_bams(target_group, group_map_tab = group_map_tab):
    target_samples_list = list(group_map_tab[group_map_tab.group.eq(target_group)]['sample'].values)
    return(expand(mapping_dir + "{ID}-" + group + ".bam", ID = target_samples_list))

for group in group_ID_list:
    group_sample_bams_list += group_sample_bams(group)
    group_bam_dict[group] = group_sample_bams(group)


group_sample_coverage_files_list, group_cov_dict = [], {}
def group_sample_coverages_files(target_group, group_map_tab = group_map_tab):
    target_samples_list = list(group_map_tab[group_map_tab.group.eq(target_group)]['sample'].values)
    return(expand(mapping_dir + "{ID}-" + group + "-gene-coverages.tsv", ID = target_samples_list))

for group in group_ID_list:
    group_sample_coverage_files_list += group_sample_coverages_files(group)
    group_cov_dict[group] = group_sample_coverages_files(group)

########################################
######## Setting up directories ########
########################################

for dir in dirs_to_create:
	try:
		os.mkdir(dir)
	except:
		pass

for dir in expand(bins_dir + "{group}", group = group_ID_list):
    try:
        os.mkdir(dir)
    except:
        pass

########################################
############# Rules start ##############
########################################

rule all:
    input:
        combined_output_dir + "All-combined-KO-function-coverages.tsv",
        combined_output_dir + "All-combined-KO-function-CPM-normalized-coverages.tsv",
        combined_output_dir + "All-combined-taxonomy-coverages.tsv",
        combined_output_dir + "All-combined-taxonomy-CPM-normalized-coverages.tsv",
        MAGs_dir + "gtdbtk-out",
        MAGs_dir + "MAGs-overview.tsv",
        expand(mapping_dir + "{group}-metabat-assembly-depth.tsv", group = group_ID_list),
        expand(combined_output_dir + "{group}-gene-coverages-annotations-and-tax.tsv", group = group_ID_list),
        expand(combined_output_dir + "{group}-CPM-normalized-gene-coverages-annotations-and-tax.tsv", group = group_ID_list),
        group_sample_coverage_files_list,
        group_sample_bams_list,
        expand(mapping_dir + "{group}-bowtie2-build.log", group = group_ID_list),
        expand(genes_dir + "{group}-genes.faa", group = group_ID_list),
        assemblies_dir + "assembly-summaries.tsv",
        expand(assemblies_dir + "{group}-coassembly.fa", group = group_ID_list),
        fastqc_out_dir + "raw_multiqc.html",
        fastqc_out_dir + "raw_multiqc_data.zip",
        fastqc_out_dir + "filtered_multiqc.html",
        fastqc_out_dir + "filtered_multiqc_data.zip"
    shell:
        """
        # copying log file to store with processing info
        cp .snakemake/log/$(ls -t .snakemake/log/ | head -n 1) snakemake-run.log
        """


rule generate_MAGs_overview_table:
    input:
        assembly_summaries = MAGs_dir + "MAG-assembly-summaries.tsv",
        checkm_results = expand(MAGs_dir + "{group}-MAGs-checkm-out.tsv", group = group_ID_list),
        gtdb_done_trigger = MAGs_dir + "gtdbtk-out"
    params:
        gtdb_results = MAGs_dir + "gtdbtk-out/gtdbtk.*.summary.tsv",
        checkm_tmp = MAGs_dir + "checkm-estimates.tmp",
        gtdb_tmp = MAGs_dir + "gtdb-taxonomies.tmp",
        checkm_w_header_tmp = MAGs_dir + "checkm-estimates-with-headers.tmp",
        gtdb_w_header_tmp = MAGs_dir + "gtdb-taxonomies-with-headers.tmp",
        overview_tmp = MAGs_dir + "MAGs-overview.tmp",
        overview_header_tmp = MAGs_dir + "MAGs-overview-header.tmp",
        overview_sorted_tmp = MAGs_dir + "MAGs-overview-sorted.tmp"
    output:
        MAGs_dir + "MAGs-overview.tsv"
    shell:
        """
        # making sure none of the intermediate files exist already
        rm -rf {params.checkm_tmp} {params.gtdb_tmp} {params.checkm_w_header_tmp} {params.gtdb_w_header_tmp} {params.overview_tmp} {params.overview_header_tmp} {params.overview_sorted_tmp}

        for MAG in $(cut -f 1 {input.assembly_summaries} | tail -n +2)
        do

            MAG_base=$(echo ${{MAG}} | cut -f 2- -d "-")

            grep -w -m 1 "^${{MAG_base}}" {input.checkm_results} | cut -f 12,13 >> {params.checkm_tmp}

            # note: below won't work with darwin sed due to tab character
            grep -w "^${{MAG}}" {params.gtdb_results} | cut -f 2 | sed 's/^.__//' | sed 's/;.__/\\t/g' | awk -F $'\\t' ' BEGIN {{ OFS=FS }} {{ for (i=1; i<=NF; i++) if ( $i ~/^ *$/) $i = "NA" }}; 1 ' >> {params.gtdb_tmp}

        done

        # adding headers
        cat <(printf "est. completeness\\test. redundancy\\n") {params.checkm_tmp} > {params.checkm_w_header_tmp}
        cat <(printf "domain\\tphylum\\tclass\\torder\\tfamily\\tgenus\\tspecies\\n") {params.gtdb_tmp} > {params.gtdb_w_header_tmp}

        paste {input.assembly_summaries} {params.checkm_w_header_tmp} {params.gtdb_w_header_tmp} > {params.overview_tmp}

        # ordering by taxonomy
        head -n 1 {params.overview_tmp} > {params.overview_header_tmp}
        tail -n +2 {params.overview_tmp} | sort -t $'\\t' -k 13,19 > {params.overview_sorted_tmp}

        cat {params.overview_header_tmp} {params.overview_sorted_tmp} > {output}

        rm -rf {params.checkm_tmp} {params.gtdb_tmp} {params.checkm_w_header_tmp} {params.gtdb_w_header_tmp} {params.overview_tmp} {params.overview_header_tmp} {params.overview_sorted_tmp}
        """


rule summarize_MAG_assemblies:
   """ summarize MAG assemblies """

    conda:
        "envs/bit.yaml"
    input:
        trigger = expand(MAGs_dir + "{group}-MAGs-checkm-out.tsv", group = group_ID_list)
    params:
        intermediate_file = MAGs_dir + "MAGs-summaries.tmp"
    output:
        MAGs_dir + "MAG-assembly-summaries.tsv"
    shell:
        """
        bit-summarize-assembly {MAGs_dir}*/*.fa -o {params.intermediate_file} -t

        # slimming down the output
        cut -f 1,2,3,5,6,8,11,18,19,20 {params.intermediate_file} > {output}
        rm {params.intermediate_file}
        """


rule run_gtdbtk_on_MAGs:
    """ assign taxonomy to MAGs with gtdb-tk """

    conda:
        "envs/gtdb-tk.yaml"
    input:
        trigger = expand(MAGs_dir + "{group}-MAGs-checkm-out.tsv", group = group_ID_list),
        db_trigger = gtdbtk_db_trigger
    params:
        temp_all_MAGs_dir = MAGs_dir + "all-for-gtdb/"
    output:
        directory(MAGs_dir + "gtdbtk-out")
    log:
        MAGs_dir + "gtdbtk-out/gtdbtk-run.log"
    shell:
        """
        mkdir -p {params.temp_all_MAGs_dir}

        cp {MAGs_dir}*/*.fa {params.temp_all_MAGs_dir}

        gtdbtk classify_wf --genome_dir {params.temp_all_MAGs_dir} -x fa --out_dir {output} --cpus {gtdb_tk_num_cpus} > {log} 2>&1

        rm -rf {params.temp_all_MAGs_dir}
        """


rule filtering_checkm_results_and_copying_MAGs:
    input:
        bins_dir + "{group}-bins-checkm-out.tsv"
    output:
        MAGs_dir + "{group}-MAGs-checkm-out.tsv"
    shell:
        """
        mkdir -p {MAGs_dir}{wildcards.group}/

        cat <( head -n 1 {input} ) <( awk -F $'\\t' ' $12 >= 90 && $13 <= 10 && $14 == 0 ' {input} ) > {wildcards.group}-MAG-info.tmp

        sed 's/bin./{wildcards.group}-MAG-/' {wildcards.group}-MAG-info.tmp > {output}

        for MAG in $(cut -f 1 {wildcards.group}-MAG-info.tmp | tail -n +2)
        do
            new_ID=$(echo $MAG | sed 's/bin./MAG-/')
            cp {bins_dir}{wildcards.group}/${{MAG}}.fa {MAGs_dir}{wildcards.group}/{wildcards.group}-${{new_ID}}.fa
        done

        rm {wildcards.group}-MAG-info.tmp
        """


rule run_checkm_on_bins:
    """ runs checkm on recovered bins """

    conda:
        "envs/checkm.yaml"
    input:
        curr_bins_dir = bins_dir + "{group}/",
        trigger = mapping_dir + "{group}-metabat-assembly-depth.tsv",
        db_trigger = checkm_db_trigger
    params:
        tmp_dir = bins_dir + "{group}/checkm-out-tmp/"
    output:
        bins_dir + "{group}-bins-checkm-out.tsv"
    log:
        bins_dir + "{group}-checkm.log"
    shell:
        """
        checkm lineage_wf -f {output} --tab_table -t {num_cpus} --pplacer_threads {num_threads} -x fa {input.curr_bins_dir} {params.tmp_dir} > {log} 2>&1
        rm -rf {params.tmp_dir}
        """


rule metabat_binning:
    """
    This rule runs metabat2 for binning contigs in each of the groups
    """

    conda:
        "envs/metabat.yaml"
    input:
        assembly = assemblies_dir + "{group}-coassembly.fa",
        bams = lambda wildcards: group_bam_dict[wildcards.group]
    params:
        prefix = bins_dir + "{group}"
    output:
        depth_file = mapping_dir + "{group}-metabat-assembly-depth.tsv"
    shell:
        """
        jgi_summarize_bam_contig_depths --outputDepth {output.depth_file} --percentIdentity 97 --minContigLength 1000 --minContigDepth 1.0  --referenceFasta {input.assembly} {input.bams} > /dev/null 2>&1
        metabat2  --inFile {input.assembly} --outFile {params.prefix}/bin --abdFile {output.depth_file} -t {num_threads} > /dev/null 2>&1
        """

rule combine_all_gene_coverages_collapsed_by_KO_and_tax:
    input:
        expand(combined_output_dir + "{group}-gene-coverages-annotations-and-tax.tsv", group = group_ID_list)
    output:
        combined_output_dir + "All-combined-KO-function-coverages.tsv",
        combined_output_dir + "All-combined-KO-function-CPM-normalized-coverages.tsv",
        combined_output_dir + "All-combined-taxonomy-coverages.tsv",
        combined_output_dir + "All-combined-taxonomy-CPM-normalized-coverages.tsv"
    shell:
        """
        python scripts/combine-all-gene-tables.py -o {combined_output_dir} {input}
        """

rule combine_group_gene_coverage_annots_and_tax:
    input:
        coverage_files = lambda wildcards: group_cov_dict[wildcards.group],
        ko_tab = annotations_and_tax_dir + "{group}-annotations.tsv",
        tax_tab = annotations_and_tax_dir + "{group}-gene-tax.tsv"
    output:
        combined_output_dir + "{group}-gene-coverages-annotations-and-tax.tsv",
        combined_output_dir + "{group}-CPM-normalized-gene-coverages-annotations-and-tax.tsv"
    shell:
        """
        python scripts/combine-gene-level-coverages-annots-and-tax-per-group.py -g {wildcards.group} -a {input.ko_tab} -t {input.tax_tab} -o {combined_output_dir} {input.coverage_files}
        """


rule get_cov_and_det:
    """
    This rule pulls out coverage and detection information for each sample, gene-level and contig-level,
    and filters the gene-level coverage information based on requiring at least 50% detection.
    """

    conda:
        "envs/mapping.yaml"
    input:
        bam = mapping_dir + "{ID}-{group}.bam",
        nt = genes_dir + "{group}-genes.fa"
    params:
        gene_cov_and_det_tmp = mapping_dir + "{ID}-{group}-gene-cov-and-det.tmp",
        contig_cov_and_det_tmp = mapping_dir + "{ID}-{group}-contig-cov-and-det.tmp",
        gene_cov_tmp = mapping_dir + "{ID}-{group}-gene-cov.tmp",
        contig_cov_tmp = mapping_dir + "{ID}-{group}-contig-cov.tmp"
    output:
        gene_covs = mapping_dir + "{ID}-{group}-gene-coverages.tsv",
        contig_covs = mapping_dir + "{ID}-{group}-contig-coverages.tsv"
    shell:
        """
        pileup.sh -in {input.bam} fastaorf={input.nt} outorf={params.gene_cov_and_det_tmp} out={params.contig_cov_and_det_tmp} > /dev/null 2>&1

        # filtering coverages based on detection
          # genes
        grep -v "#" {params.gene_cov_and_det_tmp} | awk -F $'\t' ' BEGIN {{OFS=FS}} {{ if ( $10 <= 0.5 ) $4 = 0 }} {{ print $1,$4 }} ' > {params.gene_cov_tmp}
        cat <( printf "gene_ID\tcoverage\n" ) {params.gene_cov_tmp} > {output.gene_covs}

          # contigs
        grep -v "#" {params.contig_cov_and_det_tmp} | awk -F $'\t' ' BEGIN {{OFS=FS}} {{ if ( $5 <= 50 ) $2 = 0 }} {{ print $1,$2 }} ' > {params.contig_cov_tmp}
        cat <( printf "contig_ID\tcoverage\n" ) {params.contig_cov_tmp} > {output.contig_covs}

          # removing intermediate files
        rm {params}
        """


rule run_mapping:
    """
    This rule runs the mapping for each sample to its appropriate group coassembly.
    """
    
    conda:
        "envs/mapping.yaml"
    input:
        mapping_trigger = mapping_dir + "{group}-bowtie2-build.log",
        R1 = filtered_reads_dir + "{ID}" + trimmed_R1_suffix,
        R2 = filtered_reads_dir + "{ID}" + trimmed_R2_suffix
    params:
        index = mapping_dir + "{group}-index"
    output:
        mapping_dir + "{ID}-{group}.bam"
    log:
        mapping_dir + "{ID}-{group}-mapping-info.txt"
    shell:
        """
        bowtie2 --mm -q --threads {num_threads} -x {params.index} -1 {input.R1} -2 {input.R2} --no-unal 2> {log} | samtools view -b | samtools sort -@ {num_threads} > {output} 2> /dev/null
        samtools index -@ {num_threads} {output}
        """


rule build_bowtie2_index:
    """ Builds bowtie2 databases """

    conda:
        "envs/mapping.yaml"
    input:
        assemblies_dir + "{group}-coassembly.fa"
    params:
        basename = mapping_dir + "{group}-index"
    output:
        mapping_trigger = mapping_dir + "{group}-bowtie2-build.log"
    shell:
        """
        bowtie2-build {input} {params.basename} > {output.mapping_trigger} 2>&1
        """


rule run_tax_classification:
    """
    This rule runs the gene- and contig-level taxonomic classifications for each coassembly.
    """

    conda:
        "envs/cat.yaml"
    input:
        assembly = assemblies_dir + "{group}-coassembly.fa",
        AA = genes_dir + "{group}-genes.faa",
        db_trigger = cat_db_trigger
    output:
        gene_tax_out = annotations_and_tax_dir + "{group}-gene-tax.tsv",
        contig_tax_out = annotations_and_tax_dir + "{group}-contig-tax.tsv"
    params:
        tmp_out_prefix = annotations_and_tax_dir + "{group}-tax-out.tmp",
        tmp_genes = annotations_and_tax_dir + "{group}-gene-tax.tmp",
        tmp_contigs = annotations_and_tax_dir + "{group}-contig-tax.tmp",
        cat_db = REF_DB_ROOT_DIR + config["CAT_DIR"] + config["CAT_DB"],
        cat_tax = REF_DB_ROOT_DIR + config["CAT_DIR"] + config["CAT_TAX"]
    log:
        annotations_and_tax_dir + "{group}-CAT.log"
    shell:
        """
        CAT contigs -d {params.cat_db} -t {params.cat_tax} -n {num_cpus} -r 3 --top 4 --I_know_what_Im_doing -c {input.assembly} -p {input.AA} -o {params.tmp_out_prefix} --force > {log} 2>&1

        # adding names to gene classifications
        CAT add_names -i {params.tmp_out_prefix}.ORF2LCA.txt -o {params.tmp_genes} -t {params.cat_tax} --only_official > {log} 2>&1

        # formatting gene classifications
        awk -F $'\t' ' BEGIN {{ OFS=FS }} {{ if ( $2 == "lineage" ) {{ print $1,$2,$4,$5,$6,$7,$8,$9,$10 }} \
        else if ( $2 == "ORF has no hit to database" || $2 ~ /^no taxid found/ ) {{ print $1,"NA","NA","NA","NA","NA","NA","NA","NA" }} \
        else {{ n=split($2,lineage,";"); print $1,lineage[n],$4,$5,$6,$7,$8,$9,$10 }} }} ' {params.tmp_genes} | \
        sed 's/not classified/NA/g' | sed 's/superkingdom/domain/' | sed 's/^# ORF/gene_ID/' | sed 's/lineage/taxid/' | \
        sed 's/\*//g' > {output.gene_tax_out}

        # adding names to contig classifications
        CAT add_names -i {params.tmp_out_prefix}.contig2classification.txt -o {params.tmp_contigs} -t {params.cat_tax} --only_official > {log} 2>&1

        # formatting contig classifications
        awk -F $'\t' ' BEGIN {{ OFS=FS }} {{ if ( $2 == "classification" ) {{ print $1,$4,$6,$7,$8,$9,$10,$11,$12 }} \
        else if ( $2 == "unclassified" ) {{ print $1,"NA","NA","NA","NA","NA","NA","NA","NA" }} \
        else {{ n=split($4,lineage,";"); print $1,lineage[n],$6,$7,$8,$9,$10,$11,$12 }} }} ' {params.tmp_contigs} | \
        sed 's/not classified/NA/g' | sed 's/superkingdom/domain/' | sed 's/: [0-9\.]*//g' | sed 's/^# contig/contig_ID/' | \
        sed 's/lineage/taxid/' | sed 's/\*//g' > {output.contig_tax_out}

#        rm -rf {annotations_and_tax_dir}{wildcards.group}*tmp*
        rm -rf {params}
        """


rule run_KO_annotation:
    """
    This rule runs the gene-level (KO) functional annotation for each coassembly.
    """

    conda:
        "envs/kofamscan.yaml"
    input:
        AAs = genes_dir + "{group}-genes.faa",
        db_trigger = kofamscan_db_trigger
    output:
        annotations_and_tax_dir + "{group}-annotations.tsv"
    params:
        ko_db_dir = REF_DB_ROOT_DIR + config["KOFAMSCAN_DIR"],
        tmp_out = annotations_and_tax_dir + "{group}-KO-tab.tmp",
        tmp_dir = annotations_and_tax_dir + "{group}-tmp-KO-dir"
    shell:
        """
        exec_annotation -p {params.ko_db_dir}/profiles/ -k {params.ko_db_dir}/ko_list --cpu {num_cpus} -f detail-tsv -o {params.tmp_out} --tmp-dir {params.tmp_dir} --report-unannotated {input.AAs}

        bit-filter-KOFamScan-results -i {params.tmp_out} -o {output}

        rm -rf {params.tmp_out} {params.tmp_dir}
        """


rule call_genes:
    """
    This rule calls genes on each group's coassembly file.
    """

    conda:
        "envs/prodigal.yaml"
    input:
        assemblies_dir + "{group}-coassembly.fa"
    output:
        AA = genes_dir + "{group}-genes.faa",
        nt = genes_dir + "{group}-genes.fa",
        gff = genes_dir + "{group}-genes.gff"
    shell:
        """
        prodigal -q -c -p meta -a {output.AA} -d {output.nt} -f gff -o {output.gff} -i {input}
        """


rule summarize_assemblies:
    """
    This rule summarizes and reports general stats for all individual sample assemblies in one table.
    """

    conda:
        "envs/bit.yaml"
    input:
        expand(assemblies_dir + "{group}-coassembly.fa", group = group_ID_list),
    output:
        assemblies_dir + "assembly-summaries.tsv"
    shell:
        """
        bit-summarize-assembly -o {output} {input}
        """


rule assembly:
    """
    This rule handles running the co-assembly for each individual group of samples.
    """

    conda:
        "envs/megahit.yaml"
    input:
        R1_reads = lambda wildcards: R1_group_dict[wildcards.group],
        R2_reads = lambda wildcards: R2_group_dict[wildcards.group]
    output:
        assemblies_dir + "{group}-coassembly.fa"
    log:
        assemblies_dir + "{group}-coassembly.log"
    shell:
        """
        # removing output directory if exists already but rule still needs to be run (because there is no --force option to megahit i dont't think):
        rm -rf {assemblies_dir}{wildcards.group}-megahit-out/

        R1_reads=$(echo {input.R1_reads} | tr " " ",")
        R2_reads=$(echo {input.R2_reads} | tr " " ",")

        megahit -1 ${{R1_reads}} -2 ${{R2_reads}} -m {max_mem} -t {num_threads} --min-contig-len 500 -o {assemblies_dir}{wildcards.group}-megahit-out > {log} 2>&1
        bit-rename-fasta-headers -i {assemblies_dir}{wildcards.group}-megahit-out/final.contigs.fa -w c_{wildcards.group} -o {output}

        # getting rid of assembly directory to save space
        rm -rf {assemblies_dir}{wildcards.group}-megahit-out/
        """

rule filtered_fastqc:
    """
    This rule runs fastqc on all trimmed/filtered input fastq files.
    """

    conda:
        "envs/qc.yaml"
    input:
        filtered_reads_dir + "{ID}" + trimmed_R1_suffix,
        filtered_reads_dir + "{ID}" + trimmed_R2_suffix
    output:
        filtered_reads_dir + "{ID}-R1-trimmed_fastqc.zip",
        filtered_reads_dir + "{ID}-R2-trimmed_fastqc.zip"
    shell:
        """
		fastqc {input} -t {num_threads} -q
		"""


rule filtered_multiqc:
    """
    This rule collates all trimmed/filtered fastqc outputs.
    """

    conda:
        "envs/qc.yaml"
    input:
        expand(filtered_reads_dir + "{ID}-R1-trimmed_fastqc.zip", ID = sample_ID_list),
        expand(filtered_reads_dir + "{ID}-R2-trimmed_fastqc.zip", ID = sample_ID_list)
    output:
        fastqc_out_dir + "filtered_multiqc.html",
        fastqc_out_dir + "filtered_multiqc_data.zip"
    shell:
        """
        multiqc -z -q -o {fastqc_out_dir} -n filtered_multiqc  {filtered_reads_dir} > /dev/null 2>&1
          # removing the individual fastqc files and temp locations
        rm -rf {filtered_reads_dir}*fastqc*
        """


rule bbduk:
    """
    This rule runs quality filtering/trimming on raw input fastq files for each individual sample.
    """

    conda:
        "envs/qc.yaml"
    input:
        in1 = raw_reads_dir + "{ID}" + raw_R1_suffix,
        in2 = raw_reads_dir + "{ID}" + raw_R2_suffix
    output:
        out1 = filtered_reads_dir + "{ID}" + trimmed_R1_suffix,
        out2 = filtered_reads_dir + "{ID}" + trimmed_R2_suffix
    log:
        filtered_reads_dir + "bbduk-{ID}.log"
    shell:
        """
        bbduk.sh in={input.in1} in2={input.in2} out1={output.out1} out2={output.out2} \
                ref=${{CONDA_PREFIX}}/opt/bbmap-38.86-0/resources/adapters.fa ktrim=l k=17 ftm=5 qtrim=rl \
                trimq=10 mlf=0.5 maxns=0 > {log} 2>&1
        """


rule raw_multiqc:
    """
    This rule collates all raw fastqc outputs.
    """

    conda:
        "envs/qc.yaml"
    input:
        expand(raw_reads_dir + "{ID}-R1_fastqc.zip", ID = sample_ID_list),
        expand(raw_reads_dir + "{ID}-R2_fastqc.zip", ID = sample_ID_list)
    output:
        fastqc_out_dir + "raw_multiqc.html",
        fastqc_out_dir + "raw_multiqc_data.zip"
    shell:
        """
        multiqc -z -q -o {fastqc_out_dir} -n raw_multiqc {raw_reads_dir} > /dev/null 2>&1
          # removing the individual fastqc files
        rm -rf {raw_reads_dir}*fastqc*
        """


rule raw_fastqc:
    """
    This rule runs fastqc on all raw input fastq files.
    """

    conda:
        "envs/qc.yaml"
    input:
        raw_reads_dir + "{ID}" + raw_R1_suffix,
        raw_reads_dir + "{ID}" + raw_R2_suffix        
    output:
        raw_reads_dir + "{ID}-R1_fastqc.zip",
        raw_reads_dir + "{ID}-R2_fastqc.zip"
    shell:
        """
		fastqc {input} -t {num_threads} -q
		"""

### database checking and setup rules ###
rule setup_CAT_db:
    """
    This rule checks for the CAT reference database, and downloads if needed.
    """

    conda:
        "envs/cat.yaml"
    output:
        trigger = cat_db_trigger
    params:
        cat_db_dir = REF_DB_ROOT_DIR + config["CAT_DIR"],
        compressed_cat = REF_DB_ROOT_DIR + "CAT_prepare_20200618.tar.gz"
    shell:
        """
        mkdir -p {REF_DB_ROOT_DIR}

        curl -L -o {params.compressed_cat} https://tbb.bio.uu.nl/bastiaan/CAT_prepare/CAT_prepare_20201123.tar.gz > /dev/null 2>&1
        tar -xvzf {params.compressed_cat} -C {REF_DB_ROOT_DIR} > /dev/null 2>&1
        rm {params.compressed_cat}

        touch {output.trigger}
        """


rule setup_KoFamScan_db:
    """
    This rule checks for the KoFamScan db (minimally currently) and downloads if needed.
    """

    conda:
        "envs/kofamscan.yaml"
    output:
        trigger = kofamscan_db_trigger
    params:
        ko_db_dir = REF_DB_ROOT_DIR + config["KOFAMSCAN_DIR"],
        compressed_ko_list = REF_DB_ROOT_DIR + config["KOFAMSCAN_DIR"] + "/ko_list.gz",
        compressed_profiles = REF_DB_ROOT_DIR + config["KOFAMSCAN_DIR"] + "/profiles.tar.gz"
    shell:
        """
        mkdir -p {params.ko_db_dir}

        curl -L -o {params.compressed_ko_list} ftp://ftp.genome.jp/pub/db/kofam/ko_list.gz > /dev/null 2>&1
        curl -L -o {params.compressed_profiles} ftp://ftp.genome.jp/pub/db/kofam/profiles.tar.gz > /dev/null 2>&1
        tar -xzf {params.compressed_profiles} -C {params.ko_db_dir} > /dev/null 2>&1
        gunzip {params.compressed_ko_list}

        touch {output.trigger}
        """


rule setup_checkm_db:
    """
    This rule checks for the checkm db (minimally currently) and downloads if needed and sets location for program.
    """

    conda:
        "envs/checkm.yaml"
    output:
        trigger = checkm_db_trigger
    params:
        checkm_db_dir = REF_DB_ROOT_DIR + config["CHECKM_DIR"]
    shell:
        """
        mkdir -p {params.checkm_db_dir}
        cd {params.checkm_db_dir}

        curl -LO https://data.ace.uq.edu.au/public/CheckM_databases/checkm_data_2015_01_16.tar.gz > /dev/null 2>&1
        tar -xzvf checkm_data_2015_01_16.tar.gz > /dev/null 2>&1

        checkm data setRoot $PWD > /dev/null 2>&1

        cd - > /dev/null 2>&1

        touch {output.trigger}
        """


rule setup_gtdbtk_db:
    """
    This rule checks for the gtdb-tk db (minimally currently) and downloads if needed.
    """

    conda:
        "envs/gtdb-tk.yaml"
    output:
        trigger = gtdbtk_db_trigger
    params:
        gtdbtk_db_dir = REF_DB_ROOT_DIR + config["GTDB_DATA_PATH"]
    shell:
        """
        mkdir -p {params.gtdbtk_db_dir}
        cd {params.gtdbtk_db_dir}

        # adding wanted location to this conda env PATH (gtdb-tk looks in the GTDBTK_DATA_PATH variable),
            # so will be set when the conda environment is started from now on
        mkdir -p ${{CONDA_PREFIX}}/etc/conda/activate.d/
        echo 'export GTDBTK_DATA_PATH={params.gtdbtk_db_dir}' >> ${{CONDA_PREFIX}}/etc/conda/activate.d/set_env_vars.sh

        # but still needs to be set for this particular session that is downloading the db
        GTDBTK_DATA_PATH={params.gtdbtk_db_dir}

        # now downloading
        download-db.sh > /dev/null 2>&1

        cd - > /dev/null 2>&1

        touch {output.trigger}
        """


rule clean_some:
    shell:
        "rm -rf {fastqc_out_dir} {filtered_reads_dir} {assemblies_dir} {genes_dir} {mapping_dir} {combined_output_dir}"


rule clean_all:
    shell:
        "rm -rf {dirs_to_create} .snakemake/ snakemake-run.log"