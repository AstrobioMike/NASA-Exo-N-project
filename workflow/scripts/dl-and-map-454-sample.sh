#!/usr/bin/env bash
set -u

# takes 5 positional arguments:
  # 1. the mapping directory holding the bowtie2 index and where the output bam files will be stored
  # 2. the bowtie2 index basename to map to
  # 3. number of threads to tell bowtie2 to use
  # 4. the unique sample ID
  # 5. a comma-delimited list of the associated run accessions
    # e.g.:
# bash dl-and-map-sample.sh mapping-directory syn-bowtie2-db 5 ANE_004_05M ERR598955,ERR599003

mapping_dir=$1
bowtie2_db=$2
num_threads=$3
sample=$4
list_of_run_accs=$5

timepoint=$(date)
printf "Starting on ${sample} at timepoint:  ${timepoint}\n\n" > ${mapping_dir}/${sample}-dl-and-mapping.log

printf "#################\n## Downloading ##\n#################\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log

IFS=',' read -r -a run_accs <<< "$list_of_run_accs"

for acc in "${run_accs[@]}"
do
    printf "  Working on $acc\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log

    # keeping track of number of attempts to download sra object (this seems to be where the failure can happen, but only sometimes)
    try=0

    # starting a while loop to keep trying until the sra object downloads successfully (because sometimes they just wanna fail)
    while [ ! -f ${acc}-is-good-to-go.tmp ]
    do

        try=$((${try} + 1))

        printf "    Downloading sra object for $acc, attempt # ${try}\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log

        # writing stdout and stderr to a file specific to this ERR so can track it
        prefetch ${acc} --max-size 500G > ${acc}-dl-out.tmp 2>&1

        # after that command finishes, it either worked or didn't work, so checking here
        if grep -q "was downloaded successfully" ${acc}-dl-out.tmp # if it finished successfully...
        then
            echo "got it" > ${acc}-is-good-to-go.tmp # then writing something to this file so it exists and we move on from this while loop
            rm ${acc}-dl-out.tmp
        else
            printf "    Failed, trying again...\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log # if it didn't finish successfully, then this while loop will restart and try to download the sra object again and print to the screen that it is trying again
            rm -rf "$acc" # removing to try again fresh
        fi

    done

    rm ${acc}-is-good-to-go.tmp # removing the tmp file needed for the while loop now that moving on

    printf "    ${acc} downloaded successfully. Parsing reads now...\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
    fasterq-dump ${acc}/${acc}.sra >> ${mapping_dir}/${sample}-dl-and-mapping.log 2>&1

    # removing SRA object
    rm -rf ${acc}

done

printf "\n#############\n## Mapping ##\n#############\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log

# setting variables that hold all reads if stored in multiple ERRs/files, comma-delimited so can be passed to bowtie2 together
READS=$(echo $list_of_run_accs | sed 's/,/.sra.fastq,/g' | sed 's/$/.sra.fastq/')

# mapping command
map_command="bowtie2 --mm -q -x ${bowtie2_db} -U ${READS} -p ${num_threads} --no-unal 2> ${mapping_dir}/${sample}-mapping-info.txt | samtools view -b | samtools sort -@ ${num_threads} > ${mapping_dir}/${sample}.bam 2> /dev/null"
  # writing to log file
printf "Commands as called:\n${map_command}\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
  # executing
eval ${map_command}

# remaining index command
index_command="samtools index -@ ${num_threads} ${mapping_dir}/${sample}.bam"
  # writing to log file
printf "${index_command}\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
  # executing
eval ${index_command}

# removing files
rm $(echo $READS | tr "," " ")

printf "  Fun is fun, but this sample is done.\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
timepoint=$(date)
printf "    ${timepoint}\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
