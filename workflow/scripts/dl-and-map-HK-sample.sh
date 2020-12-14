#!/usr/bin/env bash
set -u

# takes 6 positional arguments:
  # 1. the mapping directory holding the bowtie2 index and where the output bam files will be stored
  # 2. the bowtie2 index basename to map to
  # 3. number of threads to tell bowtie2 to use
  # 4. the unique sample ID
  # 5. forward reads link
  # 6. reverse reads link
    # e.g.:
# bash dl-and-map-HK-sample.sh mapping-directory syn-bowtie2-db 5 AT04-0m-8Ar ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/037/SRR10912837/SRR10912837.fastq.gz ftp.sra.ebi.ac.uk/vol1/fastq/SRR109/036/SRR10912836/SRR10912836.fastq.gz

mapping_dir=$1
bowtie2_db=$2
num_threads=$3
sample=$4
fwd_reads_link=$5
rev_reads_link=$6

timepoint=$(date)
printf "Starting on ${sample} at timepoint:  ${timepoint}\n\n" > ${mapping_dir}/${sample}-dl-and-mapping.log

printf "#################\n## Downloading ##\n#################\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log

curl -s -L -o ${sample}-R1.fq.gz ${fwd_reads_link} > /dev/null
curl -s -L -o ${sample}-R2.fq.gz ${rev_reads_link} > /dev/null

printf "\n#############\n## Mapping ##\n#############\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log

# mapping command
map_command="bowtie2 --mm -q -x ${bowtie2_db} -1 ${sample}-R1.fq.gz -2 ${sample}-R2.fq.gz -p ${num_threads} --no-unal 2> ${mapping_dir}/${sample}-mapping-info.txt | samtools view -b | samtools sort -@ ${num_threads} > ${mapping_dir}/${sample}.bam 2> /dev/null"
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
rm ${sample}-R1.fq.gz ${sample}-R2.fq.gz

printf "  Fun is fun, but this sample is done.\n\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
timepoint=$(date)
printf "    ${timepoint}\n" >> ${mapping_dir}/${sample}-dl-and-mapping.log
