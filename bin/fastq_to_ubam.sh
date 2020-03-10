#!/bin/bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# generates uBAM from R1 and R2 reads for a given sample / lane

r1_fastq=$1
r2_fastq=$2
sampleId=$3
worklistId=$4
panel=$5
expectedInsertSize=$6
general_Id=$(echo "${r1_fastq//_R1.fastq}")
r1_prefix=$(echo "${r1_fastq//.fastq}")
r2_prefix=$(echo "${r2_fastq//.fastq}")



picard \
    FastqToSam \
    F1="$r1_fastq" \
    F2="$r2_fastq"\
    O="$general_Id"_unaligned.bam \
    READ_GROUP_NAME="$general_Id" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId"_"$panel" 
    PLATFORM_UNIT="$general_Id" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="IMG" \
    PREDICTED_INSERT_SIZE="$expectedInsertSize" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/var/folders/59/6qr118qs6hv7gk0nz4lb7_6m0000gv/T/  \
    RUN_DATE=`date +%s` \
    QUIET=true \
    VERBOSITY=ERROR \
    QUALITY_FORMAT=Standard