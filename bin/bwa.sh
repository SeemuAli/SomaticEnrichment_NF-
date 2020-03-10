#!/bin/bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# runs BWA-mem over given sample: input unaligned BAM output aligned BAM

unaligned_bam=$1
bam_prefix=$(echo "${unaligned_bam//_unaligned.bam}")
reference_genome=$2


    picard \
    SamToFastq \
    I="$unaligned_bam" \
    FASTQ=fastqstdout \
    INTERLEAVE=true \
    NON_PF=true \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/var/folders/59/6qr118qs6hv7gk0nz4lb7_6m0000gv/T/ \
    COMPRESSION_LEVEL=0 \
    QUIET=true \
    VERBOSITY=ERROR \
    | \
    bwa mem \
    -M \
    -t 12 \
    -p \
    -v 1 \
    "$reference_genome"\
    /dev/stdin \
    | \
    picard \
    MergeBamAlignment \
    EXPECTED_ORIENTATIONS=FR \
    ALIGNED_BAM=/dev/stdin \
    UNMAPPED_BAM="$unaligned_bam"\
    OUTPUT="$bam_prefix"_aligned.bam \
    REFERENCE_SEQUENCE="$reference_genome" \
    PAIRED_RUN=true \
    SORT_ORDER="coordinate" \
    CLIP_ADAPTERS=false \
    MAX_RECORDS_IN_RAM=2000000 \
    MAX_INSERTIONS_OR_DELETIONS=-1 \
    PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
    CREATE_INDEX=true \
    QUIET=true \
    VERBOSITY=ERROR \
    TMP_DIR=/var/folders/59/6qr118qs6hv7gk0nz4lb7_6m0000gv/T/ 

        
