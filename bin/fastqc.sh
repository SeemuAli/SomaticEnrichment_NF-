#!/usr/bin/env bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# runs FASTQC for R1 & R2 for given sample/lane

r1_fastq=$1
r2_fastq=$2
r1_prefix=$(echo "${r1_fastq//.fastq}")
r2_prefix=$(echo "${r2_fastq//.fastq}")

mkdir -p FASTQC

# consider adding --adapter to command
fastqc \
    --dir . \
    --threads 12 \
    --extract \
    --quiet \
    $r1_fastq \
    $r2_fastq

mv "$r1_prefix"_fastqc/summary.txt "$r1_prefix"_fastqc.txt
mv "$r2_prefix"_fastqc/summary.txt "$r2_prefix"_fastqc.txt
mv *fastqc* ./FASTQC/