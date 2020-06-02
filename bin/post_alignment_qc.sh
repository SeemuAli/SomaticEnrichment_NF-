#!/bin/bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# generation of PICARD metrics for each sample 

seqId=$1
sampleId=$2
panel=$3
minimumCoverage=$4
vendorCaptureBed=$5
vendorPrimaryBed=$6
padding=$7
minBQS=$8
minMQS=$9
reference_genome=${10}
refdict=${11}

echo $reference_genome


#Convert capture BED to interval_list for later
picard BedToIntervalList \
    I=$vendorCaptureBed \
    O="$panel"_capture.interval_list \
    SD=$refdict

#Convert primary BED to interval_list for later
picard BedToIntervalList \
    I=$vendorPrimaryBed \
    O="$panel"_primary.interval_list \
    SD=$refdict

#Alignment metrics: library sequence similarity
picard CollectAlignmentSummaryMetrics \
    R=$reference_genome \
    I="$seqId"_"$sampleId".bam \
    O="$seqId"_"$sampleId"_AlignmentSummaryMetrics.txt \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR= . 

#Calculate insert size: fragmentation performance
picard CollectInsertSizeMetrics \
    I="$seqId"_"$sampleId".bam \
    O="$seqId"_"$sampleId"_InsertMetrics.txt \
    H="$seqId"_"$sampleId"_InsertMetrics.pdf \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR= .

#HsMetrics: capture & pooling performance
picard CollectHsMetrics \
     I="$seqId"_"$sampleId".bam \
     O="$seqId"_"$sampleId"_HsMetrics.txt \
     R=$reference_genome \
     BAIT_INTERVALS="$panel"_capture.interval_list \
     TARGET_INTERVALS="$panel"_primary.interval_list \
     MAX_RECORDS_IN_RAM=2000000 \
     TMP_DIR= . \
     MINIMUM_MAPPING_QUALITY=$minMQS \
     MINIMUM_BASE_QUALITY=$minBQS \
     CLIP_OVERLAPPING_READS=false



