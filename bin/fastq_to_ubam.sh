#!/bin/bash
set -euo pipefail

# Christopher Medway AWMGS
# generates uBAM from R1 and R2 reads for a given sample / lane

r1_fastq=$1
r2_fastq=$2
worklistId=$3
panel=$4
expectedInsertSize=$5

r1_prefix=$(echo "${r1_fastq//.fastq}")
r2_prefix=$(echo "${r2_fastq//.fastq}")




echo "converting fastq to ubam"

/share/apps/jre-distros/jre1.8.0_131/bin/java \
    -jar /share/apps/picard-tools-distros/picard-tools-2.18.5/picard.jar \
    FastqToSam \
    F1=/data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_"$laneId"_R1.fastq \
    F2=/data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_"$laneId"_R2.fastq \
    O=/data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_"$laneId"_unaligned.bam \
    QUALITY_FORMAT=Standard \
    READ_GROUP_NAME="$seqId"_"$sampleId"_"$laneId" \
    SAMPLE_NAME="$sampleId" \
    LIBRARY_NAME="$worklistId"_"$sampleId"_"$panel" \
    PLATFORM_UNIT="$seqId"_"$sampleId"_"$laneId" \
    PLATFORM="ILLUMINA" \
    SEQUENCING_CENTER="IMG" \
    PREDICTED_INSERT_SIZE="$expectedInsertSize" \
    SORT_ORDER=queryname \
    MAX_RECORDS_IN_RAM=2000000 \
    TMP_DIR=/state/partition1/tmpdir \
    RUN_DATE=`date +%s` \
    QUIET=true \
    VERBOSITY=ERROR

rm /data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_"$laneId"_R1.fastq
rm /data/results/$seqId/$panel/$sampleId/"$seqId"_"$sampleId"_"$laneId"_R2.fastq
