#!/bin/bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# marks duplicated for all BAM files that have been generated for a sample
# merges across lanes


seqId=$1
sampleId=$2

picard \
    MarkDuplicates \
    $(ls "$seqId"_"$sampleId"_*_aligned.bam | \sed 's/^/I=/' | tr '\n' ' ') \
    OUTPUT="$seqId"_"$sampleId"_rmdup.bam \
    METRICS_FILE="$seqId"_"$sampleId"_markDuplicatesMetrics.txt \
    CREATE_INDEX=true \
    MAX_RECORDS_IN_RAM=2000000 \
    VALIDATION_STRINGENCY=SILENT \
    TMP_DIR=/state/partition1/tmpdir \
    QUIET=true \
    VERBOSITY=ERROR
