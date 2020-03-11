
#!/bin/bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# runs BWA-mem over given sample: input unaligned BAM output aligned BAM

seqId=$1
laneId=$2
sequencing_centre=$3
sampleId=$4
reference_genome=$5
read1=$6
read2=$7

bwa mem \
    -v 1 \
    -t 12 \
    -p \
    -M \
    -R "@RG\\tID:"$seqId"."$laneId"\\tCN:"$sequencing_centre"\\tSM:"$sampleId"\\tLB:"$seqId"\\tPL:ILLUMINA" \
    $reference_genome \
    $read1 \
    $read2 | samtools view -Sb - | \
    samtools sort -T . -O bam > ""$seqId"_"$sampleId"_"$laneId".bam" 
    samtools index ""$seqId"_"$sampleId"_"$laneId".bam"
