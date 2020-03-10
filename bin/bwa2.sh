
#!/bin/bash
set -euo pipefail

# Christopher Medway, Seemu Ali AWMGS
# runs BWA-mem over given sample: input unaligned BAM output aligned BAM

r1_fastq=$1
r2_fastq=$2
reference_genome=$3
seqId=$4
laneId=$5
sequencing_centre=$6
sample_id=$7

general_Id=$(echo "${r1_fastq//_R1.fastq}")

	bwa mem \
    -t 12 \
    -p \
    -v 1 \
    -M \
    $reference_genome \
    $r1_fastq \
    $r2_fastq \
    > tes.sam 



    #| samtools view -Sb - | \
    #samtools sort -T /var/folders/59/6qr118qs6hv7gk0nz4lb7_6m0000gv/T/ -O bam > "$general_Id".bam \
    #samtools index "$general_Id".bam




