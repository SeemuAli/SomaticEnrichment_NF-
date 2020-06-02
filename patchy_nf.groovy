
process calculating_coverage {

    input: 

    file(bam) from bam_cov
    file(bai) from bai_cov

    output: 
    vendorCaptureBed_100pad.bed

"$seqId"_"$sampleId"_DepthOfCoverage
"""
hotspot_coverage.sh \
    ${params.seqId} \
    ${params.sampleId} \
    ${params.panel} \
    ${params.pipelinename} \
    ${params.pipelineversion} \
    ${params.minimiumCoverage} \
    $capture_bed \
    ${params.padding} \
    ${params.minBQS} \
    ${params.minMQS} \
    $bedtools_genome \
    $reference_genome 
"""

}



