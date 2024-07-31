include { BWAMEM2_MEM }                                 from '../../../modules/nf-core/bwamem2/mem/main'
include { SENTIEON_BWAMEM }                             from '../../../modules/nf-core/bwamem2/mem/main'
include { PICARD_SORTSAM as PRE_SORTSAM }               from '../../../modules/nf-core/picard/sortsam/main'
include { PICARD_SORTSAM as SORT_CONSENSUS }            from '../../../modules/nf-core/picard/sortsam/main'

workflow AlignReads {
    take:
    reads                   // channel: [mandatory] [meta, reads]
    fasta                   // channel: [mandatory] [meta, fasta]
    index                   // path: bwamem2/*
    sort_bam                // bool: [mandatory] true -> sort, false -> don't sort
    aligner

    main:
    // Initialize output channels
    dupmarked_bams = Channel.empty()
    reports = Channel.empty()
    versions = Channel.empty()

    // Align
    reads_to_align = reads

    switch (aligner) {
        case "bwa2":
            BWAMEM2_MEM(reads,  index.map{ it -> [ [ id:'index' ], it ] }, sort_bam)
            aligned_reads = BWAMEM2_MEM.out.bam
            break
        case "sentieon-bwamem":
            // SENTIEON_BWAMEM()
            // aligned = SENTIEON_BWAMEM.out.bam
            break
        default:
            error "Unknown aligner: ${aligner}"
    }


    // MarkDuplicates - must be coordinate-sorted
    PICARD_MARKDUPLICATES(PRE_SORTSAM(aligned_reads, "coordinate").out.bam)

    // Gather outputs
    dupmarked_bams = dupmarked_bams.mix(PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, failOnDuplicate: true, failOnMismatch: true))
    reports = reports.mix(PICARD_MARKDUPLICATES.out.metrics)

    versions = versions.mix(BWAMEM2_MEM.out.versions)
    versions = versions.mix(PRE_SORTSAM.out.versions)
    versions = versions.mix(PICARD_MARKDUPLICATES.out.versions)

    emit:
    aligned_bams    =   aligned_reads       // channel: [ [meta], bam ]
    dupmarked_bams                          // channel: [ [meta], bam, bai ]
    reports
    versions                                // channel: [ versions.yml ]
}
