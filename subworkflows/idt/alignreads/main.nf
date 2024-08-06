include { BWAMEM2_MEM }                                 from '../../../modules/idt/bwamem2/mem/main'
include { SENTIEON_BWAMEM }                             from '../../../modules/idt/sentieon/bwamem/main'
include { PICARD_SORTSAM as PRE_SORTSAM }               from '../../../modules/idt/picard/sortsam/main'
include { PICARD_MARKDUPLICATES }                       from '../../../modules/idt/picard/markduplicates/main'

workflow ALIGNREADS {
    take:
    reads                   // channel: [mandatory] [meta, reads]
    fasta                   // channel: [mandatory] [meta, fasta]
    fasta_fai
    index                   // path: bwamem2/*
    sort_bam                // bool: [mandatory] true -> sort, false -> don't sort
    aligner                 // str: "bwa2" or "sentieon-bwamem"

    main:
    // Initialize output channels
    dupmarked_bams = Channel.empty()
    reports = Channel.empty()
    versions = Channel.empty()

    // Align
    reads_to_align = reads

    switch (aligner) {
        case "bwa2":
            BWAMEM2_MEM(reads,  index, fasta, sort_bam)
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
    PICARD_MARKDUPLICATES(aligned_reads, fasta, fasta_fai)

    // Gather outputs
    dupmarked_bams = dupmarked_bams.mix(PICARD_MARKDUPLICATES.out.bam.join(PICARD_MARKDUPLICATES.out.bai, failOnDuplicate: true, failOnMismatch: true))
    reports = reports.mix(PICARD_MARKDUPLICATES.out.metrics)

    versions = versions.mix(BWAMEM2_MEM.out.versions)
    // versions = versions.mix(PRE_SORTSAM.out.versions)
    versions = versions.mix(PICARD_MARKDUPLICATES.out.versions)

    emit:
    aligned_bams    =   aligned_reads       // channel: [ [meta], bam ]
    dupmarked_bams                          // channel: [ [meta], bam, bai ]
    reports
    versions                                // channel: [ versions.yml ]
}
