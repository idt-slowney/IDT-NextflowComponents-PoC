include { FASTQC as PreprocessFASTQC }  from '../../../modules/nf-core/fastqc/main'
include { SEQTK_SAMPLE }                from '../../../modules/nf-core/seqtk/sample/main'
include { CUTADAPT_EXTRACTUMI }         from '/mnt/archive/work/slowney/IDT-Nextflow/modules/IDT/cutadapt/extractumi/main'
include { TRIMGALORE }                  from '../../../modules/nf-core/trimgalore/main'
include { FASTQC as PostprocessFASTQC } from '../../../modules/nf-core/fastqc/main'

workflow PROCESSFASTQ {

    take:
    reads    // channel: [ val(meta), [ reads ] ]
    skip_fastqc
    skip_downsample
    umi_flag
    skip_extractumi
    read_structure
    skip_trimming

    main:
    // Instantiate empty channels
    ch_versions = Channel.empty()
    fastqc_html = Channel.empty()
    fastqc_zip = Channel.empty()

    // Preprocessing QC -- skip_fastqc (DEFAULT: FALSE)
    if (!skip_fastqc) {
        PreprocessFASTQC(reads)
        ch_versions = ch_versions.mix(PreprocessFASTQC.out.versions)
        fastqc_html = fastqc_html.mix(PreprocessFASTQC.out.html)
        fastqc_zip = fastqc_zip.mix(PreprocessFASTQC.out.zip)
    }

    // Downsampling -- skip_downsample (DEFAULT: FALSE)
    // For now - we branch and mix here to account for situations with partial downsampling
    // #TODO - decide if we want to keep it like this, where downsampling is done per sample in the samplesheet
    // or set one downsampling level as a parameter for all samples - this will probably impact the conditional
    sample_reads = reads
    if (!skip_downsample) {
        downsample_reads = Channel.of(sample_reads).branch{
                    Sample: it[0]["sample_size"] && it[0]["sample_size"] != null
                    No_Sample: it[0]["sample_size"] == null | !("sample_size" in it[0]) }
                    
        SEQTK_SAMPLE(
            downsample_reads.Sample
            .map{ meta, reads -> [meta, reads, meta.sample_size] } )
        ch_versions = ch_versions.mix(SEQTK_SAMPLE.out.versions)
        sample_reads = downsample_reads.No_Sample.mix(SEQTK_SAMPLE.out.reads)
    }

    // ExtractUMI -- umi_flag (FALSE), skip_extractumi (FALSE)
    extractumi_reads = sample_reads
    if (umi_flag && !skip_extractumi) {
        CUTADAPT_EXTRACTUMI(extractumi_reads, read_structure)
        ch_versions = ch_versions.mix(CUTADAPT_EXTRACTUMI.out.versions)
        extractumi_reads = CUTADAPT_EXTRACTUMI.out.reads
    }

    // Trim Adapters -- skip_trimming (FALSE)
    trim_reads = extractumi_reads
    trimgalore_log          = Channel.empty()
    trimgalore_unpaired     = Channel.empty()
    trimgalore_html         = Channel.empty()
    trimgalore_zip          = Channel.empty()
    if (!skip_trimming) {
        TRIMGALORE(trim_reads)
        // Here we can include other trimming tools and control which get run with the when directive
        // i.e. TRIMMOMATIC(trim_reads) 
        // ---  { ext.when = trimmer = trimmomatic }
        // Then we mix results of all tools with a previously instantiated empty channel and yield only the output from the process that was run
        // i.e. trimmed_out = channel.empty()
        // ---  trimmed_out = trimmed_out.mix(TRIMGALORE.out.reads)
        // ---  trimmed_out = trimmed_out.mix(TRIMMOMATIC.out.reads)
        trimgalore_log          = trimgalore_log.mix(TRIMGALORE.out.log)
        trimgalore_unpaired     = trimgalore_unpaired.mix(TRIMGALORE.out.unpaired)
        trimgalore_html         = trimgalore_html.mix(TRIMGALORE.out.html)
        trimgalore_zip          = trimgalore_zip.mix(TRIMGALORE.out.zip)
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions)
        trim_reads = TRIMGALORE.out.reads
    } 
    // Postprocessing QC
    PostprocessFASTQC(trim_reads)
    ch_versions = ch_versions.mix(PostprocessFASTQC.out.versions)
    fastqc_html = fastqc_html.mix(PostprocessFASTQC.out.html)
    fastqc_zip = fastqc_zip.mix(PostprocessFASTQC.out.zip)

    emit:
    // TODO nf-core: edit emitted channels
    reads               = trim_reads           // channel: [ val(meta), [ reads ] ]
    fastqc_html         = fastqc_html          // channel: [ val(meta), path(html) ]
    fastqc_zip          = fastqc_zip           // channel: [ val(meta), path(zip) ]
    trimgalore_log                             // channel: [ val(meta), path(log) ]
    trimgalore_unpaired                        // channel: [ val(meta), [ unpaired_reads ] ]
    trimgalore_html                            // channel: [ val(meta), path(html) ]
    trimgalore_zip                             // channel: [ val(meta), path(zip) ]
    versions = ch_versions                     // channel: [ versions.yml ]
}

