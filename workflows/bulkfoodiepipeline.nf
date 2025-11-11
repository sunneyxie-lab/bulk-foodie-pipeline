/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { FASTQC                                 } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE                             } from '../modules/nf-core/trimgalore/main'
include { BISMARK_ALIGN                          } from '../modules/nf-core/bismark/align/main'
include { SAMTOOLS_SORT                          } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_INDEX                         } from '../modules/nf-core/samtools/index/main'
include { BISMARK_DEDUPLICATE                    } from '../modules/nf-core/bismark/deduplicate/main'
include { METHYLEXTRACT                          } from '../modules/local/methylextract/main'
include { CALLTRACK                              } from '../modules/local/calltrack/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORTBYNAME   } from '../modules/nf-core/samtools/sort/main'
include { BEDTOOLS_BAMTOBED                      } from '../modules/nf-core/bedtools/bamtobed/main'
include { POSTPROCESSBED                         } from '../modules/local/postprocessbed/main'
include { BEDTOOLS_GENOMECOV                     } from '../modules/nf-core/bedtools/genomecov/main'
include { UCSC_BEDGRAPHTOBIGWIG                  } from '../modules/nf-core/ucsc/bedgraphtobigwig/main'
include { MACS3_CALLPEAK as MACS_CALLPEAK        } from '../modules/nf-core/macs3/callpeak/main'
include { CALLRATIOANDDEPTH                      } from '../modules/local/callratioanddepth/main'
include { IGVTOOLS_TOTDF as IGVTOOLS_TOTDF_RATIO } from '../modules/local/igvtools/totdf/main'
include { IGVTOOLS_TOTDF as IGVTOOLS_TOTDF_DEPTH } from '../modules/local/igvtools/totdf/main'
include { HIGHSCOREPEAKS                         } from '../modules/local/highscorepeaks/main'
include { PLOT                                   } from '../modules/local/plot/main'
include { BEDTOOLS_INTERSECT as TSS_INTERSECT    } from '../modules/nf-core/bedtools/intersect/main'
include { PLOTTSS                                } from '../modules/local/plottss/main'
include { MULTIQC                                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMultiqc                   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML                 } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText                 } from '../subworkflows/local/utils_nfcore_bulkfoodiepipeline_pipeline'
include { FOOTPRINTING                           } from '../subworkflows/local/footprinting'
include { BEDTOOLS_SPLIT as PEAKS_SPLIT          } from '../modules/nf-core/bedtools/split/main'
include { CAT_CAT as CONCAT_FOOTPRINTS           } from '../modules/nf-core/cat/cat/main'
include { BEDTOOLS_SORT as SORT_FOOTPRINTS       } from '../modules/nf-core/bedtools/sort/main'
include { ADJUSTRATIO                            } from '../modules/local/adjustratio/main'
include { GAWK as REFORMATRATIOADJUST            } from '../modules/nf-core/gawk/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow BULKFOODIEPIPELINE {
    take:
    ch_samplesheet // channel: samplesheet read in from --input
    ch_fasta // channel: fasta file
    ch_sizes // channel: sizes file
    ch_bismark_index // channel: [ path(bismark index)   ]
    genome_id
    macs_gsize
    tss
    depth
    scripts_dir
    expected_ratio_file

    main:

    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()

    // SEQKIT_SPLIT2(ch_samplesheet)>
    // ch_versions = ch_versions.mix(SEQKIT_SPLIT2.out.versions)

    // ch_fastq_pieces = SEQKIT_SPLIT2.out.reads.flatMap { meta, chunk_list ->
    //     // Group chunks by their part number
    //     def grouped_chunks = [:]
    //     chunk_list.each { chunk ->
    //         // Extract part number from filename (handle both fastq.gz and fq.gz)
    //         def part_match = chunk.name =~ /.*_(\d+)\.(fastq|fq)\.gz$/
    //         if (part_match) {
    //             def part_num = part_match[0][1]
    //             if (!grouped_chunks.containsKey(part_num)) {
    //                 grouped_chunks[part_num] = []
    //             }
    //             grouped_chunks[part_num] << chunk
    //         }
    //     }

    //     // Return paired chunks as tuples
    //     return grouped_chunks.collect { _part_num, chunks ->
    //         tuple(meta, chunks)
    //     }
    // }

    FASTQC(ch_samplesheet)
    ch_versions = ch_versions.mix(FASTQC.out.versions)

    TRIMGALORE(ch_samplesheet)
    ch_versions = ch_versions.mix(TRIMGALORE.out.versions)

    BISMARK_ALIGN(
        TRIMGALORE.out.reads,
        [[:], file(ch_fasta, checkIfExists: true)],
        [[:], file(ch_bismark_index, checkIfExists: true)],
    )
    ch_alignments = BISMARK_ALIGN.out.bam
    ch_alignment_reports = BISMARK_ALIGN.out.report.map { meta, report -> [meta, report, []] }
    ch_versions = ch_versions.mix(BISMARK_ALIGN.out.versions)

    /*
    * Run deduplicate_bismark
    */
    BISMARK_DEDUPLICATE(
        BISMARK_ALIGN.out.bam
    )
    ch_alignments = BISMARK_DEDUPLICATE.out.bam
    ch_alignment_reports = BISMARK_ALIGN.out.report.join(BISMARK_DEDUPLICATE.out.report)
    ch_versions = ch_versions.mix(BISMARK_DEDUPLICATE.out.versions)

    /*
     * MODULE: Run samtools sort on aligned or deduplicated bam
     */
    SAMTOOLS_SORT(
        ch_alignments,
        [[:], []],
    )
    ch_versions = ch_versions.mix(SAMTOOLS_SORT.out.versions)

    /*
     * MODULE: Run samtools index on aligned or deduplicated bam
     */
    SAMTOOLS_INDEX(
        SAMTOOLS_SORT.out.bam
    )
    ch_versions = ch_versions.mix(SAMTOOLS_INDEX.out.versions)

    ch_bam_sorted_bai = SAMTOOLS_SORT.out.bam.join(SAMTOOLS_INDEX.out.bai)

    METHYLEXTRACT(
        ch_bam_sorted_bai
    )
    ch_versions = ch_versions.mix(METHYLEXTRACT.out.versions)

    CALLTRACK(
        METHYLEXTRACT.out.bed
    )
    ch_versions = ch_versions.mix(CALLTRACK.out.versions)

    SAMTOOLS_SORTBYNAME(
        ch_alignments,
        [[:], []],
    )

    BEDTOOLS_BAMTOBED(
        SAMTOOLS_SORTBYNAME.out.bam
    )
    ch_versions = ch_versions.mix(BEDTOOLS_BAMTOBED.out.versions)

    POSTPROCESSBED(
        BEDTOOLS_BAMTOBED.out.bed
    )

    BEDTOOLS_GENOMECOV(
        POSTPROCESSBED.out.frag_bed.combine(Channel.from(1)),
        ch_sizes,
        'bedGraph',
        false,
    )
    ch_versions = ch_versions.mix(BEDTOOLS_GENOMECOV.out.versions)

    UCSC_BEDGRAPHTOBIGWIG(
        BEDTOOLS_GENOMECOV.out.genomecov,
        ch_sizes,
    )
    ch_versions = ch_versions.mix(UCSC_BEDGRAPHTOBIGWIG.out.versions)

    ch_bedpe = POSTPROCESSBED.out.frag_bed.map { meta, bedpe ->
        [meta, bedpe, []]
    }
    MACS_CALLPEAK(
        ch_bedpe,
        macs_gsize,
    )
    ch_versions = ch_versions.mix(MACS_CALLPEAK.out.versions)

    HIGHSCOREPEAKS(
        MACS_CALLPEAK.out.peak
    )

    PLOT(
        METHYLEXTRACT.out.txt
    )

    CALLRATIOANDDEPTH(
        CALLTRACK.out.summary_score_bed
    )

    IGVTOOLS_TOTDF_DEPTH(
        CALLRATIOANDDEPTH.out.depth,
        genome_id,
        'igv',
    )

    IGVTOOLS_TOTDF_RATIO(
        CALLRATIOANDDEPTH.out.ratio,
        genome_id,
        'igv',
    )

    TSS_INTERSECT(
        CALLRATIOANDDEPTH.out.sites.map { meta, sites ->
            [meta, sites, file(tss, checkIfExists: true)]
        },
        [[:], ch_sizes],
    )

    PLOTTSS(
        TSS_INTERSECT.out.intersect
    )

    ch_peak = HIGHSCOREPEAKS.out.bed

    PEAKS_SPLIT(ch_peak.map { meta, bed -> [meta, bed, 10] })

    ch_adjusted_sites = CALLRATIOANDDEPTH.out.sites
    if (expected_ratio_file) {
        ADJUSTRATIO(
            CALLRATIOANDDEPTH.out.sites,
            expected_ratio_file,
            ch_fasta,
        )
        REFORMATRATIOADJUST(
            ADJUSTRATIO.out.txt,
            [],
            false,
        )
        ch_adjusted_sites = REFORMATRATIOADJUST.out.output
    }

    ch_footprinting_input = ch_adjusted_sites
        .combine(PEAKS_SPLIT.out.beds.transpose(), by: 0)
        .map { meta, sites_file, bed_file ->
            def chunk_num = bed_file.baseName.tokenize('.')[-1] as Integer
            def new_meta = meta + [chunk: chunk_num]
            [new_meta, sites_file, bed_file]
        }

    FOOTPRINTING(
        ch_footprinting_input.map { meta, sites, _bed -> [meta, sites] },
        ch_footprinting_input.map { meta, _sites, bed -> [meta, bed] },
        depth,
        scripts_dir,
    )

    CONCAT_FOOTPRINTS(
        FOOTPRINTING.out.ftprts.map { meta, files ->
            [meta.findAll { !it.key.equals('chunk') }, files]
        }.groupTuple()
    )

    SORT_FOOTPRINTS(CONCAT_FOOTPRINTS.out.file_out, [])

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'bulkfoodiepipeline_software_' + 'mqc_' + 'versions.yml',
            sort: true,
            newLine: true,
        )
        .set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config = Channel.fromPath(
        "${projectDir}/assets/multiqc_config.yml",
        checkIfExists: true
    )
    ch_multiqc_custom_config = params.multiqc_config
        ? Channel.fromPath(params.multiqc_config, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_logo = params.multiqc_logo
        ? Channel.fromPath(params.multiqc_logo, checkIfExists: true)
        : Channel.empty()
    ch_multiqc_custom_methods_description = params.multiqc_methods_description
        ? file(params.multiqc_methods_description, checkIfExists: true)
        : file("${projectDir}/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description)
    )

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true,
        )
    )

    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect { it[1] }.ifEmpty([]))

    MULTIQC(
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        [],
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions // channel: [ path(versions.yml) ]
}
