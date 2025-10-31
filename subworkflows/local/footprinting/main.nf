// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { BEDTOOLS_INTERSECT                   } from '../../../modules/nf-core/bedtools/intersect/main'
include { GAWK as FILTERBYDEPTH                } from '../../../modules/nf-core/gawk/main'
include { GAWK as FILTERBYDEPTHDYN             } from '../../../modules/nf-core/gawk/main'
include { GAWK as GETPERCENTILEVAL             } from '../../../modules/nf-core/gawk/main'
include { TABIX_BGZIPTABIX                     } from '../../../modules/nf-core/tabix/bgziptabix/main'
include { GETOPENREGION                        } from '../../../modules/local/getopenregion/main'
include { GAWK as SLOPBOTH                     } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_MERGE                       } from '../../../modules/nf-core/bedtools/merge/main'
include { GAWK as DNASENUM                     } from '../../../modules/nf-core/gawk/main'
include { BEDTOOLS_INTERSECT as PEAKINTERSECT  } from '../../../modules/nf-core/bedtools/intersect/main'
include { GNU_CUT                              } from '../../../modules/local/gnu/cut/main'
include { GNU_SORT                             } from '../../../modules/nf-core/gnu/sort/main'
include { GETFPCOUNTS                          } from '../../../modules/local/getfpcounts/main'
include { GETUNINFOSITES                       } from '../../../modules/local/getuninfosites/main'
include { FOOTPRINTEXT                         } from '../../../modules/local/footprintext/main'
include { GNU_TAIL                             } from '../../../modules/local/gnu/tail/main'
include { BEDTOOLS_INTERSECT as FINALINTERSECT } from '../../../modules/nf-core/bedtools/intersect/main'
include { GAWK as FILTERP                      } from '../../../modules/nf-core/gawk/main'

workflow FOOTPRINTING {
    take:
    ch_sites    // channel: [ val(meta), [ site ] ]
    ch_peak     // channel: [ val(meta), [ peak ] ]
    depth       // integer
    scripts_dir
    percentile  // float

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    BEDTOOLS_INTERSECT(ch_sites.join(ch_peak), [[:], []])
    ch_versions = ch_versions.mix(BEDTOOLS_INTERSECT.out.versions.first())

    ch_bed = BEDTOOLS_INTERSECT.out.intersect.map { meta, bed -> [meta + ['depth': depth], bed] }

    if (percentile) {
        FILTERBYDEPTH(ch_bed, [], false)
        ch_filtered = FILTERBYDEPTH.out.output
        ch_percentile = Channel.value(percentile)
    }
    else {
        FILTERBYDEPTHDYN(ch_bed, [], false)
        ch_filtered = FILTERBYDEPTHDYN.out.output
        GETPERCENTILEVAL(ch_filtered, [], false)
        ch_percentile = GETPERCENTILEVAL.out.output
            .map { _meta, output_file ->
                def percentile_value = output_file.text.trim()
                return percentile_value
            }
            .first()
    }

    TABIX_BGZIPTABIX(ch_filtered)
    ch_versions = ch_versions.mix(TABIX_BGZIPTABIX.out.versions.first())

    ch_aligned_inputs = TABIX_BGZIPTABIX.out.gz_tbi.join(ch_peak.map { meta, bed -> [meta + ['depth': depth], bed] })
    GETOPENREGION(
        ch_aligned_inputs.map { meta, gz, tbi, _bed -> [meta, gz, tbi] },
        ch_aligned_inputs.map { meta, _gz, _tbi, bed -> [meta, bed] },
        ch_percentile,
    )

    ch_openreg = GETOPENREGION.out.bed

    SLOPBOTH(ch_openreg, [], false)

    BEDTOOLS_MERGE(SLOPBOTH.out.output)
    ch_versions = ch_versions.mix(BEDTOOLS_MERGE.out.versions.first())

    DNASENUM(BEDTOOLS_MERGE.out.bed, [], false)

    PEAKINTERSECT(ch_filtered.join(DNASENUM.out.output), [[:], []])
    ch_versions = ch_versions.mix(PEAKINTERSECT.out.versions.first())

    GNU_CUT(PEAKINTERSECT.out.intersect)
    ch_versions = ch_versions.mix(GNU_CUT.out.versions.first())

    GNU_SORT(GNU_CUT.out.cut)

    GETFPCOUNTS(GNU_SORT.out.sorted)
    ch_versions = ch_versions.mix(GETFPCOUNTS.out.versions.first())

    GETUNINFOSITES(GNU_SORT.out.sorted.join(DNASENUM.out.output))
    ch_versions = ch_versions.mix(GETUNINFOSITES.out.versions.first())

    ch_ftp_inputs = GETFPCOUNTS.out.bed.join(GETUNINFOSITES.out.bed).join(DNASENUM.out.output)
    FOOTPRINTEXT(ch_ftp_inputs, scripts_dir)
    ch_versions = ch_versions.mix(FOOTPRINTEXT.out.versions.first())

    GNU_TAIL(FOOTPRINTEXT.out.bed)

    FINALINTERSECT(GNU_TAIL.out.out_file.join(ch_openreg), [[:], []])
    ch_versions = ch_versions.mix(FINALINTERSECT.out.versions.first())

    FILTERP(FINALINTERSECT.out.intersect, [], false)
    ch_versions = ch_versions.mix(FILTERP.out.versions.first())

    emit:
    ftprts   = FILTERP.out.output // channel: [ val(meta), [ bed ] ]
    versions = ch_versions // channel: [ versions.yml ]
}
