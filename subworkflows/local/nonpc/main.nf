// TODO nf-core: If in doubt look at other nf-core/subworkflows to see how we are doing things! :)
//               https://github.com/nf-core/modules/tree/master/subworkflows
//               You can also ask for help via your pull request or on the #subworkflows channel on the nf-core Slack workspace:
//               https://nf-co.re/join
// TODO nf-core: A subworkflow SHOULD import at least two modules

include { CENTERINGQNORM                 } from '../../../modules/local/centeringqnorm/main'
include { GAWK as FILTERSITESBYDEPTH     } from '../../../modules/nf-core/gawk/main'
include { GNU_SORT as SORTORIBG          } from '../../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as SORTQNORM          } from '../../../modules/nf-core/gnu/sort/main'
include { GNU_SORT as SORTQNORMFINAL     } from '../../../modules/nf-core/gnu/sort/main'
include { IGVTOOLS_TOTDF as ORICRTOTDF   } from '../../../modules/local/igvtools/totdf/main'
include { IGVTOOLS_TOTDF as QNORMCRTOTDF } from '../../../modules/local/igvtools/totdf/main'

workflow NONPC {
    take:
    ch_sites // channel: [ val(meta), [ sites ] ]
    species  // string

    main:

    ch_versions = Channel.empty()

    // TODO nf-core: substitute modules here for the modules of your subworkflow

    FILTERSITESBYDEPTH(ch_sites, [], false)

    CENTERINGQNORM(FILTERSITESBYDEPTH.out.output, species)
    ch_versions = ch_versions.mix(CENTERINGQNORM.out.versions.first())

    SORTORIBG(CENTERINGQNORM.out.cr_ori_bg)
    SORTQNORM(CENTERINGQNORM.out.cr_qnorm_bg)
    SORTQNORMFINAL(CENTERINGQNORM.out.cr_qnorm_final)

    ORICRTOTDF(SORTORIBG.out.sorted, species, 'bg')
    QNORMCRTOTDF(SORTQNORM.out.sorted, species, 'bg')

    emit:
    ridges_plot  = CENTERINGQNORM.out.plot_file // channel: [ val(meta), [ pdf ] ]
    dropout_plot = CENTERINGQNORM.out.dropout_plot // channel: [ val(meta), [ png ] ]
    ori_bg       = SORTORIBG.out.sorted // channel: [ val(meta, [ bedgraph ] ]
    qnorm_final  = SORTQNORMFINAL.out.sorted // channel: [ val(meta), [ bed ] ]
    ori_tdf      = ORICRTOTDF.out.tdf // channel: [ val(meta), [ tdf ] ]
    qnorm_tdf    = QNORMCRTOTDF.out.tdf // channel: [ val(meta), [ tdf ] ]
    versions     = ch_versions // channel: [ versions.yml ]
}
