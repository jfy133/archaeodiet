
//
// Run SAMtools stats, flagstat and idxstats - from https://github.com/nf-core/viralrecon/blob/master/subworkflows/nf-core/bam_stats_samtools.nf
//

params.options = [:]

include { SAMTOOLS_STATS    } from '../../modules/nf-core/modules/samtools/stats/main'    addParams( options: params.options )
include { SAMTOOLS_IDXSTATS } from '../../modules/nf-core/modules/samtools/idxstats/main' addParams( options: params.options )
include { SAMTOOLS_FLAGSTAT } from '../../modules/nf-core/modules/samtools/flagstat/main' addParams( options: params.options )

workflow BAM_STATS_SAMTOOLS {
    take:
    bam_bai // channel: [ val(meta), [ bam ], [bai] ]

    main:
    SAMTOOLS_STATS    ( bam_bai )
    SAMTOOLS_FLAGSTAT ( bam_bai )
    SAMTOOLS_IDXSTATS ( bam_bai )

    emit:
    stats    = SAMTOOLS_STATS.out.stats       // channel: [ val(meta), [ stats ] ]
    flagstat = SAMTOOLS_FLAGSTAT.out.flagstat // channel: [ val(meta), [ flagstat ] ]
    idxstats = SAMTOOLS_IDXSTATS.out.idxstats // channel: [ val(meta), [ idxstats ] ]
    version  = SAMTOOLS_STATS.out.version     //    path: *.version.txt
}
