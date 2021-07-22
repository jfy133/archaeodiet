//
// Check input samplesheet and get read channels
//

params.options = [:]

include { DAMAGEPROFILER } from '../../modules/nf-core/modules/damageprofiler/main' addParams( options: params.options )

workflow RUN_DAMAGEPROFILING {
    take:
    bam_and_reflist

    main:
    bam_and_reflist
        .multiMap {
            bam: [it[0], it[1]]
            reflist: [it0, it[2]]
        }
        .set { reads }

    emit:
    reads // channel: [ val(meta), [ reads ] ]
}
