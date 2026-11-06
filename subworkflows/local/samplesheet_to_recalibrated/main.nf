include { resolveRelativeToProjectDir } from '../samplesheet_to_samples/main'

workflow SAMPLESHEET_TO_RECALIBRATED {
    take:
    samplesheet // path to samplesheet CSV (sample,alignment,alignment_index)

    main:
    rows = Channel.fromPath(samplesheet, checkIfExists: true)
        .splitCsv(header: true)
        .map { row -> validateRecalibratedRow(row) }

    unique_ids_ch = rows
        .map { meta, _alignment, _index -> meta.sample }
        .toList()
        .map { ids ->
            def dupes = ids.countBy { it }.findAll { _id, n -> n > 1 }.keySet().sort()
            if (dupes) {
                error "Duplicate sample ID(s) in samplesheet: ${dupes.join(', ')}"
            }
            [ ok: true ]
        }

    samples = rows
        .combine(unique_ids_ch)
        .map { meta, alignment, index, _ok -> [ meta, alignment, index ] }
        .ifEmpty { error "Samplesheet produced zero samples -- check that it has at least one data row below the header" }

    emit:
    samples // [ meta(sample), alignment, alignment_index ]
}

def validateRecalibratedRow(row) {
    if (!row.sample) error "Samplesheet row missing 'sample': ${row}"
    ['alignment', 'alignment_index'].each { col ->
        if (!row[col]) error "Samplesheet row '${row.sample}' missing '${col}'"
    }

    def meta      = [ sample: row.sample ]
    def alignment = file(resolveRelativeToProjectDir(row.alignment), checkIfExists: true)
    def index     = file(resolveRelativeToProjectDir(row.alignment_index), checkIfExists: true)
    [ meta, alignment, index ]
}
