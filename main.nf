#!/usr/bin/env nextflow

process DOWNLOAD_METADATA{
    scratch true
    input:
        val task_id

    output:
        path "*.tsv"

    script:
    """
    # Download and unzip the metadata.
    curl -X POST "https://massive.ucsd.edu/ProteoSAFe/DownloadResult?view=candidate_library_spectra&task=${task_id}" -o candidate_library_spectra_${task_id}.zip
    unzip candidate_library_spectra_${task_id}.zip -x params.xml
    """
}

process GROUP_TSV{
    scratch true

    input:
        path psms_tsv
        val group_column

    output:
        path "group_*.tsv"

    script:
    """
    group_tsv.py $psms_tsv $group_column
    """
}

process MZML_GROUP_TO_MGF{
    maxForks 5
    scratch true

    input:
        path mzml_group_tsv

    output:
        path "*.mgf"

    script:
    """
    mzml_group_to_mgf.py $mzml_group_tsv
    """
}

process MERGE_MGFS {
    scratch true

    input:
        path mgf_files
        val task_id

    output:
        path "massiveKB_${task_id}.mgf"

    publishDir "results", mode: 'move'

    script:
    """
    cat $mgf_files > "massiveKB_${task_id}.mgf"
    """
}

workflow extract_psms{
    metadata_tsv = DOWNLOAD_METADATA(params.task_id)
    mzml_groups = GROUP_TSV(metadata_tsv, "filename").flatten()
    mgf_files = MZML_GROUP_TO_MGF(mzml_groups)
    merged_mgf = MERGE_MGFS(mgf_files.collect(), params.task_id)
}