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
    input:
        path psms_tsv
        val group_column

    output:
        path "groups"

    script:
    """
    group_tsv.py $psms_tsv $group_column "groups"
    """
}

process MZML_GROUP_TO_MGF{
    maxForks 5

    input:
        path mzml_group_tsv

    output:
        path "*.mgf"
        path "*.failed" optional true

    errorStrategy 'ignore'

    script:
    """
    python ${projectDir}/mzml_group_to_mgf.py $mzml_group_tsv
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
    find . -maxdepth 1 -name "*.mgf" -print0 | xargs -0 cat > massiveKB_${task_id}.mgf
    """
}

process MERGE_FAILED_LOGS {
    tag "merge_failed"

    input:
        path failed_logs

    output:
        path "all_failed_logs.txt"

    publishDir "results", mode: 'copy'  // or 'move' if you want to save space

    script:
    """
    cat $failed_logs > all_failed_logs.txt
    """
}


workflow extract_psms{
    metadata_tsv = DOWNLOAD_METADATA(params.task_id)
    mzml_groups_dir = GROUP_TSV(metadata_tsv, "filename")

    mzml_groups = mzml_groups_dir.flatMap { dir -> file("${dir}/*") }

    mgf_and_fails = MZML_GROUP_TO_MGF(mzml_groups)
    mgf_files = mgf_and_fails[0]
    failed_logs = mgf_and_fails[1]

    merged_mgf = MERGE_MGFS(mgf_files.collect(), params.task_id)
    failed_summary = MERGE_FAILED_LOGS(failed_logs.collect())
    failed_summary.view()
}