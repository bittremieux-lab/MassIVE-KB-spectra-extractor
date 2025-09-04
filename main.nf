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
        path "group_*.tsv"

    script:
    """
    group_tsv.py $psms_tsv $group_column
    """
}

process MZML_GROUP_TO_MGF{
    input:
        path mzml_group_tsv

    output:
        stdout

    script:
    """
    process_mzml_group.py $mzml_group_tsv
    """
}

// process COLLECT_MZTAB_MSGF {
//     scratch true
//
//     input:
//         val task_info
//
//     output:
//        path "*.tsv"
//
//     script:
//     """
//     curl -X POST "https://proteomics2.ucsd.edu/ProteoSAFe/DownloadResult?view=group_by_spectrum&task=${task_info['search_task_id']}&file=u.batch%252F${task_info['search_task_id']}%252FmzTab%252FTop1_results.mzTab" -o top1.zip
//     unzip top1.zip -x params.xml
//     """
// }

// process COLLECT_MZTAB_SYNTHETIC {
//     scratch true
//
//     input:
//         val task_info
//
//     output:
//         path "MSGF-PLUS-SYNTHETIC-*-group_by_spectrum-main.tsv"
//
//     script:
//     """
//     curl -X POST "https://proteomics2.ucsd.edu/ProteoSAFe/DownloadResult?view=view_result_list&task=${task_info['search_task_id']}" -o search_details.zip
//     unzip search_details.zip -x params.xml
//
//     tsv=\$(ls MSGF-PLUS-SYNTHETIC-*-view_result_list-main.tsv)
//     # Extract File_descriptor column
//     file_descriptor=\$(awk -F'\\t' 'NR==1 {for (i=1;i<=NF;i++) if (\$i=="File_descriptor") col=i} NR==2 {print \$col}' "\$tsv")
//     curl -X POST "https://proteomics2.ucsd.edu/ProteoSAFe/DownloadResult?view=group_by_spectrum&task=${task_info['search_task_id']}&file=\$file_descriptor" -o search_results.zip
//     unzip search_results.zip -x params.xml
//     """
// }

// process PROCESS_MZTAB{
//     input:
//         path tsv_file
//
//     output:
//         stdout
//
//     script:
//     """
//     process_mztab.py $tsv_file
//     """
// }

workflow extract_psms{
    metadata_tsv=DOWNLOAD_METADATA(params.task_id)
    mzml_groups=GROUP_TSV(metadata_tsv, "filename")
    MZML_GROUP_TO_MGF(mzml_groups)
//     search_tasks = tasks_tsv.splitCsv(header:true, sep:'\t').take(10)

    // Branch based on 'search_workflow' value
//     msgf_tasks = search_tasks.filter { it['search_workflow'] == 'MSGF-PLUS-AMBIGUITY' || it['search_workflow'] == 'MULTIPASS_MSGF_PLUS_DB_SEARCH' }
//     synthetic_tasks = search_tasks.filter { it['search_workflow'] == 'MSGF-PLUS-SYNTHETIC' }

    // Process each branch
//     msgf_tsvs = COLLECT_MZTAB_MSGF(msgf_tasks)
//     synthetic_tsvs = COLLECT_MZTAB_SYNTHETIC(synthetic_tasks)
//     all_tsvs = msgf_tsvs.mix(synthetic_tsvs)
//     PROCESS_MZTAB(all_tsvs)
}