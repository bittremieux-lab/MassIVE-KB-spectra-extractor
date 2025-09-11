#!/usr/bin/env nextflow

process DOWNLOAD_METADATA{
    tag "download_metadata"
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
    tag "group_tsv"

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
    tag "mgf_conversion:${mzml_group_tsv.baseName}"
    maxForks 5
    errorStrategy 'ignore'

    input:
        path mzml_group_tsv
        val cleaned_signal

    output:
        path "*.mgf", optional: true

    script:
    """
    # Set environment variable so Python script knows where to write failed logs
    export PIPELINE_DIR="${projectDir}"

    # Run the Python script - let it succeed or fail naturally
    # Failed logs are written to persistent failed_logs/ directory by the Python script
    python ${projectDir}/mzml_group_to_mgf.py $mzml_group_tsv
    """
}

process MERGE_MGFS {
    input:
        path mgf_files
        val task_id

    output:
        path "massiveKB_${task_id}.mgf", optional: true

    publishDir "results", mode: 'symlink'

    when:
    mgf_files.size() > 0

    script:
    """
    if [ \$(ls *.mgf 2>/dev/null | wc -l) -gt 0 ]; then
        find . -maxdepth 1 -name "*.mgf" ! -name "massiveKB_${task_id}.mgf" -print0 \
        | xargs -0 cat > massiveKB_${task_id}.mgf
        echo "Merged \$(ls *.mgf | wc -l) MGF files"
    else
        echo "No MGF files to merge"
        touch massiveKB_${task_id}.mgf
    fi
    """
}

process CLEAN_FAILED_LOGS_DIR {
    tag "clean_failed_logs"
    cache false

    input:
        path groups_dir

    output:
        val "CLEANED"

    script:
    """
    # Clean the failed_logs directory before starting MGF processing
    # This ensures that previously failed files that now succeed are properly removed
    if [ -d "${projectDir}/failed_logs" ]; then
        echo "Cleaning existing failed_logs directory..."
        rm -rf ${projectDir}/failed_logs/*.failed 2>/dev/null || true
        echo "Failed logs directory cleaned"
    else
        echo "No failed_logs directory to clean"
    fi
    """
}

process COLLECT_FAILED_LOGS {
    tag "collect_failed"
    cache false  // Never cache this process

    input:
        val mgf_processing_done  // Dependency to ensure this runs after MGF processing

    output:
        path "all_failed_logs.txt", optional: true

    publishDir "results", mode: 'copy'

    script:
    """
    echo "=== FAILED PROCESSES SUMMARY ===" > all_failed_logs.txt
    echo "Generated on: \$(date)" >> all_failed_logs.txt
    echo "" >> all_failed_logs.txt

    # Check if failed_logs directory exists and has files in the pipeline directory
    if [ -d "${projectDir}/failed_logs" ] && [ \$(ls ${projectDir}/failed_logs/*.failed 2>/dev/null | wc -l) -gt 0 ]; then
        failed_count=\$(ls ${projectDir}/failed_logs/*.failed | wc -l)
        echo "Total failed files: \$failed_count" >> all_failed_logs.txt
        echo "" >> all_failed_logs.txt

        for failed_file in ${projectDir}/failed_logs/*.failed; do
            if [ -f "\$failed_file" ]; then
                echo "=== FAILURE: \$(basename \$failed_file) ===" >> all_failed_logs.txt
                cat "\$failed_file" >> all_failed_logs.txt
                echo "" >> all_failed_logs.txt
                echo "---" >> all_failed_logs.txt
                echo "" >> all_failed_logs.txt
            fi
        done
    else
        echo "No failed files found" >> all_failed_logs.txt
        # Remove the file if no failures to avoid confusion
        rm all_failed_logs.txt
    fi
    """
}

process CREATE_PROCESSING_SUMMARY {
    cache false  // Never cache this process

    input:
        path mgf_files
        val total_inputs
        val mgf_processing_done  // Dependency to ensure this runs after MGF processing

    output:
        path "processing_summary.txt"

    publishDir "results", mode: 'copy'

    script:
    """
    echo "=== PROCESSING SUMMARY ===" > processing_summary.txt
    echo "Generated on: \$(date)" >> processing_summary.txt
    echo "" >> processing_summary.txt

    # Count successful MGF files
    successful_files=\$(ls *.mgf 2>/dev/null | wc -l)

    # Count failed files from persistent directory in pipeline directory
    if [ -d "${projectDir}/failed_logs" ]; then
        failed_files=\$(ls ${projectDir}/failed_logs/*.failed 2>/dev/null | wc -l)
    else
        failed_files=0
    fi

    # Total should be the number of input files
    total_files=$total_inputs

    echo "Total files processed: \$total_files" >> processing_summary.txt
    echo "Successful: \$successful_files" >> processing_summary.txt
    echo "Failed: \$failed_files" >> processing_summary.txt

    if [ \$total_files -gt 0 ]; then
        success_rate=\$(echo "scale=2; \$successful_files * 100 / \$total_files" | bc -l)
        echo "Success rate: \${success_rate}%" >> processing_summary.txt
    else
        echo "Success rate: 0%" >> processing_summary.txt
    fi

    echo "" >> processing_summary.txt

    if [ \$failed_files -gt 0 ]; then
        echo "=== FAILED FILES ===" >> processing_summary.txt
        for failed_file in ${projectDir}/failed_logs/*.failed; do
            if [ -f "\$failed_file" ]; then
                basename "\$failed_file" .failed >> processing_summary.txt
            fi
        done
    fi
    """
}


workflow extract_psms{
    // Check if task_id is provided
    if (!params.task_id) {
        error "Please provide a task_id parameter: --task_id YOUR_TASK_ID"
    }

    metadata_tsv = DOWNLOAD_METADATA(params.task_id)
    mzml_groups_dir = GROUP_TSV(metadata_tsv, "filename")

    mzml_groups = mzml_groups_dir.flatMap { dir -> file("${dir}/*") }

    // Count total input files for summary
    total_count = mzml_groups.count()

    // Clean failed_logs directory before starting MGF processing
    cleaned_signal = CLEAN_FAILED_LOGS_DIR(mzml_groups_dir)

    // Process all mzML groups - let them succeed or fail naturally
    // Failed logs are written to persistent failed_logs/ directory
    mgf_files = MZML_GROUP_TO_MGF(mzml_groups, cleaned_signal)

    // Collect successful MGF files and merge them
    mgf_files_collected = mgf_files.filter { it.size() > 0 }.collect()
    merged_mgf = MERGE_MGFS(mgf_files_collected, params.task_id)

    // Create a completion signal after all MGF processing is done
    // This ensures collection processes run after all MGF processes (successful and failed)
    mgf_processing_complete = mgf_files.collect().map { "MGF_PROCESSING_DONE" }

    // Collect failed logs - never cached, always reflects current state
    failed_summary = COLLECT_FAILED_LOGS(mgf_processing_complete)

    // Create processing summary - never cached, always reflects current state
    processing_summary = CREATE_PROCESSING_SUMMARY(
        mgf_files_collected,
        total_count,
        mgf_processing_complete
    )
}

// Default workflow
workflow {
    extract_psms()
}