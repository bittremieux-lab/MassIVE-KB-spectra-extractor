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
        val task_id

    output:
        path "*.mgf", optional: true

    script:
    """
    # Set environment variable so Python script knows where to write failed logs
    export PIPELINE_DIR="${projectDir}"
    export TASK_ID="${task_id}"

    # Run the Python script - let it succeed or fail naturally
    # Failed logs are written to persistent failed_logs_${task_id}/ directory by the Python script
    python ${projectDir}/mzml_group_to_mgf.py $mzml_group_tsv
    """
}

process MERGE_MGFS {
    input:
        path mgf_files
        val task_id

    output:
        path "massiveKB_${task_id}.mgf", optional: true

    publishDir "results_${task_id}", mode: 'symlink'

    when:
    mgf_files.size() > 0

    script:
    """
    # Use find to handle large numbers of MGF files (avoiding ls *.mgf which fails with too many files)
    mgf_count=\$(find . -maxdepth 1 -name "*.mgf" ! -name "massiveKB_${task_id}.mgf" | wc -l)

    if [ \$mgf_count -gt 0 ]; then
        # Use find with xargs to safely handle large numbers of files
        find . -maxdepth 1 -name "*.mgf" ! -name "massiveKB_${task_id}.mgf" -print0 | xargs -0 cat > massiveKB_${task_id}.mgf
        echo "Merged \$mgf_count MGF files into massiveKB_${task_id}.mgf"
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
        val task_id

    output:
        val "CLEANED"

    script:
    """
    # Clean the failed_logs directory before starting MGF processing
    # This ensures that previously failed files that now succeed are properly removed
    if [ -d "${projectDir}/failed_logs_${task_id}" ]; then
        echo "Cleaning existing failed_logs_${task_id} directory..."
        rm -rf ${projectDir}/failed_logs_${task_id}/*.csv 2>/dev/null || true
        echo "Failed logs directory cleaned"
    else
        echo "No failed_logs_${task_id} directory to clean"
    fi
    """
}

process COLLECT_FAILED_LOGS {
    tag "collect_failed"
    cache false  // Never cache this process

    input:
        val mgf_processing_done  // Dependency to ensure this runs after MGF processing
        val task_id

    output:
        path "failed_processes.csv", optional: true

    publishDir "results_${task_id}", mode: 'copy'

    script:
    """
    # Check if failed_logs directory exists and has files in the pipeline directory
    if [ -d "${projectDir}/failed_logs_${task_id}" ] && [ \$(ls ${projectDir}/failed_logs_${task_id}/*.csv 2>/dev/null | wc -l) -gt 0 ]; then
        # Create CSV header
        echo "mzml_file,error_message,spectra_count" > failed_processes.csv

        # Merge all CSV files
        for failed_file in ${projectDir}/failed_logs_${task_id}/*.csv; do
            if [ -f "\$failed_file" ]; then
                cat "\$failed_file" >> failed_processes.csv
            fi
        done

        failed_count=\$(ls ${projectDir}/failed_logs_${task_id}/*.csv | wc -l)
        echo "Merged \$failed_count failed process CSV files"
    else
        echo "No failed files found"
        # Don't create the file if no failures
    fi
    """
}

process CREATE_PROCESSING_SUMMARY {
    cache false  // Never cache this process

    input:
        path mgf_files
        val total_inputs
        val mgf_processing_done  // Dependency to ensure this runs after MGF processing
        val task_id

    output:
        path "processing_summary.txt"

    publishDir "results_${task_id}", mode: 'copy'

    script:
    """
    echo "=== PROCESSING SUMMARY ===" > processing_summary.txt
    echo "Generated on: \$(date)" >> processing_summary.txt
    echo "Task ID: ${task_id}" >> processing_summary.txt
    echo "" >> processing_summary.txt

    # Count successful MGF files
    successful_files=\$(find . -maxdepth 1 -name "*.mgf" | wc -l)

    # Count failed files from persistent directory in pipeline directory
    if [ -d "${projectDir}/failed_logs_${task_id}" ]; then
        failed_files=\$(ls ${projectDir}/failed_logs_${task_id}/*.csv 2>/dev/null | wc -l)
    else
        failed_files=0
    fi

    # Total should be the number of input files
    total_files=$total_inputs

    echo "Total files processed: \$total_files" >> processing_summary.txt
    echo "Successful: \$successful_files" >> processing_summary.txt
    echo "Failed: \$failed_files" >> processing_summary.txt

    if [ \$total_files -gt 0 ]; then
        # Calculate success rate using shell arithmetic (avoiding bc dependency)
        success_rate=\$(( successful_files * 100 / total_files ))
        echo "Success rate: \${success_rate}%" >> processing_summary.txt
    else
        echo "Success rate: 0%" >> processing_summary.txt
    fi

    echo "" >> processing_summary.txt

    if [ \$failed_files -gt 0 ]; then
        echo "=== FAILED FILES ===" >> processing_summary.txt
        for failed_file in ${projectDir}/failed_logs_${task_id}/*.csv; do
            if [ -f "\$failed_file" ]; then
                basename "\$failed_file" .csv >> processing_summary.txt
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
    cleaned_signal = CLEAN_FAILED_LOGS_DIR(mzml_groups_dir, params.task_id)

    // Process all mzML groups - let them succeed or fail naturally
    // Failed logs are written to persistent failed_logs_${task_id}/ directory
    mgf_files = MZML_GROUP_TO_MGF(mzml_groups, cleaned_signal, params.task_id)

    // Collect successful MGF files and merge them
    mgf_files_collected = mgf_files.filter { it.size() > 0 }.collect()
    merged_mgf = MERGE_MGFS(mgf_files_collected, params.task_id)

    // Create a completion signal after all MGF processing is done
    // This ensures collection processes run after all MGF processes (successful and failed)
    mgf_processing_complete = mgf_files.collect().map { "MGF_PROCESSING_DONE" }

    // Collect failed logs - never cached, always reflects current state
    failed_summary = COLLECT_FAILED_LOGS(mgf_processing_complete, params.task_id)

    // Create processing summary - never cached, always reflects current state
    processing_summary = CREATE_PROCESSING_SUMMARY(
        mgf_files_collected,
        total_count,
        mgf_processing_complete,
        params.task_id
    )
}

// Default workflow
workflow {
    extract_psms()
}