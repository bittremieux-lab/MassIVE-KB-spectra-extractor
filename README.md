# MassIVE-KB-spectra-extractor

A Nextflow pipeline for extracting spectra from MassIVE-KB with robust error handling and resume functionality.

## Quick Start

### Prerequisites

Nextflow and required python packages see `environment.yml`

### Installation

```bash
# Clone the repository
git clone https://github.com/bittremieux-lab/MassIVE-KB-spectra-extractor.git
cd MassIVE-KB-spectra-extractor

# Install dependencies (using conda/mamba)
conda env create -f environment.yml
conda activate nf-mkb
```

### Usage

```bash
# Basic usage
nextflow run main.nf --task_id YOUR_MASSIVE_TASK_ID

# Resume from cache
nextflow run main.nf --task_id YOUR_MASSIVE_TASK_ID -resume
```

## Pipeline Overview

The pipeline consists of several key processes:

1. **DOWNLOAD_METADATA**: Downloads a MassIVE-KB metadata file from the task ID
2. **GROUP_TSV**: Groups metadata by mzML/mzXML filename for parallel processing
3. **MZML_GROUP_TO_MGF**: Downloads mzML/mzXML file and converts to MGF format
4. **MERGE_MGFS**: Combines all successful MGF files into a single MGF file
5. **COLLECT_FAILED_LOGS**: Aggregates failure information for analysis
6. **CREATE_PROCESSING_SUMMARY**: Generates comprehensive processing report

## Error Handling & Resume Functionality
Due to inconsistencies in mzML/mzXML files and MassIVE-KB **MZML_GROUP_TO_MGF** might fail unexpectedly. 
This pipeline is designed so all successful **MZML_GROUP_TO_MGF** processes are reused from cache when 
rerunning with the same task_id and the -resume flag. At the same time, an overview of failed processes is generated
in `results_YOUR_TASK_ID/failed_processes.csv`. These issues can be solved by editing `mzml_group_to_mgf.py`. 
Edits to any other files might invalidate the cached versions of successful **MZML_GROUP_TO_MGF** processes.


### Workflow for Large Processing Jobs

1. **Initial run**: Process all files
   ```bash
   nextflow run main.nf --task_id YOUR_TASK_ID
   ```

2. **Check failures**: Review processing summary and failed logs
   ```bash
   # Check overall results (replace YOUR_TASK_ID with actual task ID)
   cat results_YOUR_TASK_ID/processing_summary.txt
   
   # Review specific failures
   cat results_YOUR_TASK_ID/failed_processes.csv
   ```

3. **Fix issues**: Update the Python script to handle specific error cases

4. **Resume processing**: Only failed processes will retry
   ```bash
   nextflow run main.nf --task_id YOUR_TASK_ID -resume
   ```

5. **Repeat**: Continue until all files process successfully

## License

This project is licensed under an Apache 2.0 license - see the LICENSE file for details.