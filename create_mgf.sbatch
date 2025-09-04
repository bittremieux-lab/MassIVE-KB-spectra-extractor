#!/bin/bash
#SBATCH --time=168:00:00
#SBATCH --partition=seven_days
#SBATCH --mem=50G
#SBATCH --cpus-per-task=10

source activate nf-env

nextflow main.nf -entry extract_psms --task_id 3cac03860ff7453a821332ab4cff20f4 -resume -with-trace trace-$(date +%Y%m%d-%H%M%S).tsv -with-report report-$(date +%Y%m%d-%H%M%S).html