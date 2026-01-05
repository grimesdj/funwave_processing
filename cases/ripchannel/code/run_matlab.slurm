#!/bin/bash
# Array of run directories to submit
dirs=("matlab2")
# job_name="$dir"
# output_file="$dir.out"

module load matlab

for dir in "${dirs[@]}"; do
    job_name="$dir"
    output_file="$dir.out"
    sbatch --time=14400 \
       --job-name="$job_name" \
       --output="$output_file" \
       --mem=300G \
       --partition=highmem\
       --mail-user=coastesthydro@gmail.com \
       --mail-type=BEGIN,END,FAIL       \
        matlab_bg.sh main_compare_statistics_across_runs.m matlab2.out
done
