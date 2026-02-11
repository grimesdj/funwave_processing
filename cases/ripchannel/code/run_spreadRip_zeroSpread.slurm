#!/bin/bash
# Array of run directories to submit
dirs=("spreadRip0")
# job_name="$dir"
# output_file="$dir.out"

module load matlab

for dir in "${dirs[@]}"; do
    job_name="$dir"
    output_file="$dir.out"
    sbatch --time=10080 \
       --job-name="$job_name" \
       --output="$output_file" \
       --mem=300G \
       --partition=general\
       --mail-user=coastesthydro@gmail.com \
       --mail-type=BEGIN,END,FAIL          \
        matlab_bg.sh main_process_spreadRip_zeroSpread.m spreadRip_zeroSpread.out
done
