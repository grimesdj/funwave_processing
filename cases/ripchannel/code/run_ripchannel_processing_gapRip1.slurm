#!/bin/bash
# Array of run directories to submit
dirs=("gapRip1")
# job_name="$dir"
# output_file="$dir.out"

module load matlab

for dir in "${dirs[@]}"; do
    job_name="$dir"
    output_file="$dir.out"
    sbatch --time=10400 \
       --job-name="$job_name" \
       --output="$output_file" \
       --mem=500G \
       --partition=general\
       --mail-user=coastesthydro@gmail.com \
       --mail-type=BEGIN,END,FAIL          \
        matlab_bg.sh process_ripchannel_gapRip1.m gapRip1.out
done
