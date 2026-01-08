#!/bin/bash
# Array of run directories to submit
dirs=("barRip0_s00")
# job_name="$dir"
# output_file="$dir.out"

module load matlab

for dir in "${dirs[@]}"; do
    job_name="$dir"
    output_file="$dir.out"
    sbatch --time=5760 \
       --job-name="$job_name" \
       --output="$output_file" \
       --mem=300G \
       --partition=highmem\
       --mail-user=coastesthydro@gmail.com \
       --mail-type=BEGIN,END,FAIL          \
        matlab_bg.sh main_process_ripchannel_WaveAvg_output.m barRip0_s00.out
done
