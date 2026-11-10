#!/bin/bash
# Array of run directories to submit
dirs=("compStats")
# job_name="$dir"
# output_file="$dir.out"

module load matlab

for dir in "${dirs[@]}"; do
    job_name="$dir"
    output_file="$dir.out"
    sbatch --time=1000 \
       --job-name="$job_name" \
       --output="$output_file" \
       --mem=100G \
       --partition=general\
       --mail-user=coastesthydro@gmail.com \
       --mail-type=BEGIN,END,FAIL          \
        matlab_bg.sh compile_WaveAvgVelocity_stats_cmshpc.m compStats.out
done
