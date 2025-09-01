#!/bin/bash
# Array of run directories to submit
dirs=("matlab")
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
       --mail-user=coastesthydro@gmail.com \
       --mail-type=BEGIN,END,FAIL       \
        /storage/cms/grimesdj_lab/grimesdj/matlab_bg.sh main_process_ripchannel_test2D_NetCDF.m test2D.out
done
