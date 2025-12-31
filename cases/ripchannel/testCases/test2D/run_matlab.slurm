#!/bin/bash
# Array of run directories to submit
dirs=("matlab0")
# job_name="$dir"
# output_file="$dir.out"
#        --partition=highmem \

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
        /storage/cms/grimesdj_lab/grimesdj/matlab_bg.sh plot_all_test2Davg_1D_momentum.m matlab0.out &
done
