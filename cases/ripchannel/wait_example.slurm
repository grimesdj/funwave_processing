## RUN Specific Header Material ##
## #!/bin/bash
## # Array of run directories
## run_dirs=(
##     "terRip1_h10t10s00d00"
##     "terRip1_h10t10s02d00"
##     "terRip1_h10t10s04d00"
##     "terRip1_h10t10s10d00"
##     "terRip1_h10t10s20d00"
## )

ulimit -s unlimited

prev_jobid=""

for dir in "${run_dirs[@]}"; do

    job_name="${dir##*_}"
    output_file="$dir.out"

    if [ -z "$prev_jobid" ]; then
        submit_output=$(sbatch \
            --time=11520 \
            --job-name="$job_name" \
            --output="$output_file" \
            --chdir="$dir" \
            --ntasks=128 \
            --cpus-per-task=1 \
            --mail-user=coastesthydro@gmail.com \
            --mail-type=BEGIN,END,FAIL \
            mpirun -n 128 ../../../software/funwave/funwave--mpiifort-parallel-double_vWaveAvg_vWMnew_v1 input.txt)
    else
        submit_output=$(sbatch \
            --dependency=afterany:$prev_jobid \
            --time=11520 \
            --job-name="$job_name" \
            --output="$output_file" \
            --chdir="$dir" \
            --ntasks=128 \
            --cpus-per-task=1 \
            --mail-user=coastesthydro@gmail.com \
            --mail-type=BEGIN,END,FAIL \
            mpirun -n 128 ../../../software/funwave/funwave--mpiifort-parallel-double_vWaveAvg_vWMnew_v1 input.txt)
    fi

    prev_jobid=$(echo "$submit_output" | awk '{print $4}')

    echo "Submitted $dir as JobID $prev_jobid"

done
