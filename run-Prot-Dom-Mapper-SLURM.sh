#!/bin/bash

#SBATCH --nodes=1
#SBATCH --job-name="DSE-Nord3"
#SBATCH --output=logs/slurm-%j.out
#SBATCH --error=logs/slurm-%j.err
#SBATCH --qos=bsc_ls
#SBATCH --ntasks-per-node=16
#SBATCH --cpus-per-task=1
#SBATCH --time "24:00:00"
 
# Run AS_Functional_Evaluator nextflow pipeline

source_dir="$(pwd)"
report_dir="${source_dir}/reports"

## Nextflow runtime runs on top of the Java virtual machine which, by design, tries to allocate as much memory as is available.
## To avoid this, specify the maximum amount of memory that can be used by the Java VM using the -Xms and -Xmx Java flags (https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html

# Unload all modules
# module purge # This is failing. If we unload everything, we cannot only load R, we require additional modules
# Load nextflow module
module load nextflow/21.04.1

export NXF_OPTS="-Xms500M -Xmx2G"

nextflow run "${source_dir}/main.nf" \
	-with-report "${report_dir}/report.html" \
        -with-trace "${report_dir}/trace.txt" \
        -with-timeline "${report_dir}/timeline.html" \
	-with-dag "${report_dir}/flowchart.png" \
	-profile slurm \
	#-resume awesome_banach
