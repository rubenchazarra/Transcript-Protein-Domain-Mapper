#!/bin/sh

#SBATCH --partition=interactive  # Partition to submit to
#SBATCH --qos=debug
#SBATCH --workdir=$(pwd)
#SBATCH --nodes=2
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=2

source_dir="$(pwd)"
report_dir="${source_dir}/reports"
 

## Nextflow runtime runs on top of the Java virtual machine which, by design, tries to allocate as much memory as is available.
## To avoid this, specify the maximum amount of memory that can be used by the Java VM using the -Xms and -Xmx Java flags (https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html

export NXF_OPTS="-Xms500M -Xmx2G"

/home/bsc83/bsc83930/miniconda3/bin/nextflow run main.nf "${source_dir}/main.nf" \
	-with-report "${report_dir}/report.html" \
        -with-trace "${report_dir}/trace.txt" \
        -with-timeline "${report_dir}/timeline.html" \
	-with-dag "${report_dir}/flowchart.png" \
	-profile slurm \
	-resume

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"	
