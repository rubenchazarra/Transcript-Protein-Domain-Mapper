#!/bin/sh

#SBATCH --partition=interactive  # Partition to submit to
#SBATCH --qos=debug
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=2                  # Run a single task		
#SBATCH --cpus-per-task=4            # Number of CPU cores per task
#SBATCH --mem=32gb 
#SBATCH --output=logs/slurm_%j.out
#SBATCH --error=logs/slurm_%j.err

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
