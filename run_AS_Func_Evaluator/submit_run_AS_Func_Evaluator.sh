#!/bin/sh

#SBATCH --partition=interactive  # Partition to submit to
#SBATCH --workdir=/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/
#SBATCH --qos=debug


source_dir="/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator"
output_dir="${source_dir}/results"
report_dir="${output_dir}/reports"
 
/home/bsc83/bsc83930/miniconda3/bin/nextflow run main.nf "${source_dir}/main.nf" \
	-with-report "${report_dir}/report.html" \
        -with-trace "${report_dir}/trace.txt" \
        -with-timeline "${report_dir}/timeline.html" \
	-with-dag "${report_dir}/flowchart.png" \
	-profile slurm \
	#-resume

echo "Starting at `date`"
echo "Running on hosts: $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."
echo "Current working directory is `pwd`"	
