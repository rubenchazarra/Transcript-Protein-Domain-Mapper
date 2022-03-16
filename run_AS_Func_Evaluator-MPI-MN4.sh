#!/bin/bash

# Run AS_Functional_Evaluator nextflow pipeline

source_dir="$(pwd)"
report_dir="${source_dir}/reports"

## Nextflow runtime runs on top of the Java virtual machine which, by design, tries to allocate as much memory as is available.
## To avoid this, specify the maximum amount of memory that can be used by the Java VM using the -Xms and -Xmx Java flags (https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html

export NXF_OPTS="-Xms500M -Xmx2G"

# For execution as MPI (for many low-intensive processes)
export NXF_CLUSTER_SEED=$(shuf -i 0-16777216 -n 1)

srun /apps/NEXTFLOW/21.04.1/nextflow run "${source_dir}/main.nf" \
	-with-report "${report_dir}/report.html" \
        -with-trace "${report_dir}/trace.txt" \
        -with-timeline "${report_dir}/timeline.html" \
	-with-dag "${report_dir}/flowchart.png" \
	-with-mpi \
	-profile slurm \
	-resume modest_hawking
