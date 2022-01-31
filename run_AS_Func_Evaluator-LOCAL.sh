#/bin/bash

# Run AS_Functional_Evaluator nextflow pipeline

source_dir="$(pwd)"
report_dir="${source_dir}/reports"

## Nextflow runtime runs on top of the Java virtual machine which, by design, tries to allocate as much memory as is available.
## To avoid this, specify the maximum amount of memory that can be used by the Java VM using the -Xms and -Xmx Java flags (https://www.nextflow.io/blog/2021/5_tips_for_hpc_users.html

export NXF_OPTS="-Xms500M -Xmx2G"

module load nextflow
nextflow run "${source_dir}/main.nf" \
	-with-report "${report_dir}/report.html" \
        -with-trace "${report_dir}/trace.txt" \
        -with-timeline "${report_dir}/timeline.html" \
	-with-dag "${report_dir}/flowchart.png" \
	-profile local \
	-resume pedantic_torvalds
