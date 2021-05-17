#/bin/bash

# Run AS_Functional_Evaluator nextflow pipeline

source_dir="/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator"
report_dir="${source_dir}/reports"

#  export NXF_OPTS="-Xms50m -Xmx500m"
 
/home/bsc83/bsc83930/miniconda3/bin/nextflow run main.nf "${source_dir}/main.nf" \
	-with-report "${report_dir}/report.html" \
        -with-trace "${report_dir}/trace.txt" \
        -with-timeline "${report_dir}/timeline.html" \
	-with-dag "${report_dir}/flowchart.png" \
	-profile local
	-resume	
