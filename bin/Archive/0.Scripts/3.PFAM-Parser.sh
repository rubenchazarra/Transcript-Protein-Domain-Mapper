#!/bin/bash

# Params
transcript_list=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Transcript_lists/6.DS-Events-All-Tissues-Approach-1/Transcript-List-Prot-Coding-Approach-1-25Oct2021.txt
#transcript_list=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/test_list.txt
filedir=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/1.Output/10-2021-22/2.PFAM-Query/
pfam_db=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/local-install-PFAM/Pfam-A.hmm
outdir=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/1.Output/10-2021-22/3.PFAM-Parsed/

# Get rid of already loaded modules
# module purge # This is failing

# Load required modules
# module load gcc/6.1.0
# module load CURL/7.49.0
# module load BZIP2/1.0.6
module load R/3.6.1

while read -r transcript_id; 
do
echo $transcript_id
pfam_tbl=${filedir}/${transcript_id}-pfam-domtblout.txt

/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/bin/Hmmscan-dmtblout-Parser.R --hmmscan_tbl ${pfam_tbl} --rhmmer_path /gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/bin/rhmmer.R --transcript_id $transcript_id --output_table $outdir/${transcript_id}-pfam-alignment.txt 
done < ${transcript_list}
