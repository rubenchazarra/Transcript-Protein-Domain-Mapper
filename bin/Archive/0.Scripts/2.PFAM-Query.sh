#!/bin/bash

# Params
transcript_list=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Transcript_lists/6.DS-Events-All-Tissues-Approach-1/Transcript-List-Prot-Coding-Approach-1-25Oct2021.txt
# transcript_list=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/test_list.txt
filedir=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/1.Output/10-2021-22/1.Protein-FASTA/
pfam_db=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/local-install-PFAM/Pfam-A.hmm
outdir=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/1.Output/10-2021-22/2.PFAM-Query/

while read -r transcript_id; 
do
protein_fasta=${filedir}/${transcript_id}.protein.fasta
echo $transcript_id
/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/local-install-HMMER/hmmer-3.3.2/bin/hmmscan  --domtblout $outdir/${transcript_id}-pfam-domtblout.txt -E 1e-5 --cpu 4 ${pfam_db} ${protein_fasta}
done < ${transcript_list}
