#!/bin/bash

# Params
#transcript_list=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/test_list.txt
transcript_list=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Transcript_lists/6.DS-Events-All-Tissues-Approach-1/Transcript-List-Prot-Coding-Approach-1-25Oct2021.txt
genome_fasta=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Data/Genome-fasta/GRCh38.p13.genome.fasta
gtf=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Data/GTF-Files/gencode.v26.annotation.gtf
outdir=/gpfs/projects/bsc83/Projects/RubenChazarra/1.Projects/1.Alternative_Splicing/1.AS_Function_Evaluator/Manual-Execution/1.Output/10-2021-22/1.Protein-FASTA-MN4/

while read transcript_id;
do
echo "${transcript_id}"
grep ${transcript_id} ${gtf} > ${outdir}/${transcript_id}.gtf
/home/bsc83/bsc83930/miniconda3/bin/gffread -g ${genome_fasta} ${outdir}/${transcript_id}.gtf -x ${outdir}/${transcript_id}.fasta -y ${outdir}/${transcript_id}.protein.fasta;
done < ${transcript_list}
