#!/bin/sh

#SBATCH --partition=interactive  # Partition to submit to
#SBATCH --workdir=/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/
#SBATCH --qos=debug
#SBATCH --ntasks=1

SUPPA_input="/home/bsc83/bsc83930/TFM-UOC-BSC/Data/Ribosomal_proteins.AS_events.ioe.rds"
TPM_counts="/home/bsc83/bsc83930/TFM-UOC-BSC/Data/Ribosomal_proteins.Transcripts.tpm.rds"
GTF_annot="/home/bsc83/bsc83930/TFM-UOC-BSC/Data/gencode.v26.annotation.gtf"

Rscript /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Process-SUPPA-output.R \
                --input_suppa ${SUPPA_input} \
                --count_matrix ${TPM_counts} \
                --annot_gtf ${GTF_annot} \
                --concordant "/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/results/0.Process-SUPPA/AS_events.Concordant_trancript_function.tsv" \
                --non_concordant "/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/results/0.Process-SUPPA/AS_events.NON_Concordant_trancript_function.tsv" \
                --transcript_ids "/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/results/0.Process-SUPPA/Transcript_list.Func_concordant.Coding.txt"
