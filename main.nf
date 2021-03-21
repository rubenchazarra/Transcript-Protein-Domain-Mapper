#!/home/bsc83/bsc83930/miniconda3/bin nextflow

// This is a pipeline for systematic evaluatio of the functional consequences of AS events
// The idea is that the input will SUPPA's tool output (https://github.com/comprna/SUPPA) 
//process process_SUPPA_1{ 
//	
//	"""
//	module load R
//	Rscript bin/Process-SUPPA-output.R \
//		--input_suppa SUPPA_input \
//		--count_matrix TPM_counts \
//		--annot_gtf GTF_annot \
//		--concordant "AS_events.Concordant_trancript_function.tsv" \
//		--non_concordant "AS_events.NON_Concordant_trancript_function.tsv" \
//		--transcrip_ids "Transcript_list.Func_concordant.Coding.txt" 
//	"""
//
//}
ch_SUPPA_input = Channel
        .fromPath("${params.input_files.input_suppa}", checkIfExists: true)
        .ifEmpty { exit 1, "SUPPA input file found. Required!" }

ch_TPM_counts = Channel
        .fromPath("${params.input_files.tpm_counts}", checkIfExists: true)
        .ifEmpty { exit 1, "TPM counts file found. Required!" }

ch_GTF_annot = Channel
        .fromPath("${params.input_files.annot_gtf}", checkIfExists: true)
        .ifEmpty { exit 1, "Gencode GTF annotation file found. Required!" }

//ch_SUPPA_input.view()
//ch_TPM_counts.view()
//ch_GTF_annot.view()

// Process SUPPA output --> OUTCOMMENTING NOW TO RUN PIPELIN WITH PREVIOUSLY GENERATED TRANSCRIP LIST
//process process_SUPPA {
//	tag "Process SUPPA output"
//	publishDir "${params.outdir}/1.SUPPA_process", mode: 'copy',
//	    saveAs: {filename ->
//	    	if (filename.indexOf(".tsv") > 0) "$filename"
//	    	if (filename.indexOf(".txt") > 0) "$filename"
//	    }
//	
//	input:
//	set file(SUPPA_input) from ch_SUPPA_input
//	set file(TPM_counts) from ch_TPM_counts
//	set file(GTF_annot) from ch_GTF_annot
//	
//	output:
//	file("Transcript_list.Func_concordant.Coding.txt") into ch_transcript_list
//	
//	script:
//	//def SUPPA_input = params.input_files.input_suppa
//	//def TPM_counts = params.input_files.tpm_counts
//	//def GTF_annot = params.input_files.annot_gtf
//	"""
//	module load R
//	Rscript /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Process-SUPPA-output.R \
//		--input_suppa $SUPPA_input \\
//		--count_matrix $TPM_counts \\
//		--annot_gtf $GTF_annot \\
//		--concordant "AS_events.Concordant_trancript_function.tsv" \\
//		--non_concordant "AS_events.NON_Concordant_trancript_function.tsv" \\
//		--transcript_ids "Transcript_list.Func_concordant.Coding.txt"
//	"""
//
//}

// Read transcript ID list obtained from previous process 
//ch_transcriptID = ch_transcript_list.flatMap{ it.readLines() }
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }

//Get CDS sequences from trancript IDs from ENSEMBL REST API

process get_CDS {
	tag "get CDS $transcript_ID"
	publishDir "${params.outdir}/1.CDS_fasta", mode: 'copy'
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	memory = { 16.GB + 20.GB * (task.attempt) }		
	//when:
	
	input:
	val transcript_ID from ch_transcriptID	
	
	output:
	set val(transcript_ID), file("${transcript_ID}.fasta") into ch_translate_CDS	
	
	script:
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/get.CDS.Ensembl.py \
		-transcript_id $transcript_ID \
		-seq_type cds > "${transcript_ID}.fasta"
	"""
}

// Translate CDS sequence
process translate_CDS{
	tag "translate CDS $transcript_ID"
	publishDir "${params.outdir}/2.Protein_fasta", mode: 'copy'
	
	//when:
	
	input:
	set val(transcript_ID), file(fasta) from ch_translate_CDS	
	
	output:
	set val(transcript_ID), file("${transcript_ID}.protein.fasta") into ch_query_PFAM	
	
	script:
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/translate.CDS.py\
		-fasta $fasta > "${transcript_ID}.protein.fasta" 
	"""
}

//Query PFAM database
process query_PFAM{
	tag "query PFAM $transcript_ID"
	publishDir "${params.outdir}/3.PFAM_query", mode: 'copy'
	
	//when:
	
	input:
	set val(transcript_ID), file(fasta) from ch_query_PFAM	
	
	output:
	set val(transcript_ID), file("${transcript_ID}.out.txt") into ch_PFAM_output	
	file("${transcript_ID}.sequence.txt")
	file("${transcript_ID}.submission.params")
		
	script:
	def evalue_thres = params.query_PFAM.evalue_thres	
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/pfamscan-EBI.py\
		--email rc845@cam.ac.uk \
		--database pfam-a \
		--sequence $fasta \
		--evalue $evalue_thres \
		--format json \
		--outfile $transcript_ID
	"""
}

// Read PFAM output
//process read_PFAM_output{
//	tag "read PFAM output $transcript_ID"
//	//publishDir "${params.outdir}/3.PFAM_query", mode: 'copy'
//	
//	//when:
//	
//	input:
//	set val(transcript_ID), file(json_pfam) from ch_PFAM_output	
//	
//	output:
//		
//	script:
//	def evalue_thres = params.query_PFAM.evalue_thres	
//	"""
//        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/read.json.py\
//		-json_file $json_pfam 
//	"""
//}
