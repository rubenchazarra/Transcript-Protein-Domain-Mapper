#!/home/bsc83/bsc83930/miniconda3/bin nextflow

// This is a pipeline for systematic evaluatio of the functional consequences of AS events
// The idea is that the input will SUPPA's tool output (https://github.com/comprna/SUPPA) 

// Read transcript ID list [Note: no version specified ]
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }

//Get CDS sequences from trancript IDs from ENSEMBL REST API

process get_CDS {
	tag "$samplename"
	publishDir "${params.outdir}/1.CDS_fasta", mode: 'copy'
	
	//when:
	
	input:
	set val(transcript_ID) from ch_transcriptID	
	
	output:
	set val(transcript_ID), file("${transcript_ID}.fasta") into ch_translate_CDS	
	
	script:
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/get.CDS.Ensembl-2.py \
		-transcript_id $transcript_ID \
		 -seq_type cds > "${transcript_ID}.fasta"	
	"""
}
