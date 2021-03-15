#!/home/bsc83/bsc83930/miniconda3/bin nextflow

// This is a pipeline for systematic evaluatio of the functional consequences of AS events
// The idea is that the input will SUPPA's tool output (https://github.com/comprna/SUPPA) 

// Read transcript ID list [Note: no version specified ]
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }

//Get CDS sequences from trancript IDs from ENSEMBL REST API

process get_CDS {
	tag "get CDS $transcript_ID"
	publishDir "${params.outdir}/1.CDS_fasta", mode: 'copy'
	
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
process read_PFAM_output{
	tag "read PFAM output $transcript_ID"
	//publishDir "${params.outdir}/3.PFAM_query", mode: 'copy'
	
	//when:
	
	input:
	set val(transcript_ID), file(json_pfam) from ch_PFAM_output	
	
	output:
		
	script:
	def evalue_thres = params.query_PFAM.evalue_thres	
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/read.json.py\
		-json_file $json_pfam 
	"""
}
