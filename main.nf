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

// Input from SUPPA
ch_SUPPA_input = Channel
        .fromPath("${params.input_files.input_suppa}", checkIfExists: true)
        .ifEmpty { exit 1, "SUPPA input file NOT found. Required!" }

// Input TPM matrix
ch_TPM_counts = Channel
        .fromPath("${params.input_files.tpm_counts}", checkIfExists: true)
        .ifEmpty { exit 1, "TPM counts file NOT found. Required!" }

// Input GTF Annotation file
ch_GTF_annot = Channel.fromPath("${params.input_files.annot_gtf}")
        .ifEmpty { exit 1, "Gencode GTF annotation file NOT found. Required!" }

// Genome fasta
ch_genome_fasta = Channel.fromPath("${params.input_files.genome}")
        .ifEmpty { exit 1, "Genome fasta file NOT found. Required!" }

// Process SUPPA output --> OUTCOMMENTING NOW TO RUN PIPELIN WITH PREVIOUSLY GENERATED TRANSCRIP LIST
//process process_SUPPA {
//	tag "Process SUPPA output"
//	publishDir "${params.outdir}/results-${params.run_tag}/1.SUPPA_process", mode: 'copy',
//	    saveAs: {filename ->
//	    	if (filename.indexOf(".tsv") > 0) "$filename"
//	    	if (filename.indexOf(".txt") > 0) "$filename"
//	    }
//	MAX = 4
//	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
//	memory = { 16.GB + 20.GB * (task.attempt) }		
//	
//	input:
//	file(SUPPA_input) from ch_SUPPA_input
//	file(TPM_counts) from ch_TPM_counts
//	file(GTF_annot) from ch_GTF_annot
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
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }

//Get CDS sequences from trancript IDs from ENSEMBL REST API
if(params.approach == "interactive"){

process get_CDS_Ensembl_REST_API {
	tag "get CDS $transcript_ID Ensembl REST API"
	publishDir "${params.outdir}/results-${params.run_tag}/1.CDS_fasta-Ensembl-REST-API", mode: 'copy'
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	memory = { 6.GB + 2.GB * (task.attempt) }		
	maxForks 1	
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
	publishDir "${params.outdir}/results-${params.run_tag}/2.Protein_fasta-BioPython", mode: 'copy'
	memory = { 6.GB + 2.GB * (task.attempt) }		
	
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
	publishDir "${params.outdir}/results-${params.run_tag}/3.PFAM_query-REST-API", mode: 'copy'
	maxForks 4	
	memory = { 4.GB + 2.GB * (task.attempt) }		
	
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
	publishDir "${params.outdir}/results-${params.run_tag}/4.PFAM_output_CSV",  mode: 'copy'
	maxForks 1	
	memory = { 2.GB + 2.GB * (task.attempt) }		
	//when:
	
	input:
	set val(transcript_ID), file(json_pfam) from ch_PFAM_output	
	
	output:
	file("${transcript_ID}-pfam.alignment.txt") into ch_merge_PFAM_output		
	script:
	"""
       	module load R 
	Rscript /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/PFAM-output-JSON-to-csv.R \
		--input_json $json_pfam \
		--transcript_id $transcript_ID \
		--output_table "${transcript_ID}-pfam.alignment.txt"
	"""
	}


} else { // closing bracket from approach consition


// We have to do a bit of channel engineering here: 
// Combine each transcript ID with the genome_fasta and GTF_file
ch_local_transcript_ID = ch_transcriptID. combine( ch_genome_fasta ) .combine ( ch_GTF_annot )

// Retrieve CDS sequence and protein sequence locally (no REST API)
process get_CDS_and_Protein_local {
	tag "get CDS $transcript_ID Local"
	
	// TODO: FIX: All files are going to 1.CDS_fasta-Local since .fasta suffix is present in the 2 files	
	publishDir "${params.outdir}/results-${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".fasta") > 0) "1.CDS_fasta-Local/$filename"
	    	else if (filename.indexOf(".protein.fasta") > 0) "2.Protein_fasta-Local/$filename" 
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	memory = { 6.GB + 2.GB * (task.attempt) }		
	maxForks 1	
	
	input:
	set val(transcript_ID), file(genome_fasta), file(GTF_file) from ch_local_transcript_ID	
	
	output:
	set val(transcript_ID), val(protein_ID), file("${transcript_ID}.protein.fasta") into ch_query_PFAM_local	
	set val(transcript_ID), file("${transcript_ID}.fasta") 
	
	script:
	"""
	# 1. Subset GTF file
	grep ${transcript_ID} ${GTF_file} > ${transcript_ID}.gtf
	## TODO --> Improve protein_ID extraction procedure from GTF
	# 2. Extract protein ID
	# TODO --> Change this AWK approachj to avoid token recognition error at: '('	
	protein_ID=$(awk '$3 == "transcript" { print $24 }' ${transcript_ID}.gtf)
	protein_ID=$(echo $protein_id |  grep -o -P '(?<=").*(?=\.)')	 
	# 2. Extract CDS and protein sequence from $genome_fasta 
	/home/bsc83/bsc83930/miniconda3/bin/gffread -g ${genome_fasta} ${transcript_ID}.gtf \
		-x ${transcript_ID}.fasta \
		-y ${transcript_ID}.protein.fasta \
	"""
	}

process query_PFAM_local {
	tag "Query PFAM $transcript_ID Local"
	
	publishDir "${params.outdir}/results-${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "3.PFAM_query-Local/$filename"
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	memory = { 6.GB + 2.GB * (task.attempt) }		
	maxForks 1	
	
	input:
	set val(transcript_ID), val(protein_ID), file(protein_fasta) from ch_query_PFAM_local	
	
	output:
	set val(transcript_ID), val(protein_ID), file("${transcript_ID}.pfam.out.txt") into ch_PFAM_output_local	
	
	script:
	def local_PFAM_DB = params.query_PFAM.local_PFAM_DB	
	"""
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/local-install-HMMER/hmmer-3.3.2/bin/hmmscan --domtblout "${transcript_ID}.pfam.out.txt" \
		-E 1e-5 \
		--cpu 4 \
		${local_PFAM_DB} \
		${protein_fasta} \
	"""
	}

process read_PFAM_local {
	tag "Read PFAM $transcript_ID Local"
	
	publishDir "${params.outdir}/results-${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "4.PFAM-output-CSV/$filename"
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	memory = { 6.GB + 2.GB * (task.attempt) }		
	maxForks 1	
	
	input:
	set val(transcript_ID), val(protein_ID), file(pfam_alignment) from ch_PFAM_output_local	
	
	output:
	set val(transcript_ID), val(protein_ID), file("${transcript_ID}-pfam.alignment.txt") into ch_merge_PFAM_output
	
	script:
	"""
       	module load R 
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/PFAM-output-dmtblout-to-csv.R \
		--input_file  ${pfam_alignment} \
		--transcript_id ${transcript_ID} \
		--output_table  "${transcript_ID}-pfam.alignment.txt" \
	"""
	}

// Channel duplication 
ch_merge_PFAM_output.into{ ch_genomic_coord_PFAM; ch_merge_PFAM_alignments }

// Extract PFAM alignment genomic coordinates
process extract_genomic_coord {
	tag "Extract genomic coordinates $transcript_ID"
	
	publishDir "${params.outdir}/results-${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "5.Genomic-coord-PFAM/$filename"
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	memory = { 6.GB + 2.GB * (task.attempt) }		
	maxForks 1	
	
	input:
	set val(transcript_ID), val(protein_ID), file(pfam_alignment) from ch_genomic_coord_PFAM	
	
	output:
	set val(transcript_ID), val(protein_ID), file("${transcript_ID}-PFAM.genomic.coordinates.txt") into ch_merge_PFAM_coord
	
	script:
	"""
       	module load R 
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Map-PFAM-genomic-coord.R \
		--pfam_alignment ${pfam_alignment} \
		--protein_id  ${protein_id} \
		--output_coord  "${transcript_ID}-PFAM.genomic.coordinates.txt" \
	"""
	}
}

//// Merge PFAM output
process merge_PFAM_output{
	tag "Merge PFAM outputs" 
	publishDir "${params.outdir}/results-${params.run_tag}/5.Merged_PFAM_output/",  mode: 'copy'
	maxForks 1	
	memory = { 1.GB + 2.GB * (task.attempt) }		
	
	
	input:
	file("pfam/") from ch_merge_PFAM_alignments.collect()
		
	output:
	file("Merged-PFAM-output.txt") 
	script:
	"""
       	module load R 
	Rscript /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Merge-PFAM-alignments.R \
		--input_pfam pfam/ \
		--output_pfam_df "Merged-PFAM-output.txt"
	"""
}
