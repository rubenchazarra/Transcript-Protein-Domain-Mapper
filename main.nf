#!/home/bsc83/bsc83930/miniconda3/bin nextflow

// Aternative Splicing Functional Consequence Evaluator
// This pipeline takes as input a list of transcript IDs, extracts its protein sequence and aligns it to a local installation of the PFAM database. Finally it generates a visualization of the transcript Model and it's PFAM alignment
// After this, downstream processing requires to group transcript outputs by their AS event ID to inspect the consequences of the AS event at the PFAM-domain level


// Input GTF Annotation file
ch_GTF_annot = Channel.fromPath("${params.input_files.annot_gtf}")
        .ifEmpty { exit 1, "Gencode GTF annotation file NOT found. Required!" }

// Genome fasta
ch_genome_fasta = Channel.fromPath("${params.input_files.genome}")
        .ifEmpty { exit 1, "Genome fasta file NOT found. Required!" }


// Read transcript ID list obtained from previous process 
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }

//Get CDS sequences from trancript IDs from ENSEMBL REST API
if(params.approach == "interactive"){

process get_CDS_Ensembl_REST_API {
	tag "get CDS $transcript_id Ensembl REST API"
	publishDir "${params.outdir}/${params.run_tag}/1.CDS_fasta-Ensembl-REST-API", mode: 'copy'
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	// maxForks 8	
	//when:
	
	input:
	val transcript_id from ch_transcriptID	
	
	output:
	set val(transcript_id), file("${transcript_id}.fasta") into ch_translate_CDS	
	
	script:
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/get.CDS.Ensembl.py \
		-transcript_id $transcript_id \
		-seq_type cds > "${transcript_id}.fasta"
	"""
}

// Translate CDS sequence
process translate_CDS{
	tag "translate CDS $transcript_id"
	publishDir "${params.outdir}/${params.run_tag}/2.Protein_fasta-BioPython", mode: 'copy'
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	
	//when:
	
	input:
	set val(transcript_id), file(fasta) from ch_translate_CDS	
	
	output:
	set val(transcript_id), file("${transcript_id}.protein.fasta") into ch_query_PFAM	
	
	script:
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/translate.CDS.py\
		-fasta $fasta > "${transcript_id}.protein.fasta" 
	"""
}

//Query PFAM database
process query_PFAM{
	tag "query PFAM $transcript_id"
	publishDir "${params.outdir}/${params.run_tag}/3.PFAM_query-REST-API", mode: 'copy'
	// maxForks 8	
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	
	//when:
	
	input:
	set val(transcript_id), file(fasta) from ch_query_PFAM	
	
	output:
	set val(transcript_id), file("${transcript_id}.out.txt") into ch_PFAM_output	
	file("${transcript_id}.sequence.txt")
	file("${transcript_id}.submission.params")
		
	script:
	def evalue_thres = params.query_PFAM.evalue_thres	
	"""
        /home/bsc83/bsc83930/miniconda3/bin/python3 /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/pfamscan-EBI.py\
		--email rc845@cam.ac.uk \
		--database pfam-a \
		--sequence $fasta \
		--evalue $evalue_thres \
		--format json \
		--outfile $transcript_id
	"""
	}
// Read PFAM output
process read_PFAM_output{
	tag "read PFAM output $transcript_id"
	publishDir "${params.outdir}/${params.run_tag}/4.PFAM_output_CSV",  mode: 'copy'
	// maxForks 8	
	// memory = { 4.GB + 2.GB * (task.attempt) }		
	//when:
	
	input:
	set val(transcript_id), file(json_pfam) from ch_PFAM_output	
	
	output:
	file("${transcript_id}-pfam.alignment.txt") into ch_merge_PFAM_output		
	script:
	"""
       	module load R 
	Rscript /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/PFAM-output-JSON-to-csv.R \
		--input_json $json_pfam \
		--transcript_id $transcript_id \
		--output_table "${transcript_id}-pfam.alignment.txt"
	"""
	}


} else { // closing bracket from approach consition


// We have to do a bit of channel engineering here: 
// Combine each transcript ID with the genome_fasta and GTF_file
ch_local_transcript_id = ch_transcriptID. combine( ch_genome_fasta ) .combine ( ch_GTF_annot )

// Retrieve CDS sequence and protein sequence locally (no REST API)
process get_CDS_and_Protein_local {
	tag "get CDS $transcript_id Local"
	
	// TODO: FIX: All files are going to 1.CDS_fasta-Local since .fasta suffix is present in the 2 files	
	publishDir "${params.outdir}/${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".protein.fasta") > 0) "2.Protein_fasta-Local/$filename"
	    	else  "1.CDS_fasta-Local/$filename" 
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	// maxForks 8	
	
	input:
	set val(transcript_id), file(genome_fasta), file(GTF_file) from ch_local_transcript_id	
	
	output:
	set val(transcript_id), file("${transcript_id}.gtf"), file("${transcript_id}.protein.fasta") into ch_query_PFAM_local	
	set val(transcript_id), file("${transcript_id}.fasta") 
	
	script:
	"""
	# 1. Subset GTF file
	grep ${transcript_id} ${GTF_file} > ${transcript_id}.gtf
	## TODO --> Improve protein_ID extraction procedure from GTF
	# 2. Extract protein ID
	# TODO --> Change this AWK approach since it has give us a lot of headaches 
	#protein_ID=\$(awk '\$3 == "transcript" { print \$24 }' ${transcript_id}.gtf)   
	# 2. Extract CDS and protein sequence from $genome_fasta 
	/home/bsc83/bsc83930/miniconda3/bin/gffread -g ${genome_fasta} ${transcript_id}.gtf \
		-x ${transcript_id}.fasta \
		-y ${transcript_id}.protein.fasta \
	"""
	}

process query_PFAM_local {
	tag "Query PFAM $transcript_id Local"
	
	publishDir "${params.outdir}/${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "3.PFAM_query-Local/$filename"
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	// maxForks 8	
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(protein_fasta) from ch_query_PFAM_local	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}.pfam.out.txt") into ch_PFAM_output_local	
	
	script:
	def local_PFAM_DB = params.query_PFAM.local_PFAM_DB	
	"""
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/local-install-HMMER/hmmer-3.3.2/bin/hmmscan --domtblout "${transcript_id}.pfam.out.txt" \
		-E 1e-5 \
		--cpu 4 \
		${local_PFAM_DB} \
		${protein_fasta} \
	"""
	}

process read_PFAM_local {
	tag "Read PFAM $transcript_id Local"
	
	publishDir "${params.outdir}/${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "4.PFAM-output-CSV/$filename"
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	// maxForks 8	
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment) from ch_PFAM_output_local	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-pfam.alignment.txt") into ch_merge_PFAM_output
	
	script:
	"""
       	module load R 
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/PFAM-output-dmtblout-to-csv.R \
		--input_file  ${pfam_alignment} \
		--transcript_id ${transcript_id} \
		--output_table  "${transcript_id}-pfam.alignment.txt" \
	"""
	}

// Channel duplication 
ch_merge_PFAM_output.into{ ch_genomic_coord_PFAM; ch_merge_PFAM_alignments }

//// Merge PFAM output
process merge_PFAM_output{
	tag "Merge PFAM outputs" 
	publishDir "${params.outdir}/${params.run_tag}/4.Merged_PFAM_output/",  mode: 'copy'
	// maxForks 8	
	// memory = { 1.GB + 2.GB * (task.attempt) }		
	
	
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

// Extract PFAM alignment genomic coordinates
process map_genomic_coord {
	tag "Extract genomic coordinates $transcript_id"
	
	publishDir "${params.outdir}/${params.run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "5.Genomic-coord-PFAM/$filename"
	    }
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	// maxForks 8	
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment) from ch_genomic_coord_PFAM	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-PFAM.genomic.coordinates.txt") into ch_merge_genomic_coord_PFAM, ch_visualization
	
	script:
	"""
       	module load R 
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Map-PFAM-genomic-coord.R \
		--pfam_alignment ${pfam_alignment} \
		--gtf ${transcript_GTF_file} \
		--transcript_id  ${transcript_id} \
		--output_coord  "${transcript_id}-PFAM.genomic.coordinates.txt" \
	"""
	}
}

// Visualize PFAM alignment
process visualization {
	tag "Visualize PFAM alignment $transcript_id"
	
	publishDir "${params.outdir}/${params.run_tag}/6.Visualization", mode: 'copy'
	    //saveAs: {filename ->
	    //	if (filename.indexOf(".rds") > 0) "6.Visualization/$filename"
	    //	else if (filename.indexOf(".txt") > 0) "6.Visualization/$filename"
	    //}
	
	MAX = 4
	errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	// maxforks 8
	
	input:
	set val(transcript_id), file(transcript_gtf_file), file(pfam_alignment) from ch_visualization
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	"""
       	module load R 
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Visualization-Tracks-GTF.R \
		--transcript_id  ${transcript_id} \
		--pfam_genomic_coord ${pfam_alignment} \
		--gtf ${transcript_gtf_file} \
		--cytoBand ${params.visualization.cytoBand_table} \
		--viz_track_plot  "${transcript_id}-Visualization-Gviz-Trackplot.pdf" \
		--viz_track_list  "${transcript_id}-Visualization-Gviz-Trackplot.rds" \
	"""
	}

//// Merge Genomic coordinates output
process merge_genomic_coord_PFAM {
	tag "Merge PFAM outputs" 
	publishDir "${params.outdir}/${params.run_tag}/5.Merged-Genomic-coord-PFAM/",  mode: 'copy'
	// maxForks 8	
	// memory = { 6.GB + 2.GB * (task.attempt) }		
	
	
	input:
	file("genomic-coord-pfam/") from ch_merge_genomic_coord_PFAM.collect()
		
	output:
	file("Merged-Genomic-coord-PFAM-output.txt") 
	script:
	"""
       	module load R 
	Rscript /home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/bin/Merge-PFAM-alignments.R \
		--input_pfam genomic-coord-pfam/ \
		--output_pfam_df "Merged-Genomic-coord-PFAM-output.txt"
	"""
}
