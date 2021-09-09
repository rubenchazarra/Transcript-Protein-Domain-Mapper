#!/home/bsc83/bsc83930/miniconda3/bin nextflow

// Aternative Splicing Functional Consequence Evaluator
// This pipeline takes as input a list of transcript IDs, extracts its protein sequence and aligns it to a local installation of the PFAM database. Finally it generates a visualization of the transcript Model and it's PFAM alignment
// After this, downstream processing requires to group transcript outputs by their AS event ID to inspect the consequences of the AS event at the PFAM-domain level
// Developed by Ruben Chazarra-Gil (https://github.com/rubenchazarra)

// Params
gtf = params.input_files.annot_gtf
genome_fasta = params.input_files.genome

// Input GTF Annotation file
ch_GTF_annot = Channel.fromPath("${gtf}").ifEmpty { exit 1, "Gencode GTF annotation file NOT found. Required!" }
        
// Genome fasta
ch_genome_fasta = Channel.fromPath("${genome_fasta}").ifEmpty { exit 1, "Genome fasta file NOT found. Required!" }


// Read transcript ID list obtained from previous process 
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }


// Combine each transcript ID with the genome_fasta and GTF_file
ch_local_transcript_id = ch_transcriptID. combine( ch_genome_fasta ) .combine ( ch_GTF_annot )

// Retrieve CDS sequence and protein sequence locally (no REST API)
process get_CDS_and_Protein_local {
	tag "get CDS $transcript_id Local"
	
	// TODO: FIX:  Sub-GTF files are saved in 1.CDS_fasta-Local due to the "else" // If indexing ".indexOf(".fasta")" both transcript and protein fasta are saved in the same dir.	
	publishDir "${params.outdir}/${params.run_tag}/2.Protein_fasta-Local", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".protein.fasta") > 0) "$filename"
	    	else  "1.CDS_fasta-Local/$filename" 
	    }
	
	// Note: leaving this for potential future usage	
	// MAX = 4
	// errorStrategy { (task.exitStatus == 130 || task.exitStatus == 137) && task.attempt - 1 <= MAX ? 'retry' : 'ignore' }
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
	# 2. Extract CDS and protein sequence from $genome_fasta 
	/home/bsc83/bsc83930/miniconda3/bin/gffread -g ${genome_fasta} ${transcript_id}.gtf \
		-x ${transcript_id}.fasta \
		-y ${transcript_id}.protein.fasta \
	# Usage: -g // Genome FASTA, -x // Write a FASTA file with spliced CDS for each GFF transcript, -y // Write a protein FASTA file with the translation of CDS for each record 	
	"""
	}

process query_PFAM_local {
	tag "Query PFAM $transcript_id Local"
	
	publishDir "${params.outdir}/${params.run_tag}/3.PFAM_query-Local", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
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
	# Usage: --domtbout (one line per domain); -E (report models <= this E-value threshold in output); --domE ( report domains <= this E-value threshold in output  [10.0])	
	"""
	}

process read_PFAM_local {
	tag "Read PFAM $transcript_id Local"
	
	publishDir "${params.outdir}/${params.run_tag}/4.PFAM-output-CSV", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment) from ch_PFAM_output_local	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-pfam.alignment.txt") into ch_merge_PFAM_output
	
	script:
	"""
       	module load R 
	${baseDir}/bin/PFAM-domtblout-to-CSV.R \
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
	${baseDir}/bin/Merge-PFAM-alignments.R \
		--input_pfam pfam/ \
		--output_pfam_df "Merged-PFAM-output.txt"
	"""
}

// Extract PFAM alignment genomic coordinates
process map_genomic_coord {
	tag "Genomic coordinates $transcript_id"
	
	publishDir "${params.outdir}/${params.run_tag}/5.Genomic-coord-PFAM", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment) from ch_genomic_coord_PFAM	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-PFAM.genomic.coordinates.txt") into ch_merge_genomic_coord_PFAM, ch_visualization_transcript, ch_visualization_event 
	
	script:
	"""
       	module load R 
	${baseDir}/bin/Map-PFAM-genomic-coord.R \
		--pfam_alignment ${pfam_alignment} \
		--gtf ${transcript_GTF_file} \
		--transcript_id  ${transcript_id} \
		--output_coord  "${transcript_id}-PFAM.genomic.coordinates.txt" \
	"""
	}


// Visualize Transcript Model + PFAM alignment
process visualization_transcript {
	tag "Visualize PFAM alignment $transcript_id"
	
	publishDir "${params.outdir}/${params.run_tag}/6.Visualization-Transcript", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	input:
	set val(transcript_id), file(transcript_gtf_file), file(pfam_alignment) from ch_visualization_transcript
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	"""
       	module load R 
	${baseDir}/bin/Visualization-Transcript.R \
		--transcript_id  ${transcript_id} \
		--pfam_genomic_coord ${pfam_alignment} \
		--gtf ${transcript_gtf_file} \
		--cytoBand ${params.visualization.cytoBand_table} \
		--viz_track_plot  "${transcript_id}-Visualization-Gviz-Trackplot.pdf" \
		--viz_track_list  "${transcript_id}-Visualization-Gviz-Trackplot.rds" \
	"""
	}


// Visualization of Alternative Splicing Event

// Aggregation Visualizaiton Ch (from CSV). First element is Event_ID, next are Transcript_IDs participating in the event
ch_viz_aggr = Channel.fromPath(params.visualization.aggregation_csv).splitCsv(header: false).map { tuple ( it[0], it[1..-1]) }

// Duplicate viz_ch to select GTF files and PFAM outputs independently
ch_visualization_event.into { ch_viz_event_gtf; ch_viz_event_pfam }

// GTF Channel // Collect GTF Files
ch_viz_gtf = ch_viz_event_gtf.map{ it[1] }.collect()

// PFAM Gen-Coord Channel // Collect PFAM Genomic Coordinate Files 
ch_viz_pfam = ch_viz_event_pfam.map{ it[2] }.collect()

// Visualize AS Event
process visualization_event {
	tag "Visualization AS Event $event_id"
	
	publishDir "${params.outdir}/${params.run_tag}/6.Visualization-Events", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	input:
	tuple val(event_id), val(transcript_ids) from ch_viz_aggr 
	file ('gtf_path/*') from ch_viz_gtf	
	file ('pfam_path/*') from ch_viz_pfam	
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	"""
       	module load R 
	${baseDir}/bin/Visualization-AS-Event.R \
		--transcript_ids  "${transcript_ids}" \
		--event_id "${event_id}" \
		--pfam_path 'pfam_path/' \
		--gtf_path 'gtf_path/' \
		--cytoBand ${params.visualization.cytoBand_table} \
		--viz_track_list  "${event_id}-Visualization-Gviz-Trackplot.rds" \
		--viz_track_plot  "${event_id}-Visualization-Gviz-Trackplot.pdf" \
	"""
	}


//// Merge Genomic coordinates output
process merge_genomic_coord_PFAM {
	tag "Merge PFAM outputs" 
	publishDir "${params.outdir}/${params.run_tag}/5.Merged-Genomic-coord-PFAM/",  mode: 'copy'
	
	input:
	file("genomic-coord-pfam/") from ch_merge_genomic_coord_PFAM.collect()
		
	output:
	file("Merged-Genomic-coord-PFAM-output.txt") 
	script:
	"""
       	module load R 
	${baseDir}/bin/Merge-PFAM-alignments.R \
		--input_pfam genomic-coord-pfam/ \
		--output_pfam_df "Merged-Genomic-coord-PFAM-output.txt"
	"""
}
