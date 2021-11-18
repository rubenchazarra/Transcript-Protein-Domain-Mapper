#!/home/bsc83/bsc83930/miniconda3/bin nextflow

// Aternative Splicing Functional Consequence Evaluatori ||  Developed by Ruben Chazarra-Gil (https://github.com/rubenchazarra)

// Inputs of the pipeline: i) list of protein coding transcript ids, ii) AS event visualization aggregation file || Additinal files required: i) Annotation GTF file, ii) Genome.fasta, iii) Local installation of PFAM database, iv) UCSC Cytoband file for visualization

// Pipeline: for each input transcript, we extract its coding sequence and translated it. With the protein sequence we query the PFAM protein domain database. The PFAM alignment relative coordinates are converted to genomic coordinates to enable visualization of the transcript models together with their PFAM hits.

// Params
outdir = params.outdir
run_tag = params.run_tag

// Input files
// 1) Input GTF Annotation file
ch_GTF_annot = Channel.fromPath(params.input_files.annot_gtf).ifEmpty { exit 1, "Gencode GTF annotation file NOT found. Required!" }
// 2) Genome fasta
ch_genome_fasta = Channel.fromPath(params.input_files.genome).ifEmpty { exit 1, "Genome fasta file NOT found. Required!" }

// Read input transcript ID list
ch_transcriptID = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }.unique()

// Combine each transcript ID with the genome_fasta and GTF_file
ch_local_transcript_id = ch_transcriptID. combine( ch_genome_fasta ) .combine ( ch_GTF_annot )

// Get CDS sequence of each transcript and translate
process get_CDS_and_Protein_local {
	tag "Get CDS $transcript_id"
	
	publishDir "${outdir}/${run_tag}/", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf("${transcript_id}.fasta") > 0) "1.CDS-Fasta/$filename"
	    	else if (filename.indexOf("${transcript_id}.protein.fasta") > 0) "2.Protein-Fasta/$filename"
	    }
	
	errorStrategy 'ignore'
	
	input:
	set val(transcript_id), file(genome_fasta), file(GTF_file) from ch_local_transcript_id	
	
	output:
	set val(transcript_id), file("${transcript_id}.gtf"), file("${transcript_id}.protein.fasta") into ch_query_PFAM_local
	file("${transcript_id}.fasta")	
	
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
	tag "Query PFAM $transcript_id"
	
	publishDir "${outdir}/${run_tag}/3.PFAM-DomTblOut", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy 'ignore'
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(protein_fasta) from ch_query_PFAM_local	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-pfam-domtblout.txt") into ch_PFAM_output_local	
	set val(transcript_id), file(protein_fasta) into ch_protein_fasta 
	
	script:
	def local_PFAM_DB = params.query_PFAM.local_PFAM_DB	
	"""
	/home/bsc83/bsc83930/TFM-UOC-BSC/AS_Function_Evaluator/local-install-HMMER/hmmer-3.3.2/bin/hmmscan --domtblout "${transcript_id}-pfam-domtblout.txt" \
		-E 1e-5 \
		--cpu 4 \
		${local_PFAM_DB} \
		${protein_fasta} \
	# Usage: --domtbout (one line per domain); -E (report models <= this E-value threshold in output); --domE ( report domains <= this E-value threshold in output  [10.0])	
	"""
	}

process parse_PFAM_out {
	tag "Read PFAM $transcript_id"
	
	publishDir "${outdir}/${run_tag}/4.PFAM-Out-CSV", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy 'ignore'
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment) from ch_PFAM_output_local	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-pfam-alignment.txt") into ch_merge_PFAM_output, ch_genomic_coord_PFAM, ch_genomic_coord_PFAM_GTF
	
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Hmmscan-dmtblout-Parser.R \
		--hmmscan_tbl ${pfam_alignment} \
		--rhmmer_path ${baseDir}/bin/rhmmer.R \
		--transcript_id ${transcript_id} \
		--output_table "${transcript_id}-pfam-alignment.txt" \
	"""
	}

// Collect all PFAM alignments
ch_merge_PFAM_alignments = ch_merge_PFAM_output.map { it -> it[2] }.collect()

// Merge PFAM output
process merge_PFAM_output {
	tag "Merge PFAM outputs" 
	publishDir "${outdir}/${run_tag}/4.Merged-PFAM-output/",  mode: 'copy'
	
	errorStrategy 'ignore'
	
	input:
	//file("pfam/") from ch_merge_PFAM_alignments.collect()
	file("pfam/") from ch_merge_PFAM_alignments
		
	output:
	file("Merged-PFAM-output.txt") 
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Merge-tables.R \
		--input_tables pfam/ \
		--output_df "Merged-PFAM-output.txt"
	"""
}

// Extract PFAM alignment genomic coordinates [EnsemblDB]
process map_genomic_coord_ens {
	tag "Genomic coords EnsemblDB: $transcript_id"
	
	publishDir "${outdir}/${run_tag}/5.PFAM-GenCoords-EnsemblDB", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy 'ignore'
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment) from ch_genomic_coord_PFAM	
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-PFAM-GenCoords-EnsemblDB.txt") into ch_merge_gcoords_EnsemblDB
	
	script:
	"""
       	module load R
	${baseDir}/bin/Map-PFAM-to-genomic-coord-EnsemblDB.R \
		--pfam_alignment ${pfam_alignment} \
		--gtf ${transcript_GTF_file} \
		--transcript_id  ${transcript_id} \
		--output_coord  "${transcript_id}-PFAM-GenCoords-EnsemblDB.txt" \
	"""
}

// Add protein fasta to channel
ch_genomic_coord_PFAM_protein = ch_genomic_coord_PFAM_GTF.combine(ch_protein_fasta, by: 0) 

// Extract PFAM alignment genomic coordinates [GTF]
process map_genomic_coord_gtf {
	tag "Genomic coords GTF: $transcript_id"
	
	publishDir "${outdir}/${run_tag}/5.PFAM-GenCoords-GTF", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy 'ignore'
	
	input:
	set val(transcript_id), file(transcript_GTF_file), file(pfam_alignment), file(protein_fasta) from ch_genomic_coord_PFAM_protein
	
	output:
	set val(transcript_id), file(transcript_GTF_file), file("${transcript_id}-PFAM-GenCoords-GTF.txt") into ch_merge_gcoords_PFAM, ch_visualization_transcript, ch_visualization_event 
	
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Map-PFAM-to-genomic-coord-GTF.R \
		--pfam_alignment ${pfam_alignment} \
		--gtf ${transcript_GTF_file} \
		--transcript_id  ${transcript_id} \
		--protein_fasta  ${protein_fasta} \
		--output_coord  "${transcript_id}-PFAM-GenCoords-GTF.txt" \
	"""
	}

// Visualize Transcript Model + PFAM alignment
process visualization_transcript {
	tag "Visualize PFAM alignment $transcript_id"
	
	publishDir "${outdir}/${run_tag}/6.Visualization-Transcript", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	errorStrategy 'ignore'
	
	when:
	params.viz & params.viz_transcript
		
	input:
	set val(transcript_id), file(transcript_gtf_file), file(pfam_alignment) from ch_visualization_transcript
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Visualization-Transcript.R \
		--transcript_id  ${transcript_id} \
		--pfam_genomic_coord ${pfam_alignment} \
		--gtf ${transcript_gtf_file} \
		--cytoBand ${params.visualization.cytoBand_table} \
		--viz_track_plot  "${transcript_id}-Gviz-Trackplot.pdf" \
		--viz_track_list  "${transcript_id}-Gviz-Trackplot.rds" \
	"""
	}


// Visualization of Alternative Splicing Event

// Aggregation Visualizaiton Ch (from CSV). First element is Event_ID, next are Transcript_IDs participating in the event
ch_viz_aggr = Channel.fromPath(params.visualization.aggregation_csv).splitCsv(header: false).map { tuple ( it[0], it[1], it[2..-1]) }

// Duplicate viz_ch to select GTF files and PFAM outputs independently
ch_visualization_event.into { ch_viz_event_gtf; ch_viz_event_pfam }

// GTF Channel // Collect GTF Files
ch_viz_gtf = ch_viz_event_gtf.map{ it[1] }.collect()

// PFAM Gen-Coord Channel // Collect PFAM Genomic Coordinate Files 
ch_viz_pfam = ch_viz_event_pfam.map{ it[2] }.collect()

// Visualize AS Event // TODO: Control this process from config
process visualization_event {
	tag "Visualization AS Event $event_id"
	
	publishDir "${outdir}/${run_tag}/6.Visualization-Events", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	errorStrategy 'ignore'
	
	when: 
	params.viz & params.viz_event
		
	input:
	set val(event_id), val(gene_id), val(transcript_ids) from ch_viz_aggr 
	file ('gtf_path/*') from ch_viz_gtf	
	file ('pfam_path/*') from ch_viz_pfam	
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Visualization-AS-Event.R \
		--transcript_ids "${transcript_ids}" \
		--event_id "${event_id}" \
		--gene_id "${gene_id}" \
		--pfam_path "pfam_path/" \
		--gtf_path "gtf_path" \
		--cytoBand ${params.visualization.cytoBand_table} \
		--viz_track_list "${gene_id}-${event_id}-Gviz-Trackplot.rds" \
		--viz_track_plot "${gene_id}-${event_id}-Gviz-Trackplot.pdf" \
	"""
}

// Merge Genomic coordinates from GTF Approach
process merge_gencoords_GTF {
	tag "Merge PFAM GenCoords GTF" 
	publishDir "${outdir}/${run_tag}/5.PFAM-GenCoords-GTF-Merged/",  mode: 'copy'
	
	errorStrategy 'ignore'
	
	input:
	file("genomic-coord-pfam/") from ch_merge_gcoords_PFAM.collect()
		
	output:
	file("PFAM-GenCoords-GTF-Merged.txt") 
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Merge-tables.R \
		--input_tables genomic-coord-pfam/ \
		--output_df "PFAM-GenCoords-GTF-Merged.txt"
	"""
}

// Merge Genomic coordinates from EnsemblDB Approach
process merge_gencoords_EnsemblDB {
	tag "Merge PFAM GenCoords EnsemblDB" 
	publishDir "${outdir}/${run_tag}/5.PFAM-GenCoords-EnsemblDB-Merged/",  mode: 'copy'
	
	errorStrategy 'ignore'
	
	input:
	file("genomic-coord-pfam/") from ch_merge_gcoords_EnsemblDB.collect()
		
	output:
	file("PFAM-GenCoords-EnsemblDB-Merged.txt") 
	
	script:
	"""
       	module load R/3.6.1
	${baseDir}/bin/Merge-tables.R \
		--input_tables genomic-coord-pfam/ \
		--output_df "PFAM-GenCoords-EnsemblDB-Merged.txt"
	"""
}
