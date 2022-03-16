#!/apps/NEXTFLOW/21.04.1/ nextflow

// Transcript Protein Domain Mapper ||  Developed by Ruben Chazarra-Gil (https://github.com/rubenchazarra)

// Inputs of the pipeline: i) list of protein coding transcript ids, ii) AS event visualisation aggregation file || Additinal files required: i) Annotation GTF file, ii) Genome.fasta, iii) Local installation of PFAM database, iv) UCSC Cytoband file for visualisation

// Pipeline: for each input transcript, we extract its coding sequence and translated it. With the protein sequence we query the PFAM protein domain database. The PFAM alignment relative coordinates are converted to genomic coordinates to enable visualisation of the transcript models together with their PFAM hits.

// Params
outdir = params.outdir
tag = params.tag
results_dir = "${outdir}/${tag}/"

// Input files
// 1) GTF Annotation file
gtf_ch = Channel.fromPath(params.gtf).ifEmpty { exit 1, "Gencode GTF annotation file NOT found. Required!" }
// 2) Genome fasta
genome_ch = Channel.fromPath(params.genome).ifEmpty { exit 1, "Genome fasta file NOT found. Required!" }

// Read input transcript ID list
transcripts_ch = Channel.fromPath(params.transcript_list).flatMap{ it.readLines() }.unique()

// Combine each transcript ID with the genome and gtf
data_ch = transcripts_ch. combine( genome_ch ) .combine ( gtf_ch )

// Extract CDS sequence and translate
process get_cds_translate {
	tag "Get CDS $transcript_id"
	
	publishDir "${results_dir}", mode: 'copy',
	    saveAs: { filename ->
	    	if (filename == "${transcript_id}.fasta") "1.CDS-Fasta/$filename"
	    	else if (filename == "${transcript_id}.protein.fasta") "2.Protein-Fasta/$filename"
	    	else if (filename.indexOf(".gtf") > 0) "0.Transcript-GTF/$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	set val(transcript_id), file(genome), file(gtf) from data_ch	
	
	output:
	set val(transcript_id), file("${transcript_id}.gtf"), file("${transcript_id}.protein.fasta") into pfam_ch
	file("${transcript_id}.fasta")	
	
	script:
	def gffread = params.get_cds_translate.gffread	
	"""
	# 1. Subset GTF file
	grep -w ${transcript_id} ${gtf} > ${transcript_id}.gtf
	# 2. Extract CDS and protein sequence from $genome 
	${gffread} -g ${genome} ${transcript_id}.gtf \
		-x ${transcript_id}.fasta \
		-y ${transcript_id}.protein.fasta \
	# Usage: -g // Genome FASTA, -x // Write a FASTA file with spliced CDS for each GFF transcript, -y // Write a protein FASTA file with the translation of CDS for each record 	
	"""
	}

process pfam {
	tag "PFAM $transcript_id"
	
	publishDir "${results_dir}/3.PFAM-DomTblOut", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	set val(transcript_id), file(transcript_gtf), file(protein_fasta) from pfam_ch	
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-pfam-domtblout.txt") into pfam_parse_ch
	set val(transcript_id), file(protein_fasta) into protein_fasta_ch
	file("${transcript_id}-pfam-domtblout-with-alignment.txt")
	
	script:
	def hmmscan = params.pfam.hmmscan
	def pfam_db = params.pfam.pfam_db
	def seq_e_val_thres = params.pfam.seq_e_val_thres
	def dom_e_val_thres = params.pfam.dom_e_val_thres
	"""
	${hmmscan} --domtblout "${transcript_id}-pfam-domtblout.txt" \
		-E ${seq_e_val_thres} \
		--domE ${dom_e_val_thres} \
		--cpu 2 \
		${pfam_db} \
		${protein_fasta} > "${transcript_id}-pfam-domtblout-with-alignment.txt"
	# Usage: --domtbout (one line per domain); -E (report models <= this E-value threshold in output); --domE ( report domains <= this E-value threshold in output  [10.0])	
	"""
	}

process parse_pfam {
	tag "Parse PFAM $transcript_id"
	
	publishDir "${results_dir}/4.PFAM-DomTblOut-CSV", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al) from pfam_parse_ch	
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-pfam-alignment.txt") into merge_pfam_ch, pfam_genomic_coord_ch, pfam_genomic_coord_gtf_ch, pfam_genomic_coord_ens_ch
	
	script:
	"""
	module load R/3.6.1
	Hmmscan-dmtblout-Parser.R \
		--hmmscan_tbl ${pfam_al} \
		--rhmmer_path ${baseDir}/bin/rhmmer.R \
		--transcript_id ${transcript_id} \
		--output_table "${transcript_id}-pfam-alignment.txt" \
	"""
	}

// Collect all PFAM alignments
merge_pfam_al_ch = merge_pfam_ch.map { it -> it[2] }.collect()

// Merge PFAM output
process merge_pfam {
	tag "Merge PFAM outputs" 
	publishDir "${results_dir}/4.PFAM-DomTblOut-CSV-Merged/",  mode: 'copy'
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	file("pfam/") from merge_pfam_al_ch
		
	output:
	file("PFAM-DomTblOut-CSV-Merged.txt") 
	script:
	"""
       	module load R/3.6.1
	Merge-tables.R \
		--input_tables pfam/ \
		--output_df "PFAM-DomTblOut-CSV-Merged.txt"
	"""
}

// Extract PFAM alignment genomic coordinates [EnsemblDB]
process map_genomic_coord_ens {
	tag "Genomic Coords EnsDB: $transcript_id"
	
	publishDir "${results_dir}/5.PFAM-GenCoords-EnsemblDB", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al) from pfam_genomic_coord_ens_ch
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-PFAM-GenCoords-EnsemblDB.txt") into merge_gcoords_ens_ch
	
	script:
	"""
       	module load R/3.6.1
	Map-PFAM-to-genomic-coord-EnsemblDB.R \
		--pfam_al ${pfam_al} \
		--gtf ${transcript_gtf} \
		--transcript_id  ${transcript_id} \
		--output_coord  "${transcript_id}-PFAM-GenCoords-EnsemblDB.txt" \
	"""
}

// Add protein fasta to channel
pfam_genomic_coord_gtf_ch =  pfam_genomic_coord_ch. combine ( protein_fasta_ch, by: 0 ) 

// Extract PFAM alignment genomic coordinates [GTF]
process map_genomic_coord_gtf {
	tag "Genomic coords GTF: $transcript_id"
	
	publishDir "${results_dir}/5.PFAM-GenCoords-GTF", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al), file(protein_fasta) from pfam_genomic_coord_gtf_ch
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-PFAM-GenCoords-GTF.txt") into merge_gcoords_pfam_ch, visualisation_transcript_ch, visualisation_aggr_ch 
	
	script:
	"""
       	module load R/3.6.1
	Map-PFAM-to-genomic-coord-GTF.R \
		--pfam_al ${pfam_al} \
		--gtf ${transcript_gtf} \
		--transcript_id  ${transcript_id} \
		--protein_fasta  ${protein_fasta} \
		--output_coord  "${transcript_id}-PFAM-GenCoords-GTF.txt" \
	"""
	}

// Visualize Transcript Model + PFAM alignment
process visualisation_transcript {
	tag "Visualize PFAM alignment $transcript_id"
	
	publishDir "${results_dir}/6.Visualisation-Transcript", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	when:
	params.vis & params.vis_transcript
		
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al) from visualisation_transcript_ch
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	def cyto_band = params.visualisation.cyto_band	
	"""
       	module load R/3.6.1
	Visualisation-Transcript-v2.R \
		--transcript_id  ${transcript_id} \
		--pfam_genomic_coord ${pfam_al} \
		--gtf ${transcript_gtf} \
		--cytoBand ${params.visualisation.cytoBand_table} \
		--vis_track_plot  "${transcript_id}-Gviz-Trackplot.pdf" \
		--vis_track_list  "${transcript_id}-Gviz-Trackplot.rds" \
	"""
	}


// Agrgegated visualisations (from various transcripts)

// Aggregation Visualizaiton Ch (from CSV). First element is Event_ID, next are Transcript_IDs participating in the event
vis_aggr_ch = Channel.fromPath(params.visualisation.aggregation_csv).splitCsv(header: false).map { tuple ( it[0], it[1], it[2..-1]) }

// Duplicate vis_ch to select GTF files and PFAM outputs independently
visualisation_aggr_ch.into { vis_gtf_ch ; vis_pfam_ch }

// GTF Channel // Collect GTF Files
vis_gtf_all_ch = vis_gtf_ch.map{ it[1] }.collect()

// PFAM Gen-Coord Channel // Collect PFAM Genomic Coordinate Files 
vis_pfam_all_ch = vis_pfam_ch.map{ it[2] }.collect()

// Visualize AS Event // TODO: Control this process from config
process visualisation_event {
	tag "Visualization AS Event $event_id"
	
	publishDir "${results_dir}/6.Visualisation-Events", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	when: 
	params.vis & params.vis_event
		
	input:
	set val(event_id), val(gene_id), val(transcript_ids) from vis_aggr_ch 
	file ('gtf_path/*') from vis_gtf_all_ch	
	file ('pfam_path/*') from vis_pfam_all_ch	
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	"""
       	module load R/3.6.1
	Visualisation-AS-Event.R \
		--transcript_ids "${transcript_ids}" \
		--event_id "${event_id}" \
		--gene_id "${gene_id}" \
		--pfam_path "pfam_path/" \
		--gtf_path "gtf_path" \
		--cytoBand ${params.visualisation.cytoBand_table} \
		--vis_track_list "${gene_id}-${event_id}-Gviz-Trackplot.rds" \
		--vis_track_plot "${gene_id}-${event_id}-Gviz-Trackplot.pdf" \
	"""
}

// Collect GenCoords-GTF

merge_gcoords_gtf_ch = merge_gcoords_pfam_ch.map { it -> it[2] }.collect()

// Merge Genomic coordinates from GTF Approach
process merge_gencoords_GTF {
	tag "Merge PFAM GenCoords GTF" 
	publishDir "${results_dir}/5.PFAM-GenCoords-GTF-Merged/",  mode: 'copy'
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	file("genomic-coord-pfam/") from merge_gcoords_gtf_ch
		
	output:
	file("PFAM-GenCoords-GTF-Merged.txt") 
	script:
	"""
       	module load R/3.6.1
	Merge-tables.R \
		--input_tables genomic-coord-pfam/ \
		--output_df "PFAM-GenCoords-GTF-Merged.txt"
	"""
}

// Collect GenCoords-EnsemblDB
merge_gcoords_ens_ch = merge_gcoords_ens_ch.map { it -> it[2] }.collect()

// Merge Genomic coordinates from EnsemblDB Approach
process merge_gencoords_EnsDB {
	tag "Merge PFAM GenCoords EnsemblDB" 
	publishDir "${results_dir}/5.PFAM-GenCoords-EnsemblDB-Merged/",  mode: 'copy'
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	// memory = 2.GB	
	
	input:
	file("genomic-coord-pfam/") from merge_gcoords_ens_ch
	
	output:
	file("PFAM-GenCoords-EnsemblDB-Merged.txt") 
	
	script:
	"""
       	module load R/3.6.1
	Merge-tables.R \
		--input_tables genomic-coord-pfam/ \
		--output_df "PFAM-GenCoords-EnsemblDB-Merged.txt"
	"""
}
