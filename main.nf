//#!/apps/NEXTFLOW/21.04.1/ nextflow

// Transcript Protein Domain Mapper ||  Developed by Ruben Chazarra-Gil (https://github.com/rubenchazarra)
// Inputs of the pipeline: i) list of protein coding transcript ids, ii) AS event visualisation aggregation file || Additinal files required: i) Annotation GTF file, ii) Genome.fasta, iii) Local installation of PFAM database, iv) UCSC Cytoband file for visualisation
// Pipeline: for each input transcript, we extract its coding sequence and translated it. With the protein sequence we query the PFAM protein domain database. The PFAM alignment relative coordinates are converted to genomic coordinates to enable visualisation of the transcript models together with their PFAM hits.

// TODO future: i) take output file name of aggr and visualisation processes, out of the script (clearer)

// Params
outdir = params.outdir
tag = params.tag
results_dir = "${outdir}/${tag}/"
r_module_load = params.r_module_load

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
	
	input:
	set val(transcript_id), file(genome), file(gtf) from data_ch	
	
	output:
	set val(transcript_id), file("${transcript_id}.gtf"), file("${transcript_id}.protein.fasta") into pfam_ch
	file("${transcript_id}.gtf") into collect_gtf_ch
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
	
	publishDir "${results_dir}/3.0.PFAM-Not-Parsed", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	input:
	set val(transcript_id), file(transcript_gtf), file(protein_fasta) from pfam_ch	
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-pfam-domtblout.txt") into pfam_parse_ch
	set val(transcript_id), file(protein_fasta) into protein_fasta_ch
	file("${transcript_id}-pfam-domtblout-with-alignment.txt")
	
	script:
	def hmmscan = params.pfam.hmmscan
	def pfam_db = params.pfam.pfam_db
	def sequence_evalue_thres = params.pfam.sequence_evalue_thres
	def pfam_cpu = params.pfam.pfam_cpu
	"""
	${hmmscan} --domtblout "${transcript_id}-pfam-domtblout.txt" \
		-E ${sequence_evalue_thres} \
		--cpu ${pfam_cpu} \
		${pfam_db} \
		${protein_fasta} > "${transcript_id}-pfam-domtblout-with-alignment.txt"
	# Usage: --domtbout (one line per domain); -E (report models <= this E-value threshold in output); --domE ( report domains <= this E-value threshold in output  [10.0])	
	"""
}

process parse_pfam {
	tag "Parse PFAM $transcript_id"
	
	publishDir "${results_dir}/3.1.PFAM-Not-Filtered", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al) from pfam_parse_ch	
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-pfam-not-filtered.txt") into pfam_no_filter_merge_ch, pfam_filter_ch
	
	script:
	
	"""
	${r_module_load}
	Hmmscan-dmtblout-Parser.R \
		--hmmscan_tbl ${pfam_al} \
		--rhmmer_path ${baseDir}/bin/rhmmer.R \
		--transcript_id ${transcript_id} \
		--output_table "${transcript_id}-pfam-not-filtered.txt" \
	"""
}

// Collect all PFAM alignments
merge_pfam_not_filtered = pfam_no_filter_merge_ch.map { it -> it[2] }.collect()

process filter_pfam {
	tag "Filter PFAM $transcript_id"
	
	publishDir "${results_dir}/3.2.PFAM-Filtered", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al) from pfam_filter_ch	
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-pfam-filtered.txt") into pfam_filtered_merge_ch, pfam_genomic_coord_ch, pfam_genomic_coord_gtf_ch, pfam_genomic_coord_ens_ch
	
	script:
	def domain_evalue = params.pfam.domain_evalue	
	def domain_score = params.pfam.domain_score	
	def accuracy = params.pfam.accuracy	
	def partiality = params.pfam.partiality
	"""
	${r_module_load}
	PFAM-Filter.R \
		--pfam_alignment ${pfam_al} \
		--domain_evalue $domain_evalue \
		--domain_score $domain_score \
		--accuracy $accuracy \
		--partiality $partiality \
		--pfam_filtered "${transcript_id}-pfam-filtered.txt"
	"""
}

// Collect all PFAM alignments
merge_pfam_filtered = pfam_filtered_merge_ch.map { it -> it[2] }.collect()

// Merge PFAM output
process merge_pfam {
	tag "Merge PFAM outputs" 
	publishDir "${results_dir}/3.3.PFAM-Filtered-Merged/",  mode: 'copy'
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	input:
	file("pfam/") from merge_pfam_filtered
		
	output:
	file("*txt") 
	script:
	"""
	${r_module_load}
	Merge-tables.R \
		--input_tables pfam/ \
		--output_df "PFAM-Filtered-Merged.txt"
	"""
}

// Add protein fasta to channel
pfam_genomic_coord_gtf_ch =  pfam_genomic_coord_ch. combine ( protein_fasta_ch, by: 0 ) 

// Extract PFAM alignment genomic coordinates [GTF]
process map_genomic_coord_gtf {
	tag "Genomic coords GTF: $transcript_id"
	
	publishDir "${results_dir}/4.1.PFAM-Genomic-Coords-GTF", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".txt") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al), file(protein_fasta) from pfam_genomic_coord_gtf_ch
	
	output:
	set val(transcript_id), file(transcript_gtf), file("${transcript_id}-PFAM-GenCoords-GTF.txt") into merge_gcoords_pfam_ch, vis_transcript_ch
	file("${transcript_id}-PFAM-GenCoords-GTF.txt") into collect_pfam_gc_ch
	
	script:
	"""
	${r_module_load}
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
	
	publishDir "${results_dir}/5.1.Visualisation-Transcript", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	when:
	params.vis & params.vis_transcript
		
	input:
	set val(transcript_id), file(transcript_gtf), file(pfam_al) from vis_transcript_ch
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	def cyto_band = params.visualisation.cyto_band	
	def genome_id = params.visualisation.genome_id	
	"""
	${r_module_load}
	Visualisation-Transcript-v2.R \
		--transcript_id ${transcript_id} \
		--pfam_genomic_coord ${pfam_al} \
		--gtf ${transcript_gtf} \
		--genome_id ${genome_id} \
		--cytoBand ${cyto_band} \
		--vis_track_plot "${transcript_id}_Gviz_Trackplot.pdf" \
		--vis_track_list "${transcript_id}_Gviz_Trackplot.rds" \
	"""
}



// Agrgegated visualisations (from various transcripts)

// 1) Read Event Coordinates. Columns are: event_id, event_coordinates
event_coord_ch = Channel.fromPath(params.visualisation.event_coords).splitCsv(header: false). map { tuple ( it[0], it[1..-1]) }

// 2) Read Event Aggregation file. Columns are: event_id, gene_name, transcript_ids
event_tr_ch = Channel.fromPath(params.visualisation.event_transcripts) . splitCsv(header: false)
// 3) Group transcript ids in tuple (N of transcripts can vary)
event_tr_reord_ch = event_tr_ch . map { tuple ( it[0], it[1], it[2..-1]) }

// 4) Merge Event transcripts + Event coordinates ch
event_aggr_ch = event_tr_reord_ch.combine ( event_coord_ch, by: 0 )

// 4) Collect GTF files
event_aggr_gtf_ch = collect_gtf_ch.collect()
// 5) Collect  PFAM Genomic Coordinate files
collect_pfam_gc_ch.collect().into{ merge_gcoords_gtf_ch; event_aggr_pfam_ch}

// ============ Discarded logic :(( Goal: enable event_aggregation without the need to wait for all transcripts to be run; Issue: file name collisions in aggregation step if duplicated event_ids, with different transcripts
// 4) i) Transpose ch (to have one ch per transcript); ii) Rearrange ch order (transcript_id, event_id, gene_id); and iii) Duplicate
// event_tr_reord_ch . transpose () . map { tuple( it[2], it[0], it[1]) } into { event_tr_1_ch; event_tr_2_ch }

// 5.1) Combine channels with GTF file. Cols: transcript_i, event_id, gene_id, gtf_i
// event_tr_gtf_ch = event_tr_1_ch . combine ( aggr_event_gtf_ch, by : 0)
// 5.2) Combine channels with PFAM file
// event_tr_pfam_ch = event_tr_2_ch . combine ( aggr_event_pfam_ch, by : 0)

// 6) Group transcript instances back by Event ID and gene_id, and reorder 
// event_tr_gtf_2_ch = event_tr_gtf_ch . groupTuple( by: [1, 2] ) . map{ tuple( it[1], it[2], it[0], it[3]) } 
// single_tr_pfam_2_ch = event_tr_pfam_ch . groupTuple( by: [1, 2] ) . map{ tuple( it[1], it[2], it[0], it[3]) }

// 7)  Integrate GTF and PFAM channels --> ISSUE HERE, groupTuple does not respect the order, how to ensure channel element order ?
// event_tr_gtf_pfam_ch = event_tr_gtf_2_ch . mix( single_tr_pfam_2_ch ). map {tuple( it[0], it[1], val(it[2]), it[3])} view()
//groupTuple( by : [0, 1]).view()

// 8) Rearrange ch. Now: event_id, gene_id, transcript_ids, gtfs, pfams
// event_tr_gtf_pfam_2_ch = event_tr_gtf_pfam_ch . map { tuple ( it[0], it[1], it[2][0], it[3][0], it[3][1] ) }

// 9) Merge file ch with event coordinates ch
// event_tr_coords_ch = event_tr_gtf_pfam_2_ch . combine( event_coord_ch, by: 0)
// =================

// Aggregate by Event
process aggr_event {
	tag "Aggr by Event $event_id"
	
	publishDir "${results_dir}/6.Domains-per-Event", mode: 'copy', pattern: "*_DF.rds"
	publishDir "${results_dir}/6.Domains-per-Event", mode: 'copy', pattern: "*_DF.txt"
	publishDir "${results_dir}/5.2.Visualisation-Events", mode: 'copy', pattern: "*_Visualisation.rds"
	        
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	when: 
	params.vis & params.vis_event
		
	input:
	set val(event_id), val(gene_id), val(transcript_ids), val(event_coords) from event_aggr_ch
	file("gtf_path/") from event_aggr_gtf_ch
	file("pfam_path/") from event_aggr_pfam_ch

	output:
	set val(event_id), val(gene_id), file("*_GTF_list.rds"), file("*_Visualisation.rds") into vis_event_ch
	file("*_DF.rds")
	file("*_DF.txt")
	
	script:
	"""
	${r_module_load}
	Aggregate-PFAM-by-Event.R \
		--transcript_ids "${transcript_ids}" \
		--event_id "${event_id}" \
		--gene_id ${gene_id} \
		--event_coords "${event_coords}" \
		--pfam_path "pfam_path/" \
		--gtf_path "gtf_path/"
	"""
}

// Visualize AS Event 
process visualisation_event {
	tag "Vis Event $event_id"
	
	publishDir "${results_dir}/5.2.Visualisation-Events", mode: 'copy',
	    saveAs: {filename ->
	    	if (filename.indexOf(".rds") > 0) "$filename"
	    	else if (filename.indexOf(".pdf") > 0) "$filename"
	    }
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	when: 
	params.vis & params.vis_event
		
	input:
	set val(event_id), val(gene_id), val(gtf_list), val(event_pfam_list) from vis_event_ch
	
	output:
	file("*.rds")
	file("*.pdf")
	
	script:
	def cyto_band = params.visualisation.cyto_band
	def genome_id = params.visualisation.genome_id
	def show_non_overlapping = params.visualisation.show_non_overlapping

	"""
	${r_module_load}
	Visualisation-AS-Event-v3.R \
		--event_pfam_list "${event_pfam_list}" \
		--gtf_list "${gtf_list}" \
		--event_id "${event_id}" \
		--gene_id ${gene_id} \
		--genome_id ${genome_id} \
		--cyto_band ${cyto_band} \
		--show_non_overlapping ${show_non_overlapping}
	"""
}

// Merge Genomic coordinates (from GTF)
process merge_gencoords_GTF {
	tag "Merge PFAM GenCoords GTF" 
	publishDir "${results_dir}/4.2.PFAM-Genomic-Coords-GTF-Merged/",  mode: 'copy'
	
	errorStrategy { task.exitStatus == 1 ? 'ignore' : 'ignore' }
	
	input:
	file("genomic-coord-pfam/") from merge_gcoords_gtf_ch
		
	output:
	file("*.txt") 
	script:
	"""
	${r_module_load}
	Merge-tables.R \
		--input_tables genomic-coord-pfam/ \
		--output_df "PFAM-Alignment-Genomic-Coordinates-GTF-Merged.txt"
	"""
}