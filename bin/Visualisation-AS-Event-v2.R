#!/usr/bin/env Rscript

# This Script Generates Visualization of Transcript Gene Model + the PFAM alignments of the corresponding Transcript
# The Goal is to represent the Alternative Splicing Event. The transcripts involved in each AS event are defined in an external Aggregation CSV
# TODO: Extract transcript genomic coordinates from GTF and not from EnsmblDB for consistency
## Vulnerability 1): Maybe using list.files when there are a ton of files in the path is not optimal
## Vul 2) : mapply() works for integrating 2 lists, if we have > 2 transcripts, this will fail. 

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--transcript_ids"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Transcript IDs corresponding to the AS event to be represented.'
  ),
  make_option(
    c("-e", "--event_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'The ID of the AS Event to be represented. This will appear in the PDF Plot Title.'
  ),
  make_option(
    c("-g", "--gene_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Gene ID (in HGNC Format) of the Alternative Splicing Event'
  ),
  make_option(
    c("-s", "--genomic_start"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = 'Genomic Start coordinates of the Event to represent (Will not be represented if input is NULL)'
  ), 
  make_option(
    c("-n", "--genomic_end"),
    action = "store",
    default = NULL,
    type = 'numeric',
    help = 'Genomic End coordinates of the Event to represent (Will not be represented if input is NULL)'
  ), 
  make_option(
    c("-a", "--genome_id"),
    action = "store",
    default = "hg38",
    type = 'character',
    help = 'UCSC Genome Browser assembly ID. Latest IDs for Human and Mouse genomes are "hg38", and "mm39" respectively'
  ),
  make_option(
    c("-p", "--pfam_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the directory where the PFAM alignments Genomic Coordinates of the differents Transcripts to represent are located.'
  ), 
  make_option(
    c("-t", "--gtf_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the directory where the GTF files of the transcripts to be represented are located (These are a subset of the GTF annotation file containing only the query transcript.'
  ), 
  make_option(
    c("-c", "--cytoBand"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the UCSC Browser Cyto Band table. To download cytoband files for genomes hosted at UCSC, see the UCSC Table Browser, and select Group="All Tables" and Table="cytoBand".'
  ), 
  make_option(
    c("-o", "--vis_track_list"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output RDS object storing lists of GViz Tracks for vizualization.'
  ), 
  make_option(
    c("-v", "--vis_track_plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the PDF file with the Track Plots.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Gviz))
suppressPackageStartupMessages(require(rtracklayer))
suppressPackageStartupMessages(require(dplyr))

## Functions ## 
gene.track.v2 <- function(gen.coords, which, genome, group.id, plot.label){
  suppressPackageStartupMessages(require(Gviz))
  ## Generate gene region track (either the full model of transcript or its PFAM domains)
  if(!which %in% c("transcript", "pfam")) stop(paste0("which arg must be one of: transcript or pfam", which, " not valid!" ))
  if(which == "transcript"){
    feat <- as.character(tolower(gen.coords$type))
  } else if(which == "pfam"){
    feat <- group.id
  }
  chr <- unique(as.character(gen.coords[["seqnames"]]))
  strand <- unique(as.character(gen.coords[["strand"]]))
  #Track
  gene.track <- GeneRegionTrack(rstarts = gen.coords[["start"]],
                                rends = gen.coords[["end"]],
                                genome = genome,
                                chromosome = chr,
                                strand = strand,
                                name = group.id, # grouping var (I think)
                                symbol = plot.label, # plot label (I think)
                                feature = feat, # distinguishes CDS from UTRs
                                cex.axis = 15)
  # Edit dpars (to enable distinction btw CDS and UTR)
  gene.track@dp@pars[["collapse"]] <- FALSE
  l <- list(gene.track)
  names(l) <- plot.label
  return(l)
}

parse.gtf.v2 <- function(gtf, transcript.id){
  ## Subset GTF and perform pertinent modifications for visualization by Gviz
  sub.gtf <- gtf[grep(transcript.id, gtf$transcript_id), ]
  # Create dataframe --
  transcript.model <- as.data.frame(sub.gtf)
  # Note: 'exon' type in gencode includes UTR+CDS
  transcript.model <- transcript.model[transcript.model$type %in% c("UTR", "CDS"),]
  # Rename columns --
  transcript.model$gene <- transcript.model$gene_name
  transcript.model$exon <- transcript.model$exon_id
  transcript.model$transcript <- transcript.model$transcript.id
  transcript.model$gene <- transcript.model$gene_id
  direction <- ifelse(transcript.model$strand == "+", "<", ">")
  transcript.model$symbol <- paste(transcript.model$transcript.id, direction) # label in plot
  return(transcript.model)
}


# 5. Gene Region Tracks
generate.transcript.track <- function(coord_df, gen, chr, label_name, id_name, color){
  # Generate Transcript Model Track
  
  ## Iterate across 'id_name' in coord_df. This is for Multi-Hit PFAM alignments
  track_list <- list()
  for(id in unique(coord_df[[id_name]])){
    track  <- GeneRegionTrack(coord_df[coord_df[[id_name]] == id, ], 
                              genome = gen, 
                              chromosome = chr, 
                              name = id, 
                              cex.axis = 7, 
                              cex.group= 1,
                              symbol = id
    )
    
    displayPars(track) <- list(showId=TRUE, background.panel=color)
    track_list[[id]] <- track
  }
  return(track_list)
}

parse.pfam <- function(pfam.gc){
  ## Coerce PFAM alignment output df to format readable by GViz package.
  direction <- ifelse(pfam.gc$strand == "+", "<", ">")
  # Basically, create a df [almost] equivalent to the 'transcript_model'
  df = data.frame("seqnames" = pfam.gc$seqnames,
                  "start" = pfam.gc$start,
                  "end" = pfam.gc$end,
                  "width" = pfam.gc$width,
                  "strand" = pfam.gc$strand,
                  "transcript" = pfam.gc$transcript_id,
                  "exon" = pfam.gc$exon_id,
                  "protein" = pfam.gc$protein_id,
                  "symbol" = paste(pfam.gc$transcript_id, direction),
                  #"pfam_alignment_id" = paste(pfam.gc$pfam_alignment_id, direction)
                  "pfam_alignment_id" = pfam.gc$pfam_alignment_id
  )
  # convert alignment ID to character  Required for visualization
  df[["symbol"]] = as.character(df[["pfam_alignment_id"]])
  return(df)
}

chr.track.fun <- function(genome, chr, cyto.band ){
  suppressPackageStartupMessages(require(Gviz))
  ## Generate chromosome track ## TODO: what if cyto.band is NULL (is not available ?)
  itrack <- IdeogramTrack(genome = genome,
                          chromosome = chr,
                          bands = cyto.band,
                          name = chr)
  list("Chr.track" = itrack)
}

genome.track.fun <- function(){
  suppressPackageStartupMessages(require(Gviz))
  # Generate genome track
  gtrack <- GenomeAxisTrack(name = "genome.track")
  list("Genome.Track" = gtrack)
}

event.track.fun <- function(gen.start, gen.end, chr, gen){
  # Generate Track with Coordinates of Aggregation Event
  ev.track <- GeneRegionTrack(rstarts = gen.start,
                              rends = gen.end,
                              genome = gen,
                              chromosome = chr,
                              name = "event.track", 
                              symbol = "Event Track",
                              feature = "event.track",
                              cex.axis = 15)
  
  list("Event.Track" = ev.track)
}

save.plot.pdf <- function(track.list, plot.title, file.name ){
  ## Save PDF plot
  pdf(file.name)
  # Plot
  Gviz::plotTracks(track.list,
                   event.track = "darkgreen",
                   cds = "blue", utr = "blue",
                   pfam = "orange", 
                   transcriptAnnotation = "symbol", 
                   main = plot.title, cex.main = 0.75)
  dev.off()
}

## Execute
# 0. Params
## 0.1. Input params
transcript_ids <- opt$transcript_ids
genome_id <- opt$genome_id
gene_id <- opt$gene_id
event_id <- opt$event_id
genomic_start <- as.numeric(opt$genomic_start) # If not numeric(no event coords provided), will return NAs
genomic_end <- as.numeric(opt$genomic_end)
## 0.2 Input file paths
pfam_path <- opt$pfam_path
gtf_path <- opt$gtf_path
## 0.3. Output file names
out.rds.file.name <- opt$vis_track_list
out.pdf.file.name <- opt$vis_track_plot
pfam.gen.coords <- opt$pfam_genomic_coord
## 0.3. Read Cytoband information from UCSC, required to generate the chromosome Ideogram without connection to the UCSC server
if(!(is.null(opt$cytoBand))) { cytoBand <- read.table(opt$cytoBand, header = T, sep = "\t") }

# 1) Transcript IDs to represent together
transcript_ids <-  opt$transcript_ids
replace_chars <- c("[][]", " ") # Replace the square brackets --> nextflow inputs transcript_ids tuple as a value --> "[tr_id_1, tr_id_2]"
for (char in replace_chars){ transcript_ids <-  gsub(transcript_ids, pattern = char, replacement = "") }
transcript_ids <- unlist(strsplit(transcript_ids , ","))# split by comma

# 2. Transcript Model Track 
## File paths
transcript.gc.paths <- list.files(path = gtf_path, full.names = T)
## Run
transcript.track.list <- lapply(transcript_ids, function(transcript_id){
  
  # 2.0 Transcript GTF file path
  transcript.gc.path <- grep(transcript_id, transcript.gc.paths, value = T)
  if(length(transcript.gc.path) > 1) stop(paste0("More than 1 file matching transcript ", transcript_id, " \n In ", gtf_path))
  
  ## 2.1. Read GTF
  gtf <- rtracklayer::import(transcript.gc.path)
  ## 1.2. Parse GTF
  transcript_model <- parse.gtf.v2(gtf = gtf, transcript.id = transcript_id )
  ## 1.3. Read GTF
  transcript.track <- gene.track.v2(gen.coords =  transcript_model, 
                                    which = "transcript",
                                    genome = genome_id, # this is hard-coded to human genome
                                    group.id = "transcript", # for track labeling
                                    plot.label = transcript_id ) # for track labeling
})
## Add names
names(transcript.track.list) <- transcript_ids


# 3. PFAM alignment Tracks
## File paths
pfam.gc.paths <- list.files(path = pfam_path, full.names = T)
## Run
pfam.track.list <- lapply(transcript_ids, function(transcript_id){
  
  pfam.gc.path <- grep(transcript_id, pfam.gc.paths, value = T)
  
  if(length(pfam.gc.path) > 1) stop(paste0("More than 1 file matching transcript ", transcript_id, " \n In ", pfam_path))

  ## 3.1. PFAM genomic coordinates (GC)
  pfam.gc <- read.table(file = pfam.gc.path, 
                        header = T, 
                        sep = "\t", 
                        quote = "\"")
  # TODO FUTURE: Subset by proximity to Event Coordinates! Here we do not have event coordinate (Next step)
  ## 3.2. Subset by domain partiality --> To deal with transcripts mapping to a ton of domains
  pfam.gc <- pfam.gc[order(pfam.gc[["partiality"]], decreasing = T), ]
  if(length(unique(pfam.gc[["pfam_alignment_id"]])) > 7 ) { pfam.gc <- pfam.gc[1:7, ] }

  ## 2.3. Generate alignment tracks
  pfam.track.list <- if (all(pfam.gc[["pfam_match"]]) == T ) {
    ### Parse
    pfam.gc.gviz <- parse.pfam(pfam.gc = pfam.gc)
    ## Add domain name & partiality
    pfam.gc.gviz[["domain_name"]] <- pfam.gc[["pfam_domain"]]
    pfam.gc.gviz[["partiality"]] <- pfam.gc[["partiality"]]
    ## 2.3. Generate Alignment Tracks
    al.tracks <- lapply(unique(pfam.gc.gviz[["pfam_alignment_id"]]), function(al.id) {
      
      coords.df <- pfam.gc.gviz[pfam.gc.gviz[["pfam_alignment_id"]] == al.id, ]
      # Extract partiality labels
      partiality <- pfam.gc.gviz[pfam.gc.gviz[["pfam_alignment_id"]] == al.id, "partiality"]
      # Label
      plot.label <- unique(paste0(coords.df[["domain_name"]], "(P=", round(partiality, 2), ")"))
      ## Tracks
      gene.track.v2(gen.coords = coords.df, 
                    which = "pfam", 
                    genome = genome_id, 
                    group.id = "pfam", 
                    plot.label = plot.label)
    })
    ## Unlist (function already outputs a list))
    al.tracks <- unlist(al.tracks, recursive = F)
  } else if(all(pfam_al[["pfam_match"]]) == F ){ NULL }
})

## Add names
names(pfam.track.list) <- transcript_ids

## 4. Integrate lists
transcript.pfam.track.list <- mapply(function(list.1, list.2) c(list.1, list.2), transcript.track.list, pfam.track.list)
transcript.pfam.track.list <- unlist(transcript.pfam.track.list, recursive = F)

# 5. Other tracks 
## 5.1. Genome Track
genome.track <- genome.track.fun()
## 5.2. Chromosome Track
chr <- unique(unlist(lapply(transcript.pfam.track.list, function(track) track@chromosome)))
if(length(chr) > 1 | is.na(chr)) { warning(paste0("The transcripts provided ", transcript_ids, " belong to more than one chromosome"))}
chr.track <- chr.track.fun(genome = genome_id, chr = chr, cyto.band = cytoBand)
## 5.3.  Event Track (if coordinates available)
if(is.na(genomic_start) | is.na(genomic_end)) {
	event.track <- list(NULL)
	}else{
	event.track <- event.track.fun(gen.start = genomic_start, gen.end = genomic_end, chr = chr, gen = genome_id)
}

# 6. Collect Tracks 
track.list <- c(chr.track,
                genome.track,
                event.track,
                transcript.pfam.track.list
	)
## 6.1 Rm NULL tracks 
track.list <- track.list[!unlist(lapply(track.list, is.null))]

# 7. Save
## 7.1. Save track list
saveRDS(track.list, file = out.rds.file.name)
## 7.2. Save plot PDF
plot_title <- paste0(gene_id, " - ", event_id, "\n Transcripts: ", paste(transcript_ids, collapse = " - ") )

save.plot.pdf(track.list = track.list, 
              plot.title = plot_title, 
              file.name = out.pdf.file.name)
