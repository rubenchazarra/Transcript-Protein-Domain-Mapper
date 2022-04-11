#!/usr/bin/env Rscript

# This Script Generates Visualization of Transcript Gene Model + the PFAM alignments of the corresponding Transcript
# Input file is a DF containing the PFAM mappings from each Transcript associated to the Event

## TODO: Deal with possibility that no Event Coordinates are provided

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--event_pfam_gc"),
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
    c("-a", "--genome_id"),
    action = "store",
    default = "hg38",
    type = 'character',
    help = 'UCSC Genome Browser assembly ID. Latest IDs for Human and Mouse genomes are "hg38", and "mm39" respectively'
  ),
  make_option(
    c("-t", "--gtf_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the directory where the GTF files of the transcripts to be represented are located (These are a subset of the GTF annotation file containing only the query transcript.'
  ), 
  make_option(
    c("-c", "--cyto_band"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the UCSC Browser Cyto Band table. To download cytoband files for genomes hosted at UCSC, see the UCSC Table Browser, and select Group="All Tables" and Table="cytoBand".'
  ),
  make_option(
    c("-s", "--show_non_overlapping"),
    action = "store",
    default = NA,
    type = 'logical',
    help = 'If True, PFAM domains not overlapping with AS Event coordinates will be displayed'
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

event.track.fun <- function(event_coords, chr, gen){
  # Generate Track with Coordinates of Aggregation Event
  ## Sort increasingly
  event_coords <- sort(event_coords, decreasing = F)
  ## Start
  rstarts <- event_coords[seq_along(event_coords) %%2 == 1]
  ## End
  rends <- event_coords[seq_along(event_coords) %%2 == 0]
  
  ev.track <- GeneRegionTrack(rstarts = rstarts,
                              rends = rends,
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
event_id <- opt$event_id
gene_id <- opt$gene_id
genome_id <- opt$genome_id 
show_non_overlapping <- opt$show_non_overlapping

## 0.2 Input file paths
event.pfam.gc <- opt$event.pfam.gc
gtf_path <- opt$gtf_path

# 1. Read files

## 1.1 PFAM genomic coordinates
event.pfam.gc <- readRDS(opt$event_pfam_gc)
# Extract values
#event_id <- unique(event.pfam.gc[["Event.ID"]])
transcript_ids <- unique(event.pfam.gc[["transcript_id"]])
#gene_id <- unique(event.pfam.gc[["gene.name"]])

## 1.2 Cyto Band information: from UCSC ( required to generate the chromosome Ideogram )
if(!(is.null(opt$cyto_band))) { cyto_band <- read.table(opt$cyto_band, header = T, sep = "\t") }

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
                                    genome = genome_id, # I think this could be rm from here
                                    group.id = "transcript", # for track labeling
                                    plot.label = transcript_id ) # for track labeling
})
## Add names
names(transcript.track.list) <- transcript_ids

# 3. PFAM alignment Tracks
## Split pfam genomic Coord DF by transcript
event.pfam.gc.list <- split(event.pfam.gc, event.pfam.gc[["transcript_id"]])

## Filter Domains that do not overlap
event.pfam.gc.list <-  if(show_non_overlapping){
  event.pfam.gc.list
  } else {
  # Filter to only overlapping
  lapply(event.pfam.gc.list, function(df){
  df <- df[df[["overlap.domain.event"]] == T, ]
  df })
}

## Run
pfam.track.list <- lapply(event.pfam.gc.list, function(pfam.gc){
  # Not empty alignments
  if(!is.null(pfam.gc) | all(is.na(pfam.gc[["overlap.domain.event"]]))){
    
  if(length(unique(pfam.gc[["pfam_alignment_id"]])) > 10 ) {
    uniq.alignments <- unique( pfam.gc[ , c("pfam_alignment_id", "partiality")] )
    ord.uniq.alignment.ids <- uniq.alignments[order(uniq.alignments[["partiality"]], decreasing = T), "pfam_alignment_id"]
    ord.uniq.alignment.ids <- ord.uniq.alignment.ids[1:10, ]
    pfam.gc <- pfam.gc[ pfam.gc[["pfam_alignment_id"]] %in% ord.uniq.alignment.ids, ]
  }
  
  ## 2.3. Generate alignment tracks
  pfam.track.list <- if (all(!is.na(pfam.gc[["pfam_domain"]]))) {
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
  }} else { NULL }
})

## Add names <-- Names already there
#names(pfam.track.list) <- transcript_ids

## 4. Integrate lists
transcript.pfam.track.list <- mapply(function(list.1, list.2) c(list.1, list.2), transcript.track.list, pfam.track.list, SIMPLIFY = F)
transcript.pfam.track.list <- unlist(transcript.pfam.track.list, recursive = F)

# 5. Other tracks 
## 5.1. Genome Track
genome.track <- genome.track.fun()
## 5.2. Chromosome Track
chr <- unique(unlist(lapply(transcript.pfam.track.list, function(track) track@chromosome)))
if(length(chr) > 1 | is.na(chr)) { warning(paste0("The transcripts provided ", transcript_ids, " belong to more than one chromosome"))}
chr.track <- chr.track.fun(genome = genome_id, chr = chr, cyto.band = cyto_band)

## 5.3.  Event Track (if coordinates available)
event_coords <- as.character(unique(event.pfam.gc[["event.coords"]]))
event_coords <- as.numeric(unlist(strsplit(event_coords, split = "; ")))
print(event_coords)
## Event Track
if(any(is.na(event_coords))){
  event_coords <- event_coords[!is.na(event_coords)]
  print(event_coords)
  event.track <- event.track.fun(event_coords = event_coords, chr = chr, gen = genome_id)
} else if( all(is.na(event_coords))){
  event.track <- NULL
} else {
  event.track <- event.track.fun(event_coords = event_coords, chr = chr, gen = genome_id)
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
file.name.rds <- paste(gene_id, event_id, paste0(unlist(transcript_ids), collapse = "_"), "Visualisation.rds", sep = "_")
saveRDS(track.list, file = file.name.rds)

## 7.2. Save plot PDF
plot_title <- paste0(gene_id, " - ", event_id, "\n Transcripts: ", paste(transcript_ids, collapse = " - ") )
file.name.pdf <- paste(gene_id, event_id, paste0(unlist(transcript_ids), collapse = "_"), "Visualisation.pdf", sep = "_")

save.plot.pdf(track.list = track.list, 
              plot.title = plot_title, 
              file.name = file.name.pdf)