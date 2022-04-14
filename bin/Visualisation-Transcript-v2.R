#!/usr/bin/env Rscript

# Generate Visualizations of the Transcript Exon Model + PFAM alignment output with the Bioconductor GViz package
# TODO: Extract transcript genomic coordinates from GTF and not from EnsmblDB for consistency
## Vulnerabilities: This script works with Gencode GTF annotation. With the Ensembl GTF it will break

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-t", "--transcript_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Ensembl Transcript ID.'
  ), 
  make_option(
    c("-p", "--pfam_genomic_coord"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to PFAM alignment table of the transcript_id query.'
  ), 
  make_option(
    c("-g", "--gtf"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to GTF annotation file.'
  ), 
  make_option(
    c("-a", "--genome_id"),
    action = "store",
    default = "hg38",
    type = 'character',
    help = 'UCSC Genome Browser assembly ID. Latest IDs for Human and Mouse genomes are "hg38", and "mm39" respectively'
  ), 
  make_option(
    c("-c", "--cytoBand"),
    action = "store",
    default = NULL,
    type = 'character',
    help = 'Path to the UCSC Browser Cyto Band table. To download cytoband files for genomes hosted at UCSC, see the UCSC Table Browser, and select Group="All Tables" and Table="cytoBand".'
  ), 
  make_option(
    c("-o", "--vis_track_list"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the RDS object storing lists of tracks for vizualization.'
  ), 
  make_option(
    c("-v", "--vis_track_plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the PDF file with the track plots.'
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

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

chr.track.fun <- function(genome, chr, cyto.band ){
  suppressPackageStartupMessages(require(Gviz))
  ## Generate chromosome track ## TODO: what if cyto.band is NULL (is not available ?)
  itrack <- IdeogramTrack(genome = genome,
                          chromosome = chr,
                          bands = cyto.band,
                          name = "chr.track")
  list("Chr.track" = itrack)
}

genome.track.fun <- function(){
  suppressPackageStartupMessages(require(Gviz))
  # Generate genome track
  gtrack <- GenomeAxisTrack(name = "genome.track")
  list("Genome.Track" = gtrack)
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
                          name = "chr.track")
  list("Chr.track" = itrack)
}

genome.track.fun <- function(){
  suppressPackageStartupMessages(require(Gviz))
  # Generate genome track
  gtrack <- GenomeAxisTrack(name = "genome.track")
  list("Genome.Track" = gtrack)
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

# 0. Params
## 0.1. 
transcript_id <- opt$transcript_id
genome_id <- opt$genome_id
## 0.2. Output file names
out.rds.file.name <- opt$vis_track_list
out.pdf.file.name <- opt$vis_track_plot
pfam.gen.coords <- opt$pfam_genomic_coord
## 0.3. Read Cytoband information from UCSC, required to generate the chromosome Ideogram without connection to the UCSC server
if(!(is.null(opt$cytoBand))) { cytoBand <- read.table(opt$cytoBand, header = T, sep = "\t") }

# 1. Transcript Model Track 
## 1.1. Read GTF
gtf <- rtracklayer::import(opt$gtf)
## 1.2. Parse GTF
transcript_model <- parse.gtf.v2(gtf = gtf, transcript.id = transcript_id )
## 1.3. Read GTF
transcript.track <- gene.track.v2(gen.coords =  transcript_model, 
                                       which = "transcript",
                                       genome = genome_id, # this is hard-coded to human genome
                                       group.id = "transcript", # for track labeling
                                       plot.label = transcript_id ) # for track labeling


# 2. PFAM alignment Tracks
## 2.1. PFAM genomic coordinates (GC)
pfam.gc <- read.table(file = pfam.gen.coords, header = T, sep = "\t",  quote = "\"")
## 2.2. Subset by domain partiality --> To deal with transcripts mapping to a ton of domains
if(length(unique(pfam.gc[["pfam_alignment_id"]])) > 10 ) {
    uniq.alignments <- unique( pfam.gc[ , c("pfam_alignment_id", "partiality")] )
    ord.uniq.alignment.ids <- uniq.alignments[order(uniq.alignments[["partiality"]], decreasing = T), "pfam_alignment_id"]
    ord.uniq.alignment.ids <- ord.uniq.alignment.ids[1:10]
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
    gene.track.v2(gen.coords = coords.df, which = "pfam", genome = genome_id, group.id = "pfam", plot.label = plot.label)
  })
  ## Unlist (function already outputs a list))
  al.tracks <- unlist(al.tracks, recursive = F)
  } else { NULL }

# 3. Other tracks
## 3.1. Genome Track
genome.track <- genome.track.fun()
## 3.2. Chromosome Track
chr <- as.character(unique(transcript_model$seqnames)[!is.na(unique(transcript_model$seqnames))]) # skip NAs if there are
chr.track <- chr.track.fun(genome = genome_id, chr = chr, cyto.band = cytoBand)

# 4. Collect Tracks 
track.list <- c(chr.track,
                genome.track,
                transcript.track,
                pfam.track.list
                )

# 3. Save
## 3.1. Save track list
saveRDS(track.list, file = out.rds.file.name)
## 3.2. Save plot PDF
save.plot.pdf(track.list = track.list, 
              plot.title = transcript_id, 
              file.name = out.pdf.file.name)
