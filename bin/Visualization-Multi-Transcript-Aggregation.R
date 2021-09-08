#!/usr/bin/env Rscript

# This Script Generates Visualization of Transcript Gene Model + the PFAM alignments of the corresponding Transcript
# The Goal is to represent the Alternative Splicing Event. The transcripts involved in each AS event are defined in an external Aggregation CSV
# TODO: Extract transcript genomic coordinates from GTF and not from EnsmblDB for consistency

suppressPackageStartupMessages(require(optparse))

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
    c("-p", "--pfam_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the directory where the PFAM alignments of the differents Transcripts to represent are located.'
  ), 
  make_option(
    c("-g", "--gtf_path"),
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
    c("-o", "--viz_track_list"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the output RDS object storing lists of GViz Tracks for vizualization.'
  ), 
  make_option(
    c("-v", "--viz_track_plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the PDF file with the Track Plots.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Gviz))
suppressPackageStartupMessages(require(rtracklayer))

## Functions ## 

# 1. Create list of files
create.file.path.list <- function(file.path, file.patt, transcript_ids){
  # Create list of files where element name is the pertinent transcript ID
  file.vec <- list.files(path = file.path, pattern = file.patt, full.names = T)
  list.order <- lapply(transcript_ids, function(tr) grep(tr, file.vec))
  # Order list to match: TranscriptID - file
  path.list <- as.list(file.vec[unlist(list.order)])
  names(path.list) <- transcript_ids  # TODO: PROBLEM HERE
  return(path.list)
}

# 2. Extract chromosome number
extract.chr.num <- function(gtf.list){
  ## Extract chromosome number from the GTF list of the different transcripts
  chr.vec <- unique(as.character(unlist(lapply(gtf.list, function(gtf) unique(gtf$seqnames)))))
  if( length(chr.vec) == 1 ) { # There can only be one chr_n value among the transcripts
    return(chr.vec)
  } else {
    print (chr.vec)
    stop ("There can't be more than 1 chr number. \n 
          CHECK! Some transcripts may not be from the same gene.")
  }
}

# 3. GTF processing for GViz Visualization
gtf_processing <- function(gtf, transcript_id){
  ## Subset GTF and perform pertinent modifications for visualization
  sub_gtf <- gtf[grep(transcript_id, gtf$transcript_id), ]
  # Create dataframe --
  transcriptModel <- as.data.frame(sub_gtf)
  transcriptModel <- transcriptModel[transcriptModel$type =="exon",]
  # Rename columns --
  transcriptModel$gene <- transcriptModel$gene_name
  transcriptModel$exon <- transcriptModel$exon_id
  transcriptModel$transcript <- transcriptModel$transcript_id
  transcriptModel$gene <- transcriptModel$gene_id
  direction <- ifelse(transcriptModel$strand == "+", ">", "<")
  transcriptModel$symbol <- paste(transcriptModel$transcript_id, direction) # label in plot
  
  return(transcriptModel)
}

# 4. Check Gene IDs are consistent across the different transcripts of the AS Event
check.unique.gene.ids <- function(gtf.list){
  # Check that gene IDs across the different Transcript_IDs are the same. All transcripts must come from the same Gene_id
  gene_ids <- unlist(lapply(gtf.list, function(gtf) unique(gtf$gene_id)))
  transcript_ids <- unlist(lapply(gtf.list, function(gtf) unique(gtf$transcript_id)))
  
  if(length(unique(gene_ids)) > 1){
    stop(
      paste0("Genes IDs are: ", paste(gene_ids, collapse = ", "), 
             "\nFor Transcript IDs: ", paste(transcript_ids, collapse = ", "))
    )
  }
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

# 6. Coerce PFAM df to GViz format
coerce.pfam.df.to.gviz.format <- function(pfam.genomic.coord){
  ## Coerce PFAM alignment output df to format readable by GViz package.
  
  direction <- ifelse(pfam.genomic.coord$strand == "+", ">", "<")
  # Basically, create a df [almost] equivalent to the 'transcript_model' 
  df = data.frame("seqnames" = pfam.genomic.coord$seqnames, 
                  "start" = pfam.genomic.coord$start,
                  "end" = pfam.genomic.coord$end, 
                  "width" = pfam.genomic.coord$width, 
                  "strand" = pfam.genomic.coord$strand, 
                  "transcript" = pfam.genomic.coord$tx_id,
                  "exon" = pfam.genomic.coord$exon_id,
                  "protein" = pfam.genomic.coord$protein_id,
                  "symbol" = paste(pfam.genomic.coord$transcript, direction),
                  "PFAM.Alignment.ID" = paste(pfam.genomic.coord$PFAM.Alignment.ID, direction), 
                  "PFAM.Domain" = pfam.genomic.coord$PFAM.Domain, 
                  "PFAM.Domain.Description" = pfam.genomic.coord$PFAM.Domain.Description
  )
  # Convert alignment ID to character  Required for visualization
  df[["symbol"]] = as.character(df[["PFAM.Alignment.ID"]])
  return(df)
}

## 7. Integrate Tracks into a Single Alternating List
integrate.tracks <- function(itrack, gtrack, transcript.track.list, pfam.track.list){
  ## Integrate Tracks into a Single Alternating List. As: Transcript_1, PFAM_A_T1, Transcript_2, PFAM_A_T2, PFAM_B_T2 
  keys <- unique(c(names(transcript.track.list), names(pfam.track.list)))
  ordered.track.list <- setNames(mapply(c, transcript.track.list[keys], pfam.track.list[keys]), keys)
  # Add Ideogram and Genome Track
  track.list <- c(itrack, 
                  gtrack, 
                  unlist(ordered.track.list, recursive = T))
  return(track.list)
}

## 8. Plot track list
plot.tracks <- function(track.list, event.id, out.pdf.path){
  ## Save the Track Plot to PDF file in 'out_pdf_path'
  pdf(out.pdf.path)
  Gviz::plotTracks(main = event.id, track.list)
  dev.off()
}

## Execute

# 1) Transcript IDs of the AS event to represent
transcript_ids <-  gsub(opt$transcript_ids, pattern = " ", replacement = "") # remove blank spaces 
transcript_ids <- unlist(strsplit(transcript_ids , ","))# split by comma

# 2) Read PFAM File List
## Read 
pfam.list <- create.file.path.list(file.path = opt$pfam_path, 
                                   file.patt = ".txt", 
                                   transcript_ids = transcript_ids)
pfam.list <- lapply(pfam.list, function(path) read.table(path, header = T))

## 3) Read and process transcript-GTF File List
gtf.list <- create.file.path.list(file.path = opt$gtf_path,
                                  file.patt = ".gtf",
                                  transcript_ids = transcript_ids)
gtf.list <- lapply(gtf.list, function(path) data.frame(rtracklayer::import(path)))

# 4) Process Files
## PFAM File List
### Coerce PFAM Df to GViz Format
pfam.gviz.list <- lapply(pfam.list, function(pfam_df)  coerce.pfam.df.to.gviz.format(pfam.genomic.coord = pfam_df))

## GTF File List
gtf.processed.list <- lapply(names(gtf.list), function(tr_name) gtf_processing(gtf = gtf.list[[tr_name]], transcript_id = tr_name))
### Add names 
names(gtf.processed.list) <- names(gtf.list)
## Check the transcripts come from unique Gene ID
check.unique.gene.ids(gtf.list = gtf.processed.list)

# 5) Generate Tracks 

# 5.1) Genome Track
## Genome version
gen <- "hg38" # TODO: <-- This value is hard-coded. Extract from somewhere.
## Genome Track
gtrack <- GenomeAxisTrack()

# 5.2) Chromosome Track
## Chromosome Number
chr <- extract.chr.num(gtf.list = gtf.list)
## Ideogram (Chromosome) Track
cytoBand <- read.table(opt$cytoBand, header = T, sep = "\t") # 2. Read cytoband information from UCSC, required to generate the chromosome Ideogram without connection to the UCSC server
itrack <- IdeogramTrack(genome = gen, chromosome = chr, bands = cytoBand)



# 5.3) Transcript Tracks
transcript.track.list <- lapply(names(gtf.processed.list), function(tr_name) generate.transcript.track(coord_df = gtf.processed.list[[tr_name]],
                                                                                                       id_name = "transcript_id",
                                                                                                       gen = gen, chr = chr,
                                                                                                       color = "palegreen3"))
## Add names
names(transcript.track.list) <- names(gtf.processed.list)

# 5.4) PFAM Tracks
pfam.track.list <- lapply(pfam.gviz.list, function(gviz_df) generate.transcript.track(coord_df = gviz_df, 
                                                                                      gen = gen, 
                                                                                      chr = chr, 
                                                                                      id_name = "PFAM.Domain", 
                                                                                      color = "cadetblue1"))

# 6) Integrate Tracks
##  Merge Transcript + PFAM Tracks in an alternate order.
track.list <- integrate.tracks(itrack = itrack, 
                 gtrack = gtrack, 
                 transcript.track.list = transcript.track.list, 
                 pfam.track.list = pfam.track.list)

# 7) Save Complete Track List
saveRDS(track.list, file = opt$viz_track_list)

# 8) Generate and Save Track Plot
plot.tracks(track.list = track.list, event.id = opt$event_id, out.pdf.path = opt$viz_track_plot)
