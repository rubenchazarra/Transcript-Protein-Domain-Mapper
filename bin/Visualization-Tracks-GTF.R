#!/usr/bin/env Rscript

# Generate Visualizations of the Transcript Exon Model + PFAM alignment output with the Bioconductor GViz package
# TODO: Extract transcript genomic coordinates from GTF and not from EnsmblDB for consistency

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
    help = 'Path to the RDS object storing lists of tracks for vizualization.'
  ), 
  make_option(
    c("-v", "--viz_track_plot"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the PDF file with the track plots.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(Gviz))
suppressPackageStartupMessages(require(rtracklayer))

## Functions ## 

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

generate_transcript_track <- function(transcriptModel, genome, transcript_id, UCSC_cytoBand){
  chr <- unique(as.character(transcriptModel$seqnames))
  # Loop across 'transcript_ids' we have more than one in PFAM alignments
  atrack_list <- list()
  for(id in transcript_id){ # Note: this is looped because later we are using the same fubction for the PFAM alignments (where there may be more than one) 
    sub_transcriptModel <- transcriptModel[transcriptModel$symbol == id, ]
    track_name <- unique(sub_transcriptModel$symbol)
    atrack_list[[id]] <- GeneRegionTrack(sub_transcriptModel, genome = genome, chr = chr, name = track_name, cex.axis = 15)
  }
   
  gtrack <- GenomeAxisTrack()
  itrack <- IdeogramTrack(genome = genome, chromosome = chr, bands = UCSC_cytoBand )
  
  # merge track lists
  track_list <- c(list("IdeogramTrack" = itrack, 
                     "GenomeAxisTrack" = gtrack),
                  atrack_list)

  return(track_list)
}

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
             "PFAM.alignment.ID" = paste(pfam.genomic.coord$PFAM.alignment.ID, direction)
             )
  # convert alignment ID to character  Required for visualization
  df[["symbol"]] = as.character(df[["PFAM.alignment.ID"]])
  return(df)
}

merge_track_lists <- function(transcript_track_list, pfam_track_list){
  ## Merge 'transcript_track_list' + 'pfam_track_list' avoiding duplication of chromosome and genome tracks
  merged_list <- c(transcript_track_list, # Ideogram, GenomeAxis and TranscriptTrack
                   pfam_track_list[3:length(pfam_track_list)] # PFAM Alignment Tracks without Ideogram and GenomeAxis
  ) 
  return(merged_list)
}

# Plot track list
plot_track_list <- function(merged_track_list, out_pdf_path){
  ## Save the Track Plot to PDF file in 'out_pdf_path'
  
  pdf(out_pdf_path)
  Gviz::plotTracks(merged_track_list)
  dev.off()
}


# 0. Params
transcript_id <- opt$transcript_id
# 1. Read GTF 
gtf <- rtracklayer::import(opt$gtf)
# 2. Read cytoband information from UCSC, required to generate the chromosome Ideogram without connection to the UCSC server
cytoBand <- read.table(opt$cytoBand, header = T, sep = "\t") 
# 3. Generate transcript GTF
transcript_model <- gtf_processing(gtf = gtf, transcript_id = transcript_id)
# 4. Generate Transcript Track
transcript_track_list <-  generate_transcript_track(transcriptModel = transcript_model, genome = "hg38", transcript_id =  unique(transcript_model$symbol), UCSC_cytoBand = cytoBand)

# 4. Read PFAM alignment genomic coordinates 
# TODO: Find a way to map protein coordinates to genomic coordinates with GTF
pfam.genomic.coord <- read.table(opt$pfam_genomic_coord, header = T, row.names = 1)
# 5. Coere PFAM df to required visualization format
pfam_model <- coerce.pfam.df.to.gviz.format(pfam.genomic.coord)

# If NOT EMPTY PFAM alignment: 
if(!(all(is.na(pfam_model[["start"]])))){
  # 6. Generate PFAM TrackList
  pfam_track_list <- generate_transcript_track(transcriptModel = pfam_model, genome = "hg38", transcript_id = unique(pfam_model$symbol), UCSC_cytoBand = cytoBand)
  # 7. Merge Transcript and PFAM Tracks
  merged_track_list <- merge_track_lists(transcript_track_list, pfam_track_list)
  # 8. Save TrackList file
  saveRDS(merged_track_list, file = opt$viz_track_list)
  # 9. Save Plot
  plot_track_list(merged_track_list, opt$viz_track_plot)
  # If EMPTY PFAM alignment
} else {
  # 8. Save TrackList file
  saveRDS(transcript_track_list, file = opt$viz_track_list)
  # 9. Save Plot
  plot_track_list(transcript_track_list,  opt$viz_track_plot)
}
