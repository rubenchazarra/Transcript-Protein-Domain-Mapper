#!/usr/bin/env Rscript

# Generate Visualizations of the Transcript Exon Model + PFAM alignment output with the Bioconductor GViz package

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

## Functions ## 

# Extract genomic coordinates from transcript ID with biomaRt
get.exonic.coords <- function(transcript_id, mart){
  ## Function to retrieve exon genomic coordinates from transcript ID.
  # Get ensembl mart 
  suppressPackageStartupMessages(require(biomaRt))
  ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl") # <-- This is hard-coded
  # Get 'transcript_id' coordinates
  getBM(
    attributes=c('ensembl_transcript_id', 'chromosome_name', 'exon_chrom_start', 'exon_chrom_end', 'strand', 'ensembl_exon_id', 'ensembl_peptide_id',  'ensembl_gene_id', 'description' ),  
    filters = 'ensembl_transcript_id', 
    values = transcript_id, 
    mart = ensembl, 
    useCache = F)
}

# Coerce biomRt df to GViz readable format
coerce.biomart.to.gviz.format <- function(biomart.df){
  ## Coerce dataframe of transcript exonic coordinates (from bioMart) to format readable by GViz package
  
  # 1. Create 'chromosome' column (to 'chrX' format)
  biomart.df[["chromosome"]] = paste0("chr", biomart.df[["chromosome_name"]])
  # 2. Rename coordinate columns
  colnames(biomart.df)[colnames(biomart.df) == "exon_chrom_start"] = "start"
  colnames(biomart.df)[colnames(biomart.df) == "exon_chrom_end"] = "end"
  # 3. Add 'width' colum (end-start+1)
  biomart.df[["width"]] = as.numeric(biomart.df[["end"]] - biomart.df[["start"]] + 1 )
  # 4. Add exon numbers
  #biomart.df[["exon_rank"]] = c(1:nrow(biomart.df)) [WOULD THIS BE CORRECT?]
  return(biomart.df)
}

# Coerce PFAM alingment output df to GViz readable format
coerce.pfam.df.to.gviz.format <- function(pfam.df){
  ## Coerce PFAM alignment output df to format readable by GViz package.
  ## This requires less changes since the PFAM df comes from a IRanges object
  
  # 1. Create 'seqnames' (depicting chromosome name) column (tp 'chrX' format)
  pfam.df[["chromosome"]] = paste0("chr", pfam.df[["seqnames"]])
  return(pfam.df)
}


# 3. Generate track list
generate.track.list <- function(viz.df, viz_name){
  ## Function to generate a Track list for Visualization from 'viz.df'
  suppressPackageStartupMessages(require(Gviz))
  # 0. Params
  ## Chr 
  chr <- as.character(unique(viz.df[["chromosome"]]))
  ## Gen
  gen <- c(chr="hg19") # <-- Careful, this is hardcoded
  
  # 1. Genome Track 
  gtrack <- GenomeAxisTrack()
  # 2. Chromosome track
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  # 3. Annotation track
  if(viz_name == "PFAM"){ 
    group_exons = viz.df[["tx_id"]]
    exon_ids = viz.df[["exon_id"]]
  } else {
    group_exons = viz.df[["ensembl_transcript_id"]]
    exon_ids = viz.df[["ensembl_exon_id"]]
  }
  
  atrack <- AnnotationTrack(viz.df, name = viz_name,  
                            group = group_exons, 
                            showId = T)
  
  # 4. Generate Region Gene Track
  grtrack <- GeneRegionTrack(viz.df, genome = gen,
                             chromosome = chr, name = viz_name, 
                             symbol  = exon_ids, 
                             #group = group_exons, 
                             showId = T
  ) # maybe change this!
  
  
  track.list <- list("Chr.ideogram" = itrack, # Chr ideogram
                     "Genome.Track" = gtrack, # Genome Annotation, 
                     "Annot.Track" = atrack, # Annotation Track
                     "Gene.Region.Track"= grtrack # Gene Region Track [Actual PFAM alignment]
  )
  
  return(track.list)
}

# Merge track lists from Transcript + PFAM alignment
merge.track.lists <- function(transcript.track.list, pfam.track.list){
  ## Merge track lists coming from: i) Transcript and ii) PFAM alignment
  merged.track.list <- list("Chr.ideogram" = transcript.track.list$Chr.ideogram, 
                            "Genome.Track" = transcript.track.list$Genome.Track, 
                            # Transcript
                            "Transcript.Annot.Track" = transcript.track.list$Annot.Track, 
                            "Transcript.Gene.Region.Track" = transcript.track.list$Gene.Region.Track,
                            # PFAM 
                            "PFAM.Annot.Track" = pfam.track.list$Annot.Track,
                            "PFAM.Gene.Region.Track" = pfam.track.list$Gene.Region.Track)
  return(merged.track.list)
}

# Plot track list
plot.track.list <- function(merged.track.list, out_pdf_path){
  ## Save the Track Plot to PDF file in 'out_pdf_path'
  
  pdf(out_pdf_path)
  Gviz::plotTracks(merged.track.list)
  dev.off()
}


## Run 

transcript_id = opt$transcript_id
# 0. Read PFAM alignment output
pfam.df <- read.table(opt$pfam_genomic_coord, header = T, row.names = 1)
# 1. Coerce PFAM out df to GViz readable format
pfam.df <- coerce.pfam.df.to.gviz.format(pfam.df = pfam.df)
# 2. Get exonic coordinates from 'transcript_id' 
transcript.df <- get.exonic.coords(transcript_id = transcript_id , mart = ensembl)


## Condition A: Existing transcript_ID exonic coords & Existing PFAM output: 
if(nrow(transcript.df) > 0 && pfam.df[["chromosome"]] != "chrNA"){
  
  # 3. Coerce transcript exonic coordinates to GViz readable format
  transcript.viz.df <- coerce.biomart.to.gviz.format(biomart.df = transcript.df)
  # 4. Generate Transcript Track List 
  transcript.track.list <- generate.track.list(viz.df = transcript.viz.df, viz_name = "Transcript")
  ## 5. Generate PFAM Track List 
  pfam.track.list <-  generate.track.list(viz.df = pfam.df, viz_name = "PFAM")
  ## 6. Merge Trancript + PFAM Track Lists
  merged.track.list <- merge.track.lists(transcript.track.list = transcript.track.list, pfam.track.list = pfam.track.list)
  ## 7. Save Merged Track List (RDS)
  saveRDS(merged.track.list, file = opt$viz_track_list)
  ## 8. Save  Merged Track List Plot (PDF)
  track.plot <- plot.track.list(merged.track.list = merged.track.list, out_pdf_path = opt$viz_track_plot)
  
  ## Condition B: Existing transcript ID exonic coord & EMPTY PFAM output:
}else if (nrow(transcript.df) > 0 && pfam.df[["chromosome"]] == "chrNA"){
  
  # 3. Coerce transcript exonic coordinates to GViz readable format
  transcript.viz.df <- coerce.biomart.to.gviz.format(biomart.df = transcript.df)
  # 4. Generate Transcript Track List 
  transcript.track.list <- generate.track.list(viz.df = transcript.viz.df, viz_name = "Transcript")
  ## 5. Generate PFAM Track List --> NO
  ## 6. Merge Trancript + PFAM Track Lists --> NO
  ## 7. Save TRANSCRIPT Track List (RDS)
  saveRDS(transcript.track.list, file = opt$viz_track_list)
  ## 8. Save TRANSCRIPT Track List Plot (PDF)
  track.plot <- plot.track.list(merged.track.list = transcript.track.list, out_pdf_path = opt$viz_track_plot)
  
  
  ## Condition C: EMPTY transcript ID exonic coord & Existing PFAM output:
} else if (nrow(transcript.df) == 0 && pfam.df[["chromosome"]] != "chrNA"){
  
  # 3. Coerce transcript exonic coordinates to GViz readable format --> NO
  # 4. Generate Transcript Track List --> NO
  ## 5. Generate PFAM Track List 
  pfam.track.list <-  generate.track.list(viz.df = pfam.df, viz_name = "PFAM")
  ## 6. Merge Trancript + PFAM Track Lists --> NO
  ## 7. Save PFAM Track List (RDS)
  saveRDS(pfam.track.list, file = opt$viz_track_list)
  ## 8. Save PFAM Track List Plot (PDF)
  track.plot <- plot.track.list(merged.track.list = pfam.track.list, out_pdf_path = opt$viz_track_plot)
  
  ## Condition D: EMPTY transcript ID exonic cood & EMPTY PFAM output
} else if(nrow(transcript.df) == 0 && pfam.df[["chromosome"]] == "chrNA") {
  
  # 3. Coerce transcript exonic coordinates to GViz readable format --> NO
  # 4. Generate Transcript Track List --> NO
  ## 5. Generate PFAM Track List --> NO
  ## 6. Merge Trancript + PFAM Track Lists --> NO
  ## 7. Save 'NULL' Track List (RDS)
  saveRDS(list(NULL), file = opt$viz_track_list)
  ## 8. Save EMPTY PDF Plot
  pdf(opt$viz_track_plot)
  print("Nothing")
  dev.off
}
