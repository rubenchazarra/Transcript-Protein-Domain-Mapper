#!/usr/bin/env Rscript

## Script to map PFAM domain alignment coordinates to genomic coordinates.
## Note: edbx object is hardcoded to Human database

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-a", "--pfam_alignment"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'PFAM alignment file'
  ), 
  make_option(
    c("-t", "--transcript_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Ensembl Transcript ID'
  ), 
  make_option(
    c("-g", "--gtf"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to GTF annotation file.'
  ), 
  make_option(
    c("-o", "--output_coord"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output file with genomic coordinates of the PFAM domain hits to the transcript protein sequence.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

suppressPackageStartupMessages(require(ensembldb))
suppressPackageStartupMessages(require(EnsDb.Hsapiens.v86))
suppressPackageStartupMessages(require(rtracklayer))

## Functions ##
create.iranges.pfam <- function(pfam_alignment, protein_id){
  ## Create IRanges object from PFAM alignment data
  pfam.al.start <- as.numeric(pfam_alignment[["ali_from"]])
  pfam.al.stop <- as.numeric(pfam_alignment[["ali_to"]])
  n.alignments <- length(pfam.al.start)
  protein_names <- rep(protein_id, times = n.alignments)
  # Create IRanges object
  ir.prot <- IRanges::IRanges(start = pfam.al.start, end = pfam.al.stop, names = protein_names)
  return(ir.prot)
}

extract.prot.genomic.coord.iranges <- function(iranges_prot, edbx){
  ## Extract genomic coordinates from PFAM alignment IRanges object
  prot.genome.coord <- ensembldb::proteinToGenome(x = iranges_prot, db = edbx)
  return(prot.genome.coord)
}

create.empty.iranges.pfam <- function(transcript_id, chrN){
  ## Create empty Iranges-like data frame for empty PFAM alignment outputs
  table.names <- c("pfam_alignment_id", 
                   "pfam_domain", 
                   "pfam_domain_description", 
                   
                   "seqnames", 
                   "start",
                   "end", 
                   "width", 
                   "strand", 
                   "protein_id", 
                   "transcript_id", 
                   "exon_id", 
                   "exon_number",  
                   "cds_ok",  
                   "protein_start",  
                   "protein_end" )
  
  empty.align <- rep(c(NA), times = length(table.names))
  names(empty.align) <- table.names
  
  # Fill in the fields we can ('PFAM.alignment.ID', 'Transcript.ID', 'seqnames', 'tx_id')
  empty.align[["pfam_alignment_id"]] = paste0(transcript_id, "-NA", "-from-NA-to-NA-", chrN)
  empty.align[["seqnames"]] = chrN
  empty.align[["transcript_id"]] = transcript_id
  
  empty.pfam.df = data.frame(t(empty.align))
  return(empty.pfam.df)
}

# Function that runs the 2 previous ones
extract.genomic.coord.pfam.alignment <- function(pfam_alignment, protein_id, edbx){
  ## Extract genomic coord of single or multiple PFAM alignments
  ## Also Add PFAM alignment ID 
  
  # 1. Create PFAM alignment identifier string
  pfam.alignment.id <- paste0(pfam_alignment[["query_name"]], "-",  pfam_alignment[["domain_name"]], "-from-", as.numeric(pfam_alignment[["ali_from"]]),
                       "-to-", as.numeric(pfam_alignment[["ali_to"]]), "-", chrN ) 
  # 2. Create IRanges object
  ir.prot <- create.iranges.pfam(pfam_alignment, protein_id)
  # 3. Extract genomic coordinates
  prot.genome.coord.list <- extract.prot.genomic.coord.iranges(iranges_prot = ir.prot, edbx = edbx)
  prot.genome.coord.df <- data.frame(unlist(prot.genome.coord.list, recursive = F, use.names = F))
  ## 4. Add alignment and transcript ID
  n.exons <- nrow(prot.genome.coord.df)
  
  pfam.alignment.id.df <- data.frame(
    "pfam_alignment_id" = rep(pfam.alignment.id, times = n.exons), 
    "pfam_domain" = rep(pfam_alignment[["domain_name"]], times = n.exons), 
    "pfam_domain_description" = rep(pfam_alignment[["description"]], times = n.exons)
    )
  # 5. Merge 2 dfs
  cbind(pfam.alignment.id.df, prot.genome.coord.df)
}

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

# 0. Params
transcript_id <- opt$transcript_id

# 1. Read GTF 
gtf <- rtracklayer::import(opt$gtf)
# 2. Generate transcript GTF
transcript_model <- gtf_processing(gtf = gtf, transcript_id = transcript_id)
## Protein ID
protein_id <- unique(transcript_model[["protein_id"]])
## Strip gtf id version from protein_id
protein_id <- gsub(pattern = "\\..*", replacement = "",  x = protein_id)
## Chromosome name
chrN <- as.character(unique(transcript_model[["seqnames"]]))

# 2. Read files
pfam_alignment <- read.table(opt$pfam_alignment, header = T, sep = "\t", quote = "\"")

# Condition: if PFAM alignment contains some output
if(all(pfam_alignment[["pfam_match"]] == TRUE)){
  # NOTE: This is hard coded. Restricted to the Human ENS DB
  edbx <- EnsDb.Hsapiens.v86
  # 3. Extract genomic coordinates from PFAM alignments
  genomic.coord.list <- apply(pfam_alignment, 1, function(x) 
    extract.genomic.coord.pfam.alignment(pfam_alignment = x, protein_id = protein_id, edbx = edbx))
  ## Merge
  genomic.coord.df <- Reduce(rbind, genomic.coord.list)
  # If No CDS found for: 'protein_id' --> DO SOMETHING HERE -- I am not sure this is a problem: if no CDS are found, the cds_ok values are set to FALSE, but no additional proble is seen
  #if(nrow(genomic.coord.df) == 0 ){
  #  genomic.coord.df <- create.empty.iranges.pfam(transcript_id = transcript_id)
  #}
  } else if(all(pfam_alignment[["pfam_match"]] == FALSE)){
  genomic.coord.df <- create.empty.iranges.pfam(transcript_id = transcript_id, chrN)
}

# 6. Save 
write.table(genomic.coord.df, file = opt$output_coord, sep = "\t", quote = F)
