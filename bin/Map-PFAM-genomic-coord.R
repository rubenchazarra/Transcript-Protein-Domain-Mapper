#!/usr/bin/env Rscript

## Script to map PFAM domain alignment coordinates to genomic coordinates.

suppressPackageStartupMessages(require(optparse))
## TODO --> Implement the condition where the PFAM alignment is empty. We still want to get the genomic coordinates of the protein, but we need the length of the protein sequence.
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
  pfam.al.start <- as.numeric(pfam_alignment[["alignment.from"]])
  pfam.al.stop <- as.numeric(pfam_alignment[["alignment.to"]])
  # Create IRanges object
  ir.prot <- IRanges(start = pfam.al.start, end = pfam.al.stop, names = protein_id)
  return(ir.prot)
}

extract.prot.genomic.coord.iranges <- function(iranges_prot, edbx){
  ## Extract genomic coordinates from PFAM alignment IRanges object
  prot.genome.coord <- ensembldb::proteinToGenome(x = iranges_prot, db = edbx)
  return(prot.genome.coord)
}

# Function that runs the 2 previous ones
extract.genomic.coord.pfam.alignment <- function(pfam_alignment, protein_id, edbx){
  ## Extract genomic coord of single or multiple pfam alignments
  ## Also Add PFAM alignment ID 
  
  # 1. Create PFAM alignment identifier string
  pfam.alignment.id <- paste0(pfam_alignment[["transcript.id"]], "-", pfam_alignment[["pfam.name"]], "-domain-", pfam_alignment[["n.domain"]], "-from-", pfam_alignment[["alignment.from"]], "-to-", pfam_alignment[["alignment.to"]] )
  
  # 2. Create IRanges object
  ir.prot <- create.iranges.pfam(pfam_alignment, protein_id)
  # 3. Extract genomic coordinates
  prot.genome.coord <- extract.prot.genomic.coord.iranges(iranges_prot = ir.prot, edbx = edbx)
  prot.genome.coord.df <- lapply(prot.genome.coord, as.data.frame)[[1]]
  ## 4. Add alignment and transcript ID
  n.exons <- nrow(prot.genome.coord.df)
  pfam.alignment.id.df <-data.frame(
    "PFAM.alignment.ID" = rep(pfam.alignment.id, times = n.exons), 
    "Transcript.ID" = rep(pfam_alignment[["transcript.id"]], times = n.exons))
  ### Cbind pfam.alignment ID
  prot.genome.coord.df <- cbind(pfam.alignment.id.df, prot.genome.coord.df)
  
  return(prot.genome.coord.df)
}

get.protein_id.from.transcript_id <- function(transcript_id, ensembl_mart){
  ## Retrieve Ensembl protein ID from Ensembl transcript ID
  prot.id <- getBM(attributes=c('ensembl_transcript_id','ensembl_peptide_id'),   
                  filters = c('ensembl_transcript_id'), 
                  values = transcript_id,
                  mart = ensembl, 
                  useCache = F)
  if(nrow(prot.id) == 0){
    protein_id <- "ENSP.Not.Available"
  }else{
    protein_id <- prot.id[["ensembl_peptide_id"]]
  }
  return(protein_id)
}

# 0. Params
transcript_id <- opt$transcript_id
## Get protein ID from Ensembl based on Transcript ID
suppressPackageStartupMessages(require(biomaRt))
ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
protein_id <- get.protein_id.from.transcript_id(transcript_id = transcript_id, ensembl_mart = ensembl)

## 0.1. Database object
edbx <- EnsDb.Hsapiens.v86

# 1. Read files
pfam.alignment <- read.table(opt$pfam_alignment, header = T)
# 2. Extract genomic coordinates from PFAM alignments
genomic.coord.list <- apply(pfam.alignment, 1, function(x) 
  extract.genomic.coord.pfam.alignment(pfam_alignment = x, protein_id = protein_id, edbx = edbx))
## Merge
genomic.coord.df <- Reduce(rbind, genomic.coord.list)
# 3. Save
write.table(genomic.coord.df, file = opt$output_coord, quote = F)
