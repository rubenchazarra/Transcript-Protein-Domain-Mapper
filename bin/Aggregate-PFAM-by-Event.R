#!/usr/bin/env Rscript

## Merge the PFAM domain alignments of the different transcripts involved in an Event 

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
    c("-s", "--event_coords"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Event Coordinates to represent in plot (Will not be represented if input is NULL)'
  ), 
  make_option(
    c("-p", "--pfam_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Paths of the PFAM alignments Genomic Coordinate files of the differents Transcripts to represent are located.'
  ), 
  make_option(
    c("-f", "--gtf_path"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Paths of the GTF files of the differents Transcripts to represent are located.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## Functions ## 
collapse.exons.domain <- function(df){
  suppressPackageStartupMessages(require(dplyr))
  ## Collapse individual exonic coordinates of pfam-alignment-domain to a single row
  
  # 1) Obtain overall domain-alignment coordinates (max and min)
  if(all(!is.na(df[["pfam_domain"]]))){
  
  if(unique(df[["strand"]]) == "+") {
    df[["dom.start"]] <- df[["start"]] [ which.min(df[["start"]]) ]
    df[["dom.end"]] <- df[["end"]] [ which.max(df[["end"]]) ]
    
  } else if(unique(df[["strand"]]) == "-") {
    df[["dom.start"]] <- df[["end"]] [ which.max(df[["end"]]) ]
    df[["dom.end"]] <- df[["start"]] [ which.min(df[["start"]]) ]
    
  } else {
    stop("Strand must be one of the characters '+' or '-' ")
  }
  # 3) Collapse dataframe
  variable.cols <- c("exon_number", "exon_id", "width", "start", "end")
  group.cols <- names(df)[!names(df) %in% variable.cols]
  df <- df %>% group_by_at(group.cols) %>% summarise(exon_number = paste(exon_number, collapse = "; "), 
                                                     exon_id = paste(exon_id, collapse = "; "), 
                                                     start = paste(start, collapse = "; "), 
                                                     end = paste(end, collapse = "; "), 
                                                     width = paste(width, collapse = "; "), .groups = 'drop') # irrelevant arg
  as.data.frame(df)
  } else{
    df <- tibble::add_column(df, .after = "event.coords", "dom.start" = NA)
    df <- tibble::add_column(df, .after = "dom.start", "dom.end" = NA)
    as.data.frame(df)
  }
}

add.event.info <- function(df, event_id, gene_id, event_coords) {
  df[["Event.ID"]] <- event_id
  df[["gene.name"]] <- gene_id
  df[["event.coords"]] <- paste(event_coords, collapse = "; ")
  # Order
  ord.cols <- c("Event.ID", "transcript_id", "protein_id", "gene.name", "strand", "seqnames", # event info
                "pfam_alignment_id", "pfam_domain", "pfam_domain_description", "partiality", "protein_start", "protein_end", "cds_ok", # alignment info
                "event.coords", # coordinatess
                "exon_number", "exon_id", "start", "end", "width") # exon collapsed info
  df[, ord.cols]
}

find.coord.overlap <- function(df){
  ## Returns a boolean vector indicating if input ranges overlap
  # 1) Order coordinate pairs (remember we have + and - strands)
  event_coords <- as.numeric(unlist(strsplit(as.character(unique(df[["event.coords"]])), split = "; ")))
  rng = cbind(
    # Event Coordinates
    min(event_coords, na.rm = T), 
    max(event_coords, na.rm = T),
    # Domain Coordinates
    pmin(df[["dom.start"]], df[["dom.end"]]), 
    pmax(df[["dom.start"]], df[["dom.end"]]))
  # 2) Compare ranges
  overlap <- (rng[,1] <= rng[,4]) & (rng[,2] >= rng[,3])
  
  df[["overlap.domain.event"]] <- overlap
  
  return(df)
}


## Execute
# 0. Params
## 0.1. Input params
transcript_ids <- opt$transcript_ids
gene_id <- opt$gene_id
event_id <- opt$event_id
event_coords <- opt$event_coords
## 0.2 Input file paths
pfam_path <- opt$pfam_path
gtf_path <- opt$gtf_path

# 1. Transcript IDs to represent together
replace_chars <- c("[][]", " ") # Replace the square brackets --> nextflow inputs transcript_ids tuple as a value --> "[tr_id_1, tr_id_2]"
for (char in replace_chars){ 
  transcript_ids <-  gsub(transcript_ids, pattern = char, replacement = "")
  event_coords <- gsub(event_coords, pattern = char, replacement = "")
  }
transcript_ids <- unlist(strsplit(transcript_ids , ","))# split by comma
event_coords <- as.numeric(unlist(strsplit(event_coords , ",")))# split by comma
event_coords <- event_coords[!is.na(event_coords)]

# 2. Read and save GTFs
gtf.paths <- list.files(path = gtf_path, pattern = ".gtf", full.names = T)
## Run
gtf.list <- lapply(transcript_ids, function(transcript_id) {
  # Find transcript file
  transcript.gtf <- grep(transcript_id, gtf.paths, value = T)
  # Read
  gtf <- rtracklayer::import(transcript.gtf)
  gtf <- as.data.frame(gtf)
  gtf
})
## Add names
names(gtf.list) <- transcript_ids
## Save GTFs for following processes
saveRDS(gtf.list, file = paste0(paste(transcript_ids,collapse = "-"), "_GTF_list.rds"))

# 3. Read PFAM Genomic Coordinates 
pfam.paths <- list.files(path = pfam_path, pattern = ".txt", full.names = T)
## Run
pfam.list <- lapply(transcript_ids, function(transcript_id) {
  # Find transcript file
  transcript.pfam <- grep(transcript_id, pfam.paths, value = T)
  pfam.gc <- read.table(file = transcript.pfam, header = T, sep = "\t", quote = "\"")
  # Order by exon number
  pfam.gc <- pfam.gc [ order(pfam.gc[["exon_number"]], decreasing = F), ]
  pfam.gc
})
## Add names
names(pfam.list) <- transcript_ids

# 4. Add Event Information 
pfam.list <- lapply(pfam.list, function(df) add.event.info(df = df, event_id, gene_id, event_coords))

# 5. Collapse PFAM Domain Genomic Coordinate DF to one row per Exon
pfam.collapse.list <- lapply(pfam.list, collapse.exons.domain)

# 6. Domain-Event Coordinate overlap
pfam.collapse.list <- lapply(pfam.collapse.list, find.coord.overlap)

## PER DOMAIN COORDINATES
# 7. Collapse per Exon
pfam.collapse.df <- as.data.frame(data.table::rbindlist(pfam.collapse.list))
# 8. Save Event Domains (Collapsed, 1 row per domain) 
## 8.1 Txt
file.name <- paste(gene_id, event_id, paste0(unlist(transcript_ids), collapse = "_"), "DF.txt", sep = "_")
write.table(pfam.collapse.df, file = file.name, sep = "\t", row.names = F, col.names = T, quote = F)
## 8.2 RDS
file.name <- paste(gene_id, event_id, paste0(unlist(transcript_ids), collapse = "_"), "DF.rds", sep = "_")
saveRDS(pfam.collapse.df, file = file.name)

## PER EXON COORDINATES
# 9. Add Domain-Event Coordinate overlap to Per-Exon coord data.frames
pfam.list <- mapply(function(gc.df, collapse.df){
  merge(gc.df, collapse.df[, c("pfam_alignment_id", "overlap.domain.event")], by = "pfam_alignment_id", all.x = T, all.y = T)
}, pfam.list, pfam.collapse.list, SIMPLIFY = FALSE)

# 10. Save Event Domains (per-Exon) for Visualisation step
pfam.df <- as.data.frame(data.table::rbindlist(pfam.list))
file.name <- paste(event_id, paste0(unlist(transcript_ids), collapse = "_"), "Visualisation.rds", sep = "_")

saveRDS(pfam.df, file = file.name)