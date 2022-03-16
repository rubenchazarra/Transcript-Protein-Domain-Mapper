#!/usr/bin/env Rscript

## Script to map PFAM domain alignment coordinates to genomic coordinates
## This is done manually (with a GTF annotation file), without the need of Ens Annotation package

## Note: The 'cds_ok' column in the output depicts if there is a coincidence between the pfam genomic coordinates and the length of the transcript CDS 
######## This is the only thing for which we need the protein CDS sequence so we could potentially remove this

suppressPackageStartupMessages(require(optparse))

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
    c("-p", "--protein_fasta"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to the protein fasta file corresponding to transcript_id.'
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

## Functions ##

# 1) Read fasta file protein sequence to variable
read_fasta_prot <- function(file.path){
  ## Read fasta file protein sequence to variable
  if (missing(file.path)) stop("Argument 'file.path' is required, and must be a path to a PFAM dmtable output.")
  
  line.list = readLines(file.path)
  uncommented.lines <- grep(">", line.list)
  all.lines <- 1:length(line.list)
  filled.lines <- all.lines[!all.lines %in% uncommented.lines]
  seq <- paste(line.list[filled.lines], collapse = "")
  return(seq)
}

# 2) Create IRanges object from PFAM alignment coordinates
create.iranges.pfam <- function(pfam_alignment, protein_id){
  ## Create IRanges object from PFAM alignment coordinates
  if (missing(pfam_alignment)) stop("Argument 'pfam_alignment' is required.")
  if (missing(protein_id)) stop("Argument 'protein_id' is required.")
  
  pfam.al.start <- as.numeric(pfam_alignment[["ali_from"]])
  pfam.al.stop <- as.numeric(pfam_alignment[["ali_to"]])
  n.alignments <- length(pfam.al.start)
  protein_names <- rep(protein_id, times = n.alignments)
  # Create IRanges object
  ir.prot <- IRanges::IRanges(start = pfam.al.start, end = pfam.al.stop, names = protein_names) # Package
  return(ir.prot)
}

# 3)  Map protein to genomic coordinates
protein_to_genome <- function(ir.pfam, gtf, prot_seq, transcript_id) {
  ## Map protein-relative coordinates to genomic coordinates
  ## Code partially modified from ensembldb:proteinToGenome()
  suppressPackageStartupMessages(require(ensembldb))
  suppressPackageStartupMessages(require(GenomicRanges))
  
  if (missing(ir.pfam) || !is(ir.pfam, "IRanges"))
    stop("Argument 'ir.pfam' is required and has to be an 'IRanges' object")
  if (missing(gtf) || !is(gtf, "GRanges"))
    stop("Argument 'gtf' is required and has to be an 'GRanges' object, imported with rtracklayer::import")
  if (missing(prot_seq) || !is(prot_seq, "character"))
    stop("Argument 'prot_seq' is required and has to be an 'character' object")
  if (missing(transcript_id)) stop("Argument 'transcript_id' is required")
  
  # 1)  Convert aa coordinates from 'ir.pfam' to nucleotide (*3) 
  coords_cds <- protein_coords_to_tx(ir.pfam)
  
  # 2) Retrieve CDS for each protein in ir.pfam (there's only one)
  message("Fetching CDS for ", length(ir.pfam), " proteins ... ",
          appendLF = FALSE)
  
  # 3) Extract CDS coordinates of 'transcript_id'
  ## Original function: cds_genome <- .cds_for_id_range(db, x, id = id, idType = idType)
  cds_genome <- extract_cds_coords_transcript(ir.pfam = ir.pfam, 
                                              transcript_id = transcript_id, 
                                              gtf = gtf)
  miss <- lengths(cds_genome) == 0
  if (any(miss)) {
    warning("No CDS found for: ", paste0(names(cds_genome)[miss],
                                         collapse = ", "))
  }
  message(sum(!miss), " found")
  
  ### 4) Ensure that the CDS matches the AA sequence length
  message("Checking CDS and protein sequence lengths ... ", appendLF = FALSE)
  ## Original function: cds_genome <- .cds_matching_protein(db, cds_genome)
  cds_genome <- cds_matching_protein(cds = cds_genome, prot_seq = prot_seq)
  
  are_ok <- unlist(lapply(cds_genome, function(z) {
    if (is(z, "GRangesList"))
      all(z[[1]]$cds_ok)
    else NA
  }))
  are_ok <- are_ok[!is.na(are_ok)]
  ## We've got now a list of GRanges
  message(sum(are_ok), "/", length(are_ok), " OK")
  
  ### 5) Perform the mapping for each input range with each mapped cds
  res <- mapply(
    cds_genome, as(coords_cds, "IRangesList"), as(ir.pfam, "IRangesList"),
    FUN = function(gnm, cds, prt) {
      if (is.null(gnm)) {
        GRanges()
      } else {
        ## Unlist because we'd like to have a GRanges here. Will split
        ## again later.
        
        ## New function here
        maps <- unlist(ensembldb:::.to_genome(gnm, cds))
        ## Don't want to have GRanges names!
        names(maps) <- NULL
        mcols(maps)$protein_start <- start(prt)
        mcols(maps)$protein_end <- end(prt)
        maps[order(maps$exon_number)]
      }
    })
  ## Split each element again, if there are more than one protein_id. Names
  ## of the elements are then the protein_id.
  lapply(res, function(z) {
    if (length(unique(z$protein_id)) > 1)
      split(z, f = z$protein_id)
    else z
  })
}

# 4) Convert Protein to Transcript coords
protein_coords_to_tx <- function(x) {
  ## Convert protein Coordinates in x to Transcript
  # From ensembldb:::.proteinCoordsToTx()
  if (missing(x) || !is(x, "IRanges")) stop("Argument 'x' is required and has to be an 'IRanges' object")
  IRanges::end(x) <- IRanges::end(x) * 3
  IRanges::start(x) <- 1 + (IRanges::start(x) - 1) * 3
  x
}

# 5) Extract Transcript CDS coordinates with Check
extract_cds_coords_transcript <- function(ir.pfam, gtf, transcript_id){
  ## Extract CDS from a transcript, and check the returned CDS match the provided protein ranges (from PFAM)
  
  if (missing(ir.pfam) || !is(ir.pfam, "IRanges")) stop("Argument 'ir.pfam' is required and has to be an 'IRanges' object")
  if (missing(gtf) || !is(gtf, "GRanges")) stop("Argument 'gtf' is required and has to be an 'GRanges' object")
  if (missing(transcript_id)) stop("Argument 'transcript_id' is required.")
  
  
  # 1. Input transcript_id
  n.alignments <- nrow(data.frame(ir.pfam))
  transcript_ids <- rep(transcript_id, times = n.alignments)
  
  # 2. Coerce pfam.al to IRanges
  protein_ids <- names(ir.pfam)
  
  # 3. Extract CDS of transcript_ids
  cds.granges.list <- lapply(transcript_ids, function(tr_id) extract_cds_coords_granges(gtf = gtf, transcript_id = tr_id))
  names(cds.granges.list) <- protein_ids
  # 4. Check returned CDS if their width matches the provided protein ranges 
  # TODO: TBH I don't know what they are doing here (extracted from .cds_for_id_range() function)
  mapply(cds.granges.list, end(ir.pfam) * 3, # This is multiplying aa coord to genomic
         FUN = function(x, x_end) {
           if (!is.null(x))
             x[sum(width(x)) >= x_end]
         })
  
}

# 6) Extract Transcript CDS coords as GRangesList object
extract_cds_coords_granges <- function(gtf, transcript_id){
  ## Extract Transcript CDS coords as GRangesList object from 'gtf'
  if (missing(gtf) || !is(gtf, "GRanges")) stop("Argument 'gtf' is required and has to be an 'GRanges' object")
  if (missing(transcript_id) ) stop("Argument 'transcript_id' is required")
  
  # Desired output fields
  desired.fileds <- c("protein_id", "transcript_id", "exon_id", "exon_number")
  # 1. Subset gtf to 'transcript_id' and to elements being 'CDS'
  transcript.gtf <- gtf[grep( transcript_id, gtf$transcript_id ),  ]
  
  # 2. CDS GTF
  cds.gtf <- transcript.gtf[transcript.gtf$type == "CDS", ]
  coding.exon.numbers <- cds.gtf$exon_number
  
  # 3. Extract Exon_IDs from exon_numbers
  exon.gtf <- transcript.gtf[transcript.gtf$type == "exon", ]
  
  exon.id.vec <- c()
  for (exon.num in coding.exon.numbers){
    exon.df <- data.frame(exon.gtf[exon.gtf$exon_number == exon.num, ])
    exon.id <- exon.df[["exon_id"]]
    exon.id.vec <- c(exon.id.vec, exon.id)
  }
  # 4. Append Exon IDS to CDS GTF
  cds.gtf$exon_id <- exon.id.vec
  
  # 5. Coerce to GRangesList class object - This is to make equivalent to the output of the .cds_for_id_range() function
  cds.granges <- GRangesList(cds.gtf[, desired.fileds], compress = T) 
  cds.granges@partitioning@NAMES <- transcript_id
  return(cds.granges)
}

# 7) Ensure that the CDS matches the AA sequence length
cds_matching_protein <- function(cds, prot_seq){
  ## Ensure that the CDS matches the AA sequence length
  ## Based on the ensembldb:::..cds_matching_protein() function 
  if (missing(cds)) stop("Argument 'cds' is required.")
  if (missing(prot_seq)) stop("Argument 'prot_seq' is required.")
  
  
  ## 1) Fetch the protein sequences, all in one go.
  prot_ids <- unique(unlist(lapply(cds, function(z) { # Recyclable Line
    lapply(z, function(y) y$protein_id) }), use.names = FALSE))
  if (length(prot_ids) == 0)
    return(lapply(cds, function(z) NULL))
  
  ## 2) Generate df with protein sequence
  prot_seqs <- data.frame("protein_id" = prot_ids, "protein_sequence" = prot_seq)
  
  ## 3) Calculate the expected CDS length (add +3 +3 to add the start and stop codon).
  ## This is due to the way we are extracting the CDS in the extract.cds.as.granges function. Start and Stop Codons are included in the CDS, this is why the +3 `+3
  exp_cds_len <- nchar(as.character(prot_seqs$protein_sequence)) * 3 
  names(exp_cds_len) <- prot_seqs$protein_id
  
  # 4) Determine if the CDSs are OK
  mapply(cds, names(cds), FUN = function(z, nm) {
    if (!is.null(z)) {
      cds_lens <- sum(width(z))
      prt_ids <- unlist(lapply(z, function(y) unique(y$protein_id)))
      diffs <- cds_lens - exp_cds_len[prt_ids]
      ## Return all for which diff is 0, or otherwise the one with the
      ## smallest diff.
      if (any(diffs == 0)) {
        z <- lapply(z, function(grng) {
          mcols(grng)$cds_ok <- TRUE
          grng
        })
        GRangesList(z[diffs == 0])
      } else {
        warning("Could not find a CDS whith the expected length for ",
                "protein: '", nm, "'. The returned genomic ",
                "coordinates might thus not be correct for this ",
                "protein.", call. = FALSE)
        ## Alternatively we could align the RNA and AA sequences and
        ## trim the protein sequence...
        z <- lapply(z, function(grng) {
          mcols(grng)$cds_ok <- FALSE
          grng
        })
        GRangesList(z[which.min(abs(diffs))])
      }
    }
  })
  
}

# 8) Integrate PFAM alignment data with PFAM genomic coordinates
integrate.al.pfam.with.coords <- function(al.pfam, pfam_gcords){
  # Integrate PFAM alignment data with PFAM genomic coordinates
  # 1. Select chrN
  chrN <- unique(seqnames(pfam_gcords))
  # 2. Create PFAM alignment identifier string
  al.pfam.id <- paste0(al.pfam[["query_name"]], "-",  al.pfam[["domain_name"]],
                       "-from-", as.numeric(al.pfam[["ali_from"]]),
                       "-to-", as.numeric(al.pfam[["ali_to"]]), "-", chrN ) 
  # 3. Number of exons in pfam_gcords
  n.exons <- pfam_gcords@seqnames@lengths 
  # 4. PFAM Metadata DF
  al.pfam.id.df <- data.frame("pfam_alignment_id" = rep(al.pfam.id, times = n.exons), 
                              "pfam_domain" = rep(al.pfam[["domain_name"]], times = n.exons), 
                              "pfam_domain_description" = rep(al.pfam[["description"]], times = n.exons)
  )
  # 5. Merge 2 dfs
  cbind(al.pfam.id.df, pfam_gcords)
}

# 9) Create empty PFAM data frame
create.empty.iranges.pfam <- function(transcript_id, chrN) {
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
  empty.align[["pfam_alignment_id"]] = paste0(transcript_id, "-NA", "-from-NA-to-NA-", chrN )
  empty.align[["seqnames"]] = chrN
  empty.align[["transcript_id"]] = transcript_id
  
  empty.pfam.df = data.frame(t(empty.align))
  return(empty.pfam.df)
}


### Run
suppressPackageStartupMessages(require(ensembldb)) # Note: using some functionalities of ensembldb, but not the EnsDb.Hsapiens.v86 annotation package
suppressPackageStartupMessages(require(rtracklayer))

# 0. Params
transcript_id <- opt$transcript_id

# 1. Read Files
## 1.1. Transcript GTF 
gtf <- rtracklayer::import(opt$gtf)
# Protein ID ( remove suffix )
protein_id <- unlist(strsplit(unique(gtf$protein_id), "[.]"))[1]

## 1.2. Protein sequence
prot_seq <- read_fasta_prot(file.path = opt$protein_fasta)

## 1.3. Pfam Alignment
al.pfam <- read.table(opt$pfam_alignment, header = T, sep = "\t", quote="\"" ) # Note: this quoting symbol is to avoid quotes in string of pfam alignment description (5' for instance) to be understood as quoting character

# 2 Protein --> Genomic Coordinates
## Condition: if PFAM alignment contains some output
if(all(al.pfam[["pfam_match"]] == TRUE)){
    # 1. Create IRanges object 
    ir.pfam <-  create.iranges.pfam(al.pfam, protein_id) # create IRanges object from PFAM alignment
    # 2. Convert protein relative (PFAM alignment) coordinates to genomic
    pfam_gcords <- protein_to_genome(ir.pfam = ir.pfam, gtf = gtf, prot_seq = prot_seq, transcript_id) # this is a list
    # 3. Split PFAM alignment to list
    al.pfam.list <- split(al.pfam, seq(nrow(al.pfam))) # each row is an element of the list now
    # 4. Integrate PFAM domain data with PFAM alignment genomic coordinates
    pfam.coords.list <- mapply(function(al.pfam, pfam_gcords_i) integrate.al.pfam.with.coords(al.pfam, pfam_gcords_i), al.pfam.list, pfam_gcords, SIMPLIFY = F)
    # 5. Merge dfs to single df
    pfam.coords.df <- data.table::rbindlist(pfam.coords.list)
    # 6. Order df bu 'exon_number'
    #pfam.coords.df <-  pfam.coords.df[ order(as.numeric(pfam.coords.df$exon_number, decreasing = F)), ]
    # 7. Add partiality
    chrN <- as.character(unique(pfam.coords.df[["seqnames"]]))
    al.pfam[["pfam_alignment_id"]] <- paste0(al.pfam[["query_name"]], "-", al.pfam[["domain_name"]], "-from-", al.pfam[["ali_from"]], "-to-", al.pfam$ali_to, "-", rep(chrN, times = nrow(al.pfam)))
    cols.interest <- c("pfam_alignment_id", "partiality")
    pfam.coords.df <- merge(x = pfam.coords.df, y = al.pfam[, cols.interest], by = "pfam_alignment_id")
    # If PFAM alignment is empty
  } else if (all(al.pfam[["pfam_match"]] == FALSE)) {
    chrN <- unique(seqnames(gtf))
    # Create empty data.frame if no PFAM hit was retrieved
    pfam.coords.df <- create.empty.iranges.pfam(transcript_id = transcript_id, chrN = chrN)
  }

# 3. Save 
write.table(pfam.coords.df, file = opt$output_coord, sep = "\t", quote = F )
