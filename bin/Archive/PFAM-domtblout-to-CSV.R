#!/usr/bin/env Rscript

# Script to read the PFAM output (in domtblout format) into a CSV table.
## Note: Local PFAM Db output is different fromt he PFAM REST API output.

suppressPackageStartupMessages(require("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_file"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Json object output by PFAM'
  ),
  make_option(
    c("-t", "--transcript_id"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Transcript ID of the sequence that was aligned to PFAM.'
  ),
  make_option(
    c("-o", "--output_table"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output table containing json info'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

## Functions ##
read.PFAM.domtblout <- function(file.path, transcript.id){
  ## Read the PFAM output (in domtblout format) 
  
  # 1. Names to add to the PFAM alignment table
  table.names <- c("transcript.id", # Transcript ID as added by the transcript_id argument
                   "pfam.name", # The name of the target profile
                   "pfam.accesion", # Accession of the target profile, or ’-’ if none is available
                   "hmm.length", #  Length of the target profile, in residues. This (together with the query length) is useful for interpreting where the domain coordinates (in subsequent columns) lie in the sequence.
                   "transctipt.id", # Name of the query sequence or profile.
                   "accession.2", #  Accession of the target sequence or profile, or ’-’ if none is available
                   "query.length", #  Length of the query sequence or profile, in residues.
                   "E.value.Overall", # E-value of the overall sequence/profile comparison (including ALL domains)
                   "Bit.score.Overall", # Bit score of the overall sequence/profile comparison (including all domains), inclusive of a null2 bias composition correction to the score.
                   "Bias.Overall", # The biased composition score correction that was applied to the bit score.
                   "domain.number", # This domain’s number 
                   "n.domain", # The total number of domains reported in the sequence, ndom.
                   "c.Evalue", # The “conditional E-value”, a permissive measure of how reliable this particular domain may be.
                   "i.Evalue", # The “independent E-value”, the E-value that the sequence/profile comparison would have received if this were the only domain envelope found in it, excluding any others.
                   "Bit.score.Domain", # The bit score for this domain.
                   "bias.Domain", # The biased composition (null2) score correction that was applied to the domain bit score.
                   "hmm.from", # The start of the MEA alignment of this domain with respect to the PROFILE
                   "hmm.to", # The end of the MEA alignment of this domain with respect to the PROFILE
                   "alignment.from", # The start of the MEA alignment of this domain with respect to the SEQUENCE
                   "alignment.to", # The end of the MEA alignment of this domain with respect to the SEQUENCE
                   "envelope.from", # The start of the domain ENVELOPE on the SEQUENCE
                   "envelope.to", # The end of the domain ENVELOPE on the SEQUENCE
                   "accuracy", # The mean posterior probability of aligned residues in the MEA alignment;  a measure of how reliable the overall alignment is (from 0-1)
                   "target.description", # The remainder of the line is the target’s description line, as free text
                   "pfam.match")
  # 2. Read 
  line.list = readLines(file.path)
  uncommented.lines <- grep("#", line.list)
  
  ## 2.1. PFAM Empty Alignment (if all lines are uncommented "#")
  if(length(uncommented.lines) == length(line.list)){
    
    empty.align <- rep(c(NA), times = length(table.names))
    names(empty.align) <- table.names
    ## Add transcript.id
    empty.align[["transcript.id"]] = transcript.id
    ## Add PFAM match logical value
    pfam.match = FALSE
    empty.align[["pfam.match"]] = pfam.match
    pfam.df = data.frame(t(empty.align))
    
    } else {
    
    ## 2.2. Alignment with content.
    all.lines <- 1:length(line.list)
    filled.lines <- all.lines[!all.lines %in% uncommented.lines]
    pfam.align.df <- read.table(file.path, header = F, row.names = NULL, comment.char = "#", sep = "", fill = T) # Note: fill = T to deal with description of protin domain being larger than other
    ## Edits: last column ('pfam.description') is separated since it contains >1 words
    rm.desc.pfam.align.df <- pfam.align.df[, 1:22]
    # This condition is for when PFAM Domain description has only one word. # TODO: IMPROVE THIS 
    if(ncol(pfam.align.df) == 23) {pfam.desc.vec = pfam.align.df[, 23]  
    } else {
    
	pfam.desc.vec <- apply(pfam.align.df[, 23:ncol(pfam.align.df)], 1, paste , collapse = "_" )
    } 
    ### Concatenate back
    pfam.desc.df <- rm.desc.pfam.align.df
    pfam.desc.df[, 23] <- pfam.desc.vec
    rm(rm.desc.pfam.align.df)
    ## Add transcript.id
    transcript.id.df <- data.frame(Transcript.id = rep(transcript.id, times = nrow(pfam.align.df)))
    ## Add PFAM match logical value
    pfam.match = TRUE
    pfam.match.df <- data.frame(pfam.match = rep(pfam.match, times = nrow(pfam.align.df)))
    ## Merge dfs
    pfam.df <- cbind(transcript.id.df, pfam.desc.df, pfam.match.df)
    names(pfam.df) = table.names
    }
    # Close connection 
    return(pfam.df)
}

calculate.frac.alignment.domain <- function(al.to, al.from, hmm.length ){
  ## Calculate the fraction of AA of the HMM model spanned by the PFAM alignment
  ## This is to address Partiality in PFAM alignments. Based on the idea that an alignment of the query sequence (or a part of the query sequence) against 100% of the PFAM domain can be considered a functional domain. Whereas an alignment of 30% of the domain, can't.
  frac.dom = ( as.numeric(al.to) + 1 - as.numeric(al.from)) / as.numeric(hmm.length)
  frac.dom
}

# 1. Coerce domtblout to table
pfam.table <- read.PFAM.domtblout(file.path = opt$input_file, transcript.id = opt$transcript_id)
# 2. Add Fraction of the PFAM Domain (HMM Model) spanned by the Alignment
pfam.table[["Fraction.Al.Domain"]] <- apply(pfam.table[,c('alignment.to', 'alignment.from', 'hmm.length')], 1, function(x) calculate.frac.alignment.domain(al.to = x[1],  al.from = x[2],  hmm.length = x[3] ) )

# 3. Save table
write.table(pfam.table, file = opt$output_table,  quote = F, row.names = F, sep = "\t") # Important to save with a separating character that won't be present in the character strings output of the alignment, if not htis incurs problems in  the following Merging step
