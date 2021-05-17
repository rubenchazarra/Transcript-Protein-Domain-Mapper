#!/usr/bin/env Rscript

# Pre process Alternative Splicing (AS) input files
## Output various things: i) df with functional concordance, ii) df with no-functional concordance, iii) list of unique transcript ids for PFAM query

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_suppa"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to SUPPAs input in rds fomrat. This contains the Alternative Splicing event ids and metadata'
  ),
  make_option(
    c("-x", "--count_matrix"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to TPM count matrix where rows are the event_ids contained in SUPPA file and columns are samples. In rds format.'
  ),
  make_option(
    c("-g", "--annot_gtf"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to Annotation GTF file.'
  ), 
  make_option(
    c("-c", "--concordant"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output Path to dataframe with data from AS events where main transcripts involved have SAME function, according to GTF.'
  ), 
  make_option(
    c("-n", "--non_concordant"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output Path to dataframe with data from AS events where main transcripts involved have DIFFERENT function, according to GTF.'
  ), 
  make_option(
    c("-t", "--transcript_ids"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Output Path to list of unique transcripts involved in CONCORDANT AS events, and which are CODING.'
  )
)


opt <- parse_args(OptionParser(option_list=option_list))

## Functions ##

# 1. Extract AS features from IOE
extract.ioe.feature <- function(ioe.df, feature, transcript){
  # Extract feature from SUPPA's output file (IOE)
  if (transcript == TRUE){ # this line is for extracting transcript vectors
    # separate transcripts vector
    unlist(strsplit(ioe.df[[feature]], ","))
  } else {
    ioe.df[[feature]]
  }
}
# 2. Extract features from AS event id 
extract.event.id.feature <- function(event.id, feat){
  ## Extract features from AS event ID
  if(!(feat %in% c("event.type", "strand"))) stop ("Feat must be one of: 'event.type', 'strand'.")
  if(feat == "event.type"){feat.extract = strsplit(strsplit(event.id, split = ";")[[1]][[2]], ":")[[1]][1] }
  if(feat == "strand"){feat.extract = unlist(lapply(strsplit(event.id, ":"), function(x) tail(x, 1))) }
  return(feat.extract)
}

# 3. Extract expression 
median.exp <- function(transcript.vec, counts.mat){
  ## Retrieve a list of median expression values for a given vector of transcripts
  l <- list()
  for (t in transcript.vec){
    if(!(t %in% rownames(counts.mat))) stop(paste0("Transcript : ", t, " not present in csount matrix provided."))
    # extract median value
    l[[t]] =  median(as.numeric(counts.mat[t, ]))
  }
  l
}
# 4. Select max expression
select.max.exp <- function(transcript.list){
  ## Select isoform from list with higher expression 
  transcript.list[which.max(transcript.list)]
}
# 5. Query GTF --------------------------------> THIS CAN BE IMPROVED! Subset by transcript_ID first
query.gtf <- function(gtf, transcript.id, features.extract){
  ## Extract values from GTF file
  t.gtf <- gtf[gtf$type == "transcript" , ] # this line is required
  if(!(transcript.id %in% t.gtf$transcript_id)) stop (paste0(transcript.id, "not present in 'transcript_id' col of provided GTF file"))
  t.gtf <- t.gtf[t.gtf$transcript_id == transcript.id, ]
  # coerce to df and extract features
  data.frame(t.gtf[, features.extract])
}


extract.event.id.transcripts <- function(ioe.df, tpm.transcript){
  ## Extract metadata from most abundant transcripts composing an AS event
  
  # 1. Extract AS features from IOE
  gene.id <- extract.ioe.feature(ioe.df = ioe.df, feature = "gene_id", transcript = F)
  event.id <- extract.ioe.feature(ioe.df = ioe.df, feature = "event_id", transcript = F)
  ## 1.2 Extract features from event-id
  event.type <- extract.event.id.feature(event.id = event.id, feat = "event.type")
  strand <- extract.event.id.feature(event.id = event.id, feat = "strand")        
  
  # 2. Extract transcripts from IOE
  ## Alternative T
  alt.transcripts <- extract.ioe.feature(ioe.df = ioe.df, feature = "alternative_transcripts", transcript = T)
  ## Total T
  tot.transcripts <- extract.ioe.feature(ioe.df = ioe.df, feature = "total_transcripts", transcript = T)
  ## Filter Total T
  tot.transcripts <- tot.transcripts[!(tot.transcripts %in% alt.transcripts)]
  
  # 3. Extract expression 
  ## Alternative T med expression
  alt.median.exp <- median.exp(transcript.vec = alt.transcripts, counts.mat = tpm.transcript)
  ## Total T med expression
  tot.median.exp <- median.exp(transcript.vec = tot.transcripts, counts.mat = tpm.transcript)
  
  # 4. Select max expression
  ## Alternative T
  alt.max = select.max.exp(alt.median.exp)
  ## Total T
  total.max = select.max.exp(tot.median.exp)
  
  # 5. Query GTF
  ## Alternative T
  alt.gtf.data = query.gtf(gtf, transcript.id = names(alt.max), features.extract = c("transcript_id", "transcript_type", "phase"))
  ## Total T
  total.gtf.data = query.gtf(gtf, transcript.id = names(total.max), features.extract = c("transcript_id", "transcript_type", "phase"))
  
  # 6. Output df
  out.df <- data.frame("gene_id" = gene.id,
                       "event_id" = event.id, 
                       "event_type"= event.type, 
                       "strand" = strand,
                       
                       "tr_1_ID" = names(alt.max),
                       "tr_1_TPM" = unlist(alt.max),
                       "tr_2_ID" = names(total.max), 
                       "tr_2_TPM" = unlist(total.max), 
                       
                       "tr_1_Func" = alt.gtf.data$transcript_type, 
                       "tr_2_Func" = total.gtf.data$transcript_type, 
                       
                       "Func_concord" = alt.gtf.data$transcript_type == total.gtf.data$transcript_type
                       , row.names = 1)
  out.df
}


# Read input files
## 1. SUPPA's output
ioe.AS <- readRDS(opt$input_suppa)
## 2. TPM count matrix
tpm.AS <- readRDS(opt$count_matrix)
## 3. Annotation GTF -------------------> TODO: subset GTF to transcripts present in event IDS (This will speed up the code, since GTF object won't be that big!)
suppressPackageStartupMessages(require(rtracklayer))
gtf <- import(opt$annot_gtf)


# Run (apply by row, this is, by event_id)
df.list = apply(ioe.AS, MARGIN = 1, function(AS.event) extract.event.id.transcripts(ioe.df = AS.event, tpm.transcript = tpm.AS))
for(i in 1:nrow(ioe.AS)){
  extract.event.id.transcripts(ioe.df = ioe.AS[i, ], tpm.transcript = tpm.AS)
}
# Merge dfs
results.df = Reduce(rbind, df.list)
rownames(results.df) <- c(1:nrow(results.df))

# Split df into : Funcitonal concordant and Non.concordant
func.concordance = results.df[results.df$Func_concord == TRUE, ]
func.no.concordance = results.df[results.df$Func_concord == FALSE, ]

# Save 1. Concordance
write.table(func.concordance, opt$concordant, sep = "\t", quote = F, row.names = T)
# Save 2. Non-concordance
write.table(func.no.concordance,opt$non_concordant, sep = "\t", quote = F, row.names = T)

# Process 1. Concordance
## Subset to Coding-Coding events
func.concordance.coding = func.concordance[func.concordance$tr_1_Func == "protein_coding" & func.concordance$tr_2_Func == "protein_coding", ]
## Extract list of unique transcripts from the functional concordant and coding
transcripts = c(as.character(func.concordance.coding$tr_1_ID), as.character(func.concordance.coding$tr_2_ID))
## Remove Gencode dot after transcript ID
transcripts = unlist(lapply(strsplit(transcripts, "[.]"), function(x) x[1]))
uniq.transcripts = unique(transcripts)
# Save transcript List
write.table(uniq.transcripts, opt$transcript_ids, quote = F, row.names = F,  col.names=FALSE)
