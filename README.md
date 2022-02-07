# Transcript Protein Domain Mapper 

Method for retrieving and visualising the protein domains of any protein coding isoform. 

This method extracts the coding sequence (CDS) of any protein coding isoform, and maps it to the protein domain family database PFAM. In addition, it outputs a visualisation of the transcript genomic structure together with its corresponding domains. The pipeline is based in [Nextflow](https://www.nextflow.io/), a bioinformatics workflowflow schedueler.

Developed by [**Ruben Chazarra Gil**](https://github.com/rubenchazarra).


## Pipeline inputs

The input files to the method must be specified in the configuration file (`nextflow.config`). These are: 
* GTF annotation file  (`params.inputs.annot_gtf`)
* List of transcript IDs (`params.inputs.transcript_list`): Isoform nomenclature must match that one present in the GTF annotation file
* Genome fasta file `params.inputs.genome`)

## Description of the method

For each transcript in the input list, the pipeline performs the following steps: 
1. Extraction of the transcript CDS and translation (performed with [gffread](https://github.com/gpertea/gffread)
2. Query the [PFAM](http://pfam.xfam.org/) protein domain family database
3. Map the PFAM alignment coordinates to genomic coordinates
4. Visualization of the transcript model together with its PFAM alignments


## Running the pipeline

The pipeline can be run directly from the commandline as: 

`nextflow run main.nf`

Preferably, you can launch the **Local or Slurm executables**. These contain additional execution parameters (which enale generation of nextflow reports, collection of execution metadata, or resuming previous executions).

Local: 
`./run_AS_Func_Evaluator-LOCAL.sh`

or

`sbatch ./run_AS_Func_Evaluator-SLURM.sh`


Additionally, if your infrastructure uses [Slurm](https://slurm.schedmd.com/) workload manager, you can launch the pipeline as: 
```
sbatch run_AS_Func_Evaluator/run_AS_Func_Evaluator/submit_run_AS_Func_Evaluator.sh
```

## Requirements [check below, some of it already described ]
- Installation of HMMER
- Local installation of PFAM
- Install GFFread

--> R packages
- ensembldb
- GViz



## Installation Notes

The pipeline requires the following sotfware to be installed: 

### 1. Nextflow

The pipeline is built in [Nextflow](https://www.nextflow.io/), a bioinformatics workflow tool.

### 2. gffread tool
The program [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) can be used to validate, filter, convert and perform various other operations on GFF files. In this pipeline, we use gffread to extract the protein sequence of a protein coding transcript given: i) a GTF annotation file, ii) genome fasta sequence file and iii) a valid transcript ID for the provided GTF.

You can __install gffread__ from its [GitHub repository](https://github.com/gpertea/gffread) or via [conda](https://anaconda.org/bioconda/gffread)

### 3. PFAM database

In order to perform a local installation of the PFAM database we first have to install [HMMER](http://hmmer.org/). 

[__How to install HMMER__](http://hmmer.org/documentation.html). 

Then add hmmer to your $PATH variable.

Next you will have to download one of the releases of the PFAM database file from the [EBI FTP server](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/). In order to download and install PFAM directly from your terminal you can issue the following commands: 

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz 
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```


### 4. R Packages 

The following CRAN and Bioconductor packages must be installed for the pipeline to be functional: 
* [optpartse](https://cran.r-project.org/web/packages/optparse/index.html)
* [Ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html)
* [EnsDb.Hsapiens.v86](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html)
* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html)



## Pipeline diagram

This is the flow chart of the pipeline:

![](new_flowchart.png)
