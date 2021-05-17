# AS_Function_Evaluator

This pipeline evaluates the functional consequences of alternative splicing events. 

- ADD MORE. 

## Goal of the project

- ADD. 

## Description of the method

The input files to this pipeline must be specified in the `nextflow.config` file. These are: 
* List of transcript IDs
* GTF annotation file
* Genome fasta file

For each input transcript, the pipeline performs the following steps: 
1. Extraction of the transcript protein sequence
2. Query the PFAM database
3. Map the PFAM alignment coordinates to genomic coordinates
4. Visualization of the transcript model + the PFAM alignment

- TO ADD --> Mention the SUPPA script. 

## Running the pipeline

The pipeline can be run directly from the commandline as: 

`nextflow run main.nf`

Preferably, you can execute: 

`bash run_AS_Func_Evaluator/run_AS_Func_Evaluator.sh`

Which contains additional execution parameters for generation of nextflow reports, and collection of execution metadata. 

Additionally, if your infrastructure uses [Slurm](https://slurm.schedmd.com/) workload manager, you can launch the pipeline as: 
```
sbatch run_AS_Func_Evaluator/run_AS_Func_Evaluator/submit_run_AS_Func_Evaluator.sh
```

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
