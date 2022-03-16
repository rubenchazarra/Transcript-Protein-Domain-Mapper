## Transcript Protein Domain Mapper 

Method for retrieving and visualising the protein domains of any protein coding isoform. 

This method extracts the coding sequence (CDS) of any protein coding isoform, and maps it to the protein domain family database PFAM. In addition, it outputs a visualisation of the transcript genomic structure together with its corresponding domains. The pipeline is based in [Nextflow](https://www.nextflow.io/), a bioinformatics workflowflow schedueler.

Developed by [**Ruben Chazarra Gil**](https://github.com/rubenchazarra).


### Pipeline inputs

The input files to the method must be specified in the configuration file (`nextflow.config`). These are: 
* GTF annotation file  (`params.inputs.annot_gtf`)
* List of transcript IDs (`params.inputs.transcript_list`): Isoform nomenclature must match that one present in the GTF annotation file
* Genome fasta file `params.inputs.genome`)


### Description of the method

The pipeline runs once per transcript, in a parallel manner. For each transcript in the input list the following steps are conducted: 
1. Extraction of the transcript CDS and translation (performed with [gffread](https://github.com/gpertea/gffread))
2. Query the [PFAM](http://pfam.xfam.org/) protein domain family database
3. Map the PFAM alignment coordinates to genomic coordinates
4. Visualization of the transcript model together with its PFAM alignments ([GViz](https://bioconductor.org/packages/release/bioc/html/Gviz.html))


### Running the pipeline

The pipeline can be run directly from the commandline as: 

`nextflow run main.nf`

Preferably, you can launch the **Local or Slurm executables**. These contain additional execution parameters (which enable generation of nextflow reports, collection of execution metadata, or resuming previous executions).

Local: 
`./run_AS_Func_Evaluator-LOCAL.sh`

If your infrastructure uses [Slurm](https://slurm.schedmd.com/) workload manager, you can launch the pipeline as: 

`sbatch ./run_AS_Func_Evaluator-SLURM.sh`


#### Fine tunning the pipeline

A pipeline can be tuned by editting aseries of parameters in the `nextflow.config` file 

- The **significance of the sequence alignment to the PFAM domain** can be controlled via `params.query_pfam.evalue_thres` parameter. In principle, a sequence E-value << 1 is recommended for significant alignments, although this value may be modified depending on the analysis purposes. A higher E-value threshold will result in more, less significant alignments

- To **disable the visualisation step** (generates a PDF with the transcript model together with its mapped PFAM domains):  `params.viz == False`. 

- If the **visualisation step is enabled**, and in order to include a chromosome diagram in the plot, you should add to the `params.visualization.cytoBand_table` parameter, the table describing the positions of cytogenetic bands wfor each chromosomes. This file can be downloaded [here](https://genome.ucsc.edu/cgi-bin/hgTables?db=hg38&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_doSchema=describe+table+schema).

- The `params.run_tag` parameter allows to custom the results directory. By default, the method saves its output in a `results/` directory. If `params.run_tag = 2022-01-01`, the output directory will become `results/2022-01-01`.



### Software Requirements

For the successful implementation of the method, a series of tools that have to be installed previously: 

1) **Nextflow**

The pipeline is built in [Nextflow](https://www.nextflow.io/), which can be installed in various ways described [here](https://www.nextflow.io/docs/latest/getstarted.html)

2) **gffread** 

The program [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) is used for GFF file processing. We use gffread to extract and translate the CDS from each transcript. 

You can __install gffread__ from its [GitHub repository](https://github.com/gpertea/gffread) or via [conda](https://anaconda.org/bioconda/gffread)

3) **PFAM database**

In order to perform a local installation of the PFAM database we first have to install [HMMER](http://hmmer.org/). 

[__How to install HMMER__](http://hmmer.org/documentation.html). 

Then add hmmer to your $PATH variable as `ADD`

Then you can download the latest release of the PFAM DB from the [EBI FTP server](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/). This can be done directly from your terminal as: 

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz 
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```
Note that the ftp path will have to be changed, to download a different version of the PFAM database 

4) **R Packages** 

The following CRAN and Bioconductor packages must be installed for the pipeline to be functional: 
* [optpartse](https://cran.r-project.org/web/packages/optparse/index.html)
* [Ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html)
* [EnsDb.Hsapiens.v86](http://bioconductor.org/packages/release/data/annotation/html/EnsDb.Hsapiens.v86.html) [**No longer required**]
* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html)


### Visualisation example
 
 - TO ADD

### Pipeline diagram

This is the flow chart of the pipeline:

![](new_flowchart.png)
