## Transcript Protein Domain Mapper 

Method for retrieving and visualising the protein domains of any protein coding isoform and of different genomes. Enables joint visualisations representing alternative splicing events.

**In brief:** 

This method extracts the coding sequence (CDS) of any protein coding isoform, maps it to the protein domain family database [PFAM](http://pfam.xfam.org/), and retrieves high-confidence domain alignments. In addition, it outputs a visualisation of the transcript genomic structure together with its corresponding domains. It also enables a joint visualisation of different transcripts which may represent an alternative splicing event. This pipeline is writtens in [Nextflow](https://www.nextflow.io/).

Developed by [**Ruben Chazarra Gil**](https://github.com/rubenchazarra).


### Pipeline inputs

The input files to the method must be specified in the configuration file `nextflow.config`. These are: 
* List of transcript IDs (`params.transcript_list`). These are expected to be transcripts annotated as protein coding. An example of input transcript list can be seen in `test_data/Transcript-list-test.csv`.
* GTF annotation file  (`params.gtf`)
* Genome fasta file (`params.genome`)

**Note:** Transcript IDs must match those in the GTF annotation file. Also, the versions of the GTF and the Genome Fasta file should be the same.

### Description of the method

The pipeline runs once per transcript, in a parallel manner. For each transcript in the input list the following steps are conducted: 
1. Extraction of the transcript CDS and translation (performed with [gffread](https://github.com/gpertea/gffread))
2. Query the [PFAM](http://pfam.xfam.org/) protein domain family database
3. Filter PFAM alignments to retain high-confidance domains only
4. Map the PFAM alignment coordinates to genomic coordinates
5. Visualisation of the transcript model together with its PFAM alignments ([GViz](https://bioconductor.org/packages/release/bioc/html/Gviz.html))
6. Aggregated visualisation including more than one transcript, which may represent an alternative splicing event 

### Dependencies

For the successful implementation of the method, a series of tools that have to be installed previously: 

1) **Nextflow**

The pipeline is built in [Nextflow](https://www.nextflow.io/), which can be installed in various ways as described [here](https://www.nextflow.io/docs/latest/getstarted.html).

2) **gffread** 

The tool [gffread](http://ccb.jhu.edu/software/stringtie/gff.shtml#gffread) is used for GFF file processing. We use gffread to extract and translate the CDS from each transcript. You can __install gffread__ from its [GitHub repository](https://github.com/gpertea/gffread) or via [conda](https://anaconda.org/bioconda/gffread). The path to the installation must be specified in the `params.get_cds_translate.gffread` param.

3) **PFAM database**

In order to perform a local installation of the PFAM database we first have to install [HMMER](http://hmmer.org/). 

[__How to install HMMER__](http://hmmer.org/documentation.html)

You must add the hmmscan path to `params.pfam.hmmscan`. Then you can download the latest release of the PFAM DB from the [EBI FTP server](http://ftp.ebi.ac.uk/pub/databases/Pfam/releases/). This can be done directly from your terminal as: 

```
wget ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz 
gunzip Pfam-A.hmm.gz
hmmpress Pfam-A.hmm
```
The pfam database path must then be added to `params.pfam.pfam_db`. Note that the ftp path will have to be changed, to download a different version of the PFAM database.

4) **R Packages** 

The following CRAN and Bioconductor packages must be installed for the pipeline to be functional: 
* [optpartse](https://cran.r-project.org/web/packages/optparse/index.html)
* [Ensembldb](https://bioconductor.org/packages/release/bioc/html/ensembldb.html)
* [rtracklayer](https://bioconductor.org/packages/release/bioc/html/rtracklayer.html)
* [Gviz](https://bioconductor.org/packages/release/bioc/html/Gviz.html)
* [readr](https://readr.tidyverse.org/) (Included in the [tidyverse](https://www.tidyverse.org/) package)

```
install. packages(c("optpartse", "tidyverse")); BiocManager::install("ensembldb", "rtracklayer","Gviz" )
```

### Fine tunning the pipeline

The pipeline can be tuned by editting a series of parameters in the `nextflow.config` file.

- The `params.run_tag` parameter allows to custom the results directory. By default, the method saves its output in a `results/` directory. If `params.run_tag = 2022-01-01`, the output directory will become `results/2022-01-01`.

#### PFAM 

- The **significance** of the **global sequence alignment to the PFAM domain** can be controlled via `params.pfam.sequence_evalue_thres` parameter. In principle, a sequence E-value << 1 is recommended for significant alignments, higher E-value threshold will result in a higher number of less significant alignments. Nonetheless, we recommend to set a permissive sequence E-value threshold (e.g. 1) and tune the PFAM filtering parameters: 
 
    - Domain Score (`params.pfam.domain_score`):  measure of how well a sequence matches the profile HMM (i.e. if the sequence is homologous to the model). 
    - Domain E-value (`params.pfam.domain_evalue`): measures how statistically significant the bit-score is.
    - Accuracy  (`params.pfam.accuracy`): a measure of how reliable the overall alignment is (from 0 to 1, with 1.00 indicating a completely reliable alignment according to the model).
    - Partiality (`params.pfam.partiality`): measures how partial the alignment is with respect to the domain model (e.g. a partiality = 0.9 will only retrieve alignments spanning >= 90% of the domain).


#### Visualisation step

- To **enable the visualisation steps**, we must specify `params.vis = true`. There are two visualisation steps: the first one generates the transcript model with the corresponding mapped domains (`params.vis_transcript = true`); the second one generates an aggregated visualisation including multiple transcripts ( `params.vis_transcript = true`) with their domains.

The **aggregated visualisation step** requires the addition of two additional files: 

1) __Event transcripts__ file (`params.vis.event_transcripts`): specifying the transcripts involved in each alternative splicing event. An example aggregation file can be seen in `test_data/Event-transcripts-test.csv` and should look like:
```
ENSG00000158122.11;SE:chr9:96641886-96645893:96646024-96651390:-,PRXL2C,ENST00000375234.7,ENST00000411939.5
```

- Components: 
    1. First element is the Alternative Splicing Event identifier (can be any id really)
    2. Second element is the gene name
    3. Last elements are the transcripts implicated in the Alternative Splicing Event (these must be present in the input transcript list)

2) __Event coordinates__ file (`params.vis.event_coords`): containing the alternative splicing event coordinates. An example of coordinate file can be seen in `test_data/Event-coords-test.csv` and should look like:
```
ENSG00000000938.12;AF:chr1:27625151-27626046:27626240:27625151-27635065:27635277:-,27635065,27626046,27635277,27626240
```
- Components
    1. Alternative Splicing Event identifier
    2. Genomic coordinates representing the alternative splicing event (Notes AS events might be represented by 2 or 4 coordinates depending on the type of event)
 
 **Note:** You can chose whether you want to see the exclusively the domains that overlap with the event coordinates or all the domains via the parameter: `params.visualisation.show_non_overlapping`.

- In order to include a chromosome diagram in the plot, you should add to the `params.vis.cyto_band` parameter, the table describing the positions of cytogenetic bands for each chromosomes. This file can be downloaded [here](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1311491459_VyO06ty4dMBWFIL998ysQ9Q4AJld&clade=mammal&org=Mouse&db=mm10&hgta_group=map&hgta_track=cytoBand&hgta_table=cytoBand&hgta_regionType=genome&position=&hgta_outputType=primaryTable&hgta_outFileName=) for different species and genome versions.

- The genome version id used must be specified in the `params.visualisation.genome_id` parameter. The genome version ids can be checked [here](https://genome.ucsc.edu/FAQ/FAQreleases.html) 

### Running the pipeline

The pipeline can be run directly from the commandline as: 

`nextflow run main.nf`

Preferably, you can launch the **Local or Slurm executables**. These contain additional execution parameters (which enable generation of nextflow reports, collection of execution metadata, or resuming previous executions).

__Local execution:__ 
`./run-Prot-Dom-Mapper-LOCAL.sh`

__Slurm execution:__ 
`sbatch run-Prot-Dom-Mapper-SLURM.sh`


### Visualisation examples

#### A) Transcript visualisation

One of the visual outputs of the pipeline is a visualisation of the transcript model with the corresponding mapped protein domains from the pfam database:

![](test_output/ENST00000003912.7_Gviz_Trackplot.png)

#### B) Event (Aggregated) visualisation

The 2nd pipeline output is a visualisation representing various transcripts and their pfam mappings. This may represent an Alternative Splicing Event, as the example below. Here the Event Track represents the Alternative Splicing Event coordinates given in the aggregation file (`params.visualisation.aggregation_csv)`, and we can see the structure and mappings of the `ENST00000375234.7` and `ENST00000411939.5` transcripts. Here the `params.visualisation.show_non_overlapping` parameter has been set to TRUE.

![](test_output/PRXL2C_ENSG00000158122.11%3BSE:chr9:96641886-96645893:96646024-96651390:-_Gviz-Trackplot.png)

### Pipeline diagram

This is the flow chart describind the different processes of the pipeline:

**UPDATE**

![](new_flowchart.png)
