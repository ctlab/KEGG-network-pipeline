# Pipeline for the creation of KEGG metabolic network

This pipeline produces network-objects and corresponding metabolites annotation based on
[KEGG database](https://www.genome.jp/kegg) which are used for analysis by 
[Shiny GATOM](https://artyomovlab.wustl.edu/shiny/gatom). The pipeline uses Snakemake and Singularity 
for execution and uses the provided scripts while running.


## Requirements

* Snakemake 
* Singularity v3.5.3
* At least 3 GB of disk space 


## Quick Start

* Clone the repository:

```git clone git@github.com:ctlab/KEGG-network-pipeline.git```

* Activate the environment with snakemake and singularity:

```conda activate snakemake```

* Execute the pipeline:

```snakemake --use-singularity --cores 8```

Execution will need at least 3 GB of empty space and can take from 2 to 10 
days to run depending on the setup.


## Pipeline structure

The pipeline is designed to work in specifically constructed docker image which can
be found [here](https://hub.docker.com/repository/docker/mariaembio/python3.9_r4.0_java11).
It contains all necessary dependencies for Python and R code as well as Reaction Decoder Tool,
and Snakemake will use it automatically when executing.

Some steps of the pipeline will execute included scripts that are written in Python or R, 
while simpler steps will only run bash commands from Snakefile. 

Pipeline takes as input `mets2collapse.tsv` file that are stored in `pre_data` folder. This file 
contains manually created list of anomers for the KEGG database.

Network `rds` object and metabolites annotation `rds` object are considered to be the 
output of the pipeline. The files will be stored in `network` folder.


## Snakemake pipeline steps

1. List of KEGG reactions is downloaded;
2. In order to construct atom network, we need to perform atom mapping with Reaction Decoder Tool.
The Reaction Decoder Tool takes as input RXN files, however, the KEGG database does not provide them. 
To compile RXN files the following steps are done:

    a. A table of reactions in EQUATION format is downloaded via KEGG Rest;

    b. Coefficients for reactions in EQUATION format that contain n, m or x are replaced with numeric values;

    c. MOL files for compounds of the reactions are downloaded;

    d. RXN files for reactions are compiled from those MOL files with the use of the table of reactions in EQUATION format;
3. The [Reaction Decoder Tool](https://github.com/asad/ReactionDecoder) is used for atom mapping of the RXN files;
4. Atom mapping tables are created with the use of 
[ChemmineR](https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html)
from RXN files processed with the Reaction Decoder Tool;
5. The supplementary annotation files for network and metabolites are created:
    - Network annotation includes reactions' KEGG ID & link obtained from KEGG, reaction-enzyme mapping 
   obtained from KEGG & processed atom mapping data.
    - Metabolites annotation includes metabolites' KEGG ID & link extracted from KEGG, HMDB to KEGG mapping 
   obtained via metaboliteIDmapping R-package, ChEBI to KEGG mapping extracted from KEGG & anomers list.
6. The network & metabolites annotation objects are created.

    a. Metabolites annotation object for KEGG network is created;

    b. Network-object for KEGG network is created.

