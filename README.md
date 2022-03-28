## Pipeline for the creation of KEGG metabolic network

This pipeline produces network-objects and corresponding metabolite annotations based on
KEGG database (https://www.genome.jp/kegg) which are used for analysis by Shiny GATOM
(https://artyomovlab.wustl.edu/shiny/gatom)

## Docker

The pipeline is designed to work in specifically constructed docker container which can
be found here: https://hub.docker.com/repository/docker/mariaembio/python3.9_r4.0_java11.

## Execution of the pipeline 

The code is written in Python and R, and designed to be executed with Snakemake and 
singularity. Execution will need at least 2 GB of empty space and can take from 2 to 10 
days to run.

## Data

Pipeline takes as input manually created anomers list for KEGG database.

## Pipeline description

1. List of KEGG reactions is downloaded via KEGG Rest;
2. Reactions in EQUATION format are downloaded via KEGG Rest and then coefficients 
that contain n, m or x are replaced;
3. MOL-files for compounds of the reactions are downloaded and then RXN files for 
reactions are compiled from those MOL-fies;
4. Then Reaction Decoder Tool (https://github.com/asad/ReactionDecoder) is used for atom 
mapping of the RXN files;
5. Atom mapping tables are created with the use of ChemmineR
(https://www.bioconductor.org/packages/release/bioc/html/ChemmineR.html);
6. Various supplementary files networks are created, including
KEGG reactions' and metabolites' annotation;
7. Network-object and corresponding metabolites' annotation object for KEGG network 
is created .

