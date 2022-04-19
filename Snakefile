import os


rule all:
    input:
        "network/network.rda",
        "network/met.kegg.db.rda"


## 1. Downloading list of KEGG reactions
rule getting_reactions_from_KEGG:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "1_downloading_reactions_list.py"
    output:
        "data/kegg_reactions.txt"
    message:
        "...Downloading reactions list..."
    shell:
        "python3.9 1_downloading_reactions_list.py"


## 2. Compiling RXN files
# 2a. Downloading a table of reactions in EQUATION format
rule downloading_KEGG_reactions_content:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg_reactions.txt",

        "2a_downloading_EQUATION_reactions_table.py"
    output:
        "data/kegg_reactions.csv"
    message:
        "...Downloading reactions in EQUATION format for KEGG..."
    shell:
        "python3.9 2a_downloading_EQUATION_reactions_table.py"


# 2b. Coefficients for reactions in EQUATION format that contain n, m or x are replaced with numeric values
rule replacing_coefficients:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg_reactions.csv",

        "2b_replacing_coefficients.py"
    output:
        "data/kegg_reactions_redone.txt",
        "data/kegg_reactions_redone.csv"
    message:
        "...Replacing coefficients for KEGG reactions..."
    shell:
        "python3.9 2b_replacing_coefficients.py"


# 2c. MOL files for compounds of the reactions are downloaded
rule get_compounds_MOL_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg_reactions_redone.csv",

        "2c_downloading_KEGG_MOL_files.py"
    output:
        "data/no_mol_files.txt"
    message:
        "...Downloading MOL files for KEGG..."
    shell:
        "mkdir -p data/MOL_files;"
        "python3.9 2c_downloading_KEGG_MOL_files.py"


## Extra step for proper RDT execution
rule pre_RDT:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/no_mol_files.txt"
    output:
        "data/rdt_test_done.txt"
    message:
        "...Preparing for RDT step..."
    shell:
        "cwd=$(pwd);"
        "mkdir -p data/pre_rdt_output;"
        "cd data/pre_rdt_output;"
        "java -jar /rdt-2.4.1-jar-with-dependencies.jar -Q SMI -q 'CC(O)CC(=O)OC(C)CC(O)=O.O[H]>>[H]OC(=O)CC(C)O.CC(O)CC(O)=O' -g -c -j AAM -f TEXT || true;"
        "touch $cwd/{output}"


## 2d. Compiling RXN files from MOL files
checkpoint creating_RXN_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/rdt_test_done.txt",

        "2d_compiling_RXN_files.R"
    output:
        reactions=directory("data/rxn")
    message:
        "...Compiling RXN files from MOL files..."
    shell:
        "mkdir -p data/rxn;"
        "Rscript 2d_compiling_RXN_files.R"


## 3. Running RDT
rule performing_RDT:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        rxn = "data/rxn/{reaction}.rxn"
    output:
        "rdt_output/{reaction}/ECBLAST_{reaction}_AAM.rxn"
    message:
        "...Performing RDT analysis..."
    shell:
        "cwd=$(pwd);"
        "mkdir -p rdt_output; mkdir -p failed_rdt;"
        "cd rdt_output/{wildcards.reaction};"
        "java -jar -Xmx8G /rdt-2.4.1-jar-with-dependencies.jar -Q RXN -q $cwd/{input.rxn} -g -j AAM -f TEXT || true;"
        "if [[ ! -f $cwd/rdt_output/{wildcards.reaction}/ECBLAST_{wildcards.reaction}_AAM.rxn ]]; then touch $cwd/failed_rdt/ECBLAST_{wildcards.reaction}_AAM.rxn; fi;"
        "touch $cwd/rdt_output/{wildcards.reaction}/ECBLAST_{wildcards.reaction}_AAM.rxn"


def aggregate_rxns(wildcards):
    checkpoint_output = checkpoints.creating_RXN_files.get(**wildcards).output[0]
    file_names = expand("rdt_output/{reaction}/ECBLAST_{reaction}_AAM.rxn",
                        reaction = glob_wildcards(os.path.join(checkpoint_output, "{reaction}.rxn")).reaction)
    return file_names


rule aggregate:
    input:
        aggregate_rxns
    output:
        "data/rdt_done.txt"
    shell:
        "touch {output};"
        'if [ "$(ls -A failed_rdt 2> /dev/null)" != "" ]; then cwd=$(pwd); cd rdt_output;'
        'for file in $cwd/failed_rdt/*; do filename=$(basename "$file");'
        'prefolder=$(basename "$file" "_AAM.rxn"); folder=${prefolder/ECBLAST_/};'
        'rm $folder/$filename; done;'
        'fi'


## 4. Creating atom mapping tables from RDT output
rule performing_RDT_results_analysis:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/rdt_done.txt",

        "4_RDT_output_analysis.R"
    output:
        "rdt_analysis/atom_mapping.csv",
        "rdt_analysis/atom_mapping.Rda"
    message:
        "...Creating atom mapping tables..."
    shell:
        "mkdir -p rdt_analysis;"
        "Rscript 4_RDT_output_analysis.R"


## 5. Creating supplementary files for network object & metabolites annotation object
# 5a. Creating supplementary files for network object
rule creating_network_supplementary_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "rdt_analysis/atom_mapping.csv",

        "5a_creating_network_supplementary_files.py"
    output:
        "data/atom_mapping_C_atoms.csv",
        "data/rpairs_preprocessing.csv",
        "network_files/reaction2align.csv",
        "network_files/atoms.csv",
        "network_files/metabolite2atom.csv",
        "network_files/rpairs_2C.csv"
    message:
        "...Creating supplementary files for network object..."
    shell:
        "mkdir -p network_files;"
        "mkdir -p network;"
        "python3.9 5a_creating_network_supplementary_files.py"


# 5b. Downloading KEGG to ChEBI mapping for metabolites annotation object
rule supp_files_for_met_db_kegg:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "network_files/atoms.csv"
    output:
        "data/kegg2chebi.tsv"
    message:
        "...Downloading supplementary files for metabolites annotation..."
    shell:
        "wget http://rest.kegg.jp/conv/chebi/compound -O {output}"


# 5c. Downloading KEGG to HMDB mapping for metabolites annotation object
rule getting_mapping_tables_for_met_db_kegg:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg2chebi.tsv",

        "5c_getting_HMDB_mappings.R"
    output:
        "data/HMDB2metabolite.csv"
    message:
        "...Creating supplementary files for metabolites annotation..."
    shell:
        "Rscript 5c_getting_HMDB_mappings.R"


# 5d. Downloading KEGG to HMDB mapping for metabolites annotation object
rule converting_mappings_for_met_db_kegg:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg2chebi.tsv",
        "data/HMDB2metabolite.csv",

        "5d_converting_db_IDs.py"
    output:
        "data/HMDB2metabolite_old_ids.csv",
        "data/kegg2chebi_upd.csv"
    message:
        "...Creating supplementary files for metabolites annotation..."
    shell:
        "python3.9 5d_converting_db_IDs.py"


## 6. Creating network & metabolites annotation objects
# 6a. Creating metabolites annotation object
rule creating_met_db_rhea_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/HMDB2metabolite_old_ids.csv",
        "data/kegg2chebi_upd.csv",

        "6a_creating_met.kegg.db.R"
    output:
        "network/met.kegg.db.rda"
    message:
        "...Creating metabolites annotation object..."
    shell:
        "Rscript 6a_creating_met.kegg.db.R"


# 6b. Creating network object
rule creating_network_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "network_files/rpairs_2C.csv",
        "network_files/reaction2align.csv",
        "network_files/atoms.csv",

        "6b_creating_network_object.R"
    output:
        "network/network.rda"
    message:
        "...Creating network object..."
    shell:
        "Rscript 6b_creating_network_object.R"
