import os


rule all:
    input:
        "network/network.rda",
        "network/met.kegg.db.rda"


# Downloading reactions, preprocessing them, downloading compounds
rule getting_reactions_from_KEGG:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "1_downloading_reactions_list.py"
    output:
        "data/kegg_reactions.txt"
    message:
        "...Downloading reactions' list..."
    shell:
        "python3.9 1_downloading_reactions_list.py"


rule downloading_KEGG_reactions_content:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg_reactions.txt",

        "2_creating_reactions_table.py"
    output:
        "data/kegg_reactions.csv"
    message:
        "...Downloading reactions in EQUATION format for KEGG..."
    shell:
        "python3.9 2_creating_reactions_table.py"


rule replacing_coefficients:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg_reactions.csv",

        "3_Replacing_coefficients.py"
    output:
        "data/kegg_reactions_redone.txt",
        "data/kegg_reactions_redone.csv"
    message:
        "...Replacing coefficients for KEGG reactions..."
    shell:
        "python3.9 3_Replacing_coefficients.py"


rule get_compounds_MOL_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg_reactions_redone.csv",

        "4_downloading_KEGG_MOL_files.py"
    output:
        "data/no_mol_files.txt"
    message:
        "...Downloading MOL-files for KEGG..."
    shell:
        "mkdir -p data/MOL_files;"
        "python3.9 4_downloading_KEGG_MOL_files.py"


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


checkpoint creating_RXN_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/rdt_test_done.txt",

        "5_Make_RXN_files.R"
    output:
        reactions=directory("data/rxn")
    message:
        "...Creating RXN files from MOL-files..."
    shell:
        "mkdir -p data/rxn;"
        "Rscript 5_Make_RXN_files.R"


## Running RDT
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


## Analyse RDT
rule performing_RDT_results_analysis:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/rdt_done.txt",

        "7_Rdt_Output_Analysis.R"
    output:
        "rdt_analysis/atom_mapping.csv",
        "rdt_analysis/atom_mapping.Rda"
    message:
        "...Analysing RDT results and creating atom mapping table..."
    shell:
        "mkdir -p rdt_analysis;"
        "Rscript 7_Rdt_Output_Analysis.R"


## Making RPAIRS and supplementary files
rule creating_RPAIRS_and_supplementary_files:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "rdt_analysis/atom_mapping.csv",
        "8_Making_RPAIRS_and_supplementary_files.py"
    output:
        "data/atom_mapping_C_atoms.csv",
        "data/rpairs_preprocessing.csv",
        "network_files/reaction2align.csv",
        "network_files/atoms.csv",
        "network_files/metabolite2atom.csv",
        "network_files/rpairs_2C.csv"
    message:
        "...Creating RPAIRS and supplementary files..."
    shell:
        "mkdir -p network_files;"
        "mkdir -p network;"
        "python3.9 8_Making_RPAIRS_and_supplementary_files.py"


## Creating met.kegg.db
rule supp_files_for_met_db_kegg:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "network_files/atoms.csv",

        "9_Creating_met_db_supplementary_files.py"
    output:
        "data/kegg2chebi.tsv"
    message:
        "...9. Creating met.db supplementary files..."
    shell:
        "wget http://rest.kegg.jp/conv/chebi/compound -O data/kegg2chebi.tsv"


rule getting_mapping_tables_for_met_db_kegg:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg2chebi.tsv",

        "9_1_Getting_HMDB_mappings.R"
    output:
        "data/HMDB2metabolite.csv"
    message:
        "...9. Creating met.db supplementary files..."
    shell:
        "Rscript 9_1_Getting_HMDB_mappings.R"


rule converting_mappings_for_met_db_kegg:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/kegg2chebi.tsv",
        "data/HMDB2metabolite.csv",

        "10_2_Converting_db_IDs.py"
    output:
        "data/HMDB2metabolite_old_ids.csv",
        "data/kegg2chebi_upd.csv"
    message:
        "...9. Creating met.db supplementary files..."
    shell:
        "python3.9 9_2_Converting_db_IDs.py"


rule creating_met_db_rhea_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "data/HMDB2metabolite_old_ids.csv",
        "data/kegg2chebi_upd.csv",

        "9_3_creating_met.kegg.db.R"
    output:
        "network/met.kegg.db.rda"
    message:
        "...Creating met.kegg.db..."
    shell:
        "Rscript 9_3_creating_met.kegg.db.R"


## Creating network object
rule creating_network_object:
    singularity:
        "docker://mariaembio/python3.9_r4.0_java11"
    input:
        "network_files/rpairs_2C.csv",
        "network_files/reaction2align.csv",
        "network_files/atoms.csv",

        "10_Creating_network_object.R"
    output:
        "network/network.rda"
    message:
        "...Creating network object..."
    shell:
        "Rscript 10_Creating_network_object.R"
