## Libraries
import pandas as pd
import re
import requests
import time


kegg_with_enzymes = pd.read_csv('data/kegg_reactions.csv',
                                sep=';')
final_table = []
compounds = []

for i in range(len(kegg_with_enzymes)):
    reactants = re.findall(r'([CG]\d+)', str(kegg_with_enzymes.Reaction[i]))
    for reactant in reactants:
        if reactant not in compounds:
            compounds.append(reactant)

mol_files_counter = 0
no_mol_files_counter = 0
no_mol_files_list = []

count = 1

for compound_id in compounds:
    ## Example:
    url = "http://rest.kegg.jp/get/" + compound_id + "/mol"  # http://rest.kegg.jp/get/C00001/mol

    try:
        r = requests.post(url)
    except ConnectionResetError:
        time.sleep(3)
        r = requests.post(url)
    except URLError:
        time.sleep(6)
        r = requests.post(url)
    except:
        time.sleep(10)
        r = requests.post(url)

    if len(r.text) > 0:
        mol_files_counter += 1  # 'x' means that it creates a new file
        with open("data/MOL_files/" + str(compound_id) + ".txt", "x", encoding="utf8") as file:
            file.write(str(r.text))

    elif len(r.text) == 0:
        no_mol_files_counter += 1
        no_mol_files_list.append(compound_id)

    count += 1

print("\n\nTotal amount of compounds with MOL-files: ", mol_files_counter)
print("Total amount of compounds without MOL-files: ", no_mol_files_counter)

with open("data/no_mol_files.txt", "x") as file:
    print(*no_mol_files_list, file=file, sep="\n")
