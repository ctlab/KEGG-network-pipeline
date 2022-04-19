import pandas as pd


data = pd.read_csv("data/HMDB2metabolite.csv", sep=",")

old = data[data['HMDB'].str.contains('HMDB00')]
old['HMDB'] = old['HMDB'].str.replace('HMDB00', 'HMDB')

mapping = pd.concat([data, old])
mapping.to_csv(r'data/HMDB2metabolite_old_ids.csv', index=False)


data = pd.read_csv("data/kegg2chebi.tsv", sep="\t", header=None)
data.columns = ["KEGG", "ChEBI"]

data['KEGG'] = data['KEGG'].str.replace('cpd:', '')
data['ChEBI'] = data['ChEBI'].str.replace('chebi:', 'CHEBI:')

data.to_csv(r'data/kegg2chebi_upd.csv', index=False)
