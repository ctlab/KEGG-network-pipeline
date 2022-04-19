import pandas as pd


# Leaving only C atoms
atom_mapping = pd.read_csv("rdt_analysis/atom_mapping.csv", sep=',')

atom_mapping = atom_mapping[["reaction", "reactant", "num", "product"]]

cond = atom_mapping[atom_mapping['reactant'].str.contains('_C_') & atom_mapping['product'].str.contains('_C_')]
cond2 = cond[cond['reactant'].isnull() | cond['product'].isnull()]

# print("Amount of NA in C columns is:", len(cond2), "out of", len(cond),
#       "which is", round(100 * len(cond2) / len(cond), 2), "%.")

cond.to_csv(r'data/atom_mapping_C_atoms.csv', index=False)

# Making RPAIRS
atom_mapping = pd.read_csv("data/atom_mapping_C_atoms.csv", sep=",")
atom_mapping = atom_mapping.dropna()

tmp = atom_mapping["reactant"].str.split("_", n=2, expand=True)
tmp["atom.x"] = tmp[[0, 2]].agg("_".join, axis=1)
atom_mapping["metabolite.x"] = tmp[0]
atom_mapping["element.x"] = tmp[1]
atom_mapping["atom.x"] = tmp["atom.x"]

tmp = atom_mapping["product"].str.split("_", n=2, expand=True)
tmp["atom.y"] = tmp[[0, 2]].agg("_".join, axis=1)
atom_mapping["metabolite.y"] = tmp[0]
atom_mapping["element.y"] = tmp[1]
atom_mapping["atom.y"] = tmp["atom.y"]

atom_mapping = atom_mapping[["reaction", "metabolite.x", "element.x", "atom.x", "metabolite.y", "element.y", "atom.y"]]
atom_mapping = atom_mapping.drop_duplicates()
atom_mapping = atom_mapping.reset_index(drop=True)

full_atom_mapping = atom_mapping
atom_mapping = atom_mapping[["reaction", "metabolite.x", "metabolite.y"]]

atom_mapping = atom_mapping.drop_duplicates()
atom_mapping = atom_mapping.reset_index(drop=True)

atom_mapping["index"] = atom_mapping.index
atom_mapping["index"] = atom_mapping["index"] + 1
atom_mapping["index"] = atom_mapping["index"].astype(str)

atom_mapping["length"] = 5 - atom_mapping["index"].str.len()

atom_mapping["rpair"] = atom_mapping.length.apply(lambda x: ("RP" + "0" * x))
atom_mapping["rpair"] = atom_mapping["rpair"] + atom_mapping["index"]

atom_mapping = atom_mapping[["reaction", "metabolite.x", "metabolite.y", "rpair"]]
# atom_mapping["rtype"] = "-"

# atom_mapping.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]
atom_mapping.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair"]

atom_mapping = pd.DataFrame(atom_mapping)
atom_mapping.to_csv(r"data/rpairs_preprocessing.csv", index=False)

atom_mapping.columns = ["reaction", "metabolite.x", "metabolite.y", "rpair"]
full_atom_mapping_m = pd.merge(full_atom_mapping, atom_mapping, how="left",
                               left_on=["reaction", "metabolite.x", "metabolite.y"],
                               right_on=["reaction", "metabolite.x", "metabolite.y"])

# rpair2align = full_atom_mapping_m[["rpair", "atom.x", "atom.y"]]
# rpair2align = rpair2align.drop_duplicates()
# rpair2align = rpair2align.reset_index(drop=True)
# rpair2align.to_csv(r"network_files/rpair2align.csv", index=False)

reaction2align = full_atom_mapping_m[["reaction", "atom.x", "atom.y"]]
reaction2align = reaction2align.drop_duplicates()
reaction2align = reaction2align.reset_index(drop=True)
reaction2align.to_csv(r"network_files/reaction2align.csv", index=False)

# reaction2pair = full_atom_mapping_m[["rpair", "reaction"]]
# reaction2pair = reaction2pair.drop_duplicates()
# reaction2pair = reaction2pair.reset_index(drop=True)
# reaction2pair.to_csv(r"network_files/reaction2pair.csv", index=False)


atoms_x = full_atom_mapping_m[["atom.x", "metabolite.x", "element.x"]]
atoms_x.columns = ["atom", "metabolite", "element"]
atoms_y = full_atom_mapping_m[["atom.y", "metabolite.y", "element.y"]]
atoms_y.columns = ["atom", "metabolite", "element"]

atoms = pd.concat([atoms_x, atoms_y])
atoms = atoms.drop_duplicates()
atoms = atoms.reset_index(drop=True)
atoms.to_csv(r"network_files/atoms.csv", index=False)


metabolite2atom_x = full_atom_mapping_m[["atom.x", "metabolite.x"]]
metabolite2atom_x.columns = ["atom", "metabolite"]
metabolite2atom_y = full_atom_mapping_m[["atom.y", "metabolite.y"]]
metabolite2atom_y.columns = ["atom", "metabolite"]

metabolite2atom = pd.concat([metabolite2atom_x, metabolite2atom_y])
metabolite2atom = metabolite2atom.drop_duplicates()
metabolite2atom = metabolite2atom.reset_index(drop=True)
metabolite2atom.to_csv(r"network_files/metabolite2atom.csv", index=False)


## RTYPE
full_atom_mapping_m["el.check.x"] = 1
full_atom_mapping_m["el.check.y"] = 1

full_atom_mapping_m.loc[full_atom_mapping_m["element.y"] != "C", "el.check.y"] = 0
# full_atom_mapping_m[full_atom_mapping_m["element.y"] != "C"]

full_atom_mapping_2 = full_atom_mapping_m[["reaction", "metabolite.x", "metabolite.y", "rpair", "el.check.x", "el.check.y"]]

full_atom_mapping_3 = full_atom_mapping_2.groupby(["reaction", "metabolite.x", "metabolite.y", "rpair"],
                                                  as_index=False).agg({"el.check.x": "sum",
                                                                       "el.check.y": "sum"})

# full_atom_mapping_3[full_atom_mapping_3["el.check.x"] < 2 ]
full_atom_mapping_3.loc[full_atom_mapping_3["el.check.x"] >= 2, "check"] = 1
full_atom_mapping_3.loc[full_atom_mapping_3["el.check.y"] >= 2, "check"] = 1
full_atom_mapping_3.loc[full_atom_mapping_3["el.check.x"] < 2, "check"] = 0
full_atom_mapping_3.loc[full_atom_mapping_3["el.check.y"] < 2, "check"] = 0
# full_atom_mapping_3[full_atom_mapping_3["check"] != 1]

full_atom_mapping_3["rtype"] = "not_main"
full_atom_mapping_3.loc[full_atom_mapping_3["check"] == 1, "rtype"] = "main"

full_atom_mapping_3 = full_atom_mapping_3[["reaction", "metabolite.x", "metabolite.y", "rpair", "rtype"]]
full_atom_mapping_3.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair", "rtype"]

atom_mapping.columns = ["rxn", "metabolite.x", "metabolite.y", "rpair"]

merged = pd.merge(atom_mapping, full_atom_mapping_3, how="left",
                  left_on=["rxn", "metabolite.x", "metabolite.y", "rpair"],
                  right_on=["rxn", "metabolite.x", "metabolite.y", "rpair"])

merged.to_csv(r"network_files/rpairs_2C.csv", index=False)
