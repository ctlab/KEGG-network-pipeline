library(KEGGREST)
library(data.table)

reaction2enzyme <- keggLink("enzyme", "reaction")
reaction2enzyme <- data.table(
  reaction=gsub("rn:", "", names(reaction2enzyme), fixed = T),
  enzyme=gsub("ec:", "", reaction2enzyme, fixed = T))

reactions <- keggList("reaction")
reactions <- data.table(
  reaction=gsub("rn:", "", names(reactions), fixed = T),
  reaction_name=gsub(";.*$", "", reactions),
  reaction_equation=gsub("^.*; ", "", reactions))

print(head(reactions))
reactions[ , reaction_url := sprintf("http://www.kegg.jp/entry/%s", reaction)]

enzyme2reaction <- copy(reaction2enzyme)
setkey(enzyme2reaction, enzyme)

network <- list()
network$reactions <- reactions
network$enzyme2reaction <- enzyme2reaction

network$reaction2align <- fread("network_files/reaction2align.csv")
setkey(network$reaction2align, reaction)
network$atoms <- data.table(atom=union(network$reaction2align$atom.x,
                                           network$reaction2align$atom.y))

network$atoms <- fread("network_files/atoms.csv")
setkey(network$atoms, atom)

network$metabolite2atom <- network$atoms[, list(metabolite, atom)]
setkey(network$metabolite2atom, metabolite)

save(network, file="network/network.rda")

print("Pipeline finished")