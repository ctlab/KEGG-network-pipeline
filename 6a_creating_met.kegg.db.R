library(KEGGREST)
library(data.table)


# metabolites
metabolites <- keggList("compound")
metabolites <- data.table(
  metabolite=gsub("cpd:", "", names(metabolites), fixed = T),
  metabolite_name=gsub(";.*$", "", metabolites))
metabolites[, metabolite_url := sprintf("http://www.kegg.jp/entry/%s", metabolite)]
setkey(metabolites, metabolite)


# mapFrom
HMDB2metabolite <- fread("data/HMDB2metabolite_old_ids.csv")
HMDB2metabolite <- as.data.table(HMDB2metabolite)
HMDB2metabolite <- HMDB2metabolite[, list(HMDB=HMDB, metabolite=KEGG)]
setkey(HMDB2metabolite, HMDB)

ChEBI2metabolite <-fread("data/kegg2chebi_upd.csv")
ChEBI2metabolite <- as.data.frame(ChEBI2metabolite)
ChEBI2metabolite <- unique(ChEBI2metabolite)
ChEBI2metabolite <- as.data.table(ChEBI2metabolite)

ChEBI2metabolite <- ChEBI2metabolite[, list(ChEBI=ChEBI, metabolite=KEGG)]
setkey(ChEBI2metabolite, ChEBI)


# anomers
anomers <- fread("pre_data/mets2collapse.tsv")[, list(metabolite=from, base_metabolite=to)]
anomers <- unique(anomers)

met.kegg.db <- list()
met.kegg.db$metabolites <- metabolites
met.kegg.db$mapFrom <- list("HMDB"=HMDB2metabolite, "ChEBI"=ChEBI2metabolite)
met.kegg.db$baseId <- "KEGG"
met.kegg.db$anomers <- list(
  metabolite2base_metabolite=copy(anomers),
  base_metabolite2metabolite=copy(anomers))
setkey(met.kegg.db$anomers$metabolite2base_metabolite, metabolite)
setkey(met.kegg.db$anomers$base_metabolite2metabolite, base_metabolite)


save(met.kegg.db, file="network/met.kegg.db.rda")