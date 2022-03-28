library(data.table)
library(metaboliteIDmapping)

HMDB2metabolite <- metabolitesMapping[(!is.na(metabolitesMapping$KEGG)),
                                      c("KEGG", "HMDB")]
HMDB2metabolite <- na.omit(HMDB2metabolite)
HMDB2metabolite <- as.data.frame(HMDB2metabolite)
HMDB2metabolite <- unique(HMDB2metabolite)
HMDB2metabolite <- as.data.table(HMDB2metabolite)
HMDB2metabolite$ChEBI <- paste0("CHEBI:", HMDB2metabolite$ChEBI)

write.csv(HMDB2metabolite, "data/HMDB2metabolite.csv", row.names = F)
