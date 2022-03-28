library(tidyverse)

lines <- readr::read_lines("data/kegg_reactions_redone.txt")
equations <- lines %>%
    str_split_fixed("\t", 2) %>% {
        setNames(.[,1], .[,2])
    }

readr::write_rds(equations, "data/equations.rds")

eqs <- readRDS("data/equations.rds")
pattern <- "([:digit:]*)\\s*([A-za-z][:digit:]{5})"
pattern_mult <- "([:digit:]+)\\s*(C[:digit:]{5})"
reaction_error <- c()

reaction <- names(eqs[2])
reaction_id <- unname(eqs[2])

for(rr in seq_along(eqs)){
    reaction <- names(eqs[rr])
    reaction_id <- unname(eqs[rr])

    # ------------ dealing with multiple compounds -----------------
    m <- stringr::str_extract_all(reaction, pattern) #; m
    # checking <- stringr::str_detect(m, "([:digit:]+)\\s*(C[:digit:]{5})")
    # checking <- stringr::str_detect(unlist(m), "([:digit:]+)\\s*(C[:digit:]{5})")
    need_to_unmult <- stringr::str_extract(unlist(m), pattern_mult) #; need_to_mult
    
    if(stringr::str_detect(m, pattern)){
        if(any(!is.na(need_to_unmult))){
            for(not_na in need_to_unmult[!is.na(need_to_unmult)]){
                get_elemens <- stringr::str_match(not_na, pattern_mult); get_elemens
                reaction <- gsub(not_na,
                                 paste(rep(get_elemens[3], get_elemens[2]), collapse = " + "),
                                 reaction)
            }
        }
    }

    # ------ extract compound IDs -------
    comp_ids <- unlist(stringr::str_extract_all(reaction, "[A-za-z][:digit:]{5}")) #; comp_ids
    
    if(all(file.exists(
        paste0("data/MOL_files/", comp_ids, ".txt")))){
        
        # ------ detect how many reactants and products -------
        arrow <- unlist(strsplit(reaction, "<=>"))
        n <- stringr::str_count(arrow[1], pattern = "C[:digit:]{5}")
        m <- stringr::str_count(arrow[2], pattern = "C[:digit:]{5}")
        
        # ------ translate reaction in RXN format and save -----
        content_of_output <- c("$RXN", "", "  WHATEVER     blabla", "",
                               paste0("  ", n, "  ", m))
        
        for(comp in comp_ids){

            # try to extract all mol files, if not --> next reaction
            comp_mol <- readr::read_lines(
                paste0("data/MOL_files/", comp, ".txt"))
            
            where_stop <- which(grepl("> <ENTRY>", comp_mol)) - 1
            id <- gsub("cpd:", "", comp_mol[grepl("cpd:C", comp_mol)])
            
            content_of_output <- c(content_of_output,
                                   "$MOL", 
                                   id, "  WHATEVER  0000000000", 
                                   comp_mol[-c(1, 2, where_stop:length(comp_mol))])
        }
        # ------ save file ------ 
        readr::write_lines(content_of_output, 
                           path = paste0("data/rxn/", 
                                         reaction_id, ".rxn"))  }
    else {
        reaction_error <<- c(reaction_error, reaction_id)
    }
}


readr::write_lines(reaction_error, path = paste0("data/reactions_KEGG.no_rxn.txt"))
