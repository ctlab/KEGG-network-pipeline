library(devtools)
library("ChemmineR")

RDT_outputs <- list.dirs("rdt_output")
RDT_outputs <- RDT_outputs[-1]

final_mapping <- dplyr::tibble()

ncolnot15 <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(ncolnot15) <- c("react_id", "number_of_not")

withHydros <- data.frame(matrix(ncol = 2, nrow = 0))
colnames(withHydros) <- c("react_id", "number_of_pos")

hydro_true <- FALSE
problem_list <- list()

for(dir in RDT_outputs){
    #print(dir)
    rxn <- list.files(dir, pattern = ".rxn$", full.names = TRUE)
    #print(rxn)
    rdt_lines <- readr::read_lines(rxn)
    rdt <- ChemmineR::read.SDFindex(rxn,
                                    index = data.frame(
                                        "A" = which(stringr::str_detect(rdt_lines, "[CMG][:digit:]{5}")), 
                                        "B" = which(stringr::str_detect(rdt_lines, "END"))
                                    ))

    react_id <- stringr::str_extract(rdt_lines[3], "R[:digit:]{5}")
    # print(react_id)
    if(is.na(react_id)){
        print("rxn file is empty:")
        print(dir)
    } else {
        if(any(sapply(ChemmineR::atomblock(rdt), dim)[2, ] == 2)){
            hydro_pos <- which(sapply(ChemmineR::atomblock(rdt), dim)[2, ] == 2)
            hydro_true <- TRUE
            withHydros <- rbind(withHydros, 
                            data.frame("react_id" = react_id, "number_of_pos" = length(hydro_pos)))
            rdt <- ChemmineR::read.SDFindex(rxn,
                                            index = data.frame(
                                                "A" = which(stringr::str_detect(rdt_lines, "[CMG][:digit:]{5}"))[-hydro_pos],
                                                "B" = which(stringr::str_detect(rdt_lines, "END"))[-hydro_pos]
                                            ))
        }

        if(any(sapply(ChemmineR::atomblock(rdt), dim)[2, ] != 15)){
            ncolnot15 <- rbind(ncolnot15, 
                               data.frame("react_id" = react_id, 
                                          "number_of_not" = length(which(sapply(ChemmineR::atomblock(rdt), dim)[2, ] != 15))))
            hydro_true <- FALSE
        } else {
            nn <- stringr::str_match_all(rdt_lines[5], "[:digit:]+"); nn 
            n_reactants <- nn[[1]][1]
            list_of_reactants <- data.frame()
            list_of_products <- data.frame()

            if(hydro_true){
                hydro_true <- FALSE
                pos_reactants <- seq(as.numeric(n_reactants))
                pos_to_remove <- c()
                for(pos in hydro_pos){

                    if(pos <= n_reactants){
                        pos_to_remove <- c(pos_to_remove, pos)
                    }
                }

                if(length(pos_to_remove) > 0){
                    pos_reactants <- pos_reactants[-pos_to_remove]
                    n_reactants <- length(pos_reactants)
                }
            }

            if(n_reactants == 0 | n_reactants == length(rdt)){
                withHydros <- rbind(withHydros, 
                                    data.frame("react_id" = react_id, "number_of_pos" = "turned zero"))
                next
            }
        
            for(el in seq_along(rdt)){
                curr <- ChemmineR::atomblock(rdt)[[el]]
                doing_cid <- ChemmineR::sdfid(rdt)[[el]]
                
                doing_aid <- paste(doing_cid,
                                   gsub("_.*", "", rownames(curr)), curr[ , 1], curr[ , 2], 
                                   sep = "_")
                
                if(el <= n_reactants){
                    list_of_reactants <- rbind(list_of_reactants,
                                               data.frame("reaction" = react_id,
                                                          "reactant" = doing_aid,
                                                          "num" = unname(curr[ , 13])))
                } else {
                    list_of_products <- rbind(list_of_products,
                                              data.frame("reaction" = react_id,
                                                         "product" = doing_aid,
                                                         "num" = unname(curr[ , 13])))
                }
            }

            if(length(list_of_reactants) != 0 && length(list_of_products) != 0){
                final_mapping <- dplyr::bind_rows(final_mapping,
                                                  dplyr::full_join(list_of_reactants, 
                                                                   list_of_products, 
                                                                   by = c("reaction", "num")))
            } else {
                print("smth went wrong")
	        print(react_id)
                problem_list[[length(problem_list) + 1]] <- react_id    # Append new list element
            }
        }
    }
}

write.csv(final_mapping, file = "rdt_analysis/atom_mapping.csv")
save(final_mapping, file = "rdt_analysis/atom_mapping.Rda")
