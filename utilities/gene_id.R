# library
library(annotables)
library(dplyr)

# map the id
id_map <- function(id_vec, from_ref, to_ref, annotable = grch38) {
# specific operations
#     if (from_ref == "symbol" & to_ref = "ensembl") {
        
#     } else if (from_ref == "symbol" & to_ref = "ensembl") {
        
#     } else if (from_ref == "symbol" & to_ref = "entrez") {
        
#     } else if (from_ref == "ensembl" & to_ref = "symbol") {
        
#     } else if (from_ref == "ensembl" & to_ref = "entrez") {
        
#     } else if (from_ref == "entrez" & to_ref = "symbol") {
        
#     } else if (from_ref == "extrez" & to_ref = "ensembl") {
        
#     } 

    id_df <- data.frame(raw_id = id_vec)
    colnames(id_df) <- from_ref
    id_df[, from_ref] <- as.character(id_df[, from_ref])
    by_str <- rlang::set_names(quo_name(from_ref), quo_name(from_ref))
    output_df <- id_df %>%
        left_join(annotable, by = by_str) %>%
        select_(to_ref)
    return(output_df[, to_ref])                              
    
} 
                          
        
    
    
