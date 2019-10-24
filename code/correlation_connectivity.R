# main function ks score
correlation_score <- function(up_signature, down_signature, rank_matrix, compound_name) {
    n = 1
    connectivity_score_vec <- c()
    gene_rank_bing <- rank_matrix
    for (i in 1:ncol(gene_rank_bing)) {
        if (i %% 500 == 0)  { print(paste0("INFO: ", i, " Instances."))}
        pert_vec <- gene_rank_bing[, i]
        up_signature <- up_signature[up_signature %in% pert_vec]   
        down_signature <- down_signature[down_signature %in% pert_vec]
        
        # get teh rank
        up_v <- match(up_signature, pert_vec)
        down_v <- match(down_signature, pert_vec)
        
        up_score <- cor((1:length(up_v))/length(up_v), up_v/nrow(gene_rank_bing))
        down_score <- cor((1:length(down_v))/length(down_v), down_v/nrow(gene_rank_bing))
        tmp_s <- up_score + down_score
        connectivity_score_vec <- c(connectivity_score_vec, tmp_s)   
    
        
    }
    output_df <- data.frame(compound_name = compound_name, connectivity_score = connectivity_score_vec)
    return(output_df)
    
}



correlation_score_merged <- function(up_signature, down_signature, rank_matrix, compound_name) {
    n = 1
    connectivity_score_vec <- c()
    gene_rank_bing <- rank_matrix
    for (i in 1:ncol(gene_rank_bing)) {
        if (i %% 500 == 0)  { print(paste0("INFO: ", i, " Instances."))}
        pert_vec <- gene_rank_bing[, i]
        up_signature <- up_signature[up_signature %in% pert_vec]   
        down_signature <- down_signature[down_signature %in% pert_vec]
        # get teh rank
        merged_signature <- c(up_signature, down_signature)
        rank_v <- match(merged_signature, pert_vec)
        tmp_s <- cor((1:length(merged_signature))/length(merged_signature), rank_v/nrow(gene_rank_bing))
        connectivity_score_vec <- c(connectivity_score_vec, tmp_s) 
        
    }
    output_df <- data.frame(compound_name = compound_name, connectivity_score = connectivity_score_vec)
    return(output_df)
    
}


# data loading 
#correlation_score_df <- correlation_score(up_signature_df$entrezgene_id, down_signature_df$entrezgene_id, gene_rank_bing, compound_name)
#correlation_score_merged_df <- correlation_score(up_signature_df$entrezgene_id, down_signature_df$entrezgene_id, gene_rank_bing, compound_name)

#correlation_score_df <- correlation_score_df %>% 
#    arrange(desc(connectivity_score))
#correlation_score_merged_df <- correlation_score_merged_df %>% 
#    arrange(desc(connectivity_score))




