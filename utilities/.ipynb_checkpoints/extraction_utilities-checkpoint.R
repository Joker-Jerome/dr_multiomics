# library
library(annotables)
library(dplyr)

# extract top idx
extract_top_idx <- function(score_df, top_n) {
    idx <- order(score_df$connectivity_score, decreasing = T)[1:top_n]
    #tmp_df <- score_df[order(score_df$connectivity_score, decreasing = T),][1:top_n,]
    return(idx)
    
}

# extract the drug gene that could reverse
gene_extract_topn <- function(up_signature, down_signature, gene_rank_bing, compound_name, topn) {
    pert_vec <- gene_rank_bing[, 1]
    up_signature <- up_signature[up_signature %in% pert_vec]   
    down_signature <- down_signature[down_signature %in% pert_vec]
    n_genes <- nrow(gene_rank_bing)
    p_list <- list()
    for (i in 1:ncol(gene_rank_bing)) {
        if (i %% 500 == 0)  { IRdisplay::display_html(paste0("INFO: ", i, " Instances."))}
        pert_vec <- gene_rank_bing[, i]
        
        # get the rank
        up_v <- match(up_signature, pert_vec)
        #up_geneset <- up_signature[up_v[na.omit(which(up_v < round(n_genes * topn)))]]
        up_geneset <- up_signature[na.omit(which(up_v < round(n_genes * topn)))]
        down_v <- match(down_signature, pert_vec)
        #down_geneset <- down_signature[down_v[na.omit(which(down_v > round(n_genes * (1-topn))))]]
        down_geneset <- down_signature[na.omit(which(down_v > round(n_genes * (1-topn))))]
        p_list[[i]] <- list(up_geneset = up_geneset, down_geneset = down_geneset)
        
    }
    return(p_list)
    
}

# pair plots
pair_plot_gene <- function(tmp_df, title_str) {
    p <- ggpaired(tmp_df, x = "group", y = "rank_value",
         color = "group", line.color = "gray", line.size = 0.4,
         palette = "jco", title = title_str) 
    #stat_compare_means(paired = TRUE)
    p
}

# visualize plots
visualize_enrichment_gene <- function(up_signature, down_signature, gene_rank_matrix, compound_name, idx_vec) {
    # sample a vector
    pert_vec <- gene_rank_matrix[, 1]
    up_signature <- up_signature[up_signature %in% pert_vec]   
    down_signature <- down_signature[down_signature %in% pert_vec]
    p_list <- list()
    # get the matching results for all the drugs in idx_vec
    for (i in 1:ncol(gene_rank_matrix[,idx_vec])) {
        if (i %% 1000 == 0)  { IRdisplay::display_html(paste0("INFO: ", i, " Instances."))}
        pert_vec <- gene_rank_bing[, i]
        # get the rank
        up_v <- match(up_signature, pert_vec)
        down_v <- match(down_signature, pert_vec)
        #print(length(up_v))
        #print(length(down_v))
        plot_df_up <- data.frame(rank_value = c(1-1:length(up_signature)/length(up_signature),1-(up_v/length(pert_vec))),
                                  group = c(rep("Rank of Disease Signature", length(up_signature)), 
                                            rep("Rank of Drug Signature", length(up_signature))))
        plot_df_down <- data.frame(rank_value = c(1-1:length(down_signature)/length(down_signature),1-(down_v/length(pert_vec))),
                                  group = c(rep("Rank of Disease Signature", length(down_signature)), 
                                            rep("Rank of Drug Signature", length(down_signature))))
        p_list[[i]] <- list() 
        p_list[[i]][["up"]] <- pair_plot(plot_df_up, "Distribution of Up Signature ")
        p_list[[i]][["down"]] <- pair_plot(plot_df_down, "Distribution of Down Signature ")
    }
    #output_list <- list(plot_df_up = plot_df_up, plot_df_down = plot_df_down)
    return(p_list)
}



# extract idx based on the drug 
drug_idx <- function(drug_name, compound_name) {
    return(which(conpound_name == drug_name))
}

# extract top idx
extract_top_idx <- function(score_df, top_n) {
    idx <- order(score_df$connectivity_score, decreasing = T)[1:top_n]
    #tmp_df <- score_df[order(score_df$connectivity_score, decreasing = T),][1:top_n,]
    return(idx)
    
}
                          
        
    
    
