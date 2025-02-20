## Read arguments
args = commandArgs(trailingOnly=TRUE)
idx = as.numeric(args[1])
data_file = args[2]

# library 
library(dplyr)
library(data.table)

# "./data/gsea_res_2019_09_24.RData"
load(data_file)


# GSEA connectivity score
gsea_drug <- function(num_gene, rank_matrix, up_signature, down_signature, complete_cor_vec, compound_name) {
    row_mtx <- nrow(rank_matrix)
    drug_top_mtx <- rank_matrix[1:num_gene, ]
    drug_button_mtx <- rank_matrix[(row_mtx-num_gene+1):row_mtx,]
    es_vec <- c()
    output_list <- list()
    for (i in 1:ncol(rank_matrix)) {
        if (i %% 200 == 0) { print(paste0("INFO: ", i))}
        top_vec <- drug_top_mtx[, i]
        button_vec <- drug_button_mtx[, i]
        top_cor_vec <- complete_cor_vec[complete_cor_vec_df$entrezgene_id %in% top_vec]
        button_cor_vec <- complete_cor_vec[complete_cor_vec_df$entrezgene_id %in% button_vec]
        up_pathway <- list(deg = as.character(up_signature))
        down_pathway <- list(deg = as.character(down_signature))
        up_res <- fgsea::fgsea(up_pathway, top_cor_vec, nperm=1000, maxSize=500)
        down_res <- fgsea::fgsea(down_pathway, button_cor_vec, nperm=1000, maxSize=500)
        if (nrow(up_res) == 0) { up_score <- 0 }
        else {up_score <- up_res$ES} 
        if (nrow(down_res) == 0) { down_score <- 0 }
        else {down_score <- down_res$ES}
        
        es <-  up_score - down_score 
        es_vec <- c(es_vec, es)
        output_list[[i]] <- list(up_gsea_results = up_res, down_gsea_results = down_res) 
    }
    gsea_df<- data.frame(drug_name = compound_name, gsea_enrichment_score = es_vec)
    gsea_df_ordered <- gsea_df %>% 
        dplyr::arrange(desc(gsea_enrichment_score))
    output_list[["gsea_df_ordered"]] <- gsea_df_ordered
    return(output_list)
    
}


# check overlap
check_overlap <- function(gsea_df, nash_drugs) {
    lincs_drugs <- tolower(gsea_df$drug_name)
    
    midx <- na.omit(match(nash_drugs, lincs_drugs))
    top5 <- lincs_drugs[(5774-289+1):5774]
    overlapped_drugs <- lincs_drugs[midx]
    sum_overlap <- ifelse(lincs_drugs %in% overlapped_drugs, 1, 0)
    #sum_overlap_top5 <- ifelse(lincs_drugs[(5774-289+1):5774] %in% overlapped_drugs, 1, 0)
    sum_overlap_top5 <- ifelse(lincs_drugs[1:289] %in% overlapped_drugs, 1, 0)
    #sum_overlap_top10 <- ifelse(lincs_drugs[(5774-578+1):5774] %in% overlapped_drugs, 1, 0)
    sum_overlap_top10 <- ifelse(lincs_drugs[1:578] %in% overlapped_drugs, 1, 0)
    sum_overlap_top_all <- ifelse(lincs_drugs %in% overlapped_drugs, 1, 0)
    return(c(sum(sum_overlap_top5), sum(sum_overlap_top10)))
}


#save(rank_matrix, row_mtx, up_signature, down_signature, complete_cor_vec_df, complete_cor_vec, gsea_df, compound_name, nash_drugs, file = "./data/gsea_res_2019_09_24.RData")




gsea_list <- gsea_drug(idx, rank_matrix, up_signature, down_signature, complete_cor_vec, compound_name)
save(gsea_list, file = paste0("../data/gsea_connectivity_score_list_", idx, ".RData"))
