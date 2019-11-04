# library
library(dplyr)

# KS help function
.ks <- function( V, n ) {
  t <- length( V )
  
  if( t == 0 )  {
    return( 0 )
  } else {
    
    if ( is.unsorted( V ) )
      V <- sort( V )
    d <- (1:t) / t - V / n
    a <- max( d )
    b <- -min( d ) + 1 / t
    ifelse( a > b, a, -b )
  }
}

.ks_mean <- function( V, n ) {
  t <- length( V )
  
  if( t == 0 )  {
    return( 0 )
  } else {
    
    #if ( is.unsorted( V ) )
    #  V <- sort( V )
    d <- (1:t) / t - V / n
    a <- mean( d )
    b <- -mean( d ) + 1 / t
    ifelse( a > b, a, -b )
  }
}

.s <- function( V_up, V_down, n ) {
  ks_up <- .ks( V_up, n )
  ks_down <- .ks( V_down, n )
  ifelse( sign( ks_up ) == sign( ks_down ), 0, ks_up - ks_down )
}

.s_mean <- function( V_up, V_down, n ) {
  ks_up <- .ks_mean( V_up, n )
  ks_down <- .ks_mean( V_down, n )
  ifelse( sign( ks_up ) == sign( ks_down ), 0, ks_up - ks_down )
}

.S <- function( scores ) {
  p <- max( scores )
  q <- min( scores )
  ifelse(
         scores == 0,
         0,
         ifelse( scores > 0, scores / p, -scores / q )
         )
}



# main function ks score
ks_score <- function(up_signature, down_signature, rank_matrix, compound_name) {
    n = 1
    connectivity_score_vec <- c()
    gene_rank_bing <- rank_matrix
    pert_vec <- gene_rank_bing[, 1]
    up_signature <- up_signature[up_signature %in% pert_vec]   
    down_signature <- down_signature[down_signature %in% pert_vec]
    for (i in 1:ncol(gene_rank_bing)) {
        if (i %% 500 == 0)  { print(paste0("INFO: ", i, " Instances."))}
        
        # get teh rank
        up_v <- match(up_signature, pert_vec)
        down_v <- match(down_signature, pert_vec)
        tmp_s <- .s(up_v, down_v, nrow(gene_rank_bing))
        connectivity_score_vec <- c(connectivity_score_vec, tmp_s)
        
    }
    output_df <- data.frame(compound_name = compound_name, connectivity_score = connectivity_score_vec)
    return(output_df)
    
}

# main function ks score
ks_score_mean <- function(up_signature, down_signature, rank_matrix, compound_name) {
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
        tmp_s <- .s_mean(up_v, down_v, nrow(gene_rank_bing))
        connectivity_score_vec <- c(connectivity_score_vec, tmp_s)
        
    }
    output_df <- data.frame(compound_name = compound_name, connectivity_score = connectivity_score_vec)
    return(output_df)
    
}

# main function ks score
# signature_df should include logfc
ks_score_mean_weighted <- function(up_signature, down_signature, rank_matrix, compound_name) {
    n = 1
    connectivity_score_vec <- c()
    gene_rank_bing <- rank_matrix
    pert_vec <- gene_rank_bing[, 1]
    up_signature <- up_signature[up_signature %in% pert_vec]   
    down_signature <- down_signature[down_signature %in% pert_vec]
    for (i in 1:ncol(gene_rank_bing)) {
        if (i %% 500 == 0)  { print(paste0("INFO: ", i, " Instances."))}
        # get teh rank
        up_v <- match(up_signature, pert_vec)
        down_v <- match(down_signature, pert_vec)
        tmp_s <- .s_mean(up_v, down_v, nrow(gene_rank_bing))
        connectivity_score_vec <- c(connectivity_score_vec, tmp_s)
        
    }
    output_df <- data.frame(compound_name = compound_name, connectivity_score = connectivity_score_vec)
    return(output_df)
    
}

# data loading 
#ks_score_vec <- ks_score(up_signature_df$entrezgene_id, down_signature_df$entrezgene_id, gene_rank_bing, compound_name)
#ks_score_df <- ks_score_vec %>% 
#    arrange(desc(connectivity_score))





