# check overlap
drug_str <- "ACARBOSE AMLODIPINE ATORVASTATIN BETAINE CAFFEINE CENICRIVIROC CHOLECALCIFEROL DABIGATRAN DAPAGLIFLOZIN EMRICASAN ERGOCALCIFEROL EXENATIDE EZETIMIBE FENOFIBRATE GLIMEPIRIDE INSULIN GLARGINE INSULIN HUMAN PRADIGASTAT K-877 PX-102 LIRAGLUTIDE QUINIDINE LOSARTAN ROFLUMILAST METFORMIN ROSIGLITAZONE METRELEPTIN ROSUVASTATIN MGL-3196 SELADELPAR MIDAZOLAM SELONSERTIB SEMAGLUTIDE MT-3995 MUROMONAB-CD3 SIMTUZUMAB NAMODENOSON SITAGUPTIN NIVOCASAN SPIRONOLACTONE OBETICHOLIC ACID TELMISARTAN PENTOXIFYLUNE TOFOGLIFLOZIN PERINDOPRIL URSODIOL PIOGLITAZONE"

nash_drugs <- tolower(unlist(stringr::str_split(drug_str, " ")))

library(dplyr)
check_overlap <- function(gsea_df, nash_drugs) {
    gsea_df <- gsea_df %>% arrange(desc(connectivity_score))
    lincs_drugs <- tolower(gsea_df$compound_name)

    midx <- na.omit(match(nash_drugs, lincs_drugs))
    overlapped_drugs <- lincs_drugs[midx]
    sum_overlap <- ifelse(lincs_drugs %in% overlapped_drugs, 1, 0)
    #sum_overlap_top5 <- ifelse(lincs_drugs[(5774-289+1):5774] %in% overlapped_drugs, 1, 0)
    sum_overlap_top5 <- ifelse(lincs_drugs[1:289] %in% overlapped_drugs, 1, 0)
    #sum_overlap_top10 <- ifelse(lincs_drugs[(5774-578+1):5774] %in% overlapped_drugs, 1, 0)
    sum_overlap_top10 <- ifelse(lincs_drugs[1:578] %in% overlapped_drugs, 1, 0)
    sum_overlap_top_all <- ifelse(lincs_drugs %in% overlapped_drugs, 1, 0)
    return(c(sum(sum_overlap_top5), sum(sum_overlap_top10)))
}


