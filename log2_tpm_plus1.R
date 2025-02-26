library(tidyverse)
library(readr)
library(dplyr)
all_data <-  read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/NCI_Thomas_439_protein_coding_only_TPM_normalized.tsv")
df_subset <- all_data%>%
            dplyr::select(Gene_Id, Sample_23_SB_20_788, Sample_49_NC_20_182_A1, Sample_48_NC_20_182_B1,
                          sb_19_6002, Sample_47_WS_20_4104, Sample_48_ASP19_05374, Sample_10_SS_18_5516, 
                          Sample_2_WS_18_051587, Sample_36_SB_18_6369, Sample_24_AU_18_162_AA, Sample_30_512707,
                          Sample_31_512714, Sample_32_512715, Sample_12_SS_18_4268, Sample_13_SS_18_4254, SB19_368_RNA,
                          Sample_25_AU_19_79_S, Sample_26_AU_19_79_B, Sample_43_RA_18_14, SS_19_2158,
                          SB_19_2276, SB_19_3892, `27_sb_19_5484`, `25_sb_19_6623`, CL0261_T1R_T,
                          Sample_45_RA_023_11, SS19_503_RNA, Sample_26_SB_20_5744, Sample_41_RS_20_1704, Sample_25_SB_20_5907,
                          Sample_41_CN20_33_A1_B1, Sample_39_SF_20_216_A1, Sample_40_SF_20_216_B1, Sample_23_SB_20_788, Sample_9_SS_18_763, 
                          Sample_1_WS20_3053)

colnames(df_subset) <- gsub(pattern = "Sample_", replacement = "", colnames(df_subset))
write.csv(df_subset, file="/Users/kumarr9/Downloads/israel_dataset.csv", sep=",", row.names = FALSE)


#### log2 (TPM+1) normalization ####
tpm_df_subset <- read_tsv("/Users/kumarr9/Downloads/israel_dataset.tsv")
tpm_df_subset1 <- tpm_df_subset[, -1]
# TPM+1 transformation
tpm_plus1 <- log2(tpm_df_subset1 + 1)

# Log2 normalization
normalized_data <- scale(tpm_plus1, center = TRUE, scale = TRUE)
normalized_data_df <- data.frame(normalized_data)
rownames(normalized_data_df) <- tpm_df_subset$Gene_Id
# 'normalized_data' will contain the log2 normalized values
write.csv(df_subset, file="/Users/kumarr9/Downloads/israel_dataset_log2_tpm_plus1.csv", sep=",", row.names = FALSE)


### plotting the results of this as a heatmap ####
library(dplyr)
library(ComplexHeatmap)
tt <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/israel_dataset/ne_score_israel_dataset.tsv")
tt1 <- tt[, -1]
tt1 <- data.matrix(tt)
Heatmap(tt1)

library(circlize)
col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
col_fun(seq(-3, 3))

Heatmap(tt1, name = "mat", col = col_fun)


Heatmap(tt1, name = "mat", col = rev(rainbow(10)), 
        column_title = "set a color vector for a continuous matrix")
