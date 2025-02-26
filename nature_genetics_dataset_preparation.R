library(readr)
library(dplyr)
library(tidyverse)
library(tibble)
all_tpm <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_495_protein_coding_only_TPM_normalized.tsv")
tpm_validation <- all_tpm %>%
  dplyr::select(Gene_Id, Sample_1_WS20_3053, Sample_10_SS_18_5516, Sample_11_SS_18_5261,Sample_12_SS_18_4268,  Sample_13_SS_18_4254,
                SB19_368_RNA, Sample_25_AU_19_79_S, Sample_26_AU_19_79_B, Sample_43_RA_18_14,  Sample_14_SS_18_2698, Sample_20_AU_18_47_V, 
                Sample_2_WS_18_051587, Sample_36_SB_18_6369, Sample_24_AU_18_162_AA, Sample_30_512707, Sample_31_512714, Sample_32_512715,
                Sample_26_SB_20_5744, X27_sb_19_5484, Sample_28_SB_20_399,Sample_34_SB_18_7364, Sample_39_S18_11654, Sample_33_SB_19_3400, Sample_39_SF_20_216_A1,
                Sample_40_SF_20_216_B1, Sample_41_CN20_33_A1_B1, Sample_45_RA_023_11, Sample_47_WS_20_4104,Sample_48_ASP19_05374, Sample_49_NC_20_182_A1,  Sample_48_NC_20_182_B1, CL0124_T1R_T,
                CL0170_T1R_T, CL0191_T1R_T, CL0191_T2R_T, CL0191_T3R_T, CL0191_T4R_T, Sample_33_512718, Sample_34_512719, Sample_35_512720, Sample_36_512721,CL0193_T1R_T,
                CL0261_T1R_T, SB_19_3892, sb_19_4537, sb_19_4900,sb_19_6002)
write.table(tpm_validation, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/final_validation_tpm.tsv", sep="\t")



### for raw count data table
all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_Thomas_raw_count_h.txt")
#all_rna_seq <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/two_sample.tsv")

#colnames(all_rna_seq) = gsub(pattern = "Sample_", replacement = "", colnames(all_rna_seq)) ### removing this un-necessary sample_ stuff
### Reading our protein coding gene dataset (19586 genes)
only_protein <- read.csv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/gene_v_27.csv", sep = "\t", header = T)
### collecting only protein coding genes from all_rna_seq data file
all_rna_seq.protein.coding <- merge(all_rna_seq, only_protein, by=c("Gene_Id"))
write.table(all_rna_seq.protein.coding, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_thomas_raw_count_433_protein_coding.tsv", sep = "\t", row.names = FALSE)

all_raw <- read_tsv("/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/NCI_thomas_raw_count_433_protein_coding.tsv")
raw_validation <- all_raw %>%
  dplyr::select(Gene_Id, Sample_1_WS20_3053, Sample_10_SS_18_5516, Sample_11_SS_18_5261,Sample_12_SS_18_4268,  Sample_13_SS_18_4254,
                SB19_368_RNA, Sample_25_AU_19_79_S, Sample_43_RA_18_14, Sample_26_AU_19_79_B, Sample_14_SS_18_2698, Sample_20_AU_18_47_V, 
                Sample_2_WS_18_051587, Sample_36_SB_18_6369, Sample_24_AU_18_162_AA, Sample_30_512707, Sample_31_512714, Sample_32_512715,
                Sample_26_SB_20_5744, `27_sb_19_5484`, Sample_28_SB_20_399,Sample_34_SB_18_7364, Sample_39_S18_11654, Sample_33_SB_19_3400, Sample_39_SF_20_216_A1,
                Sample_40_SF_20_216_B1, Sample_41_CN20_33_A1_B1, Sample_45_RA_023_11, Sample_47_WS_20_4104,Sample_48_ASP19_05374, Sample_49_NC_20_182_A1,  Sample_48_NC_20_182_B1, CL0124_T1R_T,
                CL0170_T1R_T, CL0191_T1R_T, CL0191_T2R_T, CL0191_T3R_T, CL0191_T4R_T, Sample_33_512718, Sample_34_512719, Sample_35_512720, Sample_36_512721,CL0193_T1R_T,
                CL0261_T1R_T, SB_19_3892, sb_19_4537, sb_19_4900,sb_19_6002)
write.table(raw_validation, file="/Users/kumarr9/Desktop/rajesh_projects/second_project/home_visit_data/NE_nonNE_calculation_all_RNA/raw_validation.tsv", sep="\t")