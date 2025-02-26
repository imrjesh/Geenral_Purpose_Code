# link - https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
library(tidyverse)
data <-  read.csv("/Users/kumarr9/Downloads/cluster1_target_genes.tsv", sep = "\t", check.names=F, row.names = 1)
#norm_exp_df.z = as.data.frame(t(scale(t(data))))

heatmap_anno <- read_tsv("/Users/kumarr9/Downloads/anno.tsv")

#ann <- data.frame(heatmap_anno$biopsy, heatmap_anno$batch, heatmap_anno$cluster, heatmap_anno$origin)
#colnames(ann) <- c('biopsy', 'batch', 'cluster', 'origin')
#colours <- list('biopsy' = c('Liver' = '#F5F5DC', 'Lymph' = '#00FFFF', 'Adrenal' = '#000000', 'Other' = '#0000FF', 'Lung' = '#FF7F50'),
#                'cluster' = c('cluster_1' = 'limegreen', 'cluster_2' = 'gold', 'cluster_3' = '#0000FF', 'cluster_4' = '#FF7F50', 'cluster_5' = '#006400'),
#                'origin' = c('PDX' = '#8B0000', 'Patient' = '#008000'),
#                'batch' = c('batch1' = '#F0E68C', 'batch2' = '#FF00FF'))
#colAnn <- HeatmapAnnotation(df = ann,
#                            which = 'col',
#                            col = colours,
#                            annotation_width = unit(c(1, 4), 'cm'),
#                            gap = unit(1, 'mm'))
library(ComplexHeatmap)
norm_exp_df.z <- as.matrix(data)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/target_gene/all_NAPY_MYCLN.png", res=300, width=7000, height=7000)
Heatmap(norm_exp_df.z, name = "z score normalized", column_split = heatmap_anno$cluster)
dev.off()



#### 
df <- read_tsv("/Users/kumarr9/Downloads/cluster3_all.tsv")
target_genes <- c("ABI3BP","ACER2","ACTA2","ACTG2","ACVR2A","ADAMTS1","ADAMTS2","ADARB1","AEBP1","AGPAT4","AHI1","AKAP2","AKT3","ALDH1L2","AMOTL1","ANGPTL2","ANKRD1","ANTXR1","APOE","AQP9","ARC","ARHGAP20","ARHGAP24","ARHGAP25","ARHGEF25","ARNTL","ARSI","ASPHD2","ATP2A2","AXL","BACH1","BAG3","BCAR1","BCL2L11","BCOR","BMP1","BMP7","BNC1","BTBD11","BVES","C10orf54","C12orf24","C13orf33","C14orf149","C14orf37","C17orf81","C19orf12","C1QTNF4","C21orf7","C2orf40","C3orf54","C3orf58","C3orf64","C4orf49","C5orf25","C6orf145","C6orf228","C8orf84","C9orf3","CACNB4","CADM1","CALD1","CALML3","CALU","CARD10","CAV1","CBLB","CCDC3","CCDC85B","CCND2","CD36","CD70","CDC42EP2","CDH13","CDH3","CDH4","CDKN1A","CHST3","CHST7","CHSY3","CLIP3","CLMP","CNN1","CNP","CNRIP1","COL12A1","COL14A1","COL16A1","COL17A1","COL18A1","COL23A1","COL4A1","COL4A2","COL5A1","COL5A2","COL7A1","COL9A2","CPNE8","CPXM1","CRISPLD1","CRLF1","CRYAB","CSDC2","CSPG4","CSRP2","CTGF","CTNNAL1","CXCL14","CYGB","D4S234E","DCBLD2","DCHS1","DCUN1D3","DKK3","DLK2","DLL1","DMWD","DOCK10","DPYSL3","DST","DUSP6","DUSP7","DZIP1L","EBF3","EDARADD","EDNRB","EEPD1","EFCAB1","EFNB1","EGFR","EGR2","EGR3","EID3","ELK3","ELOVL4","ENC1","ENPP2","EPAS1","EPDR1","EPHB1","ERF","ETS1","ETV5","EVC","EXT1","FABP5","FAM101B","FAM132A","FAM176A","FAM184A","FAM70B","FAS","FBLN1","FBLN7","FBXO30","FERMT2","FEZ1","FGFRL1","FGL2","FHL1","FHOD3","FJX1","FLNC","FLRT2","FMOD","FOXP1","FST","FXYD1","FZD8","GEM","GJA1","GJC1","GNAI1","GNB4","GNG11","GOLIM4","GPC3","GPR124","GPR176","GPR3","GPR87","GPSM1","GRASP","GSN","GYLTL1B","GYPC","HAS2","HDAC4","HEG1","HGFAC","HRAS","HS3ST3A1","HSPB2","HSPG2","HTRA1","ICAM1","ID4","IGFBP2","IGFBP3","IGFBP4","IGFBP6","IL17B","IL17RD","IL1B","IL24","IL6","IL6ST","IRX4","ISM1","ITGA1","ITGA6","ITGA9","ITGB1","ITGB4","ITM2A","JAG1","JAM2","JAM3","KANK4","KCNIP3","KCNMA1","KCNMB1","KDELC1","KIAA0889","KLHDC5","KLHL21","KLHL29","KRT14","KRT16","KRT5","KRT75","LAG3","LAMA1","LAMA3","LAMB1","LAMB3","LAMC1","LBH","LCA5","LCAT","LEP","LEPRE1","LEPREL1","LGALS1","LGALS7","LGR6","LHFP","LIFR","LIMA1","LIMS2","LMOD1","LPHN1","LRCH2","LRP1","LRP4","LRRC8C","LRRN1","LTBP4","LUZP1","MALT1","MAMDC2","MAOB","MATN2","MBNL1","MCAM","MEF2C","MEG3","MEST","MFNG","MIA","MICAL2","MME","MMP2","MPDZ","MRGPRF","MRVI1","MSRB3","MSX1","MTSS1","MXRA7","MYC","MYH11","MYL9","MYLK","MYOCD","NBL1","NDN","NETO2","NGF","NGFR","NLGN2","NNAT","NNMT","NPTX2","NRCAM","NRG1","NRP1","NRP2","NT5E","NTF3","NTRK2","NUDT10","NUDT11","NXN","ODZ3","OSBPL6","OSR1","OXTR","PAMR1","PARD6G","PCBP4","PCDH18","PCDH19","PCDH7","PCDHGC3","PCOLCE","PDGFA","PDLIM4","PDLIM7","PDPN","PEG3","PELO","PGF","PHLDA3","PHLDB1","PKD1","PKD2","PKNOX2","PKP1","PLA2G7","PLCH2","PLEKHA4","PLS3","PLXNA2","PODN","POPDC2","POSTN","POU3F1","PPAP2A","PPAP2B","PPP1R14A","PPP1R16B","PPP1R18","PPP1R3C","PPP2R2B","PRDM1","PRICKLE1","PRICKLE2","PRNP","PROS1","PRRX1","PRX","PSD2","PTGS2","PTPLA","PTPRE","PTPRT","PVRL3","PXN","QKI","QRICH2","RAB34","RAPGEF1","RARB","RARRES2","RASIP1","RASL12","RBPMS","RCN3","RCSD1","RECK","RELN","RFX2","RGNEF","RHOJ","RND3","RNF165","RUSC2","SCARF2","SCHIP1","SCML2","SCN4B","SDK2","SDPR","SEC24D","SEMA3C","SEMA5A","SERPINF1","SERPING1","SERPINH1","SGCB","SGIP1","SH2D5","SH3TC1","SHE","SIAH2","SKI","SLC12A4","SLC1A3","SLC1A5","SLC25A4","SLC27A3","SLC27A6","SLC2A3","SLC38A5","SLC4A3","SLC6A8","SLCO3A1","SLIT2","SLIT3","SMTN","SNAI2","SNCA","SNTB2","SOBP","SORBS1","SORCS1","SOX11","SPARC","SPHK1","SPRED1","SRGN","SRPX","SSBP2","SSH1","STAC","STAC2","STARD8","STXBP4","SULF1","SVEP1","SYDE1","SYNM","TACC1","TAGLN","TBX2","TCF4","TCF7L1","TCOF1","TES","TGFB1I1","TGFBR3","THBS1","THSD1","THY1","TIE1","TIMP3","TINAGL1","TM7SF3","TMEM121","TMEM178","TMEM201","TMEM204","TMEM47","TMEM64","TNS1","TNS4","TOX","TP63","TPM2","TPST1","TRIM29","TRIM9","TRO","TRPC1","TSHZ2","TSHZ3","TSKU","TSPY26P","TSPYL2","TTL","TTYH2","TUBB6","TWIST2","UCN2","UNC45A","UPP1","VCAN","VGLL3","VIM","VIT","VSNL1","WIF1","WIPF1","WTIP","YAF2","ZC3H12B","ZNF219","ZNF423")

# Use logical indexing to extract rows where any of the target genes are present
selected_rows <- df[df$Gene %in% target_genes, ]
