##### clustering of the data ######, for ATAC, use the code implemented in consensus_cluster_best.R ####
BiocManager::install("cola")
library(cola)
library(readr)
raw <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_batch_corrected_TMM_normalized.tsv")
raw11 <- raw[,c(-1)]
rownames(raw11) <- raw$coordinate
matrx <- as.matrix(raw11)
matrx = cola::adjust_matrix(matrx)
# rl = cola::run_all_consensus_partition_methods(mat,
#           top_value_method = c("SD", "MAD"),
#           partition_method = c("skmeans", "pam"),
#                           cores = 1) 
# cola::cola_report(rl, output_dir = "/Users/kumarr9/Downloads/ATAC_data/",
#                   title = ("cola Report for Hierarchical Partitioning"),
#                   env = parent.frame())
# cola::select_partition_number(res)
# cola::consensus_heatmap(res, k = 3)

res = cola::consensus_partition(matrx,
                          top_value_method = "MAD",
                          top_n = c(2690, 5000),
                          partition_method = "hclust",
                          max_k = 10,
                          p_sampling = 0.9,
                          partition_repeat = 1000,
                          anno = NULL)
cola::membership_heatmap(res, k = 6)
cola::consensus_heatmap(res, k = 6)
cola::select_partition_number(res)
get_matrix(res, full = FALSE, include_all_rows = FALSE)
get_classes(res)
dimension_reduction(res, k = 5)

##############
library(ggplot2)
library(umap)
df <- read_tsv("/Users/kumarr9/Documents/raw_TMM_mad_5000.tsv")
df1 <- df[, c(-1)]
rownames(df1) <- df$coordinate
df_t <- t(df1)
umap.df = as.data.frame(df_t)
umap.df <- tibble::rownames_to_column(umap.df, "coordinate")
data1  <- read_tsv("/Users/kumarr9/Documents/df3.tsv")

allData <- cbind(umap.df, biopsy = data1$biopsy)
allData <- cbind(allData, batch = data1$batch)
allData <- cbind(allData, patient = data1$patient)
#patient <- c("p10", "p10", "p8", "p9", "p14", "p15", "p16", "p7", "p8", "p10", "p17", 
#             "p18", "p15", "p7", "p7", "p8", "p4", "p17", "p18", "p4", 
#             "p9", "p10", "p16", "p7", "p2", "p3", "p18", "p18", "p1", "p15", "p12", "p19", "p3", "p11", "p13", "p11", "p11", "p2", "p19", "p12", "p9", "p11", "p6", "p9", "p13", "p15", "p3", "p5", "p11")
# batch <- c("batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch1","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2","batch2")
# biopsy <- c("Liver","Liver","SLN","Liver","SLN","Liver","Liver","Liver","SLN","Liver","RLN","Liver","Liver","Liver","Liver","SLN","Liver","RLN","Liver","Liver","Liver","Liver","Liver","Liver","Right supraclavicularlesion","Liver","Liver","adrenal","Liver","LN","Liver","Mediastinum","Right cervical LN","Liver","Right liver lobe mass","Lung","Lung","Left adrenal mass","Right most inf liver","Right supraclavicular lesion","Liver","Right cervical LN","Liver","Right posterior infe Liver lobe","Liver","Splenic mass","Lung","Liver","Lung","Axilla","Right supra node")
# 
# umap.df <- cbind(umap.df, batch, biopsy)
# library(umap)
# penguins <- umap.df %>% 
#   mutate(ID=row_number()) 
# penguins_meta <- umap.df %>%
#   select(batch, biopsy)

penguins <- allDataaa %>% 
  mutate(ID=row_number()) 

penguins_meta <- penguins %>%
  select(ID, coordinate, patient, biopsy, batch)
set.seed(142)

umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

umap_df <- umap_fit$layout %>%
  as.data.frame()%>%
  rename(UMAP1="V1",
         UMAP2="V2") %>%
  mutate(ID=row_number())%>%
  inner_join(penguins_meta, by="ID")

umap_df %>%
  ggplot(aes(x = UMAP1, 
             y = UMAP2,
             color = coordinate)) +
  geom_point(size=3, alpha=0.5)+
  facet_wrap(~biopsy)+
  labs(x = "UMAP1",
       y = "UMAP2",
       subtitle="UMAP plot")+
  theme(legend.position="bottom")

dt <- subset(allDataaa, select = -c(biopsy,batch, patient) )
dt.data <- dt[, grep("chr", colnames(dt))]
dt.labl <- dt[, "coordinate"]

dt.umap <- umap(dt.data)
head(dt.umap$layout, 3)
plot(dt.umap, dt.labl)


spellman <- read_csv("https://github.com/Bio723-class/example-datasets/raw/master/spellman-wide.csv")
################

install.packages("healthyR")
library(healthyR)
library(healthyR.data)
library(dplyr)
library(broom)
library(ggplot2)
data_tbl <- healthyR_data %>%
  filter(ip_op_flag == "I") %>%
  filter(payer_grouping != "Medicare B") %>%
  filter(payer_grouping != "?") %>%
  select(service_line, payer_grouping) %>%
  mutate(record = 1) %>%
  as_tibble()
uit_tbl <- kmeans_user_item_tbl(data_tbl, service_line, payer_grouping, record)

uit_tbl
glimpse(allData)
column_names <- names(allData)
target_col <- "coordinate"
predictor_cols <- setdiff(column_names, target_col)

healthyR.ai::h2o.init()

output <- hai_kmeans_automl(
  .data = allData,
  .predictors = predictor_cols,
  .standardize = FALSE
)

h2o.shutdown(prompt = FALSE)
x <- rbind(cbind(rnorm(10,0,0.5), rnorm(10,0,0.5)), 
           cbind(rnorm(15,5,0.5), rnorm(15,5,0.5))) 
dtt <- subset(allData, select = -c(coordinate,biopsy,batch, patient) )
clusplot(pam(dtt, 2))
library("cluster") 
clusplot(pam(allData, 4)) 
uit_tbl <- kmeans_user_item_tbl(data_tbl, service_line, payer_grouping, record)

library(healthyR)
#library(healthyR.data)
library(dplyr)
library(broom)
library(ggplot2)
library(readr)
library(tidyverse)

data_tbl <- healthyR_data %>%
  filter(ip_op_flag == "I") %>%
  filter(payer_grouping != "Medicare B") %>%
  filter(payer_grouping != "?") %>%
  select(service_line, payer_grouping) %>%
  mutate(record = 1) %>%
  as_tibble()
# uit_tbl <- kmeans_user_item_tbl(data_tbl, service_line, payer_grouping, record)
# kmm_tbl <- kmeans_mapped_tbl(uit_tbl)
tabb <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/raw_TMM_mad_5000.tsv")
# kmm_tbl <- kmeans_mapped_tbl(tabb)
# kmeans_scree_plt(.data = kmm_tbl)
# record <- c(1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1) 
# test$record = record
# ump_lst <- umap_list(.data = test, kmm_tbl, 3)
# ump_lst <- umap_list(.data = uit_tbl, kmm_tbl, 3)
# test1 <- test[, c(-3, -4)]
# uitt_tbl <- kmeans_user_item_tbl(test1, sample, biopsy, record)
tabb1 <- tabb[, c(-1)]
rownames(tabb1) <- tabb$coordinates  
library(factoextra)
km.res <- kmeans(tabb1, 5, nstart = 25)
dd <- cbind(tabb1, cluster = km.res$cluster)
rkmm_tbl <- kmeans_mapped_tbl(dd)
ump_lst <- umap_list(.data = tabb, rkmm_tbl, 5)
umap_plt(.data = ump_lst, .point_size = 3, TRUE)
tr_tabb <- as.data.frame(t(tabb[,-1]))
colnames(tr_tabb) <- tabb$coordinates
info <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df5.tsv")
tr_tabb$biopsy <- info$biopsy
tr_tabb$sample <- info$sample
tr_tabb$patient <- info$patient
tr_tabb$batch <- info$batch
tr_tabb$cluster <- info$cluster

penguins <- tr_tabb %>% 
  drop_na() %>%
  mutate(ID=row_number())
penguins_meta <- penguins %>%
  select(ID, biopsy, sample, patient, batch, cluster)
set.seed(142)
library(umap)
umap_fit <- penguins %>%
  select(where(is.numeric)) %>%
  column_to_rownames("ID") %>%
  scale() %>% 
  umap()

###########################################################
######################## peak annotation ##################
###########################################################
library(tidyverse)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(clusterProfiler)
#BiocManager::install("ReactomePA")
library(ReactomePA)
# files <- getSampleFiles()
# print(files)
# 
# all_peaks <- list.files(path = "/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/", pattern = "*.narrowPeak")
# colnames(all_peaks) <- c('chr','start','end','name','score','strand','signal','pval','qval','peak')
# peaks.gr = makeGRangesFromDataFrame(all_peaks,keep.extra.columns=TRUE)
# for (f in all_peaks) {
#   colnames(f) <- c('chr','start','end','name','score','strand','signal','pval','qval','peak')
#   
# }



Peaks <- list.files(path = "/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/",
                    pattern = glob2rx("*.narrowPeak*"),
                    full.names = TRUE)
mPeak = GRanges()
for(hist in Peaks){
  peakRes = read.table(hist, header = FALSE, fill = TRUE)
  mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
}


peakAnno <- annotatePeak(mPeak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

peak_df <- as.data.frame(peakAnno@anno)
peak_df$coordinate = paste0(peak_df$seqnames, ":", peak_df$start, "-", peak_df$end)
peak_df = peak_df[!duplicated(peak_df$coordinate),]
rownames(peak_df) = peak_df$coordinate
#tt <- data.frame(peakAnno)
upsetplot(peakAnno, vennpie=TRUE)
st_bar <- data.frame(peakAnno@annoStat)
# ggplot(st_bar, aes(fill = Feature,
#                       y = Frequency, x = Feature))+
#   geom_bar(position = "fill", stat = "identity")
# 
# 
# BLUE <- "#076fa2"
# plt <- ggplot(st_bar) +
#   geom_col(aes(Feature, Frequency), fill = BLUE, width = 0.6) +coord_flip()
# 
# plt
#############################################################
######## bar/radar/circular plot for all the ATAC data ######
#############################################################
### link alternative to bar plot -- https://r-craft.org/r-news/bar-plots-and-modern-alternatives/
library(grid)
library(tidyverse)
library(shadowtext)
library(stringr)
library(ggpubr)
library(hrbrthemes)
st_bar <- data.frame(peakAnno@annoStat)
#### Link for code ---  https://www.data-to-viz.com/caveat/circular_barplot_accordeon.html
# dataa <- read.table("https://raw.githubusercontent.com/holtzy/data_to_viz/master/Example_dataset/7_OneCatOneNum.csv", header=TRUE, sep=",")
# dataa %>%
# filter(!is.na(Value)) %>%
#   arrange(Value) %>%
#   tail(6) %>%
#   mutate(Country=factor(Country, Country)) %>%
#   ggplot( aes(x=Country, y=Value) ) +
#   geom_bar(fill="#69b3a2", stat="identity") +
#   geom_text(hjust = 1, size = 3, aes( y = 0, label = paste(Country," "))) +
#   theme_ipsum() +
#   theme(
#     panel.grid.minor.y = element_blank(),
#     panel.grid.major.y = element_blank(),
#     legend.position="none",
#     axis.text = element_blank()
#   ) +
#   xlab("") +
#   ylab("") +
#   coord_polar(theta = "y") +
#   ylim(0,15000)



st_bar <- st_bar %>%       ### to round off the data upto one digit
  mutate_if(is.numeric,round, digits = 1)
png(filename = "/Users/kumarr9/Downloads/ATAC_data/ATAC_distributionn", res=300, width=2700, height=2700)
st_bar %>%
  arrange(Frequency) %>%
  mutate(Feature=factor(Feature, Feature)) %>%
  ggplot( aes(x=Feature, y=Frequency) ) +
  geom_bar(fill="#FA8072", stat="identity") + theme(axis.title.x=element_blank(),
                                                    axis.text.x=element_blank(),
                                                    axis.ticks.x=element_blank()) +
  geom_text(hjust = 1, size = 3, aes( y = 0, label = paste(Feature," ", Frequency,"%"))) +
  theme_ipsum() +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    legend.position="none",
    axis.text = element_blank()
  ) + theme(axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank())+
  xlab("") +
  ylab("") +
  coord_polar(theta = "y") +    #### replace "y" by "x" to have different plot 
  ylim(0,33) 
dev.off()
### when want to add title just un-comment the lines below
#+ ggtitle("ATAC feature distribution") + theme(
#    plot.title = element_text(color="black", size=10, face="plain"))
                     



##### calculate the ATAC features of each sample ##############

ttt <- annotatePeak("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/2811630_S35_L003_peaks.narrowPeak", tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

tttt <- data.frame(ttt@annoStat)

###################################################################

#####################################################
############# stacked bar plot ######################
#####################################################
## link -   https://nanx.me/ggsci/articles/ggsci.html  ## when have to change color schema
## lin -- https://rpubs.com/TX-YXL/752329  ### when have to change color schema
library(ggplot2)
library(dplyr)
library(scales)
library(readr)
library(ggsci)
library(viridis)

data3 <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_features.tsv")

##when have to make proprtional stacked bar
png(filename = "/Users/kumarr9/Downloads/ATAC_data/ATAC_features", res=300, width=2500, height=2500)
data4 <- data3 %>% group_by(sample) %>% mutate(percent = value / sum(value) * 100)
ggplot(data4, aes(x = sample, y = percent, fill = Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="plain")) + scale_fill_igv() + theme_minimal() + labs(title="ATAC features distribution")+
  # move the title text to the middle
  theme(plot.title=element_text(hjust=0.5))
dev.off()
# ATAC_new <- TMSCLC %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
# ggplot(ATAC_new, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
#   geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()
# 
# 
# Chip_new <- chip %>% group_by(sample) %>% mutate(percent_value = value / sum(value) * 100)
# ggplot(Chip_new, aes(x = sample, y = percent_value, fill = ATAC_Feature)) + coord_flip() +
#   geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv()

ggplot(data4[order(data4$Feature,decreasing=T),], aes(x = sample, y = percent, fill = Feature)) + coord_flip() +
  geom_col() + theme(axis.text = element_text(face="bold")) + scale_fill_igv() + ggtitle("Feature Distribution")

####Proportional stack bar end here ###


##############################################################
################## getting enrichment plot for RA22 samples #######
##############################################################
RA22_1_S217 <- readPeakFile("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/RA22-1_S217_peaks.narrowPeak")
RA22_1_S217_gene <- seq2gene(RA22_1_S217, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
RA22_1_S217_pathway <- enrichPathway(RA22_1_S217_gene)
dotplot(RA22_1_S217_pathway)


RA22_12_S210 <- readPeakFile("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/RA22-12_S210_peaks.narrowPeak")
RA22_12_S210_gene <- seq2gene(RA22_12_S210, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
RA22_12_S210_pathway <- enrichPathway(RA22_12_S210_gene)
dotplot(RA22_12_S210_pathway)


RA22_18_S218 <- readPeakFile("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/RA22-18_S218_peaks.narrowPeak")
RA22_18_S218_gene <- seq2gene(RA22_18_S218, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
RA22_18_S218_pathway <- enrichPathway(RA22_18_S218_gene)
dotplot(RA22_18_S218_pathway)


RA22_4_S203 <- readPeakFile("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/RA22-4_S203_peaks.narrowPeak")
RA22_4_S203_gene <- seq2gene(RA22_4_S203, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
RA22_4_S203_pathway <- enrichPathway(RA22_4_S203_gene)
dotplot(RA22_4_S203_pathway)

RA22_5_S216 <- readPeakFile("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/RA22-5_S216_peaks.narrowPeak")
RA22_5_S216_gene <- seq2gene(RA22_5_S216, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
RA22_5_S216_pathway <- enrichPathway(RA22_5_S216_gene)
dotplot(RA22_5_S216_pathway)

RA22_6_S215 <- readPeakFile("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/RA22-6_S215_peaks.narrowPeak")
RA22_6_S215_gene <- seq2gene(RA22_6_S215, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
RA22_6_S215_pathway <- enrichPathway(RA22_6_S215_gene)
dotplot(RA22_6_S215_pathway)

# this geneeee stores the gene name for RA22 primary lung sample i.e. RA22_S210
geneee <- c("101101776","100874261","100653046","100616209","100507412","100506990","100500862","100500815","100463486","100462981","100288687","100288527","100288520","100288142","100287898","100287404","100287178","100287102","100272216","100170219","100133121","100132948","100132406","100132403","100132352","100130581","729873","729515","729447","729444","728936","728841","728730","728411","728410","728393","728262","728137","727897","727758","653548","653545","653544","653188","649159","574029","554226","548645","441328","441317","441058","441056","391622","389831","388697","388581","387628","387036","375719","347735","344558","285441","284802","284800","284565","284412","284124","283788","254958","221935","200159","164045","162699","154761","153579","151987","127833","122481","120114","85443","85318","85002","84498","84221","84105","83481","83473","81704","81037","80728","80173","79034","65217","64795","64094","57597","57521","57491","57118","55755","55733","55672","54039","51078","50808","26583","26080","23784","23370","22947","11235","11039","11033","10438","10352","10267","9696","9612","9501","9373","9270","9201","8912","5799","5578","5522","4987","4588","4583","4018","3782","3607","2821","2543","2188","2069","2009","1719","1501","775","417","382","321","105","88","29")
yy = enrichPathway(geneee, pvalueCutoff=0.05)
#head(summary(yy))
#xx <- data.frame(yy@gene)
dotplot(yy)

#RA22_s210_gene <- data.frame(RA22_12_S210_pathway@gene)
#write.csv(RA22_s210_gene, file = "/Users/kumarr9/Downloads/RA22_gene.csv")
########################################################################################
##### code for cluster significance in terms of peaks #######
##############################################################
# link -- http://www.sthda.com/english/articles/32-r-graphics-essentials/132-plot-grouped-data-box-plot-bar-plot-and-more/
clustpeakdata <- read_tsv("/Users/kumarr9/Downloads/ATAC_data/ATAC_batch1_2_peaks/cluster.tsv")

clustpeakdata$cluster <- as.factor(clustpeakdata$cluster)

compare_means(peaks ~ cluster,  data = clustpeakdata)

my_comparisons <- list( c("1", "2"), c("1", "3"), c("1", "4"), c("1", "5"), c("2", "3"), c("2", "4"), c("2", "5"), c("3", "4"), c("3", "5"), c("4", "5") )
ggboxplot(clustpeakdata, x = "cluster", y = "peaks",
          color = "cluster", palette = "aaas")+ 
  stat_compare_means(comparisons = my_comparisons)+ 
  stat_compare_means(label.y = 50)   



###############################################################################################
############ ATAC Differential Analysis and peak annotation ###################################
##############################################################################################
library(rtracklayer)
library(ComplexHeatmap)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
library(org.Hs.eg.db)
library(clusterProfiler)
library(DESeq2)
library(edgeR)


# genes.gtf = genes(txdb)
# gene_map = bitr(geneID = genes.gtf$gene_id, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = org.Hs.eg.db)
# gene_map = gene_map[!duplicated(gene_map$ENTREZID),]
# rownames(gene_map) = gene_map$ENTREZID

# genes.gtf$SYMBOL=gene_map[genes.gtf$gene_id, "SYMBOL"]
# genes.gtf = genes.gtf[!is.na(genes.gtf$SYMBOL)]
# genes.gtf = sort(genes.gtf)

sample_infor = read.table(file= "/Users/kumarr9/Downloads/ATAC_data/ATAC_final_files/df4.tsv", header = T, sep = "\t")
atac_df = read.table(file= "/Users/kumarr9/Downloads/ATAC_data/raw_coverage.tsv", header = T, sep = "\t", row.names = 1, check.names = F)
colnames(atac_df) = gsub(pattern = ".sorted.dedup.bam", replacement = "", colnames(atac_df)) ####removing pattern from the colnames

##### making column names more beautiful #######
sample_names = colnames(atac_df)
for(i in 1:length(sample_names)) {
  sample_names[i] = strsplit(x = sample_names[i], split = "_")[[1]][1]
}
colnames(atac_df) = sample_names

###### cleaning and applying sample metadata #####
sample_names = sample_infor$sample
for(i in 1:length(sample_names)) {
  sample_names[i] = strsplit(x = sample_names[i], split = "_")[[1]][1]
}
sample_infor$sample = sample_names

rownames(sample_infor) = sample_infor$sample #### attaching sample names as in row also

sample_infor = sample_infor[colnames(atac_df),] #### change column names of the atac dataframe

####  getting gene length #################
gene_lengths = data.frame(GeneID = rownames(atac_df),
                          Length = width(GRanges(rownames(atac_df))))

rownames(gene_lengths) = gene_lengths$GeneID

GeneDF_EdgeR <- edgeR::DGEList(counts = atac_df, genes = gene_lengths)
GeneDF_Norm  <- edgeR::calcNormFactors(GeneDF_EdgeR, method = 'TMM')
tmm_rpkm <- as.data.frame(edgeR::rpkm(GeneDF_Norm, normalized.lib.sizes = TRUE, log = FALSE))
tmm_rpkm = round(tmm_rpkm,1)
#tmm_fpkm.z = as.data.frame(t(scale(t(tmm_fpkm))))
tmm_rpkm_log2 = log2(tmm_rpkm + 1)
tmm_rpkm_log2 = round(tmm_rpkm_log2,1)

library("DESeq2")

dds <- DESeqDataSetFromMatrix(countData = round(atac_df), ### this should be raw count, atac_df have the raw count
                              colData = sample_infor,
                              design = ~ biopsy + cluster) #### when use raw count should include batch also
dds$biopsy <- relevel(dds$biopsy, ref = "Lung")
dds <- DESeq(dds)
resultsNames(dds)   #to see the coefficients 
#res <- results(dds)  

######### when using Lymph as relevel ######
lymph_dds <- DESeqDataSetFromMatrix(countData = round(atac_df), ### this should be raw count, atac_df have the raw count
                              colData = sample_infor,
                              design = ~ biopsy + cluster) #### when use raw count should include batch also
lymph_dds$biopsy <- relevel(lymph_dds$biopsy, ref = "Lymph")
lymph_dds <- DESeq(lymph_dds)
resultsNames(lymph_dds)


######### when using Adrenal as relevel ######
adrenal_dds <- DESeqDataSetFromMatrix(countData = round(atac_df), ### this should be raw count, atac_df have the raw count
                                    colData = sample_infor,
                                    design = ~ biopsy + cluster) #### when use raw count should include batch also
adrenal_dds$biopsy <- relevel(adrenal_dds$biopsy, ref = "Adrenal")
adrenal_dds <- DESeq(adrenal_dds)
resultsNames(adrenal_dds)

###### using cluster 1 as relevl #############
cluster1_dds <- DESeqDataSetFromMatrix(countData = round(atac_df), ### this should be raw count, atac_df have the raw count
                                      colData = sample_infor,
                                      design = ~ cluster) #### when use raw count should include batch also
cluster1_dds$cluster <- relevel(cluster1_dds$cluster, ref = "cluster_1")
cluster1_dds <- DESeq(cluster1_dds)
resultsNames(cluster1_dds)

### use raw count, without normalize
### no need to normalize
### first use the dds object and then the relevel step, and then use main DEseq2 function


##############  peak annotation for dataframe used in the study ##############################################
peakAnno <- annotatePeak(GRanges(rownames(atac_df)), tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

peak_df <- as.data.frame(peakAnno@anno)
peak_df$coordinate = paste0(peak_df$seqnames, ":", peak_df$start, "-", peak_df$end)
peak_df = peak_df[!duplicated(peak_df$coordinate),]
rownames(peak_df) = peak_df$coordinate

### merging two dataframe, this is done to get NE/Non-NE score using ATAC data only, instead of RNA seq data --

peak_df_test = peak_df
atac_df_test <- atac_df
atac_df_test$Gene <- peak_df_test$SYMBOL
atac_df_test$annotation <- peak_df_test$annotation
write.csv(atac_df_test, "/Users/kumarr9/Downloads/ATAC_data/Gene_mapped_ATAC_data.csv", row.names=FALSE)
promoter <- subset(atac_df_test, annotation %in% c('Promoter (<=1kb)'))
write.csv(promoter, "/Users/kumarr9/Downloads/ATAC_data/Gene_mapped_ATAC_data_only_Promoter.csv", row.names=FALSE)
ditsal_region <- subset(atac_df_test, annotation %in% c('Distal Intergenic'))
write.csv(ditsal_region, "/Users/kumarr9/Downloads/ATAC_data/Gene_mapped_ATAC_data_distal_intergenic.csv", row.names=FALSE)

###########################################################
### this will make plot with ATAC peak coordinates  #######
library(apeglm)
res2 <- lfcShrink(dds, coef="biopsy_Lymph_vs_Lung", type="apeglm")

library("EnhancedVolcano")
EnhancedVolcano(toptable = res2,              # We use the shrunken log2 fold change as noise associated with low count genes is removed 
                x = "log2FoldChange",           # Name of the column in resLFC that contains the log2 fold changes
                y = "padj",                     # Name of the column in resLFC that contains the p-value
                lab = rownames(res2),
                title = 'Differential peaks between Liver vs Lung',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 2.0)


######### plots for using Lung as relevel #####
res2 <- lfcShrink(cluster1_dds, coef="cluster_cluster_5_vs_cluster_1", type="apeglm")
########## when want to add gene name to the peaks ###########
res2$gene = as.character(peak_df[rownames(res2), "SYMBOL"])
res2$annotation = as.character(peak_df[rownames(res2), "annotation"])
res2$annotation = sapply(1:nrow(res2), function(i) {
  strsplit(x = res2$annotation[i], split = " ")[[1]][1]
})
png(filename = "/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/cluster_5_vs_cluster_1.png", res=300, width=3000, height=3000)
#EnhancedVolcano(res2, x="log2FoldChange", y= "pvalue", lab = res2$gene)
EnhancedVolcano(res2,
                lab = res2$gene,
                x = 'log2FoldChange',
                y = 'pvalue',
                title = 'Differential peaks associated genes in cluster5 compared to cluster1',
                pCutoff = 0.01,
                FCcutoff = 2,
                pointSize = 2.0,
                labSize = 3.0,
                ylim = c(0, -log10(10e-20)))
dev.off()

# png(filename = "/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/biopsy_Adrenal_vs_Lung.png", res=300, width=3000, height=3000)
# selected <- c("ZNF407", "CCND1")
# EnhancedVolcano(res2,
#                 lab = res2$gene,
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 title = 'Differential peaks associated genes in Adrenal compared to Lung',
#                 pCutoff = 0.01,
#                 FCcutoff = 2,
#                 pointSize = 2.0,
#                 labSize = 3.0, selectLab = selected, 
#                 ylim = c(0, -log10(10e-20)))
# dev.off()

write.csv(as.data.frame(res2), 
        file="/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/cluster_5_vs_cluster_1.csv")

##### save DGE results positive ####
# library(dplyr)
# alpha <- 0.01
# diffexp <- as.data.frame(res2) %>%
#   mutate(padj = p.adjust(pvalue, method="BH")) %>%
#   filter(log2FoldChange > 2, padj < alpha)
tt <- res2[which(res2$log2FoldChange > 2 & res2$pvalue < 0.01),]
pos_sig <- as.data.frame(tt)

ttt <- res2[which(res2$log2FoldChange < - 2 & res2$pvalue < 0.01),]
neg_sig <- as.data.frame(ttt)

result <-rbind(pos_sig, neg_sig)
write.csv(result, 
          file="/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/cluster_5_vs_cluster_1_sig.csv")
##### save DGE results negative ####
## Link - https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
# Link -- https://support.bioconductor.org/p/111631/
# alpha <- 0.01
# diffexp1 <- as.data.frame(res2) %>%
#   mutate(padj = p.adjust(pvalue, method="BH")) %>%
#   filter(log2FoldChange < -2, padj < alpha)

####combine both positive and negative #####
# google search - combine two dataframe in r one below another
# Link -- https://stackoverflow.com/questions/7229608/combine-dataframe-in-the-bottom-of-another-dataframe
# result<-rbind(diffexp, diffexp1)
# write.csv(result, 
#           file="/Users/kumarr9/Downloads/ATAC_data/Liver_vs_Lung_sig.csv")

### Link -- https://hbctraining.github.io/DGE_workshop/lessons/05_DGE_DESeq2_analysis2.html
# padj.cutoff <- 0.01
# lfc.cutoff <- 1
# sigOE <- as.data.frame(res2) %>%
#   filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)
# sigOE <- liver_gene %>%
#   filter(padj < padj.cutoff & abs(log2FoldChange) > lfc.cutoff)

# mask <- liver_gene$padj < 0.01 &
#   abs(liver_gene$log2FoldChange) > log2(2)
# deGenes <- rownames(mask[mask, ])
# resLFC$gene = as.character(peak_df[rownames(resLFC), "SYMBOL"])
# resLFC$annotation = as.character(peak_df[rownames(resLFC), "annotation"])
# resLFC$annotation = sapply(1:nrow(resLFC), function(i) {
#   strsplit(x = resLFC$annotation[i], split = " ")[[1]][1]
# })
# EnhancedVolcano(resLFC,
#                 lab = resLFC$gene,
#                 x = 'log2FoldChange',
#                 y = 'pvalue',
#                 title = 'batch 2 vs batch 1 plot',
#                 pCutoff = 0.01,
#                 FCcutoff = 2,
#                 pointSize = 3.0,
#                 labSize = 2.0)
# 

############ GSEA Plots of the Dataset ######
#### Link -- https://www.biostars.org/p/445059/
### Link -- https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# BiocManager::install("clusterProfiler", version = "3.8")
# BiocManager::install("pathview")
# BiocManager::install("enrichplot")
library(clusterProfiler)
library(enrichplot)
# we use ggplot2 to add x axis labels (ex: ridgeplot)
library(ggplot2)
organism = "org.Hs.eg.db"   # SET THE DESIRED ORGANISM HERE
library(organism, character.only = TRUE)
library(clusterProfiler)
library(enrichplot)

########### GO terms GSEA enrichment analysis from DESeq2 #################
### Link - https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
# reading in data from deseq2
df = read.csv("/Users/kumarr9/Downloads/ATAC_data/Supraclavicular_sig_vs_Lung_copy.csv", header=TRUE)

# we want the log2 fold change 
original_gene_list <- df$log2FoldChange

# name the vector
names(original_gene_list) <- df$X

# omit any NA values 
gene_list<-na.omit(original_gene_list)
write.csv(df, file ="/Users/kumarr9/Downloads/ATAC_data/dfff.csv")


# sort the list in decreasing order (required for clusterProfiler)
gene_list = sort(gene_list, decreasing = TRUE)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             nPerm = 10000, 
             minGSSize = 3, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = organism, 
             pAdjustMethod = "none")
dev.off()
require(DOSE)
dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
gseaplot(gse, by = "all", title = gse$Description[1], geneSetID = 1)
ridgeplot(gse) + labs(x = "enrichment distribution")
#####################################################################
############ KEGG terms GSEA enrichment analysis from DESeq2 ########
#####################################################################
cl4 <- read.csv("/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/cluster_4_vs_cluster_1_sig.csv", sep=",")
dedup_cl4 = cl4[!duplicated(cl4[c("gene")]),]
# we want the log2 fold change 
gene_list <- dedup_cl4$log2FoldChange

# name the vector
names(gene_list) <- dedup_cl4$gene
# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted

ids<-bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=organism)

# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = dedup_cl4[dedup_cl4$gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$log2FoldChange

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

# create/obtain KEGG organism 
# Full list is here -- https://www.genome.jp/kegg/catalog/org_list.html
kegg_organism = "hsa"

kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 1000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

dotplot(kk2, showCategory = 10, title = "KEGG Enriched Pathways for significant genes" , split=".sign") + facet_grid(.~.sign)
dev.off()
# categorySize can be either 'pvalue' or 'geneNum'
# install.packages(ggnewscale)
# library("ggnewscale")
cnetplot(kk2, categorySize="pvalue", foldChange=gene_list)
dev.off()

ridgeplot(kk2) + labs(x = "KEGG pathway enrichment")
dev.off()

gseaplot(kk2, by = "all", title = kk2$Description[1], geneSetID = 1)

#########################################################
############ GSEA with hallmark gene sets ###############
##########################################################
## Link --- https://sbc.shef.ac.uk/workshops/2019-01-14-rna-seq-r/rna-seq-gene-set-testing.nb.html
## Link - https://www.biostars.org/p/467197/
## Link of above in github -- https://github.com/hamidghaedi?tab=overview&from=2022-11-01&to=2022-11-30
## Link --- https://github.com/hamidghaedi/Enrichment-Analysis
## Link -- https://hbctraining.github.io/DGE_workshop_salmon_online/lessons/11_FA_functional_class_scoring.html
## Link -- https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html ### when adding more than one/say top3-10 pathway to plot
#ress <- read.csv("/Users/kumarr9/Downloads/Liver_vs_lung_copy.csv", header = T, sep = ",")
require(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
library(fgsea)
library(dplyr)
gseaInput <- filter(ress, !is.na(SYMBOL)) %>% 
  arrange(log2FoldChange)
ranks <- pull(gseaInput,log2FoldChange)
names(ranks) <- gseaInput$SYMBOL
barplot(ranks)
pathways.hallmark <- gmtPathways("/Users/kumarr9/Downloads/h.all.v2022.1.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways.hallmark, ranks, minSize=15, maxSize = 500, nperm=1000)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy #### save it as dataframe

library(ggplot2)
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways NES from GSEA")
plotEnrichment(pathways.hallmark[["HALLMARK_MYC_TARGETS_V1"]],
               ranks)
rr <- plotEnrichment(pathways.hallmark[["HALLMARK_MYC_TARGETS_V1"]], ### we want the name to be aded in the plot
               ranks)
rr + ggtitle("HALLMARK_MYC_TARGETS_V1") + theme(
  plot.title=element_text( hjust=1, vjust=0.5,margin=margin(t=40,b=-30))
)




########################## on the other data ############################
llc <- read.csv("/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/cluster_5_vs_cluster_1.csv", header = T, sep = ",")
gseaInput <- filter(llc, !is.na(gene)) %>% 
  arrange(log2FoldChange)
ranks <- pull(gseaInput,log2FoldChange)
names(ranks) <- gseaInput$gene
#barplot(ranks)
pathways.hallmark <- gmtPathways("/Users/kumarr9/Downloads/h.all.v2022.1.Hs.symbols.gmt")
fgseaRes <- fgsea(pathways.hallmark, ranks, minSize=10, maxSize = 500, nperm=100)
fgseaResTidy <- fgseaRes %>%
  as_tibble() %>%
  arrange(desc(NES))
# Show in a nice table:
fgseaResTidy      #### save it as dataframe
dff <- apply(fgseaResTidy,2,as.character)
write.csv(dff, file = "/Users/kumarr9/Downloads/ATAC_data/ATAC_DGE_parth_suggestion/cluster_5_vs_cluster_1_GSEA.csv", row.names = FALSE)
###### Bar plot for enrichment analysis #####
fgseaResTidy$Pvalue <- ifelse(fgseaResTidy$pval <= 0.05, "significant", "non-significant")
cols <- c("non-significant" = "grey", "significant" = "red")
ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES, fill = Pvalue)) +
  geom_col() +
  scale_fill_manual(values = cols) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark pathways Enrichment Score from GSEA")
# library(ggplot2)
# ggplot(fgseaResTidy, aes(reorder(pathway, NES), NES)) +
#   geom_col(aes(fill=padj<0.01)) +
#   coord_flip() +
#   labs(x="Pathway", y="Normalized Enrichment Score",
#        title="Hallmark pathways NES from GSEA")

#plotEnrichment(pathways.hallmark[["HALLMARK_MITOTIC_SPINDLE"]],
#               ranks)
rr <- plotEnrichment(pathways.hallmark[["HALLMARK_ESTROGEN_RESPONSE_EARLY"]], ### we want the name to be aded in the plot
                     ranks)
rr + ggtitle("HALLMARK_ESTROGEN_RESPONSE_EARLY") + theme(
  plot.title=element_text( hjust=1, vjust=0.5,margin=margin(t=40,b=-30))
)

### use this link when have to add more than one pathway name in the plot --
## link -- https://yulab-smu.top/biomedical-knowledge-mining-book/enrichplot.html

