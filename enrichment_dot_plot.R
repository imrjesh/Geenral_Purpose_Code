#usefule code while making this ###
### https://www.biostars.org/p/456294/
#### https://www.google.com/search?q=reorder+dot+plot+in+r&rlz=1C5GCEM_enUS964US964&ei=-hGdY9rMCuOk5NoP_cuvkAE&oq=reorder+dotplot+in+r&gs_lcp=Cgxnd3Mtd2l6LXNlcnAQAxgAMgoIIRDDBBAKEKABOgoIABBHENYEELADOgcIABCABBANOgUIABCGA0oECEEYAEoECEYYAFCuA1iSDmD6IWgBcAF4AIABZYgBkAKSAQMyLjGYAQCgAQHIAQjAAQE&sclient=gws-wiz-serp
### https://stackoverflow.com/questions/69883865/change-order-y-axis-of-dotplot-in-ggplot2

library(DOSE)
library(ggplot2)
library(enrichplot)
library(readr)
library(dplyr)

df <- read_tsv("/Users/kumarr9/Downloads/dotplot2.tsv")



ggplot(data = df, aes(x = GeneRatio, y = reorder(Pathway_name, Count), color = `p.adjust`, size = Count)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("GeneRatio") + 
  ggtitle("")

##### for fusion data #######

df1 <- read_tsv("/Users/kumarr9/Downloads/dotplot3.tsv")
DF_filtered <- filter(df1, p < 0.02)

ggplot(data = DF_filtered, aes(x = NES, y = reorder(Pathway_name, NES), color = `p`, size = NES)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("GeneRatio") + 
  ggtitle("Positively enriched GO fusion+")

DF_tail <- filter(df1, p > 0.992188)

ggplot(data = DF_tail, aes(x = NES, y = reorder(Pathway_name, NES), color = `p`, size = ES)) + 
  geom_point() +
  scale_color_gradient(low = "red", high = "blue") +
  theme_bw() + 
  ylab("") + 
  xlab("GeneRatio") + 
  ggtitle("Negatively enriched terms fusion positive patients")

grid.arrange(pos, neg)
##### code for parth to find common #####
library(dplyr)
setwd("/Users/kumarr9/Downloads")
df1 <- read.csv("TME_proteome.csv", header = T)

df2 <- read.csv("HALLMARK_KRAS_SIGNALING_DN.csv", header = T)
#new_data_frame <- merge(df1, df2) #### when same column name
common_tumor_proteome <- merge(df1, df2, by.x=c("TME_Proteome"), by.y=c("HALLMARK_KRAS_SIGNALING_DN"))
            



##### for chris
install.packages("ggdendro")
library(ggdendro)
library(cowplot)

gene_cluster <- read_tsv('https://github.com/davemcg/davemcg.github.io/raw/master/content/post/scRNA_dotplot_data.tsv.gz')
markers <- gene_cluster$Gene %>% unique()

gene_cluster %>% filter(Gene %in% markers) %>% 
  mutate(`% Expressing` = (cell_exp_ct/cell_ct) * 100) %>% 
  ggplot(aes(x=cluster, y = Gene, color = count, size = `% Expressing`)) + 
  geom_point() 

df3 <- read_tsv("/Users/kumarr9/Downloads/chris_gsea.tsv")
ggplot(df3, aes(x = Difference, y = GO, color = Ttest, size = Difference)) + 
  geom_point(stat = 'identity') + 
  xlab("") + ylab("path") + ggtitle("") + 
  theme_bw()

