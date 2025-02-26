update.packages("rlang")
install.packages("rlang")
devtools::install_github("jtextor/epitope-prediction")
library(EpitopePrediction)
fusion <- paste("MADQLYLENIDEFVTDQNKIVTYKWLSYTLGVHVNQAKQMLYDYVERKRKENSGAQLHVTYLVSGSLIQNGHSCHKVAVVREDKLEAVKSKLAVTASIHVYSIQKAMLKDSGPLFNTDYDILKSNLQNCSKFSAIQCAAAVPRAPAESSSSSKKFEQSHLHMSSETQANNELTTNGHGPPASKQVSQQPKGIMGMFASKAAAKTQETNKETKTEAKEVTNASAAGNKAPGKGNMMSNFFGKAAMNKFKVNLDSEQAVKEEKIVEQPTVSVTEPKLATPAGLKKSSKKAEPVKVLQKEKKRGKRVALSDDETKETENMRKKRRRIKLPESDSSEDEVFPDSPGAYEAESPSPPPPPSPPLEPVPKTEPEPPSVKSSSGENKRKRKRVLKSKTYLDGEGCIFIFRAPIKLSKPGELREEYESLRKLREEKLQEEKPSEDQIHKLLPEDTETGKRKMDEQKKRDEPLVLKTNLERCPARLSDSENEEPSRGQMTQTHRSAFVSKNNSYSLAFLAGKLNSKVERSQSCSDTAQERAKSRVRAVPGNKAKVTTMTPASNPIIGVLLSTQNNRCVSAPDLTIEKRLPFSSLSSLASLHKPERSVSPESNDSISEELNHFKPIVCSPCTPPKRLPDGRVLSPLIIKSTPRNLNRSLQKQTSYEASPRILKKWEQIFQERQIKKTLSKATLTSLAPEMGEELLGSEGIHSSKEKPLVAVNTRLSGGQVLSEYTGPTSADLDHFPSVSQTKAEQDSDNKSSTEIPLETCCSSELKGGGSGTSLEREQFEGLGSTPDAKLDKTCISRAMKITTVNSVLPQNSVLGGVLKTKQQLKTLNHFDLTNGVLVESLSEEPLPSLRRGRKRHCKTKHLEQNGSLKKLRQTSGEVGLAPTDPVLREMEQKLQQEEEDRQLALQLQRMFDNERRTVSRRKGSVDQYLLRSSNMAGAK",sep="")
fusion_new <- paste ("ESYSHTAEKTSKTQGAMRGISQENKEVMNPKTTIGVHELPDVEWQELEAPPKRSVYVAFATVTNGATCKSLILSLKPSMEGEKTNSHKALLQVGFVAGLDNAEVNASLDNTSRIRRIIEQEPIKARKQIIVSKRSLTPEKQGLQAVLNLILSEPSPTFETEHPSGESNDLELKRFSEREFTDPIEMSLTVKNAPKKKELTASNTKMRSDVASEPEESKREGVKARKQLNPCLTPCNAQTDELSSMSRSVSKARQLDGAKMKDQLKDDGTTSTKNDAPLKKTHAAMLRALPVVPPKLVQSSSSRTPQVRVSELLHVKQNLLEGSESKLPPSESNCKDMKLKRNLEVNAPDTSTREEKKLKDKGTVSGRPANRACLKEDSLCIPPSAKTSPHEVGHMKGVTYSNLQGRKGCEMKEVQLLMFNQDGSSVPSITNQLQQFTSNRYLVAPFSRRGKLKKKNSPNAKSNSQVTLTISTKQQFSKQLIHPFTFVKESKILSAPLEYLESQCPETYVSEVAHGIGLSSSQAVLASTGDFVGTEKEKKKPFLTLDCEAHTKSRKSLSYEKENPLSKEAGRLKRKTERDEHKCPGDELSKDSGPVAHLTTNNLLPLLLFDAEKLEPQESKHECIVGAQFLRRELKQKHNRNLSEYIEDENESSSASPFKGQADTGPKVASQTSPSEDSLGGASLVTETRKGTSGTSKVMQEGRSLPEMAGPSAKRKLLLRSWISELRTGKVQHAYPRTKSILRALTLVFKKGLVRENKIEKIPRPSMYLGDVHKQQKGKNFLDGPREGREMQKSPKLRLASECKRLYVVQEVSKNVYDRKIPIIAEQGTAGSESSNSPLQQQSQPNSDGPTEKKQDRLKEEATEMDQSLKAENPRVEERENVKYRVGLQSEKVSYEADNASSRPSFDAILPPMRAKLIPPCSVATSADSNQFPSYKAVN",sep="")
fusion_rand <- paste("MINCQQQCRRECYGMSPFHRCPDWKNLGSAWKVYKMIWSMFCAEYCFKEEVHVDYDNAQGNNYWPQVPGTVELMHHTCQYVWPSAFGALWAPVHLDELSICVNFAIQMYEVFFCWPGQKDRPWVFVRDHGERICMKSNWICNDIGGMTDNQICFNVERFDESKKHNYHAWQWGYKRVQFFHSGYAGGQFCSSAPVYYMWCTTCNTKKNAAWVTNLQWVFAEDNLTQWMDPCVSVGITSKSGFFEGSEGETIIKRFQGTGIWEKCHRVHANGSWVMHHPLCDDYRANYGAMTQFEDESANICFIKEGTEFELKGSYPGAFRQCITAKVKQYMIHGPPFTCCPNHPCFENTHTDINLLDDVLCWAATCCGEFECMDATAGCNIWWCCTRQYWMDDYKQGRYIIDYYREWDKIAYVTMSNIPNALVKGDLGICHFPHWWITHCDVGCFMIQINTKTTRHMEFVAADPHTMYWDMYTEHQFFQKLCESQWYQPDPWEKDQMVWIWRGAAFLHICSLWDIRHVWDPQGWPSRTEYFFEKYIQKLMHKGISIILKMLIPPAGMNRHKAMSDCWFMSEGCMFLYKKQFIVNLSHHLKTDWKGMIRAKHGYHSEEVPTERKQSSTLCHSTICADDALPEANVCIGYDVGEMTETMGVLMKIDGIICNTGANTYWHTACNYTEWYEHACWFQFPTYHVNTFQITRDFKPAHEMGPPMKKIWNRISRFNWRFQWSNTYRCDEWCQEQNCLERHLGDDSMVRNCEKNFWQGFQTKVMGDMHGDKLMDTPPRMCPWIMGPEDGSMTPQVLMGSMDEHFHLMTNVPDNRNVQIMWWWMCKKQWHTALTNRCYKMEWKRQMMGVYIEMAVDPFDPINDPDGEIWWCEYNIRTKIYKQPDFLLSGNSSMLSHKICDIRWRAQWHIEGDRLIQKSMTPGGDADDPWTRACNVPCDI", sep ="")
#fusion_junction <- paste("VLKSKTYLDGEGCIFIFRAPIKLSKPG",sep="")
binders( fusion, "HLA-A-02:01", 9 )
binders( fusion_new, "HLA-A-02:01", 9)
binders( fusion_rand, "HLA-A-02:01", 9)
#binders( fusion_junction, "HLA-A-02:01", 9 )
write.table(binders, file ="/Users/kumarr9/Downloads/epitope_binder.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE)
write.table(binders, file ="/Users/kumarr9/Downloads/epitope_binder.tsv", row.names = FALSE, col.names = TRUE, sep = '\t', append = FALSE);
### Epitope prediction for sequence given by Ajit ###
seq1 <- paste("SSGENKRKRKRVLKSKTYLDGEGCIFIFRAPIKLSKPGELREEYESLRKL",sep="")
seq2 <- paste("VKSSSGENKRKRKRVLKSKTYLDGEGCIFIFRAPIKLSKPGELREEYESLRKLRE",sep="")
seq3 <- paste("SVKSSSGENKRKRKRVLKSKTYLDGEGCIFIFRAPIKLSKPGELREEYESLRKLREE",sep="")
seq4 <- paste("QEEKPSEDQIHKLLPEDTETGKRKMDEQKKRDEPLVLKTNLERCPARLSDS",sep="")
seq5 <- paste("ASKAAAKTQETNKETKTEAKEVTNASAAGNKAPGKGNMMSNFFGKAAMNKFKVNLDSEQAVKEEKIVEQPTVSVTEPKLATPAGLKKSSKKAEPVKVLQKEKKRGKRVALSDDETKETENMRKKRRRIKLPESDSSEDEVFPDSPGAYEAESPSPPPPPSPPLEPVPKTEPEPPSVKSSSGENKRKRKRVLKSKTYLDGEGCIVTEKVYESESCTDSEEELNMKTSSVHRPPAMTVKKEPREERKGPKKGTAALGKANRQVSITGFFQRK",sep="")

binders( seq1, "HLA-A-02:01", 9 )
binders( seq2, "HLA-A-02:01", 9 )
binders( seq3, "HLA-A-02:01", 9 )
binders( seq4, "HLA-A-02:01", 9 )
binders( seq5, "HLA-A-02:01", 9 )

####
BiocManager::install("SMITE")
data("hg19_genes_bed")


##### FTD #####
setwd("/Users/kumarr9/Downloads")
library("dplyr")
#test <- read_tsv("FTD.hg38.neotype_selectvar.tsv")
test <- read.table(file = 'FTD.hg38.neotype_selectvar.tsv', sep = '\t', header = TRUE)
test1 <- test[,c(1,7:ncol(test))]
# test1 <- test[-2]
# test2 <- test1[-2]
# test3 <- test2[-2]
# test4 <- test3[-2]
# test5 <- test4[-2]
FTD <- subset(test1[,c(1,7)],test1[,7]==1)

###### ChIP-seq data analysis #####
##### CHIPseeker #####
BiocManager::install("ChIPseeker")
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
BiocManager::install("clusterProfiler")
library(clusterProfiler)
setwd("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks")

#### PDX samples ####
sample_2705000_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2705000_S0_L001_summits.bed")
sample_2705010_S40_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2705010_S40_L003_summits.bed")
sample_2705020_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2705020_S0_L001_summits.bed")
sample_2781020_S37_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2781020_S37_L003_summits.bed")
sample_2781080_S30_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2781080_S30_L003_summits.bed")
sample_2781090_S41_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2781090_S41_L003_summits.bed")
sample_2792010_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2792010_S0_L001_summits.bed")
sample_2792020_S32_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2792020_S32_L003_summits.bed")
sample_2792040_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2792040_S0_L001_summits.bed")
sample_2811610_S33_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2811610_S33_L003_summits.bed")
sample_2811630_S35_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2811630_S35_L003_summits.bed")
sample_2858490_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2858490_S0_L001_summits.bed")
sample_2858800_S38_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2858800_S38_L003_summits.bed")
sample_2858810_S31_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2858810_S31_L003_summits.bed")
sample_3013500_S34_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/3013500_S34_L003_summits.bed")
sample_3049380_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/3049380_S0_L001_summits.bed")
sample_3049400_S39_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/3049400_S39_L003_summits.bed")
sample_3049430_S42_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/3049430_S42_L003_summits.bed")
sample_3049450_S36_L003 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/3049450_S36_L003_summits.bed")

#### NCI-24-mouse-Chipseq ####

sample_RA17_L2a_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/RA17_L2a_1_summits.bed")
sample_RA19Li1a_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/RA19Li1a_1_summits.bed")
sample_RA21_16_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/RA21_16_1_summits.bed")
sample_RA22_12_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/RA22_12_1_summits.bed")
sample_RA23_3_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/RA23_3_1_summits.bed")
sample_s270500_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s270500_1_summits.bed")
sample_s270502_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s270502_1_summits.bed")
sample_s278102_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s278102_1_summits.bed")
sample_s278108_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s278108_1_summits.bed")
sample_s278109_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s278109_1_summits.bed")
sample_s279201_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s279201_1_summits.bed")
sample_s279204_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s279204_1_summits.bed")
sample_s281161_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s281161_1_summits.bed")
sample_s285881_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s285881_1_summits.bed")
sample_s301350_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s301350_1_summits.bed")
sample_s304940_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s304940_1_summits.bed")
sample_s304944_1 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/s304944_1_summits.bed")

##### NCI-8-mouse-ATAC #####

sample_270502_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/270502_S0_L001_summits.bed")
sample_278106_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/278106_S0_L001_summits.bed")
sample_304944_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/304944_S0_L001_summits.bed")
sample_304955_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/304955_S0_L001_summits.bed")
sample_313308_S0_L001 <- readPeakFile("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/313308_S0_L001_summits.bed")


#peak coverage plot ###
covplot(sample_304955_S0_L001, weightCol="V5")

#peakHeatmap("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2705000_S0_L001/2705000_S0_L001_summits.bed", TxDb=txdb, upstream=3000, downstream=3000, color="red")
peakHeatmap(sample_2705000_S0_L001, TxDb=txdb, upstream=3000, downstream=3000, color="red")

plotAvgProf2("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2705000_S0_L001/2705000_S0_L001_summits.bed", TxDb=txdb, upstream=3000, downstream=3000,conf = 0.95,
             xlab="Genomic Region (5'->3')", ylab = "Read Count Frequency")


plotPeakProf2(peak = peak, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body", nbin = 800,
              TxDb = txdb, weightCol = "V5",ignore_strand = F)

#peakAnno <- annotatePeak("/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks/2705000_S0_L001/2705000_S0_L001_summits.bed", tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")

ssample_313308_S0_L001 <- annotatePeak(sample_313308_S0_L001, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(sample_313308_S0_L001)
plotAnnoBar(peakAnno)
vennpie(peakAnno)
BiocManager::install("ggupset")
upsetplot(peakAnno)

##### 
files <- list.files(path="/Users/kumarr9/Downloads/PDX_atac_narrow/narrow_peaks")
print(files)
peak <- readPeakFile(files[[40]])

covplot(peak, weightCol="V5")
peakHeatmap(peak, TxDb=txdb, upstream=3000, downstream=3000, color="red")

peakAnno <- annotatePeak(peak, tssRegion=c(-3000, 3000),
                         TxDb=txdb, annoDb="org.Hs.eg.db")
#install.packages("ggimage")
library(ggimage)
upsetplot(peakAnno, vennpie=TRUE)

BiocManager::install("ReactomePA")
library(ReactomePA)

pathway1 <- enrichPathway(as.data.frame(peakAnno)$geneId)
head(pathway1, 2)

gene <- seq2gene(peak, tssRegion = c(-1000, 1000), flankDistance = 3000, TxDb=txdb)
pathway2 <- enrichPathway(gene)
head(pathway2, 2)
dotplot(pathway2)

plotPeakProf2(files, upstream = 3000, downstream = 3000, conf = 0.95,
              by = "gene", type = "start_site", TxDb = txdb,
              facet = "row", nbin = 800)

plotPeakProf2(peak, upstream = rel(0.2), downstream = rel(0.2),
              conf = 0.95, by = "gene", type = "body",
              TxDb = txdb, facet = "row", nbin = 800)

peakAnnoList <- lapply(files, annotatePeak, TxDb=txdb, tssRegion=c(-3000, 3000), verbose=FALSE)
#### sample dataset for checking merging function ####
x <- data.frame(k1=c(NA,NA,3,4,5), k2=c(1,NA,NA,4,5), data=1:5)
y <- data.frame(k1=c(NA,2,NA,4,5), k2=c(NA,NA,3,4,5), data=1:5)
t <- merge(x, y, by=c("k1","k2"))


##### 
promoter <- getPromoters(TxDb=txdb, upstream=3000, downstream=3000)
tagMatrixList <- lapply(files, getTagMatrix, windows=promoter)


########### merging two datasets based on column name "snp_id" #####
setwd ("/Users/kumarr9/Downloads")
library(dplyr)
HET <- read.table(file="HET_genotype_CDCFY18_19_21.tsv", sep = "\t", header=T)
Path <- read.table(file="Pathogenic.tsv", sep = "\t", header = T)
df <- HET %>% inner_join(Path, by=c("snp_id"))
write.table(df, file = "/Users/kumarr9/Downloads/rajesh.tsv", append = TRUE,sep = "\t",)

####### stacked bar ####
library(tidyverse)
ranks <- read_tsv("/Users/kumarr9/Downloads/ATAC.tsv", col_names = TRUE)
library(ggplot2)
ggplot(ranks, aes(fill=Splice_type, y=percent_value, x=Chromosome)) + 
  geom_bar(position="stack", stat="identity") + coord_flip () +theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                                                                                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))