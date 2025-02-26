library(tidyverse)
setwd("/Users/kumarr9/Downloads")

#data <- data.frame(
#  id=seq(1,60),
#  individual=paste( "Mister ", seq(1,60), sep=""),
#  value=sample( seq(10,100), 60, replace=T)
#)

# Create dataset
#data <- data.frame(
#  id=seq(1,88),
#  individual=paste( "CL0083_T1R_T2" ,"CL0106_T2R_T" ,"CL0110_T2R_T" ,"CL0114_T1R_T" ,"CL0125_T4R_T" ,"CL0146_T2R_T" ,"CL0147_T1R_T" ,"CL0147_T3R_T" ,"CL0150_T1R_T" ,"CL0157_T1R_T" ,"CL0176_T1R_T" ,"CL0193_T1R_T" ,"fs_19_7842" ,"GS19_0000505B_RNA" ,"sb_18_3486" ,"SB19_1645_RNA" ,"sb_19_5365" ,"sb_19_6270" ,"SB19_858_1A_RNA" ,"SP_18_04378" ,"urmc_13" ,"urmc_15AAFS" ,"urmc_25B" ,"urmc_33" ,"urmc_7" ,"urmc_8" ,"CL0108_T2R_T" ,"CL0111_T2R_T2" ,"CL0116_T1R_T" ,"CL0116_T2R_T" ,"CL0124_T1R_T" ,"CL0191_T4R_T" ,"CL0214_T1R_T" ,"CL0267_T1R_T" ,"GS19_0004262_RNA" ,"NCI0422_T2R_T" ,"RS_19_5880" ,"sb_18_3090" ,"SB19_1389_RNA" ,"SB19_368_RNA" ,"sb_19_4537" ,"sb_19_4900" ,"SS19_503_RNA" ,"urmc_16A" ,"urmc_18B" ,"urmc_25" ,"urmc_31B" ,"urmc_35" ,"urmc_9" ,"CL0116_T4R_T" ,"CL0164_T1R_T" ,"CL0172_T1R_T" ,"NCI0402_T1R_T" ,"SB_11_6310" ,"sb_19_4812" ,"CL0108_T1R_T" ,"CL0152_T1R_T" ,"CL0191_T2R_T" ,"GS19_0000505A_RNA" ,"NCI0402_T2R_T" ,"SB19_2197_RNA" ,"SB_19_2276" ,"sb_19_3435" ,"SB19_858_1B_RNA" ,"urmc_31A" ,"urmc_5B" ,"CL0124_T2R_T" ,"CL0186_T4R_T" ,"ss_19_2398" ,"CL0107_T2R_T" ,"CL0170_T1R_T" ,"CL0191_T3R_T" ,"SB_19_3892" ,"CL0109_T2R_T" ,"CL0186_T3R_T" ,"CL0191_T1R_T" ,"SS19_1542_RNA" ,"urmc_18A" ,"urmc_24B" ,"CL0126_T1R_T" ,"CL0126_T2R_T" ,"CL0169_T1R_T" ,"sb_19_6623" ,"CL0174_T1R_T" ,"sb_19_5931" ,"CL0173_T1R_T" ,"NCI0422_T3R_T" ,"CL0196_T1R_T", seq(1,88), sep=","),
#  value=sample( seq(4,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,4,5,5,5,6,6,6,6,7,7,7,7,7,8,9,10,10,10,11,11,13,18,27), 88, replace=T)
#)
data <- read_tsv("test.tsv")
### test.csv in for_stacked_bar_graph is used here
# ----- This section prepare a dataframe for labels ---- #
# Get the name and the y position of each label
label_data <- data

# calculate the ANGLE of the labels
number_of_bar <- nrow(label_data)
angle <-  90 - 360 * (label_data$id-0.5) /number_of_bar     # I substract 0.5 because the letter must have the angle of the center of the bars. Not extreme right(1) or extreme left (0)

# calculate the alignment of labels: right or left
# If I am on the left part of the plot, my labels have currently an angle < -90
label_data$hjust<-ifelse( angle < -90, 1, 0)

# flip angle BY to make them readable
label_data$angle<-ifelse(angle < -90, angle+180, angle)
# ----- ------------------------------------------- ---- #


# Start the plot
p <- ggplot(data, aes(x=as.factor(id), y=value)) +       # Note that id is a factor. If x is numeric, there is some space between the first bar
  
  # This add the bars with a blue color
  geom_bar(stat="identity", fill=alpha("red", 0.7)) +
  
  # Limits of the plot = very important. The negative value controls the size of the inner circle, the positive one is useful to add size over each bar
  ylim(-100,120) +
  
 #Custom the theme: no axis title and no cartesian grid
 theme_minimal() +
 theme(
   axis.text = element_blank(),
  axis.title = element_blank(),
   panel.grid = element_blank(),
   plot.margin = unit(rep(-1,4), "cm")      # Adjust the margin to make in sort labels are not truncated!
 ) +
  
  # This makes the coordinate polar instead of cartesian.
  coord_polar(start = 0) +
  
  # Add the labels, using the label_data dataframe that we have created before
  geom_text(data=label_data, aes(x=id, y=value+10, label=individual, hjust=hjust), color="black", fontface="bold",alpha=0.6, size=2.5, angle= label_data$angle, inherit.aes = FALSE ) 

p

