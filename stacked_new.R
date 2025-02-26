setwd("/Users/kumarr9/Desktop/rajesh_projects/for_stacked_bar_graph/new_stack_with_reorder")
TTM <- read.delim("stacked_ttt.txt")
library(ggplot2)

ggplot(data = TTM, aes(x = reorder(Sample, +Number), y = Number, fill = Disease)) +
  geom_bar(stat = "identity")
#theme(axis.text.x = element_text(angle = 110))

### to flip the graph
ggplot(data = TTM, aes(x = reorder(Sample, +Number), y = Number, fill = Number)) +
  geom_bar(stat = "identity", position="dodge") + coord_flip() + facet_wrap(~Disease, scales = "free")+ theme(axis.text.y = element_text(angle = 00)) +
  theme(axis.text.y = element_text(size =  06, face = "italic"))
### to change color
#ggplot(data = TTM, aes(x = Sample_ID, y = Fusion, fill = factor(Fusion_Type))) +
geom_bar(stat = "identity") + coord_flip() +  scale_fill_brewer(palette = 10) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                                                      panel.background = element_blank(), axis.line = element_line(colour = "black"))
