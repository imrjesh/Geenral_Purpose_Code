install.packages("VennDiagram")
library("VennDiagram")
grid.newpage()
draw.triple.venn(area1  = 6460,
                 area2 = 64328,
                 area3 = 2065,
                 n12  = 847,
                 n23 = 573,
                 n13  = 401,
                 n123  = 376,
                 category =
                   rep("", 3),  rotation = 1, fill = c("Yellow", "Purple", "Green"),lty = "blank",cex = 2.0,
                 cat.cex = 2.0,)
