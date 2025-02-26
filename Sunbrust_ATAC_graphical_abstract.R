#https://plotly.com/r/sunburst-charts/
#https://medium.com/r-evolution/visualizing-hierarchical-data-with-sunburst-charts-in-r-22b101f0ebfc
library(plotly)
fig <- plot_ly(
  labels = c("SCLC", "Neuroendocrine", "Non-neuroendocrine", "Stem Like Subtype", "NonNE"),
  parents = c("", "SCLC", "SCLC", "Non-neuroendocrine", "Non-neuroendocrine"),
  values = c(39, 24, 15, 11, 4),
  type = 'sunburst',
  branchvalues = 'total'
)

fig

