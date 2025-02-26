library(tidyverse)

TP <- read_tsv("tt.csv")
ggplot(TP, aes(x = Cancer_type, y = SCLC)) +
  geom_col()


TP <- mutate(TP, Extrapulmonary.0001 = Extrapulmonary/0001)
TP <- mutate(TP, Extrapulmonary.0001 = -Extrapulmonary.0001)
ggplot(TP, aes(x = Cancer_type)) +
  geom_col(aes(y = SCLC), fill = "#5d8402") +
  geom_col(aes(y = Extrapulmonary.0001), fill = "#817d79") 

ggplot(TP, aes(x = Cancer_type)) +
  geom_col(aes(y = SCLC), fill = "#5d8402") +
  geom_col(aes(y = Extrapulmonary.0001), fill = "#817d79") +
  coord_polar()

ggplot(TP, aes(x = reorder(Cancer_type, SCLC))) +
  geom_col(aes(y = SCLC), fill = "#5d8402") +
  geom_col(aes(y = Extrapulmonary.0001), fill = "#817d79") +
  coord_polar()

