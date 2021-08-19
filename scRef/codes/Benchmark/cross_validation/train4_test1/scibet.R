suppressMessages(library(ggplot2))
suppressMessages(library(tidyverse))
suppressMessages(library(scibet))
suppressMessages(library(viridis))
suppressMessages(library(ggsci))

path_da <- "/home/drizzle_zhang/scRef/summary/scibet/test.rds.gz"
expr <- readr::read_rds(path = path_da) 

tibble(
    ID = 1:nrow(expr),
    label = expr$label
) %>%
    dplyr::sample_frac(0.7) %>%
    dplyr::pull(ID) -> ID

train_set <- expr[ID,]      #construct reference set
test_set <- expr[-ID,]      #construct query set

prd <- SciBet(train_set, test_set[,-ncol(test_set)])