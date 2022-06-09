library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

read_and_combine_data <- function(args){
    out <- vector('list', length(args))
    for (num in seq_along(args)){
        these_data <- read_tsv(args[num], col_types ="ci")
        these_data <- these_data %>% mutate(sample_name = basename(args[num]) %>% str_replace("_frags_per_contig.tsv", ""))
        out[[num]] <- these_data
    }
    bind_rows(out)
}

d <- read_and_combine_data(args)
write_tsv(d, "results/quality_control/frags_per_contig/all_samples.tsv")
