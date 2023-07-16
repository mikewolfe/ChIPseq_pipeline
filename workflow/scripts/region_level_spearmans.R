#' load in required libraries
suppressMessages(library(tidyverse))

# Check for correspondance between biological replicates across a set of
# regions

spearman_per_region <- function(d){
    samples <- colnames(d %>% select(-c(region,coord, chrm)))
    out <- list()
    k <- 1
    # do every combination of samples only once
    for(i in seq_along(samples)){
        for(j in seq_along(samples)){
            # skip the upper diagonal
            if(j <= i){
                next
            }
            # use syntax to get past tidy_eval issues
            samp1 <- samples[[i]]
            samp2 <- samples[[j]]
            out[[k]] <- d %>% group_by(region) %>%
                    summarize(cor = cor(!!sym(samp1), !!sym(samp2), method = "sp"),
                              samp1 = samples[[i]],
                              samp2 = samples[[j]],
                              avg_samp1 = mean(!!sym(samp1), na.rm = TRUE),
                              avg_samp2 = mean(!!sym(samp2), na.rm = TRUE))
        # increment our combo counter
        k <- k + 1
        }
    }
    # bind into one dataframe
    bind_rows(out)
} 

args <- commandArgs(trailingOnly = TRUE)

ind <- read_tsv(args[1])

out <- spearman_per_region(ind)

write_tsv(out, args[2])
