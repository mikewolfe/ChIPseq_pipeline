library(tidyverse)
d <- read_csv("log_tpm_full.csv")
d <- d %>% rename(bnumber = '...1')
gene_info <- read_tsv("../../results/alignment/combine_bed/U00096.3_AL009126.3/U00096.3_AL009126.3_annotations.tsv")
d2 <- d %>% inner_join(gene_info %>% filter(chr == "U00096.3") %>% rename(bnumber = locus_tag), by = "bnumber")
dfinal <- d2 %>% select(chr, start, end, strand, unique_name, bnumber, control__wt_glc__1, control__wt_glc__2, oxidative__wt_pq__1, oxidative__wt_pq__2)

# Make b.sub genes have a similar distro of weights
avg_std <- dfinal %>% select(bnumber, control__wt_glc__1, control__wt_glc__2, oxidative__wt_pq__1, oxidative__wt_pq__2) %>% 
    pivot_longer(-bnumber, values_to = "exp", names_to = "sample") %>%
    group_by(bnumber) %>% summarize(avg = mean(exp)) %>% ungroup() %>%
    summarize(overall_avg = mean(avg), overall_std = sd(avg))
# simulate on log scale then transform to normal scale
bsub <- gene_info %>% filter(chr == "AL009126.3") %>% rename(bnumber = locus_tag) %>% select(chr, start, end, strand, unique_name, bnumber) %>%
    mutate(score = rnorm(n(), avg_std$overall_avg[[1]], avg_std$overall_std[[1]]))%>%
    mutate(score = if_else(score < 0, 0, score)) %>%
    mutate(weight = score/sum(score))

# Simulate IP and inputs for each sample

simulate <- function(df, bsub, outpre, frac_spike = 0.01, epitope_strength = 1){
    df <- df %>% mutate(strand = ".")
    total <- df %>% mutate(cov = (end-start)*score) %>% summarize(total = sum(cov))
    df %>% 
        # We want the spike-in to be 1% of the total stuff so we take the total value of the
        # of the enrichment + 1 for every possible position in E coli - 1 for every 
        # possible position in b sub. Multiply that by 1 % and distribute it equally
        # across the genome for B. sub
        bind_rows(bsub %>% transmute(chr, start, end, unique_name, 
                                     score = (total[[1]]*(frac_spike*(1/epitope_strength)))*weight/(end-start),
                                     strand = ".")) %>%
        write_tsv(str_c(outpre, ".bed"), col_names = FALSE)
    
    tibble(chr = c("U00096.3", "AL009126.3"), start = c(0, 0), end = c(4641652,4215606), unique_name = c("nonspike","spike"), 
        value = c(1/4641652, frac_spike/4215606), strand = c(".", ".")) %>% 
        write_tsv(str_c(outpre, "_inp.bed"), col_names = FALSE)
}

## sample 1
simulate(dfinal %>% select(chr, start, end, unique_name, score = control__wt_glc__1, strand), bsub, "wt_glc_rep1")

## sample 2
simulate(dfinal %>% select(chr, start, end, unique_name, score = control__wt_glc__2, strand), bsub, "wt_glc_rep2")

## sample 3
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1")

## sample 4
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__2, strand), bsub, "wt_pq_rep2")


## sample 5
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1_2fold", epitope_strength = 2)

## sample 6
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__2, strand), bsub, "wt_pq_rep1_halffold", epitope_strength = 0.5)

## sample 7
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1_spikehalf", frac_spike = 0.005)

## sample 8
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1_spikedouble", frac_spike = 0.02)

## sample 9
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1_spikefived", frac_spike = 0.05)

## sample 10
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1_spikehalf_2fold", frac_spike = 0.005, epitope_strength = 2)

## sample 11
simulate(dfinal %>% select(chr, start, end, unique_name, score = oxidative__wt_pq__1, strand), bsub, "wt_pq_rep1_spikedouble_2fold", frac_spike = 0.02, epitope_strength = 2)

