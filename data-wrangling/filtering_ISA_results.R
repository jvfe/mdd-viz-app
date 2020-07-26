library(dplyr)
library(magrittr)
library(stringr)

# First load all the objects from ISA and rename them according to each region.
# Couldn't figure out a way to do that in a clear and automated way, so I did it manually.

# Contains all the switchAnalyzeR objects
region_list <- list(ains, cg25, dlpfc, nac, ofc, sub)

topswitch <- lapply(region_list, function(x) {
  extractTopSwitches(x, n = Inf)
}) %>%
  bind_rows()

# The feature list contains some info we may want to use, such isoform ID and isoform fractions
features_for_each <- lapply(region_list, function(x) {
  x[["isoformFeatures"]]
}) %>%
  bind_rows() %>%
  mutate(comparison = paste(condition_1, condition_2)) %>%
  filter(isoform_switch_q_value < 0.05) %>%
  select(isoform_id, gene_id, IF1, IF2, comparison)


topswitch %>%
  tidyr::separate(condition_1, c("region", "cond1_phen", "cond1_sex"), remove = FALSE) %>%
  tidyr::separate(condition_2, c(NA, "cond2_phen", "cond2_sex"), remove = FALSE) %>%
  group_by(region) %>%
  filter(cond1_phen != cond2_phen & cond1_sex == cond2_sex) %>%
  ungroup() %>%
  select(-c(cond1_phen, cond2_phen, cond2_sex)) %>%
  mutate(comparison = paste(condition_1, condition_2), sex = cond1_sex) %>%
  left_join(features_for_each, by = c("gene_id", "comparison")) %>%
  select(-c(cond1_sex, comparison)) %>%
  readr::write_csv(path = "./data-wrangling/total_results_isa.csv")
