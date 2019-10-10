library(tidyverse)
library(Biostrings)
library(stringr)
library(cowplot)
library(broom)
theme_set(theme_cowplot())

# import alignment
fasta_file <- readAAStringSet("data/real/5PTP_A_rp35.mafft.processed.afa", format = "fasta")
seq_name = names(fasta_file)
aa_seq = paste(fasta_file)
seq_num = seq(1:length(seq_name))
real_align <- data.frame(seq_num, seq_name, aa_seq, 
                         stringsAsFactors = F)

# change sequences from strings into character vectors
as.data.frame(strsplit(real_align$aa_seq, split = "")) -> int_df
as.data.frame(t(int_df)) -> int_df
row.names(int_df) <- NULL
names(int_df) <- substring(names(int_df),2)

# tidy data frame
int_df %>%
  mutate(num_seq = seq(1:nrow(int_df))) -> int_df
int_df %>%
  gather(key = "site", value = "aa", 1:221) -> tidy_align
tidy_align$site <- as.numeric(tidy_align$site)

#---------- ACTUAL DISTRIBUTION ----------
# determine the number of each AA at each site
tidy_align %>% 
  group_by(site, aa) %>% 
  summarize(count = n()) -> site_aa_count

empty_counts <- data_frame(aa = c('A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W', 'Y'), 
                           count = rep(0, 20))

# 1. First write a function that can take a data frame with amino acids at one site, and return a data frame with counts, including zeros.
count_by_site <- function(site_df) {
  site_df %>% 
    group_by(aa) %>% 
    summarize(count = n()) -> counts
  full_join(counts, anti_join(empty_counts, counts, by = 'aa'))
}

# 2. Then call that function for all sites and combine the resulting data frames.
tidy_align %>% nest(-site) %>%
  mutate(counts = map(data, count_by_site)) %>%
  dplyr::select(-data) %>%   # use sessionInfo() to find what packages are loaded, conflicted package to give warnings
  unnest() -> all_count

#----------- ESTIMATED DISTRIBUTION ------------
# USING ONLY OBSERVED AA
# order aa by relative frequency
all_count %>%
  arrange(all_count$site, desc(count)) -> ordered_count

# make the categorical variable (aa) numerical (k) by numbering aa's 1-20 in order of freq
# EXCLUDE all '-' entries/gaps in alignment
ordered_count %>%
  filter(aa != "-") -> no_gaps

no_gaps %>%
  group_by(site) %>%
  mutate(k = as.numeric(1:20)) -> no_gaps

# remove aa with zero values
no_gaps %>%
  filter(count > 0) -> observed_aa

# count relative to most frequent - rescale to 1 for highly conserved sites
observed_aa %>%
  group_by(site) %>%
  mutate(rel_count = count/max(count)) -> observed_aa

# log-transformed data
observed_aa %>%
  mutate(ln_count = log(rel_count)) -> observed_aa 

# fit to a linear function, using map and BROOM 
# SET INTERCEPT TO 0
observed_aa %>% nest(-site) %>%
  mutate(fit = map(data, ~ lm(ln_count ~ 0 + k, data = .)),
         #intercept = map_dbl(fit, ~ (.)$coefficients[1]),
         slope = map_dbl(fit, ~ (.)$coefficients[1])) %>% 
  select(-data, -fit) -> observed_fits

# number of observed amino acids at each site in the alignment
observed_aa %>% 
  group_by(site) %>%
  mutate(obs_count = sum(count)) %>% 
  select(site, obs_count) %>% distinct() -> align_count_bysite

align_count_bysite %>%
  ggplot(aes(x = site, y = obs_count)) +
  geom_point() +
  ylab("observed count") +
  theme(text = element_text(size = 35),
        axis.text = element_text(size = 35)) 
ggsave("gaps.jpg")

# estimate count
no_gaps %>% 
  left_join(observed_fits, by = "site") %>%
  mutate(est_dist = exp(k*slope)) -> no_gaps

# rescale to compare w raw count
no_gaps %>%
  group_by(site) %>%
  mutate(est_rel = est_dist/sum(est_dist)) %>%
  mutate(est_count = sum(count)*est_rel) -> no_gaps

#------- CHI-SQUARED: ACTUAL VS. ESTIMATED DIST -------
# use raw count
no_gaps %>% nest(-site) %>%
  mutate(chisq = map(data, ~ sum(((.$count - .$est_count)^2)/.$est_count)),
         p_value = map_dbl(chisq, ~ pchisq(., 18, lower.tail=FALSE))) %>% 
  select(-data) %>%
  mutate(p.adjusted = p.adjust(p_value, method= "fdr")) -> chisq_results

chisq_results %>%
  mutate(result = "pass") -> chisq_results

chisq_results$result[chisq_results$p_value < 0.05] <- "fail"

chisq_results %>% filter(result == "fail") %>% nrow()/nrow(chisq_results)


no_gaps %>% 
  filter(count != 0) %>% 
  group_by(site) %>%
  do(tidy(chisq.test(.$count, p = .$est_count/sum(.$est_count)))) %>%
  ungroup() %>%
  mutate(p.adjusted = p.adjust(p.value, method= "fdr")) -> observed_chi

observed_chi %>%
  mutate(result = "pass") -> observed_chi
observed_chi$result[observed_chi$p.value < 0.05] <- "fail"
observed_chi %>% filter(result == "fail") %>% nrow()/nrow(observed_chi)




#visual site that failed chi-squared
no_gaps %>%
  filter(site == 61) %>%
  ggplot(aes(x = k)) +
  geom_bar(aes(y = count), fill = "grey", stat='identity') +
  geom_point(aes(y = est_count, color = "violetred")) +
  geom_line(aes(y = est_count, color = "violetred"), size = 1) +
  scale_color_manual(name = "distribution",
                     values = c("grey" = "grey", "violetred" = "violetred"),
                     labels = c("actual", "estimated"))  +
  labs(x = "amino acids ranked, k", y = "count") +
  theme(legend.position = "none",
        text = element_text(size = 35),
        axis.text = element_text(size = 35))
ggsave("real_fail.jpg")

#------- EFFECTIVE NUMBER OF AMINO ACIDS --------
observed_aa %>%
  group_by(site) %>%
  mutate(frequency = (count)/sum(count)) -> site_aa_freq

site_aa_freq %>% 
  mutate(flnf = frequency*log(frequency)) %>%
  group_by(site) %>% 
  summarize(entropy = -sum(flnf),
            eff_aa = exp(entropy),
            n = length(unique(aa))) -> eff_aa_all


#------- RELATIONSHIP BETWEEN GAMMA & EFF AA -------
observed_fits %>%
  left_join(eff_aa_all) %>% 
  select(site, slope, eff_aa)-> gam_v_effaa

gam_v_effaa %>%
  left_join(chisq_results) %>% 
  select(site, slope, eff_aa, result) -> gam_v_effaa

# theoretical prediction
neff2 <- function(lambda_inv) {
  lambda <- 1/lambda_inv
  pi <- exp(-lambda*(0:19))
  C <- sum(pi)
  pi <- pi/C
  exp(-sum(pi*log(pi)))
}
#neff2(1/3)

df_theory2 <- data.frame(lambda_inv = (0:-100)/10, ne = vapply((0:-100)/10, neff2, numeric(1)))

ggplot(gam_v_effaa, aes(x = eff_aa, y = 1/slope)) + 
  geom_point(aes(color = result)) +
  #geom_line(data = df_theory, aes(ne, -1/lambda), inherit.aes = FALSE) +
  geom_line(data = df_theory2, aes(ne, lambda_inv), inherit.aes = FALSE, size = 1) +
  labs(x = "effective number of amino acids", y = "1/slope", title = "5PTP_A_rp35") +
  theme(legend.position = "none",
        text = element_text(size = 35),
        axis.text = element_text(size = 35),
        plot.title = element_text(size = 35)) #-> plot10
ggsave("real_fail_dist.jpg")
