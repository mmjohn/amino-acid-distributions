library(tidyverse)
library(stringr)
library(cowplot)
library(broom)

# read in data
alignment <-read_csv("data/simulated/results_1B4T_A_evolved_split.csv", col_types = cols(.default = "c")) 

# tidy data used for null
alignment %>% 
  mutate(sequence = as.numeric(1:nrow(alignment)) ) %>% 
  gather(key = "site", value = "aa", 1:153) -> tidy_align  

tidy_align$site = as.numeric(tidy_align$site) 

#tidy_align %>% filter(sequence <= 100) -> tidy_align

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
ordered_count %>%
  group_by(site) %>%
  mutate(k = as.numeric(1:20)) -> ordered_count

# remove aa with zero values
ordered_count %>%
  filter(count > 0) -> observed_aa

# # count relative to most frequent - rescale to 1 for highly conserved sites
# observed_aa %>%
#   group_by(site) %>%
#   mutate(rel_count = count/max(count)) -> observed_aa
# 
# # log-transformed data
# observed_aa %>%
#   mutate(ln_count = log(rel_count)) -> observed_aa 

observed_aa %>% filter(site == 1) %>% ggplot(aes(x=k, y=count)) + geom_point()

# fit with exponential directly 
# non linear fitting, use param nls
# log and non log, and zeros or not
# fitting vs distributions
observed_aa %>% nest(-site) %>%
  mutate(fit = map(data, ~ lm(log(count) ~ k, data = .)),
         slope = map_dbl(fit, ~ (.)$coefficients[2]),
         intcp = map_dbl(fit, ~ (.)$coefficients[1])) %>% 
  select(-data, -fit) -> fits_expo

# approximate
ordered_count %>% 
  left_join(fits_expo, by = "site") %>% 
  mutate(est_dist = slope*exp(k*slope)) -> ordered_count

# rescale to compare w raw count
ordered_count %>%
  group_by(site) %>%
  mutate(est_rel = est_dist/sum(est_dist)) %>%
  mutate(est_count = sum(count)*est_rel) -> ordered_count

# visual comparison
ordered_count %>%
  filter(site == 1) %>%
  ggplot(aes(x = k)) +
  geom_bar(aes(y = count), fill = "grey", stat='identity') +
  geom_point(aes(y = est_count, color = "violetred")) +
  geom_line(aes(y = est_count, color = "violetred"), size = 1) +
  scale_color_manual(name = "distribution",
                     values = c("grey" = "grey", "violetred" = "violetred"),
                     labels = c("actual", "estimated"))  +
  labs(x = "amino acids ranked, k", y = "count") +
  theme(legend.position = "none")

# add small constant to est count for sites that have estimated zero counts
# ordered_count %>%
#   mutate(est_count_mod = est_count+0.00001) -> ordered_count



#------- CHI-SQUARED: ACTUAL VS. ESTIMATED DIST -------
# use raw count
ordered_count %>% nest(-site) %>%
  mutate(chisq = map(data, ~ sum(((.$count - .$est_count)^2)/.$est_count)),
         p_value = map_dbl(chisq, ~ pchisq(., 17, lower.tail=FALSE))) %>% 
  select(-data) %>%
  mutate(p.adjusted = p.adjust(p_value, method= "fdr")) -> chisq_results

chisq_results %>%
  mutate(result = "pass") -> chisq_results

chisq_results$result[chisq_results$p.adjusted < 0.05] <- "fail"

# REMOVE conserved sites with a single aa
chisq_results %>% 
  filter(chisq != "NA") -> chisq_results

chisq_results %>% filter(result == "fail") %>% nrow()/nrow(chisq_results)

#------- EFFECTIVE NUMBER OF AMINO ACIDS --------
# site_aa_count is actual dist. DOES NOT CONTAIN AA AT ZERO

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
fits_expo %>%
  left_join(eff_aa_all) %>% 
  select(site, slope, intcp, eff_aa)-> gam_v_effaa

gam_v_effaa %>%
  left_join(chisq_results) %>% 
  select(site, slope, intcp, eff_aa, result) -> gam_v_effaa


neff2 <- function(lambda_inv) {
  lambda <- 1/lambda_inv
  pi <- lambda*exp(-lambda*(0:19))
  C <- sum(pi)
  pi <- pi/C
  exp(-sum(pi*log(pi)))
}

df_theory2 <- data.frame(lambda_inv = (0:-100)/10, ne = vapply((0:-100)/10, neff2, numeric(1)))

ggplot(gam_v_effaa, aes(x = eff_aa, y = 1/slope)) + 
  geom_point(aes(color = result)) +
  #geom_line(data = df_theory, aes(ne, -1/lambda), inherit.aes = FALSE) +
  geom_line(data = df_theory2, aes(ne, lambda_inv), inherit.aes = FALSE, size = 1) +
  #ggtitle("1B4T_A") +
  #ylim(-10,1000) ++
  labs(x = "effective number of amino acids", y = "1/slope") +
  theme(legend.position = "right",
        text = element_text(size = 20),
        axis.text = element_text(size = 20),
        plot.title = element_text(size = 20))  #-> plot1


