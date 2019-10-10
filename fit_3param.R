# fit quadratic model to normalized, log-transformed data
# fitting 3 parameter: slope, slope^2, intercept
# use non-linear model, fit to all aa (add small value to zero counts)

library(tidyverse)
library(stringr)
library(cowplot)
library(broom)
theme_set(theme_cowplot())

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

# # remove aa with zero values
# ordered_count %>%
#   filter(count > 0) -> observed_aa

# add a small constant to include aa with zero values
ordered_count %>%
  mutate(count = count + 0.000001) -> observed_aa
ordered_count %>%
  mutate(count_mod = count + 0.000001) -> ordered_count

# count relative to most frequent - rescale to 1 for highly conserved sites
observed_aa %>%
  group_by(site) %>%
  mutate(rel_count = count/max(count)) -> observed_aa

# log-transformed data
observed_aa %>%
  mutate(ln_count = log(rel_count)) -> observed_aa 

#observed_aa %>% filter(site == 75) %>% ggplot(aes(x=k, y=rel_count)) + geom_point()

# fit with two parameters 
#observed_aa %>% mutate(k2 = k*k) -> two_param

#ak+bk^2
#a(k-k0)^2; k0<1, a<= 0
# poly(k,2) not working
observed_aa %>% nest(-site) %>%
  mutate(fit = map(data, ~ lm(ln_count ~ k + I(k^2), data = .)),
         intercept = map_dbl(fit, ~ (.)$coefficients[1]),
         slope = map_dbl(fit, ~ (.)$coefficients[2]),
         sqr_slp = map_dbl(fit, ~ (.)$coefficients[3])) %>% 
  select(-data, -fit) -> fits_2param

# visualize the distribution (at a given site) with linear fits
ordered_count %>%
  group_by(site) %>%
  mutate(rel_count = count/max(count)) -> ordered_count
ordered_count %>%
  mutate(ln_count = log(rel_count)) -> ordered_count

ordered_count %>% 
  left_join(fits_2param, by = "site") %>% 
  mutate(est_dist = slope*exp(k*slope + k*k*sqr_slp + intercept)) -> ordered_count

# rescale to compare w raw count
ordered_count %>%
  group_by(site) %>%
  mutate(est_rel = est_dist/sum(est_dist)) %>%
  mutate(est_count = sum(count)*est_rel) -> ordered_count

# visual comparison
ordered_count %>%
  filter(site == 4) %>%
  ggplot(aes(x = k)) +
  geom_point(aes(y = count), fill = "grey", stat='identity') +
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
  mutate(chisq = map(data, ~ sum(((.$count_mod - .$est_count)^2)/.$est_count)),
         p_value = map_dbl(chisq, ~ pchisq(., 16, lower.tail=FALSE))) %>% 
  select(-data) %>%
  mutate(p.adjusted = p.adjust(p_value, method= "fdr")) -> chisq_results

chisq_results %>%
  mutate(result = "pass") -> chisq_results

chisq_results$result[chisq_results$p.adjusted < 0.05] <- "fail"

# REMOVE 
chisq_results %>% 
  filter(chisq != "NA") -> chisq_results
chisq_results %>%
  filter(chisq != "NaN") -> chisq_results

chisq_results %>% filter(result == "fail") %>% nrow()/nrow(chisq_results)

chisq_results$chisq <- as.numeric(chisq_results$chisq)
ggplot(chisq_results, aes(chisq)) + geom_density() 

# chisq_results %>% mutate(method = "quadratic") -> chisq_results
# chisq_results_reg %>% mutate(method = "linear") -> chisq_results_reg
# full_join(chisq_results, chisq_results_reg) -> compare_chi
# 
# ggplot(compare_chi, aes(chisq, fill = method)) + geom_density(alpha = 0.5) + xlim(0,600)
# 
# full_join(compare_chi,eff_aa_all) -> compare_chi_eff
# compare_chi_eff %>% filter(result != "NA") -> compare_chi_eff
# 
# ggplot(compare_chi_eff, aes(x=result, y=eff_aa, fill=method)) + geom_violin(alpha=0.5) #+ coord_flip()
# 
# chisq_results %>% filter(result == "fail") %>% nrow()/nrow(chisq_results)
# chisq_results_reg %>% filter(result == "fail") %>% nrow()/nrow(chisq_results_reg)

ordered_count %>%
  filter(site == 1) %>%
  ggplot(aes(x = k)) +
  geom_bar(aes(y = count, fill = "actual"), fill = "grey", stat='identity') +
  geom_point(aes(y = est_count, color = "Estimated")) +
  geom_line(aes(y = est_count, color = "Estimated")) +
  scale_color_manual(" ", values = c("actual" = "grey", "Estimated" = "violetred"))  +
  scale_fill_manual(" ", values = "grey")  +
  theme(legend.position = "none") +
  labs(x = "k", y = "Count") +
  theme(text = element_text(size = 20),
        axis.text = element_text(size = 20)) +
  ggtitle("site = 1") #+ scale_y_log10()

#------- EFFECTIVE NUMBER OF AMINO ACIDS --------
# site_aa_count is actual dist. DOES NOT CONTAIN AA AT ZERO

# eff aa for actual distribution
observed_aa %>%
  group_by(site) %>%
  mutate(frequency = (count)/sum(count)) -> site_aa_freq

site_aa_freq %>% 
  mutate(flnf = frequency*log(frequency)) %>%
  group_by(site) %>% 
  summarize(entropy = -sum(flnf),
            eff_aa = exp(entropy),
            n = length(unique(aa))) -> eff_aa_all

# eff aa for estimated distribution
ordered_count %>% 
  group_by(site) %>%
  mutate(freq = (est_count)/sum(est_count)) -> site_aa_freq_2
site_aa_freq_2 %>% 
  mutate(flnf = freq*log(freq)) %>%
  group_by(site) %>%
  summarize(entr = -sum(flnf),
            eff_aa_est = exp(entr)) -> eff_aa_all_est

left_join(eff_aa_all, eff_aa_all_est) -> compare_eff_aa

ggplot(compare_eff_aa, aes(x=eff_aa, y=eff_aa_est)) + geom_point() +
  geom_abline(slope = 1, intercept = 0) + xlim(0,16) + ylim(0,16) + 
  labs(title = "3 parameter fit", x = "eff_aa actual", y = "eff_aa fit") #-> plot3



#------- RELATIONSHIP BETWEEN GAMMA & EFF AA -------
# ISSUE: need to visualize this relationship as a surface, because 

fits_2param %>%
  left_join(eff_aa_all) %>% 
  select(site, slope, sqr_slp, eff_aa)-> gam_v_effaa

gam_v_effaa %>%
  left_join(chisq_results) %>% 
  select(site, slope, sqr_slp, eff_aa, result) -> gam_v_effaa


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


