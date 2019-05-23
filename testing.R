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

# count relative to most frequent - rescale to 1 for highly conserved sites
observed_aa %>%
  group_by(site) %>%
  mutate(rel_count = count/max(count)) -> observed_aa

# log-transformed data
observed_aa %>%
  mutate(ln_count = log(rel_count)) -> observed_aa 


observed_aa %>% 
  filter(site == 10) %>%
  ggplot(aes(x = k, y = ln_count)) + geom_point()


# fit to a linear function, using map and BROOM 
# SET INTERCEPT TO 0
observed_aa %>% group_by(site) %>% summarize(unique = nrow(.)) %>% View()

observed_aa %>% nest(-site) %>%
  mutate(fit = map(data, ~ lm(ln_count ~ 0 + poly(k,2), data = .)),
         slope = map_dbl(fit, ~ (.)$coefficients[1]),
         slope2 = map_dbl(fit, ~ (.)$coefficients[2])) %>% 
  select(-data, -fit) -> observed_fits


observed_aa %>%
  filter(site == 10) %>%
  ggplot() +
  geom_point(aes(k, ln_count)) +
  geom_smooth(aes(k, ln_count), formula = ln_count ~ poly(k,2), data = observed_aa)
  #geom_abline(intercept = 0, slope = -0.2821500)



observed_aa %>% filter(site == 10) -> test_df

lm(test_df$ln_count ~ 0 + poly(test_df$k,2))

test_df %>%
  mutate(poly_fit = exp(-4.824*k + 1.051*k*k),
         est_rel = poly_fit/sum(poly_fit),
         est_count = sum(count)*est_rel) -> test_df


test_df %>% ggplot() +
  #geom_point(aes(k, ln_count), color = "red") +
  geom_point(aes(k, est_count), color = "blue") #+
  ylim(-5,0)






lm(log(test_df$count) ~ poly(test_df$k,2))

test_df %>%
  mutate(poly_fit = exp(2.464-4.824*k + 1.051*k*k)) -> test_df


test_df %>% ggplot() +
  #geom_point(aes(k, count), color = "red") +
  geom_point(aes(k, poly_fit), color = "blue") 




# examp_site <- observed_aa %>% filter(site==10)
# examp_fit <- observed_fits %>% filter(site==10)
# 
# plot(log(examp_site$k),log(examp_site$count))
# abline(examp_fit$intercept, examp_fit$slope)
# 
# plot(examp_site$k,examp_site$count)
# #lines(examp_site$k, examp_site$k^examp_fit$slope, col = "red")
# lines(examp_site$k, 205*exp(examp_fit$slope*examp_site$k), col = "red")
# 
# examp_exp <- data_frame(k = seq(1:20))
# examp_exp %>% mutate(est = 405*exp(k*-1.512033)) -> examp_exp
# 
# ggplot() +
#   geom_point(aes(k, count), data = examp_site) +
#   geom_point(aes(k, est), data = examp_exp)


# visualize the distribution (at a given site) with linear fits
# NOTE: sites that are highly conserved for a single aa still look exponential in log scale (see sites 3 and 9 as example) 3, 6, 33
ordered_count %>%
  group_by(site) %>%
  mutate(rel_count = count/max(count)) -> ordered_count
ordered_count %>%
  mutate(ln_count = log(rel_count)) -> ordered_count

ordered_count %>% 
  left_join(observed_fits, by = "site") %>%
  mutate(est_dist = slope*exp(k*slope)) -> ordered_count

# rescale to compare w raw count
ordered_count %>%
  group_by(site) %>%
  mutate(est_rel = est_dist/sum(est_dist)) %>%
  mutate(est_count = sum(count)*est_rel) -> ordered_count

