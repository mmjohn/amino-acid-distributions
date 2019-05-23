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

# fit to a linear function, using map and BROOM 
# SET INTERCEPT TO 0
observed_aa %>% nest(-site) %>%
  mutate(fit = map(data, ~ lm(ln_count ~ 0 + k, data = .)),
         #intercept = map_dbl(fit, ~ (.)$coefficients[1]),
         slope = map_dbl(fit, ~ (.)$coefficients[1])) %>% 
  select(-data, -fit) -> observed_fits

# data frame with slope from linear fit to log transformed data - slope = gamma, intercept = 0
#head(observed_fits)

# visualize the distribution (at a given site) with linear fits
# NOTE: sites that are highly conserved for a single aa still look exponential in log scale (see sites 3 and 9 as example) 3, 6, 33
# library(forcats)
# observed_aa %>%
#   filter(site==3) %>%
#   ggplot() + 
#   geom_point(aes(x = fct_inorder(aa), y = ln_count)) +
#   geom_abline(intercept = 0, slope = -1.9044913) + # NOTE: look up slope from ordered_fits dataframe
#   ggtitle("site = 3") +
#   xlab("ordered aa") +
#   ylab("ln(count)") #-> p1
# plot_grid(p1, p2, p3, p4, ncol = 2) -> p
# title2 <- ggdraw() + 
#   draw_label("Linear fit to log-transformed distribution",
#              fontface = 'bold')
# plot_grid(title2, p, ncol = 1, rel_heights = c(0.1, 1))


# estimate the distribution
# only for observed
# observed_aa %>%
#   left_join(observed_fits, by = "site") %>%
#   mutate(est_dist = exp(k*slope)) -> observed_aa

ordered_count %>% 
  left_join(observed_fits, by = "site") %>%
  mutate(est_dist = exp(k*slope)) -> ordered_count

# rescale to compare w raw count
ordered_count %>%
  group_by(site) %>%
  mutate(est_rel = est_dist/sum(est_dist)) %>%
  mutate(est_count = sum(count)*est_rel) -> ordered_count

# visual comparison
ordered_count %>%
  filter(site == 3) %>%
  ggplot(aes(x = k)) +
  geom_point(aes(y = count, color = "steelblue")) +
  geom_point(aes(y = est_count, color = "violetred")) +
  scale_color_manual(name = "distribution",  
                     values = c("steelblue" = "steelblue", "violetred" = "violetred"), 
                     labels = c("actual", "estimated")) +
  ggtitle("site = 3") #-> q1

# ordered_count %>%
#   filter(site == 113) %>%
#   ggplot(aes(x = k)) +
#   geom_point(aes(y = count, color = "steelblue")) +
#   geom_point(aes(y = est_count, color = "violetred")) +
#   scale_color_manual(name = "distribution",  
#                      values = c("steelblue" = "steelblue", "violetred" = "violetred"), 
#                      labels = c("actual", "estimated")) +
#   ggtitle("site = 113") #-> q2
# 
# ordered_count %>%
#   filter(site == 33) %>%
#   ggplot(aes(x = k)) +
#   geom_point(aes(y = count, color = "steelblue")) +
#   geom_point(aes(y = est_count, color = "violetred")) +
#   scale_color_manual(name = "distribution",  
#                      values = c("steelblue" = "steelblue", "violetred" = "violetred"), 
#                      labels = c("actual", "estimated")) +
#   ggtitle("site = 33") #-> q3
# 
# ordered_count %>%
#   filter(site ==6) %>%
#   ggplot(aes(x = k)) +
#   geom_point(aes(y = count, color = "steelblue")) +
#   geom_point(aes(y = est_count, color = "violetred")) +
#   scale_color_manual(name = "distribution",  
#                      values = c("steelblue" = "steelblue", "violetred" = "violetred"), 
#                      labels = c("actual", "estimated")) +
#   ggtitle("site = 6") #-> q4

# q <- plot_grid(q1 + theme(legend.position = "none"),
#                q2 + theme(legend.position = "none"), 
#                q3 + theme(legend.position = "none"), 
#                q4 + theme(legend.position = "none"), 
#                ncol = 2)
# legend <- get_legend(q1)
# title <- ggdraw() + 
#   draw_label("Actual vs. estimated distributions",
#              fontface = 'bold')
# plot_grid(title, q, ncol = 1, rel_heights = c(0.1, 1)) -> tt
# plot_grid(tt, legend, rel_widths = c(2, .3))


#------- CHI-SQUARED: ACTUAL VS. ESTIMATED DIST -------
# SOMETHING WRONG HERE - SEE _diff.R
# use raw count
ordered_count %>% nest(-site) %>%
  mutate(chisq = map(data, ~ sum((.$count - .$est_count)/.$est_count)),
         p_value = map_dbl(chisq, ~ pchisq(., 18, lower.tail=FALSE))) %>% 
  select(-data) -> chisq_results

chisq_results %>%
  mutate(result = "pass") -> chisq_results

chisq_results$result[chisq_results$p_value < 0.05] <- "fail"

chisq_results %>% filter(result == "fail") %>% nrow()/nrow(chisq_results)

# test manually for a few sites
df1 <- filter(ordered_count, site == 1)
chisq1 <- sum((df1$count - df1$est_count)/df1$est_count)
pchisq(chisq1, 18, lower.tail = FALSE)


#------- EFFECTIVE NUMBER OF AMINO ACIDS --------
# site_aa_count is actual dist. DOES NOT CONTAIN AA AT ZERO

# site_aa_count %>% 
#   group_by(site) %>%
#   mutate(frequency = (count)/sum(count)) -> site_aa_freq

observed_aa %>%
  group_by(site) %>%
  mutate(frequency = (count)/sum(count)) -> site_aa_freq

site_aa_freq %>% 
  mutate(flnf = frequency*log(frequency)) %>%
  group_by(site) %>% 
  summarize(entropy = -sum(flnf),
            eff_aa = exp(entropy),
            n = length(unique(aa))) -> eff_aa_all

head(eff_aa_all)

eff_aa_all %>% 
  ggplot(aes(x=site, y=eff_aa)) + geom_point() + 
  #geom_smooth(se=F) +
  ylab("effective number of amino acids") +
  scale_x_continuous(limits = c(0, 154), expand = c(0, 0)) + 
  scale_y_continuous(limits = c(0,13), expand = c(0, 0)) + 
  theme_classic()

#------- RELATIONSHIP BETWEEN GAMMA & EFF AA -------
observed_fits %>%
  left_join(eff_aa_all) %>% 
  select(site, slope, eff_aa)-> gam_v_effaa

gam_v_effaa %>%
  left_join(chisq_results) %>% 
  select(site, slope, eff_aa, result) -> gam_v_effaa

# theoretical prediction
neff <- function(lambda) {
  pi <- exp(-lambda*(0:19))
  C <- sum(pi)
  pi <- pi/C
  exp(-sum(pi*log(pi)))
}


df_theory <- data.frame(lambda = (0:100)/10, ne = vapply((0:100)/10, neff, numeric(1)))

ggplot(gam_v_effaa, aes(x = eff_aa, y = 1/slope)) + 
  geom_point(aes(color = result)) +
  geom_line(data = df_theory, aes(ne, -1/lambda), inherit.aes = FALSE) 

# eff_aa_all %>%
#   left_join(gam_v_effaa, by = "site") -> df_plot
# 
# ggplot(gam_v_effaa, aes(x = eff_aa, y = 1/slope)) + 
#   geom_point(aes(color = result)) +
#   geom_point(data = data.frame(x = 1, y = 0), aes(x = x, y = y)) +
#   geom_line(data = df_plot, aes(x = (exp(-n)/eff_aa.x), y = 1/slope))

