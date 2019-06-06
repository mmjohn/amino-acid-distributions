# code to test the null model with a quadratic fit. 
# extending to fit all 20 aa, even if they have a zero count observed.
# 
library(tidyverse)
library(stringr)
library(cowplot)
library(broom)

# list of conditions: case_when()
# cut() to make intervals - instead of 19 ifelse statements

# generate data into tidy_align df
# correct method: multinomial with expo prob function
# p(k)=ce^(-lambda*k)
# rmultinom(10, 500, exp(-1*(0:19)))
# range of lambda, range of eff n
null0 <- data.frame(t(rmultinom(10, 500, exp(-1*(0:19))))) # 10 sites with 500 aa counts
null1 <- data.frame(t(rmultinom(10, 500, exp(-0.83*(0:19)))))
null2 <- data.frame(t(rmultinom(10, 500, exp(-0.75*(0:19)))))
null3 <- data.frame(t(rmultinom(10, 500, exp(-0.63*(0:19)))))
null4 <- data.frame(t(rmultinom(10, 500, exp(-0.5*(0:19)))))
null5 <- data.frame(t(rmultinom(10, 500, exp(-0.42*(0:19)))))
null6 <- data.frame(t(rmultinom(10, 500, exp(-0.35*(0:19)))))
null7 <- data.frame(t(rmultinom(10, 500, exp(-0.3*(0:19)))))
null8 <- data.frame(t(rmultinom(10, 500, exp(-0.25*(0:19)))))
null9 <- data.frame(t(rmultinom(10, 500, exp(-0.2*(0:19)))))
null10 <- data.frame(t(rmultinom(10, 500, exp(-0.15*(0:19)))))

full_join(null0, null1) %>% full_join(., null2) %>% full_join(., null3) %>%
  full_join(., null4) %>% full_join(., null5) %>% full_join(., null6) %>%
  full_join(., null7) %>% full_join(., null8) %>% full_join(., null9) %>%
  full_join(., null10) -> null

#---------- NULL DISTRIBUTION ----------
# recreate all_count df (site, aa, count)
null %>% mutate(site = seq(1:nrow(null))) %>% gather(aa, count, X1:X20) -> all_count


#----------- ESTIMATED DISTRIBUTION ------------
# USING ONLY OBSERVED AA
# order aa by relative frequency
all_count %>%
  arrange(all_count$site, desc(count)) -> ordered_count

# make the categorical variable (aa) numerical (k) by numbering aa's 1-20 in order of freq
ordered_count %>%
  group_by(site) %>%
  mutate(k = as.numeric(1:20)) -> ordered_count

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

# fit with two parameters 
observed_aa %>% nest(-site) %>%
  mutate(fit = map(data, ~ lm(ln_count ~ 0 + k + I(k^2), data = .)),
         slope = map_dbl(fit, ~ (.)$coefficients[1]),
         sqr_slp = map_dbl(fit, ~ (.)$coefficients[2])) %>% 
  select(-data, -fit) -> fits_2param

# visualize the distribution (at a given site) with linear fits
ordered_count %>%
  group_by(site) %>%
  mutate(rel_count = count_mod/max(count_mod)) -> ordered_count
ordered_count %>%
  mutate(ln_count = log(rel_count)) -> ordered_count

ordered_count %>% 
  left_join(fits_2param, by = "site") %>% 
  mutate(est_dist = slope*exp(k*slope + k*k*sqr_slp)) -> ordered_count

# rescale to compare w raw count
ordered_count %>%
  group_by(site) %>%
  mutate(est_rel = est_dist/sum(est_dist)) %>%
  mutate(est_count = sum(count_mod)*est_rel) -> ordered_count

# visual comparison
ordered_count %>%
  filter(site == 1) %>%
  ggplot(aes(x = k)) +
  geom_bar(aes(y = count_mod), fill = "grey", stat='identity') +
  geom_point(aes(y = est_count, color = "violetred")) +
  geom_line(aes(y = est_count, color = "violetred"), size = 1) +
  scale_color_manual(name = "distribution",
                     values = c("grey" = "grey", "violetred" = "violetred"),
                     labels = c("actual", "estimated"))  +
  labs(x = "amino acids ranked, k", y = "count") +
  theme(legend.position = "none")

# ordered_count %>%
#   group_by(site) %>%
#   filter(count != 0) %>%
#   mutate(num = n()) %>% 
#   filter(num != 1) -> ordered_count_variable

#------- CHI-SQUARED: ACTUAL VS. ESTIMATED DIST -------
ordered_count %>% nest(-site) %>%
  mutate(chisq = map(data, ~ sum(((.$count_mod - .$est_count)^2)/.$est_count)),
         p_value = map_dbl(chisq, ~ pchisq(., 18, lower.tail=FALSE))) %>% 
  select(-data) %>%
  mutate(p.adjusted = p.adjust(p_value, method= "fdr")) -> chisq_results

chisq_results %>%
  mutate(result = "pass") -> chisq_results

chisq_results$result[chisq_results$p_value < 0.05] <- "fail"

chisq_results %>% filter(result == "fail") %>% nrow()/nrow(chisq_results)


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
  labs(title = "Null model quadratic fit", x = "eff_aa actual", y = "eff_aa fit")

# eff_aa_all %>% 
#   ggplot(aes(x=site, y=eff_aa)) + geom_point() + 
#   #geom_smooth(se=F) +
#   ylab("effective number of amino acids") +
#   scale_x_continuous(limits = c(0, 154), expand = c(0, 0)) + 
#   scale_y_continuous(limits = c(0,13), expand = c(0, 0)) + 
#   theme_classic()

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
  pi <- lambda*exp(-lambda*(0:19))
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
  #ggtitle("1B4T_A") +
  #ylim(-10,1000) ++
  labs(x = "effective number of amino acids", y = "1/lambda") +
  theme(legend.position = "right",
        text = element_text(size = 40),
        axis.text = element_text(size = 40))  #-> plot1
