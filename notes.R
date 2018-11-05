?pchisq
rnorm(10)
hist(rnorm(100))
dnorm(0)
x <- seq(-2, 2, length.out = 100)
plot(x, dnorm(x))
plot(x, pnorm(x))
x2 <- seq(0, 1, length.out = 100)
plot(x2, qnorm(x2))
qnorm(.95)
?rchisq
pchisq(10, 8, lower.tail = FALSE)
pchisq(30, 8, lower.tail = FALSE)
m <- cbind(c(3, 6), c(5, 4))
chisq.test(m, correct = FALSE)
chisq.test(m, correct = TRUE)
pchisq(.225, 1, lower.tail = FALSE)
pchisq(.9, 1, lower.tail = FALSE)
o <- c(3, 5, 6, 4)
e <- c(8*9/18, 8*9/18, 10*9/18, 10*9/18)
(o-e)^2/e
sum((o-e)^2/e)


df <- filter(ordered_count, site == 55)
chisq <- sum((df$count - df$est_count)/df$est_count)
pchisq(chisq, 18, lower.tail = FALSE)


df %>% ungroup() %>% select(k, count, est_count) %>%
  gather(distribution, value, -k) %>%
  ggplot(aes(k, value, color = distribution)) + geom_point()
