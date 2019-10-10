library(tidyverse)

neff <- function(lambda) {
  pi <- exp(-lambda*(0:19))
  C <- sum(pi)
  pi <- pi/C
  exp(-sum(pi*log(pi)))
}

df <- data.frame(lambda = (0:30)/10, ne = vapply((0:30)/10, neff, numeric(1)))
ggplot(df, aes(ne, -1/lambda)) + geom_point() + xlim(0, 12.5)
