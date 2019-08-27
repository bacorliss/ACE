
# Load packages
library(readr)
library(dplyr)
library(ggplot2)
library(tibble)


# Sample two null distributions
null_dists = tibble(x1 = rnorm(10, 0, 1), x2 = rnorm(10, 0, 1))



boxplot(null_dists,x1,x2)



# Basic box plot
p <- ggplot(null_dists, aes(x=x1, y=x2)) + 
  geom_boxplot() + geom_jitter()
p

# Box plot with dot plot
p + geom_dotplot(binaxis='y', stackdir='center', dotsize=1)
