library(tidyverse, help, pos = 2, lib.loc = NULL)
a <- rnorm(1000)
p <- ggplot(data = mpg, aes(x = cyl, y = hwy)) + geom_point()
show(p)
