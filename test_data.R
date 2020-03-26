library(tidyverse)
library(formula.tools)
set.seed(10)
n = 1000
p = 10
k = 1:p
beta = c(1, -1, rep(0, p-2))
x = matrix(rnorm(n*p), ncol = p)
colnames(x) = paste0("X", 1:p)
y = rbinom(n = n, size = 1, prob = expit(x %*% beta))

data = data.frame(x, y)

test = run_everything(data = data,
                      glmFamily = binomial(link = "logit"),
                      saveName = "test")

# Param

test = model_space(data = data,
                   B = 100,
                   m = 85)

MPP = mosaic_plot(test)

mppPlot = ggplot(MPP$plotData, aes(x = `Variable X`,
                          y = `Variable Y`)) +
  geom_tile(aes(fill = Probability)) +
  ggtitle(label = "MPP using subtractive measures") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_fill_continuous(limits=c(0, 1), breaks=seq(0,1,by=0.25),low = "white", high = "black")

VIP = variable_inclusion_plot(test)

vipPlot = ggplot(VIP$plotData, aes(x = Dimension,
                             y = Probability,
                             color = Variable)) +
  geom_line(aes(group = Variable)) +
  geom_point(aes(color = Variable)) +
  ggtitle(label = "VIP using subtractive measures") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
