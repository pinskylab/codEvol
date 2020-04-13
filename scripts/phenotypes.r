# analyze phenotype changes

# functions
require(data.table)
require(ggplot2)

calc50p <- function(x,y){ # fit logistic regression and return x value for y=50%
  mod <- glm(y ~ x, family = binomial)
  betas <- coef(mod)
  return(-betas[1]/betas[2])
}

# read in data
canmat <- fread('data/phenotypes/Brattey et al_2018_Table18_proportion_mature.csv')

# long format
canmatl <- melt(canmat, id.vars = c('Year'), value.name = 'propmat')
canmatl[, age := as.numeric(gsub('Age', '', variable))]

# plot maturity curves by year
ggplot(canmatl, aes(age, propmat, group = Year, color = Year)) +
  geom_line(alpha = 0.5)

# calc 50% maturity from logistic regression models
mat50 <- canmatl[, .(age50 = calc50p(age, cbind(propmat*1000, (1-propmat)*1000))), by = Year]

# plot age of 50% maturity by year
ggplot(mat50, aes(Year, age50)) +
  geom_line()

# write out
write.csv(mat50, 'output/age_50percmature_can.csv', row.names = FALSE)
