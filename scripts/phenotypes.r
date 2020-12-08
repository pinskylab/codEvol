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
lofmat <- fread('data/phenotypes/AFWG_Table 3.11.csv', skip = 3)
lofmat <- lofmat[!is.na(Year_age), ] # trim off an NA row

# long format
canmatl <- melt(canmat, id.vars = c('Year'), value.name = 'propmat')
canmatl[, age := as.numeric(gsub('Age', '', variable))]

setnames(lofmat, c('Year_age', '+gp'), c('Year', '15')) # set older ages to 15 for simplicity
lofmatl <- melt(lofmat, id.vars = c('Year'), value.name = 'propmat', variable.name = 'age') # gives an error, but I don't understand why. looks fine
lofmatl[, age := as.numeric(age)]

# plot maturity curves by year
ggplot(canmatl, aes(age, propmat, group = Year, color = Year)) +
  geom_line(alpha = 0.5)

ggplot(lofmatl, aes(age, propmat, group = Year, color = Year)) +
  geom_line(alpha = 0.5)


# calc 50% maturity from logistic regression models
canmat50 <- canmatl[, .(age50 = calc50p(age, cbind(propmat*1000, (1-propmat)*1000)), pop = 'can'), by = Year]
lofmat50 <- lofmatl[, .(age50 = calc50p(age, cbind(propmat*1000, (1-propmat)*1000)), pop = 'lof'), by = Year]

# combine
mat50 <- rbind(canmat50, lofmat50)

# plot age of 50% maturity by year
ggplot(mat50, aes(Year, age50, group = pop, color = pop)) +
  geom_line() + 
  geom_smooth()

# write out
write.csv(mat50, 'output/age_50percmature.csv', row.names = FALSE)
