##############
# Functions
##############

require(data.table)
require(mgcv)
require(Hmisc)

# weighted mean for use with summarize in Hmisc
# x in col 1, w in col 2
wmean <- function(x){
	return(weighted.mean(x=x[,1], w=x[,2]))
}

############
# Data prep
############

# read in
age <- fread('data/kjesbu_fecundity/Rollesfsen1953_Fig1.csv') # mature age distribution (Lofoten) 1932-1953
catchage <- fread('data/kjesbu_fecundity/AFWG_Table3.6catchnumbersatage.csv') # numbers at age from catch data 1946-2015. age13 is "gp+", presumably 13+
stockage <- fread('data/kjesbu_fecundity/AFWG_Table3.16stocknumbersatage.csv') # numbers at age from stock assessment 1946-2015
maturity <- fread('data/kjesbu_fecundity/AFWG_Table3.11proportionmatureatage.csv') # proportion mature at age, for use with catchage and stockage
fec <- fread('data/kjesbu_fecundity/Fekunditet Andenes 1986_2006 Excel_2.csv') # fecundity data for making a model
liv <- fread('data/kjesbu_fecundity/Leverdata.csv') # liver index

# trim fecundity data
fec <- fec[fec$Main_series=='Y',] # trim to main series
fec <- fec[fec$Oto_Type %in% 3:5,] # trim to skrei
	
# calc liver index
fec$liver_index <- 100*fec$Fish_Liver/fec$Fish_Weight

# rescale age distribution percents to sum to 100 (Lofoften data)
yrs <- sort(unique(age$year))
age$percent <- NA
for(i in 1:length(yrs)){
	inds <- age$year == yrs[i]
	thissum <- sum(age$percent_rough[inds])
	age$percent[inds] <- age$percent_rough[inds]/thissum*100
}

# Reshape AFWG data to long format
ctc <- reshape(catchage, varying=names(catchage)[2:12], direction='long', sep='', timevar='age', v.names='catch', times=3:13, idvar='year') # catch data. age 13 is 13+
stk <- reshape(stockage, varying=names(stockage)[2:14], direction='long', sep='', timevar='age', v.names='number', times=1:13, idvar='year')
mat <- reshape(maturity, varying=names(maturity)[2:12], direction='long', sep='', timevar='age', v.names='prop', times=3:13, idvar='year') # catch data. age 13 is 13+

# Calculate proportion mature individuals in each age class (AFWG data)
ctc2 <- merge(ctc, mat, by=c('year', 'age'))
stk2 <- merge(stk, mat, by=c('year', 'age')) # drops ages 1-2

	# number mature at each age
ctc2$nummat <- ctc2$catch * ctc2$prop
stk2$nummat <- stk2$number * stk2$prop

	# scale to a percent
yrs <- sort(unique(ctc2$year))
ctc2$percent <- NA
for(i in 1:length(yrs)){
	inds <- ctc2$year == yrs[i]
	thissum <- sum(ctc2$nummat[inds])
	ctc2$percent[inds] <- ctc2$nummat[inds]/thissum*100
}

yrs <- sort(unique(stk2$year))
stk2$percent <- NA
for(i in 1:length(yrs)){
	inds <- stk2$year == yrs[i]
	thissum <- sum(stk2$nummat[inds])
	stk2$percent[inds] <- stk2$nummat[inds]/thissum*100
}

# Add liver data to age distribution data
age <- merge(age, liv, by.x='year', by.y='Year')
ctc2 <- merge(ctc2, liv, by.x='year', by.y='Year') # lose 2014 and 2015
stk2 <- merge(stk2, liv, by.x='year', by.y='Year') # lose 2014 and 2015

####################
# Fecundity model
####################

# basic plots of fecundity data
quartz(width=7, height=5)
par(mfrow=c(2,1))
with(fec, plot(Fec_Mill ~ Oto_Age))
with(fec, plot(Fec_Mill ~ liver_index))

# model

mod1 <- lm(Fec_Mill ~ liver_index + Oto_Age, data=fec) # no interactions
	summary(mod1)
	par(mfrow=c(2,2)); plot(mod1)
	plot(mod1$fitted + mod1$residuals, mod1$fitted); abline(0,1)

mod2 <- lm(Fec_Mill ~ liver_index + Oto_Age +I(Oto_Age^2), data=fec)
	summary(mod2)
	par(mfrow=c(2,2)); plot(mod2)
	plot(mod2$fitted + mod2$residuals, mod2$fitted); abline(0,1)
	plot(mod2$model$Oto_Age, mod2$fitted)
	plot(mod2$model$liver_index, mod2$fitted)
	plot(mod2$model$liver_index*mod2$model$Oto_Age, mod2$fitted)


AIC(mod1, mod2)

# GAMs	
gam1 <- gam(Fec_Mill ~ s(liver_index) + s(Oto_Age), data=fec)
	plot(gam1, pages=1, residuals=TRUE)
	plot(gam1$fitted + gam1$residuals, gam1$fitted); abline(0,1)



# test for sensitivity to oocyte diameter

mod2_g500 <- lm(Fec_Mill ~ liver_index + Oto_Age +I(Oto_Age^2), data=fec[fec$Ooc_Mean_Dia > 500,])
mod2_g600 <- lm(Fec_Mill ~ liver_index + Oto_Age +I(Oto_Age^2), data=fec[fec$Ooc_Mean_Dia > 600,])


gam1_g500 <- gam(Fec_Mill ~ s(liver_index) + s(Oto_Age), data=fec[fec$Ooc_Mean_Dia > 500,]) # only egg > 500 um, to test sensitivity
	plot(gam1_g500, pages=1, residuals=TRUE)
	plot(gam1_g500$fitted + gam1_g500$residuals, gam1_g500$fitted, xlab='observed', ylab='predicted'); abline(0,1)

gam1_g600 <- gam(Fec_Mill ~ s(liver_index) + s(Oto_Age), data=fec[fec$Ooc_Mean_Dia > 600,]) # only egg > 600 um, to test sensitivity
	plot(gam1_g600, pages=1, residuals=TRUE)
	plot(gam1_g600$fitted + gam1_g600$residuals, gam1_g600$fitted, xlab='observed', ylab='predicted'); abline(0,1)


######################
# Do calculations
######################

# calc fecundity by age in each year (from each of 6 models)
	# Lofoten data
age$fec_lm <- predict(mod2, newdata=data.frame(liver_index=age$"HSI Commercial", Oto_Age=age$age))
age$fec_lm_g500 <- predict(mod2_g500, newdata=data.frame(liver_index=age$"HSI Commercial", Oto_Age=age$age))
age$fec_lm_g600 <- predict(mod2_g600, newdata=data.frame(liver_index=age$"HSI Commercial", Oto_Age=age$age))

age$fec_gam <- predict(gam1, newdata=data.frame(liver_index=age$"HSI Commercial", Oto_Age=age$age))
age$fec_gam_g500 <- predict(gam1_g500, newdata=data.frame(liver_index=age$"HSI Commercial", Oto_Age=age$age))
age$fec_gam_g600 <- predict(gam1_g600, newdata=data.frame(liver_index=age$"HSI Commercial", Oto_Age=age$age))

	# AFWG catch data
ctc2$fec_lm <- predict(mod2, newdata=data.frame(liver_index=ctc2$"HSI Commercial", Oto_Age=ctc2$age))
ctc2$fec_lm_g500 <- predict(mod2_g500, newdata=data.frame(liver_index=ctc2$"HSI Commercial", Oto_Age=ctc2$age))
ctc2$fec_lm_g600 <- predict(mod2_g600, newdata=data.frame(liver_index=ctc2$"HSI Commercial", Oto_Age=ctc2$age))

ctc2$fec_gam <- predict(gam1, newdata=data.frame(liver_index=ctc2$"HSI Commercial", Oto_Age=ctc2$age))
ctc2$fec_gam_g500 <- predict(gam1_g500, newdata=data.frame(liver_index=ctc2$"HSI Commercial", Oto_Age=ctc2$age))
ctc2$fec_gam_g600 <- predict(gam1_g600, newdata=data.frame(liver_index=ctc2$"HSI Commercial", Oto_Age=ctc2$age))

	# turn negative fecundities to zero
for(i in grep('fec', names(ctc2))){
  thiscol <- names(ctc2)[i]
	ctc2[ctc2[,get(thiscol)]<0 & !is.na(ctc2[,get(thiscol)]), (thiscol) := 0]
}

	# AFWG stock assessment data
stk2$fec_lm <- predict(mod2, newdata=data.frame(liver_index=stk2$"HSI Commercial", Oto_Age=stk2$age))
stk2$fec_lm_g500 <- predict(mod2_g500, newdata=data.frame(liver_index=stk2$"HSI Commercial", Oto_Age=stk2$age))
stk2$fec_lm_g600 <- predict(mod2_g600, newdata=data.frame(liver_index=stk2$"HSI Commercial", Oto_Age=stk2$age))

stk2$fec_gam <- predict(gam1, newdata=data.frame(liver_index=stk2$"HSI Commercial", Oto_Age=stk2$age))
stk2$fec_gam_g500 <- predict(gam1_g500, newdata=data.frame(liver_index=stk2$"HSI Commercial", Oto_Age=stk2$age))
stk2$fec_gam_g600 <- predict(gam1_g600, newdata=data.frame(liver_index=stk2$"HSI Commercial", Oto_Age=stk2$age))

	# turn negative fecundities to zero
for(i in grep('fec', names(stk2))){
  thiscol <- names(stk2)[i]
  stk2[stk2[,get(thiscol)]<0 & !is.na(stk2[,get(thiscol)]), (thiscol) := 0]
}


# examine the 3 models
with(age, plot(fec_lm, fec_lm_g500)) # similar
with(age, plot(fec_lm, fec_lm_g600)) # similar
with(age, plot(fec_lm, fec_gam)) # different for high fecundity
with(age, plot(fec_gam, fec_gam_g500)) # similar
with(age, plot(fec_gam, fec_gam_g600)) # different for old fish

# calc generation time
	# Lofoten data
genlof <- summarize(X=cbind(age$age, age$percent*age$fec_lm), by=list(year=age$year), FUN=wmean, stat.name='gen_lm') # generation time: mean age, weighted by egg production
genlof2 <- summarize(X=cbind(age$age, age$percent*age$fec_lm_g500), by=list(year=age$year), FUN=wmean, stat.name='gen_lm_g500') # lm g500 model
genlof3 <- summarize(X=cbind(age$age, age$percent*age$fec_lm_g600), by=list(year=age$year), FUN=wmean, stat.name='gen_lm_g600') # lm g600 model
genlof4 <- summarize(X=cbind(age$age, age$percent*age$fec_gam), by=list(year=age$year), FUN=wmean, stat.name='gen_gam') # gam all model
genlof5 <- summarize(X=cbind(age$age, age$percent*age$fec_gam_g500), by=list(year=age$year), FUN=wmean, stat.name='gen_gam_g500') # for gam g500 model
genlof6 <- summarize(X=cbind(age$age, age$percent*age$fec_gam_g600), by=list(year=age$year), FUN=wmean, stat.name='gen_gam_g600') # for gam g600 model

genlof <- merge(genlof, genlof2)
genlof <- merge(genlof, genlof3)
genlof <- merge(genlof, genlof4)
genlof <- merge(genlof, genlof5)
genlof <- merge(genlof, genlof6)

	# AFWG catch data
genctc <- summarize(X=cbind(ctc2$age, ctc2$percent*ctc2$fec_lm), by=list(year=ctc2$year), FUN=wmean, stat.name='gen_lm') # generation time: mean age, weighted by egg production
genctc2 <- summarize(X=cbind(ctc2$age, ctc2$percent*ctc2$fec_lm_g500), by=list(year=ctc2$year), FUN=wmean, stat.name='gen_lm_g500') # lm g500 model
genctc3 <- summarize(X=cbind(ctc2$age, ctc2$percent*ctc2$fec_lm_g600), by=list(year=ctc2$year), FUN=wmean, stat.name='gen_lm_g600') # lm g600 model
genctc4 <- summarize(X=cbind(ctc2$age, ctc2$percent*ctc2$fec_gam), by=list(year=ctc2$year), FUN=wmean, stat.name='gen_gam') # gam all model
genctc5 <- summarize(X=cbind(ctc2$age, ctc2$percent*ctc2$fec_gam_g500), by=list(year=ctc2$year), FUN=wmean, stat.name='gen_gam_g500') # for gam g500 model
genctc6 <- summarize(X=cbind(ctc2$age, ctc2$percent*ctc2$fec_gam_g600), by=list(year=ctc2$year), FUN=wmean, stat.name='gen_gam_g600') # for gam g600 model

genctc <- merge(genctc, genctc2)
genctc <- merge(genctc, genctc3)
genctc <- merge(genctc, genctc4)
genctc <- merge(genctc, genctc5)
genctc <- merge(genctc, genctc6)

	# AFWG assessment data
genstk <- summarize(X=cbind(stk2$age, stk2$percent*stk2$fec_lm), by=list(year=stk2$year), FUN=wmean, stat.name='gen_lm') # generation time: mean age, weighted by egg production
genstk2 <- summarize(X=cbind(stk2$age, stk2$percent*stk2$fec_lm_g500), by=list(year=stk2$year), FUN=wmean, stat.name='gen_lm_g500') # lm g500 model
genstk3 <- summarize(X=cbind(stk2$age, stk2$percent*stk2$fec_lm_g600), by=list(year=stk2$year), FUN=wmean, stat.name='gen_lm_g600') # lm g600 model
genstk4 <- summarize(X=cbind(stk2$age, stk2$percent*stk2$fec_gam), by=list(year=stk2$year), FUN=wmean, stat.name='gen_gam') # gam all model
genstk5 <- summarize(X=cbind(stk2$age, stk2$percent*stk2$fec_gam_g500), by=list(year=stk2$year), FUN=wmean, stat.name='gen_gam_g500') # for gam g500 model
genstk6 <- summarize(X=cbind(stk2$age, stk2$percent*stk2$fec_gam_g600), by=list(year=stk2$year), FUN=wmean, stat.name='gen_gam_g600') # for gam g600 model

genstk <- merge(genstk, genstk2)
genstk <- merge(genstk, genstk3)
genstk <- merge(genstk, genstk4)
genstk <- merge(genstk, genstk5)
genstk <- merge(genstk, genstk6)


# Average generation time across all years
colMeans(genlof[,2:ncol(genlof)]) # generation time with each of the models
colMeans(genctc[,2:ncol(genctc)], na.rm=TRUE) # generation time with each of the models
colMeans(genstk[,2:ncol(genstk)], na.rm=TRUE) # generation time with each of the models


# Compare the overlapping years

quartz(width=6, height=4)
# pdf(width=6, height=4, file='generationtime.pdf')
par(mai=c(1,1,0.2, 0.2), las=1)
with(genlof, plot(year, gen_gam, lty=1, type='o', ylim=c(6,13), xlim=c(1932,2013), xlab='Year', ylab='Generation time (years)', pch=16, cex=0.5))
#with(genlof, lines(year, gen_lm, col='black', lty=2))
#with(genlof, lines(year, gen_lm_g500, col='black', lty=3))
#with(genlof, lines(year, gen_lm_g600, col='black', lty=4))
#with(genlof, lines(year, gen_gam_g500, col='black', lty=5))
#with(genlof, lines(year, gen_gam_g600, col='black', lty=6))

with(genctc, lines(year, gen_gam, col='blue', lty=1, pch=16, cex=0.5, type='o'))
#with(genctc, lines(year, gen_lm, col='blue', lty=2))
#with(genctc, lines(year, gen_lm_g500, col='blue', lty=3))
#with(genctc, lines(year, gen_lm_g600, col='blue', lty=4))
#with(genctc, lines(year, gen_gam_g500, col='blue', lty=5))
#with(genctc, lines(year, gen_gam_g600, col='blue', lty=6))

with(genstk, lines(year, gen_gam, col='green', lty=1, pch=16, cex=0.5, type='o'))
#with(genstk, lines(year, gen_lm, col='green', lty=2))
#with(genstk, lines(year, gen_lm_g500, col='green', lty=3))
#with(genstk, lines(year, gen_lm_g600, col='green', lty=4))
#with(genstk, lines(year, gen_gam_g500, col='green', lty=5))
#with(genstk, lines(year, gen_gam_g600, col='green', lty=6))

legend('bottomleft', legend=c('Rollefsen', 'AFWG catch', 'AFWG assessment'), lty=1, col=c('black', 'blue', 'green'), bty='n')

dev.off()


# Calculation for 1907-2014
final <- data.frame(year=c(1907:1945, genstk$year, 2014), gentime=c(rep(11.2, 39), genstk$gen_gam, genstk$gen_gam[genstk$year==2012]))
final$gentime[final$year==2013] <- final$gentime[final$year==2012] # fill in missing

with(final, plot(year, gentime, type='o'))
	# 1: 1907-1918
	# 2: 1918-1929
	# 3: 1929-1940
	# 4: 1940-1951
	# 5: 1951-1961
	# 6: 1961-1970
	# 7: 1970-1979
	# 8: 1979-1987
	# 9: 1987-1995
	# 10: 1995-2003
	# 11: 2003-2011

mean(final$gentime)