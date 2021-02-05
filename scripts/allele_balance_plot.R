# plot output from allele_balance_calc.R
require(data.table)

# read in
can13 <- fread('analysis/allele_balance_popCan13.csv.gz', header= FALSE)
can40 <- fread('analysis/allele_balance_popCan40.csv.gz', header= FALSE)
lof07 <- fread('analysis/allele_balance_popLof07.csv.gz', header= FALSE)
lof11 <- fread('analysis/allele_balance_popLof11.csv.gz', header= FALSE)
lof14 <- fread('analysis/allele_balance_popLof14.csv.gz', header= FALSE)



# summary
summary(can13$V1)
summary(can40$V1)
summary(lof07$V1)
summary(lof11$V1)
summary(lof14$V1)

# t-test
t.test(can13$V1)

# plot
hist(can13$V1)
hist(can40$V1)
hist(lof07$V1)
hist(lof11$V1)
hist(lof14$V1)

