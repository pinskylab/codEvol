# Calculate CIs from ANGSD bootstrapped fsts
# read in results from angsd_fst_boot.sh

require(data.table)

# read in the bootstrapped Fsts
can <- fread('analysis/Can_40.Can_14.fsts.boot.txt')
setnames(can, 1:2, c('unweighted', 'weighted'))
can[, pop := 'can']

lof0711 <- fread('analysis/Lof_07.Lof_11.fsts.boot.txt')
setnames(lof0711, 1:2, c('unweighted', 'weighted'))
lof0711[, pop := 'lof0711']

lof0714 <- fread('analysis/Lof_07.Lof_14.fsts.boot.txt')
setnames(lof0714, 1:2, c('unweighted', 'weighted'))
lof0714[, pop := 'lof0714']

lof1114 <- fread('analysis/Lof_11.Lof_14.fsts.boot.txt')
setnames(lof1114, 1:2, c('unweighted', 'weighted'))
lof1114[, pop := 'lof1114']

dat <- rbind(can, lof0711, lof0714, lof1114)

# range and 95% CIs
dat[, .(min = min(weighted), max = max(weighted)), by = pop]
dat[, .(l95 = quantile(weighted, 0.025), u95 = quantile(weighted, 0.975)), by = pop]
