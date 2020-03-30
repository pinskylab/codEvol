# Make SLiM soft sweep simulations
# then plot Fst from sliding-window

require(data.table)


###########
# functions
###########

runslim <- function(L=30000000, ne=50, f=0.1, ftol = 0.025, n1=22, n2=22, s=1.0, r=3.11e-8, g=11, o='tmp/test', i='tmp/test.vcf'){
  # notice complicated quoting needed to pass a string to slim
  # use -long to enable verbose output for troubleshooting
  system2("slim", c(paste0("-d ", c(paste0("L=", L), paste0("ne=", ne), paste0("f=", f), paste0("ftol=", ftol), 
                                    paste0("n1=", n1), paste0("n2=", n2), paste0("s=", s), paste0("r=", r), 
                                    paste0("g=", g), paste0("\"o='", o, "'\""), paste0("\"i='", i, "'\""))), 
                    '-long', " scripts/slim_softsweep.slim"))
}

###########
# run slim
###########

bifiles <- list.files('analysis/', pattern = 'slim_burnin_n') # get the burn-in files
nes <- as.numeric(gsub('slim_burnin_n|_[12345]\\.vcf', '', bifiles)) # find the ne values from the file names
iter <- as.numeric(gsub('slim_burnin_n[[:digit:]]*_|\\.vcf', '', bifiles)) # find the iteration #s from the file names

# Lofoten 07-11 parameters
paramtable <- expand.grid(iter = sort(unique(iter)), ne = sort(unique(nes)), f = c(0.05, 0.2), ftol = 0.025, n1 = 21, n2 = 23, 
                          s = c(0.1, 0.3, 1, 1.5), r = 3.11e-8, g = 11, L = 3e7)
dim(paramtable)

# run slim for Lofoten
for(i in 1:nrow(paramtable)){
  infile = with(paramtable[i,], paste0('analysis/slim_burnin_n', ne, '_', iter, '.vcf'))
  outfile = with(paramtable[i,], paste0('analysis/slim_sim/slim_sim_n', ne, '_s', s, '_f', f, '_i', iter))
  with(paramtable[i,], runslim(L = L, ne = ne, f = f, ftol = ftol, n1 = n1, n2 = n2, s = s, r = r, g = g, o = outfile, i = infile))
  
}

# check which sims ran successfully
outfiles <- list.files('analysis/slim_sim/', pattern = 'slim_sim.+_1\\.vcf$') # get the output files (only the first one)
length(outfiles)
outnes <- as.numeric(gsub('slim_sim_n|_s[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the ne values from the file names
outss <- as.numeric(gsub('slim_sim_n[130]+_s|_f[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the s
outfs <- as.numeric(gsub('slim_sim_n[130]+_s[.0135]+_f|_i[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the f
outiters <- as.numeric(gsub('slim_sim_n[130]+_s[.0135]+_f[.025]+_i|_[1]+\\.vcf', '', outfiles)) # find the iteration #s
outtable <- data.frame(iter = outiters, ne = outnes, f = outfs, s = outss, success = 1)

outtable <- merge(paramtable, outtable, all.x = TRUE) # keep all the initial parameter combinations
outtable$success[is.na(outtable$success)] <- 0 # if missing, they didn't run
nrow(outtable)
sum(outtable$success) # 172 of 200 simulations ran, 28 failed
outtable[outtable$success == 0, ]
aggregate(list(nfail = outtable$success == 0), by = list(s = outtable$s, f = outtable$f, ne = outtable$ne), FUN = sum) # mostly in ne=100, f=0.05

#############
# calculate fst
#############


#rsync -tz analysis/slim_sim/slim_sim_n*.vcf USER@saga.sigma2.no:~/historic_malin/analysis/slim_sim # copy to saga
# on saga:
#sbatch scripts/slim_calcfst.sh # calc fst with vcftools
#rsync -tz USER@saga.sigma2.no:~/historic_malin/analysis/slim_sim/slim_sim_n*.fst analysis/slim_sim/ # copy fsts back from saga

