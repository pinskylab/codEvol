# Make SLiM soft sweep simulations
# assumes that burn-in files have been created with slim_runburnins.sh

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
nes <- as.numeric(gsub('slim_burnin_n|_[1234567890]+\\.vcf', '', bifiles)) # find the ne values from the file names
iter <- as.numeric(gsub('slim_burnin_n[[:digit:]]*_|\\.vcf', '', bifiles)) # find the iteration #s from the file names

# Lofoten 07-11 parameters
paramtable <- expand.grid(iter = sort(unique(iter)), ne = sort(unique(nes)), f = c(0.05, 0.2), ftol = 0.025, n1 = 21, n2 = 23, 
                          s = c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5), r = 3.11e-8, g = 11, L = 3e7)
dim(paramtable)

# trim out parameter combinations already run
outfiles <- list.files('analysis/slim_sim/', pattern = 'slim_sim.+_1\\.vcf$') # get the output files (only the first one)
length(outfiles)
outnes <- as.numeric(gsub('slim_sim_n|_s[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the ne values from the file names
outss <- as.numeric(gsub('slim_sim_n[130]+_s|_f[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the s
outfs <- as.numeric(gsub('slim_sim_n[130]+_s[.01357]+_f|_i[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the f
outiters <- as.numeric(gsub('slim_sim_n[130]+_s[.01357]+_f[.025]+_i|_[1]+\\.vcf', '', outfiles)) # find the iteration #s
outtable <- data.frame(iter = outiters, ne = outnes, f = outfs, s = outss, success = 1)

paramtable <- merge(paramtable, outtable, all.x = TRUE) # keep all the initial parameter combinations
paramtable$success[is.na(paramtable$success)] <- 0 # if missing, they didn't run

with(paramtable, aggregate(list(nfail = success == 0), by = list(s = s, f = f, ne = ne), FUN = sum)) # which ones missing? mostly in ne=100, f=0.05


# run slim for Lofoten
sum(paramtable$success == 0) # number to run
for(i in which(paramtable$success == 0)){
  infile = with(paramtable[i,], paste0('analysis/slim_burnin_n', ne, '_', iter, '.vcf'))
  outfile = with(paramtable[i,], paste0('analysis/slim_sim/slim_sim_n', ne, '_s', s, '_f', f, '_i', iter))
  with(paramtable[i,], runslim(L = L, ne = ne, f = f, ftol = ftol, n1 = n1, n2 = n2, s = s, r = r, g = g, o = outfile, i = infile))
  
}

# check which sims ran successfully
outfiles <- list.files('analysis/slim_sim/', pattern = 'slim_sim.+_1\\.vcf$') # get the output files (only the first one)
length(outfiles)
outnes <- as.numeric(gsub('slim_sim_n|_s[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the ne values from the file names
outss <- as.numeric(gsub('slim_sim_n[130]+_s|_f[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the s
outfs <- as.numeric(gsub('slim_sim_n[130]+_s[.01357]+_f|_i[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the f
outiters <- as.numeric(gsub('slim_sim_n[130]+_s[.01357]+_f[.025]+_i|_[1]+\\.vcf', '', outfiles)) # find the iteration #s
outtable <- data.frame(iter = outiters, ne = outnes, f = outfs, s = outss, success = 1)

outtable <- merge(paramtable[, !grepl('success', names(paramtable))], outtable, all.x = TRUE) # keep all the initial parameter combinations
outtable$success[is.na(outtable$success)] <- 0 # if missing, they didn't run
nrow(outtable)
sum(outtable$success) # 260 of 300 simulations ran, 40 failed (first run). 790 of 900 (110 failed, second run)
outtable[outtable$success == 0, ]
aggregate(list(nfail = outtable$success == 0), by = list(s = outtable$s, f = outtable$f, ne = outtable$ne), FUN = sum) # mostly in ne=100, f=0.05

sum(outtable$success == 0) # number missing
