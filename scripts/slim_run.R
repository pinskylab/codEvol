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

bifiles <- list.files('analysis/', pattern = 'slim_burnin_n') # get the burn-in files. this will have errors if a burnin file is present both as .vcf and .vcf.zip
nes <- as.numeric(gsub('slim_burnin_n|_[1234567890]+\\.vcf|\\.zip', '', bifiles)) # find the ne values from the file names
iter <- as.numeric(gsub('slim_burnin_n[[:digit:]]*_|\\.vcf|\\.zip', '', bifiles)) # find the iteration #s from the file names

# Lofoten 07-11 parameters
paramtable <- expand.grid(iter = sort(unique(iter)), ne = sort(unique(nes)), f = c(0.05, 0.2), ftol = 0.025, n1 = 22, n2 = 24, 
                          s = c(0, 0.1, 0.3, 0.5, 0.7, 1, 1.5), r = 3.11e-8, g = 11, L = 3e7)
dim(paramtable)

# trim out parameter combinations already run
outfiles <- list.files('analysis/slim_sim/', pattern = 'slim_sim.+_1\\.vcf$') # get the output files (only the first one)
length(outfiles)
if(length(outfiles)>0){
  outnes <- as.numeric(gsub('slim_sim_n|_s[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the ne values from the file names
  outss <- as.numeric(gsub('slim_sim_n[130]+_s|_f[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the s
  outfs <- as.numeric(gsub('slim_sim_n[130]+_s[.01357]+_f|_i[[:alnum:][:punct:]]+vcf', '', outfiles)) # find the f
  outiters <- as.numeric(gsub('slim_sim_n[130]+_s[.01357]+_f[.025]+_i|_[1]+\\.vcf', '', outfiles)) # find the iteration #s
  outtable <- data.frame(iter = outiters, ne = outnes, f = outfs, s = outss, success = 1)
  
  paramtable <- merge(paramtable, outtable, all.x = TRUE) # keep all the initial parameter combinations
  paramtable$success[is.na(paramtable$success)] <- 0 # if missing, they didn't run
} else {
  paramtable$success <- 0
}
  

missing <- with(paramtable, aggregate(list(nfail = success == 0), by = list(s = s, f = f, ne = ne), FUN = sum)) # which ones missing?
missing[missing$nfail > 0,]

# run slim for Lofoten
sum(paramtable$success == 0) # number to run
if(sum(paramtable$success == 0) > 0){
  todo <- which(paramtable$success == 0)
  for(i in 1:length(todo)){
    thisiter <- paramtable$iter[todo[i]]
    infileroot <- with(paramtable[todo[i],], paste0('slim_burnin_n', ne, '_', iter, '.vcf'))
    infileindex <- grep(infileroot, bifiles)
    infile <- paste0('analysis/', bifiles[infileindex])
    if(grepl('zip', infile)){ # unzip the infile if it is zipped and an unzipped version doesn't exist
      infilezip <- infile
      infile <- gsub('\\.zip', '', infile)
      if(!file.exists(infile)){
        unzip(infilezip, unzip = Sys.which('unzip'), exdir = 'analysis/') # uses the system's unzip function to avoid file size limits
      }
    }
    outfile <- with(paramtable[todo[i],], paste0('analysis/slim_sim/slim_sim_n', ne, '_s', s, '_f', f, '_i', iter))
    with(paramtable[todo[i],], runslim(L = L, ne = ne, f = f, ftol = ftol, n1 = n1, n2 = n2, s = s, r = r, g = g, o = outfile, i = infile))
    if(i != length(todo)){
      if(thisiter != paramtable$iter[todo[i+1]]){ # if next iteration is a different burn-in file
        if(file.exists(infilezip)){ # clean up by deleting the unzipped file
          file.remove(infile)
        }
      }
    }
  }
} else {
  print('All done!')
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
missing <- aggregate(list(nfail = outtable$success == 0), by = list(s = outtable$s, f = outtable$f, ne = outtable$ne), FUN = sum)
missing[missing$nfail > 0, ] # mostly in ne=100, f=0.05


sum(outtable$success == 0) # number missing
