#Read in file

#accept arguments
args<-commandArgs(TRUE)


pdf(paste0(c(args[1]),c(".pdf")), width = 8, height = 2.5)

tbl=read.table(args[1])


barplot(t(as.matrix(tbl)), col=rainbow(5), xlab="Individual", ylab="Ancestry", border=NA)




#save file
dev.off()
#quit
q()
