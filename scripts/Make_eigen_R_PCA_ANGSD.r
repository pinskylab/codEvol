
#Read in file

#accept arguments
args<-commandArgs(TRUE)

#install package
require(RColorBrewer)
require(maps)
require(spam)
require(fields)


C <- as.matrix(read.table(args[1]))




e <- eigen(C)

e.var.per <- round(e$values/sum(e$values)*100, 2)

e_df <- as.data.frame(e$vectors, row.names = NULL)

write.table(e_df, file = paste0(c(args[1]),c(".tab")), sep = "\t", row.names = TRUE, col.names = NA)



#quit
q()

