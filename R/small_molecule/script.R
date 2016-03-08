load("R/small_molecule/bestGLM.Rdata")
args <- commandArgs(trailingOnly = FALSE)
#print(args)
table <- read.table(args[length(args)], sep="\t", header=TRUE)

result <- predict(bestGLM, table)

cat(result)

#print(table)
