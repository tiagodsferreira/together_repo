###########################################
# loading required packages
###########################################
library(pacman) 
p_load(foreign, plyr, dplyr, psych, haven, lsr, nFactors, GPArotation, psychTools, semPlot, VIM, naniar, MissMech, lavaan, semTools, FactoMineR, missMDA)

###########################################
# User defined functions
###########################################
as.numeric.factor <- function(x) {as.numeric(levels(x))[x]} #convert factors in numeric

f <- function(x) {
  freq <- table(x)
  perc <- prop.table(freq)
  return(as.matrix(rbind(freq, perc)))
} # return frequencies and percentages


###########################################
# Load data
###########################################
#do not run
# load("DF_VAL.RData")
# load("DF_ITEMS.RData")
# load("DF_CONVAL.RData")
###########################################
## Missing data analysis
###########################################
# Patterns of missing data
miss_var_summary(DF_ITEMS)
res<-summary(aggr(DF_ITEMS[5:24], sortVar=FALSE))$combinations
head(res[rev(order(res[,2])),])
matrixplot(DF_ITEMS, sortby = 2)
vis_miss(DF_ITEMS, sort_miss = TRUE) 

# Percentage of missing data
sum(is.na(DF_ITEMS))/prod(dim(DF_ITEMS))

# MCAR test
TestMCARNormality(DF_ITEMS[5:24])
