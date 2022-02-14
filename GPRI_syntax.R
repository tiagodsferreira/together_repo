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

