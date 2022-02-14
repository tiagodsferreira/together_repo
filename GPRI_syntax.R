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

###########################################
# Item descriptive statistics
###########################################
sapply(DF_ITEMS, f)

################################################
# Inter-item correlation
################################################
summary(DF_ITEMS)

# Perceived impact and personal adjustment to genetic testing (AGT: Adjustment to Genetic Testing [12 items]): 
names_AGT <- c("GPRI_ITEM4a", "GPRI_ITEM4b", "GPRI_ITEM4c", "GPRI_ITEM5",
               "GPRI_ITEM6", "GPRI_ITEM7", "GPRI_ITEM8", "GPRI_ITEM9",
               "GPRI_ITEM10", "GPRI_ITEM11", "GPRI_ITEM12", "GPRI_ITEM13")
describe(DF_ITEMS[ ,names_AGT])
corr.AGT <- corr.test(DF_ITEMS[ ,names_AGT])

# History of mental health concerns (HMH: History of Metal Health [5 items]): 
names_HMH <- c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16", "GPRI_ITEM17", "GPRI_ITEM18")
sapply(DF_ITEMS[ ,names_HMH], table)
(tetra.HMH <- tetrachoric(DF_ITEMS[ ,names_HMH], correct=0))
x <- tetra.HMH$rho
mean(x[lower.tri(x, diag = FALSE)])

# Personal history/family history/loss to cancer (PFH: Personal & Family History [3 items]): 
names_PFH <- c("GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3")
sapply(DF_ITEMS[ ,names_PFH], table)
(tetra.PFH <- tetrachoric(DF_ITEMS[ ,names_PFH], correct=0))
x <- tetra.PFH$rho
mean(x[lower.tri(x, diag = FALSE)])

#Full correlation using mixed correlation tests
names(DF_ITEMS[ ,c(names_PFH, names_AGT, names_HMH)])
mixed <- mixedCor(DF_ITEMS[ ,c(names_PFH, names_AGT, names_HMH)], c=4:15, d=c(16:20, 1:3), correct = 0)
mixed
mixed$rho

# do not run
# write.csv(mixed$rho, "iter-item correlations.csv", row.names = T)

pvalue <- corr.p(mixed$rho, n=207)
pvalue$p

# do not run
# write.csv(pvalue$p, "iter-item pvalues.csv", row.names = T)

