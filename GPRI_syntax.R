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


################################################
# Confirmatory Factor analysis
################################################
# Factors:
# AGT = Adjustment to Genetic Testing (12 items)
# HMH = History of Metal Health (5 items)
# PFH = Personal & Family History (3 items)


# Model 1: full model (ORIGINAL MODEL)
## items 17 and 14 are perfectly correlated (r=1)
## some estimated lv variances are negative
model1 <- "
AGT =~ GPRI_ITEM4a + GPRI_ITEM4b + GPRI_ITEM4c + GPRI_ITEM5 + GPRI_ITEM6 + 
GPRI_ITEM7 + GPRI_ITEM8 + GPRI_ITEM9 + GPRI_ITEM10 + GPRI_ITEM11 + 
GPRI_ITEM12 + GPRI_ITEM13

HMH =~ GPRI_ITEM14 + GPRI_ITEM15 + GPRI_ITEM16 + GPRI_ITEM17 + GPRI_ITEM18

PFH =~ GPRI_ITEM1 + GPRI_ITEM2a + GPRI_ITEM3 
"
fit_model1 <- cfa(model1, data=DF_ITEMS, estimator="WLSMV", meanstructure = TRUE, missing="pairwise",
                  ordered=c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16", "GPRI_ITEM17",
                            "GPRI_ITEM18",
                            "GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3"))
summary(fit_model1, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE)
fitMeasures(fit_model1, c("chisq.scaled", "df.scaled", "pvalue.scaled", 
                          "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", 
                          "cfi.scaled", "srmr"))
modificationindices(fit_model1, minimum.value = 10)

# Model 2: Original model without item 17
## Heywood cases involving items for PFH factor 

model2 <- "
AGT =~ GPRI_ITEM4a + GPRI_ITEM4b + GPRI_ITEM4c + GPRI_ITEM5 + GPRI_ITEM6 + 
GPRI_ITEM7 + GPRI_ITEM8 + GPRI_ITEM9 + GPRI_ITEM10 + GPRI_ITEM11 + 
GPRI_ITEM12 + GPRI_ITEM13

HMH =~ GPRI_ITEM14 + GPRI_ITEM15 + GPRI_ITEM16 + GPRI_ITEM18

PFH =~ GPRI_ITEM1 + GPRI_ITEM2a + GPRI_ITEM3 
"
fit_model2 <- cfa(model2, data=DF_ITEMS, estimator="WLSMV", meanstructure = TRUE, missing="pairwise", 
                  ordered=c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16",
                            "GPRI_ITEM18",
                            "GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3"))
summary(fit_model2, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_model2, c("chisq.scaled", "df.scaled", "pvalue.scaled", 
                          "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", 
                          "cfi.scaled", "srmr"))

modificationindices(fit_model2, minimum.value = 100)

# Model 3: Original Model without item 17 plus modification indexes 
## Heywood cases involving items for PFH factor 
model3 <- "
AGT =~ GPRI_ITEM4a + GPRI_ITEM4b + GPRI_ITEM4c + GPRI_ITEM5 + GPRI_ITEM6 + 
GPRI_ITEM7 + GPRI_ITEM8 + GPRI_ITEM9 + GPRI_ITEM10 + GPRI_ITEM11 + 
GPRI_ITEM12 + GPRI_ITEM13

HMH =~ GPRI_ITEM14 + GPRI_ITEM15 + GPRI_ITEM16 + GPRI_ITEM18

PFH =~ GPRI_ITEM1 + GPRI_ITEM2a + GPRI_ITEM3 

AGT =~ GPRI_ITEM14
GPRI_ITEM12 ~~ GPRI_ITEM13
GPRI_ITEM4c ~~  GPRI_ITEM6

"
fit_model3 <- cfa(model3, data=DF_ITEMS, estimator="WLSMV", meanstructure = TRUE, missing="pairwise", 
                  ordered=c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16",
                            "GPRI_ITEM18",
                            "GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3"))
summary(fit_model3, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_model3, c("chisq.scaled", "df.scaled", "pvalue.scaled", 
                          "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", 
                          "cfi.scaled", "srmr"))

anova(fit_model2, fit_model3) # Comparing models

# Model 4: Model propose by Maheu & Esplen, 2018 
## Heywood cases involving items for PFH factor 

model4 <- "
Exper_AGT =~ GPRI_ITEM4a + GPRI_ITEM8 + GPRI_ITEM9 + GPRI_ITEM10 + 
GPRI_ITEM12 + GPRI_ITEM13

Ant_AGT =~ GPRI_ITEM4b + GPRI_ITEM4c + GPRI_ITEM5 + GPRI_ITEM6 + 
GPRI_ITEM7 + GPRI_ITEM11

HMH =~ GPRI_ITEM14 + GPRI_ITEM15 + GPRI_ITEM16 +  GPRI_ITEM18

PFH =~ GPRI_ITEM1 + GPRI_ITEM2a + GPRI_ITEM3 

Ant_AGT =~ GPRI_ITEM14
GPRI_ITEM12 ~~ GPRI_ITEM13

PFH ~~ VAR*PFH
"
fit_model4 <- cfa(model4, data=DF_ITEMS, estimator="WLSMV", meanstructure = TRUE, 
                  missing="pairwise",
                  ordered=c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16", "GPRI_ITEM17",
                            "GPRI_ITEM18",
                            "GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3"))
summary(fit_model4, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_model4, c("chisq.scaled", "df.scaled", "pvalue.scaled", 
                          "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", 
                          "cfi.scaled", "srmr"))
modificationindices(fit_model4, minimum.value = 5)
anova(fit_model3, fit_model4)

# Model 5: Models without PFH factor
## Heywood cases removes
## This model provided adequate fit to the data

model5 <- "
AGT =~ GPRI_ITEM4a + GPRI_ITEM4b + GPRI_ITEM4c + GPRI_ITEM5 + GPRI_ITEM6 + 
GPRI_ITEM7 + GPRI_ITEM8 + GPRI_ITEM9 + GPRI_ITEM10 + GPRI_ITEM11 + 
GPRI_ITEM12 + GPRI_ITEM13

HMH =~ GPRI_ITEM14 + GPRI_ITEM15 + GPRI_ITEM16 + GPRI_ITEM18

AGT =~ GPRI_ITEM14
GPRI_ITEM12 ~~ GPRI_ITEM13 
GPRI_ITEM4c ~~  GPRI_ITEM6
"
fit_model5 <- cfa(model5, data=DF_ITEMS, estimator="WLSMV", meanstructure = TRUE, missing="pairwise",
                  ordered=c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16",
                            "GPRI_ITEM18"))
summary(fit_model5, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_model5, c("chisq.scaled", "df.scaled", "pvalue.scaled", 
                          "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", 
                          "cfi.scaled", "srmr"))

modificationindices(fit_model5, minimum.value = 3)

# Model 6: Final model
## Best fitting model
model6 <- "
Exper_AGT =~ GPRI_ITEM4a + GPRI_ITEM8 + GPRI_ITEM9 + GPRI_ITEM10 + 
GPRI_ITEM12 + GPRI_ITEM13

Ant_AGT =~ GPRI_ITEM4b + GPRI_ITEM4c + GPRI_ITEM5 + GPRI_ITEM6 + 
GPRI_ITEM7 + GPRI_ITEM11

HMH =~ GPRI_ITEM14 + GPRI_ITEM15 + GPRI_ITEM16 +  GPRI_ITEM18

Ant_AGT =~ GPRI_ITEM14
GPRI_ITEM12 ~~ GPRI_ITEM13
"
fit_model6 <- cfa(model6, data=DF_ITEMS, estimator="WLSMV", meanstructure = TRUE, 
                  missing="pairwise",
                  ordered=c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16", "GPRI_ITEM17",
                            "GPRI_ITEM18",
                            "GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3"))
summary(fit_model6, rsquare = TRUE, fit.measures = TRUE, standardized = TRUE)

fitMeasures(fit_model6, c("chisq.scaled", "df.scaled", "pvalue.scaled", 
                          "rmsea.scaled", "rmsea.ci.lower.scaled", "rmsea.ci.upper.scaled", 
                          "cfi.scaled", "srmr"))
anova(fit_model5, fit_model6)

model6_par <- parameterestimates(fit_model6, zstat=FALSE, pvalue=TRUE)

model6_par[model6_par$op!="|" & 
             model6_par$op!="~*~" &
             model6_par$op!="~1",]
# do not run
# write.csv(
#   model6_par[model6_par$op != "|" &
#                model6_par$op != "~*~" &
#                model6_par$op != "~1", ],
#   "G:\\My Drive\\FPCEUP\\I&D Projects\\Together\\Together_Rproject\\model6parameters.csv",
#   row.names = T
# )


################################################
# Reliability analysis
################################################
PIGT <- dplyr::select(DF_ITEMS, c("GPRI_ITEM4a", "GPRI_ITEM4b", "GPRI_ITEM4c", "GPRI_ITEM5", 
                                  "GPRI_ITEM6", "GPRI_ITEM7", "GPRI_ITEM8", "GPRI_ITEM9", 
                                  "GPRI_ITEM10", "GPRI_ITEM11", "GPRI_ITEM12", "GPRI_ITEM13"))

IIGT <- dplyr::select(DF_ITEMS, c("GPRI_ITEM4a", "GPRI_ITEM8", "GPRI_ITEM9", 
                                  "GPRI_ITEM10", "GPRI_ITEM12", "GPRI_ITEM13"))

EIGT <- dplyr::select(DF_ITEMS, c("GPRI_ITEM4b", "GPRI_ITEM4c", "GPRI_ITEM5", 
                                  "GPRI_ITEM6", "GPRI_ITEM7", "GPRI_ITEM11"))

HMHC <- dplyr::select(DF_ITEMS, c("GPRI_ITEM14", "GPRI_ITEM15", "GPRI_ITEM16", 
                                  "GPRI_ITEM18"))

PFH <- dplyr::select(DF_ITEMS, c("GPRI_ITEM1", "GPRI_ITEM2a", "GPRI_ITEM3"))

psych::alpha(PIGT, check.keys=T)
psych::alpha(IIGT, check.keys=T)
psych::alpha(EIGT, check.keys=T)
psych::alpha(HMHC, check.keys=T)
psych::alpha(PFH, check.keys=T)

# Calculate the Kuder-Richardson formula 20 for estimate of reliability for tests with dichotomous items
kr20 <- function(x,hit=1)
{
  x <- na.omit(x)
  k <- ncol(x)
  n <- nrow(x)
  x <- data.frame(apply(x==hit, 2, as.numeric))
  totalVar <- var(apply(x, 1, sum))
  p <- colSums(x)/n
  q <- 1-p
  r <- (k/(k-1))*(1-sum(p*q)/totalVar)
  return(r)
}

# function from: https://rdrr.io/github/DavideMassidda/testing/man/kr20.html
kr20(HMHC)



################################################
# Convergent validity
################################################
# AGT (Adjustment to Genetic Testing)
# HMHC (History of Metal Health Concerns)
# PFHL (Personal/Family History of Loss)

# Correlations

cor_names <- c("agt_exper_new", "agt_ant_new", "hmhc_new", 
               "ies_intrusion", "ies_avoidance", "ies_total",
               "hads_anxiety", "hads_depression", 
               "core_distress", 
               "age", "education", "income", "nr_children")

rcorr(as.matrix(DF_CONVAL[ ,cor_names]))
corr.test(DF_CONVAL[ ,cor_names])

# describe by sex
describeBy(DF_CONVAL[c(4, 21:23)], "sex" )

# t test
t.test(agt_exper_new ~ sex, DF_CONVAL) # 0 = female, 1 = male
cohensD(agt_exper_new ~ sex,
        data = DF_CONVAL)

t.test(agt_ant_new ~ sex, DF_CONVAL) # 0 = female, 1 = male
cohensD(agt_ant_new ~ sex,
        data = DF_CONVAL)

t.test(gpri_hmhc ~ sex, DF_CONVAL) # 0 = female, 1 = male
cohensD(gpri_hmhc ~ sex,
        data = DF_CONVAL)

# describe by cancer types
describeBy(DF_CONVAL[c(10,  21:23)], "hist_cancer" )

# t-tests
t.test(agt_exper_new ~ hist_cancer, DF_CONVAL)
cohensD(agt_exper_new ~ hist_cancer,
        data = DF_CONVAL)

t.test(agt_ant_new  ~ hist_cancer, DF_CONVAL)
cohensD(agt_ant_new  ~ hist_cancer,
        data = DF_CONVAL)

t.test(gpri_hmhc ~ hist_cancer, DF_CONVAL)
cohensD(gpri_hmhc ~ hist_cancer,
        data = DF_CONVAL)