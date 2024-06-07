#
# Code for the meta-analysis of the lambda estimates.
# Walasek, L., Mullett, L. T., Stewart, N. (2024). A meta-analysis of loss aversion in risky contexts. Journal of Economic Psychology.
# Contact: Lukasz Walasek, L.Walasek@warwick.ac.uk 
#
# Last update: 6/7/2024
#

library(data.table)
library(metafor)
library(boot)
library(ggplot2)

set.seed(1234)

fits <- read.csv("AllFitsCombined_osf.csv",as.is=TRUE,header=TRUE)

######################################################################
################### Bootstrap parameter estimates ####################
######################################################################

### boot strap medians of lambda

fits$study_name <- factor(fits$study_name)

make.boot.objects <- function(data) {
  boot.median <- function(x, i) { median(x[i]) }
  boot(data=data$lambda.best, statistic=boot.median, R=9999)
}
boot.objs <- by(fits, INDICES=fits$study_name, FUN=make.boot.objects)
boot.objs

### boot strap CIs

get.booted.lambdas <- function(boot.obj) {
  boot.CIs <- boot.ci(boot.obj, type="bca")
  c(median=boot.obj$t0, standard.error=sd(boot.obj$t), LCL.bca=boot.CIs$bca[4], UCL.bca=boot.CIs$bca[5])
}

lambdas <- data.table((t(sapply(boot.objs, get.booted.lambdas))), keep.rownames=TRUE)
setnames(lambdas, "rn", "studyname")

lambdas <- lambdas[,LCL.approx:=median-qnorm(0.975)*standard.error]
lambdas <- lambdas[,UCL.approx:=median+qnorm(0.975)*standard.error]
lambdas

######################################################################
###### Bootstrap medians of alpha ####################################
######################################################################

make.boot.objects <- function(data) {
  boot.median <- function(x, i) { median(x[i]) }
  boot(data=data$alpha.best, statistic=boot.median, R=9999)
}
boot.objs <- by(fits, INDICES=fits$study_name, FUN=make.boot.objects)
boot.objs

### boot strap CIs (alpha)

alphas <- data.table((t(sapply(boot.objs, get.booted.lambdas))), keep.rownames=TRUE)
setnames(alphas, "rn", "studyname")

alphas <- alphas[,LCL.approx:=median-qnorm(0.975)*standard.error]
alphas <- alphas[,UCL.approx:=median+qnorm(0.975)*standard.error]
alphas

#######################################################################
###### Bootstrap medians of delta #####################################
#######################################################################

fits_delta <- fits[!is.na(fits$delta.best),]

make.boot.objects <- function(data) {
  boot.median <- function(x, i) { median(x[i]) }
  boot(data=data$delta.best, statistic=boot.median, R=9999)
}

# Relevel
fits_delta$study_name <- factor(fits_delta$study_name)

boot.objs <- by(fits_delta, INDICES=fits_delta$study_name, FUN=make.boot.objects)
boot.objs

### boot strap CIs (alpha)

deltas <- data.table((t(sapply(boot.objs, get.booted.lambdas))), keep.rownames=TRUE)
setnames(deltas, "rn", "studyname")

deltas <- deltas[,LCL.approx:=median-qnorm(0.975)*standard.error]
deltas <- deltas[,UCL.approx:=median+qnorm(0.975)*standard.error]
deltas

# Create Table 2.

lambdas_t <- round(lambdas[,2:ncol(lambdas)],2)
lambdas_t$studyname <- lambdas$studyname  
lambdas_t$median <- paste0(lambdas_t$median," [",lambdas_t$LCL.bca,",",lambdas_t$UCL.bca,"]")
setnames(lambdas_t,c("median","standard.error"),new = c("lambda","standard.error_lambda"))
lambdas_t <- lambdas_t[,.(studyname,lambda,standard.error_lambda)]

alphas_t <- round(alphas[,2:ncol(alphas)],2)
alphas_t$studyname <- alphas$studyname
alphas_t$median <- paste0(alphas_t$median," [",alphas_t$LCL.bca,",",alphas_t$UCL.bca,"]")
setnames(alphas_t,c("median","standard.error"),new = c("alpha","standard.error_alpha"))
alphas_t <- alphas_t[,.(studyname,alpha,standard.error_alpha)]

deltas_t <- round(deltas[,2:ncol(deltas)],2)
deltas_t$studyname <- deltas$studyname
deltas_t$median <- paste0(deltas_t$median," [",deltas_t$LCL.bca,",",deltas_t$UCL.bca,"]")
setnames(deltas_t,c("median","standard.error"),new = c("gamma","standard.error_gamma"))
deltas_t <- deltas_t[,.(studyname,gamma,standard.error_gamma)]

table_1 <- merge(lambdas_t,alphas_t,by="studyname")
table_1 <- merge(table_1,deltas_t,by="studyname",all.x = TRUE)

#write.csv(table_1,"Table2.csv")

######################################################################

# Perform bootstrapping on log transformed medians.

fits$lambda.best.log <- log(fits$lambda.best)
fits_log <- fits[!is.na(fits$lambda.best.log),]

make.boot.objects <- function(data) {
  boot.median <- function(x, i) { median(x[i]) }
  boot(data=data$lambda.best.log, statistic=boot.median, R=9999)
}

boot.objs_log <- by(fits_log, INDICES=fits_log$study_name, FUN=make.boot.objects)
boot.objs_log

lambdas_log <- data.table((t(sapply(boot.objs_log, get.booted.lambdas))), keep.rownames=TRUE)
setnames(lambdas_log, "rn", "studyname")

lambdas_log <- lambdas_log[,LCL.approx:=median-qnorm(0.975)*standard.error]
lambdas_log <- lambdas_log[,UCL.approx:=median+qnorm(0.975)*standard.error]
lambdas_log

######################################################################
################### Random effect meta-analysis # ####################
######################################################################

# Random effect meta-analysis
(  rma.1 <- rma(yi=median, sei=standard.error, data=lambdas)  )

# Random-Effects Model (k = 19; tau^2 estimator: REML)
# 
# tau^2 (estimated amount of total heterogeneity): 0.1537 (SE = 0.0681)
# tau (square root of estimated tau^2 value):      0.3921
# I^2 (total heterogeneity / total variability):   91.60%
# H^2 (total variability / sampling variability):  11.90
# 
# Test for Heterogeneity:
#   Q(df = 18) = 226.8211, p-val < .0001
# 
# Model Results:
#   
#   estimate      se     zval    pval   ci.lb   ci.ub      
# 1.3128  0.1093  12.0076  <.0001  1.0985  1.5271  *** 
#   
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

with(lambdas,cor.test(median,standard.error,method="spearman"))

fitstats(rma.1)
confint(rma.1)
predict(rma.1)

# Visualization

authors <- c("BrooksP","CanessaC","ChibD","FrydmanC","GlocknerP","KocherP","LorainsD",
             "PahlkeS","PighinB","Rieskampp","Sokol-HessnerH_A","Sokol-HessnerH_R","Sokol-HessnerC_A",
             "Sokol-HessnerC_R","Sokol-HessnerH","TomF","ErnerK","BrooksZ","ZeisbergerV")
dates <- c("(2014)","(2013)","(2012)","(2011)","(2012)","(2013)","(2014)","(2012)","(2014)","(2008)",
           "(2009)","(2009)","(2012)","(2012)","(2014)","(2007)","(2013)","(2005)","(2012)")
study.name <- c()
names <- data.frame(authors=authors,dates=dates,study_name=levels(fits$study_name))

fits <- merge(fits,names,by=c("study_name"))
fits$full_name <- paste(fits$authors,fits$dates,sep=" ")

setnames(lambdas,old = "studyname","study_name")
lambdas <- merge(lambdas,names,by=c("study_name"))
lambdas$full_name <- paste(lambdas$authors,lambdas$dates,sep = " ")

### Reorder for forest plot
lambdas <- lambdas[order(lambdas$median,decreasing = TRUE),]

(  rma.1 <- rma(yi=median, sei=standard.error, data=lambdas)  )
rma.1$slab <- lambdas$full_name

#png(file='forest.png',height = 22,width = 22,res = 300,units = "cm")
forest(rma.1,xlim=c(-0.5,3.9),slab=lambdas$full_name,alim=c(.5,3),at=c(.5,1,1.5,2,2.25,2.5,3),
       xlab = "Median Estimated Lambda",cex=.75,mlab="RE Model for all studies",col="dark gray",
       refline=c(1,2.25))
text(-0.35,21, c("Study"), font=2, cex=1)
text(1.75,21, c("Weights"), font=2, cex=1)
text(3.5,21, c("Medians [95% CIs]"), font=2, cex=1)
#dev.off() 

# funnel plot with inverse standard error
#png(file='funnel.png') # Open PNG device with specific file name
funnel(rma.1,yaxis = "seinv",xlim = c(-1.69,4.31),xlab = "Lambda",ylim = c(.001,25))
#dev.off() 


# funnel plot with inverse standard error (trim and fill)
taf <- trimfill(rma.1)
#png(file='funnel_tf.png') # Open PNG device with specific file name
funnel(taf,yaxis = "seinv",xlim = c(-1.69,4.31),xlab = "Lambda",ylim = c(.001,25))
#dev.off() 

######################################################################
# Appendix plots for log lambda
#

setnames(lambdas_log,old = "studyname","study_name")
lambdas_log <- merge(lambdas_log,names,by=c("study_name"))
lambdas_log$full_name <- paste(lambdas_log$authors,lambdas_log$dates,sep = " ")
lambdas_log <- lambdas_log[order(lambdas_log$median,decreasing = TRUE),]

(  rma.2 <- rma(yi=median, sei=standard.error, data=lambdas_log)  )
rma.2$slab <- lambdas$full_name

#png(file='forrest_log.png',height = 22,width = 22,res = 300,units = "cm")
forest(rma.2,xlim=c(-2,3),slab=lambdas_log$full_name,alim=c(-1,2),at=c(-1,-0.5,0,0.693,1,1.5,2),
       xlab = "Median Estimated Lambda",cex=.75,mlab="RE Model for all studies",col="dark gray",
       refline=c(0,0.693))
text(-1.75,21, c("Study"), font=2, cex=1)
text(0.5,21, c("Weights"), font=2, cex=1)
text(2.35,21, c("Medians [95% CIs]"), font=2, cex=1)
#dev.off() 

# funnel plot with inverse standard error
#png(file='funnel_log.png') # Open PNG device with specific file name
funnel(rma.2,yaxis = "seinv",xlim = c(-1.69,4.31),xlab = "Lambda",ylim = c(.001,25))
#dev.off() 

# funnel plot with inverse standard error (trim and fill)
taf <- trimfill(rma.2)
#png(file='funnel_tf_log.png') # Open PNG device with specific file name
funnel(taf,yaxis = "seinv",xlim = c(-2.69,4.31),xlab = "Lambda",ylim = c(.001,25))
#dev.off()

