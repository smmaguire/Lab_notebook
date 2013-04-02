#########################################################################
#off curve Tdata
hormoneData[1,25]<-sqrt(63.286919)
hormoneData[42,25]<-sqrt(90.9433139)
load("cort.RDATA")
#########################################################################
setwd("/Users/seanmaguire/desktop/sean lab notebook/TC_correlation")
#load("alldata_updated.RDATA")
load("tcCor.RDATA")

hormoneData<-alldata[alldata$day==10,]
library(effects); library(lme4); library(xtable)
########################## Functions ######################################

get.pvalues.mcmc<-function(mcmc){
hpd = lme4::HPDinterval(mcmc)
mcmcfixef = t(mcmc@fixef)
nr <- nrow(mcmcfixef)
prop <- colSums(mcmcfixef > 0)/nr
pvalues <- 2 * pmax(0.5/nr, pmin(prop, 1 - prop))
modelInfo<-list(p.vals=pvalues,confint=hpd)
return(modelInfo)
}

makeNiceTable<-function(mcmcSummary,modSummary){
	coefs<-as.data.frame(modSummary@coefs)
	coefs$p_value<-mcmcSummary$p.vals
	coefs$HPDlower<-as.data.frame(mcmcSummary$confint$fixef)$lower
	coefs$HPDupper<-as.data.frame(mcmcSummary$confint$fixef)$upper
	coefs<-coefs[,c(1,2,5,6,3,4)]
	return(coefs)
}
###########################################################################

names(hormoneData)[24:25]<-c("cort","test")
palette(c("gray50","black"))
regT<-hormoneData$test^2
plot(test~cort,data=hormoneData,pch=16,col=as.numeric(Status),cex=1.5)
model<-lmer(test~cort+Status+(1|tank),data=hormoneData)
model2<-lmer(test~cort*Status+(1|tank),data=hormoneData)
#(hormoneData$test)^2
modelT<-lmer(test~cort+(1|tank),data=hormoneData[hormoneData$Status=="t",])
modelNT<-lmer(test~cort+(1|tank),data=hormoneData[hormoneData$Status=="nt",])

mcmcMod1<-mcmcsamp(model,10000)
mcmcMod2<-mcmcsamp(model2,10000)
mcmcModT<-mcmcsamp(modelT,10000)
mcmcModNT<-mcmcsamp(modelNT,10000)

mcmc1pval<-get.pvalues.mcmc(mcmcMod1)
mcmc2pval<-get.pvalues.mcmc(mcmcMod2)
mcmcTpval<-get.pvalues.mcmc(mcmcModT)
mcmcNTpval<-get.pvalues.mcmc(mcmcModNT)

table1<-makeNiceTable(mcmc1pval,summary(model))
table2<-makeNiceTable(mcmc2pval,summary(model2))
table3<-makeNiceTable(mcmcTpval,summary(modelT))
table4<-makeNiceTable(mcmcNTpval,summary(modelNT))

#xtable(table1)

eff<-allEffects(model2)
str(eff)
corteff<-eff$ cort
effects<-data.frame(fit=corteff$fit,lower=corteff$lower,upper=corteff$upper,status=corteff$x[,2],cort=corteff$x[,1])
lines(effects$cort[1:10],effects$fit[1:10],col="gray50",lwd=3)
lines(effects$cort[1:10],effects$upper[1:10],col="gray50",lty="dotted")
lines(effects$cort[1:10],effects$lower[1:10],col="gray50",lty="dotted")

lines(effects$cort[11:20],effects$fit[11:20],col="black",lwd=3)
lines(effects$cort[11:20],effects$upper[11:20],col="black",lty="dotted")
lines(effects$cort[11:20],effects$lower[11:20],col="black",lty="dotted")

pdf(height=5,width=5,file="t-c interaction.pdf")
plot2py<-function(){
plot(test~cort,data=hormoneData,pch=16,col=as.numeric(Status),cex=1.5,xlab="cortisol (residuals)",ylab="sqrt(testosterone) ng/mL")
lines(effects$cort[1:10],effects$fit[1:10],col="gray50",lwd=3)
lines(effects$cort[1:10],effects$upper[1:10],col="gray50",lty="dotted")
lines(effects$cort[1:10],effects$lower[1:10],col="gray50",lty="dotted")

lines(effects$cort[11:20],effects$fit[11:20],col="black",lwd=3)
lines(effects$cort[11:20],effects$upper[11:20],col="black",lty="dotted")
lines(effects$cort[11:20],effects$lower[11:20],col="black",lty="dotted")
}

save.image("tcCor.RDATA")

