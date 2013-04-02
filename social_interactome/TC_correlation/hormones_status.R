load("hormonesxStatus.RDATA")

save(list=c("hormoneData","get.pvalues.mcmc","makeNiceTable","fullTable","cort.status","test.status","cort.status.fullTable","test.status.fullTable","plot2py"),file="hormonesxStatus.RDATA")


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
fullTable<-function(model){
mod.samp<-mcmcsamp(model,10000)
mcmcSum<-get.pvalues.mcmc(mod.samp)
return(makeNiceTable(mcmcSum,summary(model)))
}

test.status<-lmer(test~Status+(1|tank),data=hormoneData)
cort.status<-lmer(cort~Status+(1|tank),data=hormoneData)
cort.status.fullTable<-fullTable(cort.status)
test.status.fullTable<-fullTable(test.status)

plot2py<-function(){
plot(hormoneData$cort~hormoneData$Status,xlab="Status",ylab="Cortisol residuals")
plot(hormoneData$test~hormoneData$Status,xlab="Status",ylab="sqrt(testosterone) ng/mL")
print(cort.status.fullTable)
print(xtable(cort.status.fullTable))
print(test.status.fullTable)
print(xtable(test.status.fullTable))
}

pdf(height=4,width=3,file="cortisolxstatus.pdf")
dev.off()


pdf(height=4,width=3,file="testxstatus.pdf")
dev.off()



