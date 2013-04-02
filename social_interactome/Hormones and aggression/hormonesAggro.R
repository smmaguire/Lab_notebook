library(lme4); library(effects)

#functions
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
cort.status.samp<-mcmcsamp(model,10000)
p.values<-get.pvalues.mcmc(cort.status.samp)
return(makeNiceTable(p.values,summary(model)))
}

########################################################################
save(list=c("hormoneData",""

palette(c("gray50","black"))
plot(cort~inFromT,col=as.numeric(Status),pch=16,data=hormoneData[hormoneData$inFromT<60,])
cortMo1<-lmer(cort~inFromT*Status+SL+(1|tank),data=hormoneData[hormoneData$inFromT<60,])
fullTable(cortMo1)

eff<-allEffects(cortMo1,xlevels=list(inFromT=seq(from=0,to=60,length.out=10)))
STeff<-eff$inFromT
effects<-data.frame(fit=STeff$fit,lower=STeff$lower,upper=STeff$upper,status=STeff$x[,2],cort=STeff$x[,1])
lines(effects$cort[effects$status=="nt"],effects$fit[effects$status=="nt"],lwd=3,col="#7F7F7F90")
lines(effects$cort[effects$status=="nt"],effects$lower[effects$status=="nt"],lwd=3,col="#7F7F7F90",lty="dotted")
lines(effects$cort[effects$status=="nt"],effects$upper[effects$status=="nt"],lwd=3,col="#7F7F7F90",lty="dotted")
lines(effects$cort[effects$status=="t"],effects$fit[effects$status=="t"],lwd=3,col="#00000090")
lines(effects$cort[effects$status=="t"],effects$lower[effects$status=="t"],lwd=3,col="#00000090",lty="dotted")
lines(effects$cort[effects$status=="t"],effects$upper[effects$status=="t"],lwd=3,col="#00000090",lty="dotted")

cortMo2<-lmer(cort~inFromT+SL+(1|tank),data=hormoneData[hormoneData$inFromT<60&hormoneData$Status=="t",])

tankStat<-vector()
for(i in 1:length(hormoneData[,1])){
tank<-as.character(hormoneData$tank[i])
status<-as.character(hormoneData$Status[i])
tankStat<-c(tankStat,paste(tank,status))
}
tankStat<-factor(tankStat)

plot(cort~tankStat,data=hormoneData)
levels(hormoneData)
lm1<-lm(cort~factor(tank)*Status,data=hormoneData)
summary(lm1)
(anova(lm1))

names(hormoneData)
#out/ in / status / SL / CF/Weight/GSI
#out2t/out2NT/out2F/ SL/CF/weight/GSI/inFromT/ in
mo1<-lmer(cort~1+(1|tank),data=hormoneData[hormoneData$inFromT<60,])
mo2<-lmer(cort~1+SL+inFromNT+Status+(1|tank),data=hormoneData[hormoneData$inFromT<60,])
summary(mo2)
anova(mo1,mo2)
plot(inFromT~inFromNT,data=hormoneData,col=as.numeric(Status))

palette(c("gray50","black"))
plot(cort~CFsquare,data=hormoneData[hormoneData$inFromT<60,],col=as.numeric(Status),pch=21,cex=2)
palette(rainbow(9))
points(cort~CFsquare,data=hormoneData[hormoneData$inFromT<60,],col=as.numeric(tank),pch=16)

x<-aov(cort~factor(tank),data=hormoneData)
plot(x$residuals,hormoneData$cort[!is.na(hormoneData$cort)])
plot(hormoneData$CFsquare[!is.na(hormoneData$cort)],x$residuals,col=as.numeric(hormoneData$tank),pch=16)
summary(lm(x$residuals~hormoneData$SL[!is.na(hormoneData$cort)]))

#lenZ<-vector()
#for(j in 1:length(hormoneData$tank)){
#subdata<-hormoneData[hormoneData$tank==hormoneData$tank[j]&hormoneData$Status==hormoneData$Status[j],]
#mean<-mean(subdata$SL,na.rm=TRUE)
#stdev<-sd(subdata$SL,na.rm=TRUE)
#z<-((hormoneData$SL[j]-mean)/stdev)
#lenZ<-c(lenZ,z)
#}
##cfVar<-vector()
###lenVar<-vector()
####weightVar<-vector()
####lenZ<-vector()
###for(j in 1:length(hormoneData$tank)){
###subdata<-hormoneData[hormoneData$tank==hormoneData$tank[j]&hormoneData$Status==hormoneData$Status[j],]
###mean<-mean(subdata$CFsquare,na.rm=TRUE)
###stdev<-sd(subdata$CFsquare,na.rm=TRUE)
###z<-((hormoneData$CFsquare[j]-mean)/stdev)
###lenZ<-c(lenZ,z)
###cfVar<-c(cfVar,var(subdata$CFsquare,na.rm=TRUE))
###lenVar<-c(lenVar,var(subdata$SL,na.rm=TRUE))
##weightVar<-c(weightVar,var(subdata$weight,na.rm=TRUE))
#}

#cfVar<-vector()
#lenVar<-vector()
#weightVar<-vector()
#lenZ<-vector()
#for(j in 1:length(hormoneData$tank)){
#subdata<-hormoneData[hormoneData$tank==hormoneData$tank[j],]
#mean<-mean(subdata$CFsquare,na.rm=TRUE)
#stdev<-sd(subdata$CFsquare,na.rm=TRUE)
#z<-((hormoneData$CFsquare[j]-mean)/stdev)
#lenZ<-c(lenZ,z)
#cfVar<-c(cfVar,var(subdata$CFsquare,na.rm=TRUE))
#lenVar<-c(lenVar,var(subdata$SL,na.rm=TRUE))
#weightVar<-c(weightVar,var(subdata$weight,na.rm=TRUE))
#}
