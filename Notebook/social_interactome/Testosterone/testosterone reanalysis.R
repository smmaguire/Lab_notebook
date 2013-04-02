#testosterone reanalysis with 4 parameter non-linear model for standard curves
######################################################################################
load("rosetta.RDATA")

txtfilename1<-"4-10-12 Testosterone assay melanocortin 1.txt"
txtfilename2<-"5-3-12 melanocortin T assay.txt"


save.image("testData.RDATA")
sampleName1<-c(rep(NA,9),21,20,27,3,8,22,30,39,18,9,2,7,4,13,31,29,23,21,16,37,1,34,11,6,5,38,17,33,28,19,35,14,25,32,24,26,40,15,36)
sampleName2<-c(rep(NA,9),21,50,43,55,52,60,49,73,75,47,59,77,48,71,56,76,44,46,74,58,68,57,64,45,72,63,61,62,51,69,70,65,42,66,67,78,54,41,79)
dil1<-rep(30,48)
dil1[c(20,21,23)]<-c(42,70,42)
dil2<-rep(30,48)
Tassay<-function(txtfilename,dil,names=rep(NA,48)){

notes<-c(rep(NA,48))
filename<-txtfilename
master<-as.data.frame(matrix(nrow=48,ncol=5))
colnames(master)<-c("mean_OD","CV","net_OD","log_concent","ng_ml")
rownames(master)[1:9]<-c("blank","ta","nsb","bo","std1","std2","std3","std4","std5")

data<-read.delim(filename,skip=3,header=FALSE)
data405<-as.matrix(data[1:8,3:14])
data595<-as.matrix(data[1:8,16:27])
data405<-dataPrePro(data405)
data595<-dataPrePro(data595)

sub<-(data405-data595)
master$mean_OD<-getMeanOD(sub)
#assumes you are following the std Hofmann lab protocol of putting the replicates one below the other
master$CV<-getCV(sub)
#master$percent_bound<-getBound(master)
master$net_OD<-getNetOD(master$mean_OD,master)
x<-log10(c(2000,500,125,31.25,7.81))
y<-master$net_OD[5:9]
#A=max assmptote, B=Slope factor,C= inflection point, D=min assymptote
model=nls(y ~ d+(a-d)/(1+(x/cc)^b),start=list(a=.7,b=7,cc=2.5,d=.1))
xplot <- seq(min(x),max(x),length=100)
yplot <- (coef(model)["d"]+(coef(model)["a"]-coef(model)["d"])/(1+(xplot/coef(model)["cc"])^coef(model)["b"]))
plot(x,y,xlab="log T concentration log(pg/mL)",ylab="net net_OD",pch=16,col="blue")
lines(xplot,yplot, lty="dotted", col="gray50",lwd=3)
master$log_concent<-coef(model)["cc"]*(((-1*coef(model)["a"]+master$net_OD)/(coef(model)["d"]-master$net_OD))^(1/coef(model)["b"]))
points(y=master$net_OD[12:48],x=master$log_concent[12:48],pch=16,col="#1A985085")
master$ng_ml<-((10^master$log_concent)*dil)/1000
master$log_concent[1:9]<-NA
master$ng_ml[1:9]<-NA
master$sampNames<-names
suppressWarnings(master$seanCode<-addSCode(master$sampNames))
notes[which(master$log_concent>max(x))]<-"off the curve (high)"
notes[which(master$log_concent<min(x))]<-"off the curve (low)"
master$notes<-notes
print(master)
print(summary(model))
fit<-1-(deviance(model)/sum((y-mean(y))^2))
#write.csv(master,paste(filename,".csv",sep=""))
text(2.5,.55,labels=paste("Quality of fit","=",fit))
return(master)
}

dataPrePro<-function(pdata){
  colCount<-1
  cdata<-matrix(nrow=length(pdata[,1]),ncol=length(pdata[1,]))
  while(colCount<=length(pdata[1,])){
    for(i in 1:length(pdata[,1])){
      cdata[i,colCount]<-as.numeric(pdata[i,colCount])
    }
    colCount<-colCount+1
  }
  cdata
}

getMeanOD<-function(subD2){
  subData<-as.vector(subD2[,])
  vecCount<-1
 stc<-1
  meanOD<-rep(NA,48)
    while(stc<=96){
      meanOD[vecCount]<-mean(subData[stc:(stc+1)])
      stc<-stc+2
      vecCount<-vecCount+1
      }
  meanOD
}

getCV<-function(subD2){
  subData<-as.vector(subD2[,])
  vecCount<-1
  stc<-1
  cv<-rep(NA,48)
  while(stc<=96){
    cv[vecCount]<-100*((sd(subData[stc:(stc+1)]))/(mean(subData[stc:(stc+1)])))
    stc<-stc+2
    vecCount<-vecCount+1
  }
  cv
}

#getBound<-function(master){
#  subBlank<-master[,1]-master[1,1]
#  net<-subBlank-subBlank[3]
#  bound<-100*(net/net[4]) 
#  bound
#}

getNetOD<-function(vector,master){
  blank<-master[1,1]
  subBlank<-(vector-blank)
  NSB<-subBlank[3]
  netOD<-(subBlank-NSB)
}

addSCode<-function(pCodeVector){
	sCode<-vector()
	pCodeV<-as.numeric(pCodeVector)
	pCodeV[pCodeV>64]<-NA
	for(i in pCodeV){
		if(!is.na(i))(sCode<-c(sCode,rosetta$sCode[rosetta$pCode==i])) else (sCode<-c(sCode,NA))
	}
	return(sCode)
}

t1<-Tassay("4-10-12 Testosterone assay melanocortin 1.txt",dil1,sampleName1)
t2<-Tassay("5-3-12 melanocortin T assay.txt",dil2,sampleName2)


intraplate<-100*(sd(c(t1[10,5],t2[10,5]))/(mean(c(t1[10,5],t2[10,5]))))
interplate<-mean(c(t1$CV[10:48],t2$CV[10:48]))

length(t1[11:48,1]) #38

t2.inSample<-t2[!is.na(t2$seanCode),]
rownames(t2.inSample) <- seq(length=nrow(t2.inSample))
length(t2.inSample[2:24,1]) #23

length(unique(c(t1$seanCode[11:48],t2.inSample$seanCode[2:24])))
length(c(t1$seanCode[11:48],t2.inSample$seanCode[2:24]))
duplicated(c(t1$seanCode[11:48],t2.inSample$seanCode[2:24]))

fullT<-as.data.frame(matrix(ncol=length(names(t1)),nrow=61))
names(fullT)<-names(t1)
fullT[1:38,]<-t1[11:48,]
fullT[39:61,]<-t2.inSample[2:24,]
fullT<-fullT[order(fullT$seanCode),]

unique(c(fullT$seanCode,1:64))
fullT<-rbind(fullT,c(NA,NA,NA,NA,NA,NA,9,NA))
fullT<-rbind(fullT,c(NA,NA,NA,NA,NA,NA,13,NA))
fullT<-rbind(fullT,c(NA,NA,NA,NA,NA,NA,56,NA))

fullT<-fullT[order(fullT$seanCode),]


setwd("/Users/seanmaguire/Desktop/sean lab notebook/Testosterone")
list.files()
load("cort.RDATA")
fullCort[29,]
setwd("/Users/seanmaguire/Desktop/Figures/fig2/hormone analysis/")

fullCort2<-fullCort[c(1:24,26:63,65),]
fullCort2$peterCode<-as.numeric(as.character(fullCort2$peterCode))
fullCort2<-rbind(fullCort2,c(47,44,NA,NA,NA))
fullCort2[order(fullCort2$peterCode),]
fullCort2<-fullCort2[order(fullCort2$peterCode),]
blt<-read.csv("bloodTime.csv",header=T)
blt$CODE==fullCort2$peterCode
fullCort2$bloodTime<-blt[,2]
qqnorm(log10(fullCort2$cort_ng_ml))
shapiro.test(log10(fullCort2$cort_ng_ml))
plot(log10(fullCort2$cort_ng_ml)~fullCort2$bloodTime,pch=16)
#plot(fullCort2$bloodTime~factor(fullCort2$tank))
#plot(fullCort2$cortResiduals~factor(fullCort2$tank))

bloodMo<-lm(log10(fullCort2$cort_ng_ml)~fullCort2$bloodTime)
abline(bloodMo)
summary(bloodMo)
residuals<-bloodMo$residuals
names(residuals)
#residuals["1"]
for(name in names(residuals)){
fullCort2$cortResiduals[fullCort2$peterCode==as.numeric(name)]<-residuals[name]
}

fullCort2<-fullCort2[order(fullCort2$seanCode),]
fullCort2$Status<-alldata[alldata$day==10,]$Status
fullCort2$tank<-alldata[alldata$day==10,]$tank
#bloodMo2<-lmer(log10(cort_ng_ml)~bloodTime*Status+(1|tank),data=fullCort2) just wanted to see if there was an interaction w/ status. there isn't


tcCor<-data.frame(cortResid=fullCort2$cortResiduals,sqrtT=sqrt(fullT$ng_ml),regT=fullT$ng_ml)
tcCor<-tcCor[is.na(fullCort2$notes)&is.na(fullT$notes),]
plot(tcCor$cortResid,tcCor$sqrtT,ylim=c(0,9),xlim=c(-.8,.8),pch=16)
abline(lm(sqrtT~cortResid,data=tcCor))


setwd("/Users/seanmaguire/Desktop/Figures/")
load(list.files()[7])
setwd("/Users/seanmaguire/Desktop/sean lab notebook/")
save(alldata,file="alldata.RDATA")
rm(list=ls())
load("alldata.RDATA")

fullCort2[!is.na(fullCort2$notes),]
fullCort2$cortResiduals[!is.na(fullCort2$notes)]<-NA
fullT$sqrtT[!is.na(fullT$notes)]<-NA

names(alldata)
alldata$ID[alldata$day==10]==fullCort2$seanCode
names(alldata)[24:25]
alldata[alldata$day==10,24]<-fullCort2$cortResiduals

alldata$ID[alldata$day==10]==fullT$seanCode
alldata[alldata$day==10,25]<-sqrt(fullT$sqrtT)
alldata[alldata$day==10,24:25]

hormoneDat<-alldata[alldata$day==10,]
hormoneDat$cortNotes<-fullCort2$notes
hormoneDat$testNotes<-fullT$notes
names(hormoneDat)[24:25]<-c("cort","test")
par(mfrow=c(1,2))
plot(test~cort,data=hormoneDat[is.na(hormoneDat$cortNotes)&is.na(hormoneDat$testNotes),],pch=16,col=as.numeric(Status))
plot(test~cort,data=hormoneDat,pch=16,col=as.numeric(Status))

library(lme4)
mo1<-lmer(test~cort*Status+(1|tank),data=hormoneDat)
mo2<-lmer(test~cort+Status+(1|tank),data=hormoneDat[is.na(hormoneDat$cortNotes)&is.na(hormoneDat$testNotes),])

library(effects)
eff<-allEffects(mo1)
str(eff)
corteff<-eff$ cort
effects<-data.frame(fit=corteff$fit,lower=corteff$lower,upper=corteff$upper,status=corteff$x[,2],cort=corteff$x[,1])
lines(effects$cort[1:10],effects$fit[1:10],col="black")
lines(effects$cort[1:10],effects$upper[1:10],col="black",lty="dotted")
lines(effects$cort[1:10],effects$lower[1:10],col="black",lty="dotted")

lines(effects$cort[11:20],effects$fit[11:20],col="red")
lines(effects$cort[11:20],effects$upper[11:20],col="red",lty="dotted")
lines(effects$cort[11:20],effects$lower[11:20],col="red",lty="dotted")

plot()