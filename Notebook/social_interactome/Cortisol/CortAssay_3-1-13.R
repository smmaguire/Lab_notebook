######################################################################################
#Cortisol Assay v1 2/20/2013 Sean Maguire 
######################################################################################
#functions are listed in the end so read them in first the first time you use it!!
#if you make changes please save as your OWN FILE do not save over this file please!
#Run this after reading in all of the functions:  cort("my file name.txt")
######################################################################################
#dilution vector (change here if necessary) samples start at row 10.
#for example: If your first sample was diluted less than 30 times (lets say 15), then you would set row 10 (first sample is row 10) to 15 rather than 30:
#dil[10]<-15
#For more than one sample:
#dil[c(10,27,22)]<-15
#If all of your dilutions are different from the standard 30, you can just change the original code:
#dil<-rep(15,length(master[,1]))
#dil<-rep(30,48)
#cort("",dil)
###############################some functions#########################################
dil<-rep(30,48); dil[12]<-42
#names(load("names.R"))
cort2<-cort("3-1-2013_cort.txt",dil,names=names)
dil2<-rep(30,48); dil2[12]<-34.8; dil2[33]<-41.6; names2<-c(rep(NA,11),53,35,14,44,54,18,26,25,20,50,46,31,64,8,55,41,42,33,11,29,49,19,51,17,28,30,15,45,13,22,32,23,56,43,27,52,21); cort1<-cort("Cortisol raw 6-20-11.txt",dil2,names2)

cort<-function(txtfilename,dil,names=rep(NA,48)){
notes<-c(rep(NA,48))
filename<-txtfilename
master<-as.data.frame(matrix(nrow=48,ncol=5))
colnames(master)<-c("mean_OD","CV","net_OD","log_concent","ng_ml")
rownames(master)[1:11]<-c("blank","ta","nsb","bo","std1","std2","std3","std4","std5","std6","std7")

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
x<-log10(c(10000,5000,2500,1250,625,313,156))
y<-master$net_OD[5:11]
#A=max assmptote, B=Slope factor,C= inflection point, D=min assymptote
model=nls(y ~ d+(a-d)/(1+(x/cc)^b),start=list(a=.6,b=6,cc=3.2,d=.2))
xplot <- seq(min(x),max(x),length=100)
yplot <- (coef(model)["d"]+(coef(model)["a"]-coef(model)["d"])/(1+(xplot/coef(model)["cc"])^coef(model)["b"]))
plot(x,y,xlab="log cort concentration log(pg/mL)",ylab="net net_OD",pch=16,col="blue")
lines(xplot,yplot, lty="dotted", col="gray50",lwd=3)
master$log_concent<-coef(model)["cc"]*(((-1*coef(model)["a"]+master$net_OD)/(coef(model)["d"]-master$net_OD))^(1/coef(model)["b"]))
notes[which(master$log_concent>max(x))]<-"off the curve (high)"
notes[which(master$log_concent<min(x))]<-"off the curve (low)"
points(y=master$net_OD[12:48],x=master$log_concent[12:48],pch=16,col="#1A985085")
master$ng_ml<-((10^master$log_concent)*dil)/1000
master$log_concent[1:11]<-NA
master$ng_ml[1:11]<-NA
master$sampNames<-names
master$notes<-notes
suppressWarnings(master$seanCode<-addSCode(master$sampNames))
print(master)
print(summary(model))
fit<-1-(deviance(model)/sum((y-mean(y))^2))
#write.csv(master,paste(filename,".csv",sep=""))
text(3.5,.55,labels=paste("Quality of fit","=",fit))
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

cort1s<-cort1[!is.na(cort1$sCode),]
cort2s<-cort2[!is.na(cort2$sCode),]
fullCort<-data.frame(peterCode=c(cort1s$sampNames,cort2s$sampNames),seanCode=c(cort1s$seanCode,cort2s$seanCode),notes=c(cort1s$notes,cort2s$notes),cort_ng_ml=c(cort1s$ng_ml,cort2s$ng_ml),cv=c(cort1s$CV,cort2s$CV))

#removes peter's sample w/ high CV
fullCort_badSamp<-fullCort[c(1:24,26:65),]

#parallelism
cort1
x<-log10(c(10000,5000,2500,1250,625,313,156))
x<-c(x,x)
y<-c(cort1$net_OD[5:11],cort2$net_OD[5:11])
exp<-factor(c(rep(1,7),rep(2,7)))
newData<-data.frame(x,y,exp)

fullMo<-nls(y ~ d+(a-d)/(1+(x/cc)^(b[exp])),start=list(a=.6,b=c(6,7),cc=3.2,d=.2),data=newData)
subMO<-nls(y ~ d+(a-d)/(1+(x/cc)^b),start=list(a=.6,b=6,cc=3.2,d=.2),data=newData)

m1<-nls(y ~ d+(a-d)/(1+(x/cc)^b),start=list(a=.6,b=6,cc=3.2,d=.2),data=(newData[exp==1,]))

m2<-nls(y ~ d+(a-d)/(1+(x/cc)^b),start=list(a=.6,b=6,cc=3.2,d=.2),data=(newData[exp==2,]))


xplot<-(seq(min(x),max(x),length=100))
yplot1 <- (coef(m1)["d"]+(coef(m1)["a"]-coef(m1)["d"])/(1+(xplot/coef(m1)["cc"])^coef(m1)["b"]))

yplot2<-(coef(m2)["d"]+(coef(m2)["a"]-coef(m2)["d"])/(1+(xplot/coef(m2)["cc"])^coef(m2)["b"]))

yplotSUB<-(coef(subMO)["d"]+(coef(subMO)["a"]-coef(subMO)["d"])/(1+(xplot/coef(subMO)["cc"])^coef(subMO)["b"]))

plot2py<-function(){
plot(y~x,pch=16,col="purple",data=newData[exp==1,],xlab="log cort concentration",ylab="net_OD")
points(y~x,pch=16,col="yellow",data=newData[exp==2,])
lines(xplot,yplot1, lty="dotted", col="purple",lwd=3)
lines(xplot,yplot2, lty="dotted", col="yellow",lwd=3)
lines(xplot,yplotSUB,lty="solid", col="gray50",lwd=3)
}
CV<-vector()
for(x in newData2$x[1:7]){
	subdata<-newData2[newData2$x==x,]
	CV<-c(CV,(sd(subdata$y)/mean(subdata$y)))
}

save.image("cort.RDATA")