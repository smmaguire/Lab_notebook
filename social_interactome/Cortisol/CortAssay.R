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
#dil<-rep(30,48); dil[12]<-42
#cort("3-1-2013_cort.txt",dil,names=names)
cort<-function(txtfilename,dil,names=rep(NA,48)){

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
points(y=master$net_OD[12:48],x=master$log_concent[12:48],pch=16,col="#1A985085")
master$ng_ml<-((10^master$log_concent)*dil)/1000
master$log_concent[1:11]<-NA
master$ng_ml[1:11]<-NA
master$sampNames<-names
print(master)
print(summary(model))
fit<-1-(deviance(model)/sum((y-mean(y))^2))
#write.csv(master,paste(filename,".csv",sep=""))
text(3.5,.55,labels=paste("Quality of fit","=",fit))
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
