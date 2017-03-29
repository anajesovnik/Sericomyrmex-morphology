# MorphoTools: a set of R functions for morphometric analysis
# 
# Plant Systematics and Evolution xxx: xxx-xxx
# 
# Petr Koutecký
# University of South Bohemia, Faculty of Science, Department of Botany, Branišovská 1760, České Budějovice, CZ-37005, Czechia
# kouta@prf.jcu.cz
# 
# Version 5 Nov 2014
#
# The function definitions (to be used as a source in R)
#
# --------------------------------------------------------------------------------
read.morfodata<-function(FILE){
  xdata<-read.delim(FILE,header=T)
  names(xdata)[1:3]<-c("ID","Population","Taxon")
  if (!is.factor(xdata$ID)) xdata$ID<-as.factor(xdata$ID)
  if (!is.factor(xdata$Population)) xdata$Population<-as.factor(xdata$Population)
  if (!is.factor(xdata$Taxon)) xdata$Taxon<-as.factor(xdata$Taxon)
  return(xdata)}
# --------------------------------------------------------------------------------
read.morfodata2<-function(FILE){
  xdata<-read.delim2(FILE,header=T)
  names(xdata)[1:3]<-c("ID","Population","Taxon")
  if (!is.factor(xdata$ID)) xdata$ID<-as.factor(xdata$ID)
  if (!is.factor(xdata$Population)) xdata$Population<-as.factor(xdata$Population)
  if (!is.factor(xdata$Taxon)) xdata$Taxon<-as.factor(xdata$Taxon)
  return(xdata)}
# --------------------------------------------------------------------------------
export.res<-function(RESULTS,file="clipboard"){
  write.table(RESULTS,file=file,sep="\t",quote=F,row.names=F,col.names=T,na="")}
# --------------------------------------------------------------------------------
descr.tax<-function(DATA){
  taxa<-levels(DATA$Taxon)
  noi<-seq(1,ncol(DATA)-3,1)
  res<-data.frame(factor(),factor(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric())
  for (i in taxa) {
    datai<-with(DATA,DATA[Taxon==i,])
    chari<-colnames(datai)[4:ncol(datai)]
    taxi<-rep(i,ncol(datai)-3)
    lengthi<-sapply(datai[,-(1:3)],function(x) length(na.omit(x)))
    meani<-sapply(datai[,-(1:3)],mean,na.rm=T)
    sdi<-sapply(datai[,-(1:3)],sd,na.rm=T)
    quantilei<-sapply(datai[,-(1:3)],function(x) quantile(x,probs=c(0.00,0.05,0.25,0.50,0.75,0.95,1.00),na.rm=T))
    resi<-data.frame(chari,taxi,lengthi,meani,sdi,t(quantilei),noi)
    res<-rbind(res,resi)}
  names(res)<-c("Character","Taxon","N","Mean","SD","Min","5%","25%","Median","75%","95%","Max","No")
  res<-with(res,res[order(No,Taxon),])
  res$Mean[which(is.nan(res$Mean))]<-NA
  res<-res[,-ncol(res)]
  return(res)}
# --------------------------------------------------------------------------------
descr.pop<-function(DATA){
  popul<-levels(DATA$Popul)
  noi<-seq(1,ncol(DATA)-3,1)
  res<-data.frame(factor(),factor(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric())
  for (i in popul) {
    datai<-with(DATA,DATA[Population==i,])
    chari<-colnames(datai)[4:ncol(datai)]
    popi<-rep(i,ncol(datai)-3)
    lengthi<-sapply(datai[,-(1:3)],function(x) length(na.omit(x)))
    meani<-sapply(datai[,-(1:3)],mean,na.rm=T)
    sdi<-sapply(datai[,-(1:3)],sd,na.rm=T)
    quantilei<-sapply(datai[,-(1:3)],function(x) quantile(x,probs=c(0.00,0.05,0.25,0.50,0.75,0.95,1.00),na.rm=T))
    resi<-data.frame(chari,popi,lengthi,meani,sdi,t(quantilei),noi)
    res<-rbind(res,resi)}
  names(res)<-c("Character","Population","N","Mean","SD","Min","5%","25%","Median","75%","95%","Max","No")
  res<-with(res,res[order(No,Population),])
  res$Mean[which(is.nan(res$Mean))]<-NA
  res<-res[,-ncol(res)]
  return(res)}
# --------------------------------------------------------------------------------
descr.all<-function(DATA){
  taxa<-levels(DATA$Taxon)
  res<-data.frame(factor(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric(),numeric())
  chari<-colnames(DATA)[4:ncol(DATA)]
  lengthi<-sapply(DATA[,-(1:3)],function(x) length(na.omit(x)))
  meani<-sapply(DATA[,-(1:3)],mean,na.rm=T)
  sdi<-sapply(DATA[,-(1:3)],sd,na.rm=T)
  quantilei<-sapply(DATA[,-(1:3)],function(x) quantile(x,probs=c(0.00,0.05,0.25,0.50,0.75,0.95,1.00),na.rm=T))
  res<-data.frame(chari,lengthi,meani,sdi,t(quantilei))
  names(res)<-c("Character","N","Mean","SD","Min","5%","25%","Median","75%","95%","Max")
  res$Mean[which(is.nan(res$Mean))]<-NA
  return(res)}
# --------------------------------------------------------------------------------
na.meansubst<-function(DATA){
  popul<-levels(DATA$Population)
  xdata<-DATA[1,]
  for (i in popul) {
    xp<-with(DATA,DATA[Population==i,])
    meansubst<-function(x){
      m<-mean(x,na.rm=T)
      if (is.nan(m)) m<-NA
      x[which(is.na(x))]<-m
      return(x)}
    xp[,-(1:3)]<-as.data.frame(sapply(xp[,-(1:3)],meansubst))
    xdata<-rbind(xdata,xp)}
  xdata<-xdata[-1,]
  return(xdata)}
# --------------------------------------------------------------------------------
popul.means<-function(DATA){
  poplist<-unique(DATA[,c(2,3)])
  xdata<-na.omit(DATA)
  xn<-aggregate(xdata[,4],by=list(Population=xdata$Population),length)
  xn<-merge(poplist,xn,all=T)
  names(xn)[3]<-"N"
  xn$N[which(is.na(xn$N))]<-0
  res<-aggregate(DATA[,-(1:3)],by=list(Population=DATA$Population),mean,na.rm=T)
  nanx<-function(x) {
    x[which(is.nan(x))]<-NA
    return(x)}
  res[,-1]<-sapply(res[,-1],nanx)
  xn<-merge(xn,res,all=T)
  return(xn)}
# --------------------------------------------------------------------------------
popul.otu<-function(DATA){
  xpop<-DATA[,(1:3)]
  xdata<-DATA[,-(1:3)]
  colnames(xpop)<-c("ID","Population","Taxon")
  xpop[,2:3]<-DATA[,1:2]
  xpop<-cbind(xpop,xdata)
  return(xpop)}
# --------------------------------------------------------------------------------
charhist<-function(DATA){
  xn<-4:ncol(DATA)
    for (i in xn) {
    xh<-as.numeric(DATA[,i])
    if (is.integer(DATA[,i]) & max(DATA[,i])==0) next
    if (is.integer(DATA[,i]) & max(DATA[,i])==1) next
    if(sum(is.na(xh))>0) xh<-xh[-which(is.na(xh))]
    xs<-seq(min(xh),max(xh),(max(xh)-min(xh))/100)
    hist(xh,freq=F,main="",col="grey",xlab=names(DATA)[i],cex.axis=0.7,cex.lab=0.8)
    lines(xs,dnorm(xs,mean(xh),sqrt(var(xh))),col="red",lty=1)
    par(ask=T)}
  par(ask=F)}
# --------------------------------------------------------------------------------
charhist.t<-function(DATA,TAXON){
  xtax<-DATA[DATA$Taxon==TAXON,]
  xn<-4:ncol(xtax)
  for (i in xn) {
    xh<-as.numeric(xtax[,i])
    if (is.integer(xtax[,i]) & max(xtax[,i])==0) next
    if (is.integer(xtax[,i]) & max(xtax[,i])==1) next
    if (sum(is.na(xh))>0) xh<-xh[-which(is.na(xh))]
    xs<-seq(min(xh),max(xh),(max(xh)-min(xh))/100)
    hist(xh,freq=F,main="",col="grey",xlab=names(DATA)[i],cex.axis=0.7,cex.lab=0.8)
    lines(xs,dnorm(xs,mean(xh),sqrt(var(xh))),col="red",lty=1)
    par(ask=T)}
  par(ask=F)}
# --------------------------------------------------------------------------------
cormat.p<-function (DATA) {
  x<-cor(DATA[,-(1:3)],use="pairwise.complete.obs",method="pearson")
  x<-data.frame(x)
  x<-data.frame(attr(x,"row.names"),x,row.names=NULL)
  names(x)[1]<-"Pearson"
  return(x)}
# --------------------------------------------------------------------------------
cormat.s<-function (DATA) {
  x<-cor(DATA[,-(1:3)],use="pairwise.complete.obs",method="spearman")
  x<-data.frame(x)
  x<-data.frame(attr(x,"row.names"),x,row.names=NULL)
  names(x)[1]<-"Spearman"
  return(x)}
# --------------------------------------------------------------------------------
clust.upgma<-function(DATA){
  xdata<-DATA[,-c(2,3)]
  xdata<-data.frame(xdata,row.names=1)
  xdata<-scale(xdata)
  xdist<-dist(xdata,method="euclidean")
  xclust<-hclust(xdist,method="average")
  return(xclust)}
# --------------------------------------------------------------------------------
clust.ward<-function(DATA){
  xdata<-DATA[,-c(2,3)]
  xdata<-data.frame(xdata,row.names=1)
  xdata<-scale(xdata)
  xdist<-dist(xdata,method="euclidean")
  xclust<-hclust(xdist,method="ward.D")
  return(xclust)}
# --------------------------------------------------------------------------------
pca.calc<-function(DATA){
  xdata<-na.omit(DATA)
  x<-prcomp(xdata[,-(1:3)],center=T,scale.=T)
  return(x)}
# --------------------------------------------------------------------------------
pca.eigen<-function(RESULTS){
  x<-as.numeric(1:length(RESULTS$sdev))
  x<-sapply(x,function(x) paste("PC",x,sep=""))
  res<-sapply(RESULTS$sdev,function(x) x^2)
  names(res)<-x
  print(res,digits=3)}
# --------------------------------------------------------------------------------
pca.cor<-function(RESULTS,N=4){
  coeff<-data.frame(RESULTS$rotation)
  coeff<-coeff[,1:N]
  coeff<-data.frame(attr(coeff,"row.names"),coeff,row.names=NULL)
  names(coeff)[1]<-"Character"
  sdev<-RESULTS$sdev
  sdev<-sdev[1:N]
  coeff[,-1]<-t(apply(coeff[,-1],1,function(x) x*sdev))
  return(coeff)}
# --------------------------------------------------------------------------------
pca.scores<-function(RESULTS,DATA,N=4){
  xdata<-na.omit(DATA)
  res<-data.frame(predict(RESULTS,xdata))
  res<-res[,1:N]
  scores<-data.frame(ID=xdata$ID,Population=xdata$Population,Taxon=xdata$Taxon,res)
  return(scores)}
# --------------------------------------------------------------------------------
recode.col<-function(DATA,LEVEL,NEWLEVEL){
  lev<-with(DATA,which(Taxon==LEVEL))
  DATA$Col[lev]<-NEWLEVEL
  return(DATA)}
# --------------------------------------------------------------------------------
recode.symb<-function(DATA,LEVEL,NEWLEVEL){
  lev<-with(DATA,which(Taxon==LEVEL))
  DATA$Symb[lev]<-NEWLEVEL
  return(DATA)}
# --------------------------------------------------------------------------------
discr.calc<-function(DATA){
  require(vegan)
  xdata<-na.omit(DATA)
  xdata<-droplevels(xdata)
  xclass<-data.frame(model.matrix(~xdata$Taxon-1),row.names=xdata$ID)
  names(xclass)<-levels(xdata$Taxon)
  xchar<-data.frame(xdata[,-(1:3)],row.names=xdata$ID)
  discr.data<<-list("class"=xclass,"char"=xchar)
  res<-cca(discr.data$class~.,data=discr.data$char)
  return(res)}
# --------------------------------------------------------------------------------
discr.sum<-function(RESULTS,perm=500){
  require(vegan)
  x<-summary(RESULTS,scaling=-2)
  eig<-eigenvals(RESULTS,constrained=T)
  eig<-sapply(eig,function(x) x/(1-x))
  evar<-x$concont
  ccoeff<-spenvcor(RESULTS)
  test<-anova(RESULTS,step=perm)
  if (length(eig)>1) atest<-anova(RESULTS,step=perm,by="axis") else atest<-NULL
  res<-list(eig,evar,ccoeff,test,atest)
  names(res)<-c("eigenvalues","explained.variation","canonical.correlation.coefficients","test","test.of.axes")
  return(res)}
# --------------------------------------------------------------------------------
discr.coef<-function(RESULTS){
  require(vegan)
  coeff<-data.frame(coef(RESULTS))
  xmean<-sapply(discr.data$char,mean)
  coeff<-data.frame(attr(coeff,"row.names"),xmean,coeff,row.names=NULL)
  colnames(coeff)[c(1,2)]<-c("Character","Mean")
  return(coeff)}
# --------------------------------------------------------------------------------
discr.taxa<-function(RESULTS){
  require(vegan)
  taxa<-data.frame(scores(RESULTS,display="species",scaling=-2,choices=1:RESULTS$tot.chi))
  taxa<-data.frame(attr(taxa,"row.names"),taxa,row.names=NULL)
  colnames(taxa)[1]<-"Taxon"
  return(taxa)}
# --------------------------------------------------------------------------------
discr.scores<-function(RESULTS,NEWDATA){
  require(vegan)
  ndata<-na.omit(NEWDATA)
  res<-predict(RESULTS,ndata,type="lc",scaling=-2)
  res<-data.frame(ID=ndata$ID,Population=ndata$Population,Taxon=ndata$Taxon,res)
  return(res)}
#--------------------------------------------------------------------------------
discr.bip<-function(RESULTS){
  require(vegan)
  eig<-eigenvals(RESULTS,constrained=T)
  xsum<-summary(RESULTS,axes=length(eig))
  res<-xsum$biplot
  res<-t(apply(res,1,function(x) x/(sqrt(1/(1-eig)))))
  f<-envfit(RESULTS,discr.data$char,choices=1:length(eig),permutations=0)
  f<-f$vectors[[2]]
  f<-sqrt(1/(1-f))  
  res<-res*f
  {if (length(eig)>1) 
    res<-data.frame("Character"=rownames(res),res,row.names=NULL)
   else
     {res<-t(res)
      res<-data.frame(res)
      res<-data.frame("Character"=row.names(res),"CCA1"=res$res,row.names=NULL)}}
  return(res)}
#--------------------------------------------------------------------------------
discr.test<-function(RESULTS,perm=500){
  require(vegan)
  cca0<-cca(discr.data$class~1,data=discr.data$char)
  mtest<-add1(cca0,scope=formula(RESULTS),test="permutation",pstep=perm,perm.max=perm)
  mtest<-data.frame("Character"=attr(mtest,"row.names"),mtest,row.names=NULL)
  utest<-anova(RESULTS,step=perm,perm.max=perm,by="margin")
  utest<-data.frame("Character"=attr(utest,"row.names"),utest,row.names=NULL)
  xtest<-list("single.characters"=mtest,"unique.contributions"=utest)
  return(xtest)}
#--------------------------------------------------------------------------------
discr.step<-function(RESULTS,perm=500,p=0.05,dir="forward"){
  require(vegan)
  cca0<-cca(discr.data$class~1,data=discr.data$char)
  res<-ordistep(cca0,scope=formula(RESULTS),direction=dir,pstep=perm,perm.max=perm,steps=100,Pin=p)
  return(res)}
#--------------------------------------------------------------------------------
classif.da<-function(DATA){
  require(MASS)
  xdata<-na.omit(DATA)
  xdata<-droplevels(xdata)
  ntax<-length(levels(xdata$Taxon))
  char<-colnames(xdata)[-c(1:3)]
  pop<-levels(xdata$Population)
  x<-data.frame(replicate(4,factor(),simplify=F))
  y<-data.frame(replicate(ntax,numeric(),simplify=F))
  res<-merge(x,y)
  names(res)<-c("ID","Population","Taxon","Classif",levels(xdata$Taxon))
  for (i in pop) {
    samp<-with(xdata,xdata[Population==i,])
    train<-with(xdata,xdata[(Population!=i),])
    lda.train<-lda(as.formula(paste("Taxon ~ ",paste(char,collapse="+"))),data=train,prior=rep(1/ntax,ntax))
    lda.samp<-predict(lda.train,samp)
    resi<-data.frame(ID=samp$ID,Population=samp$Population,Taxon=samp$Taxon,Classif=lda.samp$class,lda.samp$posterior)
    res<-rbind(res,resi)}
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
#--------------------------------------------------------------------------------
classif.da.1<-function(DATA){
  require(MASS)
  xdata<-na.omit(DATA)
  xdata<-droplevels(xdata)
  ntax<-length(levels(xdata$Taxon))
  char<-colnames(xdata)[-c(1:3)]
  lda.res<-lda(as.formula(paste("Taxon ~ ",paste(char,collapse="+"))),data=xdata,prior=rep(1/ntax,ntax),CV=TRUE)
  res<-data.frame(ID=xdata$ID,Population=xdata$Population,Taxon=xdata$Taxon,Classif=lda.res$class,lda.res$posterior)
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
# --------------------------------------------------------------------------------
classif.samp<-function(SAMPLE,TRAINING){
  require(MASS)
  samp<-na.omit(SAMPLE)
  samp<-droplevels(samp)  
  train<-na.omit(TRAINING)
  train<-droplevels(train)
  ntax<-length(levels(train$Taxon))
  char<-colnames(train)[-c(1:3)]
  lda.train<-lda(as.formula(paste("Taxon ~ ",paste(char,collapse="+"))),data=train,prior=rep(1/ntax,ntax))
  lda.samp<-predict(lda.train,samp)
  res<-data.frame(ID=samp$ID,Population=samp$Population,Taxon=samp$Taxon,Classif=lda.samp$class,lda.samp$posterior)
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
#--------------------------------------------------------------------------------
classif.matrix<-function(DATA){
  classif<-table(DATA[,3:4])
  classif<-data.frame(unclass(classif))
  classif<-data.frame(Taxon=attr(classif,"row.names"),classif,row.names=NULL)
  classif$N<-rowSums(classif[2:ncol(classif)])
  ncor<-aggregate(Correct~Taxon,data=DATA,sum)
  classif<-merge(classif,ncor)
  last<-classif[1,]
  last[1]<-"Total"
  last[2:length(last)]<-colSums(classif[2:ncol(classif)])
  classif<-rbind(classif,last)
  names(classif)[ncol(classif)]<-"percent.correct"
  classif$percent.correct<-with(classif,(percent.correct/N)*100)
  return(classif)}
#--------------------------------------------------------------------------------
classif.pmatrix<-function(DATA){
  classif<-table(DATA[,c(2,4)])
  classif<-data.frame(unclass(classif))
  classif<-data.frame(Population=attr(classif,"row.names"),classif,row.names=NULL)
  tax<-unique(DATA[,c(2,3)])
  classif<-merge(tax,classif)
  classif$N<-rowSums(classif[3:ncol(classif)])
  ncor<-aggregate(Correct~Population,data=DATA,sum)
  classif<-merge(classif,ncor)
  names(classif)[ncol(classif)]<-"percent.correct"
  classif$percent.correct<-with(classif,(percent.correct/N)*100)
  return(classif)}
# --------------------------------------------------------------------------------
knn.select<-function(DATA){
  require(class)
  k<-as.numeric(1:30)
  xdata<-na.omit(DATA)
  xdata[,-(1:3)]<-scale(xdata[,-(1:3)])
  pop<-levels(xdata$Population)
  knn.correct<-function(x){
    res<-numeric()
    for (i in pop) {
      samp<-with(xdata,xdata[Population==i,-(1:3)])      
      samp.t<-with(xdata,xdata[Population==i,3])      
      train<-with(xdata,xdata[(Population!=i),-(1:3)])
      classif<-with(xdata,xdata[(Population!=i),3])
      knn.samp<-knn(train,samp,classif,k=x)
      resi<-sum(as.character(samp.t)==as.character(knn.samp))
      res<-sum(res,resi)}
    return(res)}
  for (j in 1:10){
    kselj<-sapply(k,knn.correct)
    if (j==1) ksel<-kselj else ksel<-rbind(ksel,kselj)}
  kselmean<-apply(ksel,2,mean)
  kselmax<-apply(ksel,2,max)
  kselmin<-apply(ksel,2,min)
  plot(kselmean,type="p",pch=16,xlab="K",ylab="correct classifications",ylim=c(min(kselmin),max(kselmax)))
  sapply(k[-1],function(x) arrows(x,kselmin[x],x,kselmax[x],code=3,angle=90,length=0.07))
  cat("The optimal K is:",which(kselmean==max(kselmean)))}
# --------------------------------------------------------------------------------
knn.select.1<-function(DATA){
  require(class)
  k<-as.numeric(1:30)
  xdata<-na.omit(DATA)
  xdata[,-(1:3)]<-scale(xdata[,-(1:3)])
  knn.correct.1<-function(x){
    train<-xdata[,-(1:3)]      
    classif<-xdata[,3]
    knn.samp<-knn.cv(train,classif,k=x,prob=F,use.all=T)
    res<-sum(as.character(classif)==as.character(knn.samp))
    return(res)}
  for (j in 1:10){
    kselj<-sapply(k,knn.correct.1)
    if (j==1) ksel<-kselj else ksel<-rbind(ksel,kselj)}
  kselmean<-apply(ksel,2,mean)
  kselmax<-apply(ksel,2,max)
  kselmin<-apply(ksel,2,min)
  plot(kselmean,type="p",pch=16,xlab="K",ylab="correct classifications",ylim=c(min(kselmin),max(kselmax)))
  sapply(k[-1],function(x) arrows(x,kselmin[x],x,kselmax[x],code=3,angle=90,length=0.07))
  cat("The optimal K is:",which(kselmean==max(kselmean)))}
# --------------------------------------------------------------------------------
knn.classif<-function(DATA,K){
  require(class)
  xdata<-na.omit(DATA)
  xdata[,-(1:3)]<-scale(xdata[,-(1:3)])
  pop<-levels(xdata$Population)
  res<-data.frame(ID=factor(),Population=factor(),Taxon=factor(),Classif=factor(),Prob=numeric())
  for (i in pop) {
    samp<-with(xdata,xdata[Population==i,-(1:3)])      
    samp.a<-with(xdata,xdata[Population==i,1:3])      
    train<-with(xdata,xdata[Population!=i,-(1:3)])
    classif<-with(xdata,xdata[Population!=i,3])
    knn.samp<-knn(train,samp,classif,k=K,prob=T,use.all=T)
    resi<-data.frame(ID=samp.a$ID,Population=samp.a$Population,Taxon=samp.a$Taxon,Classif=knn.samp,Prob=attr(knn.samp,"prob"))
    res<-rbind(res,resi)}
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
# --------------------------------------------------------------------------------
knn.classif.1<-function(DATA,K){
  require(class)
  xdata<-na.omit(DATA)
  xdata[,-(1:3)]<-scale(xdata[,-(1:3)])
  train<-xdata[,-(1:3)]      
  classif<-xdata[,3]
  knn.samp<-knn.cv(train,classif,k=K,prob=T,use.all=T)
  res<-data.frame(ID=xdata$ID,Population=xdata$Population,Taxon=xdata$Taxon,Classif=knn.samp,Prob=attr(knn.samp,"prob"))
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
# --------------------------------------------------------------------------------
knn.samp<-function(SAMPLE,TRAINING,K){
  require(class)
  xsamp<-na.omit(SAMPLE)
  samp<-scale(xsamp[,-(1:3)])
  xtrain<-na.omit(TRAINING)
  train<-scale(xtrain[,-(1:3)])
  classif<-xtrain[,3]
  knn.samp<-knn(train,samp,classif,k=K,prob=T,use.all=T)
  res<-data.frame(ID=xsamp$ID,Population=xsamp$Population,Taxon=xsamp$Taxon,Classif=knn.samp,Prob=attr(knn.samp,"prob"))
  res$Correct<-as.character(res$Taxon)==as.character(res$Classif)
  return(res)}
