#   MorphoTools is a set of R functions for morphometric analysis created by Petr Koutecký
#   Koutecký, P. Plant Syst Evol (2015) 301: 1115. doi:10.1007/s00606-014-1153-2
#   
#   Customized for analysis of Sericomyrmex ants morphological dataset by Ana Jesovnik for 
#   Jesovnik & Schultz: Revision of the fungus-farming ant genus Sericomyrmex Mayr (Hymenoptera, Formicidae, Myrmicinae). ZooKeys.
# 
#  Source file:MorphoTools.R

#Set wd
setwd("~/Desktop/R_morpho")


#Read in the data of the full dataset:
indiv.orig <- read.morfodata("~/Desktop/R_morpho/Morpho_data.txt")
str(indiv.orig) #this shows me data structure
summary(indiv.orig) #this gives data summary
d.sum<-summary(indiv.orig)
export.res(d.sum, "DataSummary.txt")

#Data summary and structure outputs are good for checking the data for errors: unexpected symbols (?, /) and weird outliers

#Descriptiove stats:
d.all<-descr.all(indiv.orig)
d.tax<-descr.tax(indiv.orig)
d.pop<-descr.pop(indiv.orig)
export.res(d.pop,"dpop.txt") 
export.res(d.tax,"dtax.txt") 
export.res(d.all,"dall.txt") 

#Boxplots:
#Boxplot za HWe
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(HWe~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="HWe: Head Width (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#Boxplot za HW
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(HW~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="HWe: Head Width (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#IFW1
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(IFW1~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="IFW1: Inter Frontal Width (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#IFW2
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(IFW2~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="IFW2: Inter Frontal Width 2 (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#HL1
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(HL1~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="HL1: Head length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#HL2
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(HL2~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="HL2: Head Length 2 (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#ScL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(ScL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="SL: Scape length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#EL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(EL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="EL: Eye length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#Om
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(Om~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="Om: Ommatidia Count",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#WL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(WL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="WL: Weber length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#PL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(PL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="PL: Petiole length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#PPL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(PPL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="PPL: Postpetiole length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values


#GL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(GL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="GL: Gaster length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#HFL
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(HFL~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="HFL: Hind Femur length (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values


#PrW
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(PrW~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="PW: Pronotum Width (mm)",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#Frontal Lobe Index
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(FLI~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="FLI: Frontal Lobe Index",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#Cephalic Index
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(CI~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="CI: Cephalic Index",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#SI
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(SI~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="SI: Scapus Index",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#CEI
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(CEI~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="CEI",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#hw1
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(hw1~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="hw1: posterior head width",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#hw2
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(hw2~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="hw2: mid head width",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#hw3
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(hw3~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="hw3: anterior head width",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values

#IOD
#males only
default.par<-par(no.readonly=T) # stores the default values of par
par(mar=c(4.5,4.5,1,1)) # bottom,left,top,right margin width (as no. of lines)
boxplot(IOD~Taxon,data=indiv.orig,cex.axis=1,cex.lab=1,range=1.5,col="yellow",
        xlab="Taxon",ylab="IOD: inter ocular distance",pars=list(whisklty=1,boxwex=0.7))
par(default.par) # restores the default values


# POPULATION MEANS calculates the average value for each character in each population, 
# with the pairwise deletion of missing data. The resulting data frame has populations 
# as rows, and the first three columns are Population, Taxon, and N (the number of complete rows for each population).

popul.orig<-popul.means(indiv.orig)
export.res(popul.orig, "PopulationMeans.txt")

# TRANSFORMATIONS

#To test for normal distribution, Shapiro Wil
shapiro.test(indiv.orig$HWe)
shapiro.test(indiv.orig$HW)
shapiro.test(indiv.orig$IFW1)
shapiro.test(indiv.orig$IFW2)
shapiro.test(indiv.orig$HL1)
shapiro.test(indiv.orig$HL2)
shapiro.test(indiv.orig$ScL)
shapiro.test(indiv.orig$EL)
shapiro.test(indiv.orig$Om)
shapiro.test(indiv.orig$WL)
shapiro.test(indiv.orig$PL)
shapiro.test(indiv.orig$PPL)
shapiro.test(indiv.orig$GL)
shapiro.test(indiv.orig$HFL)
shapiro.test(indiv.orig$PrW)
shapiro.test(indiv.orig$hw1)
shapiro.test(indiv.orig$hw2)
shapiro.test(indiv.orig$hw3)
shapiro.test(indiv.orig$CI)
shapiro.test(indiv.orig$FLI)
shapiro.test(indiv.orig$SI)
shapiro.test(indiv.orig$EI)
shapiro.test(indiv.orig$CEI)


# based on the Shapiro-Wilk, all characters for Sericomyrmex except CI(p=0.06) and SI (p=0.3) are not normal, so data needs to be transformed
# in Shapiro-Wilk the nul hypothesis is that distribution is normal. if p is really small (less than 0.05 or what ever is choosen) 
# than hypothesis is rejected an our distribution is not normal. If the p is higher than 0.5, it is normal.

#Log transform data set:

indiv<-indiv.orig  #creates new data set, indiv, which we then transform, and the old indiv.orig stays the same
indiv$HWe<-log10(indiv$HWe)
indiv$HW<-log10(indiv$HW)
indiv$IFW1<-log10(indiv$IFW1)
indiv$IFW2<-log10(indiv$IFW2)
indiv$HL1<-log10(indiv$HL1)
indiv$HL2<-log10(indiv$HL2)
indiv$ScL<-log10(indiv$ScL)
indiv$EL<-log10(indiv$EL)
indiv$Om<-log10(indiv$Om)
indiv$WL<-log10(indiv$WL)
indiv$PL<-log10(indiv$PL)
indiv$PPL<-log10(indiv$PPL)
indiv$GL<-log10(indiv$GL)
indiv$HFL<-log10(indiv$HFL)
indiv$PrW<-log10(indiv$PrW)
indiv$hw1<-log10(indiv$hw1)
indiv$hw2<-log10(indiv$hw2)
indiv$hw3<-log10(indiv$hw3)
indiv$CI<-log10(indiv$CI)
indiv$FLI<-log10(indiv$FLI)
indiv$SI<-log10(indiv$SI)
indiv$EI<-log10(indiv$EI)
indiv$CEI<-log10(indiv$CEI)

#to export the data frame with all variables transformed 
export.res(indiv,"Morpho_data_log.txt") 

#The same is done for popul.orig data set
#popul.orig is a data set with just one row per population, mean value for that taxon for each character(variable)
popul<-popul.otu(popul.orig) 
popul$HWe<-log10(popul$HWe)
popul$HW<-log10(popul$HW)
popul$IFW1<-log10(popul$IFW1)
popul$IFW2<-log10(popul$IFW2)
popul$HL1<-log10(popul$HL1)
popul$HL2<-log10(popul$HL2)
popul$ScL<-log10(popul$ScL)
popul$EL<-log10(popul$EL)
popul$Om<-log10(popul$Om)
popul$WL<-log10(popul$WL)
popul$PL<-log10(popul$PL)
popul$PPL<-log10(popul$PPL)
popul$GL<-log10(popul$GL)
popul$HFL<-log10(popul$HFL)
popul$PrW<-log10(popul$PrW)
popul$hw1<-log10(popul$hw1)
popul$hw2<-log10(popul$hw2)
popul$hw3<-log10(popul$hw3)
popul$CI<-log10(popul$CI)
popul$FLI<-log10(popul$FLI)
popul$SI<-log10(popul$SI)
popul$EI<-log10(popul$EI)
popul$CEI<-log10(popul$CEI)

export.res(popul,"popul.log.txt")
#after running this my popul data set is log transformed

# CORRELATION OF CHARACTERS 
correl.p<-cormat.p(indiv)
correl.s<-cormat.s(indiv)
export.res(correl.p, "CorrelationP.txt")
export.res(correl.s, "CorrelationS.txt")

#Based on results of Spearman  i removed the following (correlated) variables: HW, HL1, HL2, SCL, WL, HFL, hw1, hw2
#I used the Spearman because: "When the variables are not normally distributed or the relationship between 
#the variables is not linear, it may be more appropriate to use the Spearman rank correlation method" 

# PRINCIPAL COMPONENT ANALYSIS (PCA)
#I run PCA on log transformed, reduced data set
#Reducing dataset, taking out all characters that are correlated:

#Removes columns 5,8, etc from the indiv dataset (indiv is the log transformed original data set)
indiv.red<-indiv[,-c(5,8,9,13,17,19,20)] 

#Removes columns 5,8, etc from the popul dataset (popul is the log transformed data set)
popul.red<-popul[,-c(5,8,9,13,17,19,20)]

# PCA 
indiv.pca<-pca.calc(indiv.red)
summary(indiv.pca)
indiv.eigen<-pca.eigen(indiv.pca) # eigenvalues
summary(indiv.eigen)
plot(indiv.pca) # scree plot 
plot(indiv.pca,type="lines") # scree plot that shows how much variation is acounted in each PC

indiv.pca.co<-pca.cor(indiv.pca, 6) # correlation of characters and ordination axies, so called character loading
#I have put 6, the default is 4

indiv.pca.sc<-pca.scores(indiv.pca,indiv.red) # computes the ordination scores for objects tj. OTU-s
#which means this compute the PCA value for each individual, for each variable

export.res(indiv.pca.sc,"indiv.pca.sc.txt") # this is the table with PCA value for each individual, for each variable
export.res(indiv.pca.co,"indiv_pca_co.txt") # this is the with character loadings, the correlation of each characters with the cach of the PCA. 
#The higher the value more important that character (variable) is for that PCA

# scatterplot of individuals
# setting symbol types and colours
indiv.pca.sc$Symb<-as.numeric(indiv.pca.sc$Taxon)
indiv.pca.sc$Col<-as.numeric(indiv.pca.sc$Taxon)
levels(indiv.pca.sc$Taxon)

# scatterplot
par(mar=c(4.5,4.5,2,2))
with(indiv.pca.sc,plot(PC1,PC2,type="p",pch=Symb,col=Col,bg="grey",xlim=c(-9,8),
                       ylim=c(-4.5,6),cex=1.5,cex.axis=1,cex.lab=1,xlab="PC1",ylab="PC2"))
abline(h=0,v=0,lty=2,col="black")
legend(-8.5,6.0,levels(indiv.pca.sc$Taxon),pch=c(1,2,3,4,5,6,7,8,9,10,11),col=c(1,2,3,4,5,6,7,8,9,10,11),
       bty="n",pt.bg="grey",cex=1,pt.cex=1.2,ncol=1)
par(default.par)


# contributions of the characters (as arrows) for individuals pca
par(mar=c(4.5,4.5,1,1))
with(indiv.pca.co,plot(PC1,PC2,type="n",asp=1,xlim=c(-1,1),ylim=c(-1,1)))
abline(h=0,v=0,lty=2,col="black")
with(indiv.pca.co,arrows(0,0,PC1,PC2,length=0.15,lty=1,col="red"))
with(indiv.pca.co[c(6,12,13,17,22,23),],text(PC1,PC2,Character,pos=4,offset=0.3,cex=0.75))
with(indiv.pca.co[c(1,11,14,16,18,20,21,24),],text(PC1,PC2,Character,pos=2,offset=0.2,cex=0.75))
with(indiv.pca.co[c(4,8,9),],text(PC1,PC2,Character,pos=1,offset=0.3,cex=0.75))
with(indiv.pca.co[c(7),],text(PC1,PC2,Character,pos=1,offset=0.7,cex=0.75))
with(indiv.pca.co[c(2,3,10,19,25),],text(PC1,PC2,Character,pos=3,offset=0.3,cex=0.75))
with(indiv.pca.co[c(5,15),],text(PC1,PC2,Character,pos=3,offset=0.15,cex=0.75))
par(default.par)

# 3D plot of individuals
popul.pca.sc$Symb<-as.numeric(popul.pca.sc$Taxon)
popul.pca.sc$Col<-as.numeric(popul.pca.sc$Taxon)
levels(indiv.pca.sc$Taxon)
library(scatterplot3d)
with(indiv.pca.sc,scatterplot3d(PC1,PC2,PC4,color=Col,pch=16,type="h",
                                cex.symbols=1,cex.axis=0.75,cex.lab=0.75,lty.hplot=1,scale.y=0.9,angle=50,box=F,grid=T,
                                y.margin.add=0.3,mar=c(3,3,0,2.5),col.grid="grey",
                                xlab="PC1",ylab="PC2",zlab="PC4"))
#I changed code to use PC4 instead of PC3, it looks more interesting (SI is the highest correlated char). Doesnt make a lot of difference though
#these 3d plots in general look pretty bad for Serico

# "spider" scatterplot using ade4 package
library(ade4)
par(oma=c(2,2,1,1))
s.class(indiv.pca.sc,indiv.pca.sc$Taxon,xax=4,yax=5,cstar=1,cellipse=0,axesell=F,clabel=0.8,cpoint=0.9,
        label=levels(indiv.pca.sc$Taxon),col=c(1,2,3,4,5,6,7,8,9,10,11),pch=16,grid=F,addaxes=T)
axis(1,pos=-4.9,at=c(-8,0,5),labels=c("-8.0","0.0","5.0"),cex.axis=1,lwd=0)
axis(2,pos=-9.5,at=c(-4,0,5),labels=c("-4.0","0.0","5.0"),cex.axis=1,lwd=0)
par(default.par)


#Runing PCA of individuals on reduced datasets:
#Creating subsets of data, with only some species

#Subset 1 is the reduced original dataset, indiv.red, without the correlated variables
#Subset 2 is S. maravalhasi and S. scrobifer 
indiv.red2<-with(indiv.red,indiv.red[Taxon=="MAR"|Taxon=="SCR",])
indiv.red2<-droplevels(indiv.red2)

#Subset 3 is S. amabilis and S. opacus
indiv.red3<-with(indiv.red,indiv.red[Taxon=="AMA"|Taxon=="OPA",])
indiv.red3<-droplevels(indiv.red3)

#Subset 4 is S. amabilis and S. saussurei
indiv.red4<-with(indiv.red,indiv.red[Taxon=="AMA"|Taxon=="SAU",])
indiv.red4<-droplevels(indiv.red4)

#Subset 5 is S. saramama and S. parvulus
indiv.red5<-with(indiv.red,indiv.red[Taxon=="AZT"|Taxon=="CAP"|Taxon=="PAR",])
indiv.red5<-droplevels(indiv.red5)

#Subset 6 is S. bondari and S. radioheadi
indiv.red6<-with(indiv.red,indiv.red[Taxon=="BON"|Taxon=="RAD",])
indiv.red6<-droplevels(indiv.red6)

#Subset 7 is S. saussurei and S. mayri
indiv.red7<-with(indiv.red,indiv.red[Taxon=="SAU"|Taxon=="MAY",])
indiv.red7<-droplevels(indiv.red7)

#Subset 8 is S. parvulus and S. opacus
indiv.red8<-with(indiv.red,indiv.red[Taxon=="PAR"|Taxon=="AZT",])
indiv.red8<-droplevels(indiv.red8)

#Subset 9 large species only
indiv.red9<-with(indiv.red,indiv.red[Taxon=="BON"|Taxon=="MAY"|Taxon=="LUT"|Taxon=="AMA"|Taxon=="SAU"|Taxon=="RAD",])
indiv.red9<-droplevels(indiv.red9)

#Subset 10 smaller spacies only
indiv.red10<-with(indiv.red,indiv.red[Taxon=="SCR"|Taxon=="MAR"|Taxon=="SAR"|Taxon=="PAR"|Taxon=="OPA",])
indiv.red10<-droplevels(indiv.red10)

#Subset 11 S. opacus populations
indiv.red11<-with(indiv.red,indiv.red[Population=="OPA_1"|Population=="OPA_2"|Population=="OPA_3",])
indiv.red11<-droplevels(indiv.red11)

#Subset 12 S. bondari populations
indiv.red12<-with(indiv.red,indiv.red[Population=="BON"|Population=="BON_RH",])
indiv.red12<-droplevels(indiv.red12)

#Subset13 S. amabilis populations
indiv.red13<-with(indiv.red,indiv.red[Population=="AMA"|Population=="AMA_SM",])
indiv.red13<-droplevels(indiv.red13)

#Subset 14 S. saussurei populations
indiv.red14<-with(indiv.red,indiv.red[Population=="SAU"|Population=="SAU_SM",])
indiv.red14<-droplevels(indiv.red14)

#Subset 15 S. mayri populations
indiv.red15<-with(indiv.red,indiv.red[Population=="MAY1"|Population=="MAY2"|Population=="MAY3"|Population=="MAY4"|Population=="MAY",])
indiv.red15<-droplevels(indiv.red15)

#This is the same code as for PCA of individuals above, but with some lines deleted
#Instead of copying the code 5 times for each subset, i just change the reduced dataset name (on three places)
# PCA of individuals for Subset 15:

indiv.pca<-pca.calc(indiv.red15)
summary(indiv.pca)
pca.eigen(indiv.pca) # eigenvalues
indiv.pca.co<-pca.cor(indiv.pca) # loadings of the characters
indiv.pca.sc<-pca.scores(indiv.pca,indiv.red15) # sample scores
export.res(indiv.pca.sc,"indiv.pca.sc_red15.txt") #pca for each individual
export.res(indiv.pca.co,"indiv.pca.co_red15.txt") #table with character loadings

# scatterplot of individuals
indiv.pca.sc$Symb<-as.numeric(indiv.pca.sc$Taxon)
indiv.pca.sc$Col<-as.numeric(indiv.pca.sc$Taxon)
levels(indiv.pca.sc$Taxon)
# scatterplot
par(mar=c(4.5,4.5,2,2))
with(indiv.pca.sc,plot(PC1,PC2,type="p",pch=Symb,col=Col,bg="grey",xlim=c(-9,8),
                       ylim=c(-4.5,6),cex=1.5,cex.axis=1,cex.lab=1,xlab="PC1",ylab="PC2"))
abline(h=0,v=0,lty=2,col="black")
legend(-8.6,6.0,levels(indiv.pca.sc$Taxon),pch=c(1,2,3,4,5),col=c(1,2,3,4,5),
       bty="n",pt.bg="grey",cex=1,pt.cex=1.2,ncol=1)
par(default.par)

