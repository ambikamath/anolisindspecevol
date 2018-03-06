#Examining the evolution of individual specialization in perch height and perch diameter in Anolis lizards. 
#Ambika Kamath, James Stroud, Michele Johnson

library(plyr)
library(nlme)
library(ggplot2)
library(ape)
library(geiger)
library(phytools)
library(caper)
library(RInSp)

#load data
fulldat=read.csv("indspecfull.csv")
ecomorph=read.csv("ecomorph.csv")
indsvl = read.csv("indsvl.csv")

#organizing datasets
fulldat$diam=as.numeric(as.character(fulldat$diam))
spsum=ddply(fulldat, "Species", summarize, h=mean(height, na.rm=TRUE), hSD=sd(height, na.rm=TRUE), d=mean(diam[diam<999], na.rm=TRUE), dSD=sd(diam[diam<999], na.rm=TRUE), hrange=(max(height, na.rm=TRUE)-min(height, na.rm=TRUE)), drange=(max(diam[diam<999], na.rm=TRUE)-min(diam[diam<999], na.rm=TRUE)))


#create distribution dataframes (for height, this also includes a column for mean height by individual, )
height=ddply(fulldat, c("Species", "ID"), summarize, meanh=mean(height, na.rm=TRUE), check=sum(height<2401, na.rm=TRUE), cat1=(sum(height<=50, na.rm=TRUE)), cat2=(sum(height>50&height<=100, na.rm=TRUE)), cat3=(sum(height>100&height<=150, na.rm=TRUE)),
             cat4=(sum(height>150&height<=200, na.rm=TRUE)), cat5=(sum(height>200&height<=250, na.rm=TRUE)), 
             cat6=(sum(height>250&height<=300, na.rm=TRUE)), cat7=(sum(height>300&height<=350, na.rm=TRUE)), cat8=(sum(height>350, na.rm=TRUE)))
diam=ddply(fulldat, c("Species", "ID"), summarize, check=sum(diam<=1001, na.rm=TRUE),cat1=sum(diam<=1, na.rm=TRUE), cat2=sum(diam>1&diam<=10, na.rm=TRUE), cat3=sum(diam>10&diam<=20, na.rm=TRUE), cat4=sum(diam>20&diam<=40, na.rm=TRUE), 
           cat5=sum(diam>40, na.rm=TRUE))

#remove any individuals with 0 observations
height=height[height$check!=0,]   
diam=diam[diam$check!=0,]
#summarize sample sizes
hsample=ddply(height, "Species", summarise, repmeanh=mean(check), repminh=min(check), repmaxh=max(check))
dsample=ddply(diam, "Species", summarise, repmeand=mean(check), repmind=min(check), repmaxd=max(check))
#remove column for checking and then merge with ecomorph to add ecomorph classifications
#keeping "check" in for now to see if weighting by sample size is possible in mixed effects models
#height$check=NULL
height=merge(ecomorph, height, by="Species")
#diam$check=NULL
diam=merge(ecomorph, diam, by="Species")

ecomorph=merge(ecomorph, spsum, by="Species")
ecomorph=merge(ecomorph, hsample, by="Species")
ecomorph=merge(ecomorph, dsample, by="Species")



#for calculating individual specialization
n=length(levels(fulldat$Species))

splist=levels(fulldat$Species)

#height
indspech=NULL
indspechSD=NULL
indspechP=NULL
Nh=NULL
indspeclisth=list()


for (i in 1:n){
  temp1=subset(height, height$Species==splist[i])
temp2 = import.RInSp(temp1, row.names = 1, info.cols = c(1:7))
PSi = PSicalc(temp2, exclude = FALSE, pop.diet="sum", replicates = 10000)
indspeclisth[[i]]=PSi$PSi
  indspech[i]=mean(PSi$PSi)
  indspechSD[i]=sd(PSi$PSi)
  indspechP[i]=PSi$IS.pvalue
  Nh[i]=nrow(temp1)
}

v1=unlist(indspeclisth)

height$indspec=v1

ecomorph$indspech=indspech
ecomorph$indspechSD=indspechSD
ecomorph$indspechP=indspechP
ecomorph$Nh=Nh

#diameter
indspecd=NULL
indspecdSD=NULL
indspecdP=NULL
Nd=NULL
indspeclistd=list()


for (i in 1:n){
  temp1=subset(diam, diam$Species==splist[i])
  temp2 = import.RInSp(temp1, row.names = 1, info.cols = c(1:6))
  PSi = PSicalc(temp2, exclude = FALSE, pop.diet="sum", replicates = 10000)
  indspeclistd[[i]]=PSi$PSi
  indspecd[i]=mean(PSi$PSi)
  indspecdSD[i]=sd(PSi$PSi)
  indspecdP[i]=PSi$IS.pvalue
  Nd[i]=nrow(temp1)
}

v2=unlist(indspeclistd)

diam$indspec=v2

ecomorph$indspecd=indspecd
ecomorph$indspecdSD=indspecdSD
ecomorph$indspecdP=indspecdP
ecomorph$Nd=Nd

#to split by sex
#merge height and diam datasets
all=merge(height, diam, by=c("ID", "Species", "Ecomorph"))

#merge height, diam with sex, SVL
indsvl$SVL=as.numeric(as.character(indsvl$SVL))
indsvl=subset(indsvl, indsvl$Sex=="M"|indsvl$Sex=="F")


full=merge(fulldat, indsvl, by=c("ID", "Species"))
fullheight=merge(indsvl, height, by=c("ID", "Species"))
fulldiam=merge(indsvl, diam, by=c("ID", "Species"))
fullsp=ddply(full, c("Sex", "Species"), summarize, svl=mean(SVL, na.rm=TRUE), svlSD=sd(SVL, na.rm=TRUE), h=mean(height, na.rm=TRUE), 
             hSD=sd(height, na.rm=TRUE), d=mean(diam[diam<999], na.rm=TRUE), dSD=sd(diam[diam<999], na.rm=TRUE),hrange=(max(height, na.rm=TRUE)-min(height, na.rm=TRUE)), drange=(max(diam[diam<999], na.rm=TRUE)-min(diam[diam<999], na.rm=TRUE)))

dathm=subset(fullheight, fullheight$Sex=="M")
dathm$indspec=NULL
dathf=subset(fullheight, fullheight$Sex=="F")
dathf$indspec=NULL

datdm=subset(fulldiam, fulldiam$Sex=="M")
datdm$indspec=NULL
datdf=subset(fulldiam, fulldiam$Sex=="F")
datdf$indspec=NULL

#redoing individual specialization calculations for males and females separately. 
indspechm=NULL
indspechSDm=NULL
indspechPm=NULL
Nhm=NULL
indspeclisthm=list()

for (i in 1:n){
  temp1=subset(dathm, dathm$Species==splist[i])
  temp2 = import.RInSp(temp1, row.names = 1, info.cols = c(1:9))
  PSi = PSicalc(temp2, exclude = FALSE, pop.diet="sum", replicates = 10000)
  indspechm[i]=mean(PSi$PSi)
  indspechSDm[i]=sd(PSi$PSi)
  indspechPm[i]=PSi$IS.pvalue
  Nhm[i]=nrow(temp1)
  indspeclisthm[[i]]=PSi$PSi
}

v3=unlist(indspeclisthm)

dathm$indspec=v3

indspechf=NULL
indspechSDf=NULL
indspechPf=NULL
Nhf=NULL
indspeclisthf=list()

for (i in 1:n){
  temp1=subset(dathf, dathf$Species==splist[i])
  temp2 = import.RInSp(temp1, row.names = 1, info.cols = c(1:9))
  PSi = PSicalc(temp2, exclude = FALSE, pop.diet="sum", replicates = 10000)
  indspechf[i]=mean(PSi$PSi)
  indspechSDf[i]=sd(PSi$PSi)
  indspechPf[i]=PSi$IS.pvalue
  Nhf[i]=nrow(temp1)
  indspeclisthf[[i]]=PSi$PSi
}

v4=unlist(indspeclisthf)

dathf$indspec=v4

fullheightbysex=rbind(dathm, dathf)


fullsp$indspech=0
fullsp$indspech[1:15]=indspechf
fullsp$indspech[16:30]=indspechm
fullsp$indspechSD=0
fullsp$indspechSD[1:15]=indspechSDf
fullsp$indspechSD[16:30]=indspechSDm
fullsp$indspechP=0
fullsp$indspechP[1:15]=indspechPf
fullsp$indspechP[16:30]=indspechPm
fullsp$Nh=0
fullsp$Nh[1:15]=Nhf
fullsp$Nh[16:30]=Nhm

#split by sex for diameter
indspecdm=NULL
indspecdSDm=NULL
indspecdPm=NULL
Ndm=NULL
indspeclistdm=list()

for (i in 1:n){
  temp1=subset(datdm, datdm$Species==splist[i])
  temp2 = import.RInSp(temp1, row.names = 1, info.cols = c(1:8))
  PSi = PSicalc(temp2, exclude = FALSE, pop.diet="sum", replicates = 10000)
  indspecdm[i]=mean(PSi$PSi)
  indspecdSDm[i]=sd(PSi$PSi)
  indspecdPm[i]=PSi$IS.pvalue
  Ndm[i]=nrow(temp1)
  indspeclistdm[[i]]=PSi$PSi
}

v5=unlist(indspeclistdm)

datdm$indspec=v5

indspecdf=NULL
indspecdSDf=NULL
indspecdPf=NULL
Ndf=NULL
indspeclistdf=list()

for (i in 1:n){
  temp1=subset(datdf, datdf$Species==splist[i])
  temp2 = import.RInSp(temp1, row.names = 1, info.cols = c(1:8))
  PSi = PSicalc(temp2, exclude = FALSE, pop.diet="sum", replicates = 10000)
  indspecdf[i]=mean(PSi$PSi)
  indspecdSDf[i]=sd(PSi$PSi)
  indspecdPf[i]=PSi$IS.pvalue
  Ndf[i]=nrow(temp1)
  indspeclistdf[[i]]=PSi$PSi
}

v6=unlist(indspeclistdf)

datdf$indspec=v6

fulldiambysex=rbind(datdm, datdf)


fullsp$indspecd=0
fullsp$indspecd[1:15]=indspecdf
fullsp$indspecd[16:30]=indspecdm
fullsp$indspecdSD=0
fullsp$indspecdSD[1:15]=indspecdSDf
fullsp$indspecdSD[16:30]=indspecdSDm
fullsp$indspecdP=0
fullsp$indspecdP[1:15]=indspecdPf
fullsp$indspecdP[16:30]=indspecdPm
fullsp$Nd=0
fullsp$Nd[1:15]=Ndf
fullsp$Nd[16:30]=Ndm

fullsphd=merge(fullsp, ecomorph[,1:2], by="Species", all.x=TRUE)

fullsphdm=subset(fullsphd, fullsphd$Sex=="M")
fullsphdf=subset(fullsphd, fullsphd$Sex=="F")
fullsex=merge(fullsphdm, fullsphdf, by=c("Species", "Ecomorph"))



#phylogenetic tree set up
#NO NEED TO RE-RUN THIS PART UNLESS SPECIES CHANGE
#Read tree from RAxML
ComparativeTree <- read.nexus("Anolis_216_species.nex.con")
ComparativeTree <- ComparativeTree[[1]]
ComparativeTree
pruned.tree<-drop.tip(ComparativeTree,ComparativeTree$tip.label[-match(splist, ComparativeTree$tip.label)])
write.tree(pruned.tree)
write.nexus(pruned.tree, file="prunedtree.nex")

#rename node labels by hand
#change all "1.00"  to 16:22


#START HERE
pruned.tree <- read.nexus("prunedtree.nex")
plot(pruned.tree)


#align.tip.label used to align all tips to the same point
#like in the previous pruned tree
plot.phylo(pruned.tree, cex=1,show.tip.label= TRUE ,type="phylogram", label.offset=.01, align.tip.label=TRUE)

#for phylogenetically corrected correlations
comp1 <- comparative.data(phy=pruned.tree, data=ecomorph, names.col = Species, vcv = TRUE)
comp2 <- comparative.data(phy=pruned.tree, data=fullsex, names.col = Species, vcv = TRUE)

#to prep for phytools
ecomorphfactor=as.factor(ecomorph$Ecomorph)
names(ecomorphfactor)=ecomorph$Species
indspech=ecomorph$indspech
names(indspech)=ecomorph$Species
indspecd=ecomorph$indspecd
names(indspecd)=ecomorph$Species

#ANALYSES by question

#1. Is there phylogenetic signal in the evolution of individual specialization?
#Run the analysis for Bloomberg's K and test for significance
phylosig(pruned.tree, indspech, method="K", nsim=10000, test=TRUE)
phylosig(pruned.tree, indspecd, method="K", nsim=10000, test=TRUE)

#plot ind spec on phylogeny 
plotBranchbyTrait(pruned.tree,mode="tips", indspech, palette="gray", xlims=c(0.25,0.75))
plotBranchbyTrait(pruned.tree,mode="tips", indspecd, palette="gray", xlims=c(0.5,1))


#2.	Does individual specialization vary by ecomorph, after accounting for phylogeny?

#height
#model looking at effect of ecomorph on individual specialization, no phylogenetic correction
mod1=lme(indspec~Ecomorph, random=~1|Species, data=height, na.action=na.omit)
summary(mod1)

#with phylogenetic correction (on species means)
aov.phylo(indspech~ecomorphfactor, pruned.tree, n=10000)

#boxplot
plot1=ggplot(data=height)
plot1+theme_classic()+geom_boxplot(aes(x=Ecomorph, y=indspec, fill=Ecomorph, color=Species))+
  scale_color_manual(values = rep("black", 15), guide=FALSE)+xlab("Ecomorph")+ylab("Degree of Individual Specialization \n More Specialized                    Less Specialized")+
  scale_x_discrete(labels=c("Generalist", "Grass Bush", "Trunk Crown", "Trunk Ground", "Twig"))+
  scale_fill_manual(values=c("grey", "#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"))


#diameter
mod2=lme(indspec~Ecomorph, random=~1|Species, data=diam, na.action=na.omit)
summary(mod2)
aov.phylo(indspecd~ecomorphfactor, pruned.tree, n=10000)

plot2=ggplot(data=diam)
plot2+theme_classic()+geom_boxplot(aes(x=Ecomorph, y=indspec, fill=Ecomorph, color=Species))+
  scale_color_manual(values = rep("black", 15), guide=FALSE)+xlab("Ecomorph")+ylab("Degree of Individual Specialization \n More Specialized                    Less Specialized")+
  scale_x_discrete(labels=c("Generalist", "Grass Bush", "Trunk Crown", "Trunk Ground", "Twig"))+
  scale_fill_manual(values=c("grey","#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"))

#correlation between height and height specialization. 
mod3=lme(indspec~meanh, random=~1|Species, data=height, na.action=na.omit)
summary(mod3)

#correcting for phylogeny
fit1 <- pgls(indspech~h, comp1, lambda="ML")
summary(fit1)

#plot
plot3=ggplot(ecomorph)
plot3+theme_classic()+geom_point(aes(x=h, y=indspech, color=Ecomorph), size=4)+
  scale_color_manual(values=c("grey","#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"))+
  xlab("Mean Perch Height (in cm)")+ylab("Mean Individual Specialization \n More Specialized                    Less Specialized")+
  ylim(c(0.25,0.75))

#3.	Is individual specialization in height related to individual specialization in diameter?
mod4=lme(indspec.x~indspec.y, random=~1|Species, data=all, na.action=na.omit)
summary(mod4)

fit2 <- pgls(indspech~indspecd, comp1, lambda="ML")
summary(fit2)


#4.	Is individual specialization of a species related to its
#   niche breadth (i.e. range of perch height or diameter)?
fit3 <- pgls(indspech~hrange, comp1, lambda="ML")
summary(fit3)

fit4 <- pgls(indspecd~drange, comp1, lambda="ML")
summary(fit4)

plot3+theme_classic()+geom_point(aes(x=hrange, y=indspech, color=Ecomorph), size=4)+
  scale_color_manual(values=c("grey","#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"))+
  xlab("Perch Height Range (in cm)")+ylab("Mean Individual Specialization \n More Specialized                    Less Specialized")
  

plot3+theme_classic()+geom_point(aes(x=drange, y=indspecd, color=Ecomorph), size=4)+
  scale_color_manual(values=c("grey","#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"))+
  xlab("Perch Diameter Range (in cm)")+ylab("Mean Individual Specialization \n More Specialized                    Less Specialized")



#5.	Are species with more congeners more likely to be specialists?

fit5 <- pgls(indspech~Congeners, comp1, lambda="ML")
summary(fit5)

fit6 <- pgls(indspecd~Congeners, comp1, lambda="ML")
summary(fit6)

#6.	Does individual specialization in height or diameter differ between the sexes?
mod5=lme(indspec~Sex, random=~1|Ecomorph/Species, data=fullheightbysex)
summary(mod5)
anova(mod5)

mod6=lme(indspec~Sex, random=~1|Ecomorph/Species, data=fulldiambysex)
summary(mod6)
anova(mod6)


plot4=ggplot(fullheightbysex)
plot4+theme_classic()+geom_boxplot(aes(x=Sex, y=indspec, fill=Ecomorph))+
  scale_fill_manual(values = c("grey", "#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"), name="Ecomorph")+
  ylab("Individual Height Specialization \n More Specialized                    Less Specialized")+xlab("Sex")

plot5=ggplot(fulldiambysex)
plot5+theme_classic()+geom_boxplot(aes(x=Sex, y=indspec, fill=Ecomorph))+
  scale_fill_manual(values = c("grey", "#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"), name="Ecomorph")+
  ylab("Individual Diameter Specialization \n More Specialized                    Less Specialized")+xlab("Sex")


#7.	Is individual specialization in height or diameter correlated between the sexes,
#after accounting for phylogeny?
#only done on species-sex means
mod7=lme(indspech.x~indspech.y, random=~1|Ecomorph, data=fullsex)
summary(mod7)
mod8=lme(indspecd.x~indspecd.y, random=~1|Ecomorph, data=fullsex)
summary(mod8)

fit7 <- pgls(indspech.x~indspech.y, comp2, lambda="ML")
summary(fit7)
fit8 <- pgls(indspecd.x~indspecd.y, comp2, lambda="ML")
summary(fit8)

 
plot6=ggplot(fullsex)
plot6+theme_classic()+geom_point(aes(x=indspech.x, y=indspech.y, colour=Ecomorph), size=4)+
  geom_abline(intercept=0, slope=1, linetype=2)+xlim(c(0.25, 1))+ylim(c(0.25,1))+
  scale_color_manual(values = c("grey", "#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"), name="Ecomorph")+
  xlab("Male Individual Height Specialization")+
  ylab("Female Individual Height Specialization")

plot6+theme_classic()+geom_point(aes(x=indspecd.x, y=indspecd.y, color=Ecomorph), size=4)+
  geom_abline(intercept=0, slope=1, linetype=2)+xlim(c(0.4, 1))+ylim(c(0.4,1))+
  scale_color_manual(values = c("grey", "#00c6c3ff", "#FFC981","#e24344ff" , "#9cd78aff"), name="Ecomorph")+
  ylab("Individual Female Diameter Specialization \n More Specialized                          Less Specialized")+
  xlab(" More Specialized                          Less Specialized \n Individual Male Diameter Specialization")




#CODE FOR RANDOMIZATION ACROSS SYMPATRY GROUPS, BUT I'M NOT SURE THIS IS BIOLOGICALLY INTERESTING. 
#pairwise diffs in height specialization within sympatry groups.

pairs=list()

for (i in 1:length(levels(ecomorph$Sympatry))){
  pairs[[i]]=if(nrow(subset(ecomorph, ecomorph$Sympatry==levels(ecomorph$Sympatry)[i]))==1)  c("NA", "NA") else combn(ecomorph$Species[ecomorph$Sympatry==levels(ecomorph$Sympatry)[i]], 2)
}

pairs2=data.frame(matrix(unlist(pairs), nrow=14, byrow=T))
pairs2$sp1=splist[as.numeric(as.character(pairs2$X1))]
pairs2$sp2=splist[as.numeric(as.character(pairs2$X2))]

#for calculating differences in individual specialization for sympatric species in original dataset
for (i in 1:nrow(pairs2)){
  pairs2$diff[i]=abs(ecomorph$indspech[ecomorph$Species==pairs2$sp1[i]]-ecomorph$indspech[ecomorph$Species==pairs2$sp2[i]])
}



temp=ecomorph
randdiff=NULL
set.seed(123)
#randomizing
for (j in 1:1000){
  temp$Sympatry=sample(temp$Sympatry)
  rpairs=list()
  for (i in 1:length(levels(temp$Sympatry))){
    rpairs[[i]]=if(nrow(subset(temp, ecomorph$Sympatry==levels(temp$Sympatry)[i]))==1)  c("NA", "NA") else combn(temp$Species[temp$Sympatry==levels(temp$Sympatry)[i]], 2)
  }
  
  rpairs2=data.frame(matrix(unlist(rpairs), nrow=14, byrow=T))
  rpairs2$sp1=splist[as.numeric(as.character(rpairs2$X1))]
  rpairs2$sp2=splist[as.numeric(as.character(rpairs2$X2))]
  for (i in 1:nrow(rpairs2)){
    rpairs2$diff[i]=abs(temp$indspech[temp$Species==rpairs2$sp1[i]]-temp$indspech[temp$Species==rpairs2$sp2[i]])
  }
 randdiff[j]=mean(rpairs2$diff, na.rm=TRUE)
}

length(randdiff[randdiff>mean(pairs2$diff, na.rm=TRUE)])/length(randdiff)
#p-value of 0.076 


