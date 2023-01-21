if (!require("pacman")) install.packages("pacman")
pacman::p_load(bipartite,reshape, dplyr, MASS, performance, pscl, car)


data<-read.csv("data/occurences.csv", sep=",", header=T, row.names = 1)
fisio<- read.csv('data/fisio.csv')


fisio$cols<- viridis::viridis(nrow(fisio), option = 'magma')
organsP<-fisio %>%
  filter(P.g>0) 

names(organsP)[1] = names(fisio)[1] = 'organ'

data2<- data[rownames(data) %in% organsP$organ, colnames(data) %in% organsP$organ]

S_feat<- data.frame(organ=rownames(data2), sMI=rowSums(data2)/sum(data2), S_k = apply(data2, 1, function(i) sum(i > 0)))
A_feat<- data.frame(organ=colnames(data2), aMI=colSums(data2)/sum(data2), A_k = apply(data2, 2, function(i) sum(i > 0)))

organs_netfeat = merge(S_feat, A_feat, by='organ')

merge(organs_netfeat, fisio[, c('organ', 'P.g')], by='organ') %>% 
  filter(organ != 'Bone', P.g>0) %>% 
  glm.nb(formula=S_k ~ P.g) %>% # the response variable can be sMI, aMI, A_k, or S_k
  summary()






overall<-sum(data2)
dataBSI2<-data2/overall


Pmatrix<- dataBSI2
for (i in 1:nrow(Pmatrix)) {
  for (j in 1:ncol(Pmatrix)) {
    nameSource<-rownames(Pmatrix)[i]
    nameAcceptor<- colnames(Pmatrix)[j] 
    Pdis<- (organsP[organsP$organ== nameSource,]$P.g - organsP[organsP$organ== nameAcceptor,]$P.g)
    Pmatrix[i,j]<-(Pdis)
  }
}

dataBSI2<-as.matrix(dataBSI2)
Pmatrix<-as.matrix(Pmatrix)

#cor.test(dataBSI2,Pmatrix, method = "kendall", exact = NULL)
#cor.test(dataBSI2,Pmatrix, method = "spearman")


BSI_melt<-reshape::melt(dataBSI2)
Pdiff_melt<-reshape::melt(Pmatrix)#; a2$Pst<- (a2$value-min(a2$value))/(max(a2$value)-min(a2$value))

BSI_melt$Pdif<- Pdiff_melt$value
colnames(BSI_melt)<-c('Source', 'Acceptor', 'Incidence', 'Pdiff')
BSI_melt$Pdiff<-BSI_melt$Pdiff*1000

mod<-merge(BSI_melt, organsP, by.x = 'Acceptor', by.y = 'organ') %>%
  filter(Source != 'Bone', Acceptor != 'Bone') %>%
  mutate(Source=factor(Source),Incidence = Incidence * overall)%>%
  zeroinfl(formula= Incidence ~ Pdiff)


sum_mod = summary(mod)
r2_zeroinflated(mod, method = c("default"))



require(quantreg)

BSI_melt_excBone = BSI_melt %>%
  filter(Acceptor != 'Bone') %>%
  filter(Source!= 'Bone')

randx = rnorm(nrow(BSI_melt_excBone[BSI_melt_excBone$Incidence>0,]), mean(BSI_melt_excBone[BSI_melt_excBone$Incidence>0,]$Pdiff), sqrt(var(BSI_melt_excBone[BSI_melt_excBone$Incidence>0,]$Pdiff)))


Dat.nls <- nls(Incidence ~ SSlogis(Pdiff, Asym, mid, scal), data=BSI_melt_excBone[BSI_melt_excBone$Incidence>0,])
Dat.nlrq.lower <- nlrq(Incidence ~ SSlogis(Pdiff, Asym, mid, scal), data=BSI_melt_excBone[BSI_melt_excBone$Incidence>0,], tau=0.1, trace=TRUE)
Dat.nlrq.upper <- nlrq(Incidence ~ SSlogis(Pdiff, Asym, mid, scal), data=BSI_melt_excBone[BSI_melt_excBone$Incidence>0,], tau=0.9, trace=TRUE)
Dat.nlrq.median <- nlrq(Incidence ~ SSlogis(Pdiff, Asym, mid, scal), data=BSI_melt_excBone[BSI_melt_excBone$Incidence>0,], tau=0.50, trace=TRUE)
