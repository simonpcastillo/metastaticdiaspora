if (!require("pacman")) install.packages("pacman")
pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)

data<-read.csv("data/occurences.csv", sep=",", header=T, row.names = 1)
fisio<- read.csv('data/fisio.csv', check.names = F)
names(fisio)[1] = 'organ'

sim<- 1000
mets<-10000
n.cores=20



#Case 1: cumscld and dP

fisiosimm<-fisio %>%
  filter(P.g>0 & cumlscd>0)%>%
  mutate(probPrim=cumlscd/sum(cumlscd))

occ<- data[rownames(data) %in% fisiosimm$organ, colnames(data) %in% fisiosimm$organ]

sMI<- data.frame(Source=rownames(occ), sMI=rowSums(occ)/sum(occ))

overall<-sum(occ)
dataBSI3<-occ/overall


Pmatrix3<- dataBSI3
for (i in 1:nrow(Pmatrix3)) {
  for (j in 1:ncol(Pmatrix3)) {
    nameSource<-rownames(Pmatrix3)[i]
    nameAcceptor<- colnames(Pmatrix3)[j] 
    Pdis<- (fisiosimm[fisiosimm$organ== nameSource,]$P.g - fisiosimm[fisiosimm$organ== nameAcceptor,]$P.g)
    
    Pmatrix3[i,j]<-Pdis
  }
}

dataBSI3<-as.matrix(dataBSI3)
Pmatrix3<-as.matrix(Pmatrix3)

Pmatrix3b<- Pmatrix3
for(prow in 1:nrow(Pmatrix3b)){  
  for(pcol in 1:ncol(Pmatrix3b)){
    if(min(Pmatrix3b[prow, ])<0){Pmatrix3b[prow, pcol]<- 1- (Pmatrix3[prow, pcol] + min(Pmatrix3[prow, ])*-1)}
  }
}
rowsums<-rowSums(Pmatrix3b)
for(prow in 1:nrow(Pmatrix3b)){  
  for(pcol in 1:ncol(Pmatrix3b)){
    Pmatrix3b[prow, pcol]<- Pmatrix3b[prow,pcol]/ rowsums[prow] 
    }
}






nulltopo_1<- data.frame()
cl <- makePSOCKcluster(n.cores)
registerDoParallel(cl)  # use multicore, set to the number of our cores

nulltopo_1<-  foreach (nullnets=1:sim, .combine = rbind) %dopar% {
#for(nullnets in 1:sim){
  pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)
  
  print(nullnets/sim*100)
  nullBSI3<- matrix(data = 0,ncol=ncol(Pmatrix3), nrow=nrow(Pmatrix3))
  rownames(nullBSI3)<-rownames(dataBSI3);colnames(nullBSI3)<-colnames(dataBSI3)
  
  #case A
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI3),1)
    rcol<-sample(1:ncol(dataBSI3),1)
    pPrim<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probPrim
    
    if(pPrim>runif(1)){nullBSI3[rrow, rcol]<-nullBSI3[rrow, rcol]+1}#& Pmatrix3b[rrow, rcol]>runif(1)
  }
  nullBSI3<- nullBSI3/sum(nullBSI3)
  
  gm<- computeModules(nullBSI3)
  
  simnesttemp1a<- nestedtemp(nullBSI3)$statistic
  simnestnodf1a<- nestednodf(nullBSI3)$statistic[3]
  simnestwine1a<-wine(nullBSI3)$wine
  simmod1a<-gm@likelihood 
  
  obstemp1a<- nestedtemp(dataBSI3)$statistic
  obsnodf1a<- nestednodf(dataBSI3)$statistic[3]
  obswine1a<- wine(dataBSI3,nreps = 999)$wine
  obsmod1a<-computeModules(dataBSI3)@likelihood

  fit<-lm(c(nullBSI3)~c(dataBSI3))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI1a<- cor.test(nullBSI3, dataBSI3,method = 'spearman')$p.value
  corBSI1a<-cor.test(nullBSI3, dataBSI3,method = 'spearman')$estimate
  pval01a<-summary(fit)$coefficients[2,4]
  
  
  
  #case B
  nullBSI3<- matrix(data = 0,ncol=ncol(Pmatrix3), nrow=nrow(Pmatrix3))
  rownames(nullBSI3)<-rownames(dataBSI3);colnames(nullBSI3)<-colnames(dataBSI3)
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI3),1)
    rcol<-sample(1:ncol(dataBSI3),1)
    pPrim<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probPrim
    
    if(Pmatrix3b[rrow, rcol]*pPrim>runif(1)){nullBSI3[rrow, rcol]<-nullBSI3[rrow, rcol]+1}#& Pmatrix3b[rrow, rcol]>runif(1)
  }
 
  nullBSI3<- nullBSI3
  
  gm<- computeModules(nullBSI3)
  
  simnesttemp1b<- nestedtemp(nullBSI3)$statistic
  simnestnodf1b<- nestednodf(nullBSI3)$statistic[3]
  simnestwine1b<-wine(nullBSI3)$wine
  simmod1b<-gm@likelihood 
  
  obstemp1b<- nestedtemp(dataBSI3)$statistic
  obsnodf1b<- nestednodf(dataBSI3)$statistic[3]
  obswine1b<- wine(dataBSI3,nreps = 999)$wine
  obsmod1b<-computeModules(dataBSI3)@likelihood

  fit<-lm(c(nullBSI3)~c(dataBSI3))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI1b<- cor.test(nullBSI3, dataBSI3,method = 'spearman')$p.value 
  corBSI1b<-cor.test(nullBSI3, dataBSI3,method = 'spearman')$estimate
  pval01b<-summary(fit)$coefficients[2,4]
  
  
  df0= data.frame(nullnets,simnesttemp1a, simnesttemp1b,
  simnestnodf1a, simnestnodf1b,
  simnestwine1a, simnestwine1b,
  simmod1a,simmod1b,
  corBSI1a,corBSI1b,
  pcorBSI1a, pcorBSI1b,
  pval01a, pval01b,
  obstemp1a, obstemp1b,
  obsnodf1a, obsnodf1b,
  obswine1a, obswine1b,
  obsmod1a, obsmod1b)
  

}

stopCluster(cl)





#Case 2 dP and blood flow

#######Simulation DP and flux

fisiosimm<-fisio %>%
  filter(P.g>0 & flow>0)%>%
  mutate(probFlow=flow/sum(flow)) #() (blood.flow.g-min(blood.flow.g))/(max(blood.flow.g)-min(blood.flow.g))

occ<- data[rownames(data) %in% fisiosimm$organ, colnames(data) %in% fisiosimm$organ]

sMI<- data.frame(Source=rownames(occ), sMI=rowSums(occ)/sum(occ))

overall<-sum(occ)

dataBSI3<-occ/overall


Pmatrix3<- dataBSI3
for (i in 1:nrow(Pmatrix3)) {
  for (j in 1:ncol(Pmatrix3)) {
    nameSource<-rownames(Pmatrix3)[i]
    nameAcceptor<- colnames(Pmatrix3)[j] 
    Pdis<- (fisiosimm[fisiosimm$organ== nameSource,]$P.g - fisiosimm[fisiosimm$organ== nameAcceptor,]$P.g)
    
    Pmatrix3[i,j]<-Pdis
  }
}

dataBSI3<-as.matrix(dataBSI3)
Pmatrix3<-as.matrix(Pmatrix3)

Pmatrix3b<- Pmatrix3
for(prow in 1:nrow(Pmatrix3b)){  
  for(pcol in 1:ncol(Pmatrix3b)){
    if(min(Pmatrix3b[prow, ])<0){Pmatrix3b[prow, pcol]<- 1- (Pmatrix3[prow, pcol] + min(Pmatrix3[prow, ])*-1)}
  }
}
rowsums<-rowSums(Pmatrix3b)
for(prow in 1:nrow(Pmatrix3b)){  
  for(pcol in 1:ncol(Pmatrix3b)){
    Pmatrix3b[prow, pcol]<- Pmatrix3b[prow,pcol]/ rowsums[prow] 
  }
}

nullBSI3<- matrix(data = 0,ncol=ncol(Pmatrix3), nrow=nrow(Pmatrix3))
rownames(nullBSI3)<-rownames(dataBSI3);colnames(nullBSI3)<-colnames(dataBSI3)

nulltopo_2<- data.frame()

cl <- makePSOCKcluster(n.cores)
registerDoParallel(cl)  # use multicore, set to the number of our cores

nulltopo_2<-  foreach (nullnets=1:sim, .combine = rbind) %dopar% {
  pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)
  print(nullnets)
  nullBSI3<- matrix(data = 0,ncol=ncol(Pmatrix3), nrow=nrow(Pmatrix3))
  rownames(nullBSI3)<-rownames(dataBSI3);colnames(nullBSI3)<-colnames(dataBSI3)
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI3),1)
    rcol<-sample(1:ncol(dataBSI3),1)
    pFlow<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probFlow
    
    if(pFlow>runif(1) ){nullBSI3[rrow, rcol]<-nullBSI3[rrow, rcol]+1} #& Pmatrix3b[rrow, rcol]>runif(1)
    
  }
  nullBSI3<- nullBSI3/sum(nullBSI3)
  
  gm<- computeModules(nullBSI3)
  
  
  simnesttemp2a<- nestedtemp(nullBSI3)$statistic
  simnestnodf2a<- nestednodf(nullBSI3)$statistic[3]
  simnestwine2a<-wine(nullBSI3, nreps = 999)$wine
  simmod2a<- gm@likelihood
  obstemp2a<- nestedtemp(dataBSI3)$statistic
  obsnodf2a<- nestednodf(dataBSI3)$statistic[3]
  obswine2a<- wine(dataBSI3,nreps = 999)$wine
  obsmod2a<-computeModules(dataBSI3)@likelihood
  fit<-lm(c(nullBSI3)~c(dataBSI3))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI2a<- cor.test(nullBSI3, dataBSI3,method = 'spearman')$p.value
  corBSI2a<-cor.test(nullBSI3, dataBSI3,method = 'spearman')$estimate
  pval02a<-sfit$coefficients[2,4]

  
  
  #Case B
  nullBSI3<- matrix(data = 0,ncol=ncol(Pmatrix3), nrow=nrow(Pmatrix3))
  rownames(nullBSI3)<-rownames(dataBSI3);colnames(nullBSI3)<-colnames(dataBSI3)
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI3),1)
    rcol<-sample(1:ncol(dataBSI3),1)
    pFlow<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probFlow
    
    if(pFlow*Pmatrix3b[rrow, rcol]>runif(1) ){nullBSI3[rrow, rcol]<-nullBSI3[rrow, rcol]+1} 
    
  }
  gm<- computeModules(nullBSI3)
  
  nullBSI3<- nullBSI3/sum(nullBSI3)
  simnesttemp2b<- nestedtemp(nullBSI3)$statistic
  simnestnodf2b<- nestednodf(nullBSI3)$statistic[3]
  simnestwine2b<-wine(nullBSI3,nreps = 999)$wine
  obstemp2b<- nestedtemp(dataBSI3)$statistic
  obsnodf2b<- nestednodf(dataBSI3)$statistic[3]
  obswine2b<- wine(dataBSI3,nreps = 999)$wine
  obsmod2b<-computeModules(dataBSI3)@likelihood
  simmod2b<-  gm@likelihood
  
  fit<-lm(c(nullBSI3)~c(dataBSI3))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)

  corBSI2b<-cor.test(nullBSI3, dataBSI3,method = 'spearman')$estimate
  pcorBSI2b<-cor.test(nullBSI3, dataBSI3,method = 'spearman')$p.val
  pval02b<-sfit$coefficients[2,4]
  
  
  df0= data.frame(nullnets, simnesttemp2a, simnesttemp2b,
                  simnestnodf2a, simnestnodf2b,
                  simnestwine2a, simnestwine2b,
                  simmod2a,simmod2b,
                  corBSI2a,corBSI2b,
                  pcorBSI2a, pcorBSI2b,
                  pval02a, pval02b,
                  obstemp2a, obstemp2b,
                  obsnodf2a, obsnodf2b,
                  obswine2a, obswine2b,
                  obsmod2a, obsmod2b)
  
  
};beepr::beep(sound = 4)

stopCluster(cl)


rm(list = setdiff(ls(), c('nulltopo_1', 'nulltopo_2')))

