if (!require("pacman")) install.packages("pacman")
pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)


load('data/occurrences.rdata')
load('data/physio_df.rdata')

set.seed(08132023)
sim<- 1000
mets<-10000

# This pipeline uses parallelisation, the number of cores needs to be defined here:
n.cores=2


#Scenarios I and II: cumscld and dP

fisiosimm<-physio_df %>%
  filter(P.g>0 & cumlscd>0)%>%
  mutate(probPrim=cumlscd/sum(cumlscd))

occ<- df0[rownames(df0) %in% fisiosimm$organ, colnames(df0) %in% fisiosimm$organ]

sMI<- data.frame(Source=rownames(occ), sMI=rowSums(occ)/sum(occ))

overall<-sum(occ)
dataBSI<-occ/overall


Pmatrix<- dataBSI
for (i in 1:nrow(Pmatrix)) {
  for (j in 1:ncol(Pmatrix)) {
    nameSource<-rownames(Pmatrix)[i]
    nameAcceptor<- colnames(Pmatrix)[j] 
    Pdis<- (fisiosimm[fisiosimm$organ== nameSource,]$P.g - fisiosimm[fisiosimm$organ== nameAcceptor,]$P.g)
    
    Pmatrix[i,j]<-Pdis
  }
}

dataBSI<-as.matrix(dataBSI)
Pmatrix<-as.matrix(Pmatrix)

Pmatrixb<- Pmatrix
for(prow in 1:nrow(Pmatrixb)){  
  for(pcol in 1:ncol(Pmatrixb)){
    if(min(Pmatrixb[prow, ])<0){Pmatrixb[prow, pcol]<- 1- (Pmatrix[prow, pcol] + min(Pmatrix[prow, ])*-1)}
  }
}
rowsums<-rowSums(Pmatrixb)
for(prow in 1:nrow(Pmatrixb)){  
  for(pcol in 1:ncol(Pmatrixb)){
    Pmatrixb[prow, pcol]<- Pmatrixb[prow,pcol]/ rowsums[prow] 
  }
}


nulltopo_1<- data.frame()
cl <- makePSOCKcluster(n.cores)
registerDoParallel(cl)  # use multicore, set to the number of our cores

nulltopo_1<-  foreach (nullnets=1:sim, .combine = rbind) %dopar% {
  pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)
  
  print(nullnets/sim*100)
  nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
  rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)
  
  sc1_probmat = nullBSI
  
  for(sourceid in fisiosimm$organ){
    sc1_probmat[sourceid,] = fisiosimm[fisiosimm$organ==sourceid,'probPrim']
  }
  
  write.csv(sc1_probmat, file='data/sc1_probmatrix.csv')
  
  
  #scenario I
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI),1)
    rcol<-sample(1:ncol(dataBSI),1)
    #pPrim<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probPrim
    prob_mets <- sc1_probmat[rrow, rcol]
    if(prob_mets>runif(1)){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1}
  }
  
  nullBSI<- nullBSI/sum(nullBSI)
  gm<- computeModules(nullBSI)
  simnesttemp1<- nestedtemp(nullBSI)$statistic
  simnestnodf1<- nestednodf(nullBSI)$statistic[3]
  simnestwine1<-wine(nullBSI)$wine
  simmod1<-gm@likelihood 
  
  obstemp1<- nestedtemp(dataBSI)$statistic
  obsnodf1<- nestednodf(dataBSI)$statistic[3]
  obswine1<- wine(dataBSI,nreps = 999)$wine
  obsmod1<-computeModules(dataBSI)@likelihood
  
  fit<-lm(c(nullBSI)~c(dataBSI))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI1<- cor.test(nullBSI, dataBSI,method = 'spearman')$p.value
  corBSI1<-cor.test(nullBSI, dataBSI,method = 'spearman')$estimate
  pval01<-summary(fit)$coefficients[2,4]
  
  write.csv(data.frame(run = nullnets, melt.matrix(nullBSI), melt.matrix(dataBSI)), paste0('data/raw_cor_trueplot/','sc1_run',nullnets,'.csv'))
  
  
  #scenario II
  nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
  rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)
  
  sc2_probmat = sc1_probmat * Pmatrixb
  write.csv(sc2_probmat, file='data/sc2_probmatrix.csv')
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI),1)
    rcol<-sample(1:ncol(dataBSI),1)
    # pPrim<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probPrim
    prob_mets = sc2_probmat[rrow, rcol]
    # if(Pmatrixb[rrow, rcol]*pPrim>runif(1)){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1}
    if(prob_mets>runif(1)){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1}
    
  }
  
  nullBSI<- nullBSI/sum(nullBSI)
  gm<- computeModules(nullBSI)
  simnesttemp2<- nestedtemp(nullBSI)$statistic
  simnestnodf2<- nestednodf(nullBSI)$statistic[3]
  simnestwine2<-wine(nullBSI)$wine
  simmod2<-gm@likelihood 
  
  obstemp2<- nestedtemp(dataBSI)$statistic
  obsnodf2<- nestednodf(dataBSI)$statistic[3]
  obswine2<- wine(dataBSI,nreps = 999)$wine
  obsmod2<-computeModules(dataBSI)@likelihood
  
  fit<-lm(c(nullBSI)~c(dataBSI))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI2<- cor.test(nullBSI, dataBSI,method = 'spearman')$p.value 
  corBSI2<-cor.test(nullBSI, dataBSI, method = 'spearman')$estimate
  pval02<-summary(fit)$coefficients[2,4]
  write.csv(data.frame(run = nullnets, melt.matrix(nullBSI), melt.matrix(dataBSI)), paste0('data/raw_cor_trueplot/','sc2_run',nullnets,'.csv'))
  
  
  df0= data.frame(nullnets,simnesttemp1, simnesttemp2,
                  simnestnodf1, simnestnodf2,
                  simnestwine1, simnestwine2,
                  simmod1,simmod2,
                  corBSI1,corBSI2,
                  pcorBSI1, pcorBSI2,
                  pval01, pval02,
                  obstemp1, obstemp2,
                  obsnodf1, obsnodf2,
                  obswine1, obswine2,
                  obsmod1, obsmod2)
  
  
}

stopCluster(cl)



#Scenarios III and IV: dP and blood flow

fisiosimm<-physio_df %>%
  filter(P.g>0 & flow>0)%>%
  mutate(probFlow=flow/sum(flow))

occ<- df0[rownames(df0) %in% fisiosimm$organ, colnames(df0) %in% fisiosimm$organ]

sMI<- data.frame(Source=rownames(occ), sMI=rowSums(occ)/sum(occ))

overall<-sum(occ)

dataBSI<-occ/overall


Pmatrix<- dataBSI
for (i in 1:nrow(Pmatrix)) {
  for (j in 1:ncol(Pmatrix)) {
    nameSource<-rownames(Pmatrix)[i]
    nameAcceptor<- colnames(Pmatrix)[j] 
    Pdis<- (fisiosimm[fisiosimm$organ== nameSource,]$P.g - fisiosimm[fisiosimm$organ== nameAcceptor,]$P.g)
    Pmatrix[i,j]<-Pdis
  }
}

dataBSI<-as.matrix(dataBSI)
Pmatrix<-as.matrix(Pmatrix)

Pmatrixb<- Pmatrix
for(prow in 1:nrow(Pmatrixb)){  
  for(pcol in 1:ncol(Pmatrixb)){
    if(min(Pmatrixb[prow, ])<0){Pmatrixb[prow, pcol]<- 1- (Pmatrix[prow, pcol] + min(Pmatrix[prow, ])*-1)}
  }
}
rowsums<-rowSums(Pmatrixb)
for(prow in 1:nrow(Pmatrixb)){  
  for(pcol in 1:ncol(Pmatrixb)){
    Pmatrixb[prow, pcol]<- Pmatrixb[prow,pcol]/ rowsums[prow] 
  }
}

nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)

nulltopo_2<- data.frame()

cl <- makePSOCKcluster(n.cores)
registerDoParallel(cl)  # use multicore, set to the number of our cores

nulltopo_2<-  foreach (nullnets=1:sim, .combine = rbind) %dopar% {
  pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)
  print(nullnets)
  nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
  rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)
  
  sc3_probmat = nullBSI
  
  for(acceptorid in fisiosimm$organ){
    sc3_probmat[,acceptorid] = fisiosimm[fisiosimm$organ==acceptorid,'probFlow']
  }
  
  write.csv(sc3_probmat, file='data/sc3_probmatrix.csv')
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI),1)
    rcol<-sample(1:ncol(dataBSI),1)
    # pFlow<- fisiosimm[fisiosimm$organ==rownames(dataBSI)[rrow],]$probFlow
    prob_mets = sc3_probmat[rrow, rcol]
    
    if(prob_mets>runif(1) ){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1}
    
  }
  nullBSI<- nullBSI/sum(nullBSI)
  
  gm<- computeModules(nullBSI)
  
  
  simnesttemp3<- nestedtemp(nullBSI)$statistic
  simnestnodf3<- nestednodf(nullBSI)$statistic[3]
  simnestwine3<-wine(nullBSI, nreps = 999)$wine
  simmod3<- gm@likelihood
  obstemp3<- nestedtemp(dataBSI)$statistic
  obsnodf3<- nestednodf(dataBSI)$statistic[3]
  obswine3<- wine(dataBSI,nreps = 999)$wine
  obsmod3<-computeModules(dataBSI)@likelihood
  fit<-lm(c(nullBSI)~c(dataBSI))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI3<- cor.test(nullBSI, dataBSI,method = 'spearman')$p.value
  corBSI3<-cor.test(nullBSI, dataBSI,method = 'spearman')$estimate
  pval03<-sfit$coefficients[2,4]
  write.csv(data.frame(run = nullnets, melt.matrix(nullBSI), melt.matrix(dataBSI)), paste0('data/raw_cor_trueplot/','sc3_run',nullnets,'.csv'))
  
  
  
  #Case B
  nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
  rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)
  
  sc4_probmat = sc3_probmat * Pmatrixb
  write.csv(sc4_probmat, file='data/sc4_probmatrix.csv')
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI),1)
    rcol<-sample(1:ncol(dataBSI),1)
    # pFlow<- fisiosimm[fisiosimm$organ==rownames(dataBSI3)[rrow],]$probFlow
    prob_mets = sc4_probmat[rrow, rcol]
    if(prob_mets>runif(1) ){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1} 
    
  }
  gm<- computeModules(nullBSI)
  
  nullBSI<- nullBSI/sum(nullBSI)
  simnesttemp4<- nestedtemp(nullBSI)$statistic
  simnestnodf4<- nestednodf(nullBSI)$statistic[3]
  simnestwine4<-wine(nullBSI,nreps = 999)$wine
  obstemp4<- nestedtemp(dataBSI)$statistic
  obsnodf4<- nestednodf(dataBSI)$statistic[3]
  obswine4<- wine(dataBSI,nreps = 999)$wine
  obsmod4<-computeModules(dataBSI)@likelihood
  simmod4<-  gm@likelihood
  
  fit<-lm(c(nullBSI)~c(dataBSI))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  corBSI4<-cor.test(nullBSI, dataBSI,method = 'spearman')$estimate
  pcorBSI4<-cor.test(nullBSI, dataBSI,method = 'spearman')$p.val
  pval04<-sfit$coefficients[2,4]
  write.csv(data.frame(run = nullnets, melt.matrix(nullBSI), melt.matrix(dataBSI)), paste0('data/raw_cor_trueplot/','sc4_run',nullnets,'.csv'))
  
  
  df0= data.frame(nullnets, 
                  simnesttemp3, simnesttemp4,
                  simnestnodf3, simnestnodf4,
                  simnestwine3, simnestwine4,
                  simmod3,simmod4,
                  corBSI3,corBSI4,
                  pcorBSI3, pcorBSI4,
                  pval03, pval04,
                  obstemp3, obstemp4,
                  obsnodf3, obsnodf4,
                  obswine3, obswine4,
                  obsmod3, obsmod4)
  
  
};beepr::beep(sound = 4)

stopCluster(cl)




#Scenarios V and VI: dP and blood flow corrected

fisiosimm<-physio_df %>%
  filter(P.g>0 & flow>0)%>%
  mutate(probFlow=flow/sum(flow))

occ<- df0[rownames(df0) %in% fisiosimm$organ, colnames(df0) %in% fisiosimm$organ]

sMI<- data.frame(Source=rownames(occ), sMI=rowSums(occ)/sum(occ))

overall<-sum(occ)

dataBSI<-occ/overall


Pmatrix<- dataBSI
for (i in 1:nrow(Pmatrix)) {
  for (j in 1:ncol(Pmatrix)) {
    nameSource<-rownames(Pmatrix)[i]
    nameAcceptor<- colnames(Pmatrix)[j] 
    Pdis<- (fisiosimm[fisiosimm$organ== nameSource,]$P.g - fisiosimm[fisiosimm$organ== nameAcceptor,]$P.g)
    Pmatrix[i,j]<-Pdis
  }
}

dataBSI<-as.matrix(dataBSI)
Pmatrix<-as.matrix(Pmatrix)

Pmatrixb<- Pmatrix
for(prow in 1:nrow(Pmatrixb)){  
  for(pcol in 1:ncol(Pmatrixb)){
    if(min(Pmatrixb[prow, ])<0){Pmatrixb[prow, pcol]<- 1- (Pmatrix[prow, pcol] + min(Pmatrix[prow, ])*-1)}
  }
}
rowsums<-rowSums(Pmatrixb)
for(prow in 1:nrow(Pmatrixb)){  
  for(pcol in 1:ncol(Pmatrixb)){
    Pmatrixb[prow, pcol]<- Pmatrixb[prow,pcol]/ rowsums[prow] 
  }
}

nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)

nulltopo_3<- data.frame()

cl <- makePSOCKcluster(n.cores)
registerDoParallel(cl)  # use multicore, set to the number of our cores

nulltopo_3<-  foreach (nullnets=1:sim, .combine = rbind) %dopar% {
  pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr, ggplot2, patchwork)
  print(nullnets)
  nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
  rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)
  
  sc5_probmat = nullBSI
  
  for(acceptorid in fisiosimm$organ){
    sc5_probmat[,acceptorid] = fisiosimm[fisiosimm$organ==acceptorid,'probFlow']
  }
  
  sc5_probmat[rownames(sc5_probmat) %in% c('Spleen')] = sc5_probmat[rownames(sc5_probmat) %in% c('Spleen')]/2
  
  
  write.csv(sc5_probmat, file='data/sc5_probmatrix.csv')
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI),1)
    rcol<-sample(1:ncol(dataBSI),1)
    # pFlow<- fisiosimm[fisiosimm$organ==rownames(dataBSI)[rrow],]$probFlow
    prob_mets = sc5_probmat[rrow, rcol]
    if(prob_mets>runif(1) ){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1} #& Pmatrix3b[rrow, rcol]>runif(1)
    
  }
  nullBSI<- nullBSI/sum(nullBSI)
  
  gm<- computeModules(nullBSI)
  
  
  simnesttemp5<- nestedtemp(nullBSI)$statistic
  simnestnodf5<- nestednodf(nullBSI)$statistic[3]
  simnestwine5<-wine(nullBSI, nreps = 999)$wine
  simmod5<- gm@likelihood
  obstemp5<- nestedtemp(dataBSI)$statistic
  obsnodf5<- nestednodf(dataBSI)$statistic[3]
  obswine5<- wine(dataBSI,nreps = 999)$wine
  obsmod5<-computeModules(dataBSI)@likelihood
  fit<-lm(c(nullBSI)~c(dataBSI))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  # pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  pcorBSI5<- cor.test(nullBSI, dataBSI,method = 'spearman')$p.value
  corBSI5<-cor.test(nullBSI, dataBSI,method = 'spearman')$estimate
  pval05<-sfit$coefficients[2,4]
  
  write.csv(data.frame(run = nullnets, melt.matrix(nullBSI), melt.matrix(dataBSI)), paste0('data/raw_cor_trueplot/','sc5_run',nullnets,'.csv'))
  
  
  #Case B
  nullBSI<- matrix(data = 0,ncol=ncol(Pmatrix), nrow=nrow(Pmatrix))
  rownames(nullBSI)<-rownames(dataBSI);colnames(nullBSI)<-colnames(dataBSI)
  
  sc6_probmat = sc5_probmat * Pmatrixb
  write.csv(sc6_probmat, file='data/sc6_probmatrix.csv')
  
  for(i in 1:mets){
    rrow<-sample(1:nrow(dataBSI),1)
    rcol<-sample(1:ncol(dataBSI),1)
    pFlow<- fisiosimm[fisiosimm$organ==rownames(dataBSI)[rrow],]$probFlow
    prob_mets = sc6_probmat[rrow, rcol]
    if(prob_mets>runif(1) ){nullBSI[rrow, rcol]<-nullBSI[rrow, rcol]+1} 
    
  }
  gm<- computeModules(nullBSI)
  
  nullBSI<- nullBSI/sum(nullBSI)
  simnesttemp6<- nestedtemp(nullBSI)$statistic
  simnestnodf6<- nestednodf(nullBSI)$statistic[3]
  simnestwine6<-wine(nullBSI,nreps = 999)$wine
  obstemp6<- nestedtemp(dataBSI)$statistic
  obsnodf6<- nestednodf(dataBSI)$statistic[3]
  obswine6<- wine(dataBSI,nreps = 999)$wine
  obsmod6<-computeModules(dataBSI)@likelihood
  simmod6<-  gm@likelihood
  
  fit<-lm(c(nullBSI)~c(dataBSI))
  sfit<- summary(fit)
  tstats <- (1-sfit$coefficients[2,1])/sfit$coefficients[2,2]
  # Calculates two tailed probability
  # pval<- 2 * pt(abs(tstats), df = df.residual(fit), lower.tail = FALSE)
  
  corBSI6<-cor.test(nullBSI, dataBSI,method = 'spearman')$estimate
  pcorBSI6<-cor.test(nullBSI, dataBSI,method = 'spearman')$p.val
  pval06<-sfit$coefficients[2,4]
  write.csv(data.frame(run = nullnets, melt.matrix(nullBSI), melt.matrix(dataBSI)), paste0('data/raw_cor_trueplot/','sc6_run',nullnets,'.csv'))
  
  
  df0= data.frame(nullnets, 
                  simnesttemp5, simnesttemp6,
                  simnestnodf5, simnestnodf6,
                  simnestwine5, simnestwine6,
                  simmod5,simmod6,
                  corBSI5,corBSI6,
                  pcorBSI5, pcorBSI6,
                  pval05, pval06,
                  obstemp5, obstemp6,
                  obsnodf5, obsnodf6,
                  obswine5, obswine6,
                  obsmod5, obsmod6)
  
  
};beepr::beep(sound = 4)

stopCluster(cl)
