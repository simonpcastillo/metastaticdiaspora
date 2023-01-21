if (!require("pacman")) install.packages("pacman")
pacman::p_load(bipartite,reshape,foreach,doParallel, dplyr)


data<-read.csv("data/occurences.csv", sep=",", header=T, row.names = 1)
    
total<-sum(data)
dataBSI<-data/total

# Degree distirbution
degreedistr(dataBSI)

# Observed nestedness
BSI_nest<-as.data.frame(bipartite::nested(dataBSI, method = c("binmatnest", 'NODF2', 'wine')))    
BSI_nest<-as.data.frame((BSI_nest)); names(BSI_nest)<-'value'

temp.obs<-nestedtemp(dataBSI)$statistic
wine.obs<- wine(dataBSI)$wine
nodf.obs<-nestednodf(dataBSI, weighted = TRUE)$statistic[3]


# Observed modularity
    mod.obs<-computeModules(as.matrix(dataBSI))@likelihood
    
    # Observed modularity
    mod.obs<-computeModules(as.matrix(dataBSI))@likelihood
    
    # Null models modularity
    rep.nulls= 1000
    n.cores = 15 # use multicore, set to the number of our cores
    
    # Null model 1
    method.nulls<- 'r0_ind'  
    nm<-vegan::nullmodel(data,method.nulls)  
    null<-simulate(nm, nsim =rep.nulls)
    
    nulls_1<- data.frame()
    registerDoParallel(n.cores)  
    
    nulls_1 <- foreach (p=1:rep.nulls, .combine = rbind) %dopar% {
      
      pacman::p_load(foreach,doParallel, vegan, bipartite)   
      nullnet= null[,,p]/sum(null[,,p])
      temp.null<-nestedtemp(nullnet)$statistic
      wine.null<- wine(nullnet)$wine
      nodf.null<-nestednodf(nullnet, weighted = TRUE, order = TRUE)$statistic[3]
      md<-computeModules(nullnet) 
      mod.null<-md@likelihood
      
      null0=data.frame(temp= temp.null, wine= wine.null, nodf=nodf.null,mod=mod.null)
      nulls_1<- rbind(nulls_1, null0)
    }  
    
    # Null model 2
    method.nulls<- 'c0_ind'  
    nm<-vegan::nullmodel(data,method.nulls)  
    null<-simulate(nm, nsim =rep.nulls)
    
    nulls_2<- data.frame()
    registerDoParallel(n.cores)  # use multicore, set to the number of our cores
    
    nulls_2<-  foreach (p=1:rep.nulls, .combine = rbind) %dopar% {
      
      pacman::p_load(foreach,doParallel, vegan, bipartite)   
      nullnet= null[,,p]/sum(null[,,p])
      temp.null<-nestedtemp(nullnet)$statistic
      wine.null<- wine(nullnet)$wine
      nodf.null<-nestednodf(nullnet, weighted = TRUE)$statistic[3]
      md<-computeModules(nullnet) 
      mod.null<-md@likelihood
      
      null0=data.frame(temp= temp.null, wine= wine.null, nodf=nodf.null,mod=mod.null)
      nulls_2<- rbind(nulls_2, null0)
    }  
    
    
    
    summ_null1 = nulls_1 %>%
      summarise_all(list(mean, var)) %>%
      mutate(
        z.temp = (temp.obs- temp_fn1)/sqrt(temp_fn2),
        z.wine =(wine.obs-wine_fn1)/sqrt(wine_fn2), 
        z.nodf=(nodf.obs-nodf_fn1)/sqrt(nodf_fn2),
        z.mod =(mod.obs- mod_fn1)/sqrt(mod_fn2)
      )
    
    summ_null2 = nulls_2 %>%
      summarise_all(list(mean, var)) %>%
      mutate(
        z.temp = (temp.obs- temp_fn1)/sqrt(temp_fn2),
        z.wine =(wine.obs-wine_fn1)/sqrt(wine_fn2), 
        z.nodf=(nodf.obs-nodf_fn1)/sqrt(nodf_fn2),
        z.mod =(mod.obs- mod_fn1)/sqrt(mod_fn2)
      )
    
    
    