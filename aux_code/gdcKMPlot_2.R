##' @title Kaplan Meier plot
##' @description Plot Kaplan Meier survival curve
##' @param gene an Ensembl gene id
##' @param rna.expr \code{\link[limma]{voom}} transformed expression data
##' @param metadata metadata parsed from \code{\link{gdcParseMetadata}}
##' @param sep a character string specifying which point should be used to 
##'     separate low-expression and high-expression groups. Possible values 
##'     are \code{'1stQu'}, \code{'mean'}, \code{'median'}, 
##'     and \code{'3rdQu'}. Default is \code{'median'} 
##' @importFrom survival survfit
##' @importFrom survival survdiff
##' @importFrom survminer ggsurvplot
##' @return A plot of Kaplan Meier survival curve
##' @export
##' @author Ruidong Li and Han Qu
##' @examples 
##' ####### KM plots #######
##' 
if (!require("pacman")) install.packages("pacman")
pacman::p_load(survival, survminer)

gdcKMPlot_2 <- function(gene, gene_id,  rna.expr, metadata, sep='median') {
  
  metadata <- metadata[metadata$sample_type=='PrimaryTumor',]
  
  samples = intersect(colnames(rna.expr), metadata$sample)
  
  exprDa=rna.expr[gene,samples]
  
  
  if (sep=='1stQu') {
    thresh <- as.numeric(summary(exprDa)[2])
  } else if (sep=='median') {
    thresh <- as.numeric(summary(exprDa)[3])
  } else if (sep=='mean') {
    thresh <- as.numeric(summary(exprDa)[4])
  } else if (sep=='3rdQu') {
    thresh <- as.numeric(summary(exprDa)[5])
  }
  
  exprGroup <- exprDa > thresh
  
  nH <- sum(exprGroup)
  nL <- sum(!exprGroup)
  
  clinicalDa=metadata[match(samples,metadata$sample),]
  
  daysToDeath <- as.numeric(clinicalDa$days_to_death)
  nonComplt <- is.na(daysToDeath)
  
  vitalStatus <- as.numeric(ifelse(nonComplt, 0, 1))
  daysToDeath[nonComplt] <- 
    as.numeric(clinicalDa$days_to_last_follow_up[nonComplt])
  
  survDa <- data.frame(daysToDeath,vitalStatus, exprGroup)
  
  sdf <- survdiff(Surv(daysToDeath, vitalStatus) ~ exprGroup)
  pValue <- format(pchisq(sdf$chisq, length(sdf$n)-1, 
                          lower.tail = FALSE),digits=3)
  #pValue <- format(1-pchisq(sdf$chisq, df=1),digits=3)
  
  HR = (sdf$obs[2]/sdf$exp[2])/(sdf$obs[1]/sdf$exp[1])
  upper95 = exp(log(HR) + qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  lower95 = exp(log(HR) - qnorm(0.975)*sqrt(1/sdf$exp[2]+1/sdf$exp[1]))
  
  HR <- format(HR, digits = 3)
  upper95 <- format(upper95, digits = 3)
  lower95 <- format(lower95, digits = 3)
  
  
  label1 <- paste('HR = ', HR, ' (', lower95, '-', upper95, ')', sep='')
  label2 <- paste('logrank P value = ', pValue, sep='')
  
  fit <- survfit(Surv(daysToDeath, vitalStatus) ~ exprGroup, data=survDa)
  
  lgdXpos <- 1/1.3
  lgdYpos = 0.25
  
  xpos = max(daysToDeath, na.rm=TRUE)/1.8
  ypos1 = 0.8
  
  ypos2 = 0.75
  
  ggsurvplot(fit, data=survDa, 
             pval = paste(label1, '\n', label2), 
             pval.coord = c(xpos, ypos1),
             pval.size=5,
             font.main = c(16, 'bold', 'black'), conf.int = FALSE, 
             title = gene_id,
             legend = c(lgdXpos, lgdYpos), 
             #legend = 'none',
             #color = c('blue', 'green'),
             palette= c('green', 'purple'),
             legend.labs = c(paste('lowExp (N = ',nL,')',sep=''),paste('HighExp (N = ',nH,')',sep='')),  
             legend.title='',
             xlab = paste('Overall survival (days)'), ylab = 'Survival probability',
             font.x = c(16), font.y = c(16), ylim=c(0,1), #16
             ggtheme = theme_minimal()) 
  
}