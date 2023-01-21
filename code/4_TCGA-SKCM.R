if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("GDCRNATools")

if (!require("pacman")) install.packages("pacman")
pacman::p_load(GDCRNATools, dplyr, tidyr, ggplot2, ggrepel)

system("wget https://gdc.cancer.gov/files/public/file/gdc-client_v1.6.1_Ubuntu_x64.zip")
unzip("gdc-client_v1.6.1_Ubuntu_x64.zip", exdir=".")


gdcGetURL_new <- function(project.id, data.type) {
  urlAPI <- 'https://api.gdc.cancer.gov/files?'
  
  if (data.type=='RNAseq') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Gene Expression Quantification'
    workflow.type <- 'STAR - Counts'
  } else if (data.type=='miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'Isoform Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  } else if (data.type=='Clinical') {
    data.category <- 'Clinical'
    data.type <- 'Clinical Supplement'
    workflow.type <- NA
  } else if (data.type=='pre-miRNAs') {
    data.category <- 'Transcriptome Profiling'
    data.type <- 'miRNA Expression Quantification'
    workflow.type <- 'BCGSC miRNA Profiling'
  }
  
  project <- paste('{"op":"in","content":{"field":"cases.',
                   'project.project_id","value":["', 
                   project.id, '"]}}', sep='')
  dataCategory <- paste('{"op":"in","content":{"field":"files.', 
                        'data_category","value":"', data.category, '"}}', sep='')
  dataType <- paste('{"op":"in","content":{"field":"files.data_type",',
                    '"value":"', data.type, '"}}', sep='')
  workflowType <- paste('{"op":"in","content":{"field":"files.',
                        'analysis.workflow_type","value":"', workflow.type, '"}}', sep='')
  
  
  if (is.na(workflow.type)) {
    dataFormat <- paste('{"op":"in","content":{"field":"files.',
                        'data_format","value":"', 'BCR XML', '"}}', sep='')
    content <- paste(project, dataCategory, dataType, dataFormat, sep=',')
  } else {
    content <- paste(project, dataCategory, dataType, 
                     workflowType, sep=',')
  }
  
  filters <- paste('filters=',URLencode(paste('{"op":"and","content":[', 
                                              content, ']}', sep='')),sep='')
  
  expand <- paste('analysis', 'analysis.input_files', 'associated_entities',
                  'cases', 'cases.diagnoses','cases.diagnoses.treatments', 
                  'cases.demographic', 'cases.project', 'cases.samples', 
                  'cases.samples.portions', 'cases.samples.portions.analytes', 
                  'cases.samples.portions.analytes.aliquots',
                  'cases.samples.portions.slides', sep=',')
  
  expand <- paste('expand=', expand, sep='')
  
  payload <- paste(filters, 'pretty=true', 'format=JSON', 
                   'size=10000', expand, sep='&')
  url <- paste(urlAPI, payload, sep='')
  
  return (url)
}

toolenv <- environment(get("gdcGetURL", envir = asNamespace("GDCRNATools")))
unlockBinding("gdcGetURL", toolenv)
assignInNamespace("gdcGetURL", gdcGetURL_new, ns="GDCRNATools", envir=toolenv)
assign("gdcGetURL", gdcGetURL_new)
lockBinding("gdcGetURL", toolenv)

# downloaded files will be stored under TCGA-SKCM_1611/RNAseq folder
project <- 'TCGA-SKCM'
rnadir <- paste(project, 'RNAseq', sep='/')

# Download RNAseq data 
gdcRNADownload(project.id     = 'TCGA-SKCM', 
               data.type      = 'RNAseq', 
               write.manifest = FALSE,
               method         = 'gdc-client',
               directory      = rnadir)


# Query metadata from GDC graph associated with RNA-seq quantification file
metaMatrix.RNA <- gdcParseMetadata(project.id = 'TCGA-SKCM',
                                   data.type  = 'RNAseq', 
                                   write.meta = FALSE)
# Filter duplicates
metaMatrix.RNA <- gdcFilterDuplicate(metaMatrix.RNA)

# Melt RNAseq files to a df
merge_rna <-function(metadata, fdir){
  filelist <- list.files(fdir, pattern="*.tsv$", 
                         recursive = TRUE, full.names=TRUE)
  for (i in 1:length(filelist)){
    iname <- basename(filelist[i])
    isamplename <- metadata[metadata$file_name==iname, "sample"]
    idf <- read.csv(filelist[i], sep="\t", skip=1, header=TRUE)
    # remove first 4 rows
    remove <- 1:4
    idf_subset <- idf[-remove, c("gene_id","unstranded")]
    rm(idf)
    names(idf_subset)[2] <- isamplename
    #print(dim(idf_subset))
    if (i==1){
      combined_df <- idf_subset
      rm(idf_subset)
    } else {
      combined_df <- merge(combined_df, idf_subset, by.x='gene_id', by.y="gene_id", all=TRUE)
      rm(idf_subset)
    }
  }
  # remove certain gene ids
  combined_df <- combined_df[!(grepl("PAR_Y", combined_df$gene_id, fixed=TRUE)),]
  # modify gene_id
  combined_df$gene_id <- sapply(strsplit(combined_df$gene_id,"\\."), `[`, 1)
  # use gene_id as row names and remove gene_id column
  rownames(combined_df) <- combined_df$gene_id
  combined_df <- combined_df[,-which(names(combined_df) %in% c("gene_id"))]
  return(combined_df)
}

rnaCounts <-  merge_rna(metaMatrix.RNA, rnadir)

# Filtering for primary melanomas and metastases
removeids <- metaMatrix.RNA[metaMatrix.RNA$sample_type %in% c('AdditionalMetastatic', 'SolidTissueNormal'),]$sample
rnaCounts_2 <- rnaCounts[, !colnames(rnaCounts) %in% removeids]
metaMatrix.RNA_2 <- metaMatrix.RNA[!metaMatrix.RNA$sample_type %in% c('AdditionalMetastatic', 'SolidTissueNormal'),]

clinical_df = read.csv("C:/Users/scastillo/Dropbox (ICR)/MetastaticNetwork/0_Science/data/tcga_clinical.csv")

removeids_2 = clinical_df[clinical_df$site_2 == 'remove',]$case_submitter_id
removeids_2 = paste(removeids_2, c('06', '01'),sep = '-')
removeids_2A = paste0(removeids_2, c('A'))

rnaCounts_3 <- rnaCounts_2[, !colnames(rnaCounts_2) %in% removeids_2]
metaMatrix.RNA_3 <- metaMatrix.RNA_2[!metaMatrix.RNA_2$submitter_id %in% removeids_2A,]
remove_temp = metaMatrix.RNA_3$sample[!metaMatrix.RNA_3$sample %in% colnames(rnaCounts_3)]
metaMatrix.RNA_3 <- metaMatrix.RNA_3[!metaMatrix.RNA_3$sample %in% remove_temp,]


# Normalization of RNAseq data 
rnaExpr <- gdcVoomNormalization(counts = rnaCounts_3, filter = FALSE)

# Differential expression analysis Metastases relative to primary tumour
DE_All_SKCM<- gdcDEAnalysis(counts = rnaCounts_3, 
                            group      = metaMatrix.RNA_3$sample_type, 
                            comparison = 'Metastatic-PrimaryTumor',
                            method     = 'DESeq2')


# All genes passed filtering (FC >= 2, pval =< 0.01) including protein-coding and non-coding
de_ALL_SKCM <- gdcDEReport(deg = DE_All_SKCM, gene.type = 'all', fc = 2, pval = 0.01)


#Univateriate survival analysis for each gene expression, cut-off median
survOutput <- gdcSurvivalAnalysis(gene     = rownames(de_ALL_SKCM), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA_3, 
                                  sep      = 'median')

# Sort the output of univariate survival analysis based on pvalue
survOutput$pValue <- as.numeric(survOutput$pValue)
sorted_survOutput<- survOutput[order(survOutput$pValue),]
sorted_survOutput[1:5,]


#Univateriate survival analysis for each gene expression, Q1 cut-off
survOutput <- gdcSurvivalAnalysis(gene     = rownames(de_ALL_SKCM), 
                                  method   = 'KM', 
                                  rna.expr = rnaExpr, 
                                  metadata = metaMatrix.RNA_3, 
                                  sep      = '1stQu')

# Sort the output of univariate survival analysis based on pvalue
survOutput$pValue <- as.numeric(survOutput$pValue)
sorted_survOutput<- survOutput[order(survOutput$pValue),]
sorted_survOutput[1:5,]





# pick the most significant gene for KM visualization
gdcKMPlot(gene     = rownames(sorted_survOutput[1,]),
          rna.expr = rnaExpr,
          metadata = metaMatrix.RNA_3,
          sep      = 'median')

