rm(ALK.1);rm(BCOR.1);rm(CTNNB1.1);rm(TP53.1)
 ALK.1 <- mergeVCFObject.variantKey %>% dplyr::filter(Gene %in% c("ALK"))
 BCOR.1 <- mergeVCFObject.variantKey %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.1 <- mergeVCFObject.variantKey %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.1 <- mergeVCFObject.variantKey %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.1);dim(BCOR.1);dim(CTNNB1.1);dim(TP53.1)

 rm(ALK.2);rm(BCOR.2);rm(CTNNB1.2);rm(TP53.2)
 ALK.2 <- mergeVCFObjectDiagnosis %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.2 <- mergeVCFObjectDiagnosis %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.2 <- mergeVCFObjectDiagnosis %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.2 <- mergeVCFObjectDiagnosis %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.2);dim(BCOR.2);dim(CTNNB1.2);dim(TP53.2)
 
 rm(ALK);rm(BCOR);rm(CTNNB1);rm(TP53)
 ALK.4 <- mergeVCFObject.variantKey.Tumor %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.4 <- mergeVCFObject.variantKey.Tumor %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.4 <- mergeVCFObject.variantKey.Tumor %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.4 <- mergeVCFObject.variantKey.Tumor %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.4);dim(BCOR.4);dim(CTNNB1.4);dim(TP53.4)
 
 
 rm(ALK.5);rm(BCOR.5);rm(CTNNB1.5);rm(TP53.5)
 ALK.5 <- TumorFiltered.Normal %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.5 <- TumorFiltered.Normal %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.5 <- TumorFiltered.Normal %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.5 <- TumorFiltered.Normal %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.5);dim(BCOR.5);dim(CTNNB1.5);dim(TP53.5)
 
 
 rm(ALK.6);rm(BCOR.6);rm(CTNNB1.6);rm(TP53.6)
 ALK.6 <- TumorFiltered.Normal.freq %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.6 <- TumorFiltered.Normal.freq %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.6 <- TumorFiltered.Normal.freq %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.6 <- TumorFiltered.Normal.freq %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.6);dim(BCOR.6);dim(CTNNB1.6);dim(TP53.6)
 
 
 rm(ALK.7);rm(BCOR.7);rm(CTNNB1.7);rm(TP53.7)
 ALK.7 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.7 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.7 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.7 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.7);dim(BCOR.7);dim(CTNNB1.7);dim(TP53.7)
 
 rm(ALK.8);rm(BCOR.8);rm(CTNNB1.8);rm(TP53.8)
 ALK.8 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.8 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.8 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.8 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.8);dim(BCOR.8);dim(CTNNB1.8);dim(TP53.8)
 
 rm(ALK.9);rm(BCOR.9);rm(CTNNB1.9);rm(TP53.9)
 ALK.9 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 %>%  dplyr::filter(Gene %in% c("ALK"))
 BCOR.9 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 %>%  dplyr::filter(Gene %in% c("BCOR"))
 CTNNB1.9 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 %>%  dplyr::filter(Gene %in% c("CTNNB1"))
 TP53.9 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 %>%  dplyr::filter(Gene %in% c("TP53"))
 dim(ALK.9);dim(BCOR.9);dim(CTNNB1.9);dim(TP53.9)
 
 
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 1.0")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 1.1")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 1.2")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 1.3")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 1.4")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 2")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 3")]$variantKey))
  length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Variant.Tier %in% c("Tier 4")]$variantKey))
 
 
 
 
 
  sampleList.v15 <- read.csv("../RNASeq.Mutation.data/SampleList_Final_V15_makePaths_caseIDCorrected.txt", sep="\t", header = T)
  oldMetaDAta <- rnaseqMutationProject$validMetaDataDF
  
  removedSamples <- oldMetaDAta[which(! oldMetaDAta$Sample.Data.ID %in% sampleList.v15$SampleID),]
 
 
 
 old.V1 <- readRDS("../RNASeq.Mutation.data/outputRDSOutput/1.All.variants.RDS")
 New.V1 <- readRDS("../RNASeq.Mutation.data/outputRDSOutput/1.All.variants.v2.RDS")
 variantsAbsentInNew <- old.V1[ which(!old.V1$variantKey %in% New.V1$variantKey),] %>% dplyr::filter(Variant.Tier %in% c("Tier 1.1")); dim(variantsAbsentInNew)
 
 
 metaData <- read.csv("../RNASeq.Mutation.data/SampleList_Final_V16.txt", sep="\t", header = T)
 customFun  = function(DF) {
   printDF <- DF[,c("Patient.ID.In.paper","CaseID")]
   printDF$rnaseq <- "rnaseq" ; printDF <- printDF %>% distinct()
   write.table(printDF,paste0("../RNASeq.Mutation.data/",unique(DF$Patient.ID.In.paper),".txt"), sep = "\t", col.names = F, row.names = F, quote = F)
   return(printDF)
 }
 
 metaData %>% 
   group_by(Patient.ID.In.paper) %>% 
   do(customFun(.))
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 