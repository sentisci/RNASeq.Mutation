rm(list=ls())
## Source all classes and packages

source("./utilityPackages.R")
source("./statisticalPackages.R")
source("./class.R")

## Project Title: Expression Analysis for Landscape paper

## ## Instantiate a new Object of type ProjectSetUp
rnaseqMutationProject <- ProjectSetUp$new(
  
  date                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[1]],
  time                    = unlist(strsplit(x = as.character(Sys.time()), "\\s+"))[[2]],
  projectName             = "RNASeq.Mutation.data",
  annotationRDS           = "C:/Users/sindiris/R Scribble/Annotation RDS/annotation_ENSEMBL_gene.RDS",
  pcRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/pc.other.HGNCTableFlat.rds",
  tfRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/TFs_no_epimachines.RDS",
  csRDS                   = "C:/Users/sindiris/R Scribble/Annotation RDS/CellSurface.RDS",
  cgaRDS                  = "C:/Users/sindiris/R Scribble/Annotation RDS/cancerGermlineAntigens.rds",
  ewsr1Fli1RDS            = "C:/Users/sindiris/R Scribble/Annotation RDS/EWSR1_FL1_DownstreamTargets.RDS",
  pax3Foxo1RDS             = "C:/Users/sindiris/R Scribble/Annotation RDS/PAX3_FOXO1_DownstreamTargets.RDS",

  outputPrefix            = "_rnaseq.annotated.vcf",
  factorName              = "PATIENT.ID",
  metaDataFileName        = "MetadataMapper.txt",
  outputdirRDSDir         = "outputRDSOutput",
  outputdirTXTDir         = "outputTXTOutput",
  plotsDir                = "Figures",
  plotsDataDir            = "FigureData",
  vcfFilesDir             = "vcfFiles"
  #factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), 
  #                               "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "PolyA"),
  #                               "Tumors"=list("LIBRARY_TYPE"="Tumor", "LibraryPrep" = "PolyA"))
  #factorsToExclude        = list("CellLine"=list("LIBRARY_TYPE"="CellLine"), "Normal.ribozero"=list("LIBRARY_TYPE"="Normal", "LibraryPrep" = "Ribozero"))                                     
)

## Add utility functions to the project
corUtilsFuncs <- CoreUtilities$new(  ProjectSetUpObject = rnaseqMutationProject )

## Generate expression matrix
rm(mergeObjectsNoDup)
mergeVCFObject <- corUtilsFuncs$getMergedMatrix(dir               = "vcfFiles", 
                                                   fileFormat        = "txt", 
                                                   fileSuffix        = "_rnaseq.annotated.vcf",
                                                   metadata          = rnaseqMutationProject$metaDataDF,
                                                   metadataFileRefCol= "Patient.ID")
## Make syntactically right names.
mergeVCFObject$Patient.ID <- make.names(mergeVCFObject$Patient.ID)
## Append Diagnosis
## Not Working
# mergeVCFObjectDiagnosis <- corUtilsFuncs$leftjoinDTs(mergeVCFObject, data.table(rnaseqMutationProject$validMetaDataDF[, c("Patient.ID","DIAGNOSIS.Alias")]) , 
#                                                     key = "Patient.ID") ; class(mergeVCFObjectDiagnosis); dim(mergeVCFObjectDiagnosis)
## STEP 1 (Assemble all the variants in the cohort)
mergeVCFObjectDiagnosis <- dplyr::left_join(mergeVCFObject, data.table(rnaseqMutationProject$validMetaDataDF[, c("Patient.ID","DIAGNOSIS.Alias")]), by="Patient.ID") %>% 
                           dplyr::filter(! Patient.ID %in% c("NA.")) %>% 
                           dplyr::filter(!is.na(MAF)) %>% 
                           data.table()
dim(mergeVCFObjectDiagnosis)
## dplyr way
#mergeVCFObjectDiagnosis <- mergeVCFObjectDiagnosis %>% dplyr::mutate(variantKey = paste(Chr,Start,End,Ref,Alt))
## data.table way
mergeVCFObject.variantKey <- mergeVCFObjectDiagnosis[, variantKey := paste(Chr,Start,End,Ref,Alt,sep="_")][order(variantKey)]
dim(mergeVCFObject.variantKey)

## Separate Normal variants into different DF
## STEP 2 (separate variants presents in Normal and Tumor)
## dplyr way
# mergeVCFObjectDiagnosis.Normal <- mergeVCFObjectDiagnosis %>% dplyr::filter(DIAGNOSIS.Alias %in% c("NS"))
# dim(mergeVCFObjectDiagnosis.Normal)
## data.table way
mergeVCFObject.variantKey.Normal <- mergeVCFObject.variantKey[ DIAGNOSIS.Alias %in% c("NS") ]
dim(mergeVCFObject.variantKey.Normal)

## dplyr way
# mergeVCFObjectDiagnosis.Tumor <- mergeVCFObjectDiagnosis %>% dplyr::filter(!DIAGNOSIS.Alias %in% c("NS"))
# dim(mergeVCFObjectDiagnosis.Tumor)
## data.table way
mergeVCFObject.variantKey.Tumor <- mergeVCFObject.variantKey[ !DIAGNOSIS.Alias %in% c("NS") ]
dim(mergeVCFObject.variantKey.Tumor)


## STEP 3 (filter the Normal variants from Tumor )
## dplyr way
# TumorFiltered.Normal <- mergeVCFObjectDiagnosis.Tumor %>% dplyr::filter(!variantKey %in% mergeVCFObjectDiagnosis.Normal$variantKey)
# dim(TumorFiltered.Normal)
## data.table way
TumorFiltered.Normal <- mergeVCFObject.variantKey.Tumor[ !variantKey %in% mergeVCFObject.variantKey.Normal$variantKey ]
setkey(TumorFiltered.Normal, variantKey)
dim(TumorFiltered.Normal)

## STEP 4 Grouping the variants and calculating their frequencies
## 1st method
variantsFrequency = TumorFiltered.Normal[, .N , by = variantKey][,prop := (N/sum(N))*100]
setkey(variantsFrequency, variantKey)
variantsFrequency <- variantsFrequency[order(prop),]
dim(variantsFrequency)

## STEP 5 Join the frequency table with the main table
TumorFiltered.Normal.freq <- merge( TumorFiltered.Normal, variantsFrequency, by="variantKey" )
TumorFiltered.Normal.freq <- TumorFiltered.Normal.freq[order(prop),]
##  make columns numeric
TumorFiltered.Normal.freq[, (c("MAF", "VAF")) := replace(.SD, is.na(.SD), 0), .SDcols = c("MAF", "VAF")]
## Testing the "fkt" columns
min(TumorFiltered.Normal.freq$Total.coverage); min(TumorFiltered.Normal.freq$Variant.coverage); 
min(TumorFiltered.Normal.freq$VAF)           ; min(TumorFiltered.Normal.freq$MAF)


## STEP 6 filter by VAF, Variant coverage
TumorFiltered.Normal.freq.VAF.TC.VC.MAF <- TumorFiltered.Normal.freq[
                                              Total.coverage >= 10 & Variant.coverage >= 3 & VAF >= 0.10 & MAF <= 0.01 & prop <= 0.01 ]
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF)

## STEP 7 filter by population filter 
TumorFiltered.GermlineTierTier <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Germline.Tier][order(Germline.Tier)]
TumorFiltered.GermlineN <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = N][order(N)]
TumorFiltered.GermlineClinvar <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Clinvar.Pathogenic][order(Clinvar.Pathogenic)]
TumorFiltered.GermlineIntervar <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Intervar][order(Intervar)]
TumorFiltered.GermlineHGMD <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = HGMD.Disease][order(HGMD.Disease)]

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Germline.Tier %in% c("Tier 1.0",
                                                                                                               "Tier 1.1",
                                                                                                               "Tier 1.2",
                                                                                                               "Tier 1.3",
                                                                                                               "Tier 1.4")]
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1[
                                      Clinvar.Pathogenic %in% c("Y") | HGMD.Disease %in% c("Y") | Intervar %in% c("Likely pathogenic",
                                                                                                                  "Pathogenic",
                                                                                                                  "Uncertain significance") ]
sort(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing$Gene))


TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Germline.Tier %in% c("Tier 2")]
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Germline.Tier %in% c("Tier 3")]
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Germline.Tier %in% c("Tier 4")]


