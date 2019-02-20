rm(list=ls())
## Source all classes and packages

source("./utilityPackages.R")
source("./statisticalPackages.R")
source("./class.R")

## Project Title: Mutation Analysis for Landscape paper

## Make files to do the API calling

# metaData <- read.csv("../RNASeq.Mutation.data/SampleList_Final_V16.txt", sep="\t", header = T)
# customFun  = function(DF) {
#   printDF <- DF[,c("Patient.ID.In.paper","CaseID")]
#   printDF$rnaseq <- "rnaseq"
#   write.table(printDF,paste0("../RNASeq.Mutation.data/",unique(DF$Patient.ID.In.paper),".txt"), sep = "\t", col.names = F, row.names = F, quote = F)
#   return(printDF)
# }
# 
# metaData[1:5,] %>% 
#   group_by(Patient.ID.In.paper) %>% 
#   do(customFun(.))

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
  pax3Foxo1RDS            = "C:/Users/sindiris/R Scribble/Annotation RDS/PAX3_FOXO1_DownstreamTargets.RDS",

  outputPrefix            = "_rnaseq.annotated.vcf",
  factorName              = "PATIENT.ID",
  metaDataFileName        = "MetadataMapper.v3.txt",
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

## Generate variant matrix
rm(mergeObjectsNoDup)
mergeVCFObject <- corUtilsFuncs$getMergedMatrix(dir               = "vcfFiles", 
                                                   fileFormat        = "txt", 
                                                   fileSuffix        = "_rnaseq.annotated.vcf",
                                                   metadata          = rnaseqMutationProject$metaDataDF,
                                                   metadataFileRefCol= "Patient.ID")

## Make syntactically right names.
mergeVCFObject$Patient.ID <- make.names(mergeVCFObject$Patient.ID)
mergeVCFObject.variantKey <- mergeVCFObject[, variantKey := paste(Chr,Start,End,Ref,Alt,sep="_")][order(variantKey)]
print(paste(" Total Variants in the cohort: ", length(mergeVCFObject.variantKey$variantKey), 
            "& Total unique Variants in the cohort: ", length(unique(mergeVCFObject.variantKey$variantKey))))

## STEP 0 Make a tier DF
TierDF <- data.table(Tier=c("Tier 1.0","Tier 1.1","Tier 1.2","Tier 1.3","Tier 1.4","Tier 2","Tier 3","Tier 4", "NA"), code=c(1,2,3,4,5,6,7,8,9))
## Making HasTable environment
hash_table <- new.env(hash = TRUE, parent = emptyenv())
for(i in seq(nrow(TierDF))){
  hash_table[[ TierDF[[i,1]] ]] <- TierDF[[i,2]]
}

## STEP 1 (Assemble all the variants in the cohort and do basic data processing)
mergeVCFObjectDiagnosis <- dplyr::left_join(mergeVCFObject.variantKey, 
                                unique( data.table(rnaseqMutationProject$validMetaDataDF[, c("DIAGNOSIS.Alias","LIBRARY_TYPE","Sample.ID")]) ), 
                                by="Sample.ID") %>% 
                                data.table()

## Sanity Check : To evaluate the join // Check if join produced any NA
df <- mergeVCFObjectDiagnosis %>% data.frame()
df$Any_NA <- apply(df[,grep(paste(c("Patient.ID","DIAGNOSIS.Alias","LIBRARY_TYPE","Sample.ID"),collapse = "|"), names(df) )],1, function(x) anyNA(x))
df_NA_TRUE <- df[which(df$Any_NA == "TRUE"),] %>% dplyr::distinct(Patient.ID,Sample.ID);dim(df_NA_TRUE); View(df_NA_TRUE)

## Print
paste0( "Total Variants in the cohort ", dim(mergeVCFObjectDiagnosis)[1], 
        "; Unique Variants in the cohort ",length(unique(mergeVCFObjectDiagnosis$variantKey)) )

## Remove samples excluded from this project
mergeVCFObjectDiagnosis <- mergeVCFObjectDiagnosis %>% dplyr::filter(!Sample.ID %in% df_NA_TRUE$Sample.ID) %>% data.table()

## Sanity Check : match the samples & patient in the vcf files with the metadata
patientsEqual = length(unique(mergeVCFObjectDiagnosis$Patient.ID)) == length(unique(rrnaseqMutationProject$validMetaDataDF$Patient.ID$Patient.ID))
samplesEqual  = length(unique(mergeVCFObjectDiagnosis$Sample.ID)) == length(unique(rrnaseqMutationProject$validMetaDataDF$Patient.ID$Sample.ID))

if(! patientsEqual | ! samplesEqual ) { stop(paste0("Mismatch of Patients and/or Samples between VCF files and Metadata file !! ",
                                                       ", Total patients from VCF: ", length(unique(mergeVCFObjectDiagnosis$Patient.ID)),
                                                       ", Total patients from Metadata: ", length(unique(rrnaseqMutationProject$validMetaDataDF$Patient.ID$Patient.ID)),
                                                       ", Total samples from VCF: ", length(unique(mergeVCFObjectDiagnosis$Sample.ID)),
                                                       ", Total samples from Metadata: ", length(unique(mergeVCFObjectDiagnosis$Sample.ID)))) }
totalPatients <- length(unique(mergeVCFObjectDiagnosis$Patient.ID))
totalSamples  <- length(unique(mergeVCFObjectDiagnosis$Sample.ID))

## Print
paste0( "Total Variants in the cohort after removing extra samples ", dim(mergeVCFObjectDiagnosis)[1], 
        "; Unique Variants in the cohort after removing extra samples  ",length(unique(mergeVCFObjectDiagnosis$variantKey)) )

## 1a
mergeVCFObjectDiagnosis[, (c("Exonic.function")) := ifelse(Region %in% c("splicing","ncRNA_exonic","ncRNA_splicing","intronic",
                                                                         "UTR3","ncRNA_intronic","-","intergenic","upstream"),
                                                                                      Region , Exonic.function)]
unique(mergeVCFObjectDiagnosis$Exonic.function); unique(mergeVCFObjectDiagnosis$Region); dim(mergeVCFObjectDiagnosis)
## 1b
mergeVCFObjectDiagnosis[, (c("MAF", "VAF")) := replace(.SD, is.na(.SD), 0), .SDcols = c("MAF", "VAF")]; dim(mergeVCFObjectDiagnosis)

## 1c
## SNVs and Indels
mergeVCFObjectDiagnosisVA <- mergeVCFObjectDiagnosis %>% dplyr::filter(! Patient.ID %in% c("NA.") ) %>% 
                                                       dplyr::filter(! Exonic.function %in% c("ncRNA_exonic","ncRNA_splicing",
                                                                                     "ncRNA_intronic","-","intergenic",
                                                                                     "synonymous SNV")) %>% data.table()
## Cosmic 
mergeVCFObjectDiagnosisCosmic <- mergeVCFObjectDiagnosis %>% dplyr::filter(! Patient.ID %in% c("NA.") ) %>% 
                                                       dplyr::filter(! Exonic.function %in% c("ncRNA_exonic","ncRNA_splicing",
                                                                                     "ncRNA_intronic","-","intergenic")) %>% data.table()
## Printing
print(paste(" Total Variants in the cohort after removing non-coding,splicing aand intronic: ", length(mergeVCFObjectDiagnosisVA$variantKey), 
            " & Total Variants in the cohort after removing non-coding,splicing aand intronic: ", length(unique(mergeVCFObjectDiagnosisVA$variantKey))))

## Print & Save
saveRDS(mergeVCFObjectDiagnosisVA,
        paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, "1.All.variants.cosmic.v3.RDS", sep="/"))
write.table(mergeVCFObjectDiagnosisVA, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "1.All.variants.cosmic.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## 1d
## Separate Normal variants into different DF
mergeVCFObject.variantKey.Normal <- mergeVCFObjectDiagnosisVA[ LIBRARY_TYPE %in% c("Normal") ]
dim(mergeVCFObject.variantKey.Normal); length(unique(mergeVCFObject.variantKey.Normal$variantKey))

## Separate Tumor variants into different DF
mergeVCFObject.variantKey.Tumor <- mergeVCFObjectDiagnosisVA[ !LIBRARY_TYPE %in% c("Normal") ]
dim(mergeVCFObject.variantKey.Tumor); length(unique(mergeVCFObject.variantKey.Tumor$variantKey))

## filter the Normal variants from Tumor DF
TumorFiltered.Normal <- mergeVCFObject.variantKey.Tumor[ !variantKey %in% mergeVCFObject.variantKey.Normal$variantKey ];setkey(TumorFiltered.Normal, variantKey)
dim(TumorFiltered.Normal); length(unique(TumorFiltered.Normal$variantKey))

## 1e
## Updating the empty values in Tier column with "NA"
TumorFiltered.Normal$Somatic.Tier[TumorFiltered.Normal$Somatic.Tier==""] <- "NA"
TumorFiltered.Normal$Germline.Tier[TumorFiltered.Normal$Germline.Tier ==""] <- "NA"

## Printing
print(paste(" Total Variants in the cohort after removing Normal variants: ", length(TumorFiltered.Normal$variantKey), 
            " & Total Variants in the cohort after removing Normal variants : ", length(unique(TumorFiltered.Normal$variantKey))))

## Step 2 Add the smallest Tier
## Using R environment ( fastest )
Variant.Tier <- list(Variant.Tier = rep(NA, nrow(TumorFiltered.Normal)) )
system.time ({
for ( i in 1:nrow(TumorFiltered.Normal)) {
  Variant.Tier[i] <- ifelse( hash_table[[ TumorFiltered.Normal[i]$Somatic.Tier ]] < hash_table[[ TumorFiltered.Normal[i]$Germline.Tier ]],
                             TumorFiltered.Normal[i]$Somatic.Tier, TumorFiltered.Normal[i]$Germline.Tier)
}
})

## Save it as it takes long time ( More effeicient solution needed )
#saveRDS(Variant.Tier, "./Variant.Tier.rds")
Variant.Tier <- readRDS("./Variant.Tier.rds")
TumorFiltered.Normal$Variant.Tier <- unlist(Variant.Tier)

## Print & Save
saveRDS(TumorFiltered.Normal, 
        paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, "2.Tumor_CellLine_Variants.cosmic.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "2.Tumor_CellLine_Variants.cosmic.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## Step 3 Calculate the frequency at Patient and Sample level fro Tumor and CellLine Cohort Only
variantsFreqInCohort = TumorFiltered.Normal[, .(variantInSamples = .N, 
                                                     propInAllPatients = (length(unique(Patient.ID))/totalPatients ),
                                                     propInAllSamples = (length(unique(Sample.ID))/totalSamples   ) ) , 
                                               by = variantKey ]

setkey(variantsFreqInCohort, variantKey)
variantsFreqInCohort <- variantsFreqInCohort[order(-propInAllPatients),]
dim(variantsFreqInCohort); head(variantsFreqInCohort)
## Merge Frequency df with main df
TumorFiltered.Normal.freq <- merge( TumorFiltered.Normal, variantsFreqInCohort, by="variantKey" ); dim(TumorFiltered.Normal.freq)

## Sanity Check: Testing the "fkt" columns
min(TumorFiltered.Normal.freq$Total.coverage); min(TumorFiltered.Normal.freq$Variant.coverage); 
min(TumorFiltered.Normal.freq$VAF)           ; min(TumorFiltered.Normal.freq$MAF)

saveRDS(TumorFiltered.Normal.freq,
         paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, "3.Tumor_CellLine_Variants.freq.cosmic.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal.freq, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "3.Tumor_CellLine_Variants.freq.cosmic.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## STEP 6 filter by VAF, Variant coverage
## Mutation analysis
TumorFiltered.Normal.freq.VAF.TC.VC.MAF <- TumorFiltered.Normal.freq[
                                              Total.coverage >= 10 & Variant.coverage >= 3 & VAF >= 0.10 & MAF <= 0.01 & propInAllSamples <= 0.10 ]
## Neoantigen analysis
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen <- TumorFiltered.Normal.freq[
                                              Total.coverage >= 10 & Variant.coverage >= 3 & VAF >= 0.10 & MAF <= 0.0001 & propInAllSamples <= 0.10 ]

dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF);length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF$variantKey))
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen);length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen$variantKey))
                    

## Print and Save files                                
saveRDS(TumorFiltered.Normal.freq.VAF.TC.VC.MAF, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, 
                  "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.cosmic.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.cosmic.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

saveRDS(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen, 
        paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, 
              "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.10e4_propInTumor.LTE.10pc.Neoantigen.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.10e4_propInTumor.LTE.10pc.Neoantigen.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## STEP 7 filter by population filter 
# TumorFiltered.GermlineTierTier <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Variant.Tier][order(Variant.Tier)]
# TumorFiltered.GermlineClinvar <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Clinvar.Pathogenic][order(Clinvar.Pathogenic)]
# TumorFiltered.GermlineIntervar <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Intervar][order(Intervar)]
# TumorFiltered.GermlineHGMD <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = HGMD.Disease][order(HGMD.Disease)]

## Type of Tiering
# Tier <- "Germline.Tier"
# Tier <- "Somatic.Tier"
# Tier <- "Variant.Tier"

## No Tier
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.NoTier <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"NA\", \"\")"))
## Tier 1
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier,
                                                                                        " %in%  c(\"Tier 1.0\",\"Tier 1.1\",\"Tier 1.2\",\"Tier 1.3\",\"Tier 1.4\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.0 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.0\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.1 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.1\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.2\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.3 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.3\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.4 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.4\")"))
## Tier 2
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 2\")"))
## Tier 3
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 3\")"))
## Tier 4
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 4\")"))

## Print summary
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.NoTier)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.NoTier$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.0)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.0$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.1)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.1$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.2)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.2$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.3)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.3$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.4)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.4$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3$variantKey)) )
paste0( "Total Variants ", dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4)[1], 
        "; Unique Variants ",length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4$variantKey)) )

## Writing files
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.NoTier, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5a.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.NoTier.v2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5b.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Tier1.v2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5c.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Tier2.v2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5d.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Tier3.v2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5e.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Tier4.v2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


## Only disease causing
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[
  Clinvar.Pathogenic %in% c("Y") | HGMD.Disease %in% c("Y") | Intervar %in% c("Likely pathogenic",
                                                                              "Pathogenic",
                                                                              "Uncertain significance") ]
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4e.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.v2.DiseaseCausing.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


### Plotting
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 %>% filter(LIBRARY_TYPE %in% c("Tumor", "CellLine"))
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor)

## Sample statistics
SampleStatsMut.Tumor <- rnaseqMutationProject$validMetaDataDF %>% dplyr::filter(LIBRARY_TYPE %in% c("Tumor", "CellLine")) %>% dplyr::select(DIAGNOSIS.Alias) %>% 
                              dplyr::group_by(DIAGNOSIS.Alias) %>% 
                              dplyr::mutate(PatientSum=n()) %>% distinct()
SampleStatsMut.Tumor

### Diagnosis Wise mutation frequency 
MutationDataFilt <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor
MutationDataFiltSelect <- MutationDataFilt %>% dplyr::select(one_of("Chr","Start","End","Ref","Alt", "Gene", "DIAGNOSIS.Alias", "Patient.ID")) %>% 
  dplyr::select("Gene","Patient.ID","DIAGNOSIS.Alias") %>% 
  dplyr::group_by_(.dots = c("Gene", "DIAGNOSIS.Alias")) %>% 
  dplyr::mutate(GeneDiagByPatient = length(unique(Patient.ID))) %>% 
  dplyr::select("Gene",  "DIAGNOSIS.Alias", "GeneDiagByPatient" ) %>% distinct()
View(MutationDataFiltSelect); dim(MutationDataFiltSelect)

MutationDataFiltPercent  <-   dplyr::left_join(MutationDataFiltSelect, SampleStatsMut.Tumor, by="DIAGNOSIS.Alias" ) %>%
  dplyr::select(one_of(c("Gene", "DIAGNOSIS.Alias", "GeneDiagByPatient", "PatientSum"))) %>%
  dplyr::mutate(Percent=( (GeneDiagByPatient/PatientSum)*100)) %>% unique() %>%
  dplyr::select(Gene, DIAGNOSIS.Alias, Percent) %>% 
  tidyr::spread(DIAGNOSIS.Alias, Percent) %>% data.frame()
MutationDataFiltPercent[is.na(MutationDataFiltPercent)] <- 0
View(MutationDataFiltPercent);dim(MutationDataFiltPercent)

#### By Gene mutation function frequency
MutationDataFilt <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor
MutationDataFuncCbio <- MutationDataFilt %>% dplyr::select(one_of("Chr","Start","End","Ref","Alt", "Gene", "Exonic.function", "Patient.ID")) %>% 
  dplyr::group_by_(.dots=c("Chr","Start","End","Ref","Alt", "Gene", "Exonic.function", "Patient.ID") ) %>% 
  unique() %>% ungroup()  %>% 
  dplyr::group_by_(.dots=c("Gene", "Exonic.function") ) %>%
  dplyr::mutate(Count=n()) 
dim(MutationDataFuncCbio)
MutationDataFunc    <- MutationDataFuncCbio %>%
  dplyr::select(one_of("Gene", "Exonic.function")) %>% 
  dplyr::group_by_(.dots=c("Gene", "Exonic.function") ) %>%
  dplyr::mutate(Count=n()) %>%  unique() %>% ungroup()  %>%
  tidyr::spread(Exonic.function, Count) %>% data.frame()
MutationDataFunc[is.na(MutationDataFunc)] <- 0
MutationDataFunc$Sum = rowSums(MutationDataFunc[,-1])
dim(MutationDataFunc);View(MutationDataFunc)

### Join the above two data frames
mutationFreqPer <- dplyr::full_join(MutationDataFiltPercent, MutationDataFunc, by="Gene") %>% dplyr::arrange(-Sum) %>%
                      dplyr::mutate(Gene = factor(Gene, ordered = TRUE))
write.table(mutationFreqPer, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "6.mutation.Freq.Per.Diagnosis.and.Gene.v2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


mutationFreqPer.RM.Single <- mutationFreqPer %>% filter(Sum >=2) %>% dplyr::select(1:14)
dim(mutationFreqPer.RM.Single);View(mutationFreqPer.RM.Single)
mutationFreqPer.RM.Single.tidy.log <- gather(mutationFreqPer.RM.Single, key="Diagnosis", value="Frequency", 2:14) %>% 
         mutate( 
                #Frequency = log2(Frequency+1)/10,
                Gene = factor(Gene, levels=rev(mutationFreqPer.RM.Single$Gene))
                )


## Main heatmap
HeatMap <- ggplot(mutationFreqPer.RM.Single.tidy.log, aes(Diagnosis,Gene)) + 
  geom_tile(aes(fill = Frequency), colour = "grey") + 
  scale_fill_gradient(low = "white",high = "darkred", name="Fraction of patients") +
  scale_x_discrete(expand = c(0, 0)) +
  #scale_y_discrete(expand = c(0, 0)) +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        axis.text=element_text(face="bold", size=10),
        axis.text.x = element_text(angle = 90),
        legend.position="bottom",
        legend.text =element_text(face="bold", size=10),
        legend.key.width =unit(3,"line"))

## Side Frequency Plot
mutationType.RM.Single <- mutationFreqPer %>% filter(Sum >=2) %>% arrange(desc(Sum)) %>% dplyr::select(c(1,16:22))
mutationType.RM.Single.tidy.log <- gather(mutationType.RM.Single, key="Type", value="Percent", 2:7) %>% 
  mutate(Percent = Percent,
         Gene = factor(Gene, levels=rev(mutationFreqPer.RM.Single$Gene)),
         Type      = factor(Type,
                            levels = c("stopgain","splicing",
                                       "frameshift.insertion","frameshift.deletion",
                                       "nonframeshift.deletion","nonframeshift.insertion",
                                       "nonsynonymous.SNV"), ordered = TRUE)) %>% 
  dplyr::arrange(desc(Gene), Type)
head(mutationType.RM.Single.tidy.log)

#myColors <- setNames( c("#31a354", "#d2a679", "#669999", "#bf4040", "#810f7c", "#7094db"), 
myColors <- setNames( c("#31a354", "#d2a679", "#7094db", "#669999", "#bf4040", "#810f7c", "#737373"), 
                      c("frameshift.deletion","frameshift.insertion","nonframeshift.deletion","nonframeshift.insertion","nonsynonymous.SNV","splicing","stopgain"))

# CbioColors <- setNames(c("#636363", "#636363", "#994d00", "#31a354", "#636363", "#636363" ), 
#                              c("frameshift.deletion","frameshift.insertion","nonframeshift.deletion","nonsynonymous.SNV","splicing","stopgain"))

customColors <- myColors

## green orange/brown violet red brinjal skyblue
StackedBar <- ggplot(mutationType.RM.Single.tidy.log, aes(x=Gene, y=Percent, fill=Type) ) +
  #ggplot(mutationType.RM.Single.tidy.log, aes(x=Gene, y=Percent, fill=Type) ) +
  geom_bar(stat="identity") + coord_flip() +
  scale_fill_manual(values = customColors, name = "Mutation Type", 
                    guide = guide_legend( 
                      direction = "vertical", 
                      reverse = TRUE))+
  scale_x_discrete(expand = c(0, 0)) +
  labs(x="SNV Counts") +
  theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), 
        panel.grid.minor =  element_line(colour="grey"),
        axis.text=element_text(face="bold", size=10),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank(),
        #legend.position="bottom",
        legend.key.width =unit(3,"line")) +
  scale_y_continuous(minor_breaks = seq(0,30,5))

## Bind both of them
ggarrange(HeatMap, StackedBar, ncol=2)



