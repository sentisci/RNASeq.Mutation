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
patientsEqual = length(unique(mergeVCFObjectDiagnosis$Patient.ID)) == length(unique(rnaseqMutationProject$validMetaDataDF$Patient.ID))
samplesEqual  = length(unique(mergeVCFObjectDiagnosis$Sample.ID)) == length(unique(rnaseqMutationProject$validMetaDataDF$Sample.ID))

if(! patientsEqual | ! samplesEqual ) { stop(paste0("Mismatch of Patients and/or Samples between VCF files and Metadata file !! ",
                                                       ", Total patients from VCF: ", length(unique(mergeVCFObjectDiagnosis$Patient.ID)),
                                                       ", Total patients from Metadata: ", length(unique(rnaseqMutationProject$validMetaDataDF$Patient.ID)),
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
## Printing VAriant analysis
print(paste(" Total Variants in the cohort after removing non-coding,splicing aand intronic: ", length(mergeVCFObjectDiagnosisVA$variantKey), 
            " & Total Variants in the cohort after removing non-coding,splicing aand intronic: ", length(unique(mergeVCFObjectDiagnosisVA$variantKey))))

## Printing cosmic
print(paste(" Total Variants in the cohort after removing non-coding,splicing aand intronic: ", length(mergeVCFObjectDiagnosisCosmic$variantKey), 
            " & Total Variants in the cohort after removing non-coding,splicing aand intronic: ", length(unique(mergeVCFObjectDiagnosisCosmic$variantKey))))

## Print & Save
saveRDS(mergeVCFObjectDiagnosisVA,
        paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, "1.All.variants.v3.RDS", sep="/"))
write.table(mergeVCFObjectDiagnosisVA, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "1.All.variants.v3.txt", sep="/"),
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
saveRDS(Variant.Tier, "./Variant.Tier.cosmic.rds")
## Variant.Tier <- readRDS("./Variant.Tier.rds")
TumorFiltered.Normal$Variant.Tier <- unlist(Variant.Tier)

## Print & Save
saveRDS(TumorFiltered.Normal, 
        paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, "2.Tumor_CellLine_Variants.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "2.Tumor_CellLine_Variants.v3.txt", sep="/"),
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
         paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, "3.Tumor_CellLine_Variants.freq.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal.freq, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "3.Tumor_CellLine_Variants.freq.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## STEP 6 filter by VAF, Variant coverage
## Mutation analysis
TumorFiltered.Normal.freq.VAF.TC.VC.MAF <- TumorFiltered.Normal.freq[
                                                    Total.coverage >= 10 & 
                                                    Variant.coverage >= 3 & 
                                                    VAF >= 0.10 & 
                                                    MAF <= 0.01 & 
                                                    propInAllSamples <= 0.10 ]

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.HighConfIndels <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[which( Exonic.function %in% 
                                                                                                             c("frameshift deletion","frameshift insertion",
                                                                                                               "nonframeshift deletion","nonframeshift insertion") & 
                                                                                                             propInAllPatients <= 0.01 ), ]
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.HighConfSNVs <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[which( Exonic.function %in% c("splicing", "nonsynonymous SNV", "stopgain", "stoploss", "intronic", "UTR3")), ]

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels <- rbind(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.HighConfSNVs, TumorFiltered.Normal.freq.VAF.TC.VC.MAF.HighConfIndels)

## Review the above step
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels)
length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels$variantKey))
unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels$Exonic.function)

## Neoantigen analysis
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen <- TumorFiltered.Normal.freq[
                                                  Total.coverage >= 10 & 
                                                  Variant.coverage >= 3 & 
                                                  VAF >= 0.10 & 
                                                  MAF <= 0.0001 & 
                                                  propInAllSamples <= 0.10 ]

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.HighConfIndels <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen[which( Exonic.function %in% 
                                                                                                           c("frameshift deletion","frameshift insertion",
                                                                                                             "nonframeshift deletion","nonframeshift insertion") & 
                                                                                                           propInAllPatients <= 0.01 ), ]
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.HighConfSNVs <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen[which( Exonic.function %in% c("splicing", "nonsynonymous SNV", "stopgain", "stoploss", "intronic", "UTR3")), ]

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.indels <- rbind(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.HighConfIndels, 
                                                                   TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.HighConfSNVs )


## Review the above step
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.indels);
length(unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.indels$variantKey))
unique(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.indels$Exonic.function)
                    

## Print and Save files                                
saveRDS(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, 
                  "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.Indels.LTE.1pc.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.Indels.LTE.1pc.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

saveRDS(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.indels, 
        paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, 
              "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.10e4_propInTumor.LTE.10pc.Indels.LTE.1pc.Neoantigen.v3.RDS", sep="/"))
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Neoantigen.indels, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.10e4_propInTumor.LTE.10pc.Indels.LTE.1pc.Neoantigen.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

## STEP 7 filter by population filter 
# TumorFiltered.GermlineTierTier <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Variant.Tier][order(Variant.Tier)]
# TumorFiltered.GermlineClinvar <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Clinvar.Pathogenic][order(Clinvar.Pathogenic)]
# TumorFiltered.GermlineIntervar <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = Intervar][order(Intervar)]
# TumorFiltered.GermlineHGMD <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[, .N , by = HGMD.Disease][order(HGMD.Disease)]

## Type of Tiering
# Tier <- "Germline.Tier"
# Tier <- "Somatic.Tier"
Tier <- "Variant.Tier"

## No Tier
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.NoTier <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"NA\", \"\")"))
## Tier 1
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier,
                                                                                        " %in%  c(\"Tier 1.0\",\"Tier 1.1\",\"Tier 1.2\",\"Tier 1.3\",\"Tier 1.4\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.0 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.0\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.1 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.1\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.2\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.3 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.3\")"))
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.4 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 1.4\")"))
## Tier 2
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 2\")"))
## Tier 3
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 3\")"))
## Tier 4
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.indels %>% filter_(.dots=paste0(Tier, " %in%  c(\"Tier 4\")"))

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
                  "5a.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Indels.LTE.1pc.NoTier.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5b.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Indels.LTE.1pc.Tier1.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5c.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Indels.LTE.1pc.Tier2.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5d.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Indels.LTE.1pc.Tier3.v3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "5e.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10.Indels.LTE.1pc.Tier4.v3.txt", sep="/"),
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

###### OR get reviewd & selected variants from a file #######
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor <- read.csv("../RNASeq.Mutation.data/Fig1_manual_filter_tier1_minusfaultygermline_addselectTier2Tier3_indelback_05_21_19_manual_indel_splice_Selected.txt",
                                                                  sep="\t", header = T)
dim( TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor)
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor %>% 
                                                              dplyr::filter(grepl("Tier 1.*", Variant.Tier))
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor)

## Sample statistics
SampleStatsMut.Tumor <- rnaseqMutationProject$validMetaDataDF %>% dplyr::filter(LIBRARY_TYPE %in% c("Tumor", "CellLine")) %>% 
                              dplyr::select(DIAGNOSIS.Alias, Patient.ID.updated) %>%
                              dplyr::group_by(DIAGNOSIS.Alias) %>% 
                              summarise(PatientSum=n_distinct(Patient.ID.updated))
SampleStatsMut.Tumor

### Diagnosis Wise mutation frequency 
MutationDataFilt <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor
MutationDataFiltSelect <- MutationDataFilt %>% dplyr::select(one_of("Chr","Start","End","Ref","Alt", "Gene", "DIAGNOSIS.Alias", "Patient.ID.updated")) %>% 
  dplyr::select("Gene","Patient.ID.updated","DIAGNOSIS.Alias") %>% 
  dplyr::group_by_(.dots = c("Gene", "DIAGNOSIS.Alias")) %>% 
  dplyr::mutate(GeneDiagByPatient = length(unique(Patient.ID.updated))) %>% 
  dplyr::select("Gene",  "DIAGNOSIS.Alias", "GeneDiagByPatient" ) %>% distinct()
View(MutationDataFiltSelect); dim(MutationDataFiltSelect)

MutationDataFiltPercent  <-   dplyr::left_join(MutationDataFiltSelect, SampleStatsMut.Tumor, by="DIAGNOSIS.Alias" ) %>%
  dplyr::select(one_of(c("Gene", "DIAGNOSIS.Alias", "GeneDiagByPatient", "PatientSum"))) %>%
  dplyr::mutate(Percent=( (GeneDiagByPatient/PatientSum)*100)) %>% unique() %>%
  dplyr::select(Gene, DIAGNOSIS.Alias, Percent) %>% 
  tidyr::spread(DIAGNOSIS.Alias, Percent) %>% data.frame()
MutationDataFiltPercent[is.na(MutationDataFiltPercent)] <- 0
View(MutationDataFiltPercent);dim(MutationDataFiltPercent)

### memo sort #####
# memoSort = function(M = NA) {
#   geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
#   scoreCol <- function(x) {
#     score <- 0;
#     for(i in 1:length(x)) {
#       if(x[i]) {
#         score <- score + 2^(length(x)-i);
#       }
#     }
#     return(score);
#   }
#   scores <- apply(M[geneOrder, ], 2, scoreCol);
#   sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
#   return(M[geneOrder, sampleOrder]);
# }
# test <- MutationDataFiltPercent
# rownames(test) <- MutationDataFiltPercent$Gene
# test <- test[,-c(1)]
# MutationDataFiltPercent.memoSort <- memoSort(test)
# View(MutationDataFiltPercent.memoSort)


#### By Gene mutation function frequency ####
MutationDataFilt <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor
MutationDataFuncCbio <- MutationDataFilt %>% dplyr::select(one_of("Chr","Start","End","Ref","Alt", "Gene", "Exonic.function", "Patient.ID.updated")) %>% 
  dplyr::group_by_(.dots=c("Chr","Start","End","Ref","Alt", "Gene", "Exonic.function", "Patient.ID.updated") ) %>% 
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


mutationFreqPer.RM.Single <- mutationFreqPer %>% filter(Sum >=3) %>% dplyr::select(1:14)
dim(mutationFreqPer.RM.Single);View(mutationFreqPer.RM.Single)
mutationFreqPer.RM.Single.tidy.log <- gather(mutationFreqPer.RM.Single, key="Diagnosis", value="Frequency", 2:14) %>% 
         mutate( 
                #Frequency = log2(Frequency+1),
                Gene = factor(Gene, levels=rev(mutationFreqPer.RM.Single$Gene))
                )

## Main heatmap
# test <- mutationFreqPer.RM.Single.tidy.log
# test$Frequency <- log2(test$Frequency + 1)
HeatMap <- ggplot(mutationFreqPer.RM.Single.tidy.log, aes(Diagnosis,Gene)) + 
  geom_tile(aes(fill = Frequency), colour = "grey") + 
  scale_fill_gradientn(#low = "white",high = "darkred", name="Fraction of patients", 
                      colours = c("white", "darkred"),
                      #values = scales::rescale(c(-0.5, -0.05, 0, 0.05, 0.5))
                      values = scales::rescale(c(-0.5, -0.05, 0, seq(1, 20, 1), seq(30,55,5)))
                      ### Generate breaks ## paste(seq(-3.3, 5.7, 0.3), collapse = ",") ##
                      # breaks=c(-3.3,-2.4,-1.5,-0.6,0.3,1.2,2.1,3,3.9,4.8,5.7),
                      ### Generate labels ## paste(as.numeric(formatC(2^seq(-3.3, 5.7, 0.3), digits = 1))-0.1, collapse = ",") ##
                      # labels=c(0,0.1,0.3,0.6,1,2,4,8,10,30,50)
                      ) +
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
mutationType.RM.Single.tidy.log <- gather(mutationType.RM.Single, key="Type", value="Counts", 2:8) %>% 
  mutate(Counts = Counts,
         Gene = factor(Gene, levels=rev(mutationFreqPer.RM.Single$Gene),ordered = TRUE),
         Type      = factor(Type,
                            levels = c("stopgain","splicing",
                                       "frameshift.insertion","frameshift.deletion",
                                       "nonframeshift.deletion","nonframeshift.insertion",
                                       "nonsynonymous.SNV"),ordered = TRUE)) %>% 
  dplyr::arrange(desc(Gene), Type)
head(mutationType.RM.Single.tidy.log)

#myColors <- setNames( c("#31a354", "#d2a679", "#669999", "#bf4040", "#810f7c", "#7094db"), 
myColors <- setNames( c("#31a354", "#d2a679", "#7094db", "#669999", "#bf4040", "#810f7c", "#737373"), 
                      c("frameshift.deletion","frameshift.insertion","nonframeshift.deletion","nonframeshift.insertion","nonsynonymous.SNV","splicing","stopgain"))

# CbioColors <- setNames(c("#636363", "#636363", "#994d00", "#31a354", "#636363", "#636363" ), 
#                              c("frameshift.deletion","frameshift.insertion","nonframeshift.deletion","nonsynonymous.SNV","splicing","stopgain"))

customColors <- myColors

## green orange/brown violet red brinjal skyblue
StackedBar <- ggplot(mutationType.RM.Single.tidy.log, aes(x=Gene, y=Counts, fill=Type) ) +
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
  scale_y_continuous(minor_breaks = seq(0,105,10))

## Bind both of them
ggarrange(HeatMap, StackedBar, ncol=2)

pdf("C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/Figures/Figure.1b.Variant_percentage_and_count.Tier 1.v8.NL.v3.pdf", height = 25, width = 10 )
ggarrange(HeatMap, StackedBar, ncol=2)
dev.off()


################### Neo-antigen analysis #######################

#### STEP 1 : Make filtered VCF files #####

## Read the selected variant dataset which has variants across all samples in the study cohort.
filtered.Final.VCF <- read.csv("C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/FigureData/4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.10e4_propInTumor.LTE.10pc.Indels.LTE.1pc.Neoantigen.v3.txt",
                                 sep="\t", header = T)

## List all samples in the study cohort
allSamples <- unique(as.character(filtered.Final.VCF$Sample.ID)) ; length(allSamples)
## Convert actual sampleIDs to Biowulf sampleIDs
sampleIDs.DF <- rnaseqMutationProject$metaDataDF %>% dplyr::filter(Sample.ID %in%  allSamples ) %>% dplyr::select(Sample.ID, Sample.Biowulf.ID)

## Custom function to subset raw vcf to new vcf with selected variants only.
filterRawVCFForVariants <- function(x){
  sampleids=x
  ## start analysis for only one sample
  sample.filtered.vcf  <- filtered.Final.VCF %>% dplyr::filter(Sample.ID == x[1] ); 
  totalVariantsFinal = dim(sample.filtered.vcf)

  ## Read the raw vcf
  rawVCF <- readVcf(paste0("K:/projects/Sivasish/pvactools/VariantVCFs/",x[2], ".HC_RNASeq.raw.vcf") )
  
  ## Method 1 Subsetting / Filtering 
  ### Filter criteria
  q=GRanges(seqnames=sample.filtered.vcf$Chr,
            ranges=IRanges(start=sample.filtered.vcf$Start, end = sample.filtered.vcf$End))
  filteredVCF <- subsetByOverlaps(rawVCF, q); rowRanges(filteredVCF)
  totalVariantsRaw = dim(filteredVCF)
  
  ## Sanity Check
  ## The range of indel will be different in vcf file vs annotated vcf file.
  ## if ( totalVariantsFinal[1] != totalVariantsRaw[1] ) stop("Variants doesn't match")
  
  print(paste(x[1], x[2], paste0(" Variants in filtered file ", totalVariantsFinal[1], collapse = ":"), 
                 paste0(" Variants in raw file ", totalVariantsRaw[1], collapse = ":")))
  ### Writing the filtered raw VCF file
  writeVcf(filteredVCF, paste0("C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/FigureData/filteredvcfs/",x[2], "filtered.raw.vcf") )
}

## Execute the above function: To create  new vcfs for all samples present in the input file
outlist <- apply(sampleIDs.DF, 1, filterRawVCFForVariants)


### Neoantigens from Variants
emptyDF <- data.frame(Gene.Name=c(0), Mutation=c(0), Protein.Position=c("NA"), 
                      HGVSc=c("NF"), HGVSp=c("NF"), HLA.Allele=c("NF"), MT.Epitope.Seq=c("NF"), 
                      MT.IC50=c("NF"),	WT.IC50=c("NF"),	Fold.Change=c("NF"),	Tumor.DNA.Depth=c("NF"),
                      Tumor.DNA.VAF=c("NF"),	Tumor.RNA.Depth=c("NF"),	Tumor.RNA.VAF=c("NF"),	
                      Gene.Expression=c("NF"),	Score=c("NF")
)

### List files and read data into a single data matrix ####
folder_name = "C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/mutation_neoantigen_files/"
fileList <- list.files(folder_name)
#fileList <- c("convert.Sample_RMS248_C14C7ACXX.clones.txt")
AllNeoantigenData             <- rbindlist( lapply(fileList, function(x){
  print(x)
  exomeData <- read.csv( paste(folder_name, x, sep=""), sep="\t", header = TRUE, fill=TRUE )
  if(nrow(exomeData)>0){
    exomeData$SampleName <- x
  } else {
    emptyDF$SampleName <- c(x)
    exomeData <- emptyDF
  }
  return(exomeData)
}) )
write.table(AllNeoantigenData,"C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/mutation_neoantigen_files/AllNeoantigenData.txt", 
            sep="\t", col.names = T, row.names = F, quote = F)


################### Mutation Signature analysis #######################
library(deconstructSigs)
mut_data <- readRDS(paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirRDSDir, 
                          "4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.Indels.LTE.1pc.v3.RDS", sep="/"))
mut_data_sample <- mut_data %>% dplyr::filter(!Sample.ID %in% c("SS019tumor_T_D1UE8ACXX", "NCIEWS5000_T_C28PUACXX", "EWS114tumor_T_C1P2WACXX", "EWS124tumor_T_D1UALACXX"))
mut_data_sample_input <- mut_data_sample %>% dplyr::select(Sample.ID, Chr, Start, End, Ref, Alt); View(mut_data_sample_input)

## Count total mutations per sample
mut_data_bySample_Count <- mut_data %>% group_by(Sample.ID) %>% 
  dplyr::mutate(TotalMutations = n()) %>% 
  dplyr::select(Sample.ID, DIAGNOSIS.Alias, TotalMutations) %>% 
  distinct(); View(mut_data_bySample_Count)

## Make input for dConstructSig
sigs.input <- mut.to.sigs.input(mut.ref = mut_data_sample_input,
                                sample.id = "Sample.ID",
                                chr = "Chr",
                                pos = "Start",
                                ref = "Ref",
                                alt = "Alt")

## Total mutations per sample identifies by dConstructSig
sigs.input.bySample_Sum <- apply(sigs.input, 1, sum) %>% data.frame() %>% 
  tibble::rownames_to_column(var="Sample.ID"); 
colnames(sigs.input.bySample_Sum)<- c("Sample.ID", "dConstruct_TotalMutations"); View(sigs.input.bySample_Sum)

## Save the input for deconstructSig
samples <- rownames(sigs.input)
write.table(sigs.input , "C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/sigs.input.allSamples.GTE.50.mutations.txt", sep="\t")

## Sanity check; To see the difference in Total mutation count between dConstructSig count vs Actual count
mut_data_dconst <- dplyr::full_join(mut_data_bySample_Count, sigs.input.bySample_Sum, by="Sample.ID"); View(mut_data_dconst)
write.table(mut_data_dconst , "C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/Actual_vs_dConst.allSamples.GTE.50.mutations.txt", sep="\t")

## DconstrucSig for all samples using Cosmic signatures as REference
cosmic_list <- sapply(samples, function(x){
  cosmic = whichSignatures(tumor.ref = sigs.input,
                           signatures.ref = signatures.cosmic,
                           sample.id = x,
                           contexts.needed = TRUE,
                           tri.counts.method = 'exome2genome')
  return(cosmic)
})

## Extract different items from the above output
cosmic.df_weights <- do.call(rbind, cosmic_list[1,]); View(cosmic.df_weights)
cosmic.df_tumor <- do.call(rbind, cosmic_list[2,]); View(cosmic.df_tumor)
cosmic.df_product <- do.call(rbind, cosmic_list[3,]); View(cosmic.df_product)
cosmic.df_diff <- do.call(rbind, cosmic_list[4,]); View(cosmic.df_diff)
cosmic.df_unknown <- do.call(rbind, cosmic_list[5,]); View(cosmic.df_unknown)

### MEthod 1
## Determine percent mutation contribution by each signature
# mutation_weight_sum <- sum(as.numeric(unlist(apply(cosmic.df_weights, 2, sum))))
# cosmic.df_weights_mut <- cosmic.df_weights %>% tibble::rownames_to_column(var = "Sample.ID")
# cosmic.df_weights_diag <- full_join(cosmic.df_weights_mut, mut_data_bySample_Count, by="Sample.ID")
# cosmic.df_weights_diag <- cosmic.df_weights_diag[complete.cases(cosmic.df_weights_diag),]
# 
# mutation_weight_sum_byDiagnosis <- cosmic.df_weights_diag %>% dplyr::select(-Sample.ID) %>%
#                                         group_by(DIAGNOSIS.Alias) %>% 
#                                         summarise_each(sum);
# mutation_weight_sum_byDiagnosis <- (mutation_weight_sum_byDiagnosis[,-1]/mutation_weight_sum)*100; View(mutation_weight_sum_byDiagnosis)

## MEthod 2
sigs.input.bySample_Sum_annot <- sigs.input.bySample_Sum %>% tibble::column_to_rownames("Sample.ID")
weight_mutationCount <- t( sapply(seq(1,784), function(x){
  return(t(sigs.input.bySample_Sum_annot[x,]*cosmic.df_weights[x,]))
  #return(t(1*cosmic.df_weights[x,]))
}) )
colnames(weight_mutationCount) <- colnames(cosmic.df_weights)
rownames(weight_mutationCount) <- rownames(cosmic.df_weights)
#View(weight_mutationCount)

weight_mutationCount_annot <- weight_mutationCount %>% data.frame() %>% tibble::rownames_to_column("Sample.ID")
weight_mutationCount_final <- dplyr::full_join(weight_mutationCount_annot,mut_data_dconst[,c("Sample.ID", "DIAGNOSIS.Alias")], by="Sample.ID"); 
View(weight_mutationCount_final)

write.table(weight_mutationCount_final, "../RNASeq.Mutation.data/Weight.txt", sep="\t", quote = FALSE, row.names = FALSE)

weight_mutationCount_final <- weight_mutationCount_final %>% tibble::column_to_rownames("Sample.ID")
weight_mutationCount_final_freq <- t( apply(weight_mutationCount_final[,-31], 1, function(x){
  sum = sum(as.numeric(x))
  fraction = x/sum
  percent = fraction*100
  return(percent)
}))

weight_mutationCount_final_percent <- weight_mutationCount_final_freq[complete.cases(weight_mutationCount_final_freq),]
View(weight_mutationCount_final_percent)

weight_mutationCount_percent_annot <- weight_mutationCount_final_percent %>% data.frame() %>% tibble::rownames_to_column("Sample.ID")
weight_mutationCount_percent_annot <- dplyr::full_join(weight_mutationCount_percent_annot,mut_data_dconst[,c("Sample.ID", "DIAGNOSIS.Alias")], by="Sample.ID"); 
weight_mutationCount_percent_annot <- weight_mutationCount_percent_annot %>% tibble::column_to_rownames("Sample.ID")
weight_mutationCount_percent_annot <- weight_mutationCount_percent_annot[complete.cases(weight_mutationCount_percent_annot),]
View(weight_mutationCount_percent_annot)

write.table(weight_mutationCount_percent_annot, "../RNASeq.Mutation.data/Weight_Total_mutations_percent.txt", sep="\t", quote = FALSE )

mutationBurden_By_Diagnosis <- weight_mutationCount_percent_annot %>% dplyr::group_by(DIAGNOSIS.Alias) %>% summarise_each(sum)
mutationBurden_By_Diagnosis <- mutationBurden_By_Diagnosis %>% tibble::column_to_rownames("DIAGNOSIS.Alias")

write.table(mutationBurden_By_Diagnosis, "../RNASeq.Mutation.data/Total_percent_(contribution)_sum_by_diagnosis.txt", sep="\t", quote = FALSE)

mutationBurden_By_Diagnosis_freq <- t( apply(mutationBurden_By_Diagnosis, 1, function(x){
  sum = sum(as.numeric(x))
  fraction = x/sum
  percent = fraction*100
  return(percent)
})) %>% data.frame()
mutationBurden_By_Diagnosis_freq <- mutationBurden_By_Diagnosis_freq %>% t() %>% data.frame() %>% tibble::rownames_to_column(var="Signatures")
View(mutationBurden_By_Diagnosis_freq)
write.table(mutationBurden_By_Diagnosis_freq, "../RNASeq.Mutation.data/Total_(contribution)_sum_by_diagnosis_percent.txt", sep="\t", quote = FALSE)

## Make a bar plot
mutationBurden_By_Diagnosis_freq[mutationBurden_By_Diagnosis_freq<4]=0
pallete <- distinctColorPalette(30)
mutationBurden_By_Diagnosis_freq_tidy <- tidyr::gather(mutationBurden_By_Diagnosis_freq, "Diagnosis", "PercentMutation", 2:15); 
View(mutationBurden_By_Diagnosis_freq_tidy)
ggplot(data = mutationBurden_By_Diagnosis_freq_tidy, aes(x = Diagnosis, y = PercentMutation, fill= Signatures )) + 
  geom_bar(stat="identity") + coord_flip() + scale_fill_manual(values=pallete)



## Like the above but only for one sample using Cosmic as reference
cosmic_exp = whichSignatures(tumor.ref = sigs.input,
                         signatures.ref = signatures.cosmic,
                         sample.id = "NCI0132tumor4_T_C28D2ACXX",
                         contexts.needed = TRUE,
                         tri.counts.method = 'exome2genome')

plotSignatures(cosmic_exp, sub='Mutational Signature Based on COSMIC')
makePie(cosmic_exp, sub='Mutational Signature Based on COSMIC')

## Like the above but only for one sample using Nature publications as reference
nature = whichSignatures(tumor.ref = sigs.input,
                         signatures.ref = signatures.nature2013,
                         sample.id = "NB2050tumor_T_D1T6TACXX",
                         contexts.needed = TRUE,
                         tri.counts.method = 'default')

plotSignatures(nature, sub='Mutational Signature based on Nature 2013--23945592')
makePie(nature, sub='Mutational Signature based on Nature 2013--23945592')

#### Retrieving the dataset
jun.data.set <- read.csv("C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/Fig1_manual_filter_tier1_minusfaultygermline_addselectTier2Tier3_indelback_05_21_19_manual_indel_splice_Selected.txt",
                         sep="\t", header = TRUE)
jun.data.set$keyNew <- paste0(jun.data.set$Patient.ID,
                           jun.data.set$Chr,
                           jun.data.set$Start,
                           jun.data.set$End,
                           jun.data.set$Ref,
                           jun.data.set$Alt,
                           jun.data.set$Gene,
                           jun.data.set$VAF,
                           jun.data.set$Total.coverage,
                           jun.data.set$Variant.coverage)

actual.data.set <- read.csv("C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/outputTXTOutput/4.Tumor_CellLine_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_propInTumor.LTE.10pc.Indels.LTE.1pc.v3.txt",
                         sep="\t", header = TRUE)
actual.data.set$keyNew <- paste0(actual.data.set$Patient.ID,
                           actual.data.set$Chr,
                           actual.data.set$Start,
                           actual.data.set$End,
                           actual.data.set$Ref,
                           actual.data.set$Alt,
                           actual.data.set$Gene,
                           actual.data.set$VAF,
                           actual.data.set$Total.coverage,
                           actual.data.set$Variant.coverage)

jun.actual <- dplyr::left_join(jun.data.set, actual.data.set, by="keyNew")
write.table(jun.actual, "../RNASeq.Mutation.data/jun.actual.retrieved.txt", sep = "\t", quote = FALSE,row.names = F, col.names = T )


#### MEMO SORT ####
memoSort = function(M = NA) {
  geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  scores <- apply(M[geneOrder, ], 2, scoreCol);
  sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
  return(M[geneOrder, sampleOrder]);
}

metaData <- read.csv("../RNASeq.Mutation.data/MetadataMapper.v3.txt", sep="\t", header = T)
primaryRelapse <- function(x) {

  patientname = x
  print(patientname)
  ## Primary and relapse
  DRCT001.meta   <- metaData %>% dplyr::filter(Patient.ID.updated == patientname)
  patients <- as.character(unique(DRCT001.meta$Patient.ID))
  
  ## Using actual dataset
  DRCT001.actual <- actual.data.set %>% dplyr::filter(Patient.ID %in% patients)
  print(dim(DRCT001.actual))
  if(nrow(DRCT001.actual) > 1 ) {
      DRCT001.actual$variantKey <- paste(DRCT001.actual$Gene, paste(DRCT001.actual$Ref, DRCT001.actual$Alt, DRCT001.actual$AAChange,
                                                                    DRCT001.actual$Start, sep = "/"), sep="#")
      DRCT001.actual$dummyVal <- rep(1, nrow(DRCT001.actual))
      DRCT001.actual.heatmap <- reshape2::dcast(DRCT001.actual, Exonic.function + variantKey ~ Sample.ID, value.var = "dummyVal",  fill = 0);
      DRCT001.actual.heatmap$Name <- paste(DRCT001.actual.heatmap$Exonic.function, DRCT001.actual.heatmap$variantKey, sep="#")
      DRCT001.actual.heatmap <- DRCT001.actual.heatmap[,c(3:ncol(DRCT001.actual.heatmap))]
      DRCT001.actual.heatmap %<>% tibble::column_to_rownames(var = "Name")
      DRCT001.actual.heatmap <- memoSort(M = DRCT001.actual.heatmap)
      DRCT001.actual.heatmap$Sum <- apply(DRCT001.actual.heatmap,1, function(x) { 
        x <- as.numeric(x)
        #return(sum(x[3:length(x)]))
        return(sum(x))
      })
      
      DRCT001.actual.heatmap %<>% tibble::rownames_to_column(var = "Name")
      DRCT001.actual.heatmap %<>% tidyr::separate("Name", c("Exonic.function", "Gene", "AAChange"), sep="#")
      #DRCT001.actual.heatmap %<>% dplyr::arrange(-Sum)
      head( DRCT001.actual.heatmap)
      write.table(DRCT001.actual.heatmap, paste0("../RNASeq.Mutation.data/primary.relapse/",patientname,".actual.heatmap.pimary.relapse.heat.txt"), sep = "\t", row.names = FALSE,
                  col.names = TRUE)
  }
  
  
  ## Using Jun's filtered dataset
  DRCT001.filtered <- jun.actual %>% dplyr::filter(Patient.ID.x %in% patients)
  print(dim(DRCT001.filtered))
  if(nrow(DRCT001.filtered) > 1 ) {
      DRCT001.filtered$variantKey.x <- paste(DRCT001.filtered$Gene.x, paste(DRCT001.filtered$Ref.x, DRCT001.filtered$Alt.x, DRCT001.filtered$AAChange.x, 
                                                                          DRCT001.filtered$Start.x, sep = "/"), sep="#")
      DRCT001.filtered$dummyVal <- rep(1, nrow(DRCT001.filtered))
      DRCT001.filtered.heatmap <- reshape2::dcast(DRCT001.filtered, Exonic.function.x + variantKey.x ~ Sample.ID, value.var = "dummyVal",  fill = 0);
      DRCT001.filtered.heatmap$Name <- paste(DRCT001.filtered.heatmap$Exonic.function, DRCT001.filtered.heatmap$variantKey, sep="#")
      DRCT001.filtered.heatmap <- DRCT001.filtered.heatmap[,c(3:ncol(DRCT001.filtered.heatmap))]
      DRCT001.filtered.heatmap %<>% tibble::column_to_rownames(var = "Name")
      DRCT001.filtered.heatmap <- memoSort(M = DRCT001.filtered.heatmap)
      DRCT001.filtered.heatmap$Sum <- apply(DRCT001.filtered.heatmap,1, function(x) { 
        x <- as.numeric(x)
        #return(sum(x[3:length(x)]))
        return(sum(x))
      })
      DRCT001.filtered.heatmap %<>% tibble::rownames_to_column(var = "Name")
      DRCT001.filtered.heatmap %<>% tidyr::separate("Name", c("Exonic.function", "Gene", "AAChange"), sep="#")
      #DRCT001.filtered.heatmap %<>% dplyr::arrange(-Sum)
      head( DRCT001.filtered.heatmap)
      write.table(DRCT001.filtered.heatmap, paste0("../RNASeq.Mutation.data/primary.relapse/",patientname,".filtered.heatmap.pimary.relapse.heat.txt"), sep = "\t", row.names = FALSE,
                  col.names = TRUE)
  }
}

primaryrelapsePatients = c(
  "DRCT001","EWS5000","EWS5201","MSKCC5000",
  "MSKCC5001","MSKCC5002","MSKCC5003","MSKCC5004",
  "MSKCC5005",# "NB2305",
  "NB2307","NB2308",
  "NB2310","NB2318","NB2319","NCI0017",
  "NCI0039","NCI0064","NCI0097","NCI0107",
  "NCI0132","NCI0135","NCI0150","NCI0152",
  "NCI0155","NCI0163","NCI0165","NCI0167",
  "NCI0203",# "NCI0229",
  "NCI0245")

lapply(primaryrelapsePatients, primaryRelapse)

## Merge Files
library(readr)
fileList <- list.files("../RNASeq.Mutation.data/primary.relapse/",full.names = TRUE)
## Jun's Files 
junFile <- fileList[grepl('filtered', fileList)] %>% 
  sapply(read.csv, sep="\t") %>% 
  bind_rows
write.table(junFile, "../RNASeq.Mutation.data/all.filtered.heatmap.pimary.relapse.heat.txt", sep = "\t", row.names = F, quote = F)
## Actual File
fileList <- list.files("../RNASeq.Mutation.data/primary.relapse/",full.names = TRUE)
## Jun's Files 
actualFile <- fileList[grepl('actual', fileList)] %>% 
  sapply(read.csv, sep="\t") %>% 
  bind_rows
write.table(actualFile, "../RNASeq.Mutation.data/all.actual.heatmap.pimary.relapse.heat.txt", sep = "\t", row.names = F, quote = F)


