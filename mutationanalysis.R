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
  metaDataFileName        = "MetadataMapper.edited.txt",
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
mergeVCFObjectDiagnosis <- dplyr::left_join(mergeVCFObject, data.table(rnaseqMutationProject$validMetaDataDF[, 
                                                        c("Patient.ID","DIAGNOSIS.Alias","LIBRARY_TYPE")]), by="Patient.ID") %>% 
                           dplyr::filter(! Patient.ID %in% c("NA.")) %>% 
                           dplyr::filter(!is.na(MAF)) %>% 
                           data.table()
dim(mergeVCFObjectDiagnosis)
## dplyr way
#mergeVCFObjectDiagnosis <- mergeVCFObjectDiagnosis %>% dplyr::mutate(variantKey = paste(Chr,Start,End,Ref,Alt))
## data.table way
mergeVCFObject.variantKey <- mergeVCFObjectDiagnosis[, variantKey := paste(Chr,Start,End,Ref,Alt,sep="_")][order(variantKey)]
dim(mergeVCFObject.variantKey)
write.table(mergeVCFObject.variantKey, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "1.MAF.LT.5pc.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

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
write.table(TumorFiltered.Normal, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, "2.MAF.LT.5pc_No.NS.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

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
                                              Total.coverage >= 10 & Variant.coverage >= 3 & VAF >= 0.10 & MAF <= 0.01 & prop <= 0.1 ]
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF)
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "3.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

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
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4a.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.Tier1.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Somatic.Tier %in% c("Tier 2")]
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier2, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4b.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.Tier2.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Somatic.Tier %in% c("Tier 3")]
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier3, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4c.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.Tier3.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[ Somatic.Tier %in% c("Tier 4")]
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier4, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4d.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.Tier4.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)


## Only disease causing
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF[
  Clinvar.Pathogenic %in% c("Y") | HGMD.Disease %in% c("Y") | Intervar %in% c("Likely pathogenic",
                                                                              "Pathogenic",
                                                                              "Uncertain significance") ]
write.table(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.DiseaseCausing, 
            paste(rnaseqMutationProject$workDir, rnaseqMutationProject$projectName, rnaseqMutationProject$outputdirTXTDir, 
                  "4e.MAF.LT.5pc_No.NS_TC.GTE.10_VC.GTE.3_VAF.GTE.10pc_MAF.LTE.1pc_prop.LTE.10.DiseaseCausing.txt", sep="/"),
            sep="\t", quote = FALSE, row.names = FALSE, col.names = TRUE)



### Plotting
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter(LIBRARY_TYPE %in% c("Tumor"))
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor)
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor.MT.2 <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF %>% filter(N >=2 )
TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.CL <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1 %>% filter(LIBRARY_TYPE %in% c("CellLine"))
dim(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.CL)

## Sample statistics
SampleStatsMut.Tumor <- rnaseqMutationProject$validMetaDataDF %>% dplyr::filter(LIBRARY_TYPE %in% c("Tumor")) %>% dplyr::select(DIAGNOSIS.Alias) %>% 
                              dplyr::group_by(DIAGNOSIS.Alias) %>% 
                              dplyr::mutate(PatientSum=n()) %>% distinct()
SampleStatsMut.Tumor

### Diagnosis Wise mutation frequency 
MutationDataFilt <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor
MutationDataFiltSelect <- MutationDataFilt %>% dplyr::select(one_of("Chr","Start","End","Ref","Alt", "Gene", "DIAGNOSIS.Alias", "Patient.ID")) %>% 
  dplyr::group_by_(.dots=c("Chr", "Start","End","Ref","Alt", "Gene", "DIAGNOSIS.Alias", "Patient.ID") ) %>% 
  dplyr::mutate(Count=n(), Number = 1) %>% unique() %>% ungroup()  %>% 
  dplyr::select("Gene","DIAGNOSIS.Alias") %>% 
  dplyr::group_by_(.dots = c("Gene", "DIAGNOSIS.Alias")) %>% 
  dplyr::mutate(GeneCount = n()) %>% distinct()
View(MutationDataFiltSelect); dim(MutationDataFiltSelect)

MutationDataFiltPercent  <-   dplyr::left_join(MutationDataFiltSelect, SampleStatsMut.Tumor, by="DIAGNOSIS.Alias" ) %>%
  dplyr::select(one_of(c("Gene", "DIAGNOSIS.Alias", "GeneCount", "PatientSum"))) %>%
  dplyr::mutate(Percent=( (GeneCount/PatientSum)*100)) %>% unique() %>%
  dplyr::select(Gene, DIAGNOSIS.Alias, Percent) %>% 
  tidyr::spread(DIAGNOSIS.Alias, Percent) %>% data.frame()
MutationDataFiltPercent[is.na(MutationDataFiltPercent)] <- 0
View(MutationDataFiltPercent);dim(MutationDataFiltPercent)

#### By Gene mutation function frequency

MutationDataFilt <- TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor
MutationDataFilt[which(TumorFiltered.Normal.freq.VAF.TC.VC.MAF.Tier1.Tumor$Exonic.function %in% c("ALDH2", "FANCD2", "ERBB3")), c("Exonic.function")] <- "nonsynonymous SNV"
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

mutationFreqPer.RM.Single <- mutationFreqPer %>% filter(Sum >=2) %>% dplyr::select(1:14) #%>% tibble::column_to_rownames(var="Gene")
mutationFreqPer.RM.Single.tidy.log <- gather(mutationFreqPer.RM.Single, key="Diagnosis", value="Frequency", 2:14) %>% 
  mutate(Frequency = log2(Frequency+1)/10,
         Gene = factor(Gene, levels=rev(mutationFreqPer.RM.Single$Gene)))

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



