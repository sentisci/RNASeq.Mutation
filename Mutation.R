library(dplyr)
library(tidyr)

##Mutation Filteration
setwd("T:/Sivasish_Sindiri/R_workspace/Mutation/")
MutationData <- read.csv("Actionable_Tier1.Tier2.Tier3_10pct.v2.txt", sep="\t", header = T,
                         stringsAsFactors = FALSE)

##Extra Genes to keep
extraGenes <- read.csv("Important.Hot.Genes.txt", sep="\t", header = T, stringsAsFactors = FALSE)

##Filter1: Filter by TC, VC, VAF, 1000g, esp, ExacnonTCGA

MutationData.filter1 <- MutationData %>% 
                        filter( X1000g2014oct_all < 0.01 &
                                X1000g2014oct_eur    < 0.01 &                           
                                X1000g2014oct_afr < 0.01 &                               
                                X1000g2014oct_amr < 0.01 &                              
                                X1000g2014oct_eas  < 0.01 &                              
                                X1000g2014oct_sas < 0.01 &                              
                                esp6500_all      < 0.01 &                                
                                esp6500_ea          < 0.01 &                            
                                esp6500_aa          < 0.01 &                             
                                ExAC_ALL_nonTCGA   < 0.01 &                             
                                ExAC_AFR_nonTCGA    < 0.01 &                             
                                ExAC_AMR_nonTCGA    < 0.01 &                            
                                ExAC_EAS_nonTCGA   < 0.01 &                              
                                ExAC_FIN_nonTCGA  < 0.01 &                              
                                ExAC_NFE_nonTCGA   < 0.01 &                              
                                ExAC_OTH_nonTCGA   < 0.01 &                             
                                ExAC_SAS_nonTCGA  < 0.01 & 
                                  TotalReads >= 10 &
                                  AltReads   >= 3 &
                                  VAF        >= 0.01 &
                                  Count_154Normal <= 0) %>% 
                          filter( !grepl("Discard",PatientID))

populationFilterStr = "X1000g2014oct_all < 0.0001 & X1000g2014oct_eur < 0.0001 & X1000g2014oct_afr < 0.0001  &                               
  X1000g2014oct_amr < 0.0001 & X1000g2014oct_eas  < 0.0001 & X1000g2014oct_sas < 0.0001 & esp6500_all < 0.0001 & 
  esp6500_ea < 0.0001 & esp6500_aa < 0.0001 & ExAC_ALL_nonTCGA < 0.0001 &  ExAC_AFR_nonTCGA < 0.0001 & 
  ExAC_AMR_nonTCGA    < 0.0001 &  ExAC_EAS_nonTCGA   < 0.0001 & ExAC_FIN_nonTCGA  < 0.0001 & ExAC_NFE_nonTCGA   < 0.0001 &                              
  ExAC_OTH_nonTCGA   < 0.0001 & ExAC_SAS_nonTCGA  < 0.0001 & TotalReads >= 10 & AltReads >= 3 & VAF >= 0.0001 & Count_154Normal <= 0 "



dim(MutationData.filter1)

##Filter2: Only Tier 1 & NA
MutationData.filter2 <- MutationData.filter1 %>% filter( Level_Germline %in% c("Tier1", "NA"))
dim(MutationData.filter2)


##Filter3: Only Tier 1.2-1.4
MutationData.filter3 <- MutationData.filter1 %>% filter( Level_Somatic %in% c("Tier1.1","Tier1.2","Tier1.3","Tier1.4") &
                                                           Frac_546 <= 0.02 &
                                                          !DIAGNOSIS %in% c("RCC", "NET", "NF1") ) 
dim(MutationData.filter3)

##Filter3: Only Tier 2 & Tier 3
MutationData.filter4 <-  MutationData.filter1 %>% filter( Level_Somatic %in% c("Tier2", "Tier3") &
                                                            Gene_refGene %in% extraGenes$Important.Hot.Gene)
dim(MutationData.filter4)

### Merge all the tables
MutationData.Final <- rbind(MutationData.filter3, MutationData.filter4); dim(MutationData.Final)

### Group By Gene
MutationData.GroupBy.Gene <- MutationData.Final %>% group_by(Gene_refGene) %>% 
                                                    summarise(SampleGoup = as.character(paste0(Sample, collapse=",")), SamplesCount= n()) %>%
                                                    data.frame()
dim(MutationData.GroupBy.Gene)
                                                                                         
### Group By SNV
MutationData.Final$Level_Germline <- gsub(NA, "No Tier" ,MutationData.Final$Level_Germline)
MutationData.GroupBy.SNV <- MutationData.Final %>% group_by(X.Chr,Start, End, Ref, Alt, Gene_refGene) %>% 
                                                    summarise(SampleGoup = as.character(paste0(Sample, collapse=",")), 
                                                              SomaticTier= unique(as.character(Level_Somatic)), 
                                                              GermlineTier= unique(as.character(Level_Germline)), SamplesCount= n()) %>%
                                                    data.frame()
dim(MutationData.GroupBy.SNV)

### save the file
write.table(MutationData.filter1, "Filter+PopulationVAF+NoNormal.Raw.v2.txt", sep="\t", row.names = F, quote = F)
write.table(MutationData.filter3, "OnlyTier1.X_Frac_546LessThan_0.02.txt.v2.txt", sep="\t", row.names = F, quote = F)


write.table(MutationData.Final, "MutationData.Final.Germline.txt", sep="\t", row.names = F, quote = F)
write.table(MutationData.GroupBy.Gene, "MutationData.GroupBy.Gene.Germline.txt", sep="\t", row.names = F, quote = F)
write.table(MutationData.Final, "MutationData.GroupBy.SNV.Germline.txt", sep="\t", row.names = F, quote = F)

##########################################################################################################################
## Frequency Plot
SampleStatsMut <- StatsFinal %>% dplyr::select(Diagnosis, PatientSum, LegendPatientSum)
colnames(SampleStatsMut) <- c("Diagnosis", "PatientSum", "LegendSampleSum")
## MutationData <- read.csv("T:/Sivasish_Sindiri/R_workspace/Mutation/LandscapePaperMutation/MutationData.Final.Somatic_03_27.txt", sep="\t", header = T,stringsAsFactors = FALSE)
MutationData <- read.csv("T:/Sivasish_Sindiri/R_workspace/Mutation/LandscapePaperMutation/MutationData.Final.Somatic_06.26.18.txt", sep="\t", header = T,stringsAsFactors = FALSE) 
dim(MutationData)
MutationDataFilt <- MutationData
MutationDataFiltSelect <- MutationDataFilt %>% dplyr::select(one_of("Start","End","Ref","Alt", "Gene_refGene", "Diagnosis", "Patient")) %>% 
                              dplyr::group_by_(.dots=c("Start","End","Ref","Alt", "Gene_refGene", "Diagnosis", "Patient") ) %>% 
                              dplyr::mutate(Count=n(), Number = 1) %>% unique() %>% ungroup()  %>% 
                              dplyr::select("Gene_refGene", "Diagnosis") %>% 
                              dplyr::group_by_(.dots = c("Gene_refGene", "Diagnosis")) %>% 
                              dplyr::mutate(GeneCount = n())
View(MutationDataFiltSelect)

MutationDataFiltPercent  <-   dplyr::left_join(MutationDataFiltSelect, SampleStatsMut, by="Diagnosis" ) %>%
                              dplyr::select(one_of(c("Gene_refGene", "Diagnosis", "GeneCount", "PatientSum"))) %>%
                              dplyr::mutate(Percent=as.numeric(specify_decimal( (GeneCount/PatientSum)*100,2))) %>% unique() %>%
                              dplyr::select(Gene_refGene, Diagnosis, Percent) %>% 
                              tidyr::spread(Diagnosis, Percent) %>% data.frame()
MutationDataFiltPercent[is.na(MutationDataFiltPercent)] <- 0
View(MutationDataFiltPercent)

write.table(MutationDataFiltPercent, "./LandscapePaperMutation/MutationDataFiltPercent_06.26.18_Patient.txt", sep="\t", row.names = F, quote = F)
dim(MutationDataFiltPercent)

MutationDataGroupByPatient <- MutationDataFilt %>% 
                                dplyr::filter_(populationFilterStr) %>% 
                                dplyr::group_by_(.dots=c("Diagnosis","Sample") ) %>% dplyr::mutate(Count=n()) %>% 
                                dplyr::select(Sample, Patient, Count) %>% distinct()
View(MutationDataGroupByPatient); dim(MutationDataFiltPercent)
## For Neoantigen
write.table(MutationDataGroupByPatient, "./Neoantigen/mutation.PopulationFreq.Count.ByPatient_06.26.18.txt", sep="\t", row.names = F, quote = F)



# MutationDataFiltOnco <- MutationDataFiltPercent[,-c(1)]
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) == 0)]  <- "N"
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) > 0 & as.numeric(x) < 10)]   <- "Y1"
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) > 10 & as.numeric(x) < 20)]   <- "Y2"
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) > 20 & as.numeric(x) < 30)]   <- "Y3"
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) > 30 & as.numeric(x) < 40)]   <- "Y4"
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) > 40 & as.numeric(x) < 50)]   <- "Y5"
# MutationDataFiltOnco[apply(MutationDataFiltOnco, 2, function(x) as.numeric(x) > 50 & as.numeric(x) < 60)]   <- "Y6"
# rownames(MutationDataFiltOnco) <- MutationDataFiltPercent[,1]
# 
# order <- getOncoprintOrder(as.matrix(MutationDataFiltOnco))
# MutationDataFiltPercentOrd <- MutationDataFiltPercent[order$rowOder, order$colOder]
# View(MutationDataFiltPercentOrd)

############################################################################################################################
MutationDataFilt <- MutationData
                            # %>% dplyr::filter(!Diagnosis %in% c("NET", "NF1", "RCC")) %>% dplyr::filter( grepl("Tier1",Level_Somatic) ) 
                            # %>% dplyr::filter(ExonicFunc_refGene != -1); dim(MutationDataFiltSelect)
MutationDataFuncCbio <- MutationDataFilt %>% dplyr::select(one_of("X.Chr","Start","End","Ref","Alt", "Gene_refGene", "ExonicFunc_refGene", "AAChange_refGene", "Patient" )) %>% 
                          dplyr::group_by_(.dots=c("X.Chr","Start","End","Ref","Alt", "Gene_refGene", "ExonicFunc_refGene", "Patient" ) ) %>% 
                          unique() %>% ungroup()  %>% 
                          dplyr::group_by_(.dots=c("Gene_refGene", "ExonicFunc_refGene") ) %>%
                          dplyr::mutate(Count=n()) 
dim(MutationDataFuncCbio)
MutationDataFunc    <- MutationDataFuncCbio %>%
                          dplyr::select(one_of("Gene_refGene", "ExonicFunc_refGene")) %>% 
                          dplyr::group_by_(.dots=c("Gene_refGene", "ExonicFunc_refGene") ) %>%
                          dplyr::mutate(Count=n()) %>%  unique() %>% ungroup()  %>%
                          tidyr::spread(ExonicFunc_refGene, Count) %>% data.frame()
dim(MutationDataFunc)

MutationDataFunc[is.na(MutationDataFunc)] <- 0
View(MutationDataFunc)

write.table(MutationDataFunc, "./LandscapePaperMutation/MutationDataFunc_06.26.18_Patient.txt", sep="\t", row.names = F, quote = F)
#write.table(MutationDataFuncCbio, "./LandscapePaperMutation/MutationDataFunc_06_07_Patient.cbio.txt", sep="\t", row.names = F, quote = F)

###################################### Draw heatmap #################################
#setwd("T:/Sivasish_Sindiri/R_workspace/Mutation/LandscapePaperMutation/")
mutationFreqPer <- read.table("T:/Sivasish_Sindiri/R_workspace/Mutation/LandscapePaperMutation/MutationDataFiltFunctionPercent_06.26.18_Patient.txt", sep="\t", header = T)

#### Working ####
mutationFreqPer.RM.Single <- mutationFreqPer %>% filter(Sum >=2) %>% dplyr::select(1:14) #%>% tibble::column_to_rownames(var="GeneNames")
mutationFreqPer.RM.Single.tidy.log <- gather(mutationFreqPer.RM.Single, key="Diagnosis", value="Frequency", 2:14) %>% 
                                                  mutate(Frequency = log2(Frequency+1)/10,
                                                         GeneNames = factor(GeneNames, levels=rev(mutationFreqPer.RM.Single$GeneNames)))

HeatMap <- ggplot(mutationFreqPer.RM.Single.tidy.log, aes(Diagnosis,GeneNames)) + 
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
#  guides(fill=guide_legend(title="Fraction of patients\n altered within a histology"))


mutationType.RM.Single <- mutationFreqPer %>% filter(Sum >=2) %>% arrange(desc(Sum)) %>% dplyr::select(c(1,16:22))
mutationType.RM.Single.tidy.log <- gather(mutationType.RM.Single, key="Type", value="Percent", 2:7) %>% 
  mutate(Percent = Percent,
         GeneNames = factor(GeneNames, levels=rev(mutationFreqPer.RM.Single$GeneNames)),
         Type      = factor(Type,
                            levels = c("stopgain","splicing",
                                       "frameshift.insertion","frameshift.deletion",
                                       "nonframeshift.deletion","nonsynonymous.SNV"), ordered = TRUE)) %>% 
        dplyr::arrange(desc(GeneNames), Type)
head(mutationType.RM.Single.tidy.log)

#myColors <- setNames( c("#31a354", "#d2a679", "#669999", "#bf4040", "#810f7c", "#7094db"), 
myColors <- setNames( c("#31a354", "#d2a679", "#7094db", "#bf4040", "#810f7c", "#737373"), 
          c("frameshift.deletion","frameshift.insertion","nonframeshift.deletion","nonsynonymous.SNV","splicing","stopgain"))

# CbioColors <- setNames(c("#636363", "#636363", "#994d00", "#31a354", "#636363", "#636363" ), 
#                              c("frameshift.deletion","frameshift.insertion","nonframeshift.deletion","nonsynonymous.SNV","splicing","stopgain"))

customColors <- myColors

## green orange/brown violet red brinjal skyblue
StackedBar <- ggplot(mutationType.RM.Single.tidy.log, aes(x=GeneNames, y=Percent, fill=Type) ) +
#ggplot(mutationType.RM.Single.tidy.log, aes(x=GeneNames, y=Percent, fill=Type) ) +
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
  
#grid.arrange(HeatMap, StackedBar, ncol=2)
pdf("./MutationDataFuncPercent_06.26.18_Patient.pdf", height = 11, width = 15)
ggarrange(HeatMap, StackedBar, ncol=2)
dev.off()

#### Teseting ####
mutationFreqPer.RM.Single <- mutationFreqPer %>% dplyr::filter(Sum >=2) %>% dplyr::select(1:14) %>% tibble::column_to_rownames(var="GeneNames")
mutationFreqPer.RM.Single.log <- log2(mutationFreqPer.RM.Single+1)
pdf("./Mutation.PlotCbio.pdf", height = 10, width = 15)
superheat(mutationFreqPer.RM.Single,  
          title = "Mutation Plot",
          #linkage.method = "ward.D2",
          #legend=FALSE,
          grid.hline = FALSE,
          grid.vline = FALSE,
          # grid.hline.size = 0.01,
          # grid.vline.size = 0.01,
          scale = TRUE,
          heat.pal = c("white", "brown"),
          heat.pal.values = c(0, 5),
          bottom.label.text.angle=90,
          title.size = 6)
dev.off()



colfunc <- colorRampPalette(c( "white", "#de2d26"))
pheatmap(mutationFreqPer.RM.Single.log,
         cluster_cols = FALSE, 
         cluster_rows = FALSE, 
         clustering_distance_cols = "euclidean" ,
         clustering_method = "mcquitty",
         color = colfunc(30), 
         cellwidth = 20,
         #cellheight = 9,
         main = "Mutation Plot mcquitty")


















