library(R6)
library(assertthat)
library(microbenchmark)

# ### Project SetUp Class ####
ProjectSetUp <- R6Class(
  classname = "ProjectSetUP",
  portable = TRUE,
  private  = list(
    
    checkDirExists = function() {
      
      projectDirPath = paste(self$workDir, self$projectName, sep="/")
      ## Check & create if the project directory exists ?? 
      if(!dir.exists(file.path(self$workDir))) {
        warning(paste(self$workDir, "No Such directory exsist !!"))
        paste0("Creating working dir for this project ", projectDirPath)
        dir.create(self$workDir)
      }
      
      ## Create all the required dirs for the project.
      dir.create(projectDirPath)
      self$fileDirs <- paste(projectDirPath, c(self$outputdirRDSDir, self$outputdirTXTDir, self$vcfFilesDir, 
                                               self$plotsDir, self$plotsDataDir), sep="/")
      lapply(self$fileDirs, function(x){ if(!dir.exists(x)) dir.create(x) })
    },
    ## Generate filters to exclude given list
    generateFiltersToExclude = function(factorsToExclude=NA){
      
      lapply(factorsToExclude, function(x){
        subList=x
        sublistNames=names(x)
        paste0("!(",
               paste0(sapply(sublistNames, function(y){
                 
                 #print(subList[y])
                 filter <- paste0(y," %in% ","\"",unlist(unname(subList[y])), "\"" )
                 return(filter)
               }), collapse = " & " ),
               ")"
        )
      })
    },
    readMetaData = function() {
      self$metaDataDF <- read.csv(paste0(self$workDir,"/",self$projectName, "/",self$metaDataFileName), sep="\t", header = T);
      print(paste0("Dimension of metadata is ", paste(dim(self$metaDataDF))))
      if(!is.na(self$factorsToExclude)){
        filters <- private$generateFiltersToExclude(factorsToExclude=self$factorsToExclude)
        test <- sapply(filters, function(x){
          self$metaDataDF <- dplyr::filter_(self$metaDataDF, .dots=list(x))
        })
      }
      print(paste0("Dimension of metadata after applying parameter 'factorsToExclude' is ", paste(dim(self$metaDataDF)[1])))
      print("Make R valid names in the metadata file and storing it as validMetaDataDF")
      tempDF <- as.data.frame(lapply(self$metaDataDF[,c("Patient.ID","Sample.Data.ID","Sample.ID","Sample.ID.Alias")], make.names)) %>% data.frame()
      self$validMetaDataDF <- cbind(tempDF, self$metaDataDF[,which(names(self$metaDataDF) != c("Patient.ID","Sample.Data.ID","Sample.ID","Sample.ID.Alias"))])
    },
    readAnnotation = function() {
      self$annotationDF <- readRDS(self$annotationRDS) %>% as.data.frame()
    },
    readProteinCoding = function(){
      self$pcDF <- as.data.frame( readRDS(self$pcRDS) )
      colnames(self$pcDF)[c(1,2)] <- c("GeneName","GeneID")
    }
    ,
    readTranscriptionFactor = function(){
      self$tfDF <- readRDS(self$tfRDS)
    },
    readCellSurface = function(){
      self$csDF <- readRDS(self$csRDS)  %>% dplyr::filter(NewCount >= 5)
    },
    readCancerGermlineAntigens = function(){
      self$cgaDF <- readRDS(self$cgaRDS) 
    },
    readpax3Foxo1Targets = function(){
      self$pax3Foxo1DF <- readRDS(self$pax3Foxo1RDS) 
    },
    readewsr1Fli1Targets = function(){
      self$ewsr1Fli1DF <- readRDS(self$ewsr1Fli1RDS) 
    }
  ),
  public = list(
    date                    = NULL,
    time                    = NULL,
    projectName             = NULL,
    annotationRDS           = NULL,
    outputPrefix            = NULL,
    factorName              = NULL,
    metaDataFileName        = NULL,
    workDir                 = NULL,
    outputdirRDSDir         = NULL,
    outputdirTXTDir         = NULL,
    vcfFilesDir             = NULL,
    plotsDir                = NULL,
    plotsDataDir            = NULL,
    metaDataDF              = NULL,
    annotationDF            = NULL,
    pcDF                    = NULL,
    pcRDS                   = NULL,
    tfDF                    = NULL,
    tfRDS                   = NULL,
    csDF                    = NULL, 
    csRDS                   = NULL,
    cgaDF                   = NULL,
    cgaRDS                  = NULL,
    pax3Foxo1DF             = NULL,
    pax3Foxo1RDS            = NULL,
    ewsr1Fli1DF             = NULL,
    ewsr1Fli1RDS            = NULL,
    
    factorsToExclude        = NULL,
    fileDirs                = NULL, 
    validMetaDataDF         = NULL,
    initialize              = function(date = NA, time =NA, projectName = NA, annotationRDS = NA, outputPrefix = NA,
                                       factorName = NA, metaDataFileName = NA, workDir = NA, outputdirRDSDir = NA, 
                                       outputdirTXTDir = NA, vcfFilesDir = NA, plotsDir = NA, plotsDataDir = NA, pcRDS = NA,tfRDS=NA,
                                       csRDS =NA, cgaRDS=NA, pax3Foxo1RDS=NA, ewsr1Fli1RDS=NA, factorsToExclude=NA) {
      
      self$date <- date
      self$time <- time
      self$projectName <- projectName
      self$annotationRDS <- annotationRDS
      self$outputPrefix <- outputPrefix
      self$factorName <- factorName
      self$metaDataFileName <- metaDataFileName
      self$workDir <- "C:/Users/sindiris/R Scribble/"
      self$outputdirRDSDir <- outputdirRDSDir
      self$outputdirTXTDir <- outputdirTXTDir
      self$vcfFilesDir        <- vcfFilesDir
      self$plotsDir <- plotsDir
      self$plotsDataDir <- plotsDataDir
      
      self$pcRDS <- pcRDS
      self$tfRDS <- tfRDS
      self$csRDS <- csRDS
      self$cgaRDS <- cgaRDS
      
      self$pax3Foxo1RDS <- pax3Foxo1RDS
      self$ewsr1Fli1RDS <- ewsr1Fli1RDS
      
      self$factorsToExclude <- factorsToExclude
      private$checkDirExists()
      private$readMetaData()
      private$readAnnotation()
      if (!is.na(pcRDS)){ private$readProteinCoding() }
      if (!is.na(tfRDS)){ private$readTranscriptionFactor() }
      if (!is.na(csRDS)){ private$readCellSurface() }
      if (!is.na(csRDS)){ private$readCancerGermlineAntigens() }
      
      if (!is.na(pax3Foxo1RDS)){ private$readpax3Foxo1Targets() }
      if (!is.na(ewsr1Fli1RDS)){ private$readewsr1Fli1Targets() }
      
      print(paste0("Changing to working directory ", self$workDir))
    }
  ) 
)

# ### Utilities Class ####

CoreUtilities <- R6Class(
  classname = "CoreUtilities",
  portable  = TRUE,
  ## Read CSV or TXT files
  private   = list(
    ## Read TXT Files
    readTXTFiles  = function(x, fileSuffix=NA ){ ##, colNameSelect=NA, primaryID=NA ){
      if(!is.na(fileSuffix)) 
      {
        sampleName <- gsub(fileSuffix, "", basename(x)) 
      } else { 
        sampleName <- basename(x) 
      }
      
      rdsObj <- readRDS("C:/Users/sindiris/R Scribble/RNASeq.Mutation.data/dummyRDS")
      
      tryCatch({ 
        rdsObj <-     fread(x, sep="\t", header = TRUE, check.names = TRUE) 
      }, error=function(e){
        print(paste0("I am in readTXTFiles ", e, " ", sampleName))
        rdsObj$Patient.ID <- sampleName
      })
      return(rdsObj)
      
    },
    ## Merge CSV or TXT files
    mergeTXTFiles = function( x, fileSuffix=NA ){
      selectedFileList       <- x[which(basename(x) %in% self$allFileList)]
      notselectedFileList    <- x[which(!basename(x) %in% self$allFileList)]
      if( length(notselectedFileList) > 1 ) { stop(paste(" couldn't retrived following files. Please check the file names or suffix "
                                                                     ,paste(basename(notselectedFileList), collapse = "\n")))}
      
      print(paste0("Selecting ", length(selectedFileList), " files out of ", length(x), " from the given folder"))
      dataMatrixLists     <- lapply(selectedFileList, private$readTXTFiles, fileSuffix=fileSuffix )
      dataMatrix          <- data.table::rbindlist(dataMatrixLists, idcol = FALSE)
      return(dataMatrix)
    },
    ## Read RDS files
    readRDSFiles  = function(x, fileFormat = NA, colInterest = NA, isRowNames= FALSE, rowNamesColInFile = NA ){
      
      print(paste(x))
      featureCountRDS <- readRDS(x)
      countDF         <- featureCountRDS$count %>% data.frame()
      return(countDF)
      
    },
    ## Merge RDS files
    mergeRDSFiles = function(x ){
      # dataMatrixLists     <- lapply(unlist(x), private$readRDSFiles)
      # rowNames            <- rownames( dataMatrixLists[[1]] )
      # dataMatrix          <- rbindlist(dataMatrixLists, use.names = TRUE)[,as.matrix(.SD)]
      # rownames(dataMatrix)<- rowNames
      # return(as.data.frame(dataMatrix))
      
      dataMatrix          <- do.call(cbind,lapply(unlist(x),private$readRDSFiles))
      return(dataMatrix)
      
    }
  ),
  public    = list(
    workDir                = NULL,
    projectName            = NULL,
    allFileList            = NULL,
    ProjectSetUpObject     = NULL,
    initialize             = function(ProjectSetUpObject = NA ){
      
      assert_that("ProjectSetUP" %in% class(ProjectSetUpObject), 
                  msg="Please setup Project using ProjectSetUp class !!\nProjectSetUpObject cannot be NA !!")
      self$workDir     <- ProjectSetUpObject$workDir
      self$projectName <- ProjectSetUpObject$projectName
      self$ProjectSetUpObject <- ProjectSetUpObject
    },
    ## Get merged matrix
    getMergedMatrix = function(dir = NA, fileFormat = NA, 
                               ## colNameSelect = NA, colIndexSelect = NA, isRowNames = FALSE, rowNamesColInFile = NA,
                               fileSuffix=NA, 
                               ## primaryID=NA, 
                               metadata=NA, metadataFileRefCol=NA){
      
      # if( !is.na(colNameSelect) & !is.na(colIndexSelect)) stopifnot("Use either colNameSelect or colIndexSelect but not both")
      # if( is.na(colNameSelect) & is.na(colIndexSelect)) stopifnot("Use either colNameSelect or colIndexSelect but not both")
      # assert_that( !is.na(fileFormat)  , msg ="Please provide file format to parse. Following are supported, \"rds\", \"csv or text\", \"feather\" ")
      #
      # if( is.na(colNameSelect) & is.na(colIndexSelect)) {
      #   print("First column of every file will be selected.Or use either \"colNameSelect\" or \"colIndexSelect\"")
      #   colIndexSelect = 1
      # }
      # if ( is.na(colNameSelect) ) { 
      #   colNameSelect  =  colIndexSelect
      # } 
      # if(isRowNames) assert_that(!is.na(rowNamesColInFile), msg= "Please provide row names index , \"rowNamesColInFile\" can't be NA ")
      
      self$allFileList = paste0(as.character(metadata[,metadataFileRefCol]),fileSuffix)
      print(paste0("Total unique entries in the metadata's reference column are ",length(unique(self$allFileList))))
      
      dirs                <- list.dirs(paste0(self$workDir,"/",self$projectName,"/",dir))[-1]
      folderNames         <- basename(dirs)
      
      switch(fileFormat,
             
             "rds"     = {
               mergedDataList       <- lapply( lapply(dirs, list.files, full.names=T), private$mergeRDSFiles )
               rowNames              <- rownames(mergedDataList[[1]])
               mergedData            <- dplyr::bind_cols( mergedDataList)
               rownames(mergedData)  <- rowNames
               return(mergedData)
             },
             "txt"     = {
               mergedDataList        <- lapply( lapply(dirs, list.files, full.names=T),
                                                private$mergeTXTFiles, fileSuffix = fileSuffix ) 
               mergedData            <- data.table::rbindlist(mergedDataList, idcol = TRUE)
               colnames(mergedData)  <- gsub("\\s+", ".", colnames(mergedData))
               return(mergedData)
             }
      )
    },
    ## matchAndChangeColNames columns
    leftjoinDTs = function( DT1=NA, key=NA,  DT2=NA ){
      metaDT <- merge(DT1, DT2, by=key, all.x=TRUE) 
      return(setDT(metaDT))
    },
    ## Factorise columns
    factorizeColumn = function( toFactor, asFactor){
      factorColumn <- factor(toFactor, levels=asFactor, ordered = TRUE )
      return(factorColumn)
    },
    ## subsetMetaData
    subsetMetaData = function(colnamesDF=NA){
      df <- dplyr::left_join(colnamesDF,self$ProjectSetUpObject$metaDataDF, by="SAMPLE_ID") %>% filter(complete.cases(.))
      if(!is.null(self$ProjectSetUpObject$metaDataDF)) {
        self$ProjectSetUpObject$metaDataDF <- df
      }
    },
    ## Add gene annotation and filter genes
    filterAnotateGenes = function(filterName= NA , folderName = "data", queryDF = NA, rdsDir=NA, txtDir=NA){
      
      print(paste0( rdsDir ,"/", folderName, "/" ))
      rdsDirPath  <- paste0( rdsDir ,"/", folderName, "/" )
      txtDirPath  <- paste0( txtDir ,"/", folderName, "/" )
      
      if(!dir.exists(rdsDirPath)) { dir.create(rdsDirPath) }
      if(!dir.exists(txtDirPath)) { dir.create(txtDirPath) }
      
      rdsfile <- paste0(rdsDirPath, paste0(folderName,"_",filterName,".rds") )
      txtFile <- paste0(txtDirPath, paste0(folderName,"_",filterName,".txt") )
      
      switch(filterName,
             
             "proteinCoding"= {
               #print("Before filtering for pritein coding genes")
               GeneDF_DiffExp <- na.omit(dplyr::left_join(self$ProjectSetUpObject$pcDF, queryDF, by="GeneID"))
               #private$GeneDF_DiffExp <- dplyr::left_join( private$GeneDF_DiffExp, rnaseqProject$csDF, by="GeneName")
               #print("Filter matched ", dim(GeneDF_DiffExp_PC)[1], " out of  ", dim(rnaseqProject$pcDF)[1], " given protein coding genes")
               #pcDF <- private$GeneDF_DiffExp %>% filter(GeneName %in% rnaseqProject$pcDF)
             },
             "cellsurface"= {
               #print("filtering for cellsurface genes")
               GeneDF_DiffExp <- queryDF %>% filter(GeneName %in% as.character(self$ProjectSetUpObject$csDF[,"GeneName"]))
               GeneDF_DiffExp <- dplyr::left_join(GeneDF_DiffExp, self$ProjectSetUpObject$pcDF, by = "GeneID")
               #print("Filter matched ", dim(GeneDF_DiffExp_csDF)[1], " out of  ", dim(rnaseqProject$csDF)[1], " given CS genes")
             },
             "transcriptionFactor"= {
               # print("filtering for transcriptionFactor genes")
               GeneDF_DiffExp <- queryDF %>% filter(GeneName %in% as.character(self$ProjectSetUpObject$tfDF[,"GeneName"]))
               GeneDF_DiffExp <- dplyr::left_join(GeneDF_DiffExp, self$ProjectSetUpObject$pcDF, by = "GeneID")
               #print("Filter matched ", dim(GeneDF_DiffExp_tfDF)[1], " out of  ", dim(rnaseqProject$tfDF)[1], " given TF genes")
             },
             "cancergermlineantigen"= {
               # print("filtering for cancergermlineantigen genes")
               GeneDF_DiffExp <- queryDF %>% filter(GeneName %in% as.character(self$ProjectSetUpObject$cgaDF[,"GeneName"]))
               GeneDF_DiffExp <- dplyr::left_join(GeneDF_DiffExp, self$ProjectSetUpObject$pcDF, by = "GeneID")
               #print("Filter matched ", dim(GeneDF_DiffExp_cgaDF)[1], " out of  ", dim(rnaseqProject$tfDF)[1], " given TF genes")
             },
             "all" = {
               GeneDF_DiffExp <- queryDF
               GeneDF_DiffExp <- dplyr::left_join(GeneDF_DiffExp, self$ProjectSetUpObject$pcDF,  by="GeneID")
               GeneDF_DiffExp <- GeneDF_DiffExp %>% mutate(
                 meanBrainExp  = ifelse(GeneName.x %in% self$ProjectSetUpObject$BrainExpDF[,"GeneName"], self$ProjectSetUpObject$BrainExpDF[,"MeanExp"], "N"),
                 meanHeartExp  = ifelse(GeneName.x %in% self$ProjectSetUpObject$HeartExpDF[,"GeneName"], self$ProjectSetUpObject$HeartExpDF[,"MeanExp"], "N"),
                 meanKidneyExp = ifelse(GeneName.x %in% self$ProjectSetUpObject$KidneyExpDF[,"GeneName"], self$ProjectSetUpObject$KidneyExpDF[,"MeanExp"], "N"),
                 meanLiverExp  = ifelse(GeneName.x %in% self$ProjectSetUpObject$LiverExpDF[,"GeneName"], self$ProjectSetUpObject$LiverExpDF[,"MeanExp"], "N"),
                 meanLungExp   = ifelse(GeneName.x %in% self$ProjectSetUpObject$LungExpDF[,"GeneName"], self$ProjectSetUpObject$LungExpDF[,"MeanExp"], "N"),
                 
                 CellSurface   = ifelse(GeneName.x %in% self$ProjectSetUpObject$csDF[,"GeneName"], "Y", "N"),
                 TranscriptionFactor     = ifelse(GeneName.x %in% self$ProjectSetUpObject$tfDF[,"GeneName"], "Y", "N"),
                 CancerGermlineAntigen   = ifelse(GeneName.x %in% self$ProjectSetUpObject$cgaDF[,"GeneName"], "Y", "N"),
                 PAX3FOXO1     = ifelse(GeneName.x %in% self$ProjectSetUpObject$pax3Foxo1DF[,"GeneName"], "Y", "N"),
                 EWSR1FL1      = ifelse(GeneName.x %in% self$ProjectSetUpObject$ewsr1Fli1DF[,"GeneName"], "Y", "N") )
             }
      )
      saveRDS(GeneDF_DiffExp, rdsfile)
      write.table(x = GeneDF_DiffExp, file = txtFile, sep="\t", row.names = FALSE, quote=FALSE)
      return(GeneDF_DiffExp)
    }  )
)

