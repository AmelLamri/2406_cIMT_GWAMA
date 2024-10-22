library(readxl)
library(R.utils)
library(data.table)

m_p<- "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
o_p <- paste0(m_p, "/Results")
gw_p <- paste0(o_p, "/GWAS") 
ma_p <- paste0(o_p, "/GWAMA") 

pops2 <- c("EUR", "ME")
imt <- c("maxLCFp", "meanLCFp", "minLCFp")
adj <- c("adjBMI", "noBMI")




d <- as.data.frame(read_excel(paste0(o_p,"/FilesTables/241017_Files.xlsx" )))
stds <- unique(d$Study)


#Fill in empty raw file paths 
  table(d$pop, useNA="ifany")
  table(is.na(d$RawFilePaths), useNA="ifany")
  table(d[is.na(d$RawFilePaths), "Study"] )
  table(is.na(d$Study) )
  table(d$AgeGroup, useNA="ifany")

  d[is.na(d$RawFilePaths), "RawFilePaths"]  <- paste0("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/", d[is.na(d$RawFilePaths), "Study"], "/Raw")

  d[d$Study=="ALSPAC", "RawFilePaths"] <- paste0("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/", d[d$Study=="ALSPAC", "Study"], "/Raw")



# Clean column names of the table 
  colnames(d) <- gsub(" ", "_", colnames(d))
# Correct AgeGroup Var

  table(d$AgeGroup, useNA="ifany")
  d[d$AgeGroup=="CLD", "AgeGroup"] <- "ADO"
  d[d$AgeGroup=="xAGE", "AgeGroup"] <- "AADxAGE"


ages <- c("OAD", "ADL", "YAD", "AAD", "ADO", "AADxAGE")


# Correct some file names 
  d[d$Study=="FAMILY", "File"] <- gsub("ALL", "MW", d[d$Study=="FAMILY", "File"])
  d[d$Study=="FAMILY", "File"] <- gsub("WOMEN", "W", d[d$Study=="FAMILY", "File"])
  d[d$Study=="FAMILY", "File"] <- gsub("MEN", "M", d[d$Study=="FAMILY", "File"])
 
  d[d$Study=="ALSPAC", "File"] <- gsub("minIMT", "minLCF", d[d$Study=="ALSPAC", "File"])
  d[d$Study=="ALSPAC", "File"] <- gsub("meanIMT", "meanLCF", d[d$Study=="ALSPAC", "File"])
  d[d$Study=="ALSPAC", "File"] <- gsub("maxIMT", "maxLCF", d[d$Study=="ALSPAC", "File"])
  d[d$Study == "UKB", "File"] <- gsub("240920", "240927", d[d$Study == "UKB", "File"])
 
  d[which(d$Study=="UKB" & d$pop=="NbEUR"),"File"] <- gsub("NbEUR", "NBEUR", d[which(d$Study=="UKB" & d$pop=="NbEUR"),"File"] )

# Correct some file paths 
  d[d$Study=="FAMILY", "RawFilePaths"] <- paste0("/genetics/Nutrigen/FAMILY/2301_FAMILY_cIMT_GWAS_PHRI/Outputs/",d[d$Study=="FAMILY", "AgeGroup"])

# Create Cleaned Raw files variable 

  d$RawFilePaths_m <- d$RawFilePaths
  d[d$Study=="CHCP", "RawFilePaths_m"] <- paste0(d[d$Study=="CHCP", "RawFilePaths"], "/Modified/")
  d[d$Study=="ABCD", "RawFilePaths_m"] <- paste0(d[d$Study=="ABCD", "RawFilePaths"], "/Modified/")
  d[d$Study=="SWS", "RawFilePaths_m"] <- paste0(d[d$Study=="SWS", "RawFilePaths"], "/Modified/")
  #d[d$Study=="ALSPAC", "RawFilePaths_m"] <- paste0(d[d$Study=="ALSPAC", "RawFilePaths"], "")


# Makesure all RawFilePaths exist
  for (i in 1:nrow(d)){
    stdi  <- unique(d[i, "Study"])
    pheni <- unique(d[i, "cIMT_measure"])
    if(stdi=="CHCP" & pheni %in% c("minRCN","meanRCN", "minRCF" , "meanRCF", "maxRCN", "maxRCF")) next()
  
    if(!dir.exists(d$RawFilePaths_m[i])) cat (d$RawFilePaths_m[i],  " path does not exist \n")
  
  }


# Delete non existing files 
  #d <- d[-which(d$Study == "ALSPAC" & d$AgeGroup %in% c("OAD", "AAD", "ADL") & d$sexGroup %in% c("M")), ] 
  # Alpac daults are only moms  
  #d <- d[-which(d$Study == "FAMILY" & d$AgeGroup %in% c("OAD", "AAD", "ADL") & d$sexGroup %in% c("M")), ] 




# Correct some file names to use the updated version 


  d$File2 <- d$File 
  d[d$Study == "CHCP", "File2"] <- gsub(".zip", "_3.txt.gz", d[d$Study == "CHCP", "File"])
  d[d$Study == "GRACE", "File2"] <- gsub("20240903.txt.gz", "20240903_tmp.txt.gz", d[d$Study == "GRACE", "File"])
  d[d$Study == "FAMILY", "File2"] <- gsub("20240825.txt.gz", "20240825_2.txt.gz", d[d$Study == "FAMILY", "File"])
  d[d$Study == "SWS", "File2"] <- gsub("240130", "240130_2", d[d$Study == "SWS", "File"])
  d[d$Study == "ALSPAC", "File2"] <- gsub("\\.20240702\\.", ".\\20240702_3\\.", d[d$Study == "ALSPAC", "File"])

# Correct Population variables 

  #table(d$pop, useNA="ifany") 
  # AFR    AMR    EUR nonEUR    SAS
  #  18     18    440    132     18
  d[grep("NbEUR", d$File), "pop"] <- "NbEUR"
  d[grep("BRIT", d$File), "pop"] <- "BRIT"
  table(d$pop, useNA="ifany") 
  # AFR    AMR   BRIT    EUR  NbEUR nonEUR    SAS
  #  18     18     18    404     18    132     181

  d[d$pop %in% c("BRIT", "EUR", "NbEUR"), "pop2"] <- "EUR"
  d[d$pop %in% c("AFR", "AMR", "nonEUR", "SAS"), "pop2"] <- "nonEUR"
  table(d$pop2, useNA="ifany") 
  # EUR nonEUR
  # 422    186



# Correct some phenos Make sure that the pheno is correct 
  d[d$cIMT_measure== "MeanLCF", "cIMT_measure"] <- "meanLCF"
  d[d$cIMT_measure== "min240", "cIMT_measure"] <- "minIMT240"
  d[d$cIMT_measure== "mean240", "cIMT_measure"] <- "meanIMT240"
  d[d$cIMT_measure== "max240", "cIMT_measure"] <- "maxIMT240"

  for (i in 1:nrow(d)){
    if(length(grep(d$cIMT_measure[i], d$File[i]) ) != 1) cat (i, "\n")
    } # all good 

# Check that sex group matches 
  table(d$sexGroup, useNA="ifany")
  #  M  MW   W
  # 180 248 210
  mw <- d[d$sexGroup=="MW",]

  for (i in 1:nrow(mw)){
    if(length(grep("xAGE.", mw$File[i]) ) == 0 & length(grep(".MW.", mw$File[i]) ) == 0 ) cat (i, "\n")
  } # ALL GOOD 
  m <- d[d$sexGroup=="M",]
  for (i in 1:nrow(m)){
    if(length(grep(".M.", m[i,]) ) == 0 ) cat (i, "\n")
    } 

  w <- d[d$sexGroup=="W",]
  for (i in 1:nrow(w)){
    if(length(grep(".W.", w[i,]) ) == 0 ) cat (i, "\n")
    } 

# Make sure that the pop is correct 
  for (i in 1:nrow(d)){
    if(length(grep( d[i,"pop"], d[i,"File"]) ) == 0 ) cat (i, "\n")
    } #all good 



# Create variable to see if raw- modified files exist 
  d$File_Exists <- file.exists(paste0(d$RawFilePaths_m, "/",d$File2 ))
  table(d$File_Exists)

  #FALSE  TRUE
  #  414   236

  # Makesure sure all modified Files exist
  table(duplicated(d$File))
  #FALSE 650

  # unique(d[!d$File_Exists & d$AgeGroup=="ADO", c("Study", "pop", "AgeGroup", "cIMT_measure", "RawFilePaths_m", "File2")])

  unique(d[!d$File_Exists & d$AgeGroup=="ADO" &  d$cIMT_measure=="maxLCF" , c("Study", "pop", "AgeGroup", "cIMT_measure","RawFilePaths_m", "File2")])
# Create a variable for QCed files
  d$QCed_fn <- paste0("CLEANED.", gsub(".txt.gz", "", d$File2),".cpaid.gz")
  d$QCed_fn <- gsub(".gz.cpaid.gz", ".cpaid.gz",d$QCed_fn  ) 
  d[d$Study=="UKB", "QCed_fn"  ] <- gsub("240920", "240927", d[d$Study=="UKB", "QCed_fn"  ] ) 
  d[d$Study=="GenR2", "QCed_fn"]  <- gsub(".MetaScore.assoc.cpaid.gz", ".cpaid.MetaScore.assoc.gz", d[d$Study=="GenR2", "QCed_fn"]) 

  adjsd <- d$AdjBMI
  table(adjsd)
  d[d$AdjBMI=="T","AdjBMI2"] <- "adjBMI"
  d[d$AdjBMI=="F","AdjBMI2"] <- "noBMI"
  table(d$AdjBMI)
  table(d$AdjBMI2)
  d$AdjBMI <- as.logical(d$AdjBMI)

   d$MetaA_p <- paste0(ma_p, "/", d$AgeGroup, "/", d$cIMT_measure, "/", d$pop2, "/", d$sexGroup, "/", d$AdjBMI2)


#Create a variable for QCed file paths
  d$QCed_p <- paste0(gw_p,"/",d$Study, "/QCed","/", d$AgeGroup,"/", d$cIMT_measure,"/",d$pop )
  
  d[d$Study=="UKB" & d$cIMT_measure == "maxIMT240", "QCed_p"] <- gsub("maxIMT240", "maxIMT240_2",d[d$Study=="UKB" & d$cIMT_measure == "maxIMT240", "QCed_p"])
  d[d$Study=="UKB" & d$pop == "NbEUR", "QCed_p"] <- gsub("NbEUR", "NBEUR",d[d$Study=="UKB" & d$pop == "NbEUR", "QCed_p"])


# See if QCed files exist 

for (i in 1:nrow(d)){
  if (!dir.exists(d$QCed_p[i])){
    dir.create(d$QCed_p[i], rec=T)
    }
  }


for (i in 1:nrow(d)){
  qcf <- paste0(d$QCed_p[i] ,"/",  d$QCed_fn[i] )
  if( file.exists(qcf)) {
    d$QCed_File_Exists[i] <- T
    } 
  if(!file.exists(qcf)) {
    d$QCed_File_Exists[i] <- F
    }
}

table(d$QCed_File_Exists)
 # FALSE  TRUE
 # 548   102




  #d[d$QCed_File_Exists & d$AgeGroup =="OAD" & d$cIMT_measure =="maxIMT240" & d$Study=="UKB",c("QCed_File_Exists" ,"pop","sexGroup")]

 # !d$QCed_File_Exists 
 # d$Study=="GenR2"
 # d$AdjBMI2=="adjBMI"
 # d$AgeGroup
 # d$pop2 =="EUR"
 # d$cIMT_measure =="maxLCF"
 # d$sexGroup%in% c("M","W")
 # c("Study", "pop", "AgeGroup", "cIMT_measure", "sexGroup", "QCed_fn", "QCed_p", "pop2", "QCed_File_Exists")


# Generate file names for easyQC 


  # paste0(d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"RawFilePaths_m"], "/", d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"File2"])

  # unique(d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"QCed_p"])





# create QCed noAFO files 
  d$QCednoAFO_fn <- gsub("cpaid", "cpaid_noDup_noAFO.tbl", d$QCed_fn)
  d$QCednoAFO_fn <- gsub("cpaid_noDup_noAFO.tbl.MetaScore.assoc", "cpaid.MetaScore.assoc_noDup_noAFO.tbl", d$QCednoAFO_fn)
  

# See which QCed noAFO files exist

  for (i in 1:nrow(d)){
    qcf <- paste0(d$QCed_p[i] ,"/",  d$QCednoAFO_fn[i] )
    if( file.exists(qcf)) {
      d$QCed_noAFO_File_Exists[i] <- T
      } 
    if(!file.exists(qcf)) {
      d$QCed_noAFO_File_Exists[i] <- F
      }
  }

  table(d$QCed_noAFO_File_Exists)
  #FALSE  TRUE
  # 584    66

  



# Rename files if needed  
  #CHCPfiles <- paste0( unique(d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"QCed_p"]) ,"/",  list.files( unique(d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"QCed_p"])))
  # CHCPrenamfiles <- gsub("230922.", "230922_2.", CHCPfiles)
  #length(CHCPfiles)==length(CHCPrenamfiles)
  #for (i in 1:length(CHCPfiles)){
  #  file.rename(CHCPfiles[i], CHCPrenamfiles[i])
  #}


# Create a variable for who will be in maxIMT MEta-A 

  # ADOS

    Adostds <- unique(d[d$AgeGroup=="ADO", "Study"])
    
    # "ABCD"   "CHCP"   "FAMILY" "GenR1"  "GenR2"  "PANIC"  "SWS"


    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADO_maxLCF_EUR"] <- T
    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADO_maxLCF_ME"] <- T
    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADO_maxLCF_EUR"] <- F
    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADO_maxLCF_ME"] <- T


    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADO_maxLCFp_EUR"] <- T
    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADO_maxLCFp_ME"] <- T
    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADO_maxLCFp_EUR"] <- F
    d[d$AgeGroup=="ADO" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADO_maxLCFp_ME"] <- T
    
    d[ d$Study == "ABCD" & d$AgeGroup=="ADO" & d$cIMT_measure == "maxCCA", "in_ADO_maxLCFp_EUR"] <- T
    d[ d$Study == "ABCD" & d$AgeGroup=="ADO" & d$cIMT_measure == "maxCCA", "in_ADO_maxLCFp_ME"] <- T
    
    d[ d$Study == "CHCP" & d$AgeGroup=="ADO" & d$cIMT_measure == "maxRCF" & d$pop2 == "EUR", "in_ADO_maxLCFp_EUR"] <- T
    d[ d$Study == "CHCP" & d$AgeGroup=="ADO" & d$cIMT_measure == "maxRCF" & d$pop2 == "EUR", "in_ADO_maxLCFp_ME"] <- T
    d[ d$Study == "CHCP" & d$AgeGroup=="ADO" & d$cIMT_measure == "maxRCF" & d$pop2 != "EUR", "in_ADO_maxLCFp_EUR"] <- F
    d[ d$Study == "CHCP" & d$AgeGroup=="ADO" & d$cIMT_measure == "maxRCF" & d$pop2 != "EUR", "in_ADO_maxLCFp_ME"] <- T
    
    d[ d$Study == "PANIC" & d$AgeGroup=="ADO" & d$cIMT_measure == "meanLCF" & d$pop2 == "EUR", "in_ADO_maxLCFp_EUR"] <- T
    d[ d$Study == "PANIC" & d$AgeGroup=="ADO" & d$cIMT_measure == "meanLCF" & d$pop2 == "EUR", "in_ADO_maxLCFp_ME"] <- T
    d[ d$Study == "PANIC" & d$AgeGroup=="ADO" & d$cIMT_measure == "meanLCF" & d$pop2 != "EUR", "in_ADO_maxLCFp_EUR"] <- F
    d[ d$Study == "PANIC" & d$AgeGroup=="ADO" & d$cIMT_measure == "meanLCF" & d$pop2 != "EUR", "in_ADO_maxLCFp_ME"] <- T
    




    #Exceptions
    # ABCD has maxCCA and meanCCA 
    # CHCP has RCF and RCN 
    # PANIC has only meanLCF
    # UKB has max240 


    #See who we got so far 
    unique(d[which(d$in_ADO_maxLCFp_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    
    #      Study AgeGroup cIMT_measure pop2
    # 1     ABCD      ADO       maxCCA  EUR
    # 29    CHCP      ADO       maxRCF  EUR
    # 279 FAMILY      ADO       maxLCF  EUR
    # 301  GenR1      ADO       maxLCF  EUR
    # 337  GenR2      ADO       maxLCF  EUR
    # 409  PANIC      ADO      meanLCF  EUR
    # 519    SWS      ADO       maxLCF  EUR
    
    unique(d[which(d$in_ADO_maxLCFp_ME==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    
    #      Study AgeGroup cIMT_measure   pop2
    # 1     ABCD      ADO       maxCCA    EUR
    # 29    CHCP      ADO       maxRCF    EUR
    # 30    CHCP      ADO       maxRCF nonEUR
    # 279 FAMILY      ADO       maxLCF    EUR
    # 301  GenR1      ADO       maxLCF    EUR
    # 337  GenR2      ADO       maxLCF    EUR
    # 409  PANIC      ADO      meanLCF    EUR
    # 519    SWS      ADO       maxLCF    EUR
    
  # YAD

    yadstds <- unique(d[d$AgeGroup=="YAD", c("Study")])   
 
    # "YFS"    "ALSPAC"
     unique(d[d$AgeGroup=="YAD", c("Study","cIMT_measure","pop2")])    
    d[d$AgeGroup=="YAD" & d$cIMT_measure == "maxLCF", "in_YAD_maxLCF_EUR"] <- T
    d[d$AgeGroup=="YAD" & d$cIMT_measure == "maxLCF", "in_YAD_maxLCFp_EUR"] <- T

    unique(d[which(d$in_YAD_maxLCF_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])


  # ADL

    ADLstds <- unique(d[d$AgeGroup=="ADL", c("Study")])   
    #"CHCP"   "FAMILY" "YFS"    "ALSPAC"
    unique(d[d$AgeGroup=="ADL", c("Study","cIMT_measure","pop2")])    
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF", "in_ADL_maxLCF_EUR"] <- T
    unique(d[which(d$in_ADL_maxLCF_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])

    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADL_maxLCF_EUR"] <- T
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADL_maxLCF_ME"] <- T
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADL_maxLCF_ME"] <- T
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADL_maxLCF_EUR"] <- F

    unique(d[which(d$in_ADL_maxLCF_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])

    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADL_maxLCFp_EUR"] <- T
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_ADL_maxLCFp_ME"] <- T
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADL_maxLCFp_ME"] <- T
    d[d$AgeGroup=="ADL" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_ADL_maxLCFp_EUR"] <- F
    d[d$AgeGroup=="ADL" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 == "EUR", "in_ADL_maxLCFp_EUR"] <- T
    d[d$AgeGroup=="ADL" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 == "EUR", "in_ADL_maxLCFp_ME"] <- T
    d[d$AgeGroup=="ADL" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 != "EUR", "in_ADL_maxLCFp_ME"] <- T
    d[d$AgeGroup=="ADL" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 != "EUR", "in_ADL_maxLCFp_EUR"] <- F

    unique(d[which(d$in_ADL_maxLCFp_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    unique(d[which(d$in_ADL_maxLCFp_ME==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])


  # OAD 


    OADstds <- unique(d[d$AgeGroup=="OAD", c("Study")])   
    # "CHCP"   "NFBC"   "GRACE"  "UKB"    "FAMILY" "ALSPAC"
    unique(d[d$AgeGroup=="OAD", c("Study","cIMT_measure","pop2")])    

    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_OAD_maxLCF_EUR"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_OAD_maxLCF_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_OAD_maxLCF_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_OAD_maxLCF_EUR"] <- F

    unique(d[which(d$in_OAD_maxLCF_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])

    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_OAD_maxLCFp_EUR"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 == "EUR", "in_OAD_maxLCFp_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_OAD_maxLCFp_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxLCF"& d$pop2 != "EUR", "in_OAD_maxLCFp_EUR"] <- F
    d[d$AgeGroup=="OAD" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 == "EUR", "in_OAD_maxLCFp_EUR"] <- T
    d[d$AgeGroup=="OAD" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 == "EUR", "in_OAD_maxLCFp_ME"] <- T
    d[d$AgeGroup=="OAD" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 != "EUR", "in_OAD_maxLCFp_ME"] <- T
    d[d$AgeGroup=="OAD" & d$Study=="CHCP" & d$cIMT_measure == "maxRCF"& d$pop2 != "EUR", "in_OAD_maxLCFp_EUR"] <- F

    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240" & d$pop2 == "EUR" & d$Study=="UKB", "in_OAD_maxLCF_EUR"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240" & d$pop2 == "EUR" & d$Study=="UKB", "in_OAD_maxLCF_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240" & d$pop2 != "EUR" & d$Study=="UKB", "in_OAD_maxLCF_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240" & d$pop2 != "EUR" & d$Study=="UKB", "in_OAD_maxLCF_EUR"] <- F

    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240"& d$pop2 == "EUR" & d$Study=="UKB", "in_OAD_maxLCFp_EUR"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240"& d$pop2 == "EUR" & d$Study=="UKB", "in_OAD_maxLCFp_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240"& d$pop2 != "EUR" & d$Study=="UKB", "in_OAD_maxLCFp_ME"] <- T
    d[d$AgeGroup=="OAD" & d$cIMT_measure == "maxIMT240"& d$pop2 != "EUR" & d$Study=="UKB", "in_OAD_maxLCFp_EUR"] <- F

    unique(d[which(d$in_OAD_maxLCF_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    unique(d[which(d$in_OAD_maxLCF_ME==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    unique(d[which(d$in_OAD_maxLCF_ME==T), c("Study", "AgeGroup", "cIMT_measure", "pop")])

    unique(d[which(d$in_OAD_maxLCFp_EUR==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    unique(d[which(d$in_OAD_maxLCFp_ME==T), c("Study", "AgeGroup", "cIMT_measure", "pop2")])
    unique(d[which(d$in_OAD_maxLCFp_ME==T), c("Study", "AgeGroup", "cIMT_measure", "pop")])


# Create table by meta-group 

  d[which(d$in_OAD_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==T), "in_OAD_maxLCFp_EUR_MW_adjBMI"] <- T
  d[which(d$in_ADO_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==T), "in_ADO_maxLCFp_EUR_MW_adjBMI"] <- T
  d[which(d$in_ADL_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==T), "in_ADL_maxLCFp_EUR_MW_adjBMI"] <- T
  d[which(d$in_YAD_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==T), "in_YAD_maxLCFp_EUR_MW_adjBMI"] <- T
  
  d[which(d$in_OAD_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==T), "in_OAD_maxLCFp_EUR_M_adjBMI"] <- T
  d[which(d$in_ADO_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==T), "in_ADO_maxLCFp_EUR_M_adjBMI"] <- T
  d[which(d$in_ADL_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==T), "in_ADL_maxLCFp_EUR_M_adjBMI"] <- T
  d[which(d$in_YAD_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==T), "in_YAD_maxLCFp_EUR_M_adjBMI"] <- T
  
  d[which(d$in_OAD_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==T), "in_OAD_maxLCFp_EUR_W_adjBMI"] <- T
  d[which(d$in_ADO_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==T), "in_ADO_maxLCFp_EUR_W_adjBMI"] <- T
  d[which(d$in_ADL_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==T), "in_ADL_maxLCFp_EUR_W_adjBMI"] <- T
  d[which(d$in_YAD_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==T), "in_YAD_maxLCFp_EUR_W_adjBMI"] <- T
  
  
  d[which(d$in_OAD_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==F), "in_OAD_maxLCFp_EUR_MW_noBMI"] <- T
  d[which(d$in_ADO_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==F), "in_ADO_maxLCFp_EUR_MW_noBMI"] <- T
  d[which(d$in_ADL_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==F), "in_ADL_maxLCFp_EUR_MW_noBMI"] <- T
  d[which(d$in_YAD_maxLCFp_EUR==T & d$sexGroup=="MW" & d$AdjBMI==F), "in_YAD_maxLCFp_EUR_MW_noBMI"] <- T
  
  d[which(d$in_OAD_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==F), "in_OAD_maxLCFp_EUR_M_noBMI"] <- T
  d[which(d$in_ADO_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==F), "in_ADO_maxLCFp_EUR_M_noBMI"] <- T
  d[which(d$in_ADL_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==F), "in_ADL_maxLCFp_EUR_M_noBMI"] <- T
  d[which(d$in_YAD_maxLCFp_EUR==T & d$sexGroup=="M" & d$AdjBMI==F), "in_YAD_maxLCFp_EUR_M_noBMI"] <- T
  
  d[which(d$in_OAD_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==F), "in_OAD_maxLCFp_EUR_W_noBMI"] <- T
  d[which(d$in_ADO_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==F), "in_ADO_maxLCFp_EUR_W_noBMI"] <- T
  d[which(d$in_ADL_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==F), "in_ADL_maxLCFp_EUR_W_noBMI"] <- T
  d[which(d$in_YAD_maxLCFp_EUR==T & d$sexGroup=="W" & d$AdjBMI==F), "in_YAD_maxLCFp_EUR_W_noBMI"] <- T






  d[which(d$Study=="UKB" & d$pop=="NbEUR"), "in_OAD_maxLCFp_EUR_MW_noBMI"]
  
  ==T & d$sexGroup=="MW" & d$AdjBMI==F), "in_OAD_maxLCFp_EUR_MW_noBMI"] <- T





# Create variable for files for AFO and no liftover 
  d$aff_fn <- gsub("CLEANED.", "", d[,"QCed_fn"])
  d$aff_fn <- gsub(".cpaid.gz", ".AFCHECK.outlier.txt", d[,"aff_fn"])
  d$aff_fn <- gsub(".cpaid.MetaScore.assoc.gz", ".AFCHECK.outlier.MetaScore.assoc.txt", d$aff_fn)

  d$noLO_fn <- gsub("CLEANED.", "", d[,"QCed_fn"])
  d$noLO_fn <- gsub(".cpaid.gz", ".notlifted.txt",  d$noLO_fn)
  d$noLO_fn <- gsub(".cpaid.MetaScore.assoc.gz", ".notlifted.MetaScore.assoc.txt",  d$noLO_fn)


# See who's missing QCed file
  d[which(!d$QCed_File_Exists & d$in_ADO_maxLCF_EUR),c("Study", "sexGroup", "AdjBMI2")]

# See who's missing the file without AFO 
   d[which(!d$QCed_noAFO_File_Exists & d$in_ADO_maxLCF_EUR) ,c("Study", "sexGroup",  "AdjBMI2", "QCednoAFO_fn")] #, "QCednoAFO_fn"




# Get list of files for meta-A 

  mtav <- "maxLCFp"
  selcv <- "in_ADO_maxLCFp_EUR_MW_adjBMI" 
  #meta-A input files 
    fns <- paste0(d[which(d[,selcv]) , "QCed_p"], "/", d[which(d[,selcv]) , "QCednoAFO_fn"])
    if(any(!file.exists(fns))) fns[which(!file.exists(fns))]

 # MetaA path + file name  
  maf <- unique(paste0(ma_p, "/",d[which(d[,selcv]) , "AgeGroup"], "/",d[which(d[,selcv]) , "pop2"], "/", mtav, "/", d[which(d[,selcv]) , "AdjBMI2"], "/241015_", mtav, "_", d[which(d[,selcv]) , "AgeGroup"], "_",d[which(d[,selcv]) , "pop2"], "_", d[which(d[,selcv]) , "sexGroup"], "_", d[which(d[,selcv]) , "AdjBMI2"], "_GWAMA_noAFO"))

 # MetaA path 
  MetaA_p <- paste0(ma_p, "/",d[which(d[,selcv]) , "AgeGroup"], "/", mtav, "/", d[which(d[,selcv]) , "AgeGroup"], "_",d[which(d[,selcv]) , "pop2"], "/", d[which(d[,selcv]) , "sexGroup"], "/", d[which(d[,selcv]) , "AdjBMI2"])

 # MetaA file name  
  paste0( "/241015_", mtav, "_", d[which(d[,selcv]) , "AgeGroup"], "_",d[which(d[,selcv]) , "pop2"], "_", d[which(d[,selcv]) , "sexGroup"], "_", d[which(d[,selcv]) , "AdjBMI2"], "_GWAMA_noAFO")

# Save output table 
  write.table(d,paste0(o_p,"/FilesTables/241022_GWAS_Files.tbl"), row.names=F, col.names=T, quote=F, sep="\t")






tbg <- d[which(d$in_YAD_maxLCFp_EUR_W_noBMI | d$in_YAD_maxLCFp_EUR_M_noBMI | d$in_YAD_maxLCFp_EUR_MW_noBMI | d$in_ADO_maxLCFp_EUR_W_noBMI | d$in_ADO_maxLCFp_EUR_M_noBMI | d$in_ADO_maxLCFp_EUR_MW_noBMI |d$in_ADL_maxLCFp_EUR_W_noBMI | d$in_ADL_maxLCFp_EUR_M_noBMI | d$in_ADL_maxLCFp_EUR_MW_noBMI | d$in_OAD_maxLCFp_EUR_W_noBMI | d$in_OAD_maxLCFp_EUR_M_noBMI | d$in_OAD_maxLCFp_EUR_MW_noBMI ), c("Study", "pop", "pop2", "AgeGroup","AdjBMI","AdjBMI2","sexGroup","cIMT_measure", "QCed_fn", "MetaA_p")]
tbg$in_maxLCFp <- T
tbg$Gr <- paste(tbg$"AgeGroup", "maxLCFp", tbg$pop2, tbg$sexGroup, tbg$AdjBMI2, sep="_")

tbg$metaA_fn <- paste("maxLCFp",tbg$"AgeGroup",  tbg$pop2, tbg$sexGroup, tbg$AdjBMI2, "_GWAMA_1.TBL", sep="_")


write.table(tbg,paste0(o_p,"/241010_GWAMA_Files.tbl"), row.names=F, col.names=T, quote=F, sep="\t")

# SelectionVariables 

selectVs <- c("in_ADO_maxLCF_EUR" , "in_YAD_maxLCFp_EUR", "in_ADO_maxLCFp_EUR", "in_ADO_maxLCFp_ME", "in_YAD_maxLCF_EUR",  "in_ADO_maxLCF_ME",   "in_ADL_maxLCF_ME", "in_ADL_maxLCF_EUR",  "in_ADL_maxLCFp_EUR", "in_ADL_maxLCFp_ME", "in_OAD_maxLCF_EUR",  "in_OAD_maxLCF_ME",   "in_OAD_maxLCFp_EUR", "in_OAD_maxLCFp_ME")

