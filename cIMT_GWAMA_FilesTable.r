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




d <- as.data.frame(read_excel(paste0(o_p,"/240925_Files.xlsx" )))
stds <- unique(d$Study)


#Fill in empty raw file paths 
  table(d$pop, useNA="ifany")
  table(is.na(d$RawFilePaths), useNA="ifany")
  table(d[is.na(d$RawFilePaths), "Study"] )
  table(is.na(d$Study) )
  table(d$AgeGroup, useNA="ifany")

  d[is.na(d$RawFilePaths), "RawFilePaths"]  <- paste0("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/", d[is.na(d$RawFilePaths), "Study"], "/Raw")

# Clean column names of the table 
colnames(d) <- gsub(" ", "_", colnames(d))
# correct AgeGroup Var

  table(d$AgeGroup, useNA="ifany")
 d[d$AgeGroup=="CLD", "AgeGroup"] <- "ADO"
 d[d$AgeGroup=="xAGE", "AgeGroup"] <- "AADxAGE"

ages <- c("OAD", "ADL", "YAD", "AAD", "ADO", "AADxAGE")


# correct some file names 
 d[d$Study=="FAMILY", "File"] <- gsub("ALL", "MW", d[d$Study=="FAMILY", "File"])
 d[d$Study=="FAMILY", "File"] <- gsub("WOMEN", "W", d[d$Study=="FAMILY", "File"])
 d[d$Study=="FAMILY", "File"] <- gsub("MEN", "M", d[d$Study=="FAMILY", "File"])

 d[d$Study=="ALSPAC", "File"] <- gsub("minIMT", "minLCF", d[d$Study=="ALSPAC", "File"])
 d[d$Study=="ALSPAC", "File"] <- gsub("meanIMT", "meanLCF", d[d$Study=="ALSPAC", "File"])
 d[d$Study=="ALSPAC", "File"] <- gsub("maxIMT", "maxLCF", d[d$Study=="ALSPAC", "File"])

# Correct some file paths 
d[d$Study=="FAMILY", "RawFilePaths"] <- gsub("2301_cIMT_GWAS_PHRI", "2301_FAMILY_cIMT_GWAS_PHRI", d[d$Study=="FAMILY", "RawFilePaths"])

#Create Cleaned Raw files variable 

d$RawFilePaths_m <- d$RawFilePaths
d[d$Study=="CHCP", "RawFilePaths_m"] <- paste0(d[d$Study=="CHCP", "RawFilePaths"], "/Modified/")
d[d$Study=="ABCD", "RawFilePaths_m"] <- paste0(d[d$Study=="ABCD", "RawFilePaths"], "/Modified/")
d[d$Study=="ABCD", "RawFilePaths_m"] <- paste0(d[d$Study=="ABCD", "RawFilePaths"], "/Modified/")


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
    d[d$Study == "CHCP", "File2"] <- gsub(".zip", "_2.txt.gz", d[d$Study == "CHCP", "File"])
    d[d$Study == "GRACE", "File2"] <- gsub("20240903.txt.gz", "20240903_tmp.txt.gz", d[d$Study == "GRACE", "File"])
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
  d[d$cIMT_measure== "min240", "cIMT_measure"] <- "minIMT240"
  d[d$cIMT_measure== "mean240", "cIMT_measure"] <- "meanIMT240"
  d[d$cIMT_measure== "max240", "cIMT_measure"] <- "maxIMT240"

  for (i in 1:nrow(d)){
    if(length(grep(d$cIMT_measure[i], d$File[i]) ) != 1) cat (i, "\n")
    } # all good 

# Check that sex group matches 
  table(d$sexGroup, useNA="ifany")
  #  M  MW   W
  # 180 248 198
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






# Makesure sure all modified Files exist
table(duplicated(d$File))
#FALSE



fs<- paste0(d$RawFilePaths_m, "/", d$File2) 

for (i in 1:length(fs)){
  #if(d$Study[i] %in% c( "GRACE", "ALSPAC", "FAMILY") & d$cIMT_measure[i] %in% c("meanLCF", "minLCF")) next()
  #if(d$Study[i] == "UKB") {
  #  if(d$cIMT_measure[i] %in% c("mean240", "min240")) next()
  #  if(d$pop[i] %in% c("AFR", "SAS")) next()
  #  if(length(grep("NbEUR", d$File[i]))==1 ) next()
  #}
  #if(d$Study[i] %in% c("ALSPAC", "FAMILY") & (d$"AgeGroup"[i] %in% c("ADL", "OAD", "AAD")) & (d$"sexGroup"[i] %in% c("MW", "M")))next()
  #if(d$Study[i] ==c("CHCP") & (d$"cIMT_measure"[i] %in% c("minRCN", "minRCF", "meanRCN", "meanRCF", "maxRCN" )))next()

  #if(d$Study[i] %in% c("GRACE") & d$pop[i] %in% c("EUR"))next()
  if( file.exists(fs[i])) {
    d$File_Exists[i] <- T
    } 
  if(!file.exists(fs[i])) {
    d$File_Exists[i] <- F
    }
}
table(d$File_Exists)

unique(d[!d$File_Exists, c("Study", "pop", "AgeGroup", "cIMT_measure", "RawFilePaths_m", "File2")])

unique(d[!d$File_Exists, c("Study", "pop", "AgeGroup", "cIMT_measure")])


# TO DO !!!  
# GENERATE THESE FILES That are missing  above !!! unique(d[!d$File_Exists, c("Study", "pop", "AgeGroup", "cIMT_measure")])

d$QCed_fn <- paste0("CLEANED.", gsub(".txt.gz", "", d$File2),".cpaid.gz")
d$QCed_fn <- gsub(".gz.cpaid.gz", ".cpaid.gz",d$QCed_fn  ) 



  adjsd <- d$AdjBMI
  table(adjsd)
  d[d$AdjBMI=="T","AdjBMI2"] <- "adjBMI"
  d[d$AdjBMI=="F","AdjBMI2"] <- "noBMI"
  table(d$AdjBMI)
  table(d$AdjBMI2)


  d$MetaA_p <- paste0(ma_p, "/", d$AgeGroup, "/", d$cIMT_measure, "/", d$pop2, "/", d$sexGroup, "/", d$AdjBMI2)







# See if QCed files exist 

d$QCed_p <- paste0(gw_p,"/",d$Study, "/QCed","/", d$AgeGroup,"/", d$cIMT_measure,"/",d$pop )
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





#see who's missing 
d[!d$QCed_File_Exists & d$AgeGroup =="OAD",c("Study", "pop", "AgeGroup", "cIMT_measure", "sexGroup")]
d[!d$QCed_File_Exists & d$AgeGroup =="OAD" & d$pop2=="EUR" & d$cIMT_measure %in% c("maxLCF", "maxIMT240") ,c("Study", "pop", "AgeGroup", "cIMT_measure", "sexGroup")]


d[ d$AgeGroup =="OAD" & d$pop2=="EUR" & d$cIMT_measure %in% c("maxLCF", "maxIMT240") & d$Study=="NFBC",c("QCed_File_Exists")]


d[d$Study=="GenR1" & d$AgeGroup =="ADO" & d$pop2 =="EUR" & d$cIMT_measure =="maxLCF" & d$sexGroup=="MW" & d$AdjBMI2=="adjBMI" ,"File2"] 


d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$cIMT_measure =="maxRCF" & d$pop =="EUR" ,c("Study", "pop", "AgeGroup", "cIMT_measure", "sexGroup", "QCed_fn", "QCed_p", "pop")]

#generate data for easyQC 


 paste0(d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"RawFilePaths_m"], "/", d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"File2"])

 unique(d[!d$QCed_File_Exists & d$Study=="CHCP" & d$AgeGroup =="OAD" & d$pop =="EUR" & d$cIMT_measure =="maxRCF" ,"QCed_p"])



## Failed to open file '/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/PANIC/QCed/ADO/meanLCF/EUR/CLEANED.PANIC_cIMT.meanLCF.ADO.MW.adjBMI.EUR.TS.231218.cpaid.gz'





#rename files if needed  
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

d$AdjBMI <- as.logical(d$AdjBMI)

#create table by meta-group 

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



write.table(d,paste0(o_p,"/241004_GWAS_Files.tbl"), row.names=F, col.names=T, quote=F, sep="\t")

tbg <- d[which(d$in_YAD_maxLCFp_EUR_W_noBMI | d$in_YAD_maxLCFp_EUR_M_noBMI | d$in_YAD_maxLCFp_EUR_MW_noBMI | d$in_ADO_maxLCFp_EUR_W_noBMI | d$in_ADO_maxLCFp_EUR_M_noBMI | d$in_ADO_maxLCFp_EUR_MW_noBMI |d$in_ADL_maxLCFp_EUR_W_noBMI | d$in_ADL_maxLCFp_EUR_M_noBMI | d$in_ADL_maxLCFp_EUR_MW_noBMI | d$in_OAD_maxLCFp_EUR_W_noBMI | d$in_OAD_maxLCFp_EUR_M_noBMI | d$in_OAD_maxLCFp_EUR_MW_noBMI ), c("Study", "pop", "pop2", "AgeGroup","AdjBMI","AdjBMI2","sexGroup","cIMT_measure", "QCed_fn", "MetaA_p")]
tbg$in_maxLCFp <- T
tbg$Gr <- paste(tbg$"AgeGroup", "maxLCFp", tbg$pop2, tbg$sexGroup, tbg$AdjBMI2, sep="_")

tbg$metaA_fn <- paste("maxLCFp",tbg$"AgeGroup",  tbg$pop2, tbg$sexGroup, tbg$AdjBMI2, "_GWAMA_1.TBL", sep="_")


write.table(tbg,paste0(o_p,"/241004_GWAMA_Files.tbl"), row.names=F, col.names=T, quote=F, sep="\t")

# SelectionVariables 

selectVs <- c("in_ADO_maxLCF_EUR" , "in_YAD_maxLCFp_EUR", "in_ADO_maxLCFp_EUR", "in_ADO_maxLCFp_ME", "in_YAD_maxLCF_EUR",  "in_ADO_maxLCF_ME",   "in_ADL_maxLCF_ME", "in_ADL_maxLCF_EUR",  "in_ADL_maxLCFp_EUR", "in_ADL_maxLCFp_ME", "in_OAD_maxLCF_EUR",  "in_OAD_maxLCF_ME",   "in_OAD_maxLCFp_EUR", "in_OAD_maxLCFp_ME")


#TO DO< update metaA path in the following table !! 

# GWAMA files 
 fsd <- read.table(paste0(o_p,"/GWAMA_Files_info.txt" ), head=T, string=F)
 fsd$IMT <- paste0( fsd$IMT , "p")
  fsd$IMT <- gsub("IMT", "LCF", fsd$IMT )

 fsd$MetaAFilePath <- paste0(fsd$MetaAFilePath, "/",fsd$AgeGroup, "/",fsd$IMT, "/",fsd$Pop,  "/",fsd$Adj)

  # d$MetaA_p <- paste0(ma_p, "/", d$AgeGroup, "/", d$cIMT_measure, "/", d$pop2, "/", d$sexGroup, "/", d$AdjBMI2)



# Create files 
for (i in ages){
  for (j in pops2){
    for (k in imt){
      for (l in adj){
        gwamap <- paste0(ma_p, "/", i, "/", j, "/",k , "/", l )
if (!dir.exists(gwamap)) dir.create(gwamap, rec=T)
      }}}}




write.table(fsd, paste0(o_p,"/240926_GWAMA_Files.tbl"), row.names=F, col.names=T, quote=F)
