
 # #Study-QC
 

    # https://www.uni-regensburg.de/medizin/epidemiologie-praeventivmedizin/genetische-epidemiologie/charge-gli
    # https://homepages.uni-regensburg.de/~wit59712/easyqc/EasyQC_9.0_Commands_140918_2.pdf



    library(EasyQC2)
    library(data.table)
    library(R.utils)

  m_p="/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/"
  qc2_p=paste0(m_p, "/Scripts/QC/ALSPAC/")
  or_p=paste0(m_p, "/Results/GWAS/ALSPAC/Raw")
  o_p=paste0(m_p, "/Results/GWAS/ALSPAC/QCed")
  ir_p <- "/genetics/Users/lamria/ALSPAC/Projects/2406_ALSPAC_cIMT_GWAS_PHRI/Outputs"


# hd<- fread("/genetics/Users/lamria/ALSPAC/Projects/2406_ALSPAC_cIMT_GWAS_PHRI/Outputs/YAD/ALSPAC_cIMT.maxLCF.YAD.MW.adjBMI.EUR.AL.20240702.txt.gz")

# hdc<- fread("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ALSPAC/QCed/YAD/maxLCF/EUR/CLEANED.ALSPAC_cIMT.maxLCF.YAD.MW.adjBMI.EUR.AL.20240702.cpaid.gz")
#table(is.na(d[d$INFO_TYPE==0,"INFO"]))
#unique(d[d$INFO_TYPE==0,"INFO"])
# Prepare files for QC 

#(summary(d[d$INFO_TYPE!=0,"INFO"]))


   # d <- fread(paste0(if2, ".gz"))
  options(warn=2)

  for (pheni in c("maxLCF")){ #
    for (sexi in c("MW", "M", "W")){ # "MW", 
      for (adji in c("noBMI", "adjBMI") ){
        for(popi in c("EUR")){ # , "nonEUR"
          for (agi in c("YAD", "ADL", "OAD", "AAD")){# , "AAD"
            filei <- paste0("/ALSPAC_cIMT.", pheni, ".", agi,".",sexi,".",adji,".",popi, ".AL.20240702")
            if1<- paste0(ir_p,"/", agi, "/", filei, ".txt")
            if2<- paste0(or_p,"/", agi, "/", filei,  "_3.txt")
            if(!dir.exists(paste0(or_p,"/", agi)))dir.create(paste0(or_p,"/", agi))

            if( agi %in% c("ADL", "OAD", "AAD") &  sexi %in% c("MW","M")) next()
            if(!file.exists(paste0(if2, ".gz"))){
              
              cat (" \n reading file", pheni, sexi, adji, popi, agi, "\n") 


              fam <- as.data.frame(fread(paste0(if1, ".gz")))
              colnames(fam)[colnames(fam)=="INFO_TYPE"] <- "INFO_Measure"
              colnames(fam)[colnames(fam)=="IMPUTED"] <- "INFO_TYPE"
              summary(fam[which(fam$INFO_TYPE==0),"INFO"])
                            table((fam$INFO_TYPE))

              fam[which(fam$INFO_TYPE==0),"INFO"] <- 1

              # numDrop.BETA.missing

              # numDrop.Imputed.MissingInformation


              table(is.na(fam$INFO))
              table(is.na(fam$INFO_TYPE))
              table((fam$INFO_TYPE))

              table(is.na(fam$BETA))
              summary(fam[is.na(fam$BETA), "INFO"])
              summary(fam[is.na(fam$BETA), "EAF"])

              cat (" \n            saving file \n") 
              if(any(is.na(fam$INFO) | is.na(fam$INFO_TYPE) ))warning( " NAs in info \n")
              write.table(fam[!is.na(fam$BETA),], if2, row.names=F, col.names=T, quote=F, sep="\t" )
              cat (" \n            zipping file \n") 
              gzip(if2, ext="gz", FUN=gzfile)
              #if(file.exists(if1))file.remove(if1)

            }
          } 
        }
      }
    }
  }




    # YAD maxLCF
    setwd(paste0(o_p, "/YAD/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_YAD.ecf"))

    # OAD maxLCF
    setwd(paste0(o_p, "/OAD/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_OAD.ecf"))


    # ADL maxLCF
    setwd(paste0(o_p, "/ADL/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/EGG_cIMT_EasyQC2_ALSPAC_maxLCF_ADL.ecf"))

    # AAD maxLCF
    setwd(paste0(o_p, "/AAD/maxLCF/EUR"))
    EasyQC2(paste0(qc2_p, "/1EGG_cIMT_EasyQC2_ALSPAC_maxLCF_AAD.ecf"))
