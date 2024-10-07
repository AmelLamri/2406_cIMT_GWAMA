library(qqman)
library(R.utils)
library(data.table)


# m<- "/home/lamria/avr/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results
#m_p <- "/home/lamria/avr/Projects/byStudy/EGGC/2211_cIMT_GWAMA/Results/GWAMA/"#
m_p <- "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI"
rgma <- paste0(m_p,"/Results/GWAMA/")
kg_p <- "/home/lamria/avr/Data/1000Genomes/20130502_Release/Genotypes/PLINK/Sample_Subsets/Europeans/SNP_Subsets/mac10/"

# Create empty list 
  eumaxp <- list() 
  for (i in 1:4){
    eumaxp[[i]] <- list() 
    for (j in 1:3 ){  
      eumaxp[[i]][[j]] <- list()  
      for (k in 1:2){  
        tmp <- list()  
        }
      }
    }


# Open files tables 

  lft <- read.delim(paste0(m_p, "/Results/241004_GWAS_Files.tbl"), head=T, string=F)
  gwi <- colnames(lft)[ grep("_maxLCFp_EUR_", colnames(lft))]




# Create a table of GWAMA files 
  fts <- as.data.frame(matrix(ncol=11, nrow=3*5*2))
  fts[,1] <- rep( c("ADO", "YAD", "ADL", "OAD", "AAD"), 6)
  fts[,2] <- rep( c("MW", "M", "W"), 10)
  fts[,3] <- rep( c("noBMI", "adjBMI"), 15)
  fts[,4] <- paste(fts[,1], fts[,2] ,fts[,3] , "eur_maxp",sep="_")
  fts[,5] <- 1:5
  fts[,6] <- 1:3
  fts[,7] <- 1:2
  fts <- fts[order(fts[,5], fts[,6], fts[,7]), ]
  fts[,8] <- paste0(rgma,"/",fts[,1], "/EUR/maxLCFp/", fts[,3],"/240926_maxIMT_",fts[,1],"_EUR_",fts[,2], "_", fts[,3] , "_GWAMA.TBL")
  fts[,9] <- gsub("_GWAMA.TBL","_GWAMA_ChrBp.TBL",fts[,8])
  colnames(fts) <- c("AgeGr", "SexGr", "Adj", "abv", "i", "j", "k", "ifn", "ifn2", "N_SNP", "N_SNP_nodup")


# Correct file names in GWAMA files table 
  fts[fts$AgeGr=="OAD","ifn"] <- gsub("_GWAMA1.TBL","_GWAMA_1.TBL", fts[fts$AgeGr=="OAD","ifn"])
  fts[fts$AgeGr=="OAD","ifn"] <- gsub("240926","241004", fts[fts$AgeGr=="OAD","ifn"])


# Open files tables 

  lft <- read.delim(paste0(m_p, "/Results/241004_GWAS_Files.tbl"), head=T, string=F)
  eumaxpfs<-lft[which(lft$in_OAD_maxLCFp_EUR_MW_noBMI |lft$in_ADL_maxLCFp_EUR_MW_noBMI | lft$in_YAD_maxLCFp_EUR_MW_noBMI |lft$in_ADO_maxLCFp_EUR_MW_noBMI | lft$in_OAD_maxLCFp_EUR_MW_adjBMI | lft$in_ADL_maxLCFp_EUR_MW_adjBMI | ft$in_YAD_maxLCFp_EUR_MW_adjBMI |lft$in_ADO_maxLCFp_EUR_MW_adjBMI ), ]

  eumaxpfs$QCedgz_fn <- paste0(eumaxpfs[,"QCed_fn"], ".gz")
  eumaxpfs$QCed_noDup <- gsub(".gz", "_noDup.tbl",  eumaxpfs$QCed_fn)
  eumaxpfs$QCed_noDupgz_fn <- paste0(eumaxpfs[,"QCed_noDup"], ".gz")
  eumaxpfs$QCed_noDupNoAFO_fn <- gsub("_noDup.tbl", "_noDup_noAFO.tbl",  eumaxpfs$QCed_noDup)
  eumaxpfs$QCed_noDupNoAFOgz_fn <-  paste0(eumaxpfs[,"QCed_noDupNoAFO_fn"], ".gz")
  eumaxpfs$aff_fn <- gsub("CLEANED.", "", eumaxpfs[,"QCed_fn"])
  eumaxpfs$aff_fn <- gsub(".cpaid.gz", ".AFCHECK.outlier.txt", eumaxpfs[,"aff_fn"])
  eumaxpfs$noLO_fn <- gsub("CLEANED.", "", eumaxpfs[n,"QCed_fn"])
  eumaxpfs$noLO_fn <- gsub(".cpaid.gz", ".notlifted.txt",  eumaxpfs$noLO_fn)
  write.table(eumaxpfs, paste0(m_p, "/Results/241004_GWAS_Files_maxp_MW_EUR_noAFO.tbl"),row.names=F, quote=F, sep="\t")


options(warn=2)

    # Open gwama data and save for locus zoom 
    for (n in 1:nrow(eumaxpfs)){ #which(eumaxpfs[,1]=="OAD" :nrow()
    
        if(!file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_fn"]))){
        cat ( "file", eumaxpfs[n,4], " not found \n")
        next()
        }else{
        if(!file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"])) & !file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"], ".gz" ))& ! file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"] ))){
            cat ( "running ", eumaxpfs[n,4],n,  " \n")
            tmp <- as.data.frame(fread(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_fn"])))
            eumaxpfs$N_SNP_QCed[n] <- nrow(tmp)
            tmp$bp_QCed<-sapply(strsplit(tmp$cpaid,      ":"), `[`, 2)
           
            # re-order chr and position before save for locus zoom
            tmpord <- tmp[which(!is.na(tmp$CHR)),]
            tmpord$chrpos <- paste0(tmpord$CHR, tmpord$bp_QCed)
            eumaxpfs$N_SNP_QCed_Dupl[n] <- sum(duplicated(tmpord$chrpos ) | duplicated(tmpord$chrpos, fromLast=T))

            tmpord <-tmpord[!duplicated(tmpord$chrpos ) & !duplicated(tmpord$chrpos, fromLast=T), ] 
            #check 
            if(length(names(table(tmpord$bp_QCed==tmpord$POS)))>1 | names(table(tmpord$bp_QCed==tmpord$POS)) != "TRUE") cat( "\n WARNING !!!!! , check BP for ", n,  eumaxpfs$QCed_fn, "\n")
            eumaxpfs[n,"N_SNP_QCed_noDup"] <- nrow(tmpord)

            eumaxpfs$N_SNP_QCed_Dupl_cpaid[n] <- sum(duplicated(tmpord$cpaid ) | duplicated(tmpord$cpaid, fromLast=T))
            tmpord <-tmpord[!duplicated(tmpord$cpaid ) & !duplicated(tmpord$cpaid, fromLast=T), ] 
            eumaxpfs[n,"N_SNP_QCed_noDup2"] <- nrow(tmpord)

            tmpord <-tmpord[order(tmpord[,"CHR"], tmpord[,"bp_QCed"]), ] 
            # Save file with chr and position for locus zoom 
            #if(!file.exists(eumaxpfs[n,"QCedgz_fn"]) ) 
            write.table(tmpord,paste0(eumaxpfs[n,"QCed_p"], "/", eumaxpfs[n,"QCed_noDup"] ), row.names=F, col.names=T, quote=F, sep="\t")
            # Open list of variants that fail AF check 


            aff <- read.table(paste0(eumaxpfs[n,"QCed_p"], "/", eumaxpfs[n,"aff_fn"]), head=T, string=F)
            eumaxpfs[n,"N_SNP_QCed_noDup_AFO"] <- length(which(tmpord$cpaid %in% aff$cpaid))
            tmpord_noafo <- tmpord[which(!tmpord$cpaid %in% aff$cpaid ),]
            eumaxpfs[n,"N_SNP_QCed_noDup_noAFO"] <- nrow(tmpord_noafo)
            write.table(tmpord_noafo, paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"]) , row.names=F, col.names=T, quote=F, sep="\t")

            #noLO <- as.data.frame(fread(paste0(eumaxpfs[n,"QCed_p"], "/", eumaxpfs[n,"noLO_fn"])))
            #table(noLO$SNPID %in% tmpord_noafo$SNPID  )
            #summary(tmpord_noafo$N)


            #if(!file.exists(paste0(eumaxpfs[n,"QCedgz_fn"], ".gz"))) 
            gzip (paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDup"]))
            gzip (paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"]))


gzip("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ABCD/QCed/ADO/maxCCA/EUR/CLEANED.EGG_cIMT.maxCCA.ADO.MW.noBMI.EUR.EM.20230626.cpaid_noDup_noAFO.tbl")
gzip("/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI//Results/GWAS/ABCD/QCed/ADO/maxCCA/EUR/CLEANED.EGG_cIMT.maxCCA.ADO.MW.noBMI.EUR.EM.20230626.cpaid_noDup.tbl")


            #tmpord$P-value <- as.numeric(tmpord$P-value)
            # if(!file.exists(gsub("_GWAMA.TBL", "_GWAMA_manh.", eumaxpfs[n,"QCed_fn"]))){
            #   jpeg(gsub("_GWAMA.TBL", "_GWAMA_manh.", eumaxpfs[n,"QCed_fn"]))
            #     manhattan(tmpord[which(tmpord$chr %in% 1:22),],p="P-value", chr="chr", bp="bp", snp="MarkerName" )
            #   dev.off()
            # } 

        }
      }
    }

