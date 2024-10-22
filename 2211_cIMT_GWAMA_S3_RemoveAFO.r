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

  lft <- read.delim(paste0(m_p, "/Results/FilesTables/241018_GWAS_Files.tbl"), head=T, string=F)

  # eumaxpfs <- lft[which(lft$in_OAD_maxLCFp_EUR_MW_noBMI |lft$in_ADL_maxLCFp_EUR_MW_noBMI | lft$in_YAD_maxLCFp_EUR_MW_noBMI |lft$in_ADO_maxLCFp_EUR_MW_noBMI | lft$in_OAD_maxLCFp_EUR_MW_adjBMI | lft$in_ADL_maxLCFp_EUR_MW_adjBMI | lft$in_YAD_maxLCFp_EUR_MW_adjBMI |lft$in_ADO_maxLCFp_EUR_MW_adjBMI | lft$in_OAD_maxLCFp_EUR_M_noBMI |lft$in_ADL_maxLCFp_EUR_M_noBMI | lft$in_YAD_maxLCFp_EUR_M_noBMI |lft$in_ADO_maxLCFp_EUR_M_noBMI | lft$in_OAD_maxLCFp_EUR_M_adjBMI | lft$in_ADL_maxLCFp_EUR_M_adjBMI | lft$in_YAD_maxLCFp_EUR_M_adjBMI |lft$in_ADO_maxLCFp_EUR_M_adjBMI | lft$in_OAD_maxLCFp_EUR_W_noBMI |lft$in_ADL_maxLCFp_EUR_W_noBMI | lft$in_YAD_maxLCFp_EUR_W_noBMI |lft$in_ADO_maxLCFp_EUR_W_noBMI | lft$in_OAD_maxLCFp_EUR_W_adjBMI | lft$in_ADL_maxLCFp_EUR_W_adjBMI | lft$in_YAD_maxLCFp_EUR_W_adjBMI |lft$in_ADO_maxLCFp_EUR_W_adjBMI ), ]
    
    eumaxpfs <- lft[which(lft$Study=="FAMILY" & lft$in_ADL_maxLCFp_EUR  & lft$sexGroup=="W"), ]
    


  #ados <- lft[which((lft$in_ADO_maxLCFp_EUR_MW_noBMI |  lft$in_ADO_maxLCFp_EUR_MW_adjBMI |  lft$in_ADO_maxLCFp_EUR_M_noBMI |  lft$in_ADO_maxLCFp_EUR_M_adjBMI |  lft$in_ADO_maxLCFp_EUR_W_noBMI |  lft$in_ADO_maxLCFp_EUR_W_adjBMI) & lft$Study=="SWS" ), ]


  #yad <- lft[which((lft$in_YAD_maxLCFp_EUR_MW_noBMI |  lft$in_YAD_maxLCFp_EUR_MW_adjBMI |  lft$in_YAD_maxLCFp_EUR_M_noBMI |  lft$in_YAD_maxLCFp_EUR_M_adjBMI |  lft$in_YAD_maxLCFp_EUR_W_noBMI |  lft$in_YAD_maxLCFp_EUR_W_adjBMI) & lft$Study=="YFS" ), ]


  #adl<- lft[which((lft$in_ADL_maxLCFp_EUR_MW_noBMI |  lft$in_ADL_maxLCFp_EUR_MW_adjBMI |  lft$in_ADL_maxLCFp_EUR_M_noBMI |  lft$in_ADL_maxLCFp_EUR_M_adjBMI |  lft$in_ADL_maxLCFp_EUR_W_noBMI |  lft$in_ADL_maxLCFp_EUR_W_adjBMI) & lft$Study=="YFS"), ]
  # eumaxpfs <- adl
  
  file.exists(paste0(eumaxpfs[,"QCed_p"]))
  file.exists(paste0(eumaxpfs[,"QCed_p"], "/",eumaxpfs[,"QCed_fn"]))


  options(warn=2)

    # Open gwama data and save for locus zoom 
    for (n in  1:nrow(eumaxpfs)){ #which(eumaxpfs[,1]=="OAD" :nrow() #


    
        if(!file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_fn"]))){
        cat ( "file", eumaxpfs[n,2], eumaxpfs[n,4],eumaxpfs[n,6],eumaxpfs[n,7],eumaxpfs[n,"AdjBMI2"], " not found \n")
        next()
        }else{
        if(!file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"])) & !file.exists(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"], ".gz" ))){

          cat ( "running ", eumaxpfs[n,2], eumaxpfs[n,4],eumaxpfs[n,6],eumaxpfs[n,7],eumaxpfs[n,"AdjBMI2"], n,  " \n")
           
           # GenR had leading NULLs and comments, so use read.table instead of fread
          if(eumaxpfs[n,3] %in% c("GenR1", "GenR2")){
            tmp <- read.table(gzfile(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_fn"])), head=T, string=F)
          }else {
            tmp <- as.data.frame(fread(paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_fn"])))
          } 
          cat ( "  data read \n")

            #eumaxpfs$N_SNP_QCed[n] <- nrow(tmp)
            tmp$bp_QCed<-sapply(strsplit(tmp$cpaid,      ":"), `[`, 2)
           
            # re-order chr and position before save for locus zoom
            tmpord <- tmp[which(!is.na(tmp$CHR)),]
            tmpord$chrpos <- paste0(tmpord$CHR, tmpord$bp_QCed)
            #eumaxpfs$N_SNP_QCed_Dupl[n] <- sum(duplicated(tmpord$chrpos ) | duplicated(tmpord$chrpos, fromLast=T))

            tmpord <-tmpord[!duplicated(tmpord$chrpos ) & !duplicated(tmpord$chrpos, fromLast=T), ] 
            #check 
            if(length(names(table(tmpord$bp_QCed==tmpord$POS)))>1 | names(table(tmpord$bp_QCed==tmpord$POS)) != "TRUE") cat( "\n WARNING !!!!! , check BP for ", n,  eumaxpfs$QCed_fn, "\n")
            #eumaxpfs[n,"N_SNP_QCed_noDup"] <- nrow(tmpord)

            #eumaxpfs$N_SNP_QCed_Dupl_cpaid[n] <- sum(duplicated(tmpord$cpaid ) | duplicated(tmpord$cpaid, fromLast=T))
            
            tmpord <-tmpord[!duplicated(tmpord$cpaid ) & !duplicated(tmpord$cpaid, fromLast=T), ] 
            #eumaxpfs[n,"N_SNP_QCed_noDup2"] <- nrow(tmpord)

            tmpord <-tmpord[order(tmpord[,"CHR"], tmpord[,"bp_QCed"]), ] 
            # Save file with chr and position for locus zoom 
            #if(!file.exists(eumaxpfs[n,"QCedgz_fn"]) ) 
            #write.table(tmpord,paste0(eumaxpfs[n,"QCed_p"], "/", eumaxpfs[n,"QCed_noDup"] ), row.names=F, col.names=T, quote=F, sep="\t")
            # Open list of variants that fail AF check 
            

            aff <- read.table(paste0(eumaxpfs[n,"QCed_p"], "/", eumaxpfs[n,"aff_fn"]), head=T, string=F)
            eumaxpfs[n,"N_SNP_QCed_noDup_AFO"] <- length(which(tmpord$cpaid %in% aff$cpaid))
            tmpord_noafo <- tmpord[which(!tmpord$cpaid %in% aff$cpaid ),]
            #eumaxpfs[n,"N_SNP_QCed_noDup_noAFO"] <- nrow(tmpord_noafo)
            cat ( "  saving output file\n")

            write.table(tmpord_noafo, paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"]) , row.names=F, col.names=T, quote=F, sep="\t")

            #noLO <- as.data.frame(fread(paste0(eumaxpfs[n,"QCed_p"], "/", eumaxpfs[n,"noLO_fn"])))
            #table(noLO$SNPID %in% tmpord_noafo$SNPID  )
            #summary(tmpord_noafo$N)

            cat ( "  zipping output file\n")

            #if(!file.exists(paste0(eumaxpfs[n,"QCedgz_fn"], ".gz"))) 
            gzip (paste0(eumaxpfs[n,"QCed_p"], "/",eumaxpfs[n,"QCed_noDupNoAFO_fn"]))

            #tmpord$P-value <- as.numeric(tmpord$P-value)
            # if(!file.exists(gsub("_GWAMA.TBL", "_GWAMA_manh.", eumaxpfs[n,"QCed_fn"]))){
            #   jpeg(gsub("_GWAMA.TBL", "_GWAMA_manh.", eumaxpfs[n,"QCed_fn"]))
            #     manhattan(tmpord[which(tmpord$chr %in% 1:22),],p="P-value", chr="chr", bp="bp", snp="MarkerName" )
            #   dev.off()
            # } 

        }
      }
    }





# Copy adult cleaned women files to MW 

  file.copy(paste0(o_p, "/ADL/FAMILY_cIMT.maxLCF.ADL.W.noBMI.EUR.AL.20240825_2.txt.gz" ), paste0(o_p, "/ADL/FAMILY_cIMT.maxLCF.ADL.MW.noBMI.EUR.AL.20240825_2.txt.gz" ))
  file.copy(paste0(o_p, "/ADL/FAMILY_cIMT.maxLCF.ADL.W.adjBMI.EUR.AL.20240825_2.txt.gz"), paste0(o_p, "/ADL/FAMILY_cIMT.maxLCF.ADL.MW.adjBMI.EUR.AL.20240825_2.txt.gz"))
  file.copy(paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.W.noBMI.EUR.AL.20240825_2.txt.gz" ), paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.MW.noBMI.EUR.AL.20240825_2.txt.gz") )
  file.copy(paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.W.adjBMI.EUR.AL.20240825_2.txt.gz"), paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.MW.adjBMI.EUR.AL.20240825_2.txt.gz"))
  file.copy(paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.W.noBMI.EUR.AL.20240825_2.txt.gz" ), paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.MW.noBMI.EUR.AL.20240825_2.txt.gz") )
  file.copy(paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.W.adjBMI.EUR.AL.20240825_2.txt.gz"), paste0(o_p, "/OAD/FAMILY_cIMT.maxLCF.OAD.MW.adjBMI.EUR.AL.20240825_2.txt.gz"))
