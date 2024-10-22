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
        eumaxp[[i]][[j]][[k]] <- list()  
        }
      }
    }


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
  fts[fts$AgeGr=="OAD","ifn"] <- gsub("GWAMA.TBL","GWAMA1.TBL", fts[fts$AgeGr=="OAD","ifn"])
  fts[fts$AgeGr=="OAD","ifn"] <- gsub("240926","241005", fts[fts$AgeGr=="OAD","ifn"])

# Keep only a few columns
  ftss <- fts[!is.na(fts$N_SNP),c(1:7, 10:11)]
  ftss[,ncol(ftss)] <- paste(ftss[,1], ftss[,2] , ftss[,3] , sep="_")

# Open files tables 

  lft <- read.delim(paste0(m_p, "/Results/241004_GWAS_Files_noAFO.tbl"), head=T, string=F)
  eumaxpfs <- unique(lft[which(lft$AdjBMI2=="noBMI" & lft$sexGroup == "MW" & lft$in_OAD_maxLCFp_EUR_MW_noBMI |lft$in_ADL_maxLCFp_EUR_MW_noBMI | lft$in_YAD_maxLCFp_EUR_MW_noBMI |lft$in_ADO_maxLCFp_EUR_MW_noBMI ), c("Study","pop", "pop2", "AgeGroup","sexGroup", "AdjBMI","AdjBMI2","cIMT_measure","MetaA_p")])
  nrow(eumaxpfs)

  unique(eumaxpfs[,c("pop2", "AgeGroup","sexGroup", "AdjBMI2")])

,
# Open gwama data and save for locus zoom 
  for (n in 7:nrow(fts)){ #which(fts[,1]=="OAD" 
    i <- fts[n,5]
    j <- fts[n,6]
    k <- fts[n,7]
    if(!file.exists(fts[n,8])){
      cat ( "file", fts[n,4], " not found \n")
      next()
    }else{
      if(file.exists(paste0(fts[n,9], ".gz" ))) {
        cat ( "running ", fts[n,4], " \n")
        eumaxp[[i]][[j]][[k]] <- as.data.frame(fread(fts[n,8]))
        fts[n,11] <- nrow(eumaxp[[i]][[j]][[k]])
        names(eumaxp)[i] <- fts[n,1]
        names(eumaxp[[i]])[j] <- fts[n,2]
        names(eumaxp[[i]][[j]])[k] <- fts[n,3]
       } else if (file.exists(fts[n,9]) & !file.exists(paste0(fts[n,9], ".gz" ))){
        cat ( "running ", fts[n,4], " \n")
        eumaxp[[i]][[j]][[k]] <- as.data.frame(fread(fts[n,9]))
        fts[n,11] <- nrow(eumaxp[[i]][[j]][[k]])
        names(eumaxp)[i] <- fts[n,1]
        names(eumaxp[[i]])[j] <- fts[n,2]
        names(eumaxp[[i]][[j]])[k] <- fts[n,3]
        gzip (fts[n,9])
       } else if(!file.exists(fts[n,9]) & !file.exists(paste0(fts[n,9], ".gz" ))){
        cat ( "running ", fts[n,4], " \n")
        eumaxp[[i]][[j]][[k]] <- as.data.frame(fread(fts[n,8]))
        names(eumaxp)[i] <- fts[n,1]
        names(eumaxp[[i]])[j] <- fts[n,2]
        names(eumaxp[[i]][[j]])[k] <- fts[n,3]
        fts[n,10] <- nrow(eumaxp[[i]][[j]][[k]])
        tmp <- eumaxp[[i]][[j]][[k]]
        #tmp$chr<-sapply(strsplit(tmp$MarkerName,":"), `[`, 1)
        #tmp$bp<-sapply(strsplit(tmp$MarkerName,      ":"), `[`, 2)
        #tmp$chr<- as.numeric(tmp$chr)
        #tmp$bp <- format(as.numeric(tmp$bp), scientific = FALSE) 
        
        # re-order chr and position before save for locus zoom
        tmpord <- tmp[which(!is.na(tmp$chr)),]
        tmpord$chrpos <- paste0(tmpord$CHR, tmpord$POS)
        tmpord <-tmpord[!duplicated(tmpord$chrpos ) & !duplicated(tmpord$chrpos, fromLast=T), ] 
        tmpord <-tmpord[order(tmpord[,"CHR"], tmpord[,"POS"]), ] 
        eumaxp[[i]][[j]][[k]] <- tmpord
        fts[n,11] <- nrow(tmpord)
        # Save file with chr and position for locus zoom 
        #if(!file.exists(fts[n,9]) )
        #write.table(tmpord, fts[n,9] , row.names=F, col.names=T, quote=F, sep="\t")
        #if(!file.exists(paste0(fts[n,9], ".gz"))) 
        #gzip (fts[n,9])
        #tmpord$P-value <- as.numeric(tmpord$P-value)
        # if(!file.exists(gsub("_GWAMA.TBL", "_GWAMA_manh.", fts[n,8]))){
        #   jpeg(gsub("_GWAMA.TBL", "_GWAMA_manh.", fts[n,8]))
        #     manhattan(tmpord[which(tmpord$chr %in% 1:22),],p="P-value", chr="chr", bp="bp", snp="MarkerName" )
        #   dev.off()
        # } 
      }
    }
  }



  for (i in 1:4){
    for (j in 1:3 ){  
      for (k in 1:2){  
        cat( i, j, k, colnames(eumaxp[[i]][[j]][[k]]), "\n\n\n\n" )
        }
      }
    }


# Keep only significant loci 
  sig <- eumaxp

  for (i in 1:4){
    for (j in 1 ){  
      for (k in 1:2){  
        cat (names(eumaxp)[i] , names(eumaxp[[i]])[j] , names(eumaxp[[i]][[j]])[k] , "\n" ) 
        colnames(sig[[i]][[j]][[k]]) <- gsub("P.value","P-value", colnames(sig[[i]][[j]][[k]]))
        sig[[i]][[j]][[k]] <- sig[[i]][[j]][[k]][sig[[i]][[j]][[k]][,"P-value",] <= 0.05,]
        sig[[i]][[j]][[k]] <- sig[[i]][[j]][[k]][order(sig[[i]][[j]][[k]][,"P-value"]),]

      }
    }
  }

  sig2 <- sig
   for (i in 1:4){
      for (j in 1 ){  
        for (k in 1:2){  
          colnames(sig2[[i]][[j]][[k]]) [! colnames(sig2[[i]][[j]][[k]]) %in% c("MarkerName","chr","bp","chrpos")] <- paste0(colnames(sig2[[i]][[j]][[k]]) [! colnames(sig2[[i]][[j]][[k]]) %in% c("MarkerName","chr","bp","chrpos")], "_",names(sig2)[i],"_",names(sig2[[i]])[[j]],"_",names (sig2[[i]][[j]])[[k]])
          sig2[[i]][[j]][[k]] <- sig2[[i]][[j]][[k]][order(sig2[[i]][[j]][[k]]$"P-value" ),]
        }
      }
   }
  
  head(sig2[[3]][[1]][[1]])
  head(sig2[[4]][[1]][[1]])
  
  min(sig2[[1]][[1]][[1]]$"P-value" )


summary(sig[[4]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    491     491     864    1295    2114    3469
summary(sig[[3]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    657    1779    2114    2395    2771    4550

summary(sig[[2]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    302    1575    1575    1383    1877    1877
summary(sig[[1]][[1]][[1]][, "N"])
#   Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
#    241    1650    2752    2729    3875    5284




# EUR MW NoBMI 



# Keep only significant loci 
  eumaxpa <- eumaxp

  for (i in 1:4){
    for (j in 1 ){  
      for (k in 1:2){  
        cat (names(eumaxpa)[i] , names(eumaxpa[[i]])[j] , names(eumaxpa[[i]][[j]])[k] , "\n" ) 
        colnames(eumaxpa[[i]][[j]][[k]]) <- gsub("P.value","P-value", colnames(sig[[i]][[j]][[k]]))
          colnames(eumaxpa[[i]][[j]][[k]]) [! colnames(eumaxpa[[i]][[j]][[k]]) %in% c("MarkerName")] <- paste0(colnames(eumaxpa[[i]][[j]][[k]]) [! colnames(eumaxpa[[i]][[j]][[k]]) %in% c("MarkerName","chr","bp","chrpos")], "_",names(eumaxpa)[i],"_",names(eumaxpa[[i]])[[j]],"_",names (eumaxpa[[i]][[j]])[[k]])


      }
    }
  }


pci <- c("P-value_ADO_MW_noBMI","P-value_YAD_MW_noBMI","P-value_ADL_MW_noBMI","P-value_OAD_MW_noBMI")
eci <- c("Effect_ADO_MW_noBMI","Effect_YAD_MW_noBMI", "Effect_ADL_MW_noBMI","Effect_OAD_MW_noBMI")
epci <-  c("MarkerName","Effect_ADO_MW_noBMI","P-value_ADO_MW_noBMI","Effect_YAD_MW_noBMI","P-value_YAD_MW_noBMI", "Effect_ADL_MW_noBMI","P-value_ADL_MW_noBMI","Effect_OAD_MW_noBMI","P-value_OAD_MW_noBMI")



  MW_noBMI_2 <- merge(eumaxpa[[1]][[1]][[1]], eumaxpa[[2]][[1]][[1]], all=T)
  MW_noBMI_2 <- merge(MW_noBMI_2           , eumaxpa[[3]][[1]][[1]], all=T)
  MW_noBMI_2 <- merge(MW_noBMI_2           , eumaxpa[[4]][[1]][[1]], all=T)

MW_noBMI_2a <- MW_noBMI_2[order(MW_noBMI_2$"P-value_OAD_MW_noBMI"),]
write.csv(MW_noBMI_2a[MW_noBMI_2a$"P-value_OAD_MW_noBMI"<0.0005,epci], paste0(rgma, "/OAD_Sig_loci_across_Ages.csv"), row.names=F, quote=F)


MW_noBMI_2b<- MW_noBMI_2[order(MW_noBMI_2$"P-value_ADL_MW_noBMI"),]
write.csv(MW_noBMI_2b[which(MW_noBMI_2$"P-value_ADL_MW_noBMI"<0.0005),epci], paste0(rgma, "/ADL_Sig_loci_across_Ages.csv"), row.names=F, quote=F)


MW_noBMI_2c<- MW_noBMI_2[order(MW_noBMI_2$"P-value_YAD_MW_noBMI"),]
write.csv(MW_noBMI_2c[which(MW_noBMI_2$"P-value_YAD_MW_noBMI"<0.005),epci], paste0(rgma, "/YAD_Sig_loci_across_Ages.csv"), row.names=F, quote=F)

MW_noBMI_2d<- MW_noBMI_2[order(MW_noBMI_2$"P-value_ADO_MW_noBMI"),]
write.csv(MW_noBMI_2d[which(MW_noBMI_2$"P-value_ADO_MW_noBMI"<0.005),epci], paste0(rgma, "/ADO_Sig_loci_across_Ages.csv"), row.names=F, quote=F)




colnames(MW_noBMI_2) <- gsub("_ADO_MW_noBMI_ADO_MW_noBMI", "_ADO_MW_noBMI", colnames(MW_noBMI_2))
colnames(MW_noBMI_2) <- gsub("_ADL_MW_noBMI_ADL_MW_noBMI", "_ADL_MW_noBMI", colnames(MW_noBMI_2))

colnames(MW_noBMI_2) <- gsub("_OAD_MW_noBMI_OAD_MW_noBMI", "_OAD_MW_noBMI", colnames(MW_noBMI_2))

colnames(MW_noBMI_2) <- gsub("_YAD_MW_noBMI_YAD_MW_noBMI", "_YAD_MW_noBMI", colnames(MW_noBMI_2))

