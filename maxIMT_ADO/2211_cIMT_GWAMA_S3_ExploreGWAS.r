library(qqman)

# m<- "/home/lamria/avr/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results

m_p <- "/home/lamria/avr/Projects/byStudy/EGGC/2211_cIMT_GWAMA/Results/GWAMA/"
#m_p <- "/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAMA/maxIMT_ADO/2211_cIMT_GWAMA_PHRI/Results/GWAMA/"
kg_p <- "/home/lamria/avr/Data/1000Genomes/20130502_Release/Genotypes/PLINK/Sample_Subsets/Europeans/SNP_Subsets/mac10/"
pheno <- "maxIMT_ADO" 

d <- read.table(paste0(m_p, pheno, "/METAL_MW_maxIMT_NoAdjBMI_240619.TBL"), head=T, string=F)
nrow(d)
d$chr<-sapply(strsplit(d$MarkerName,":"), `[`, 1)
d$bp<-sapply(strsplit(d$MarkerName,":"), `[`, 2)
d$chr<- as.numeric(d$chr)
d$bp <- as.numeric(d$bp)

jpeg(paste0(m_p, pheno, "/METAL_MW_maxIMT_NoAdjBMI_manha_240619.jpeg"))
manhattan(d[which(d$chr %in% 1:22),],p="P.value", chr="chr", bp="bp", snp="MarkerName" )
dev.off()


# get SNPnames from TOPMED 

sig <- d[d$P.value< 10^-8, ]
sig$chrpos <- paste0(bim$chr, ":", bim$bp)
sig$bpp1 <- sig$bp + 1
sig$chr2<- paste0("chr", sig$chr)
write.table(sig[,c("chr2", "bp", "bpp1", "MarkerName")], paste0(m_p, pheno, "/METAL_MW_maxIMT_NoAdjBMI_manha_240619_sig_for_liftover"), row.names=F, col.names=F, quote=F)

lo <- read.table(paste0(m_p, pheno, "/METAL_MW_maxIMT_NoAdjBMI_manha_240619_sig_hg19.bed"), head=F, string=F)[,c(1,2,4)]
colnames(lo) <- c("chr", "pos_19", "rsid_metal")
lo$chr <- gsub("chr", "",lo$chr )
lo$MarkerName <- lo$rsid_metal

sig2 <- merge(sig, lo )
nrow(sig2)
sig2$chrpos19 <- paste0(sig2$chr, ":",sig2$pos_19)


# Open 1KG SNP names
bim <- read.table(paste0(kg_p, "//European1000G_KeepAlleleOrder_chr1..22_DupVars_removed_mac10.bim"), head=F, string=F )




bim$chrpos19 <- paste0(bim$V1, ":", bim$V4)
colnames(bim)<- c("chr", "rsid_kg", "cm", "pos_19", "a1_kg", "a2_kg", "chrpos19")

m<- merge(sig2, bim)
nrow(m)

ms<- m[,c("rsid_kg", "P.value")]
colnames(ms) <- c("SNP", "P")

write.table(ms, paste0(m_p, "/maxIMT_ADO/METAL_MW_maxIMT_NoAdjBMI_manha_240619_clump"), row.names=F, col.names=T, quote=F)

plink \
--bfile /home/lamria/avr/Data/1000Genomes/20130502_Release/Genotypes/PLINK/Sample_Subsets/Europeans/SNP_Subsets/mac10/European1000G_KeepAlleleOrder_chr1..22_DupVars_removed_mac10 \
--clump /home/lamria/avr/Projects/byStudy/EGGC/2211_cIMT_GWAMA/Results/GWAMA/maxIMT_ADO/METAL_MW_maxIMT_NoAdjBMI_manha_240619_clump \
--out /home/lamria/avr/Projects/byStudy/EGGC/2211_cIMT_GWAMA/Results/GWAMA/maxIMT_ADO/METAL_MW_maxIMT_NoAdjBMI_240619_clumps


