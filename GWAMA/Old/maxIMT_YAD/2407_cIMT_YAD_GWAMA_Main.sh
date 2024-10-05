# DESCIPTION 
#This code complies the differents codes used to run the GWAS meta-analysis of GWAMA 
#path on git /y/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Scripts 
hostname=cat /proc/sys/kernel/hostname
host=rigenep1
if [ $(hostname) = rigenep1 ] | [ $(hostname) = rigenep2 ] | [ $(hostname) = ristatp25 ] | [ $(hostname) = ristatp23 ] | [ $(hostname) = ristatp21 ] | [ $(hostname) = ristatp19 ] | [ $(hostname) = ristatp17 ];then 
echo TRUE 
m_p=/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI
metal=/genetics/Programs/Linux-metal/generic-metal/metal 
fi 

m_p=/genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI
metal=/genetics/Programs/Linux-metal/generic-metal/metal 
o_p=$m_p/Results
s_p=$m_p/Scripts
eq2_p=$s_p/EasqyQC2

# Prep files for EasyQC 
  # YFS 

  
    cd /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/YFS/Raw
  
    abcdfs=("EGG_cIMT.maxLCF.YAD.MW.noBMI.YFS.EUR.LPL.231001.txt.gz")
  
    for filei in ${abcdfs[@]};do 
    
      echo ${abcdfs[@]}
    
      #cp /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/$filei.gz /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei.gz
    
      #gunzip /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei.gz
    
      #sed -i.bak 's/ /\t/g' /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei
      
      #gzip /genetics/MixedStudies/Projects/2211_cIMT_GWAMA_PHRI/Results/GWAS/ABCD/Raw/Modified/$filei
    
    done 

    $eq2_p/ABCD/ADO_CCA/EGG_cIMT_ABCD_ADO_maxCCA_QC_quant.ecf



  
# Run GWAMA for maxIMT 
  cd $o_p/GWAMA/maxIMT_ADO
  $metal $s_p/2211_cIMT_GWAMA_S2_RunMetal_MW_ADO_maxIMT_NoAdjBMI.metal > $o_p/GWAMA/maxIMT_ADO/METAL_MW_maxIMT_NoAdjBMI_240619.log
