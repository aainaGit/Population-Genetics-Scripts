#!/bin/bash
##--------------------------------------------------------------------------------------------------------------------------------------------------
##Script to manipulate and analyse VCF files from GBS for population structure and genetic diversity indices
##Get into the right conda environment where all software/tools required are already loaded
##Conda active bioinformatics

# Removing unwanted sequencing info attached to file names e.g -plate-- including all blanks
sed -e 's/-FeralHemp22Plate-.*//' -e 's/-Plate-.*//' -e '/^Blank/d' samples.txt > samples2.txt

##----------MAKE A COPY OF THE ORIGINAL VCF FILE FOR DOWNSTREAM ANALYSIS----------------------------------------------------------------------------
cp SNPs.mergedAll.vcf SNPs.mergedAll_copy.vcf

##----------VIEW SAMPLE NAMES
bcftools query -l SNPs.mergedAll.vcf 

#----------BGZIP VCF FILE FOR EASE OF USE/OPERATION 
bgzip SNPs.mergedAll.vcf
bcftools index SNPs.mergedAll.vcf.gz

##----------VIEW CHR NUMBER AND EDIT TO NUMERIC ALSO SET CONTIGS TO NUMERIC FOR EASY DELETION  ---------------------------------------------------

bcftools query -f '%CHROM\n'  SNPs.mergedAll.vcf.gz 
#correct header error if you get one like this [W::bcf_hdr_check_sanity] PL should be declared as Number=G 
##Unzip gz file file, correct sanity and gzip again
#bgzip -d SNPs.mergedAll.vcf.gz
#sed 's/PL,Number=./PL,Number=G/g' SNPs.mergedAll.vcf > SNPs.merged.vcf
#bgzip SNPs.merged.vcf
#bcftools index SNPs.merged.vcf.gz

#bcftools query -f '%CHROM\n'  SNPs.merged.vcf.gz | sort | uniq > chromosomes.txt
#bcftools annotate --rename-chrs chromosomes.txt SNPs.merged.vcf.gz -Oz -o SNP.merged2.vcf.gz
##Use bcftools to filter out chromosome 11 (renamed unmapped contigs)and create a new VCF file
#bcftools view -t ^11 -o SNP.merged3.vcf.gz -O v SNP.merged2.vcf.gz
#Check to see if unmapped contigs have been removed make ABOVE OUTPUT ONLY VCF no GZ extension
#bcftools query -f '%CHROM\n' SNP.merged3.vcf.gz | sort | uniq

##-----------------------------------------------------------------------------------------------------------------------------------------------
#Uncompress the gz.vcf file
#bgzip -d SNP.merged3.vcf.gz
##make ped and map file from vcf file
#plink --vcf SNP.merged3.vcf.gz --recode --nonfounders --allow-no-sex --out SNPs.merged

##make bed bim and fam files from plink ped and map file
#plink --file SNPs.merged --make-bed --nonfounders --allow-no-sex --allow-extra-chr --out SNPs.merged1

#Filter with QC options
##missingness per SNP: --geno
##missingness per Individual: --mind
##minor allele frequency: --maf
#run each command individual to avoid plink eliminating all data
#plink --bfile SNPs.merged1 --geno 0.1 --allow-extra-chr --make-bed --out SNPs.merged2
#plink --bfile SNPs.merged2 --maf 0.05 --allow-extra-chr --make-bed --out SNPs.merged3
#plink --bfile SNPs.merged3 --mind 0.1 --allow-extra-chr --make-bed --out SNPs.merged4

#recode SNPs.merged4.bed file back to vcf (use this bed file for admixture analysis)
#plink --bfile SNPs.merged4 --recode vcf-iid --allow-extra-chr --out SNPs.merged4

##RENAME SAMPLES IN SNPs.merged4 and RECODE TO BED FILE FOR ADMIX
# Remove the duplicate sample_name from the input VCF file and output the result to a new VCF file
#bcftools view --samples ^GWAS-435-Arl2022-GWAS-Plate4-A11,GWAS-435-Arl2022-GWAS-Plate4-A10 -o SNPs.merged5.vcf SNPs.merged4.vcf
#remove plate info from all sample names before renaming using header
#bcftools reheader -s old_new_rename.txt -o SNPs.merged6.vcf SNPs.merged5.vcf
#rename again to final names
#bcftools reheader -s sampleIDs_to_rename.txt -o SNPs.merged7.vcf SNPs.merged6.vcf
#bcftools query -l SNPs.merged7.vcf
#plink --vcf SNPs.merged7.vcf --make-bed --double-id --out SNPs.merged7
#extract sample names for metadata creation/same for bedfile from which it was created
#bcftools query -l SNPs.merged7.vcf > checkfiles.txt
#cut -c 1-2 checkfiles.txt > checkfiles2.txt

##MOVE TO R
##PCA analysis -----------------------------------------------------------------------------------------------------------------------------------------
##PCA analysis STEP1 (LD thining)
#plink --vcf SNPs.merged7.vcf --double-id --set-missing-var-ids @:# --indep-pairwise 50 2 0.2 --out SNPs.QCmerged7
##PCA analysis STEP 2
#plink --vcf SNPs.merged7.vcf --double-id --set-missing-var-ids @:# --extract SNPs.QCmerged7.prune.in --pca --make-bed --out SNPs.QCmerged7.prune_pca
##For PCA plot (move to R/Rstudio)


##Admixture analysis-------------------------------------------------------------------------------------------------------------------------------------
#for K in 1 2 3 4 5 6 7 8 9 10 11 12; \
#do admixture --cv SNPs.merged7.bed $K | tee log${K}.out; done > admixresult.txt 2> admixresult.stderr 

###Filter out admixed individuals from VCF file and rerun PCA in Plink-----------------------------------------------------------------------------------
# Ensure vcftools is installed and available
#vcftools --vcf SNPs.merged7.vcf --remove na_remove@K7.txt --recode --out SNPs.merged@K7
##PCA analysis STEP1 (LD thining)
#plink --vcf SNPs.merged@K7.recode.vcf --double-id --set-missing-var-ids @:# --indep-pairwise 50 2 0.2 --out SNPs.merged@K7
##PCA analysis STEP 2
#plink --vcf SNPs.merged@K7.recode.vcf --double-id --set-missing-var-ids @:# --extract SNPs.merged@K7.prune.in --pca --make-bed --out SNPs.merged@K7.prune_pca
##For PCA plot (move to R/Rstudio)


#CORRECTING ERROR IN SAMPLE_IDs
#bcftools reheader -s renameIDs_finalVcfs@K7.txt -o SNPs.merged@K7Final.vcf SNPs.merged@K7.recode.vcf
#bcftools query -l SNPs.merged@K7Final.vcf


#Generate distance matrix for NJ tree 
#plink --vcf SNPs.merged@K7Final.vcf --make-bed --double-id --out SNPs.merged@K7Final 
#plink --bfile SNPs.merged@K7Final --distance-matrix --out ForGen

#########
##adjusting NJ Tree
#vcftools --vcf SNPs.merged@K7Final.vcf --remove takeout.txt --recode --out SNPs.merged@K7NJ 

### For Bimbo New dist matrix for another file
#plink --vcf bimbo312RENAMED_Fil.vcf --make-bed --double-id --out bimbo312_Fil 
#plink --bfile bimbo312_Fil --distance-matrix --out Bimb
# Ensure vcftools is installed and available
vcftools --vcf bimbo312RENAMED_Fil.vcf --remove K7data_na.txt --recode --out bimboK3
##PCA analysis STEP1 (LD thining)
plink --vcf bimboK3.recode.vcf --double-id --set-missing-var-ids @:# --indep-pairwise 50 2 0.2 --out bimboK3SNP
##PCA analysis STEP 2
plink --vcf bimboK3.recode.vcf --double-id --set-missing-var-ids @:# --extract bimboK3SNP.prune.in --pca --make-bed --out bimboK3SNP.prune_pca
##For PCA plot (move to R/Rstudio)


