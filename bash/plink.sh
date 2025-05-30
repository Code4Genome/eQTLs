
vcf_file <- "VTE_IMFgenos.vcf"  # Change this to your VCF filename
vcf <- read.vcfR(vcf_file)

# STEP 2: Extract the Genotype (GT) data
geno <- extract.gt(vcf, element = "GT", as.numeric = FALSE)

head(geno_numeric)

plink --vcf VTE_IMFgenos.vcf --recode A  --make-bed  --mind 0.02 --out recoded_geno_VTE

plink --vcf VTE_IMFgenos.vcf --freq --out allele_freq

plink --bfile recoded_geno_VTE --missing --out sampleQC

plink --bfile recoded_geno_VTE --het --out sampleQC

plink --bfile recoded_geno_VTE --indep-pairwise 50 5 0.2 --out LD_VTE

plink --bfile recoded_geno_VTE --extract LD_VTE.prune.in --make-bed --out VTE_pruned

plink --bfile VTE_pruned --genome --min 0.125 --out Related

#---------------------------------------------------------------------------------

plink.exe --bfile VTE_pruned --make-bed -out fullQCsample

plink --bfile fullQCsample --missing -out SNPQC

plink --bfile fullQCsample --geno 0.02 --make-bed -out VTE_fullQC

plink--bfile VTE_fullQC --hwe 1e-6  --make-bed -out VTE_full_HWE

plink --bfile VTE_full_HWE --maf 0.01 --make-bed --out VTE_total

plink --bfile VTE_total --pca header tabs --out PCA
