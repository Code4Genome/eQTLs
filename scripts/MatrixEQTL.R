library(MatrixEQTL)
library(devtools)
library(tidyverse)


# Genotype file name
SNP_file_name = paste("geno_case_2.txt", sep="\t");
# Gene expression file name
expression_file_name = paste("exp_cases.txt", sep="\t");
# Covariates file name
# Set to character() for no covariates
covariates_file_name = paste("cov_cases.txt", sep="\t");

# Output file name
output_file_name_cis = tempfile();

pvOutputThreshold_cis = 1;

# Distance for local gene-SNP pairs

# Error covariance matrix
# Set to numeric() for identity.
errorCovariance = numeric();

# errorCovariance = read.table("Sample_Data/errorCovariance.txt");

cisDist = 1e7;

## Load genotype data

snps = SlicedData$new();
snps$fileDelimiter = "\t";      # the TAB character
snps$fileOmitCharacters = "NA"; # denote missing values;
snps$fileSkipRows = 1;          # one row of column labels
snps$fileSkipColumns = 1;       # one column of row labels
snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name);

## Load gene expression data

gene = SlicedData$new();
gene$fileDelimiter = "\t";      # the TAB character
gene$fileOmitCharacters = "NA"; # denote missing values;
gene$fileSkipRows = 1;          # one row of column labels
gene$fileSkipColumns = 1;       # one column of row labels
gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name);

## Load covariates

cvrt = SlicedData$new();
cvrt$fileDelimiter = "\t";      # the TAB character
cvrt$fileOmitCharacters = "NA"; # denote missing values;
cvrt$fileSkipRows = 1;          # one row of column labels
cvrt$fileSkipColumns = 1;       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name);
}

for( sl in 1:length(gene) ) {
  mat = gene[[sl]];
  mat = t(apply(mat, 1, rank, ties.method = "average"));
  gene[[sl]] = mat;
}
rm(sl, mat);


maf.list = vector('list', length(snps))
for(sl in 1:length(snps)) {
  slice = snps[[sl]];
  maf.list[[sl]] = rowMeans(slice,na.rm=TRUE)/2;
  maf.list[[sl]] = pmin(maf.list[[sl]],1-maf.list[[sl]]);
}
maf = unlist(maf.list)

## Look at the distribution of MAF
# hist(maf[maf<0.1],seq(0,0.1,length.out=100))

cat('SNPs before filtering:',nrow(snps))
# snps$RowReorderSimple(maf>0.15);
snps$RowReorder(maf>0.05);
cat('SNPs after filtering:',nrow(snps))

#----------------------------------------------------------------------------------#####

## Run the analysis
snpspos = read.table("SNP_Position_GRCh38.txt", header = TRUE, stringsAsFactors = FALSE);
genepos = read.table("Gene_Position_circ.txt", header = TRUE, stringsAsFactors = FALSE);

me = Matrix_eQTL_main(
  snps = snps,
  gene = gene,
  cvrt = cvrt,
  useModel = modelLINEAR_CROSS,
  errorCovariance = errorCovariance,
  verbose = TRUE,
  output_file_name.cis = output_file_name_cis,
  pvOutputThreshold.cis = pvOutputThreshold_cis, # Significance threshold
  pvOutputThreshold = 0,  # No trans-eQTL analysis
  pvalue.hist = "qqplot",
  snpspos = snpspos,  # SNP genomic positions
  genepos = genepos,
  cisDist = cisDist,  # 1 Mb for cis-eQTLs
  min.pv.by.genesnp = TRUE,
  noFDRsaveMemory = FALSE);

unlink(output_file_name_cis);

## Results:

cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');
cat('Detected eQTLs:', '\n');
show(me$cis$eqtls)


## Plot the Q-Q plot of local and distant p-values
plot(me)

#Save results as txt files
results_cis_eqtls = (me$cis$eqtls)


# Define significance threshold

write.table(results_cis_eqtls, file = "cis_eqtls_cases.txt", quote = FALSE,  sep = "\t",
            row.names = TRUE, col.names = TRUE)


head(results_cis_eqtls)


# Filter for significant cis-eQTLs (e.g., FDR < 0.05)

significant_cis_eqtls = results_cis_eqtls[results_cis_eqtls$FDR < 0.05, ]
significant_cis_eqtls

write.table(significant_cis_eqtls, file = "cis_eqtls_cases_significant.txt", quote = FALSE,  sep = "\t",
            row.names = TRUE, col.names = TRUE)
