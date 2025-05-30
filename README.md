
# 🧬 eQTL Analysis Pipeline
eQTL (expression quantitative trait loci) analysis identifies genetic variants that influence gene expression levels. 
It's a powerful tool for understanding how genetic variation contributes to differences in gene expression across individuals or tissues. 

This repository contains a reproducible pipeline and scripts for performing **eQTL (expression Quantitative Trait Loci) analysis**, 
integrating genotyping and gene expression data to identify genetic variants that regulate gene expression.

**Expression QTLs (eQTLs)** are genomic loci that explain variation in gene expression levels across individuals. A project typically involves:

- Genotype data (e.g., VCF, PLINK)
- RNA-seq gene expression data (e.g., TPM, counts)
- Covariates (e.g., population structure, batch effects)

We use tools like **Matrix eQTL** to detect **cis-** and **trans-eQTLs**.

## 📁 Project Structure

├── README.md  
├── bash/   
├── scripts/  
