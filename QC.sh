#!/bin/bash

plink=/home/admin/Desktop/./plink

# Define base file
BASE="/home/admin/Desktop/ELIAS/Cohort_PAFIP-Santander/pafip_santander_merged"

$plink --bfile $BASE --recode vcf --out $BASE


# Step 1: Missingness
$plink --pafip_santander_merged $BASE --missing  



# Step 2: Filter SNPs and individuals by missingness
$plink --pafip_santander_merged $BASE --geno 0.02 --make-bed --out qc_step1
$plink --pafip_santander_merged qc_step1 --mind 0.02 --make-bed --out qc_step2

# Step 3: Filter Indels
$plink --bfile qc_step2 --snps-only just-acgt --make-bed --out qc_step3

# Step 4.1
$plink --bfile qc_step3 --update-sex sex_info.txt --make-bed --out qc_step3_sex

# Step 4.2: Eliminate duplicated SNPs
$plink --bfile qc_step3_sex --list-duplicate-vars
$plink --bfile qc_step3_sex --exclude plink.dupvar --make-bed --out qc_step4

# Step 5: Eliminate IDs duplicados
$plink --bfile qc_step4 --check-sex --out sex_check
grep "PROBLEM" sex_check.sexcheck | awk '{print $1, $2}' > sex_discrepancy.txt
$plink --bfile qc_step4 --remove sex_discrepancy.txt --make-bed --out qc_step5

# Step 6: Assign IDs 
$plink --bfile qc_step5 --set-missing-var-ids @:# --make-bed --out qc_step6

# Step 7: Filter SNPs with low MAF
$plink --bfile qc_step6 --maf 0.01 --make-bed --out STEP_FINAL/qc_final

# Final files: qc_final.bed, qc_final.bim, qc_final.fam
echo "QC done. Output file: qc_final"

