#!/bin/bash

plink=/home/admin/Desktop/./plink

BASE="/home/admin/Desktop/ELIAS/Cohort_PAFIP-Santander/pafip_santander_no_dups"


mkdir /home/admin/Desktop/ELIAS/OUTPUT_$GWAS

cat _pst* > scores.txt

awk '{print $2, $4, $6}' scores.txt > weights.txt

$plink \
  --bfile $BASE   
  --score weights.txt 
  --out my_prs_output


#!/bin/bash

# Directorio donde están los archivos de entrada
INPUT_DIR="/media/admin/Seagate Basic/Elias/Cohort PAFIP - Santander/GWAS/prscs"


# Directorio con el archivo de genotipos (en formato bed/bim/fam)
PLINK_PREFIX="/home/admin/Desktop/ELIAS/Cohort_PAFIP-Santander/pafip_santander_no_dups" 

# Iterar sobre archivos
for FILE in ${INPUT_DIR}/*.prscs; do
    BASENAME=$(basename "$FILE" .prscs)
    mkdir "media/admin/Seagate Basic/Elias/Cohort PAFIP - Santander/GWAS/OUTPUT_PRSCS/${FILE}_out"
    OUTPUT_DIR="/media/admin/Seagate Basic/Elias/Cohort PAFIP - Santander/GWAS/OUTPUT_PRSCS/${FILE}_out"
    
    echo "Procesando $BASENAME..."

    # Ejecutar PRS-CS
    python run_prscs.py \
        --ref_dir=ldblk_1kg_eur \
        --bim_prefix=${PLINK_PREFIX} \
        --sst_file=${FILE} \
        --n_gwas=100000 \
        --out_dir=${OUTPUT_DIR} \
        --out_name=${BASENAME}
    
    # Archivo .txt con los efectos posterior a PRS-CS
    EFFECT_FILE="${OUTPUT_DIR}/${BASENAME}_pst_eff_a1_b0.5_phi1e-2.txt"

    # Calcular el PRS con plink
    plink \
        --bfile ${PLINK_PREFIX} \
        --score ${EFFECT_FILE} 1 2 4 sum \
        --out ${OUTPUT_DIR}/${BASENAME}_prs
done

