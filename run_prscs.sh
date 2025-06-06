#!/bin/bash

# Paths
PRSCSPATH="/home/admin/Desktop/PRScs"
PRSCSPY="${PRSCSPATH}/PRScs.py"

REF_DIR="/home/admin/Desktop/ELIAS/ldblk_1kg_eur"
BIM_PREFIX="/home/admin/Desktop/ELIAS/Cohort_PAFIP-Santander/pafip_santander_no_dups_nomhc"
GWAS_DIR="/media/admin/Seagate Basic/Elias/Cohort PAFIP - Santander/GWAS/prscs"
OUT_BASE="/media/admin/Seagate Basic/Elias/Cohort PAFIP - Santander/GWAS/OUTPUT_PRSCS"


# PLINK binary path (ajústalo si necesario)
PLINK=/home/admin/Desktop/./plink

# Iterar sobre todos los archivos .txt de GWAS
for SST_FILE in "${GWAS_DIR}"/*prscs; do
    BASENAME=$(basename "$SST_FILE" _prscs)
    OUT_DIR="${OUT_BASE}/${BASENAME}"
    mkdir -p "$OUT_DIR"
    
    N_GWAS=$(wc -l < "$SST_FILE")
    
    echo "Procesando $BASENAME..."

    # Ejecutar PRS-CS
    python3 "$PRSCSPY" \
        --ref_dir="$REF_DIR" \
        --bim_prefix="$BIM_PREFIX" \
        --sst_file="$SST_FILE" \
        --n_gwas="$N_GWAS" \
        --out_dir="$OUT_DIR"


    cat ${OUT_DIR}/*eff* > ${OUT_DIR}/${BASENAME}_scores.txt
    awk '{print $2, $4, $6}' ${OUT_DIR}/${BASENAME}_scores.txt > ${OUT_DIR}/${BASENAME}_weights.txt
    
    EFFECT_FILE="${OUT_DIR}/${BASENAME}_weights.txt"
    
    if [[ -f "$EFFECT_FILE" ]]; then
        echo "Calculando PRS con PLINK para $BASENAME..."

        $PLINK \
            --bfile "$BIM_PREFIX" \
            --score "$EFFECT_FILE" \
            --out "${OUT_DIR}/${BASENAME}_prs"
    else
        echo "ERROR: No se encontró archivo de efectos para $BASENAME"
    fi
done


#plink=/home/admin/Desktop/./plink
#cat *.txt > GWAS_scores.txt
#awk '{print $2, $4, $6}' GWAS_scores.txt > GWAS_weights.txt
#$plink --bfile /home/admin/Desktop/ELIAS/Cohort_PAFIP-Santander/pafip_santander_no_dups_nomhc --score GWAS_weights.txt --out GWAS_prs
