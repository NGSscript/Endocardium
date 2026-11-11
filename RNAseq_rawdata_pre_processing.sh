#!/bin/bash
# ============================================================================
# Bulk RNA-seq processing pipeline
# ============================================================================

# ------------------------- path -------------------------
RAW_DIR="/path"
WORK_DIR="/path"
SAMPLE_INFO="/path/SAMPLE_INFO.txt" #A comma-delimited .txt file with six columns: sample ID(sample identifier), group(treatment), sequencing strategy(e.g., PE for paired-end), batch ID(integer), reads1 and reads2(FASTQ paths). Each row represents one RNA-seq sample.
fastqc_BIN="/path/fastqc"
fastq_screen_BIN="/path/fastq_screen"
FASTQ_SCREEN_CONF="/path/fastq_screen.mus.conf" 
fastp_BIN="/path/fastp"
hisat2_BIN="/path/hisat2"    
samtools_BIN="/path/samtools"
sambamba_BIN="/path/sambamba"
samblaster_BIN="/path/samblaster" 
Rscript_BIN="/path/Rscript" 
HISAT2_INDEX="/path/hisat2/Mus_musculus"
GTF="/path/genes.gtf"
QUANT_R="/path/Quantification_featurecounts.R"
MERGE_PY="/path/Merge_sample_quantification.py"   
Rscript_BIN="/path/Rscript" 
python_BIN="/path/python"
THREADS=16

# -------------------------sample-------------------------
SAMPLES=(sample1 sample2 sample3) 

# ----------------------------------------------------------------------------
for SAMPLE in "${SAMPLES[@]}"; do
    echo "[$(date +%F\ %T)] ==================== $SAMPLE ===================="

    QC="$WORK_DIR/$SAMPLE/01_reads_processing"
    MAP="$WORK_DIR/$SAMPLE/02_mapping_bam_processing"
    MAT="$WORK_DIR/$SAMPLE/03_matrix_processing"
    mkdir -p "$QC/fastqc" "$QC/fastq_screen" "$QC/fastp" "$MAP" "$MAT"

    r1="$RAW_DIR/${SAMPLE}_1.fq.gz"
    r2="$RAW_DIR/${SAMPLE}_2.fq.gz"
    if [[ ! -f "$r1" || ! -f "$r2" ]]; then
        echo "ERROR: fail ($SAMPLE): $r1 / $r2" >&2
        exit 1
    fi

    # ---------------------- 1. QC ----------------------
    "$fastqc_BIN" -o "$QC/fastqc" -t "$THREADS" "$r1" "$r2" \
        &> "$QC/fastqc/fastqc.log"

    # ---------------------- 2. fastq_screen ----------------------
    "$fastq_screen_BIN" --force --conf "$FASTQ_SCREEN_CONF" \
        --threads "$THREADS" --subset 1000000 \
        "$r1" "$r2" \
        --outdir "$QC/fastq_screen" \
        2> "$QC/fastq_screen/fastq_screen.log"

    # ---------------------- 3. fastp  ----------------------
    "$fastp_BIN" --thread "$THREADS" -z 4 \
        --in1 "$r1" --in2 "$r2" \
        --trim_front1 2 --trim_tail1 0 \
        --trim_front2 2 --trim_tail2 0 \
        --qualified_quality_phred 20 --unqualified_percent_limit 50 \
        --cut_right --cut_window_size 4 --cut_mean_quality 20 \
        --low_complexity_filter --trim_poly_x --trim_poly_g \
        --overrepresentation_analysis --detect_adapter_for_pe --correction \
        --length_required 20 \
        -j "$QC/fastp/${SAMPLE}.stat.json" \
        -h "$QC/fastp/${SAMPLE}.stat.html" \
        --out1 "$QC/fastp/${SAMPLE}_R1.fastq.gz" \
        --out2 "$QC/fastp/${SAMPLE}_R2.fastq.gz" \
        2> "$QC/fastp/fastp.log"

    # ---------------------- 4. HISAT2 fixmate ----------------------
    "$hisat2_BIN" -p "$THREADS" -x "$HISAT2_INDEX" \
        -1 "$QC/fastp/${SAMPLE}_R1.fastq.gz" \
        -2 "$QC/fastp/${SAMPLE}_R2.fastq.gz" \
        --new-summary 2> "$MAP/${SAMPLE}.hisat2.alignmentStatus.log" \
        | "$samtools_BIN" fixmate -@ "$THREADS" -O bam - "$MAP/${SAMPLE}.hisat2.bam" \
        2>> "$MAP/mapping_bam_processing.log"

    # ---------------------- 5. samblaster ----------------------
    "$samtools_BIN" view -@ "$THREADS" -h "$MAP/${SAMPLE}.hisat2.bam" \
        | "$samblaster_BIN" -r 2> "$MAP/samblaster.delDup.log" \
        | "$samtools_BIN" view -@ "$THREADS" -b \
              -o "$MAP/${SAMPLE}.hisat2.dedup.bam" \
        2>> "$MAP/mapping_bam_processing.log"

    # ---------------------- 6. sort ----------------------
    ulimit -n 4096
    "$sambamba_BIN" sort --tmpdir="$MAP" -t "$THREADS" \
        -o "$MAP/${SAMPLE}.hisat2.dedup.sorted.bam" \
        "$MAP/${SAMPLE}.hisat2.dedup.bam" \
        2>> "$MAP/mapping_bam_processing.log"
    ulimit -n 1024
    rm -f "$MAP/${SAMPLE}.hisat2.dedup.bam" "$MAP/${SAMPLE}.hisat2.bam"

    # ---------------------- 7.featureCounts ----------------------
    "$Rscript_BIN" "$QUANT_R" 3 "$GTF" 0 gene_id \
        "$MAP/${SAMPLE}.hisat2.dedup.sorted.bam" "$MAT/${SAMPLE}" \
        > "$MAT/quantification.R.log" 2>&1
done

# ---------------------- 8. merging per-sample counts ----------------------
echo "[$(date +%F\ %T)] Start merging per-sample counts"
mkdir -p "${WORK_DIR}/Quantification"
"$python_BIN" "$MERGE_PY" -w "$WORK_DIR" -i "$SAMPLE_INFO" -a "$GTF" -c "gene"
echo "[$(date +%F\ %T)] Pipeline finished"
