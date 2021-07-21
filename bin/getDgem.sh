#!/bin/bash

#Load samples
home=$1
indir=$2
run=$3
sample=$4
#BCFile=cell_bc_file.txt

# Create directories
mkdir -p $home/$run/$sample
mkdir -p $home/$run/fastqc_out

# Define input and output directories
DATA=$indir/$run/ExonReadsCorrected
OUT=$home/$run/$sample
#BC=/data/sanli267/Array_loading_Ncells_pilot-project/180409_NS500340_0257_AHMWJ5BGX5/HistSeq_out/2000_human_PBMC_S1_001

TMP=/data/sharedData/tmp
code=/opt/Drop-seq_tools/1.12

### Update with;
### /home/sanli267/software/FastQC/fastqc -f bam -j /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -o fastqc_out/ /data/sharedData/SingleCellProcData/180429_NS500340_0265_AHVC2JBGX5/ExonReadsCorrected/HC2_0h_S2_001.bam
### Fix the files etc.

#####################################################
echo "Remove short reads"

/opt/samtools/1.4/bin/samtools view -h $DATA/$sample.bam | awk 'length($10) > 35 || $1 ~ /^@/' | /opt/samtools/1.4/bin/samtools view -bS - > $OUT/lengthFiltered_$sample.bam

echo "fastqc"

/home/sanli71/software/FastQC/fastqc -f bam -j /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -o $home/$run/fastqc_out/ $OUT/lengthFiltered_$sample.bam

echo "Primary output - UMIs per gene and cell"

    $code/DigitalExpression \
	I=$OUT/lengthFiltered_$sample.bam \
	O=$OUT/umi.dge.txt.gz \
	SUMMARY=$OUT/umi.dge.summary.txt \
	CELL_BARCODE_TAG=XC \
	MIN_NUM_READS_PER_CELL=10000 \
	TMP_DIR=$TMP

echo "Extra output - reads per gene and cell"

    $code/DigitalExpression \
        I=$OUT/lengthFiltered_$sample.bam \
        O=$OUT/reads.dge.txt.gz \
        SUMMARY=$OUT/reads.dge.summary.txt \
        CELL_BARCODE_TAG=XC \
        MIN_NUM_READS_PER_CELL=10000 \
	OUTPUT_READS_INSTEAD=true \
	TMP_DIR=$TMP

