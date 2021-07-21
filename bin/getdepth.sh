#!/bin/bash

#Load samples
home=$1
indir=$2
run=$3
sample=$4
#BCFile=cell_bc_file.txt

# Create directories
mkdir -p $home/$run/$sample/read_depth

# Define input and output directories
DATA=$indir/$run/ExonReadsCorrected
OUT=$home/$run/$sample/read_depth
#BC=/data/sanli267/Array_loading_Ncells_pilot-project/180409_NS500340_0257_AHMWJ5BGX5/HistSeq_out/2000_human_PBMC_S1_001

TMP=/data/sanli267/TMP
code=/opt/Drop-seq_tools/1.12

### Update with;
### /home/sanli267/software/FastQC/fastqc -f bam -j /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java -o fastqc_out/ /data/sharedData/SingleCellProcData/180429_NS500340_0265_AHVC2JBGX5/ExonReadsCorrected/HC2_0h_S2_001.bam
### Fix the files etc.

#####################################################
echo "Calculate distribution of nr of reads per cell"
    $code/BAMTagHistogram \
        I=$DATA/$sample.bam \
        O=$OUT/output.txt.gz \
        TAG=XC \
        TMP_DIR=$TMP \
	READ_QUALITY=0

gzip -d $OUT/output.txt.gz
cat $home/$run/$sample/reads.dge.summary.txt | tail -n +4 | cut -f 1 > $OUT/cells.txt
echo NUM_READS CELL_BARCODE > $OUT/depth_per_cell.txt
grep -Fwf $OUT/cells.txt $OUT/output.txt >> $OUT/depth_per_cell.txt
rm $OUT/output.txt
rm $OUT/cells.txt

#echo $run $sample >> $home/NReadsPerCell_summary.txt
#awk '{sum+=$1;a[x++]=$1;b[$1]++}b[$1]>Mode{Mode=$1}END{print "Mean: " sum/x "\nMedian: "a[int((x-1)/2)]"\nMax: " Mode}' $OUT/depth_per_cell.txt >> $home/NReadsPerCell_summary.txt
#N=$(wc -l $OUT/depth_per_cell.txt | cut -f 1 -d ' ')
#echo N: $N >> $home/NReadsPerCell_summary.txt
