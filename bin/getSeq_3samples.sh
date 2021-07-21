#!/bin/bash

## The inputs parameters are (aboslute paths):
##   - Single Cell Seq folder
##   - Genome
## nohups into a folder
## Given a sample dir SampleDir, the structure of the folder will be
## /SampleDir
##  +--/Bcl2Bam
##  +--/TagReads
##  |   +--/TagCell
##  |   +--/TagUmi
##  +--/TrimReads
##  |   +--/TrimAdapter
##  |   +--/TrimPolyA
##  |   +--/MergeTrimmed
##  |   +--/Fastq
##  +--/STAR
##  |   +--/Aligment
##  |   +--/MergeLanes


###############
## FUNCTIONS ##
###############


work_flow(){
    echo ""
    echo "Here the workflow of the precoessing data
+---------+      +--------------+      
| Bcl2Bam | ---> | Taging Reads | ---> 
+---------+      | Cell & UMIs  |
                 +--------------+

       |
       V
+-----------------+
|   Trim Reads    |
| Adapter & PolyA |
+-----------------+
        |
        V
+---------------+
| STAR aligment |
+---------------+
        |
        V

| Merge Barcodes |
| Tag Exons |
       |
       V
+-------------+
| Merge Files |
+-------------+
       |
       V
| Filter Barcodes |
| |
"
    echo ""
}


clean_dir(){
    ## this function checks if a string ends with '/' and, if it the
    ## case, remove it.
    if [[ $1 =~ .*/$ ]]; then
	echo ${1%/}
    else
	echo ${1}
    fi
}


check_dir(){
    ## just check if the dir exist. If exists choose the smallest
    ## natural number 'n' that the DIR ended with '-n' dosn't exist.
    if [ -e $1 ]; then
	i=1
	s=0
 	while [ "$s" -eq "0" ]
	do
	    if [ -e "$1-$i" ]; then
		i=`expr $i \+ 1`
	    else
		echo "${1}-${i}"
		s=1
	    fi
	done
    else
	echo $1
    fi
}


mkdirs(){
    ## This functions makes all the required dirs for the output
    ## files.
    if [ -e $1 ]; then
	# echo "Xaval, la que estas liant (${1##*/})!!"
	echo "Xaval, la que estas liant (${1})!!"
	exit 1
    else
	if ! mkdir -p $1; then
	    echo "The $1 folder was not created!"
	    exit 1
	fi
    fi
    
}


bcl2bam(){
    # Here two actions are performed:
    #  1. the bcl files are converted to fastq
    #  2. the fastq files are converted to [bs]am
    # the (2) only on identified samples.
    local outdir=$output_dir/Bcl2Bam
    local outnohup=$info_dir/Bcl2Bam
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
    	mkdirs $outnohup
    fi
    echo " +--/info"
    echo " |   +--/Bcl2Bam"
    echo " |   |   +--*.bam"
    nohup bcl2fastq \
    	  -R $input_dir \
    	  -o $outdir \
    	  -p 36 \
    	  -w 36 \
    	  --minimum-trimmed-read-length 20 \
    	  --mask-short-adapter-reads 18 \
	  --stats-dir=$outnohup/Stats \
	  --reports-dir=$outnohup/Reports \
    &> $outnohup/nohup_bcl2fastq.out
    ## finally, to [sb]am files.  Here it is assumed that the files from
    ## bcl2fastq functions have the format
    ##    {sample_name}_S{number_of_sample}_L{number of lane}_R{read}_001.fastq.gz
    ## where {number_of_sample} looks like to be the position into the
    ## SampleSheet file.
    echo " |   |   +--*.fastq.gz"
    local x=`ls $outdir`
    for v in $x; do
    	if [[ (! $v =~ .*Undetermined.*) && ($v =~ .*"_R1_".*) ]]; then
	    local v2=${v/_R1_/_R2_}
	    local outfile=${v/_R[0-9]_/_}
	    outfile=${outfile/\.fastq\.gz/.bam}
	    nohup java -Xmx4g -Djava.io.tmpdir=$tmp_dir/tmpBcl2fasq -jar $PICARD/picard.jar FastqToSam \
    		  F1=$outdir/$v \
    		  F2=$outdir/$v2 \
    		  O= $outdir/$outfile\
    		  SM=fastq_raw_files_to_bam \
	    &> $outnohup/nohup_${outfile/\.bam/\.out} &
    	fi
    done
    wait    
    # for v in $x; do
    # 	if [[ $v =~ .*\.fastq.gz ]]; then
    # 	    rm $outdir/$v
    # 	fi
    # done    
}


tag_barcode(){
    # This function tags the cell or molecular barcodes (the later
    # also knowns as umi)
    # This functions need the output from 'bcl2fastq' function
    # The arguments expected are:
    #   -
    echo " |   |   +--/TagCell"
    if (($# < 0)); then
	echo ""
	echo "I need more arguments!"
	echo ""
	exit 1
    fi
    local outdir=$output_dir/TagReads/TagCell
    local outnohup=$info_dir/TagReads/TagCell
    local indir=$output_dir/Bcl2Bam
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    local x=`ls $indir`
    for v in $x; do
	if [[ $v =~ .*\.bam$ ]]; then
	    nohup $dropseq/TagBamWithReadSequenceExtended \
    	    	  INPUT=$indir/$v \
    	    	  OUTPUT=$outdir/$v \
    	    	  SUMMARY=$outnohup/summary_${v/\.bam/\.txt} \
    	    	  BASE_RANGE=1-12 \
    	    	  BASE_QUALITY=10 \
    	    	  BARCODED_READ=1 \
    	    	  DISCARD_READ=False \
    	    	  TAG_NAME=XC \
    	    	  NUM_BASES_BELOW_QUALITY=1 \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


tag_umi(){
    # This function tags the cell or molecular barcodes (the later
    # also knowns as umi)
    # This functions need the output from 'bcl2fastq' function
    # The arguments expected are:
    #   -  
    echo " |   |   +--/TagUmi"
    if (($# < 0)); then
	echo ""
	echo "I need more arguments!"
	echo ""
	exit 1
    fi
    local outdir=$output_dir/TagReads/TagUmi
    local outnohup=$info_dir/TagReads/TagUmi
    local indir=$output_dir/TagReads/TagCell
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    local x=`ls $indir`
    for v in $x; do
	if [[ $v =~ .*\.bam$ ]]; then
	    nohup $dropseq/TagBamWithReadSequenceExtended \
		  INPUT=$indir/$v \
    	    	  OUTPUT=$outdir/$v \
		  SUMMARY=$outnohup/summary_${v/\.bam/\.txt} \
		  BASE_RANGE=13-20 \
		  BASE_QUALITY=10 \
		  BARCODED_READ=1 \
		  DISCARD_READ=True \
		  TAG_NAME=XM \
		  NUM_BASES_BELOW_QUALITY=1 \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


filter_bam(){
    ## the input files are bam files into TagUmi
    echo " |   |   +--/FilterBam"
    if (($# < 0)); then
	echo ""
	echo "I need more parameters!"
	echo ""
	exit 1
    fi
    local indir=$output_dir/TagReads/TagUmi
    local outdir=$output_dir/TagReads/FilterBam
    local outnohup=$info_dir/TagReads/FilterBam
    local x=`ls $indir`
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    for v in $x; do
	if [[ $v =~ .*\.bam ]]; then
	    nohup $dropseq/FilterBAM \
		  TAG_REJECT=XQ \
		  INPUT=$indir/$v \
		  OUTPUT=$outdir/$v \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


tag_reads(){
    # This funcion make:
    #   1. Tag cell barcodes
    #   2. Tag umis
    #   3. Filter bad quality reads
    # and the input are the bam file from Bcl2Bam
    echo " |   +--/TagReads"
    local  outdir=$output_dir/TagReads
    local  outnohup=$info_dir/TagReads
    # echo "tagging reads"
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    tag_barcode
    tag_umi
    filter_bam
}


trim_adapter(){
    # The input files are from the TrimReads
    echo " |   |   +--/TrimAdapter"
    local outdir=$output_dir/TrimReads/TrimAdapter
    local outnohup=$info_dir/TrimReads/TrimAdapter
    local indir=$output_dir/TagReads/FilterBam
    local x=`ls $indir`
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    for v in $x; do
	if [[ $v =~ .*\.bam ]];then
	    $dropseq/TrimStartingSequence \
		INPUT=$indir/$v \
		OUTPUT=$outdir/$v \
		OUTPUT_SUMMARY=$outnohup/summary_${v/\.bam/\.txt} \
		SEQUENCE=AAGCAGTGGTATCAACGCAGAGTGAATGGG \
		MISMATCHES=0 \
		NUM_BASES=5 \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


trim_polyA(){
    # The input is from trim_adapter
    echo " |   |   +--/TrimPolyA"
    local outdir=$output_dir/TrimReads/TrimPolyA
    local outnohup=$info_dir/TrimReads/TrimPolyA
    local indir=$output_dir/TrimReads/TrimAdapter
    x=`ls $indir`
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    for v in $x; do
	if [[ $v =~ .*\.bam ]]; then
	    $dropseq/PolyATrimmer \
		INPUT=$indir/$v \
		OUTPUT=$outdir/$v \
		OUTPUT_SUMMARY=$outnohup/summary_${v/\.bam/\.txt} \
		MISMATCHES=0 \
		NUM_BASES=6 \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


trim_reads(){
    # Trims the seqs:
    #  1. trim the adapter
    #  2. trims the poly A tail
    #  3. Merge files
    echo " |   +--/TrimReads"
    local outdir=$output_dir/TrimReads
    local outnohup=$info_dir/TrimReads
    local outsubdir=$outdir/Fastq
    local outsubnohup=$outnohup/Fastq
    local indir=$outdir/TrimPolyA
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    trim_adapter
    trim_polyA
    echo " |   |   +--Fastq"
    mkdirs $outsubdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outsubnohup
    fi
    local outdir=$outsubdir
    local outnohup=$outsubnohup
    local x=`ls $indir`
    for v in $x; do
    	if [[ $v =~ .*\.bam ]]; then
    	    nohup java -Xmx4g -jar $PICARD/picard.jar SamToFastq\
    	    	  INPUT=${indir}/$v \
    	    	  FASTQ=${outdir}/${v/\.bam/\.fastq} \
    	    &> ${outnohup}/nohup_${v/\.bam/\.out} &
    	fi
    done
    wait
}



sort_files0(){
    nohup java -Xmx4g -jar $PICARD/picard.jar SortSam \
    	  I=$1 \
    	  O=$2 \
    	  SO=queryname \
          TMP_DIR=/data/sharedData/tmp \
    &> $3 &
}

sort_files_from_star(){
    echo " |   |   +--/SortAligned"
    local indir1=$output_dir/STAR/Aligment
    local outdir1=$output_dir/STAR/SortAligned
    local outnohup=$info_dir/STAR/SortAligned
    mkdirs $outdir1
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    local x=`ls $indir1`
    y=()
    local i=0
    for v in $x; do
	if [[ $v =~ .*\.sam$ ]]; then
	    ((i=i+1))
	    local infile=${v%%\.*}
	    sort_files0 $indir1/$v $outdir1/${infile}.bam $outnohup/nohupt_${infile}.out
	    if [[ $(($i%$picardruns)) == 0 ]]; then
		wait
	    fi
	fi
    done
    wait
}


star_aligment(){
    # The input are the fastq files from trim_reads
    local outdir=$output_dir/STAR
    local outsubdir=$output_dir/STAR/Aligment
    local outnohup=$info_dir/STAR
    local outsubnohup=$info_dir/STAR/Aligment
    local indir=$output_dir/TrimReads/Fastq
    mkdirs $outdir
    mkdirs $outsubdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    if [[ rm_tmp=true ]]; then
	mkdirs $outsubnohup
    fi
    local outdir=$outsubdir
    local outnohup=$outsubnohup
    local x=`ls $indir`
    local tmpdir=`check_dir $tmp_dir/tmpSTAR`
    mkdirs $tmpdir
    local i=0
    for v in $x; do
    	if [[ $v =~ .*\.fastq ]]; then
    	    ((i=i+1))
    	    nohup $star/STAR --runThreadN $starprocs \
    	    	  --genomeDir $genome_dir \
    	    	  --readFilesIn $indir/$v \
    	    	  --outFileNamePrefix  $outdir/${v}_ \
    	    	  --outTmpDir $tmpdir/${v%%\.*} \
    	    &> $outnohup/nohup_${v/\.fastq/\.out} &
	    if [[ $(($i%$starruns)) == 0 ]]; then
		wait
	    fi
	fi
    done
    wait
    # move Log outputs to main folder
    local x=`ls $outdir`
    for v in $x; do
	if [[ ($v =~ .*Log.*) || ($v =~ .*\.tab$) ]]; then
	    mv $outdir/$v $outnohup/$v
	fi
    done
    echo " |   +--/STAR"
    echo " |   |   +--/Aligment"
    sort_files_from_star
    # for v in ${tmpDirs[*]}; do echo $v; done
}


merge_barcodes(){
    ## reference file eg /data/sharedData/Genomes/GRCm38.p5/Fasta/GRCm38.primary_assembly.genome.fa
    nohup java -Xmx4g -jar $PICARD/picard.jar MergeBamAlignment \
	  REFERENCE_SEQUENCE=$4 \
	  UNMAPPED_BAM=$2 \
	  ALIGNED_BAM=$1 \
	  OUTPUT=$3 \
	  INCLUDE_SECONDARY_ALIGNMENTS=false \
	  PAIRED_RUN=false \
          TMP_DIR=/data/sharedData/tmp \
    &> $5 &
}


merge_aligned_with_barcodes(){
    # Here, the aligned reads are merged with (cell & molecular)
    # barcodes.
    # The input files are from the aligment step (sorted).
    echo " |   +--/MergeAlignedTags"
    if (($# < 1)); then
	echo ""
	echo "I need the fasta reference !!!"
	echo ""
	exit 1
    fi
    local indir1=$output_dir/STAR/SortAligned
    local indir2=$output_dir/TrimReads/TrimPolyA
    local outdir=$output_dir/MergeWithBarcodes
    local outnohup=$info_dir/MergeWithBarcodes
    mkdirs $outdir
    if [[ rm_tmp=true ]]; then
	mkdirs $outnohup
    fi
    local x1=`ls $indir1`
    local x2=`ls $indir2`
    for v in $x1; do
	if [[ $v =~ .*\.bam$ ]]; then
	    if [[ -e $indir2/$v ]]; then
		local nohup_file=$outnohup/nohup_${v/\.bam/\.out}
		merge_barcodes $indir1/$v $indir2/$v $outdir/$v $1 $nohup_file
	    else
		echo ""
		echo "Merging With Barcodes. The file "
		echo "         '$indir2/$v'"
		echo " doesn't exists."
		echo ""
	    fi
	fi
    done
    wait
}


tag_exons(){
    # Here the input are the ones from merge_aligned_with_barcodes
    echo " |   +--/TagExons"
    local indir=$output_dir/MergeWithBarcodes
    local outdir=$output_dir/TagExons
    local outnohup=$info_dir/TagExons
    mkdirs $outdir
    if [[ rm_tmp=true && !($outdir = $outnohup)]]; then
    	mkdirs $outnohup
    fi
    echo " |   |   +--/TagExons"
    local x=`ls $indir`
    for v in $x; do
	if [[ $v =~ .*.bam$ ]]; then
	    nohup $dropseq/TagReadWithGeneExon \
    		  I=$indir/$v \
    		  O=$outdir/$v \
    		  ANNOTATIONS_FILE=$annot_file \
		  TAG=GE \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


merge_and_sort(){
    local outdir1=${5%/*} # tmp folder
    local outdir1=${outdir1##*/}
    local file1=${5##*/}
    local outdir2=${6%/*} # home dir
    local outdir2=${outdir2##*/}
    local file2=${6##*/}
    nohup java -Xmx4g -jar $PICARD/picard.jar MergeSamFiles \
    	  I=$1 \
    	  I=$2 \
    	  I=$3 \
    	  I=$4 \
    	  O=$5 \
    &> $info_dir/$outdir1/nohup_${file1/\.bam/\.out}
    echo " |   +--/ExonsReadsNotCorrected"
    nohup java -Xmx4g -jar $PICARD/picard.jar SortSam \
	  I=$5 \
	  O=$6 \
	  SO=coordinate \
          TMP_DIR=/data/sharedData/tmp \
    &> $info_dir/$outdir2/nohup_${file2/\.bam/\.out}
}

merge_aligned_files(){
    local indir=$output_dir/TagExons
    local outdir1=$output_dir/MergeAlignedMappedReads
    local outnohup1=$info_dir/MergeAlignedMappedReads
    local outdir2=$home_dir/ExonReadsNotCorrected
    local outnohup2=$info_dir/ExonReadsNotCorrected
    echo " |   +--/MergeAlignedMappedReads"
    mkdirs $outdir1
    mkdirs $outdir2
    if [[ rm_tmp=true ]]; then
    	mkdirs $outnohup1
    	mkdirs $outnohup2
    fi    
    local x=`ls $indir`
    local files1=()
    local files2=()
    local files3=()
    for v in $x; do
	if [[ $v =~ .*\.bam$ ]]; then
	    if [[ $v =~ .*_S1_L00.*\.bam$ ]]; then
		files1+=($v)
	    elif [[ $v =~ .*_S2_L00.*\.bam$ ]]; then
		files2+=($v)
	    elif [[ $v =~ .*_S3_L00.*\.bam$ ]]; then
		files3+=($v)
	    fi	    
	fi
    done
    if [[ ${#files1[*]} -ne 4 ]]; then
	echo ""
	echo "It is supesed to have 4 files for each sample (got ${#files1[*]})"
	echo ""
	exit 1
    fi
    if [[ ${#files2[*]} -ne 4 ]]; then
	echo ""
	echo "It is supesed to have 4 files for each sample (got ${#files2[*]})"
	echo ""
	exit 1
    fi
    if [[ ${#files3[*]} -ne 4 ]]; then
	echo ""
	echo "It is supesed to have 4 files for each sample (got ${#files3[*]})"
	echo ""
	exit 1
    fi    
    local outfile=${files1[0]/_L00[0-9]_/_}
    merge_and_sort \
	$indir/${files1[0]} \
	$indir/${files1[1]} \
	$indir/${files1[2]} \
	$indir/${files1[3]} \
	$outdir1/$outfile \
	$outdir2/$outfile &
    
    local outfile=${files2[0]/_L00[0-9]_/_}
    merge_and_sort \
	$indir/${files2[0]} \
	$indir/${files2[1]} \
	$indir/${files2[2]} \
	$indir/${files2[3]} \
	$outdir1/$outfile \
	$outdir2/$outfile &

    local outfile=${files3[0]/_L00[0-9]_/_}
    merge_and_sort \
	$indir/${files3[0]} \
	$indir/${files3[1]} \
	$indir/${files3[2]} \
	$indir/${files3[3]} \
	$outdir1/$outfile \
	$outdir2/$outfile &    
    wait
}


bead_synthesis_errors(){
    local indir=$home_dir/ExonReadsNotCorrected
    local outdir=$home_dir/ExonReadsCorrected
    local outnohup=$info_dir/ExonReadsCorrected
    echo " |   +--/ExonReadsCorrected"
    mkdirs $outdir
    if [[ rm_tmp=true && !($outdir = $outnohup) ]]; then
    	mkdirs $outnohup
    fi
    local x=`ls $indir`
    for v in $x; do
	if [[ $v =~ .*\.bam$ ]]; then
	    $dropseq/DetectBeadSynthesisErrors \
    		I=$indir/$v \
    		O=$outdir/$v \
    		OUTPUT_STATS=$outdir/my.synthesis_stats_${v/\.bam/\.txt} \
    		SUMMARY=$outdir/summary_file_${v/\.bam/\.txt} \
    		NUM_BARCODES=$num_cells \
    		PRIMER_SEQUENCE=AAGCAGTGGTATCAACGCAGAGTAC \
	    &> $outnohup/nohup_${v/\.bam/\.out} &
	fi
    done
    wait
}


##########
## MAIN ##
##########

## set rm_tmp to false if you want to perserve all the data generated
## by executing the pipeline.
rm_tmp=true
## starruns define the number of times to run STAR
## to compute the number, the size of the gnome sould be computed
## eg, if the size is 25GB, then running STAR 8 times will use 200GB
starruns=8
picardruns=9
## now number of processes for each STAR instance
starprocs=4


if [[ $# -lt 4 ]]; then
    echo ""
    echo "not all argumets are given!"
    echo ""
    exit 1
fi

if [[ $# -gt 5 ]]; then
    echo ""
    echo "Too many argumets are given!"
    echo ""
    exit 1
fi

clear


# initial variables
PICARD=/opt/picard/current
dropseq=/opt/Drop-seq_tools/1.12
star=/opt/star/2.5.3a/bin/Linux_x86_64


if [[ ! $PATH =~ .*/opt/bcl2fastq2-v2.19.1/bin.* ]]; then
    PATH="$PATH:/opt/bcl2fastq2-v2.19.1/bin"
fi

user=$(id -un)
host=$(hostname)

if [ ! -e /data/$user ]; then
    echo ""
    echo "Looks like $host doesn't know you."
    echo "Please, seek some advice from your counselor!"
    echo " I am really sad for you :("
    echo ""
    exit 1
fi


# taking input_dir in right format
# eg /data/sharedData/SingleCellRawData/170725_NS500340_0199_AHCG5FBGX2
input_dir=`clean_dir $1`
# eg /data/sharedData/Genomes/GRCm38.p5/GenomeDir
genome_dir=`clean_dir $2`
# eg /data/sharedData/Genomes/GRCm38.p5/Fasta/GRCm38.primary_assembly.genome.fa
fasta_ref=$3
# eg /data/sharedData/Genomes/GRCm38.p5/GTF/gencode.vM14.primary_assembly.annotation.gtf
annot_file=$4
# Number of expected cells into the array
# This argument if optional
if [[ $# -eq 5 ]]; then 
    num_cells0=$5
    num_cells=$((2*$num_cells0))
else
    num_cells0=false
    num_cells=false
fi


sample=${input_dir##*/}

# home_dir is tha main output folder
prefix=/data/${user}/getSeq/
home_dir=${prefix}$sample
home_dir=${home_dir}
home_dir=`check_dir $home_dir`
mkdirs $home_dir

info_dir=$home_dir/info
mkdirs $info_dir

tmp_dir=$home_dir/tmp
mkdirs $tmp_dir

# defining DIRs
if [ rm_tmp=true ]; then
    output_dir=$tmp_dir
else
    output_dir=$home_dir
fi


echo -e "$0 $@\n" > $home_dir/call
echo -e "$0 $@"
echo ""
echo "Here who you are and your options:"
echo "user:       "$user
echo "host:       "$host
echo "sample:     "$sample
echo "input:      "$input_dir
echo "output:     "$home_dir
echo "genome:     "$genome_dir
echo "fasta:      "$fasta_ref
echo "annotation: "$annot_file
echo "tmp:        "$tmp_dir
echo "remove tmp: "$rm_tmp
echo ""


echo "And here a summary of the data:"
echo "sample Ids:   "
echo "sample names: "
echo "num cells:    "$num_cells0
echo ""

echo "The following structure will be created and uncovered through
the pipeline's execution at each respective step, in which the
information files will be stored."

echo $home_dir
echo " +--/info
 |   +--/Bcl2Bam
 |   |   +--*.bam
 |   |   +--*.fastq.gz
 |   +--/TagReads
 |   |   +--/TagCell
 |   |   +--/TagUmi
 |   |   +--/FilterBam
 |   +--/TrimReads
 |   |   +--/TrimAdapter
 |   |   +--/TrimPolyA
 |   |   +--Fastq
 |   +--/STAR
 |   |   +--/Aligment
 |   |   +--/SortAligned
 |   +--/MergeAlignedTags
 |   +--/TagExons
 |   |   +--/TagExons
 +--/AlignedReads
 +--CorrectedBeads"
echo ""

echo "A 'tmp' folder will be crated and the same folders and structure
as explained above will be created. This purpose of that folder is to
store temporary files/folders or data used to get teh final aligned data."

echo $home_dir
echo " +--/tmp"
echo ""

echo "If you want to save the output data from each step, then you
must enable the 'save tmp files' (change 'rm_tmp' to false). Then, the
data will be stored into the folders into the main folder."
echo ""

echo "Oh yes, just in case, the output contains reads that are not sorted!"
echo ""

echo "So, let us start, be patient..."
echo ""
echo $home_dir


## 1. Convert the files from the sequencer
##############################################

bcl2bam

# 2. Tagging the sequences
##########################

tag_reads


# 3. Trimming sequences
#######################

trim_reads


# 4. Aligment
#############
## here:
##   i. the reads in each lane and sample are mapped.
##  ii. the aligned files are sorted
## iii. the unaligned files are sorted (from TagUmi)
## Here, we sort the unaligned because we can run  in parallel with aligned.

star_aligment

# 5. Merge with Barcodes
########################
# Now the files are sorted and we merge with (cell and molecular) barcodes

merge_aligned_with_barcodes $fasta_ref


# 6. Tag with exons
###################

tag_exons


# 7. merge files
################

merge_aligned_files

# 8. Detecting errors in barcodes
#################################
# This step is optional and is performed if the optional parameter 'number of cells' is given

if [[ num_cell != false ]]; then
    bead_synthesis_errors
fi

rm $tmp_dir -r

echo ""
echo "Process Done :)"
echo ""

exit 0

# 8. Moving the final outputs
#############################
# They are from bead_synthesis_errors
final_files(){
    local indir=$output_dir/FilterBeads
    local outdir=$home_dir
    local x=`ls $indir`
    for v in $x; do
	if [[ $v =~ .*\.bam$ ]]; then
	    mv $indir/$v $outdir/$v
	    echo " +--"$v
	fi
    done
}

final_files

echo ""
echo "Process Done :)"

exit 0

