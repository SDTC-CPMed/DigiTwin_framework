#!/bin/bash

indir=data
infile=_min200genesPerCell_sorted_expression_matrix.csv
outdir=data/knn_smoothing
outfile=_min200genesPerCell_sorted_expression_matrix

### Run
tissue=HA

for i in 14
do
  python3 knn_smooth.py --k $i --fpath $indir/$tissue$infile --saveto $outdir/$tissue$outfile.knn-smooth_k$i.csv --sep ,
done

tissue=HC

for i in 21
do
  python3 knn_smooth.py --k $i --fpath $indir/$tissue$infile --saveto $outdir\
/$tissue$outfile.knn-smooth_k$i.csv --sep ,
done
gzip $outdir/*.csv


