# Extract data
### Human data
/opt/getSeq/bin/getSeq_1sample /data/sharedData/SingleCellRawData/180810_NS500340_0296_AHK2W2BGX7 /data/sharedData/Genomes/GRCh38.p10/GenomeDir/ /data/sharedData/Genomes/GRCh38.p10/Fasta/GRCh38.primary_assembly.genome.fa /data/sharedData/Genomes/GRCh38.p10/GTF/gencode.v27.primary_assembly.annotation.gtf 20000

### get Digital Expression Matrix
/opt/getSeq/bin/getDgem /data/sanli71/RA-CIA_project/CIA_VT2019_project/scData /data/sanli71/getSeq 180810_NS500340_0296_AHK2W2BGX7 HC7_12h_D_S1_001
/opt/getSeq/bin/getDgem /data/sanli71/RA-CIA_project/CIA_VT2019_project/scData /data/sanli71/getSeq 180810_NS500340_0296_AHK2W2BGX7 HC7_2day_A_S2_001
/opt/getSeq/bin/getDgem /data/sanli71/RA-CIA_project/CIA_VT2019_project/scData /data/sanli71/getSeq 180810_NS500340_0296_AHK2W2BGX7 HC8_2day_A_S3_001

### Check read depth
### $1:home $2:indir $3:run $4:sample
/opt/getSeq/bin/getdepth /data/sanli71/RA-CIA_project/CIA_VT2019_project/scData /data/sanli71/getSeq 180810_NS500340_0296_AHK2W2BGX7 HC7_12h_D_S1_001
/opt/getSeq/bin/getdepth /data/sanli71/RA-CIA_project/CIA_VT2019_project/scData /data/sanli71/getSeq 180810_NS500340_0296_AHK2W2BGX7 HC7_2day_A_S2_001
/opt/getSeq/bin/getdepth /data/sanli71/RA-CIA_project/CIA_VT2019_project/scData /data/sanli71/getSeq 180810_NS500340_0296_AHK2W2BGX7 HC8_2day_A_S3_001
