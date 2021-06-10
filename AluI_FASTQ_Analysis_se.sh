#!/bin/bash                                                   
#PBS -l nodes=1:ppn=6                                                                                                                                    
#PBS -l walltime=3:00:00:00

#User Input    
START_DIR="/home/achialastri/Plate37_H9_scAbaSI_with_Genome_and_Mitodondria/ACSD020420-P37L2_Left-scH9-BseRI_and_AluI_and_AbaSI"

#Analyze AluI reads with BseRI/AbaSI Reads - Taken From 12/10/19 David Podorefsky and edited by Alex Chialastri2/12/2020

#Standard Usage, no input required
BARCODES="/home/achialastri/perlscripts/AluI/ALUI"
GENOME="/home/achialastri/Genomes/hg19_Zymo_LambdaPhage/hg19_Zymo_LambdaPhage.fa"
PERL_DIR="/home/achialastri/perlscripts/AluI"
ABASI_BC="/home/achialastri/perlscripts/5hmC/aba_barcodes.csv"

#Do not change
PAST_DIR=${START_DIR%/*}
#(Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq)
OUT_NAME=${START_DIR##*/}-AluI
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2
intermediateR1=${FASTQ_R1%??????}
intermediateR2=${FASTQ_R2%??????}
OUT_NAME_R1=$intermediateR1-ALUI
OUT_NAME_R2=$intermediateR2-ALUI

#Cat Fastq Files, starts as 4 lanes off of Illumina sequencing at UCSB
cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R1

# 1. #Pull out only Celseq Lines Can use the same script pretending AluI reads are MspJI because they have the same 3 bp UMI then 8 Barcode, but then AluI has a CA (set --scTHseq == 1 if done with ABASI)
perl /home/achialastri/perlscripts/MspJI/ExtractingAluIReads_UserInput_NoBarcodeCollisions.pl --FASTQ_R1 $START_DIR/$FASTQ_R1 --MSPJI_BC $BARCODES --scTHseq 1 --ABASI_BC $ABASI_BC
#rename output
mv $START_DIR/$intermediateR1-AluI.fastq $START_DIR/$OUT_NAME_R1.fastq

# 2. Map data to genome
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa aln -q 0 -n 0.04 -k 2 -l 200 -t 8 -B 13 $GENOME $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME_R1.sai

# 3. Convert .sai to .sam file
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa samse -n 100 $GENOME $START_DIR/$OUT_NAME_R1.sai $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME-se.sam
/home/cwangsanuwat/src/samtools/samtools flagstat $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se.flagstat
#3.1 Save as Bam File
/home/cwangsanuwat/src/samtools/samtools view -bS $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-se.bam


# 4. Identify AluI sites with UMI 
IN="-se.sam"
OUT="-se_MappedReads.txt"
OUT2="-se_correctpos_stringent.txt"
UMIIN=$OUT_NAME$IN
UMIOUT=$OUT_NAME$OUT
UMIOUT2=$OUT_NAME$OUT2
perl $PERL_DIR/Identifying_AluI_Sites_withUMI_SCdata_ForkedAdapters.pl $START_DIR/$UMIIN $START_DIR/$UMIOUT $START_DIR/$UMIOUT2

# 5. Simplify data
OUT3="-se_correctpos_stringent_simplified.txt"
OUTOUT3=$OUT_NAME$OUT3
perl $PERL_DIR/SimplifyData.pl $START_DIR/$UMIOUT2 $START_DIR/$OUTOUT3 $BARCODES

# 6. Sort
OUT4="-se_correctpos_stringent_simplified_sort.txt"
OUTOUT4=$OUT_NAME$OUT4
sort -k1,1n -k2,2V -k3,3n $START_DIR/$OUTOUT3 > $START_DIR/$OUTOUT4

# 7. Remove PCR duplicates
OUT5="-se_correctpos_stringent_simplified_sort_rmdup.txt"
OUTOUT5=$OUT_NAME$OUT5
perl $PERL_DIR/RemoveDup.pl $START_DIR/$OUTOUT4 $START_DIR/$OUTOUT5

#Remove excess files created
rm $START_DIR/$OUT_NAME-se.sam
rm $START_DIR/$OUT_NAME_R1.sai
rm $START_DIR/$OUT_NAME_R1.fastq 
rm $START_DIR/$OUT_NAME_R2.fastq 
rm $START_DIR/$FASTQ_R1
rm $START_DIR/$FASTQ_R2