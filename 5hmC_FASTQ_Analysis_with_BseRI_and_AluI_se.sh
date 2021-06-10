#!/bin/bash                                                                                                                                                  
                                       
#PBS -l nodes=1:ppn=6                                                                                                                                        
                                     
#PBS -l walltime=2:00:00:00


#User Input    Make OUT_NAME the same as the fastq files before _L00#_R#_001.fastq
START_DIR="/home/achialastri/Plate37_H9_scAbaSI_with_Genome_and_Mitodondria/ACSD020420-P37L2_Left-scH9-BseRI_and_AluI_and_AbaSI"




#Standard Usage, no input required
BARCODES="/home/achialastri/perlscripts/5hmC/aba_barcodes.csv"
GENOME="/home/achialastri/Genomes/hg19_Zymo_LambdaPhage/hg19_Zymo_LambdaPhage.fa"
PERL_DIR="/home/achialastri/perlscripts/5hmC"
ALUIBARCODES="/home/achialastri/perlscripts/AluI/ALUI.txt"

#Do not change
PAST_DIR=${START_DIR%/*}
OUT_NAME=${START_DIR##*/}
RUN_NAME=${PAST_DIR##*/}
R1="_L001-4_R1_001.fastq"
R2="_L001-4_R2_001.fastq"
FASTQ_R1=$OUT_NAME$R1
FASTQ_R2=$OUT_NAME$R2

intermediateR1=${FASTQ_R1%??????}
intermediateR2=${FASTQ_R2%??????}
OUT_NAME_R1=$intermediateR1-ABA
OUT_NAME_R2=$intermediateR2-ABA

#Cat Fastq Files
cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R1
#Fake R2 is made from R1 since a R2 is needed for the next perl file but is not used downstream in identifying 5hmC or BseRI sites
cat $START_DIR/*L001_R1* $START_DIR/*L002_R1* $START_DIR/*L003_R1* $START_DIR/*L004_R1* > $START_DIR/$FASTQ_R2


#Order Matters for Arguments, MSPJI Barcodes are same format as ALUI
perl $PERL_DIR/ExtractingAbaReads_96BC_SimultanousWithAluI_UserInput.pl $START_DIR $OUT_NAME_R1 $OUT_NAME_R2 $FASTQ_R1 $FASTQ_R2 $BARCODES $ALUIBARCODES

#Mapping
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa aln -q 0 -n 0.04 -k 2 -l 200 -t 6 -B 6 $GENOME $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME_R1.sai
/home/cwangsanuwat/bwa/bwa-0.7.15/bwa samse -n 100 $GENOME $START_DIR/$OUT_NAME_R1.sai $START_DIR/$OUT_NAME_R1.fastq > $START_DIR/$OUT_NAME-se.sam

#call 5hmC sites
perl $PERL_DIR/process_scaba_with_BseRI_AnyInput.pl $GENOME $START_DIR/$OUT_NAME-se.sam $PERL_DIR/aba_barcodes.csv

#Convert faba files into simpler base 1 counted text files
perl $PERL_DIR/Aba_or_BseRI_FabaToText.pl $START_DIR/$OUT_NAME-se-ABA.faba
perl $PERL_DIR/Aba_or_BseRI_FabaToText.pl $START_DIR/$OUT_NAME-se-BseRI.faba

#mapping info
/home/cwangsanuwat/src/samtools/samtools flagstat $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-2bpOverhang-se.flagstat


#Same as Bam File
/home/cwangsanuwat/src/samtools/samtools view -bS $START_DIR/$OUT_NAME-se.sam > $START_DIR/$OUT_NAME-2bpOverhang-se.bam



#Remove excess files created
rm $START_DIR/$OUT_NAME-se.sam

rm $START_DIR/$OUT_NAME_R1.sai

rm $START_DIR/$OUT_NAME_R1.fastq 
rm $START_DIR/$OUT_NAME_R2.fastq 
 
rm $START_DIR/$FASTQ_R1
rm $START_DIR/$FASTQ_R2

#copy all important files to a single location for that run
FABA="_FABA"
FLAGSTAT="_FLAGSTAT"
RABA="_RABA"


mkdir $PAST_DIR/$RUN_NAME$FABA
cp $START_DIR/$OUT_NAME-se.faba $PAST_DIR/$RUN_NAME$FABA

mkdir $PAST_DIR/$RUN_NAME$FLAGSTAT
cp $START_DIR/$OUT_NAME-se.flagstat $PAST_DIR/$RUN_NAME$FLAGSTAT

mkdir $PAST_DIR/$RUN_NAME$RABA
cp $START_DIR/$OUT_NAME-se.raba $PAST_DIR/$RUN_NAME$RABA




