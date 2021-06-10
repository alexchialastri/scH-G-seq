#!/usr/bin/perl

use strict "subs";
use strict "refs";
#use strict "vars";
 
print STDERR "
  processes sam file of single-cell aba-seq data
 
  Usage
  process_scaba.pl genome_file.fa sam_file.sam aba_barcodes.csv
 
";

# read in genome file
my $genomefile = $ARGV[0];
open(MULTIFASTA, "$genomefile") || die " FATAL ERROR:\n Unable to load '$genomefile'.\n";
while (<MULTIFASTA>) {
	chomp($_);
	if ($_=~/^>(.*)/) {
		if ($seq) {
			$sequence{$header}=$seq;
		}
		$header = $1;
		$counter++;
		$seq    = '';
	} else {
		$seq.=$_;
	}
 
}
$sequence{$header}=$seq;

# read in cell-specific barcode file
my $file2 = $ARGV[2] or die "Need to get cell specific barcode file on the commend line\n";
my @csbc;
open(my $cscodes, '<', $file2) or die "Could not open '$file2' $!\n";
while (my $line = <$cscodes>) {
  chomp $line;
    my @fields = split ";" , $line;
    push @csbc,@fields[0];
}

# read in sam file, filter reads, and count CG dinucleotides
my $file1 = $ARGV[1] or die "Need to get SAM file on the command line\n";
open(my $data, '<', $file1) or die "Could not open '$file1' $!\n";
for ($i=0;$i<105;$i++)
 {
   $N_CG[$i] = 0;
 }
for ($k=0;$k<96;$k++)
 {
   $allreads[$k]=0;
   $rawreads[$k] = 0;
   $mappedreads[$k] = 0;
   $cleanreads[$k] = 0;
 }
my $rABAoutputfile = substr($file1,0,-4)."-ABA.raba";
my $fABAoutputfile = substr($file1,0,-4)."-ABA.faba";
my $Mitooutputfile = substr($file1,0,-4)."-Mito.raba";
my $Zymooutputfile = substr($file1,0,-4)."-Zymo.raba";
my $BseRIoutputfile = substr($file1,0,-4)."-BseRI.raba";
my $fabaBseRIoutputfile = substr($file1,0,-4)."-BseRI.faba";

my @chrvec=qw(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY chrM ZymoDNAStandardSet PhageLambda);

open(my $abaout, '>', $rABAoutputfile) or die "Can't open file for writing: $!\n";
open(my $Mitoout, '>', $Mitooutputfile) or die "Can't open file for writing: $!\n";
open(my $BseRIout, '>', $BseRIoutputfile) or die "Can't open file for writing: $!\n";
open(my $Zymoout, '>', $Zymooutputfile) or die "Can't open file for writing: $!\n";

while (my $line = <$data>) {
    chomp $line;
    my @fields = split "\t" , $line;
    my $seqflag = $fields[1];
    my $chr = $fields[2];
    my $cutcoord = $fields[3];
    my $cigar = $fields[5];
    my $cellbarcode = substr($fields[11],5,6);
    my $uniqueflag = $fields[12];
    my $strand = 0;
    my $offset = 0;
    my $CGcount = 0;
 
if ($cellbarcode ~~ @csbc ){
    my ($csbcindex) = grep {$csbc[$_] ~~ $cellbarcode } 0 .. 95;
	$allreads[$csbcindex]++;}

if ($cellbarcode ~~ @csbc && $chr ~~ @chrvec){
    my ($csbcindex) = grep {$csbc[$_] ~~ $cellbarcode } 0 .. 95;
    $rawreads[$csbcindex]++;}


    #Pull out reads that map to the Mitocondrial DNA, disregaurd the Cigar string since we are looking for SNPS
    if (($seqflag == 0 || $seqflag == 16) && $uniqueflag =~ "XT:A:U" && $cellbarcode ~~ @csbc && $chr eq "chrM"){
       my ($csbcindex) = grep {$csbc[$_] ~~ $cellbarcode } 0 .. 95;
       print $Mitoout "$csbcindex \t $chr \t $cutcoord \t $seqflag \t $cigar \t $fields[9]\n";
    }
    elsif(($seqflag == 0 || $seqflag == 16) && $uniqueflag =~ "XT:A:U" && $cellbarcode ~~ @csbc && $chr eq "ZymoDNAStandardSet"){
	 print $Zymoout "$csbcindex \t $chr \t $cutcoord \t $seqflag \t $cigar \t $fields[9]\n";
    }
    elsif (($seqflag == 0 || $seqflag == 16) && ($cigar =~ /(^\d\d)M/) && (length($cigar)== 3) && $uniqueflag =~ "XT:A:U" && $cellbarcode ~~ @csbc && $chr ~~ @chrvec){
       my ($csbcindex) = grep {$csbc[$_] ~~ $cellbarcode } 0 .. 95;
       my ($chrindex) = grep {$chrvec[$_] ~~ $chr} 0 .. 26;       
       $mappedreads[$csbcindex]++;
       my $cutsite = uc(substr($sequence{$chr},$cutcoord-20,125));
       my $char = 'CG';
       my $result = index($cutsite, $char, $offset);

       my $MatchLength = substr($cigar,0,2);
       my $AdjustedLength = 70 - int($MatchLength);


    #This section will Separate out Reads from BseRI
       my $BseRIflag = 0;
       #Direct, BseRI Sequence originally On + strand - MSPJI stuff plus 6
       if ($seqflag == 16 & substr($cutsite,73-$AdjustedLength,6) eq "GAGGAG"){
	    $CGpos=$cutcoord+53-$AdjustedLength; $CGcand = uc(substr($sequence{$chr},$CGpos,6));$strand=1; $BseRIflag++;	    
        }
       #Direct, BseRI originally on - strand +6
       if ($seqflag == 0 & substr($cutsite,29,6) eq "CTCCTC"){	   
	   $CGpos=$cutcoord+9; $CGcand = uc(substr($sequence{$chr},$CGpos,6)); $seqmatch = reverse $CGcand; $seqmatch =~ tr/ACGTacgt/TGCAtgca/; $CGcand = $seqmatch; $strand=-1; $BseRIflag++;
       }
       #Indirect, BseRI originally on + Strand 
        if ($seqflag == 0 & substr($cutsite,5,6) eq "GAGGAG"){ 
	    $CGpos=$cutcoord-15; $CGcand = uc(substr($sequence{$chr},$CGpos,6));$strand=1; $BseRIflag++;	    
	}
       #Indirect, 5mC orignally on - strand & substr($cutsite,96-$AdjustedLength,1) eq 'G'
       if ($seqflag == 16 &  substr($cutsite,107,6) eq "CTCCTC"){
	   $CGpos=$cutcoord+77-$AdjustedLength; $CGcand = uc(substr($sequence{$chr},$CGpos,6));$seqmatch = reverse $CGcand; $seqmatch =~ tr/ACGTacgt/TGCAtgca/; $CGcand = $seqmatch; $strand=-1; $BseRIflag++;
       }





       # Look for AbaSI based Reads
       if ($result > -1) {$N_CG[$result]++; 
			  if ($result == 8  && $seqflag == 0){$CGpos = $cutcoord-12; $strand=1; $CGcount++}; #Indirect, 5mC originally on + Strand
			  if ($result == 9  && $seqflag == 0){$CGpos = $cutcoord-11; $strand=1; $CGcount++}; #Indirect, 5mC originally on + Strand
								
			  if ($result == 29  && $seqflag == 0){$CGpos = $cutcoord+9; $strand=-1; $CGcount++}; #Direct, 5mC originally on - strand
			  if ($result == 30  && $seqflag == 0){$CGpos = $cutcoord+10; $strand=-1; $CGcount++}; #Direct, 5mC originally on - strand
								
			  if ($result == 76-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+56-$AdjustedLength; $strand=1; $CGcount++}; #Direct, 5mC originally On + strand
			  if ($result == 77-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+57-$AdjustedLength; $strand=1; $CGcount++}; #Direct, 5mC originally On + strand
								
			  if ($result == 97-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+77-$AdjustedLength; $strand=-1; $CGcount++}; #Indirect, 5mC orignally on - strand
			  if ($result == 98-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+78-$AdjustedLength; $strand=-1; $CGcount++}; #Indirect, 5mC orignally on - strand
			  if ($CGcount > 0) {$CGcand = uc(substr($sequence{$chr},$CGpos,2));}
                                }

       while ($result != -1) { #Just slides down to account for hitting a CG early by chance

          $offset = $result + 2;
          $result = index($cutsite, $char, $offset);
          if ($result > -1) {$N_CG[$result]++;

			     if ($result == 8  && $seqflag == 0){$CGpos = $cutcoord-12; $strand=1; $CGcount++};
			     if ($result == 9  && $seqflag == 0){$CGpos = $cutcoord-11; $strand=1; $CGcount++};
								
			     if ($result == 29  && $seqflag == 0){$CGpos = $cutcoord+9; $strand=-1; $CGcount++};
			     if ($result == 30  && $seqflag == 0){$CGpos = $cutcoord+10; $strand=-1; $CGcount++};
								
			     if ($result == 76-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+56-$AdjustedLength; $strand=1; $CGcount++};
			     if ($result == 77-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+57-$AdjustedLength; $strand=1; $CGcount++};
								
			     if ($result == 97-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+77-$AdjustedLength; $strand=-1; $CGcount++};
			     if ($result == 98-$AdjustedLength  && $seqflag == 16){$CGpos = $cutcoord+78-$AdjustedLength; $strand=-1; $CGcount++};
			     if ($CGcount > 0) {$CGcand = uc(substr($sequence{$chr},$CGpos,2));}
        			}}
       #Print to AbaSI raba if only fits with AbaSI Data (Overlap could happen if detect the non direct strands!)
	if ($CGcount==1 & $BseRIflag == 0){	
		print $abaout "$csbcindex \t $chrindex \t $CGpos \t $strand \t $seqflag \t $CGcount \t $CGcand \n";
		$cleanreads[$csbcindex]++;}
       	if ($CGcount==0 & $BseRIflag == 1){	
		print $BseRIout "$csbcindex \t $chrindex \t $CGpos \t $strand \t $seqflag \t $BseRIflag \t $CGcand \n";
	}

    }
} 
# write CG count to file
my $CGoutputfile = substr($file1,0,-3)."CG";
open(my $out, '>', $CGoutputfile) or die "Can't open file for writing: $!\n";
for ($i=0;$i<105;$i++)
{
    print $out "$N_CG[$i]\n";
}
# write QC count to file
my $QCoutputfile = substr($file1,0,-3)."QC";
open(my $qcout, '>', $QCoutputfile) or die "Can't open file for writing: $!\n";
my $QCoutputfile2 = substr($file1,0,-3)."QC2";
open(my $qcout2, '>', $QCoutputfile2) or die "Can't open file for writing: $!\n";

for ($k=0;$k<96;$k++)
 {
   print $qcout "$k \t $rawreads[$k] \t $mappedreads[$k] \t $cleanreads[$k] \n";
   print $qcout2 "$allreads[$k]\n";
 }

# keep only unique cuts and sort
system("sort -u -nk2 -nk3 -nk1 $rABAoutputfile > $fABAoutputfile");
system("sort -u -nk2 -nk3 -nk1 $BseRIoutputfile > $fabaBseRIoutputfile");
