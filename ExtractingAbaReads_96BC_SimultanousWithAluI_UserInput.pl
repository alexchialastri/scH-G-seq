# Extract Aba reads based on 96 barcodes from fastq files (containing Aba, and AluI (or MspJI) reads)                                                                                                 
use warnings;

 
my ($START_DIR, $OUT_NAME_R1, $OUT_NAME_R2, $FASTQ_R1, $FASTQ_R2, $BARCODES, $MSPJIBARCODES)= @ARGV;



if (!($START_DIR && $OUT_NAME_R1 && $OUT_NAME_R2 && $FASTQ_R1 && $FASTQ_R2 && $BARCODES && $MSPJIBARCODES)){
    die "usage: All Arguments in this order are needed: -START_DIR is the starting directory && -OUT_NAME_R1 is the output name of your output R1 fastq file && -OUT_NAME_R1 is the output name of your output R2 fastq file && -FASTQ_R1 is the input name of your input R1 fastq file && -FASTQ_R2 is the input name of your input R2 fastq file && -BARCODES is a .csv file of cell specific AbaSI barcodes && -MSPJIBARCODES is a .csv file of AluI (or MspJI) barcodes)\n";
}

$FileName1=join "",$START_DIR,"/",$FASTQ_R1;
$FileName2=join "",$START_DIR,"/",$FASTQ_R2;
$FileName3=join "",$START_DIR,"/",$OUT_NAME_R1,".fastq";
$FileName4=join "",$START_DIR,"/",$OUT_NAME_R2,".fastq";



open($fastafileR1, "$FileName1");
open($fastafileR2, "$FileName2");

open($outputR1, ">$FileName3");
open($outputR2, ">$FileName4");

open($BC, "$BARCODES");
open($BCMSPJI, "$MSPJIBARCODES");




$i = 0;
while ($BCline = <$BC>)
{
    @BarC = split(";",$BCline);
    $BarCode[$i] = $BarC[0];
    $i++;
    #print "$BarCode[1]\n";                                                                                                                                                                     
                                                                                                                                                                               
}

$i = 0;
while ($BClineMSPJI = <$BCMSPJI>)
{
    chomp($BClineMSPJI);
    $BarCodeMSPJI[$i] = $BClineMSPJI;
        
   # print "$BarCodeMSPJI[$i]\n";
$i++;
   
}
                  

#$BarCode1 = "CATCACGC";

while ($header = <$fastafileR1>)
{
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;

    $headerR2   = <$fastafileR2>;
    $sequenceR2 = <$fastafileR2>;
    $plusR2     = <$fastafileR2>;
    $qualR2     = <$fastafileR2>;

$SKIPFlag=0;
    for ($i=0;$i<96;$i++)
    {
	$start = index($sequence,$BarCode[$i]);


	if ($start == 0)
	{

	    for ($j=0;$j<96;$j++){
		$MatchMSPJI = index($sequence,$BarCodeMSPJI[$j]);
		if ($MatchMSPJI == 3){
		    $SKIPFlag=1;
		    last;
		}}
        if ($SKIPFlag == 0)
	{
	    print $outputR1 ("$header$sequence$plus$qual");
	    print $outputR2 ("$headerR2$sequenceR2$plusR2$qualR2");
	    last;
	}
	    #else {
		#print ("$header$sequence$plus$qual"); }
}
$SKIPFlag=0;
    }
}
