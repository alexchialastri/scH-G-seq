# Extract Aba reads based on 96 barcodes from fastq files (containing Aba, and AluI (or MspJI) reads)                                                                                                 
use warnings;
use Getopt::Long;
 
GetOptions (   
	    "START_DIR=s" => \$START_DIR,
    "FASTQ_R1=s" => \$FASTQ_R1,
		"BARCODES=s" => \$BARCODES,
    "ALUIBARCODES=s" => \$ALUIBARCODES
    )   
 or die("usage: All Arguments in this order are needed: -START_DIR is the starting directory && -OUT_NAME_R1 is the output name of your output R1 fastq file && -FASTQ_R1 is the input name of your input R1 fastq file && -BARCODES is a .csv file of cell specific AbaSI barcodes && -ALUIBARCODES is a .csv file of AluI (or MspJI) barcodes)\n");

$MSPJIBARCODES = $ALUIBARCODES;

$FileName1=join "",$START_DIR,"/",$FASTQ_R1;

$OUT_NAME_R1 = substr($FASTQ_R1,0,-6)."-AbaSI-BseRI";
$FileName3=join "",$START_DIR,"/",$OUT_NAME_R1,".fastq";




open($fastafileR1, "$FileName1");

open($outputR1, ">$FileName3");


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
	    last;
	}
	    #else {
		#print ("$header$sequence$plus$qual"); }
}
$SKIPFlag=0;
    }
}
