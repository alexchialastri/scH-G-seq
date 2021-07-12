# Extract AluI (or 5mC) reads based on 96 barcodes from fastq files (containing 3 UMIs too)
                                                                                                                
use warnings;
use Getopt::Long;
 
GetOptions (   
	    "FASTQ_R1=s" => \$FASTQ_R1,
		"ALUI_BC=s" => \$AluIBCpath,
		"ABASI_BC=s" => \$ABASIBCpath
    )   
 or die("Error in command line arguments\n usage: All Arguments in this order are needed: --FASTQ_R1 is the input name of your input R1 fastq file && ALUI_BC is the file path to the AluI (or MspJI) barcodes && ABASI_BC is the file path to the ABASI barcodes )\n");
 
$scTHseq = 1;
open($fastafileR1, "$FASTQ_R1");

my $outputR1In = substr($FASTQ_R1,0,-6)."-AluI.fastq";

open($outputR1, '>', "$outputR1In");


# READ in AluI Barcodes 
open($AluIBCfile, "$AluIBCpath"); #opens AluI barcode file
$i = 0;
while ($BCline2 = <$AluIBCfile>)
{
chomp $BCline2;   
@AluIBarC = split("\t",$BCline2);
$AluIBarCode[$i] = $AluIBarC[0];
$i++;
}
  


# Read in ABASI Barcodes if is part of Joint sequencing
if ($scTHseq == 1) {
    open($ABASIBCfile, "$ABASIBCpath"); #opens ABASI barcode file
    $i = 0;
		while ($BCline3 = <$ABASIBCfile>)
		{
		chomp $BCline3;   
		@ABASIBarC = split(";",$BCline3);
		$ABASIBarCode[$i] = $ABASIBarC[0];
		$i++; 
		}
 }    


while ($header = <$fastafileR1>)
{
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;


    $AluIBCSection = substr($sequence,3,8); #Get AluI Barcode Position sequence
    $ABASIBCSection = substr($sequence,0,6); #Get ABASI Barcode Position sequence
	
    
	
	
	$AluICount = 0;
	$ABASICount = 0;
	
	
		 #Check if AluI Barcode is a match
	foreach (@AluIBarCode){
		if ($_ eq $AluIBCSection){$AluICount++;}
	}
	
	
	if ($scTHseq == 1) { #Check if ABASI Barcode is a match
	foreach (@ABASIBarCode){
		if ($_ eq $ABASIBCSection){$ABASICount++;}
		}
	}
	
	$NonAluISEQ = $ABASICount; #Match count of all non AluI Barcodes
	
	#Print out New Fastq File if barcode match only AluI
	if ( $AluICount > 0){
		if ( $NonAluISEQ == 0 ){
			print $outputR1 ("$header$sequence$plus$qual");
		}
	}
	

}
