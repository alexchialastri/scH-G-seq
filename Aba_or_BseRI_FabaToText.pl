# Takes in the faba file of BseRI or AbaSI and converts it to a simpler text file with 1 base numbering for cell and chromosome.



$filename = $ARGV[0];
open($file, '<', $filename);

if (!($ARGV[0])){
    die "usage: All Arguments in this order are needed: Argument 0 is the input faba file for BseRI or AbaSI\n";
}

$outputfilenameCG      = substr($filename,0,-5).".txt";
open($outputCG, '>', $outputfilenameCG);



	
while ($line = <$file>)
{
    @x = split("\t",$line);
    #$ReadsPerCell[$x[0]-1]++;
	$ReadsPerCell[$x[0]]++;
	$CellNum = $x[0]+1;
	$ChrNum = $x[1]+1;
	if($ChrNum<=25){
		print $outputCG ("$CellNum\t$ChrNum\t$x[2]\t$x[3]\n");
	}
}
