# Shorten a fastq file i.e. 150 read length to our normal 76
                                                                                                                
use warnings;
use Getopt::Long;
 
GetOptions (   
	    "Input=s" => \$Input,
		"LengthMake=i" => \$LengthMake
    )   
 or die("usage: All Arguments in this order are needed: -Input is the input fastq file  && -LengthMake is the number of base pairs to trim to\n");


open($fastafileR1, $Input);

$AutoNameGen = substr($Input,0,-6)."-Trimed-".$LengthMake.".fastq";
open($outputR1,'>',$AutoNameGen);

while ($header = <$fastafileR1>)
{
    $sequence = <$fastafileR1>;
    $plus     = <$fastafileR1>;
    $qual     = <$fastafileR1>;

    $sequence2 = substr($sequence,0,$LengthMake);
    $qual2 = substr($qual,0,$LengthMake);
    
    print $outputR1 ("$header$sequence2\n$plus$qual2\n");
	    
	
    
}
