##### This code requires sorted file input (sorted on cell barcode, chr and AluI coordinate!!!#####

#Removes PCR duplicates

use warnings;

my $fileinput = $ARGV[0];
my $fileoutput = $ARGV[1];

open($file, "$fileinput");
open($output, '>', "$fileoutput");

#open($file, "/hpc/hub_oudenaarden/s.dey/Hyd_RNA_E14_sc_20141112/E14HYDRNA_AbaReads/E14HYDRNA_AbaReads_se_half5hmc_correctpos_PseudoBulk_simplified_sort.txt");
#open($output, ">/hpc/hub_oudenaarden/s.dey/Hyd_RNA_E14_sc_20141112/E14HYDRNA_AbaReads/E14HYDRNA_AbaReads_se_half5hmc_correctpos_PseudoBulk_simplified_sort_rmdup.txt");

$line1       = <$file>;
#$count       = 0;
#$dup[$count] = 0; 
print $output ("$line1");

#for ($i=0;$i<=8;$i++)
#{
while ($line2 = <$file>)
{
    #$count++;
    #$line2 =<$file>;
    if ($line1 eq $line2)
    {
	#$dup[$count] = 1;
    }
    else
    {
	print $output ("$line2");
	#$dup[$count] = 0;
    }
    $line1 = $line2;
}

#print "@dup\n";

close ($file);
