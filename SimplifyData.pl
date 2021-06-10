#Extract only the most important info from the file that identifies AluI sites (including UMIs)

use warnings;

my $fileinput = $ARGV[0];
my $fileoutput = $ARGV[1];
my $barcode = $ARGV[2];

open($file, "$fileinput");
open($output, '>', "$fileoutput");

#open($file, "/hpc/hub_oudenaarden/s.dey/DamID_RNA_EBHyd_sc_20160523/EBDLRNALib2_DamIDReads/EBDLRNALib2_DamIDReads_se_correctpos_stringent.txt");
#open($output, ">/hpc/hub_oudenaarden/s.dey/DamID_RNA_EBHyd_sc_20160523/EBDLRNALib2_DamIDReads/EBDLRNALib2_DamIDReads_se_correctpos_stringent_simplified.txt");

open($bcfile, "$barcode");

%bc_hash = ();

$c = 0;
while ($bcline = <$bcfile>)
{
    chomp($bcline);
    $c++;
    $bc_hash{$bcline} = $c;
}

#for($i=0;$i<96;$i++)
#{
#    @p = each (%bc_hash);
#    print "@p\n";
#}

#for($i=0;$i<100;$i++)
#{
#    $line = <$file>;

while ($line = <$file>)
{
    @y = split(" ",$line);
    if (exists($bc_hash{$y[0]}))
    {
	$bc_num = $bc_hash{$y[0]};
	chomp($bc_num);
	print $output ("$bc_num $y[1] $y[6] $y[3]\n");
    }
##   if ($y[5] =~ /^-?\d+$/)     # if $y[5] is a integer      #this part of the script is introduced because of wrong order of columns in pervious script 
##   {
###print $output ("$bc_num $y[1] $y[5] $y[5] $y[0]\n");
##print $output ("$bc_num $y[1] $y[5] $y[5]\n");
##    }
##    else
##    {
##print $output ("$bc_num $y[1] $y[4] $y[4]\n");
##    }
}
