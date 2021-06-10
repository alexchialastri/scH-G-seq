# Extracting AluI sequences that mapped to the genome

use warnings;

my $samfileinput = $ARGV[0];
my $output1input = $ARGV[1];
my $output2input = $ARGV[2];

open($samfile, "$samfileinput");
open($output, '>', "$output1input");
open($output2, '>', "$output2input");

#open($samfile, "/home/dpodorefsky/sequencing_data/5mC_DamID/5mC_DamID_03/5mC_DamID_03-se.sam");
#open($output, ">/home/dpodorefsky/sequencing_data/5mC_DamID/5mC_DamID_03/5mC_DamID_03-se_MappedReads_test.txt");
#open($output2, ">/home/dpodorefsky/sequencing_data/5mC_DamID/5mC_DamID_03/5mC_DamID_03-se_correctpos_stringent_test.txt");

#open($samfile, "/hpc/hub_oudenaarden/s.dey/DamID_RNA_EBHyd_sc_20160523/EBDLRNALib2_DamIDReads/EBDLRNALib2_DamIDReads_se.sam");
#open($output, ">/hpc/hub_oudenaarden/s.dey/DamID_RNA_EBHyd_sc_20160523/EBDLRNALib2_DamIDReads/EBDLRNALib2_DamIDReads_se_MappedReads.txt");
#open($output2, ">/hpc/hub_oudenaarden/s.dey/DamID_RNA_EBHyd_sc_20160523/EBDLRNALib2_DamIDReads/EBDLRNALib2_DamIDReads_se_correctpos_stringent.txt");

##open($output3, ">/home/s.dey/GroupVanOudenaarden/s.dey/Hydroxy_RNA_E14_20140919/DylanNORTHL_AbaReads/DylanNORTHL_AbaReads_MappedReads_NucleotideFreq.csv");
##open($output4, ">/home/s.dey/GroupVanOudenaarden/s.dey/Hydroxy_RNA_E14_20140919/DylanNORTHL_AbaReads/DylanNORTHL_AbaReads_MappedReads_CGFreq.csv");

$tot = 0;
$cor = 0;

#$w      = -1;
#@freqCG = ();

#for($j=0;$j<=64;$j++)
#{
#    for($k=0;$k<=3;$k++)
#    {
#$freq[$j][$k] = 0;
#    }
#}

while ($samline = <$samfile>)
{
    if (substr($samline,0,2) eq "NS")
    {
	@R1 = split(" ",$samline);
	if (($R1[1] == 0) && ($R1[12] eq "XT:A:U"))
	{
	    $chr          = substr($R1[2],3);
	    $coord        = $R1[3];
	    $actual_coord = $coord;
	    $bc           = substr($R1[11],8,8);
	    $bc           =~ tr/ACGTacgtNn/ACGTACGTNN/;
	    $umi          = substr($R1[11],5,3);
	    $umi          =~ tr/ACGTacgtNn/ACGTACGTNN/;
	        #$start       = $coord - 65;
	        #$stop        = $coord + 65;
	    $seq          = $R1[9];
	    $tot++;
	        if (substr($seq,0,2) eq "CT")
		{
		    $cor++;
		    $cor_cut = 1;
		    print $output2 ("$bc $chr $coord $umi $R1[1] $seq $actual_coord\n");          #start gives 65bp upstream and stop gives 65bp downstream of cut site
		}
	        else
		{
		    $cor_cut = 0;
		}
	        #print ("$chr $coord $seq @poss\n");
	        #last;
	    print $output ("$bc $chr $coord $umi $R1[1] $seq $actual_coord $cor_cut\n");
	        
	        #for($j=0;$j<=64;$j++)           #Looking at frequencies of nucleotides in the first 65 bases
	        #{
	    #$nuc = substr($seq,$j,1);
	    #if ($nuc eq "A")
	    #{
	    #    $freq[$j][0] = $freq[$j][0] + 1;
	    #}
	    #if ($nuc eq "C")
        #        {
        #            $freq[$j][1] = $freq[$j][1] + 1;
        #        }
	    #if ($nuc eq "G")
        #        {
        #            $freq[$j][2] = $freq[$j][2] + 1;
        #        }
	    #if ($nuc eq "T")
        #        {
        #            $freq[$j][3] = $freq[$j][3] + 1;
        #        }
	    #    }
	}
        if (($R1[1] == 16) && ($R1[12] eq "XT:A:U"))
        {
            $chr          = substr($R1[2],3);
            $coord        = $R1[3];
	    $actual_coord = $coord + length($R1[9]) - 1;
	    $bc           = substr($R1[11],8,8);
	    $bc           =~ tr/ACGTacgtNn/ACGTACGTNN/;
	    $umi          = substr($R1[11],5,3);
	    $umi          =~ tr/ACGTacgtNn/ACGTACGTNN/;
	        #$start       = $coord;
	        #$stop        = $coord + 130;
            $seq          = $R1[9];
	    $seqrc        = reverse($seq);
            $seqrc        =~ tr/ACGTacgtN/TGCAtgcaN/;
	    $tot++;
	        if (substr($seqrc,0,2) eq "CT")
		{
		    $cor++;
		    $cor_cut = 1;
		    print $output2 ("$bc $chr $coord $umi $R1[1] $seqrc $actual_coord\n");
		}
            else
            {
                $cor_cut = 0;
            }
            #print ("$chr $coord $seqrc @poss\n");
            #last;
	    print $output ("$bc $chr $coord $umi $R1[1] $seqrc $actual_coord $cor_cut\n");
	        
            #for($j=0;$j<=64;$j++)         #Looking at frequencies of nucleotides in the first 65 bases                                                                               
            
            #{
	        #$nuc = substr($seqrc,$j,1);
	        #if ($nuc eq "A")
            #    {
            #        $freq[$j][0] = $freq[$j][0] + 1;
            #    }
            #    if ($nuc eq "C")
            #    {
            #        $freq[$j][1] = $freq[$j][1] + 1;
            #    }
            #    if ($nuc eq "G")
            #    {
            #        $freq[$j][2] = $freq[$j][2] + 1;
            #    }
            #    if ($nuc eq "T")
            #    {
            #        $freq[$j][3] = $freq[$j][3] + 1;
            #    }
            #}
	}
    }
}

#for($j=0;$j<=64;$j++)
#{
#    #print ("\n");
#    for($k=0;$k<=3;$k++)
#    {
#        #print ("$freq[$j][$k] ");
#print $output3 ("$freq[$j][$k] ");
#    }
#    print $output3 ("\n");
#}
##print ("\n");

##print ("@freqCG\n");
##print $output4 ("@freqCG");

#for ($j=0;$j<scalar(@freqCG);$j++)
#{
#    print $output4 ("$freqCG[$j]\n");
#}

print ("Total Mapped Reads = $tot , AluI-cut reads = $cor\n");
