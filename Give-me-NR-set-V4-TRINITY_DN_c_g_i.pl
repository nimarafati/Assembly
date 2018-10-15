#!/usr/bin/perl
#This script gives you the non-redundant set using blast output & a file including length of each sequence
#Usage:
#
#Nima Rafati 20131124:20131129

use Getopt::Long;
my $unaligned_flag="T";
my $query="";
my $target="";
my %lengthHash=();
my $evalue=100;
my %countQHash=();
my %countTHash=();
my $evalueCut=1e-10;
my $alignHash=();

$usage="Use Give-me-NR-set.pl to get non-redundant list of \"Trinity\" sequences based on their alignment.
Nima Rafati 2013-11-29
Required:
-blast alignment file from blastn in tabular format
-length a list of sequences with their length (sequence_name\tsequence_length)
		Please note that the sequence names in length and blast file should be the same pattern.
-output (optional) name the output file. Otherwise blast_file-NR.list
-unaligned (T/F) This flag reports sequences which have not been observed in blastn output.
			By default it is true (T).
-evalue	cutoff evalue to use for filtering the alignments (default = 1e-10).
-h prints this help\n";

&GetOptions ('h'=>\$help_flag,
			 'blast=s' =>\$alignFile,
			 'length=s' =>\$lengthFile,
			 'output=s' =>\$output_file,
			 'unaligned=s' =>\$unaligned_flag,
			 'evalue=f' =>\$evalueCut);
if($flag_hel)
{
	die $usage;
}

unless($alignFile&&$lengthFile)
{
	die $usage
}

if ($output_file eq "")
{
	$output_file="$alignFile"."-NR.list";
}

 if ($unaligned_flag eq "T")
 {
 	$out_unaligned_file="unaligned.list";
 }



open(inF1,$lengthFile) || die print "No input file";
while(<inF1>)
{
	#Save lengths in a hash
    chomp($_);
    @lineArr1=split("\t",$_);
	$lengthHash{$lineArr1[0]}=$lineArr1[1];
}
close inF1;
open(inF2,$alignFile) || die print "No input file";
while(<inF2>)
{
    chomp($_);
    @lineArr2=split("\t",$_);
    $query=$lineArr2[0];
    $target=$lineArr2[1];
    $evalue=$lineArr2[10];
	my ($comp_Query,$comp_Query_c)=NameParser($query);
	my ($comp_Target,$comp_Target_c)=NameParser($target);
#print $comp_Query,"\t",$comp_Target;<STDIN>;
    #Filter by Evalue
    if ($evalue<=$evalueCut)
    {
    	#Save self-alignments
		if($query eq $target) 
		{
			$uniqueHash{$query}=$target;
		}
    	#put multiple hits in one hash with e-value if they are different
    	elsif($query ne $target) 
	    {
	    	#Check if the sequences are from the same components
    		if ($comp_Query eq $comp_Target)
	    	{
	    		#Save alignments of the same components and determine which one (query || subject) is larger.
				if (!exists $alignHash{$query}{$target} || !exists $alignHash{$target}{$query})
				{
					if($lengthHash{$query}>$lengthHash{$target} && !exists $deletedHash{$target})
					{
						$alignHash{$query}{$target}=$query;
						$deletedHash{$target}=$target;
					}
					elsif($lengthHash{$query}<$lengthHash{$target} && !exists $deletedHash{$query} )
					{
						$alignHash{$query}{$target}=$target;
						$deletedHash{$query}=$query;
					}
					elsif($lengthHash{$target}>$lengthHash{$query} && !exists $deletedHash{$query})
					{
						$alignHash{$target}{$query}=$target;
						$deletedHash{$query}=$query;
					}
					elsif($lengthHash{$target}<$lengthHash{$query}&& !exists $deletedHash{$target})
					{
						$alignHash{$target}{$query}=$query;
						$deletedHash{$target}=$target;
					}
				}
    		}
    	}
    }
}
close inF2;

#Choose the largest sequence in cases with multiple hits:
foreach my $kQ (keys %alignHash) #Read the saved alignments of the same components
{
	foreach my $kT (keys %{$alignHash{$kQ}})
	{
		if(exists $tmpHash{$kQ}) #Check if saved sequence is larger than other sequence which it has alignment with
		{
			if($lengthHash{$tmpHash{$kQ}}>$lengthHash{$kT})
			{
				$deletedHash{$kT}=$kT; #sequences which have alignment with a larger sequence of the same comp. and need to be removed.
			}
			else
			{
				$tmpHash{$kQ}=$kT;
				$deletedHash{$kQ}=$kQ;
			}
		}
		else #Save in a temporary hash if the query is larger than target and pass the smaller sequence to $deletedHash
		{
			if($lengthHash{$kQ}>$lengthHash{$kT})
			{
				$tmpHash{$kQ}=$kQ;
				$deletedHash{$kT}=$kT;
			}
			else
			{
				$tmpHash{$kQ}=$kT;
				$deletedHash{$kQ}=$kQ;
			}
		}
	}
	$nowFinal=$tmpHash{$kQ}; 
	if (!exists $finalHash{$nowFinal}) #Pass the larger sequence into finalHash
	{
		$finalHash{$nowFinal}=$nowFinal;
	}
	else
	{
		if($finalHash{$nowFinal} ne $nowFinal)
		{
			$finalHash{$nowFinal}=$nowFinal;
		}
	}
	$size=scalar(keys(%tmpHash));
	%tmpHash=();
}

##Unique hit report:
open (outF, ">$output_file");
foreach my $k (keys %uniqueHash)
{
	if ($uniqueHash{$k} && !exists $deletedHash{$k} && !exists $finalHash{$k})
	{
		print outF "$k\n";
	}
}

##Multiple hits report:
foreach my $kout (sort keys %finalHash)
{
	if(!exists $deletedHash{$kout})
	{
		print outF  "$kout\n";
	}
}
close outF;

##Unaligned sequences:
if ($unaligned_flag eq "T")
{
	open (outFU, ">$out_unaligned_file");
	foreach my $kU (keys %lengthHash)
	{
		if(!exists $finalHash{$kU} && !exists $uniqueHash{$kU})
		{
			print outFU "$kU\n";
		}
	}
close outFU;
}
################################subroutines
sub NameParser ($)
{
	$inName=shift;
	$inName=~ m/TRINITY_(DN\d+)_c(\d+)_g(\d+)_i(\d+).*/;
	$DN=$1;
	$comp=$2;
	$g=$3;
	$iso=$4;
#	$outComp="comp".$comp;
	$outComp=$DN;
	$outComp_c=$outComp;
	return($outComp,$outComp_c);
}

sub CheckLength($$)
{
	my ($nowLength,$hashLength)=@_;

	if($nowLength>$hashLength)
	{
		$outLength=$nowLength;
	}
	else
	{
			$outLength="HASH";
	}
	return ($outLength);
}
