#!/usr/bin/perl
#This script reprots the number of components, contigs, seqs, and the length of each sequence.
#inputfile is a sequence file
#Usage= *.pl Trinity.fasta
#Nima Rafati 20120912
my $fastaFile=$ARGV[0];
my %compHash=();
$compCntr=0;
$contigCntr=0;
$seqCntr=0;

open (inFastaFile, $fastaFile) || die print "No scaffolds file\n";
while(<inFastaFile>)
{
	
	if (/^\>/)
	{
		chomp($_);
		$_=~ m/\>(.*)\_(.*)\_(.*)/;#\-.*\=(\d+)/;
		$compId=$1;
		$contigId=$2;
		$seqId=$3;
		$seqLength=$4;
#		print "$compId\t$contigId\t$seqId\t$seqLength";<STDIN>;
		$components{$compId}{$contigId}{$seqId} = $seqLength;
	}
}
close inFastaFile;
foreach my $compK (keys  %components)
{
	$compCntr++; #countign the number of components.
	foreach my $contigK (sort keys  %{$components{$compK}})
	{
		$contigCntr++; #counting number of contigs per component.
		foreach my $seqK (keys %{$components{$compK}{$contigK}})
		{
			#adding seqlengths in a hash for statistical report
			$seqValue=$components{$compK}{$contigK}{$seqK};
			$seqKey="$compK\_$contigK\_$seqK";
			$seqLengthHash{$seqKey}=$seqValue; #header=seqlength since there are some contigs with same length.
			$seqCntr++;
			$sumLength=$sumLength+$seqValue;
			#Finding max length of seqs among seqs in one contig.
			if ($seqValue>=$maxLength)
			{
				$maxLength=$seqValue;
			}
			if (!exists $thisSeqHash{$seqValue})
			{
				$thisSeqHash{$seqValue}=1;
			}
			else
			{
				$thisSeqHash{$seqValue}=$thisSeqHash{$seqValue}+1;
			}
		}
		my ($seqLengthID,$maxCountSeqLength) = each %thisSeqHash;
		while(my ($k,$v)=each %thisSeqHash)
		{
			if ($v>=$maxCountSeqLength)
			{
				$maxCountSeqLength=$v;
				$seqLengthID=$k;
			}
		}
#		"len=$seqLengthID:$maxCountSeqLength";
		%thisSeqHash=();
		$meanLength=($sumLength/$seqCntr);#calculating mean length of seqs in one contig
		$meanLength=sprintf "%.2f", $meanLength;
		$sumLength=0;
		$returnSeq="$seqCntr\:$maxLength\:$meanLength\:len=$seqLengthID-$maxCountSeqLength";
#		print $returnSeq;<STDIN>;
		$compContigKey="$compK\_$contigK";
		$seqCountHash{$compContigKey}=$returnSeq; #generating a hash of seqcount within each contig of components.
		$seqCntr=0;
		$sumLength=0;
		$meanLength=0;
		$maxLength=0;
		$seqLengthID=0;
		$maxCountSeqLength=0;
		$k=0;
		$v=0;
	}
	$compCountHash{$compK}=$contigCntr;	
	$contigCntr=0;
}

print "Number of components: $compCntr\n";
#Print number of sequences in each contig, length mean, max length, and thelenght of sequence with highest observation.
open (outSeqCount, ">>$fastaFile.seqCount") || die print "Error: was not able to create seqCount file\n";
foreach my $seqCountK (keys %seqCountHash)
{
	@seqArr=split("\:",$seqCountHash{$seqCountK});
	print  outSeqCount "$seqCountK\t";
	foreach (@seqArr)
	{
		$i++;
		if ($i==4)
		{
			$_=~ m/(.*)\-(\d+)/;
			print outSeqCount "$1\t$2";
		}
		else
		{
			print outSeqCount "$_\t";
		}
	}
	$i=0;
	print outSeqCount "\n"
}
close outSeqCount;
#Print number of contigs per component.
open (outCompCount, ">>$fastaFile.compCount") || die print "Error: was not able to create compCount file\n";
foreach my $compCountK (keys %compCountHash)
{
	print outCompCount "$compCountK\t$compCountHash{$compCountK}\n";
}
close outCompCount; 
