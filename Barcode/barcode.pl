# !/usr/bin/perl -w
#
# barcode.pl
# Written by: David B. Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: January 12, 2010
#################################################
#
#	Parameters:
#		-s sequences input
#		-b barcodes in txt file
#	Options:
#		-d directory of output
#
#################################################

if ($#ARGV < 3)
{
	print "Please enter the -s sequences and -b barcodes.\n";
	exit;
}
elsif ($#ARGV >5)
{
	print "Too many inputs\n";
	exit;
}

for($a = 0; $a < $#ARGV + 1; $a++)
{
	if ($ARGV[$a] eq "-s")
	{
		$seqIn = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-b")
	{
		$barIn = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-d")
	{
		$dirName = $ARGV[$a+1];
	}
}

if(!defined($seqIn) || !defined($barIn))
{
	print "Must enter at least -s sequences and -b barcodes.\n";
	exit;
}

print "Opening $seqIn...";
unless (open(SEQIN, $seqIn))       #tries to open file
{
	print "Unable to open $seqIn\nMake sure you entered the extension when entering the file name.";
	exit;
}
print "Successful.\nOpening $barIn...";

unless (open(BARIN, $barIn))       #tries to open file
{
	print "Unable to open $barIn\nMake sure you entered the extension when entering the file name.";
	exit;
}
@barNF = <BARIN>;
close BARIN;
print "Successful.\nCreating Files...\n";

if(defined($dirName))
{
	open NOBAR, ">$dirName/"."nobar.fas" or die $!;
	open COUNT, ">$dirName/"."barcode_count.txt" or die $!;
}
else
{
	open NOBAR, ">nobar.fas" or die $!;
	open COUNT, ">barcode_count.txt" or die $!;
}
@barCount = ();
@storedLines = ();
$count = 0;
$countNoBar = 0;
$countThrownOut = 0;
$barSize = scalar @barNF;
for($a = 0; $a < $barSize; $a++)
{
	$dex = index($barNF[$a], "\t") + 1;
	$line = $barNF[$a];
	$line =~ s/\s+$//;
	$line =  uc substr($line, $dex, length($line) - $dex);
	push(@bar, $line);
	push(@barCount, 0);
	push(@out, []);
}
while($line = <SEQIN>)
{
	if($line =~ />/)
	{
		$count++;
		$storeSize = scalar @storedLines;
		if($storeSize > 0)
		{
			find_and_store_barcode(@storedLines);
			@storedLines = ();
		}
		push(@storedLines, $line);
	}
	else
	{
		if($line ne "\n")
		{
			push(@storedLines, $line);
		}
	}
}
close SEQUIN;
find_and_store_barcode(@storedLines);

print "Found $count sequences.\nFound $countNorm barcodes at the beginning.\nFound no barcodes in $countNoBar sequences.\nThrew out $countThrownOut barcodes based on length.\nWriting Files...";
$min = $#{$out [0]} + 1;
$count = 0;
for($a = 0; $a < $barSize; $a++)
{
	$size = $#{$out [$a]}+1;
	$num = $a + 1;
	if($num < 10)
	{
		print COUNT "0"."$num"."\t"."$bar[$a]"."\t"."$barCount[$a]"."\n";
	}
	else
	{
		print COUNT "$num"."\t"."$bar[$a]"."\t"."$barCount[$a]"."\n";
	}
	if(defined($dirName))
	{
		open FILEOUT, ">$dirName/"."bar".substr($barNF[$a], 0, index($barNF[$a], "\t")).".fas" or die $!;
	}
	else
	{
		open FILEOUT, ">bar".substr($barNF[$a], 0, index($barNF[$a], "\t")).".fas" or die $!;
	}
	for($b = 0; $b < $size; $b++)
	{
		print FILEOUT $out[$a][$b];
	}
	if($barCount[$a] > 100)
	{
		if($barCount[$a] < $min)
		{
			$min = $barCount[$a];
		}
	}
	close FILEOUT;
}
open MIN, ">$dirName"."/info.txt" or die $!;
print MIN $min;
print "Successful.\nFinished!\n";

close MIN;
close NOBAR;
close COUNT;

#######################################	  SUBROUTINES	#######################################

sub find_and_store_barcode()
{
	my $first = 1;
	my $thisBar = 0;
	my $place = 0;
	my $len;
	my $bases = 0;
	my $storeSize = scalar @storedLines;
	for($b = 0; $b < $storeSize; $b++)
	{
		if($storedLines[$b] =~ />/)
		{
			$first = 1;
			$thisBar = -1;
			$bases = 0;
			$place = 0;
		}
		else
		{
			if($first == 1)
			{
				for($a = 0; $a < $barSize; $a++)
				{
					if($storedLines[$b] =~ /^$bar[$a]/)
					{
						$thisBar = $a;
						$a = $barSize;
						$place = 1;
						$bases = return_bases($thisBar, $place);
					}
				}
				$first == 0;
			}
			if($thisBar == -1)
			{
				for($a = 0; $a < $barSize; $a++)
				{
					if($storedLines[$b] =~ /$bar[$a]/)
					{
						$thisBar = $a;
						$a = $barSize;
						$place = $b;
						$bases = return_bases($thisBar, $place);
					}
				}
			}
		}
	}
	if($thisBar == -1)
	{
		$countNoBar++;
		while($temp = shift(@storedLines))
		{
			print NOBAR $temp;
		}
		print NOBAR "\n";
	}
	elsif($bases < 100)
	{
		$countThrownOut++;
		@storedLines = ();
	}
	else
	{
		$trimmed = 0;
		$barCount[$thisBar]++;
		$countNorm++;
		$temp = shift(@storedLines);
		push @{ $out[$thisBar] }, ">".substr($barNF[$thisBar], 0, index($barNF[$thisBar], "\t"))."_".substr($temp, 1);
		$place--;
		while($temp = shift(@storedLines))
		{
			if($place == 0)
			{
				if($trimmed == 0)
				{
					$len = length($bar[$thisBar]);
					push @{ $out[$thisBar] }, substr($temp, index($temp, $bar[$thisBar]) + $len);
					$trimmed = 1;
				}
				else
				{
					push @{ $out[$thisBar] }, $temp;
				}
			}	
			else
			{
				$place--;
			}
		}
		push @{ $out[$thisBar] }, "\n";
	}
}

sub return_bases
{
	my ($thisBar, $place) = @_;
	my $len = length($bar[$thisBar]);
	my $bases = 60 - (index($storedLines[$place], $bar[$thisBar]) + $len);
	for($b = $place + 1; $b < $storeSize; $b++)
	{
		$bases += length($storedLines[$b]);
	}
	return $bases;
}