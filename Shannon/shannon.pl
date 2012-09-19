# !/usr/bin/perl -w
#	
# shannon.pl
# Written by: David Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: January 14, 2010
#######################################################
#
#	Parameters:
#		-t hybrid table output file
#		-n number of samples
#		-o output name
#
#######################################################

if ($#ARGV != 5)
{
	print "Please make sure that you are entering the correct parameters.\n";
	exit;
}

#Organizes input
for($a = 0; $a < 5; $a++)
{
	if($ARGV[$a] eq "-t")
	{
		$a++;
		$inputTable = $ARGV[$a];
	}
	elsif($ARGV[$a] eq "-n")
	{
		$a++;
		$num = $ARGV[$a];
		$num++;
	}
	elsif($ARGV[$a] eq "-o")
	{
		$a++;
		$output = $ARGV[$a];
	}
}

unless (open(INPUT, $inputTable))       #tries to open file
{
	print "Unable to open $inputTable\nMake sure you entered the extension when entering the file name.";
	exit;
}
@in = <INPUT>;
$inSize = scalar @in;
open OUT, ">$output" or die $!;
for($a = 0; $a < $inSize; $a++)
{
	print OUT $in[$a];
}
print OUT "Shannon-Weaver Diversity Index H'(loge):";
for($b = 1; $b < $num; $b++)
{
	@data = ();
	for($c = 1; $c < $inSize; $c++)
	{
		@line = split(/\t/, $in[$c]);
		push(@data, $line[$b]);
		
	}
	$size = scalar @data;
	$sum = 0;
	$sum2 = 0;
	for($a = 0; $a < $size; $a++)
	{
		$sum += $data[$a];	
	}
	for($a = 0; $a < $size; $a++)
	{
		if($data[$a] != 0)
		{
			$data[$a] = $data[$a]/$sum;
			$data[$a] = $data[$a] * log($data[$a]);
			$sum2 += $data[$a];
		}
	}
	$sum2 = $sum2*(-1);
	$sum2 = sprintf("%.2f", $sum2);
	print OUT "\t$sum2";
}