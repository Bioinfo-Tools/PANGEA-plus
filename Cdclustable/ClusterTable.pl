# !/usr/bin/perl -w
#
# ClusterTable.pl
# Created by David Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: October 29, 2009
#################################################################
#
#	Parameters:	
#		-c fas.clstr.clstr input file
#		-n  number of samples
#		-o name of output file
#
#################################################################


if($#ARGV != 5)
{
	print "Please enter: -c your fas.clstr.clstr input file -n number of samples -o the name of your output file.\n";
	exit;
}

print "Opening $ARGV[1]...";

unless (open(INPUT, $ARGV[1]))       #tries to open file
{
	print "Unable to open $ARGV[1]\nMake sure you entered the extension when entering the file name.";
	exit;
}

@input = <INPUT>;
print "successful.\nCreating $ARGV[5]...";
open FILE, ">$ARGV[5]" or die $!;
print "successful.\nGenerating table...";
@samples = ();
for($a = 0; $a < $ARGV[3]; $a++)
{
	unshift(@samples, "0");
}

$size = scalar @input;
for($a = 0; $a < $size; $a++)
{
	if($input[$a] =~ /Cluster/)
	{
		if($a > 0)
		{
			for($b = 0;$b<$ARGV[3];$b++)
			{
				print FILE $samples[$b]."\t";
				$samples[$b] = 0;
			}
			print FILE "\n";
		}
		print FILE substr($input[$a], 9, length($input[$a]) - 10)."\t";
	}
	else
	{
		$samp = substr($input[$a], index($input[$a], ">") + 1, 2);
		$samp--;
		$samples[$samp]++;
		if($input[$a] =~ /\*/)
		{
			print FILE substr($input[$a], index($input[$a], ">"), index($input[$a], "\.")-index($input[$a], ">"))."\t";
		}
	}
}
for($b = 0;$b<$ARGV[3];$b++)
{
	print FILE $samples[$b]."\t";
	$samples[$b] = 0;
}
			print FILE "\n";
print "successful.\nFinished!\n";
close FILE;
close INPUT;
