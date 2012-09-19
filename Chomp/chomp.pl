#!/usr/bin/perl -w
#
# chomp.pl
# Written by David Crabb
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

if ($#ARGV !=3)
{
	print "Please enter the program with  -i inputfile -o outputfile";
	exit;
}

for($a = 0; $a < 4; $a++)
{
	if($ARGV[$a] eq "-i")
	{
		$inputFile = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-o")
	{
		$outputFile = $ARGV[$a + 1];
	}
}

print "Opening $inputFile...";

unless (open(INPUT, $inputFile))       #tries to open file
{
	print "Unable to open $inputFile\nMake sure you entered the extension when entering the file name.";
	exit;
}

print "\n";
@input = <INPUT>;
open FILE, ">$outputFile" or die $!;
$size = scalar @input;
for($loop = 0; $loop < $size; $loop++)
{
	if($input[$loop] =~ />/)
	{
		print FILE substr($input[$loop], 0, index($input[$loop], " ") + 1)."\n";
	}
	else
	{
		print FILE $input[$loop];
	}
}
print "Writing $outputFile\n";
print "Done!\n";

close INPUT;
close FILE;
