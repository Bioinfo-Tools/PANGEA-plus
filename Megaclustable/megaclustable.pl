# !/usr/bin/perl -w
#	
# megaclustable.pl
# Written by: David Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: October 29, 2009
###########################################################################################
#
#	Parameters:
#		-m list of all megaclust output files seperated by spaces
#		-t level of taxa that megaclust was run for where 0= domain and 6= species
#		-o output table name
#
###########################################################################################

if($#ARGV < 5)
{
	print "Please enter the correct parameters.\n";
	exit;
}
$i=0;
#organizes input
$mIn = "false";
for($a = 0; $a <= $#ARGV; $a++)
{
	if($ARGV[$a] eq "-m")
	{
		$mIn = "true";
	}
	elsif($ARGV[$a] eq "-o")
	{
		$mIn = "false";
		$a++;
		$output = $ARGV[$a];
	}
	elsif($ARGV[$a] eq "-t")
	{
		$mIN = "false";
		$a++;
		$taxLevel = $ARGV[$a];
		if($taxLevel > 6 || $taxLevel < 0)
		{
			print "You must enter a number between 0 and 6 for taxonomy level where 0 = domain and 6 = species.\n";
			exit;
		}
		$taxLevel = "["."$taxLevel"."]";
	}
	elsif($mIn eq "true")
	{
		push(@megaclustFile, $ARGV[$a]); $i++;
	}
}

#$numMega = scalar @megaclustFile;
$numMega = $i;
@table = ();
@taxa = ();
$size = 0;

for($b = 0; $b < $numMega; $b++)
{
	unless (open(MEGA, $megaclustFile[$b]))       #tries to open file
	{
		print "Unable to open $megaclustFile[$b]\nMake sure you entered the extension when entering the file name.\n";
		exit;
	}
	@file = ();
	for($a = 0; $a < $size; $a++)
	{
		push(@file, "0");
	}
	while($line = <MEGA>)
	{
		chomp($line);
		$loc = index($line, $taxLevel);
		if($loc > -1)
		{
			$found = "false";
			$loc += 3;
			$end = index($line, ";", $loc);
			if($end == -1)
			{
				$end = index($line, ",", $loc);
			}
			$name = substr($line, $loc, $end - $loc);
			$numStart = index($line, ",", $end) + 1;
			#$numEnd = index($line, "\n", $numStart);
			$num = substr($line, $numStart, length($line) - $numStart);
			for($a = 0; $a < $size; $a++)
			{
				if($taxa[$a] eq "$name")
				{
					$found = "true";
					$file[$a]+= $num;
				}
			}
			
			if($found eq "false")
			{
				push(@taxa, $name);
				push(@file, $num);
				$size++;
			}
		}
	}
	close MEGA;
	push(@table, [@file]);
}

unshift(@table, [@taxa]);
open OUTPUT, ">$output" or die $!;
for($a = 1; $a <= $numMega; $a++)
{
	print OUTPUT "\t$a";
}
for($a = 0; $a < $size; $a++)
{
	print OUTPUT "\n";
	for($b = 0; $b < $numMega + 1; $b++)
	{
		if($table[$b][$a] eq "")
		{
			$table[$b][$a] = 0;
		}
		print OUTPUT "$table[$b][$a]\t";
	}
}


