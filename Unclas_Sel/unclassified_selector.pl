#!/usr/bin/perl -w
#
# unclassified_selector.pl
# Written by: David Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: November 2, 2009
####################################################################################################
#
# 	Parameters: 
#		-m megablast results
# 	     	-s sequences in .fas file
#	     	-o name of output file
#	Options:
#		-t similarity threshold at which you want to reject below (DEFAULT: 95)
#		-e e-value upper threshold (DEFAULT: -20)
#		-b bitscore lower threshold (DEFAULT: 200)
#
####################################################################################################

if ($#ARGV < 5 || $#ARGV > 11)
{
	print "Please enter the -m megablast, -s sequences, -t threshold, -e e-value upper threshold, -b bitscore lower threshold, and -o output file.\n";
	exit;
}

#DEFAULT VALUES
$thres = 95;
$eThres = exp(-20);
$bThres = 200;

for($a = 0; $a < 12; $a++)
{
	if ($ARGV[$a] eq "-m")
	{
		$mega = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-s")
	{
		$sequ = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-t")
	{
		$thres = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-e")
	{
		$eThres = $ARGV[$a+1];
		$eThres = exp($eThres);
	}
	elsif($ARGV[$a] eq "-b")
	{
		$bThres = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-o")
	{
		$out = $ARGV[$a+1];
	}
}

if(!(defined($mega) && defined($sequ) && defined($out)))
{
	print "Must have at least -m megablast -s sequences -o output file.\n";
	exit;
}

print "Opening $mega...";
unless (open(MEGA, $mega))       #tries to open file
{
	print "Unable to open $mega\nMake sure you entered the extension when entering the file name.";
	exit;
}

print "successful.\nRejecting...";
$count = 0;
%rejecthash = ();
%megahash = ();
while($megLine = <MEGA>)
{
	if($megLine eq "\n")
	{
		last;
	}
	chomp($megLine);
	$start = index($megLine, "\t") + 1;
	$start = index($megLine, "\t", $start) + 1;
	$end = index($megLine, "\t", $start);
	$percent = substr($megLine, $start, $end - $start);
	$Bstart = rindex($megLine, "\t") + 1;
	$B = substr($megLine, $Bstart);
	$Estart = rindex($megLine, "\t", $Bstart - 3) + 1;
	$Eend = $Bstart - 1;
	$E = substr($megLine, $Estart, $Eend - $Estart);
	$nameEnd = index($megLine, "\t");
	$name = substr($megLine, 0, $nameEnd);
	if($percent < $thres || $E > $eThres || $B < $bThres)
	{
		if(!(exists($megahash{$name}) || exists($rejecthash{$name})))
		{
			$rejecthash{$name} = 1;
		}
	}
	else
	{
		$megahash{$name} = 1;
		if(exists($rejecthash{$name}))
		{
			delete($rejecthash{$name});
		}
	}
}
close MEGA;
$found = 0;

print "successful.\nOpening $sequ...";
unless (open(SEQ, $sequ))       #tries to open file
{
	print "Unable to open $sequ\nMake sure you entered the extension when entering the file name.";
	exit;
}

print "successful.\nCreating $out...";
open OUT, ">$out" or die $!;

print"successful.\nPrinting...";
while($seqLine = <SEQ>)
{
	if($seqLine =~ />/)
	{
		$seqLine = substr($seqLine, 1);
		$seqLine =~ s/\s+$//;
		$found = 0;
		if(exists($rejecthash{$seqLine})) 
		{
			$found = 1;
			print OUT ">$seqLine \n";
			$count++;

		}
		if($found == 0)
		{
			$found = 1;
			if (exists($megahash{$seqLine})) 
			{
				$found = 0;
				if($megahash{$seqLine} > 1)
				{
					$megahash{$seqLine}--;
				}
				else
				{
					delete($megahash{$seqLine});
				}
			}
			else
			{
				print OUT ">$seqLine \n";
				$count++;
			}
		}		
	}
	elsif($found == 1)
	{
		print OUT $seqLine;
	}
}
print "successful.\nRejected $count sequence(s).\nFinished!\n";
close SEQ;
close OUT;