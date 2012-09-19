#!/usr/bin/perl -w
#
# selector.pl
# Written by: David Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: May 11, 2010
###########################################################################################
#
#	Parameters:
#		-i input file
#		-s number of sequences
#		-o output file
#
###########################################################################################

if ($#ARGV != 5)		#Check number of inputs
{
	print "Please enter the program with  -i input file -o output file -s number of sequences";
	exit;
}

for($a = 0; $a < 6; $a++)		#Organize input
{
	if($ARGV[$a] eq "-i")		#Input file name
	{
		$inputFile = $ARGV[$a+1];
	}
	elsif($ARGV[$a] eq "-o")	#Output file name
	{
		$outputFile = $ARGV[$a + 1];
	}
	elsif($ARGV[$a] eq "-s")	#number of sequences to select
	{
		$limit = $ARGV[$a + 1];
	}
}
print "Opening $inputFile...\n";

unless (open(INPUT, $inputFile))       #tries to open file
{
	print "Unable to open $inputFile\nMake sure you entered the extension when entering the file name.";
	exit;
}
print "Selecting...\n";

$line = <INPUT>;
open OUTPUT, ">$outputFile" or die $!;

if($line =~ /\n/)
{
	$limit--;
	if($line =~ /\r/)		#Newline characters are in CRLF format
	{
		&select_Normal($line, $limit);	
	}
	else				#Newline characters are in LF format
	{
		&select_Normal($line, $limit);
	}
}
elsif($line =~ /\r/)
{
	&select_CR($line, $limit);	#Newline characters are in CR format
}
else
{
	print "\nError: Either there is only one line in the input file or there is some newline character problem.\nExiting...\n";
}
print "Finished!\n";
close INPUT;
close OUTPUT;

###############################################   SUBROUTINES   ###############################################

sub select_Normal()	#CRLF or LF format
{
	my ($line, $limit) = @_;
	if($line =~ />/)
	{
		$limit--;
		print OUTPUT "$line";
	}
	while($line = <INPUT>)		#selects the correct number of sequences or stops at the end of the file
	{
		if($limit < 0)		#end if the correct number of sequences have been found
		{
			if(&check_Line($line) eq "exit")	#subroutine prints or tells it to exit
			{
				last;
			}
			while($line = <INPUT>)		#loop to find the end of the sequence
			{
				if(&check_Line($line) eq "exit")	#subroutine prints or tells it to exit
				{
					last;
				}
			}
			last;
		}
		else
		{
			if($line =~ />/)		#if new sequence, count it and print it
			{
				$limit--;
			}
			print OUTPUT "$line";
		}
	}
}

sub select_CR()		#CR format
{
	my($file, $limit) = @_;
	@lines = split(/\r/, $file);		#splits file into array of lines
	
	$size = scalar @lines;
	for($a = 0; $a < $size; $a++)
	{
		if($limit < 1)		#end if the correct number of sequences have been found
		{
			for($a = $a; $a < $size; $a++)		#loop to find the end of the sequence
			{
				$lines[$a] = $lines[$a]."\n";
				if(&check_Line($lines[$a]) eq "exit")	#subroutine prints or tells it to exit
				{
					$a = $size;
				}
			}
		}
		else
		{
			$lines[$a] = $lines[$a]."\n";
			if($lines[$a] =~ />/)		#if new sequence, count it and print it
			{
				$limit--;
			}
			print OUTPUT "$lines[$a]";
		}
	}
}

sub check_Line()
{
	my ($line) = @_;
	if($line =~ /[a-zA-Z]/)		#if line contains letters (bases) continue
	{
		if($line =~ />/)		#if new sequence, end
		{
			return "exit";
		}
		else
		{
			print OUTPUT "$line";
		}
	}
	else
	{
		return "exit";
	}
}
