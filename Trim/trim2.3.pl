# !/usr/bin/perl -w
#
# Trim2_Illumina.pl
# Written by: David B. Crabb
# Eric Triplett's Group
# University of Florida
# Last Modified: December 6, 2011
# Trim2.2.pl
# Edited by: Raquel Dias
# Pontificia Universidade Cat—lica do Rio Grande do Sul - Brazil
# Last modified: February 26, 2012
####################################################################
#
# 	Parameters:
#		-a raw illumina input file read 1
#		-b raw illumina input file read 2 (if any) 
#		-g size of GAP between paired-ends (if any) 
#		-t truncate size (if any)
#		-q quality file (in case of FASTA input)
# 	Supported formats: FASTA, FASTQ (Ilumina) and QSEQ (pired/single-end)
####################################################################

#TODO NOW
#include lc and qc cutoffs variables inputs
#include FASTQ support Illumina

#TODO later
#include FASTQ support for Solexa sequencer
#include FASTQ support for Sanger sequencer
#FASTQ Paired-end support

#DEFAULT VALUES
$LENGTH_CUTOFF = 70;
$QUALITY_CUTOFF = 20;
$TRUNCATE = 11;
$GAPSIZE = 189;

# Example:
# If Length < 70; throw out.
# If throw out 1; then throw out both.
# Delete first 11 on both sides
# Number of reads in each barcode
# Number thrown out
# GAPSIZE is the distance, in number of bases, between the paired-ends

use Getopt::Std;
use File::Basename;

# PARSE ARGUMENTS
my %parameters;
getopts( 'a:b:g:t:q:qc:lc:', \%parameters );    #Takes parameters

unless ( $parameters{a} ) {
    print "Usage: perl trim2.pl 
	-a raw illumina input file read 1
	-b raw illumina input file read 2 (if any) 
	-g size of GAP between paired-ends (if any) 
	-t truncate size (if any)
	-q quality file (in case of FASTA input)
	-qc quality cutoff value
	-lc minimum length \n";
    print "Supported formats: FASTA, FASTQ and QSEQ.\n";
    exit;
}

$read1 = $parameters{a};
unless ( open( READ1, "$read1" ) ) {
print "Error: Unable to open $read1.\n";
exit;
}


if ( $parameters{b} ) {
    print "Paired-ends.\n";
    $paired = 1;
    $read2 = $parameters{b};
	unless ( open( READ2, "$read2" ) ) {
    	print "Error: Unable to open $read2.\n";
    	exit;
	}
}



if  ( $parameters{g} ) {
$GAPSIZE = $parameters{g};
}

if  ( $parameters{t} ) {
$TRUNCATE = $parameters{t};
}

if  ( $parameters{qc} ) {
$QUALITY_CUTOFF = $parameters{qc};
}

if  ( $parameters{lc} ) {
$LENGTH_CUTOFF = $parameters{lc};
}


#Identify file format
open( FILETEST, "$read1" );
read FILETEST, $format, 1;
#print "$format\n";

$out_dir = dirname($read1);
$prefix = basename($read1);

# Make output directory, formated to run BLAST
system("mkdir -p output_files/trim2/");

# Make output file
open RUNBLAST, ">output_files/trim2/$prefix" . "_runblast.fasta" or die $!;

if ( $format eq ">" ){
    print "FASTA file format found.\n";
    #Read FASTA quality file
    if ( $parameters{q} ) {
	print "$parameters{q}\n";
	$readq = $parameters{q};
	unless ( open( READQ, "$readq" ) ) {
    	print "Error: Unable to open $readq required for FASTA file triming.\n";
    	exit;
	}
    #FASTA triming function
    parse_fasta();

    }else {
	print "Error: Please, specify the FASTA quality file with -q option.\n";
	exit;
     }
    
}elsif ( $format eq "@" ){
    print "FASTQ file format found.\n";
    #FASTQ triming function
    parse_fastq();
    
 }else {
    $format = <FILETEST>;
    chomp( $format );
    @format = split( "\t", $format );
    
    if ( ( ($format[7] eq "1") or ( $format[7] eq "2" ) ) and ( ( $format[10] eq "0" ) or ( $format[10] eq "1" ) ) ) {
	
	#print "@line1[7]\n@line1[10]\n";
	print "QSEQ file format found.\n";
	#QSEQ triming function
	parse_qseq();
    }else {
	print "Error: file format not recognized.\n";
     }  
  }

print "Trimming complete.\n";

sub parse_qseq {

@split_out1 = split( /\./, $prefix );



# Make Singletons Directory
system("mkdir -p $out_dir/singletons/");

open SINGLETONS, ">$out_dir/singletons/$prefix" . "_single.txt" or die $!;

while ( $read1_line = <READ1> ) {
    chomp( $read1_line );
    chomp( $read2_line = <READ2> );

    @line1 = split( /\t/, $read1_line );
    @line2 = split( /\t/, $read2_line );
    
    if ( (@line1[7] eq "1") ){
    
    $line1[8] =~ s/\./N/g;
    $line2[8] =~ s/\./N/g;

    $line1[8] = trim_qseq( $line1[8], $line1[9] );
    $line2[8] = trim_qseq( $line2[8], $line2[9] );

    $line1[9] = "";
    $line2[9] = "";
    }
    
    next if ( ( $line1[8] eq "0" ) && ( $line2[8] eq "0" ) );


    # Print Singletons
    if ( $line1[8] eq "0" ) {
        #print SINGLETONS join( "\t", @line2 ) . "\n";
        #print "@line1[8]\n @line1\n";

    }
    elsif ( $line2[8] eq "0" ) {
        #print SINGLETONS join( "\t", @line1 ) . "\n";
    }
    # Or, print entire thing (interleaved QSEQ)
    else {
       # Print in QSEQ format, uncomment (and comment FASTQ line)
       # print join( "\t", @line1 ) . "\n" . join( "\t", @line2 ) . "\n";

       # Print in FASTQ (interleaved) format:
       # TODO shorten this!
       $header = join(':', @line1[0..7]);
         
       $sequence_quality = $line1[8];
       @sq = split( /\t/, $sequence_quality);
       $sequence = $sq[0];
       $quality = $sq[1];
       #FASTQ format
       #print "\@$header:A\n$sequence\n\+$header:A\n$quality\n";
       
       # Print in FASTA format
       print RUNBLAST ">$header:AB\n$sequence";
       
       #Insert the GAP size between paired-ends, fill gap with Ns
       for ( $r = 0 ; $r < $GAPSIZE ; $r++){
		    print RUNBLAST "N";
       }
       
       $header = join(':', @line2[0..7]);
       $sequence_quality = $line2[8];
       @sq = split( /\t/, $sequence_quality);
       $sequence = $sq[0];
       $quality = $sq[1];
       # Print in FASTQ format
       #print "\@$header:B\n$sequence\n\+$header:B\n$quality\n";
       
       #FASTA format
       print RUNBLAST "$sequence\n";
    } 
}

#close SINGLE1;
#close SINGLE2;
}

# trim($seq, $qual);
sub trim_qseq() {
    $seq  = shift;
    $qual = shift;


    #Trim off primer
    $seq  = substr( $seq,  $TRUNCATE );
    $qual = substr( $qual, $TRUNCATE );

    substr($seq, 0, ($TRUNCATE-1)) = '';

    $max       = 0;
    $sum       = 0;
    $first     = 0;
    $start     = 0;
    $end       = 0;
    @qualities = split( //, $qual );
    $size      = scalar @qualities;

   for ( $a = 0 ; $a < $size ; $a++ ) {
        $sum += ( ord( $qualities[$a] ) - 64 - $QUALITY_CUTOFF );
        if ( $sum > $max ) {
            $max     = $sum;
            $end     = $a;
	    #$teste = (ord( $qualities[$a] ) - 64);
	    #print "$teste\n";
            #$start   = $first;
	    
        }
	#if the quality is under cutoff, replace by N
        if ( $sum < 0 ){
	    $sum = 0;
	    $seq[$a] = "N";
        }
    }

    $seq = substr( $seq, $start, $end - $start );
    
    #print "Seq = $seq\nStart = $start\nEnd = $end\n";
    
    return 0 if ( length($seq) < $LENGTH_CUTOFF );
    
    $seq = $seq . "\t" . substr( $qual, $start, $end - $start );
    
    return $seq;
}


sub parse_fasta {

$lengthCut = LENGTH_CUTOFF;
$qualityCut = QUALITY_CUTOFF;
    
$Rejected = -1;
@FinalTrim = ();

$end = 0;
$max = 0;
$sum = 0;
$start = 0;
$first = 0;
$lineNum = 0;

    while($lineSeq = <READ1>) {
    	
	$lineQual = <READQ>;

	if($lineSeq =~ />/) {	    
	    $length = ($end + 1) - $start;
		if( $length < $lengthCut ) {
		    $Rejected = 1;
		}
		if($Rejected == 0) {
		    print "$header\n";
		    $countDown = 60;
		    for($a = $start; $a <= $end; $a++) {
			$countDown--;
			print "$FinalTrim[$a]";
			if($countDown == 0) {
			    print "\n";
			    $countDown = 60;
			}
		    }
		    print "\n";
		}
		$Rejected = 0;
		
		#Trim off heading
		$space = index($lineSeq, " ") + 1;
		$header = substr($lineSeq, 0, $space);
		@FinalTrim = ();
		$max = 0;
		$sum = 0;
		$first = 0;
		$lineNum = 0;
	}
	elsif($Rejected == 0) {
		#Store for use when trimming bases
		@bases = split(//, $lineSeq);
		$size  = scalar @bases;
		$size--;
		for($a = 0; $a < $size; $a++) {
		    push(@FinalTrim, $bases[$a]);
		}
		
		chomp($lineQual);
		@line = split(/ /, $lineQual);
		$size = scalar @line;
		for($a = 0; $a < $size; $a++) {
		    $sum += ($line[$a] - $qualityCut);
		    if($sum > $max) {
			$max = $sum;
			$end = $a + (60*$lineNum);
			$start = $first;
		    }
		    if($sum < 0) {
			$sum = 0;
			$first = $a + (60*$lineNum);
		    }
		}
		$lineNum++;
	}
    }
    

close READ1;
close READQ;
#close OUTPUT;

}

sub parse_fastq {
    

    while ($header1 = <READ1>) {
    
 
	$sequence = <READ1>;
	$header2 = <READ1>;
	$quality = <READ1>;
   
	
    $fastq1 = trim_fastq( $sequence, $quality );
    
    chomp($header1);
    chomp($fastq1);
    
    $header1 =~ s/@//g;
    $fastq1 =~ s/\s//g;
    $fastq2 =~ s/\s//g;
    
    print RUNBLAST ">$header1:AB\n$fastq1";
    
    print ">$header1:AB\n$fastq1";
    
	if ($paired eq 1){
	$header1_2 = <READ1>;
	$sequence_2 = <READ1>;
	$header2_2 = <READ1>;
	$quality_2 = <READ1>;
   
    
	$fastq2 = trim_fastq( $sequence_2, $quality_2 );
	
	#print "$fastq1\n$fastq2\n";
	#Insert the GAP size between paired-ends, fill gap with Ns
	for ( $r = 0 ; $r < $GAPSIZE ; $r++){
	    print "N";
	    print RUNBLAST "N";
	 }
    
       #FASTA format
       chomp($fastq2);
       print "$fastq2\n";
       print RUNBLAST "$fastq2\n";
       
	}else{
	    
	    print "\n";
	    print RUNBLAST "\n";
	    
	    }
    }
    

	#print "$header1|$sequence|$header2|$quality\n";
	#my $solill_regexp = qr/[JKLMNOPQRSTUVWXYZ\[\]\^\_\`abcdefgh]/;
	
}
    

 
 
sub trim_fastq {
    
    
	    my $sequence = shift;
            my $quality = shift;
	    
	    $max       = 0;
	    $sum       = 0;
	    $first     = 0;
	    $start     = 0;
	    $end       = 0;
	    @qualities = split( //, $quality );
	    $size      = scalar @qualities;

	for ( $a = 0 ; $a < $size ; $a++ ) {
	    #$sum += 10 * log(1 + 10 ** (ord($qualities[$a]) - 64) / 10.0) / log(10);
	    
	    #$sum += ( ord( $qualities[$a] ) - 64 - $QUALITY_CUTOFF );
	    $sum += ( ord( $qualities[$a] ) - 33 - $QUALITY_CUTOFF );
	    
	    #print "$sum\n";
	     if ( $sum > $max ) {
	         $max     = $sum;
	         $end     = $a;
		    #$teste = (ord( $qualities[$a] ) - 64);
		    #print "$teste\n";
	         #$start   = $first;
	    
	     }
		#if the quality is under cutoff, replace by N
	     if ( $sum < 0 ){
		    $sum = 0;
		    $seq[$a] = "N";
	     }
	    }

	$seq = substr( $sequence, $start, $end - $start );
    
	#print "Seq = $seq\nStart = $start\nEnd = $end\n";
    
	return 0 if ( length($seq) < $LENGTH_CUTOFF );
    
	$seq = $seq . "\t" . substr( $qual, $start, $end - $start );
    
	return $seq;
	


    
}
