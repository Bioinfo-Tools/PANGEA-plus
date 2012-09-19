#!/usr/bin/perl

# $Id: cluster-blast-output.pl,v 1.3 2009/05/28 02:00:09 waltsbm Exp $

use strict;
use warnings;

use Getopt::Std;

######################################################################
# Options
#
# Required options:
# -i input BLAST tabular results file (megablast or blastall -m 8)
# -o output file name
#
# Optional parameters:
# -s similarity lower threshold (percent) (default 95)
# -e e-value upper threshold (default 1e-20)
# -b bitscore lower threshold (default 200)
# -d delimiter (default to comma)
#
# Optional switches:
# -c count every query hit (if -c not given, then only count
#                           any query-genome pair as one genome hit)
# -h print usage summary
######################################################################

#### Code starts

my %options;
getopts('i:o:s:e:b:c:d:h', \%options);

### Error checking and options setup
if ($options{h}) {
    print_usage();
    exit;
}

unless ($options{i} && $options{o}) {
    print "Must specify both an input and output filename\n";
    print_usage();
    exit;
};

my $sim_threshold = 95;
if ($options{s}) {
    if (($options{s} < 0) || ($options{s} > 100)) {
	print "similarity threshold must be between 0 and 100\n";
	print_usage();
	exit;
    } 

    $sim_threshold = $options{s};
}

my $eval_threshold = 1e-20 + 0;
if ($options{e}) {
    $eval_threshold = $options{e} + 0;
}

my $bitscore_threshold = 200;
if ($options{b}) {
    $bitscore_threshold = $options{b};
}

my $delimiter = ",";
if ($options{d}) {
    $delimiter = $options{d};
}

open INFILE, $options{i} or die "couldn't open infile";
open OUTFILE, ">$options{o}" or die "couldn't open outfile";

### Read BLAST results
my $lines_processed = 0;
my $lines_beyond_threshold = 0;
my %subject_counts;
my %queries_for_subject; #used to detect double hits;
while (<INFILE>) {
    next if (/^\#/);
    $lines_processed++;
    chomp;
    my ($query_id,
	$subject_id,
	$pct_identity,
	$aln_len,
	$mismatches,
	$gaps,
	$q_start,
	$q_end,
	$s_start,
	$s_end,
	$e_value,
	$bitscore,
	$temp ) = split /\t\t|\t\s|\s\t|\t/;
    
#    if ($pct_identity =~ /\D/){
#	$pct_identity = $aln_len;
#	$aln_len = $mismatches;
#	$mismatches = $gaps;
#	$gaps = $q_start;
#	$q_start = $q_end;
#	$q_end = $s_start;
#	$s_start = $s_end;
#	$s_end = $e_value;
#	$e_value = $bitscore;
#	$bitscore = $temp;
 #  }

    
# print "$query_id,
#	$subject_id,
#	$pct_identity,
#	$aln_len,
#	$mismatches,
#	$q_end,
#	$s_start,
#	$s_end,
#	$e_value,
#	$bitscore\n";
 
# print "evalue: $e_value, $s_end\n";
    
    #toss anything outside the three thresholds
    if (($pct_identity < $sim_threshold) ||
	($e_value > $eval_threshold) ||
	($bitscore < $bitscore_threshold)) {
	$lines_beyond_threshold++;
	next;
    }

    if ($options{c}) {
	#if we're couting every hit, it's easy, just increment the hitcount...
	$subject_counts{$subject_id}++;
    } else {
	#otherwise, need to check to see if we've already counted this query
	#seq hitting the subject
	unless (defined($queries_for_subject{$subject_id}{$query_id})) {
	    $subject_counts{$subject_id}++;
	}
	$queries_for_subject{$subject_id}{$query_id}++;
    }    
}

### create output
my @header = qw(OTU times_hit);
print OUTFILE join($delimiter, @header) . "\n";

foreach my $subject (keys(%subject_counts)) {
    my @outline = ($subject, $subject_counts{$subject});
    print OUTFILE join($delimiter, @outline) . "\n";
}

### print brief run summary to screen
print "Run complete:\n";
print "$lines_processed hits examined\n";
print "$lines_beyond_threshold hits beyond thresholds and therefore not counted.\n";


######
# End main, begin subs
######

sub print_usage {
    my $usage = qq|Usage:
		   cluster-blast-output.pl -i infile -o outfile [options]
		   
		   Required options:
		   -i input BLAST tabular results file (megablast or blastall -m 8)
		   -o output file name

		   Optional parameters:
		   -s similarity lower threshold (percent, between 0-100) (default 95)
		   -e e-value upper threshold (default 1e-20)
		   -b bitscore lower threshold (default 200)
		   -d delimiter (default to comma)
		   
		   Optional switches:
		   -c count every query hit (if -c not given, then only count
					     any query-genome pair as one genome hit)
		   -h print usage summary
		   |;

    print $usage . "\n";
		   
}
