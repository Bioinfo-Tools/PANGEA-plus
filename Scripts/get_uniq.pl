# !/usr/bin/perl -w


use Getopt::Std;

use File::Basename;

my %parameters;
getopts( 'f:', \%parameters );    #Takes parameters

unless ( $parameters{f} ) {
    print "Usage: perl taxcollector_ncbi-0.01.pl 
	-f Classification results (tabular text file)\n";
    exit;
}

print "\nLoading input file...\n";

###############Check names##########################
$readnameid = $parameters{f};
unless ( open( READ, "$readnameid" ) ) {
print "Error: Unable to open database file $readnameid.\n";
exit;
}

$output = "$readnameid."."unique";

unless ( open( OUT, ">$output" ) ) {
print "Error: Unable to open output file $output.\n";
exit;
}


my %seen = ();

while ($line = <READ>) {
    chomp;
    my @columns = split ('\t', $line);
    print OUT "$line" if ! $seen{$columns[0]}++;
}