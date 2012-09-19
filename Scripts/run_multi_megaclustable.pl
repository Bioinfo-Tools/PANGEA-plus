# !/usr/bin/perl -w
#################################################
# run_multi_trim_qseq.pl
# by Raquel Dias
# Laboratorio de Alto Desempenho (LAD) - PUCRS
# e-mail: raquel.dias001@acad.pucrs.br
#################################################
# This scripts runs trimming on multiple QSEQ
# files inside a same directory
# Usage:
# perl run_multi_megaclustable.pl ../megaclust_dir/ extension
# megaclust_dir: megaclust result files directory (default ./).
# extension: tabular file (default txt).
#################################################

use Cwd;


($dir, $extension) = @ARGV;

if (!$dir and !$extension) {
print   "
        run_multi_megaclust.pl Usage:

        perl run_multi_megaclustable.pl ../megaclust_dir/ extension
        
        megaclust_dir: megaclust result files directory (default ./).
        extension: tabular file (default txt).\n\n ";
        
        exit;
}

open FILE80, ">$dir/80.txt" or die $!;
open FILE90, ">$dir/90.txt" or die $!;
open FILE95, ">$dir/95.txt" or die $!;
open FILE99, ">$dir/99.txt" or die $!;

my $wdir = getcwd;

$dir = "." unless $dir;

$extension = "txt" unless $extension;

#chdir( $dir ) or die "Cannot chdir to $dir\n";
   
#@files80 = <*80_hits.$extension>;
opendir(DIR, "$dir");
@list = readdir(DIR);

@files80 = grep( /80_hits\.$extension$/ , @list );
@files90 = grep( /90_hits\.$extension$/ , @list );
@files95 = grep( /95_hits\.$extension$/ , @list );
@files99 = grep( /99_hits\.$extension$/ , @list );

closedir( DIR );

chomp(@files80);
chomp(@files90);
chomp(@files95);
chomp(@files99);

@files80 = sort(@files80);
@files90 = sort(@files90);
@files95 = sort(@files95);
@files99 = sort(@files99);


$files80f = join( " $dir", @files80 );
$files90f = join( " $dir", @files90 );
$files95f = join( " $dir", @files95 );
$files99f = join( " $dir", @files99 );


print FILE80 "$files80f";
print FILE90 "$files90f";
print FILE95 "$files95f";
print FILE99 "$files99f";


$comands = "";

$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files80f -t 0 -o ./results/megaclustable/DomainTable.txt\n";
$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files80f -t 1 -o ./results/megaclustable/PhylumTable.txt\n";
$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files90f -t 2 -o ./results/megaclustable/ClassTable.txt\n";
$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files90f -t 3 -o ./results/megaclustable/OrderTable.txt\n";
$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files90f -t 4 -o ./results/megaclustable/FamilyTable.txt\n";
$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files95f -t 5 -o ./results/megaclustable/GenusTable.txt\n";
$comands = "$comands" . "perl ./Megaclustable/megaclustable.pl -m $dir$files99f -t 6 -o ./results/megaclustable/SpeciesTable.txt\n";

print "$comands";

my $script = `$comands`;

print $script;

system( "$script" );

exit;
 
