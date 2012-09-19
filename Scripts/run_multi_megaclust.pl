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
# perl run_multi_megaclust.pl ../qseq_dir/ extension
# qseqdir: qseq files directory
# extension ("qseq" or "txt")
#################################################

use Cwd;


($dir, $extension) = @ARGV;

if (!$dir and !$extension) {
print   "
        run_multi_megaclust.pl Usage:

        perl run_multi_megaclust.pl ../blast_dir/ extension
        
        blast_dir: blast result files directory (default ./).
        extension: tabular file (default txt).\n\n ";
        
        exit;
}

$bitscore = 50;

my $wdir = getcwd;

$dir = "." unless $dir;

$lextension = "txt" unless $lextension;

chdir( $dir ) or die "Cannot chdir to $dir\n";
   
my @files = <*.$extension>;

closedir( DIR );

chdir( $wdir );
   
my $count = @files;
   
#print "@files";

chomp( @files );

$comands = "";
  
for ($i=0; $i < $count; $i++) {
    
   $comands = "$comands" . "perl ./Megaclust/megaclust2.pl -i $dir$files[$i] -o ./results/megaclust/$files[$i].megaclust_80_hits.txt -b $bitscore -s 80 -e 1e-20\n";
   $comands = "$comands" . "perl ./Megaclust/megaclust2.pl -i $dir$files[$i] -o ./results/megaclust/$files[$i].megaclust_90_hits.txt -b $bitscore -s 90 -e 1e-20\n";
   $comands = "$comands" . "perl ./Megaclust/megaclust2.pl -i $dir$files[$i] -o ./results/megaclust/$files[$i].megaclust_95_hits.txt -b $bitscore -s 95 -e 1e-20\n";
   $comands = "$comands" . "perl ./Megaclust/megaclust2.pl -i $dir$files[$i] -o ./results/megaclust/$files[$i].megaclust_98_hits.txt -b $bitscore -s 98 -e 1e-20\n";
   $comands = "$comands" . "perl ./Megaclust/megaclust2.pl -i $dir$files[$i] -o ./results/megaclust/$files[$i].megaclust_99_hits.txt -b $bitscore -s 99 -e 1e-20\n";
    
}

print "$comands";

my $script = `$comands`;

print $script;

#system( "$script" );

#exec( "$script" );



exit;
 
