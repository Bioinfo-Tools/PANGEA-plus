#!/usr/bin/perl -w
#################################################
# run_multi_trim_qseq.pl
# by Raquel Dias
# Laboratorio de Alto Desempenho (LAD) - PUCRS
# e-mail: raquel.dias001@acad.pucrs.br
#################################################
# This scripts runs trimming on multiple QSEQ
# files inside a same directory
# Usage:
# perl run_multi_trim_qseq.pl ../qseq_dir/ extension
# qseqdir: qseq files directory
# extension ("qseq" or "txt")
#################################################

use Cwd;


($dir, $extension) = @ARGV;

if (!$dir and !$extension) {
print   "
        run_multi_join_fasta.pl Usage:

        perl run_multi_join_fasta.pl ./fasta_dir/ extension
        
        fasta_dir: fasta files directory (default ./).
        extension: fasta or txt (default txt).\n\n ";
        
        exit;
}

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
    
   $comands = "$comands" . "perl ./Trim/trim2.4.pl -a $dir$files[$i] ";
   
   $i++;
   
   $comands = "$comands" . "-b $dir$files[$i] -j \n";
    
}

print "$comands";

my $script = `$comands`;

system( "$script" );

exec( "$script" );

print $script;

exit;

