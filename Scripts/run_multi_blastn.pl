# !/usr/bin/perl -w
#################################################
# run_multi_blastn.pl
# by Raquel Dias
# Laboratorio de Alto Desempenho (LAD) - PUCRS
# e-mail: raquel.dias001@acad.pucrs.br
#################################################
# This scripts runs trimming on multiple fasta
# files inside a same directory
# Usage:
# perl run_multi_blastn.pl ../fasta_dir/ extension
# fastadir: fasta files directory
# extension ("fasta" or "txt")
#################################################

use Cwd;


($dir, $database) = @ARGV;

if (!$dir or !$database) {
print   "
        run_multi_blastn.pl Usage:

        perl run_multi_blastn.pl ../fasta_dir/ database        
        fastadir: input fasta files directory (default ./).
        database: path and name to the database file for running Blastn. Example: \$PANGEAWD/database/my_database_name\n\n ";
        
        exit;
}

my $wdir = getcwd;

$dir = "." unless $dir;

$extension = "fasta";

chdir( $dir ) or die "Cannot chdir to $dir\n";
   
my @files = <*.$extension>;

closedir( DIR );

chdir( $wdir );
   
my $count = @files;
   
#print "@files";

chomp( @files );

$comands = "";
   
for ($i=0; $i < $count; $i++) {
    
   $comands = "$comands" . "./Classify/RunBlast/blastn -query $dir$files[$i] -db $database -outfmt 6 -out blastn_output_$files[$i].txt\n";
   
    
}

print "$comands";

my $script = `$comands`;

system( "$script" );

exec( "$script" );

print $script;

exit;
 
