# !/usr/bin/perl -w

use Getopt::Std;

use File::Basename;

my %parameters;
getopts( 'b:r:s:o:', \%parameters );    #Takes parameters

unless ( $parameters{b} && $parameters{r} && $parameters{o}) {
    print "Usage: perl Consensus-1.0.pl 
	-b Classification results (Blast)
	-r Classification results (RDP)
	-s Classification results (SOAP2)
	-o Output file (txt)\n";
    exit;
}

print "\nLoading input files...\n";



###############Check files##########################
$readnameid1 = $parameters{b};
unless ( open( READ1, "$readnameid1" ) ) {
print "Error: Unable to open $readnameid1 file.\n";
exit;
}

###############Check files##########################
$readnameid2 = $parameters{r};
unless ( open( READ2, "$readnameid2" ) ) {
print "Error: Unable to open $readnameid2 file.\n";
exit;
}


###############Check files##########################

if ( $parameters{s} ) {
    $readnameid3 = $parameters{s};
    unless ( open( READ3, "$readnameid3" ) ) {
    print "Error: Unable to open $readnameid3 file.\n";
    exit;
    }
}



$output = $parameters{o};
print "$parameters{o}\n";
unless ( open( OUT, ">$output" ) ) {
print "Error: Unable to open output file $output.\n";
exit;
}


$count = 0;

$linecount = 0;

$maxrank = 7;

@rdpranks = (
	     "domain",
	     "phylum",
	     "class",
	     "order",
	     "family",
	     "genus",
	     "species"
	     );

@blastranks = (
	       "0",
	       "1",
	       "2",
	       "3",
	       "4",
	       "5",
	       "6"
	       );

$maxlastcount = 0;

$maxrankmatches = 0;

@read1 = <READ1>;

$blastelements = scalar(@read1);

$i = 0;

while( $rdpline = <READ2> ) {
    
    chomp($rdpline);
    $i = 0;
    
    GETBLAST: $blastline = $read1[$i];
    
    $i++;
    
    chomp($blastline);
    
    @blastline = split( /\t/, $blastline);
   
    
    $blasttax = $blastline[1];
    
    chomp($blasttax);
    
    @blasttax =  split( /\[|\]|\;/, $blasttax);

    #print "$blastline[1]\n";
    
    $blasttax = join (' ', @blasttax);
    
    @blasttax = split (' ', $blasttax);
    
    #print "$blasttax[0]\n";
    
    @rdpline = split( /\t\t\t\t\t/, $rdpline);
    
    #print "$rdpline[1]\n";
    
    $rdptax = $rdpline[1];
    
    @rdptax = split( /\t/, $rdptax);
    
   chomp($blastline[0]);
   chomp($rdpline[0]);

    #print "$blastline[0]\t $rdpline[0]\n";
    
    
	# If blast and rdp query sequences are the same
	if ( $blastline[0] eq $rdpline[0] ) {
    #print "$blastline[0]\t $rdpline[0]\n";
	
	#found the same query 
	$found = 1;
		
	#counts the rank matches among classifications
	$rankmatches = 0;
		
	$a = 0;
	$b = 0;
	
	    #compare results of Blast and RDP rank by rank of all taxonomic levels, if any
	    for ( $a = 0; $a < @blasttax; $a = $a+2 ) {
			    
		for ($b = 0; $b < @rdptax; $b = $b+3 ) {
			
			
		    $rdptax[$b] =~ s/"|\\//g;
		    $rdptax[$b] =~ s/[\W\d_]//g;

		    #print "comparing $blasttax[$a] and $rdptax[$b+1] = $blasttax[$a+1] and $rdptax[$b]\n";
		    
		    my %index1;
		    @index1{@blastranks} = (0..$#blastranks);
		    my $index1 = $index1{$blasttax[$a]};
		    
		    my %index2;
		    @index2{@rdpranks} = (0..$#rdpranks);
		    my $index2 = $index2{$rdptax[$b+1]};
		    
		    #print "blasttax=$blasttax[$a+1] ---- rdptax=$rdptax[$b] \n indx1=$index1 ----- indx2=$index2\n";

		    if ( $blasttax[$a+1] eq $rdptax[$b] && $index1 eq $index2 ) {
			    #|| grep /^$blasttax[$a+1]/, $rdptax[$b]
			    #print "Match: $blasttax[$a+1] and $rdptax[$b]\n";
		    	    $rankmatches++;
						
		    }
		    
		}
	    
	    
	    }
	    
	    $blastcount = @blasttax;
	    #print "matches = $rankmatches of $maxrankmatches\n";
	    #print "blastcount = $blastcount of $maxblastcount\n";
	    
	    #counts if rank match found is bigger than the stored one
	    if ($rankmatches gt $maxrankmatches) {
		
		$maxrankmatches = $rankmatches;
		$tempresult = $blastline;
		}
	    
	    #in case of draw, counts if the blast result is more complete than the stored one, and if the similarity is higher
	    if ( ($blastcount > $maxblastcount or $blastsim < $blastline[2]) and $rankmatches eq $maxrankmatches){
	    $maxblastcount = $blastcount;
	    $tempresult = $blastline;
	    $blastsim = $blastline[2];
		
	    }
		
		$rankmatches = 0;
		
	#print "$blastelements\n";
    
	#if($i lt $blastelements){
	goto GETBLAST;
	#}
	
	}else{
	
	    if( $found eq 0 ){
		goto GETBLAST;
		$rankmatches = 0;

	     }
	    
	    if( $found eq 1 ){
		
		print OUT "$tempresult\n";
		print OUT "#Matches found: $maxrankmatches\n";
		$found = 0;
		$maxblastcount = 0;
		$maxrankmatches = 0;
		$rankmatches = 0;
		$blastsim = 0;
		
	    
	    }
	    
	     
	}
    
    
    }
    
    print "\nDone!\n";
    
    close (OUT);
