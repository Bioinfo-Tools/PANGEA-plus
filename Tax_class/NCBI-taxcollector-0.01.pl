# !/usr/bin/perl -w

############################################################
#
#NCBI-taxcollector
#Attaches taxonomic information to BLAST and SOAP2 results
#Results must be in tabular format (-m 8 or -m 6)
#By: Raquel Dias and Marcelo V. Neves
#Contact: raquel.dias.001@acad.pucrs.br
#Pontifical Catholic University of Rio Grande do Sul - Brazil
#High Performance Laboratory (LAD PUCRS)
#
############################################################

use Getopt::Std;
use File::Basename;


my %parameters;
getopts( 'f:o:', \%parameters );    #Takes parameters

unless ( $parameters{f} && $parameters{o} ) {
    print "Usage: perl taxcollector_ncbi-0.01.pl 
	-f Classification results (tabular text file)
	-o Output file \n";
    exit;
}

################Check INPUT FILE########################
$read1 = $parameters{f};
unless ( open( READ1, "$read1" ) ) {
print "Error: Unable to open classification results file $read1.\n";
exit;
}
###########################################################



###################Create OUTPUT###########################
$out = $parameters{o};
unless ( open( OUT, ">$out" ) ) {
print "Error: Unable to open output file $out.\n";
exit;
}
###########################################################

chdir("./Tax_class/");


#START
read_fasta();

close READ1;
close OUT;



sub read_fasta {
        my $line;

        $seqnum = 0;
	
        while ( $line = <READ1> ) {
                
                START:

                chomp( $line );

                $seqheader[$seqnum] = $line;
                          
                #Original header
	        ####print "$seqheader[$seqnum]\n";                                        
 
                #Gets the seuqnce GI ID
                @gi = split( /\|/, $seqheader[$seqnum] );
		 
                @id = split( /\ |\t\t|\t/, $seqheader[$seqnum] );
                #Gets Taxonomic information

		#print "\n...Parsing taxonomy for GI: @gi[1]...\n\n";
	
		#Convert the TaxID to complete taxonomic information
		if ( @gi ) {
                parse_taxonomy( @gi[1] );
		}else{
		    exit;
		}
		#print "\nArray = @ClassArray\n\n";

		#Print sample sequence name
		####print "$id[0]\t"; 	
		print OUT "$id[0]\t"; 	

		
		#Prepare results array to sort results
		$ClassArray2 = join( /./, @ClassArray);
		@ClassArray2 = split( /\|/, $ClassArray2);
		#print "Splited array = $ClassArray2\n\n";
		
		#Get array size	
		$size = @ClassArray2;
		#print "Size = $size\n\n";
		
		#Sort taxonomic header results
		for ($i = $size+1; $i>=0 ; $i--) {
		    #print "@ClassArray2[$i]";
		    
		    #If finds the Species rank level
		    if ( grep /[6]/, @ClassArray2[$i] ){
			
			#Replaces the spaces by '_'
			@ClassArray2[$i] =~ s/\s/_/g;
			#Prints the result
			
			
			if ( !grep /[5]/, @ClassArray2 ){
			
			#If can't find genus, reprint the species at genus rank level
			@ClassArray2[$i] =~ s/6/5/;
			print OUT "@ClassArray2[$i]";
			@ClassArray2[$i] =~ s/5/6/;
			print OUT "@ClassArray2[$i]";
			
			}else{
			    print OUT "@ClassArray2[$i]";
			}
			
		    }else{
		    
		    if ( grep /[7]/, @ClassArray2[$i] ){
			
			@ClassArray2[$i] =~ s/7/9/;
			#In any other rank level cases
			print OUT "@ClassArray2[$i]";
		    }else{
			print OUT "@ClassArray2[$i]";
		    }
		    

		    }
		    
		    
		    
		}

		#Print the rest of blast results
		####print "\t@gi[4]\n";
		for ($i = 2; $i <= 12; $i++){
		    
		    if ( $id[$i] ne "" ){
		    print OUT "\t$id[$i]";
		    }
		}

		print OUT "\n";
		
		$ClassArray2 = "";
		@ClassArray = "";
                		    		
                $seqnum = $seqnum + 1;
                        

        }
}
        
sub parse_taxonomy {
        
    $gi = shift;
        
    my $taxline =  `./tax_class -s $gi`;
		    
    $taxline =~ s/\t//g;
    $taxline =~ s/\ //g;
    @line = split( /\|/, $taxline );
               
    print "Searching upper node for TAXID $line[0].\n";
    if ( $line[0] eq "0\n") {
	print "\n\nTAXID zero GI = $gi.\n\n";
	push ( @ClassArray, "Unidentified(GI:$gi);|" );
	return;
	}else{
        get_uptaxa( $line[0] );
	print "Done for TAXID $line[0].\n";
    }
    return;
}

sub get_name {

    $taxid = shift;

    my $namearray =  `./tax_class -n $taxid`;

    @namearray = split( /\n/, $namearray);

#print "Name array: $namearray\n\n";

    foreach $nameline ( @namearray ) {
	
	chomp($nameline);
	$nameline =~ s/\t//g;
	@nameline = split( /\|/, $nameline );
	
	#print "Nameline: $nameline[3]\n\n";
	        
	#chomp( $taxid );
      
	    if( ( grep /scientific name/, $nameline[3] ) ) {
                # or ( grep /synonym/, $nameline[3] )
	    #print taxonomic rank name
	    $nameline[1] =~ s/\s+$//; #remove spaces from end
	    $nameline[1] =~ s/^\s+//; #remove spaces from begin
	    #$nameline[1] =~s/\ //g;
	    $nameline[1] =~ s/\t//g; #remove tabs
	    
	    #print "Name: $nameline[1];\n";
	    push ( @ClassArray, "$nameline[1];|" );
	    
	    return;
	
	    }
    }

}

sub get_uptaxa {
 
    @ranklist = (
                 "superkingdom",
                 "phylum",
                 "class",
                 "order",
                 "family",
                 "genus",
                 "species",
		 "kingdom"
                 );
        
    local ( $taxid ) = shift; 
    
    my $nodeline =  `./tax_class -t $taxid`;	
    
    chomp($nodeline);
    
    #print "$nodeline\n";
    
    if(!$nodeline){print "NODE ID NOTFOUND"; return;}	
    
        
    $nodeline =~ s/\t//g;
    
    $nodeline =~ s/\ //g;
        
    @nodeline = split( /\|/, $nodeline );

    

        
    chomp( $taxid );
        
            if ( grep /^$nodeline[2]/, @ranklist ) {
                
                #get taxonomic rank number
                $ind = indexArray($nodeline[2], @ranklist);
                
                push (@ClassArray, "[$ind]");
                
		#print "[$ind]\n";
		
                #print "Type: $nodeline[2]\n";

                #Get name from taxon ID
                get_name( $taxid );
        
                        if ( $nodeline[2] eq $ranklist[0] ) {
                
                        #print "Found $nodeline[2]\n";
                        return;
                
                        }else{
                                #If still not found upmost rank, look for upper taxon ID
                                get_uptaxa( $nodeline[1] );
                        }
            }else{
		if ( grep /^$nodeline[2]/, "norank" ){
		
		#get_name( $taxid );
		
		    if($nodeline[1] eq "1"){
			push (@ClassArray, "[0]Unclassified;|");
		    }else{
			get_uptaxa( $nodeline[1] );
		    }
		
		return;
		}		
	        get_uptaxa( $nodeline[1] );
	     }
          
}


sub indexArray {

 local ($string) = shift;       
 local (@array) = @_;
 
 $i = 0;
 
        foreach $item (@array) {
                if ($string eq $item) {
                        return $i; 
                }else{
                        
                     #print "|$item|$string|\n";   
                }
                $i++;
        }
}

