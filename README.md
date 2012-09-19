#PANGEA-plus

A new implementation of PANGEA pipeline for metagenomics with multiple classification methods and consensus analysis


#Extract the files:

            tar –xvf BioinfoTools_PANGEA-plus.tar.gz


Your work dir should be set as the PANGEA-plus directory.

            cd BioinfoTools_PANGEA-plus
            export $PANGEAWD=$PWD


#Trimming your input sequences

            cd $PANGEAWD/Trim
            perl trim2.3.pl -a ../input_A.txt -b ../input_B.txt -g 100

where: perl trim2.3.pl ...
	-a raw illumina input file read 1
	-b raw illumina input file read 2 (if any) 
	-g size of GAP between paired-ends (if any) 
	-t truncate size (if any)
	-q quality file (in case of FASTA input)
	-qc quality cutoff value
	-lc minimum length 

Supported formats: FASTA, FASTQ and QSEQ.

Results will be saved in $PANGEAWD/output/trim2 folder


#Download NCBI database for classification

            cd $PANGEAWD
            wget ftp://ftp.ncbi.nih.gov/blast/db/FASTA/nt.gz
            gunzip nt.gz

#Format the database

            $PANGEAWD/Classify/Runblast/formatdb -i $PANGEAWD/nt -p F

#Classify your sequences using parallel BLAST search

            cd $PANGEAWD/Classify/Runblast/Release

Example of parallel BLAST (MPI-blastn) executed in a PBS/Torque/Maui HPC cluster:

Submit interactive job:

            qsub –I -V

Create nodes list file:

            cat $PBS_NODEFILE > nodes.txt

Run the parallel BLAST:

            mpirun -np 4 -machinefile nodes.txt mpiblastn input.fasta database.formated $PANGEAWD/parallel_output.txt 12

where: 	input.fasta refers to your sequences after trimming.
database.formated is the name of database file. 
formatted by formatdb. 
12 refers to the total number of processes to be executed by MPI-blastn (nodes x cores).


Example using your own megablast installation:

            megablast -i input.fasta -d database/nt -m 8 -a 4 -o megablast_output.txt -D 3

#Parse the taxonomic classification based on the NCBI taxonomy databases

Running NCBI-taxcollector:

            cd $PANGEAWD/Tax_class
            make all
            wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz
            wget ftp://ftp.ncbi.nih.gov/pub/taxonomy/gi_taxid_nucl.dmp.tar.gz
            tar -xvf taxdump.tar.gz
            tar -xvf gi_taxid_nucl.dmp.tar.gz
            ./tax_class -c
            perl NCBI-taxcollector-0.01.pl –f $PANGEAWD/parallel_output.txt -o $PANGEAWD/parallel_output_class.txt > report.txt

where: 	parallel_output.txt is the mpiblastn results
parallel_output_class.txt are the parsed and classified output generated by this program.

#Cluster your results by similarity:

Example for 80% similarity*:

            perl $PANGEAWD/Megaclust/megaclust2.pl -i $PANGEAWD/parallel_output_class.txt -o $PANGEAWD/parallel_output_class.txt.megaclust_80_hits.txt -b 100 -s 80 -e 1e-20

*More examples and automatic scripts at $PANGEAWD/Scripts


#Generate summary table for classified results:

Example for Domain level (80%) similarity*:

            perl $PANGEAWD/Megaclustable/megaclustable.pl -m $PANGEAWD/parallel_output_class.txt.megaclust_80_hits.txt -t 0 -o $PANGEAWD /results/megaclustable/DomainTable.txt

*Note: in the –m option you shall list all the ouput files generated by the megaclust execution for every sample. More examples and automatic scripts at $PANGEAWD/Scripts.

The classification output should be like this:

            		1	2	3	4	5	6	7	8	9	10
            Bacteria	479	4	32	7507	11977	13245	2129	11222	539	2411	
            Eukaryota	1	4	5	5	2	17	78	3	10	3	
            Archaea		1	0	0	0	0	0	0	0	0	1		

