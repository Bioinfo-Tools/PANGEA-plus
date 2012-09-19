/*
 * 127mer/main.c
 * 
 * Copyright (c) 2008-2010 BGI-Shenzhen <soap at genomics dot org dot cn>. 
 *
 * This file is part of SOAPdenovo.
 *
 * SOAPdenovo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * SOAPdenovo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with SOAPdenovo.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include "stdinc.h"
#include "newhash.h"
#include "extfunc.h"
#include "global.h"

extern int call_pregraph ( int arc, char ** argv );
extern int call_heavygraph ( int arc, char ** argv );
extern int call_map2contig ( int arc, char ** argv );
extern int call_scaffold ( int arc, char ** argv );
extern int call_align ( int arc, char ** argv );

static void display_usage();
static void display_all_usage();
static void pipeline ( int argc, char ** argv );

int main ( int argc, char ** argv )
{
	printf ( "\nVersion 1.05: testing... 2010\n\n" );
	argc--;
	argv++;

	/*
	    printf("Size of Kmer is %d\n",sizeof(Kmer));
	    printf("Size of kmer_t is %d\n",sizeof(kmer_t));
	*/
	if ( argc == 0 )
	{
		display_usage();
		return 0;
	}

	if ( strcmp ( "pregraph", argv[0] ) == 0 )
		{ call_pregraph ( argc, argv ); }
	else if ( strcmp ( "contig", argv[0] ) == 0 )
		{ call_heavygraph ( argc, argv ); }
	else if ( strcmp ( "map", argv[0] ) == 0 )
		{ call_align ( argc, argv ); }
	//call_map2contig(argc,argv);
	else if ( strcmp ( "scaff", argv[0] ) == 0 )
		{ call_scaffold ( argc, argv ); }
	else if ( strcmp ( "all", argv[0] ) == 0 )
		{ pipeline ( argc, argv ); }
	else
		{ display_usage(); }

	return 0;
}

static void display_usage()
{
	printf ( "\nUsage: SOAPdenovo <command> [option]\n" );
	printf ( "    pregraph     construction kmer-graph\n" );
	printf ( "    contig       eliminate errors and output contigs\n" );
	printf ( "    map          map reads to contigs\n" );
	printf ( "    scaff        scaffolding\n" );
	printf ( "    all          doing all the above in turn\n" );
}

static void pipeline ( int argc, char ** argv )
{
	char * options[16];
	unsigned char getK, getRfile, getOfile, getD, getDD, getL, getR, getP;
	char readfile[256], outfile[256];
	char temp[128];
	char * name;
	int kmer = 0, cutoff_len = 0, ncpu = 0;
	char kmer_s[16], len_s[16], ncpu_s[16], M_s[16];
	int i, copt, index, M = 1;
	extern char * optarg;
	time_t start_t, stop_t;
	time ( &start_t );
	getK = getRfile = getOfile = getD = getDD = getL = getR = getP = 0;

	while ( ( copt = getopt ( argc, argv, "s:o:K:M:L:p:G:dDRua:" ) ) != EOF )
	{
		switch ( copt )
		{
			case 's':
				getRfile = 1;
				sscanf ( optarg, "%s", readfile );
				continue;
			case 'o':
				getOfile = 1;
				sscanf ( optarg, "%s", outfile ); //
				continue;
			case 'K':
				getK = 1;
				sscanf ( optarg, "%s", temp ); //
				kmer = atoi ( temp );
				continue;
			case 'G':
				sscanf ( optarg, "%s", temp ); //
				GLDiff = atoi ( temp );
				continue;
			case 'M':
				sscanf ( optarg, "%s", temp ); //
				M = atoi ( temp );
				continue;
			case 'p':
				getP = 1;
				sscanf ( optarg, "%s", temp ); //
				ncpu = atoi ( temp );
				continue;
			case 'L':
				getL = 1;
				sscanf ( optarg, "%s", temp ); //
				cutoff_len = atoi ( temp );
				continue;
			case 'R':
				getR = 1;
				continue;
			case 'u':
				maskRep = 0;
				continue;
			case 'd':
				getD = 1;
				continue;
			case 'D':
				getDD = 1;
				continue;
			case 'a':
				initKmerSetSize = atoi ( optarg );
				break;
			default:

				if ( getRfile == 0 || getOfile == 0 )       //
				{
					display_all_usage();
					exit ( -1 );
				}
		}
	}

	if ( getRfile == 0 || getOfile == 0 )       //
	{
		display_all_usage();
		exit ( -1 );
	}

	if ( thrd_num < 1 )
		{ thrd_num = 1; }

	// getK = getRfile = getOfile = getD = getL = getR = 0;
	name = "pregraph";
	index = 0;
	options[index++] = name;
	options[index++] = "-s";
	options[index++] = readfile;

	if ( getK )
	{
		options[index++] = "-K";
		sprintf ( kmer_s, "%d", kmer );
		options[index++] = kmer_s;
	}

	if ( getP )
	{
		options[index++] = "-p";
		sprintf ( ncpu_s, "%d", ncpu );
		options[index++] = ncpu_s;
	}

	if ( getD )
		{ options[index++] = "-d"; }

	if ( getDD )
		{ options[index++] = "-D"; }

	if ( getR )
		{ options[index++] = "-R"; }

	options[index++] = "-o";
	options[index++] = outfile;

	for ( i = 0; i < index; i++ )
		{ printf ( "%s ", options[i] ); }

	printf ( "\n" );
	call_pregraph ( index, options );
	name = "contig";
	index = 0;
	options[index++] = name;
	options[index++] = "-g";
	options[index++] = outfile;
	options[index++] = "-M";
	sprintf ( M_s, "%d", M );
	options[index++] = M_s;

	if ( getR )
		{ options[index++] = "-R"; }

	for ( i = 0; i < index; i++ )
		{ printf ( "%s ", options[i] ); }

	printf ( "\n" );
	call_heavygraph ( index, options );
	name = "map";
	index = 0;
	options[index++] = name;
	options[index++] = "-s";
	options[index++] = readfile;
	options[index++] = "-g";
	options[index++] = outfile;

	if ( getP )
	{
		options[index++] = "-p";
		sprintf ( ncpu_s, "%d", ncpu );
		options[index++] = ncpu_s;
	}

	for ( i = 0; i < index; i++ )
		{ printf ( "%s ", options[i] ); }

	printf ( "\n" );
	call_align ( index, options );
	name = "scaff";
	index = 0;
	options[index++] = name;
	options[index++] = "-g";
	options[index++] = outfile;
	options[index++] = "-F";

	if ( getL )
	{
		options[index++] = "-L";
		sprintf ( len_s, "%d", cutoff_len );
		options[index++] = len_s;
	}

	for ( i = 0; i < index; i++ )
		{ printf ( "%s ", options[i] ); }

	printf ( "\n" );
	call_scaffold ( index, options );
	time ( &stop_t );
	printf ( "time for the whole pipeline: %dm\n", ( int ) ( stop_t - start_t ) / 60 );
}

static void display_all_usage()
{
	printf ( "\nSOAPdenovo all -s configFile [-K kmer -d -D -M mergeLevel -R -u -G gapLenDiff -L minContigLen -p n_cpu] -o Output\n" );
	printf ( "  -s ShortSeqFile: The input file name of solexa reads\n" );
	printf ( "  -K kmer(default 23): k value in kmer\n" );
	printf ( "  -p n_cpu(default 8): number of cpu for use\n" );
	printf ( "  -a initKmerSetSize: define the initial KmerSet size(unit: GB)\n" );
	printf ( "  -M mergeLevel(default 1,min 0, max 3): the strength of merging similar sequences during contiging\n" );
	printf ( "  -d (optional): delete kmers with frequency one (default no)\n" );
	printf ( "  -D (optional): delete edges with coverage one (default no)\n" );
	printf ( "  -R (optional): unsolve repeats by reads (default no)\n" );
	printf ( "  -G gapLenDiff(default 50): allowed length difference between estimated and filled gap\n" );
	printf ( "  -L minLen(default K+2): shortest contig for scaffolding\n" );
	printf ( "  -u (optional): un-mask contigs with high coverage before scaffolding (default mask)\n" );
	printf ( "  -o Output: prefix of output file name\n" );
}
