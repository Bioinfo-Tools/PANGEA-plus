/*
 * 31mer/pregraph.c
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
#include "extvab.h"

static void initenv ( int argc, char ** argv );
static char shortrdsfile[256];
static char graphfile[256];
static int cutTips = 1;

static void display_pregraph_usage();

int call_pregraph ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	initenv ( argc, argv );

	if ( overlaplen % 2 == 0 )
	{
		overlaplen++;
		printf ( "K should be an odd number\n" );
	}

	if ( overlaplen < 13 )
	{
		overlaplen = 13;
		printf ( "K should not be less than 13\n" );
	}
	else if ( overlaplen > 31 )
	{
		overlaplen = 31;
		printf ( "K should not be greater than 31\n" );
	}

	time ( &time_bef );
	prlRead2HashTable ( shortrdsfile, graphfile );
	time ( &time_aft );
	printf ( "time spent on pre-graph construction: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	printf ( "deLowKmer %d, deLowEdge %d\n", deLowKmer, deLowEdge );
	//analyzeTips(hash_table, graphfile);

	if ( !deLowKmer && cutTips )
	{
		time ( &time_bef );
		removeSingleTips();
		removeMinorTips();
		time ( &time_aft );
		printf ( "time spent on cutTipe: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	}
	else
	{
		time ( &time_bef );
		removeMinorTips();
		time ( &time_aft );
		printf ( "time spent on cutTipe: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	}

	initKmerSetSize = 0;
	//combine each linear part to an edge
	time ( &time_bef );
	kmer2edges ( graphfile );
	time ( &time_aft );
	printf ( "time spent on making edges: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	//map read to edge one by one
	time ( &time_bef );
	prlRead2edge ( shortrdsfile, graphfile );
	time ( &time_aft );
	printf ( "time spent on mapping reads: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	output_vertex ( graphfile );
	free_Sets ( KmerSets, thrd_num );
	free_Sets ( KmerSetsPatch, thrd_num );
	time ( &stop_t );
	printf ( "overall time for lightgraph: %dm\n\n", ( int ) ( stop_t - start_t ) / 60 );
	return 0;
}

/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/


void initenv ( int argc, char ** argv )
{
	int copt;
	int inpseq, outseq;
	extern char * optarg;
	char temp[100];
	optind = 1;
	inpseq = outseq = 0;

	while ( ( copt = getopt ( argc, argv, "a:s:o:K:p:d:DR" ) ) != EOF )
	{
		//printf("get option\n");
		switch ( copt )
		{
			case 's':
				inpseq = 1;
				sscanf ( optarg, "%s", shortrdsfile );
				continue;
			case 'o':
				outseq = 1;
				sscanf ( optarg, "%s", graphfile ); //
				continue;
			case 'K':
				sscanf ( optarg, "%s", temp ); //
				overlaplen = atoi ( temp );
				continue;
			case 'p':
				sscanf ( optarg, "%s", temp ); //
				thrd_num = atoi ( temp );
				continue;
			case 'R':
				repsTie = 1;
				continue;
			case 'd':
				sscanf ( optarg, "%s", temp ); //
				deLowKmer = atoi ( temp ) >= 0 ? atoi ( temp ) : 0;
				continue;
			case 'D':
				deLowEdge = 1;
				continue;
			case 'a':
				initKmerSetSize = atoi ( optarg );
				break;
			default:

				if ( inpseq == 0 || outseq == 0 )       //
				{
					display_pregraph_usage();
					exit ( -1 );
				}
		}
	}

	if ( inpseq == 0 || outseq == 0 )  //
	{
		//printf("need more\n");
		display_pregraph_usage();
		exit ( -1 );
	}
}

static void display_pregraph_usage()
{
	printf ( "\npregraph -s readsInfoFile [-d KmerFreqCutoff -R -K kmer -p n_cpu] -o OutputFile\n" );
	printf ( "  -s readsInfoFile: The file contains information of solexa reads\n" );
	printf ( "  -a initMemoryAssumption: Initiate the memory assumption to avoid further reallocation\n" );
	printf ( "  -p n_cpu(default 8): number of cpu for use\n" );
	printf ( "  -K kmer(default 21): k value in kmer\n" );
	printf ( "  -d KmerFreqCutoff(optional): delete kmers with frequency no larger than (default 0)\n" );
	printf ( "  -R (optional): unsolve repeats by reads (default no)\n" );
	printf ( "  -o OutputFile: prefix of output file name\n" );
}

