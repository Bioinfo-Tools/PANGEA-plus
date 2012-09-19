/*
 * 31mer/scaffold.c
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
static void display_scaff_usage();

static boolean LINK, SCAFF;
static char graphfile[256];

int call_scaffold ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	initenv ( argc, argv );
	loadPEgrads ( graphfile );
	time ( &time_bef );
	loadUpdatedEdges ( graphfile );
	time ( &time_aft );
	printf ( "time spent on loading edges %ds\n", ( int ) ( time_aft - time_bef ) );

	if ( !SCAFF )
	{
		time ( &time_bef );
		PE2Links ( graphfile );
		time ( &time_aft );
		printf ( "time spent on loading pair end info %ds\n", ( int ) ( time_aft - time_bef ) );
		time ( &time_bef );
		Links2Scaf ( graphfile );
		time ( &time_aft );
		printf ( "time spent on creating scaffolds %ds\n", ( int ) ( time_aft - time_bef ) );
		scaffolding ( 100, graphfile );
	}

	prlReadsCloseGap ( graphfile );
	//  locateReadOnScaf(graphfile);
	free_pe_mem();

	if ( index_array )
		{ free ( ( void * ) index_array ); }

	freeContig_array();
	destroyPreArcMem();
	destroyConnectMem();
	deleteCntLookupTable();
	time ( &stop_t );
	printf ( "time elapsed: %dm\n", ( int ) ( stop_t - start_t ) / 60 );
	return 0;
}

/*****************************************************************************
 * Parse command line switches
 *****************************************************************************/


void initenv ( int argc, char ** argv )
{
	int copt;
	int inpseq;
	extern char * optarg;
	char temp[256];
	inpseq = 0;
	LINK = 0;
	SCAFF = 0;
	optind = 1;

	while ( ( copt = getopt ( argc, argv, "g:L:p:G:FuS" ) ) != EOF )
	{
		switch ( copt )
		{
			case 'g':
				inGraph = 1;
				sscanf ( optarg, "%s", graphfile ); //
				continue;
			case 'G':
				sscanf ( optarg, "%s", temp ); //
				GLDiff = atoi ( temp );
				continue;
			case 'L':
				sscanf ( optarg, "%s", temp );
				ctg_short = atoi ( temp );
				continue;
			case 'F':
				fillGap = 1;
				continue;
			case 'S':
				SCAFF = 1;
				continue;
			case 'u':
				maskRep = 0;
				continue;
			case 'p':
				sscanf ( optarg, "%s", temp ); //
				thrd_num = atoi ( temp );
				continue;
			default:

				if ( inGraph == 0 )       //
				{
					display_scaff_usage();
					exit ( -1 );
				}
		}
	}

	if ( inGraph == 0 )   //
	{
		display_scaff_usage();
		exit ( -1 );
	}
}

static void display_scaff_usage()
{
	printf ( "\nscaff -g InputGraph [-F -u -S] [-G gapLenDiff -L minContigLen] [-p n_cpu]\n" );
	printf ( "  -g InputFile: prefix of graph file names\n" );
	printf ( "  -F (optional) fill gaps in scaffold\n" );
	printf ( "  -S (optional) scaffold structure exists(default: NO)\n" );
	printf ( "  -G gapLenDiff(default 50): allowed length difference between estimated and filled gap\n" );
	printf ( "  -u (optional): un-mask contigs with high coverage before scaffolding (default mask)\n" );
	printf ( "  -p n_cpu(default 8): number of cpu for use\n" );
	printf ( "  -L minLen(default K+2): shortest contig for scaffolding\n" );
}
