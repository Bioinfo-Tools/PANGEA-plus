/*
 * 63mer/contig.c
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
static void display_contig_usage();
char shortrdsfile[256], graphfile[256];

static boolean repeatSolve;
static int M = 1;

int call_heavygraph ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	boolean ret;
	initenv ( argc, argv );
	loadVertex ( graphfile );
	loadEdge ( graphfile );

	if ( repeatSolve )
	{
		time ( &time_bef );
		ret = loadPathBin ( graphfile );

		if ( ret )
		{
			solveReps();
		}
		else
			{ printf ( "repeat solving can't be done...\n" ); }

		time ( &time_aft );
		printf ( "time spent on solving repeat: %ds\n", ( int ) ( time_aft - time_bef ) );
	}

	//edgecvg_bar(edge_array,num_ed,graphfile,100);
	//0531
	if ( M > 0 )
	{
		time ( &time_bef );
		bubblePinch ( 0.90, graphfile, M );
		time ( &time_aft );
		printf ( "time spent on bubblePinch: %ds\n", ( int ) ( time_aft - time_bef ) );
	}

	removeWeakEdges ( 2 * overlaplen, 1 );
	removeLowCovEdges ( 2 * overlaplen, 1 );
	cutTipsInGraph ( 0, 0 );
	//output_graph(graphfile);
	output_contig ( edge_array, num_ed, graphfile, overlaplen + 1 );
	output_updated_edges ( graphfile );
	output_heavyArcs ( graphfile );

	if ( vt_array )
	{
		free ( ( void * ) vt_array );
		vt_array = NULL;
	}

	if ( edge_array )
	{
		free_edge_array ( edge_array, num_ed_limit );
		edge_array = NULL;
	}

	destroyArcMem();
	time ( &stop_t );
	printf ( "time elapsed: %dm\n\n", ( int ) ( stop_t - start_t ) / 60 );
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
	inpseq = outseq = repeatSolve = 0;
	optind = 1;

	while ( ( copt = getopt ( argc, argv, "g:M:R" ) ) != EOF )
	{
		switch ( copt )
		{
			case 'M':
				sscanf ( optarg, "%s", temp ); //
				M = atoi ( temp );
				continue;
			case 'g':
				inGraph = 1;
				sscanf ( optarg, "%s", graphfile ); //
				continue;
			case 'R':
				repeatSolve = 1;
				continue;
			default:

				if ( inGraph == 0 )       //
				{
					display_contig_usage();
					exit ( -1 );
				}
		}
	}

	if ( inGraph == 0 )   //
	{
		display_contig_usage();
		exit ( -1 );
	}
}

static void display_contig_usage()
{
	printf ( "\ncontig -g InputGraph [-M mergeLevel -R]\n" );
	printf ( "  -g InputFile: prefix of graph file names\n" );
	printf ( "  -M mergeLevel(default 1,min 0, max 3): the strength of merging similar sequences during contiging\n" );
	printf ( "  -R solve_repeats (optional): solve repeats by read paths(default: no)\n" );
}
