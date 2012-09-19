/*
 * 127mer/map.c
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

static void display_map_usage();

static int getMinOverlap ( char * gfile )
{
	char name[256], ch;
	FILE * fp;
	int num_kmer, overlaplen = 23;
	char line[1024];
	sprintf ( name, "%s.preGraphBasic", gfile );
	fp = fopen ( name, "r" );

	if ( !fp )
		{ return overlaplen; }

	while ( fgets ( line, sizeof ( line ), fp ) != NULL )
	{
		if ( line[0] == 'V' )
		{
			sscanf ( line + 6, "%d %c %d", &num_kmer, &ch, &overlaplen );
		}
		else if ( line[0] == 'M' )
		{
			sscanf ( line, "MaxReadLen %d MinReadLen %d MaxNameLen %d", &maxReadLen, &minReadLen, &maxNameLen );
		}
	}

	fclose ( fp );
	return overlaplen;
}

int call_align ( int argc, char ** argv )
{
	time_t start_t, stop_t, time_bef, time_aft;
	time ( &start_t );
	initenv ( argc, argv );
	overlaplen = getMinOverlap ( graphfile );
	printf ( "K = %d\n", overlaplen );
	time ( &time_bef );
	ctg_short = overlaplen + 2;
	printf ( "contig len cutoff: %d\n", ctg_short );
	prlContig2nodes ( graphfile, ctg_short );
	time ( &time_aft );
	printf ( "time spent on De bruijn graph construction: %ds\n\n",
	         ( int ) ( time_aft - time_bef ) );
	//map long read to edge one by one
	time ( &time_bef );
	prlLongRead2Ctg ( shortrdsfile, graphfile );
	time ( &time_aft );
	printf ( "time spent on mapping long reads: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	//map read to edge one by one
	time ( &time_bef );
	prlRead2Ctg ( shortrdsfile, graphfile );
	time ( &time_aft );
	printf ( "time spent on mapping reads: %ds\n\n", ( int ) ( time_aft - time_bef ) );
	free_Sets ( KmerSets, thrd_num );
	time ( &stop_t );
	printf ( "overall time for alignment: %dm\n\n", ( int ) ( stop_t - start_t ) / 60 );
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

	while ( ( copt = getopt ( argc, argv, "s:g:K:p:" ) ) != EOF )
	{
		//printf("get option\n");
		switch ( copt )
		{
			case 's':
				inpseq = 1;
				sscanf ( optarg, "%s", shortrdsfile );
				continue;
			case 'g':
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
			default:

				if ( inpseq == 0 || outseq == 0 )       //
				{
					display_map_usage();
					exit ( 1 );
				}
		}
	}

	if ( inpseq == 0 || outseq == 0 )  //
	{
		//printf("need more\n");
		display_map_usage();
		exit ( 1 );
	}
}

static void display_map_usage()
{
	printf ( "\nmap -s readsInfoFile -g graphfile [-p n_cpu]\n" );
	printf ( "  -s readsInfoFile: The file contains information of solexa reads\n" );
	printf ( "  -p n_cpu(default 8): number of cpu for use\n" );
	printf ( "  -g graphfile: prefix of graph files\n" );
}
