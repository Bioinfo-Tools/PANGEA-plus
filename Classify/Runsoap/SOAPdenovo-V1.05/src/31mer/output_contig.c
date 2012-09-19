/*
 * 31mer/output_contig.c
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

static char * kmerSeq;

void output_graph ( char * outfile )
{
	char name[256];
	FILE * fp;
	unsigned int i, bal_i;
	sprintf ( name, "%s.edge.gvz", outfile );
	fp = ckopen ( name, "w" );
	fprintf ( fp, "digraph G{\n" );
	fprintf ( fp, "\tsize=\"512,512\";\n" );

	for ( i = num_ed; i > 0; i-- )
	{
		if ( edge_array[i].deleted )
			{ continue; }

		/*
		        arcCount(i,&arcNum);
		        if(arcNum<1)
		            continue;
		*/
		bal_i = getTwinEdge ( i );
		/*
		        arcCount(bal_i,&arcNum);
		        if(arcNum<1)
		            continue;
		*/
		fprintf ( fp, "\tV%d -> V%d[label =\"%d(%d)\"];\n", edge_array[i].from_vt, edge_array[i].to_vt, i, edge_array[i].length );
	}

	fprintf ( fp, "}\n" );
	fclose ( fp );
}

static void output_1contig ( int id, EDGE * edge, FILE * fp, boolean tip )
{
	int i;
	Kmer kmer;
	char ch;
	fprintf ( fp, ">%d length %d cvg_%.1f_tip_%d\n", id, edge->length + overlaplen, ( double ) edge->cvg / 10, tip );
	//fprintf(fp,">%d\n",id);
	kmer = vt_array[edge->from_vt].kmer;

	for ( i = overlaplen - 1; i >= 0; i-- )
	{
		ch = kmer & 3;
		kmer >>= 2;
		kmerSeq[i] = ch;
	}

	for ( i = 0; i < overlaplen; i++ )
		{ fprintf ( fp, "%c", int2base ( ( int ) kmerSeq[i] ) ); }

	for ( i = 0; i < edge->length; i++ )
	{
		fprintf ( fp, "%c", int2base ( ( int ) getCharInTightString ( edge->seq, i ) ) );

		if ( ( i + overlaplen + 1 ) % 100 == 0 )
			{ fprintf ( fp, "\n" ); }
	}

	fprintf ( fp, "\n" );
}

int cmp_int ( const void * a, const void * b )
{
	int * A, *B;
	A = ( int * ) a;
	B = ( int * ) b;

	if ( *A > *B )
		{ return 1; }
	else if ( *A == *B )
		{ return 0; }
	else
		{ return -1; }
}

int cmp_edge ( const void * a, const void * b )
{
	EDGE * A, *B;
	A = ( EDGE * ) a;
	B = ( EDGE * ) b;

	if ( A->length > B->length )
		{ return 1; }
	else if ( A->length == B->length )
		{ return 0; }
	else
		{ return -1; }
}
void output_contig ( EDGE * ed_array, unsigned int ed_num, char * outfile, int cut_len )
{
	char    temp[256];
	FILE * fp, *fp_contig;
	int flag, count, len_c;
	int signI;
	unsigned int i;
	long long sum = 0, N90, N50;
	unsigned int * length_array;
	boolean tip;
	sprintf ( temp, "%s.contig", outfile );
	fp = ckopen ( temp, "w" );
	qsort ( &ed_array[1], ed_num, sizeof ( EDGE ), cmp_edge );
	length_array = ( unsigned int * ) ckalloc ( ed_num * sizeof ( unsigned int ) );
	kmerSeq = ( char * ) ckalloc ( overlaplen * sizeof ( char ) );
	//first scan for number counting
	count = len_c = 0;

	for ( i = 1; i <= ed_num; i++ )
	{
		if ( ( ed_array[i].length + overlaplen ) >= len_bar )
			{ length_array[len_c++] = ed_array[i].length + overlaplen; }

		if ( ed_array[i].length < 1 || ed_array[i].deleted )
			{ continue; }

		count++;

		if ( EdSmallerThanTwin ( i ) )
			{ i++; }
	}

	sum = 0;

	for ( signI = len_c - 1; signI >= 0; signI-- )
		{ sum += length_array[signI]; }

	if ( len_c > 0 )
		{ printf ( "%d ctgs longer than %d, sum up %lldbp, with average length %lld\n", len_c, len_bar, sum, sum / len_c ); }

	qsort ( length_array, len_c, sizeof ( length_array[0] ), cmp_int );
	printf ( "the longest is %dbp, ", length_array[len_c - 1] );
	N50 = sum * 0.5;
	N90 = sum * 0.9;
	sum = flag = 0;

	for ( signI = len_c - 1; signI >= 0; signI-- )
	{
		sum += length_array[signI];

		if ( !flag && sum >= N50 )
		{
			printf ( "contig N50 is %d bp,", length_array[signI] );
			flag = 1;
		}

		if ( sum >= N90 )
		{
			printf ( "contig N90 is %d bp\n", length_array[signI] );
			break;
		}
	}

	//fprintf(fp,"Number %d\n",count);

	for ( i = 1; i <= ed_num; i++ )
	{
		//if(ed_array[i].multi!=1||ed_array[i].length<1||(ed_array[i].length+overlaplen)<cut_len)
		if ( ed_array[i].deleted || ed_array[i].length < 1 )
			{ continue; }

		if ( ed_array[i].arcs && ed_array[getTwinEdge ( i )].arcs )
			{ tip = 0; }
		else
			{ tip = 1; }

		output_1contig ( i, & ( ed_array[i] ), fp, tip );

		if ( EdSmallerThanTwin ( i ) )
			{ i++; }
	}

	fclose ( fp );
	free ( ( void * ) kmerSeq );
	free ( ( void * ) length_array );
	printf ( "%d contigs longer than %d output\n", count, cut_len );
	sprintf ( temp, "%s.ContigIndex", outfile );
	fp_contig = ckopen ( temp, "w" );
	fprintf ( fp_contig, "Edge_num %d %d\n", ed_num, count );
	fprintf ( fp_contig, "index\tlength\treverseComplement\n" );

	for ( i = 1; i <= num_ed; i++ )
	{
		fprintf ( fp_contig, "%d\t%d\t", i, edge_array[i].length + overlaplen );

		if ( EdSmallerThanTwin ( i ) )
		{
			fprintf ( fp_contig, "1\n" );
			i++;
		}
		else if ( EdLargerThanTwin ( i ) )
			{ fprintf ( fp_contig, "-1\n" ); }
		else
			{ fprintf ( fp_contig, "0\n" ); }
	}

	fclose ( fp_contig );
}

void output_updated_edges ( char * outfile )
{
	FILE * fp;
	char name[256];
	unsigned int i, validCounter = 0;
	EDGE * edge;
	sprintf ( name, "%s.updated.edge", outfile );
	fp = ckopen ( name, "w" );

	for ( i = 1; i <= num_ed; i++ )
		{ validCounter++; }

	fprintf ( fp, "EDGEs %d\n", validCounter );
	validCounter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		edge = &edge_array[i];
		fprintf ( fp, ">length %d,", edge->length );
		print_kmer ( fp, vt_array[edge->from_vt].kmer, ',' );
		print_kmer ( fp, vt_array[edge->to_vt].kmer, ',' );

		if ( EdSmallerThanTwin ( i ) )
			{ fprintf ( fp, "1," ); }
		else if ( EdLargerThanTwin ( i ) )
			{ fprintf ( fp, "-1," ); }
		else
			{ fprintf ( fp, "0," ); }

		fprintf ( fp, "%d\n", edge->cvg );
	}

	fclose ( fp );
}

void output_heavyArcs ( char * outfile )
{
	unsigned int i, j;
	char name[256];
	FILE * outfp;
	ARC * parc;
	sprintf ( name, "%s.Arc", outfile );
	outfp = ckopen ( name, "w" );

	for ( i = 1; i <= num_ed; i++ )
	{
		parc = edge_array[i].arcs;

		if ( !parc )
			{ continue; }

		j = 0;
		fprintf ( outfp, "%u", i );

		while ( parc )
		{
			fprintf ( outfp, " %u %u", parc->to_ed, parc->multiplicity );

			if ( ( ++j ) % 10 == 0 )
			{
				fprintf ( outfp, "\n%u", i );
			}

			parc = parc->next;
		}

		fprintf ( outfp, "\n" );
	}

	fclose ( outfp );
}

