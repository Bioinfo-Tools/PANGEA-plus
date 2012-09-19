/*
 * 31mer/cutTip_graph.c
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


static int caseA, caseB, caseC, caseD, caseE;

void destroyEdge ( unsigned int edgeid )
{
	unsigned int bal_ed = getTwinEdge ( edgeid );
	ARC * arc;

	if ( bal_ed == edgeid )
	{
		edge_array[edgeid].length = 0;
		return;
	}

	arc = edge_array[edgeid].arcs;

	while ( arc )
	{
		arc->bal_arc->to_ed = 0;
		arc = arc->next;
	}

	arc = edge_array[bal_ed].arcs;

	while ( arc )
	{
		arc->bal_arc->to_ed = 0;
		arc = arc->next;
	}

	edge_array[edgeid].arcs = NULL;
	edge_array[bal_ed].arcs = NULL;
	edge_array[edgeid].length = 0;
	edge_array[bal_ed].length = 0;
	edge_array[edgeid].deleted = 1;
	edge_array[bal_ed].deleted = 1;
}

ARC * arcCount ( unsigned int edgeid, unsigned int * num )
{
	ARC * arc;
	ARC * firstValidArc = NULL;
	unsigned int count = 0;
	arc = edge_array[edgeid].arcs;

	while ( arc )
	{
		if ( arc->to_ed > 0 )
		{
			count++;

			if ( count == 1 )
				{ firstValidArc = arc; }
			else if ( count > 1 )
			{
				*num = count;
				return firstValidArc;
			}
		}

		arc = arc->next;
	}

	*num = count;
	return firstValidArc;
}

void removeWeakEdges ( int lenCutoff, unsigned int multiCutoff )
{
	unsigned int bal_ed;
	unsigned int arcRight_n, arcLeft_n;
	ARC * arcLeft, *arcRight;
	unsigned int i;
	int counter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].deleted || edge_array[i].length == 0
		        || edge_array[i].length > lenCutoff
		        || EdSameAsTwin ( i ) )
			{ continue; }

		bal_ed = getTwinEdge ( i );
		arcRight = arcCount ( i, &arcRight_n );

		if ( arcRight_n > 1 || !arcRight || arcRight->multiplicity > multiCutoff )
			{ continue; }

		arcLeft = arcCount ( bal_ed, &arcLeft_n );

		if ( arcLeft_n > 1 || !arcLeft || arcLeft->multiplicity > multiCutoff )
			{ continue; }

		destroyEdge ( i );
		counter++;
	}

	printf ( "Remove weakly linked edges: %d weak inner edges destroyed\n", counter );
	removeDeadArcs();
	//linearConcatenate();
	//compactEdgeArray();
}

void removeLowCovEdges ( int lenCutoff, unsigned short covCutoff )
{
	unsigned int bal_ed;
	unsigned int arcRight_n, arcLeft_n;
	ARC * arcLeft, *arcRight;
	unsigned int i;
	int counter = 0;

	for ( i = 1; i <= num_ed; i++ )
	{
		if ( edge_array[i].deleted || edge_array[i].cvg == 0
		        || edge_array[i].cvg > covCutoff * 10
		        || edge_array[i].length >= lenCutoff
		        || EdSameAsTwin ( i ) || edge_array[i].length == 0 )
			{ continue; }

		bal_ed = getTwinEdge ( i );
		arcRight = arcCount ( i, &arcRight_n );
		arcLeft = arcCount ( bal_ed, &arcLeft_n );

		if ( arcLeft_n < 1 || arcRight_n < 1 )
			{ continue; }

		destroyEdge ( i );
		counter++;
	}

	printf ( "Remove low coverage(%d): %d inner edges destroyed\n", covCutoff, counter );
	removeDeadArcs();
	linearConcatenate();
	compactEdgeArray();
}

boolean isUnreliableTip ( unsigned int edgeid, int cutLen, boolean strict )
{
	unsigned int arcRight_n, arcLeft_n;
	unsigned int bal_ed;
	unsigned int currentEd = edgeid;
	int length = 0;
	unsigned int mult = 0;
	ARC * arc, *activeArc = NULL, *tempArc;

	if ( edgeid == 0 )
		{ return 0; }

	bal_ed = getTwinEdge ( edgeid );

	if ( bal_ed == edgeid )
		{ return 0; }

	arcCount ( bal_ed, &arcLeft_n );

	if ( arcLeft_n > 0 )
		{ return 0; }

	while ( currentEd )
	{
		arcCount ( bal_ed, &arcLeft_n );
		tempArc = arcCount ( currentEd, &arcRight_n );

		if ( arcLeft_n > 1 || arcRight_n > 1 )
			{ break; }

		length += edge_array[currentEd].length;

		if ( tempArc )
		{
			activeArc = tempArc;
			currentEd = activeArc->to_ed;
			bal_ed = getTwinEdge ( currentEd );
		}
		else
			{ currentEd = 0; }
	}

	if ( length >= cutLen )
	{
		return 0;
	}

	if ( currentEd == 0 )
	{
		caseB++;
		return 1;
	}

	if ( !strict )
	{
		if ( arcLeft_n < 2 )
			{ length += edge_array[currentEd].length; }

		if ( length >= cutLen )
			{ return 0; }
		else
		{
			caseC++;
			return 1;
		}
	}

	if ( arcLeft_n < 2 )
	{
		return 0;
	}

	if ( !activeArc )
		{ printf ( "no activeArc while checking edge %d\n", edgeid ); }

	if ( activeArc->multiplicity == 1 )
	{
		caseD++;
		return 1;
	}

	for ( arc = edge_array[bal_ed].arcs; arc != NULL; arc = arc->next )
		if ( arc->multiplicity > mult )
			{ mult = arc->multiplicity; }

	if ( mult > activeArc->multiplicity )
		{ caseE++; }

	return mult > activeArc->multiplicity;
}

boolean isUnreliableTip_strict ( unsigned int edgeid, int cutLen )
{
	unsigned int arcRight_n, arcLeft_n;
	unsigned int bal_ed;
	unsigned int currentEd = edgeid;
	int length = 0;
	unsigned int mult = 0;
	ARC * arc, *activeArc = NULL, *tempArc;

	if ( edgeid == 0 )
		{ return 0; }

	bal_ed = getTwinEdge ( edgeid );

	if ( bal_ed == edgeid )
		{ return 0; }

	arcCount ( bal_ed, &arcLeft_n );

	if ( arcLeft_n > 0 )
		{ return 0; }

	while ( currentEd )
	{
		arcCount ( bal_ed, &arcLeft_n );
		tempArc = arcCount ( currentEd, &arcRight_n );

		if ( arcLeft_n > 1 || arcRight_n > 1 )
		{
			if ( arcLeft_n == 0 || length == 0 )
				{ return 0; }
			else
				{ break; }
		}

		length += edge_array[currentEd].length;

		if ( length >= cutLen )
			{ return 0; }

		if ( tempArc )
		{
			activeArc = tempArc;
			currentEd = activeArc->to_ed;
			bal_ed = getTwinEdge ( currentEd );
		}
		else
			{ currentEd = 0; }
	}

	if ( currentEd == 0 )
	{
		caseA++;
		return 1;
	}

	if ( !activeArc )
		{ printf ( "no activeArc while checking edge %d\n", edgeid ); }

	if ( activeArc->multiplicity == 1 )
	{
		caseB++;
		return 1;
	}

	for ( arc = edge_array[bal_ed].arcs; arc != NULL; arc = arc->next )
		if ( arc->multiplicity > mult )
			{ mult = arc->multiplicity; }

	if ( mult > activeArc->multiplicity )
		{ caseC++; }

	return mult > activeArc->multiplicity;
}

void removeDeadArcs()
{
	unsigned int i, count = 0;
	ARC * arc, *arc_temp;

	for ( i = 1; i <= num_ed; i++ )
	{
		arc = edge_array[i].arcs;

		while ( arc )
		{
			arc_temp = arc;
			arc = arc->next;

			if ( arc_temp->to_ed == 0 )
			{
				count++;
				edge_array[i].arcs = deleteArc ( edge_array[i].arcs, arc_temp );
			}
		}
	}

	printf ( "%d dead arcs removed\n", count );
}

void cutTipsInGraph ( int cutLen, boolean strict )
{
	int flag = 1;
	unsigned int i;

	if ( !cutLen )
		{ cutLen = 2 * overlaplen; }

	printf ( "strict %d, cutLen %d\n", strict, cutLen );

	if ( strict )
		{ linearConcatenate(); }

	caseA = caseB = caseC = caseD = caseE = 0;

	while ( flag )
	{
		flag = 0;

		for ( i = 1; i <= num_ed; i++ )
		{
			if ( edge_array[i].deleted )
				{ continue; }

			/*
			if(strict){
			    if(isUnreliableTip_strict(i,cutLen)){
			        destroyEdge(i);
			        flag++;
			    }
			}else
			*/
			if ( isUnreliableTip ( i, cutLen, strict ) )
			{
				destroyEdge ( i );
				flag++;
			}
		}

		printf ( "a cutTipsInGraph lap, %d tips cut\n", flag );
	}

	removeDeadArcs();

	if ( strict )
		{ printf ( "case A %d, B %d C %d D %d E %d\n", caseA, caseB, caseC, caseD, caseE ); }

	linearConcatenate();
	compactEdgeArray();
}
