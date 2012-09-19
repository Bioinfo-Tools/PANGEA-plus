/*
 * 63mer/kmer.c
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

__uint128_t Kmer2int128 ( Kmer seq )
{
	__uint128_t temp;
	temp = seq.high;
	temp <<= 64;
	temp |= seq.low;
	return temp;
}

boolean KmerSmaller ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high < kmer2.high )
		{ return 1; }
	else if ( kmer1.high == kmer2.high )
	{
		if ( kmer1.low < kmer2.low )
			{ return 1; }
		else
			{ return 0; }
	}
	else
		{ return 0; }
}

boolean KmerLarger ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high > kmer2.high )
		{ return 1; }
	else if ( kmer1.high == kmer2.high )
	{
		if ( kmer1.low > kmer2.low )
			{ return 1; }
		else
			{ return 0; }
	}
	else
		{ return 0; }
}

boolean KmerEqual ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high == kmer2.high && kmer1.low == kmer2.low )
		{ return 1; }
	else { return 0; }
}
// kmer1 = kmer1 & kmer2
Kmer KmerAnd ( Kmer kmer1, Kmer kmer2 )
{
	kmer1.high &= kmer2.high;
	kmer1.low &= kmer2.low;
	return kmer1;
}
// kmer <<= 2
Kmer KmerLeftBitMoveBy2 ( Kmer word )
{
	ubyte8 temp = word.low >> 62;
	word.high <<= 2;
	word.high |= temp;
	word.low <<= 2;
	return word;
}

// kmer >>= 2
Kmer KmerRightBitMoveBy2 ( Kmer word )
{
	ubyte8 temp = ( word.high & 0x3 ) << 62;
	word.high >>= 2;
	word.low >>= 2;
	word.low |= temp;
	return word;
}

Kmer KmerPlus ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word.low |= ch;
	return word;
}

Kmer nextKmer ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word = KmerAnd ( word, WORDFILTER );
	word.low |= ch;
	return word;
}

Kmer prevKmer ( Kmer next, char ch )
{
	Kmer word = KmerRightBitMoveBy2 ( next );

	if ( 2 * ( overlaplen - 1 ) < 64 )
		{ word.low |= ( ( ( ubyte8 ) ch ) << 2 * ( overlaplen - 1 ) ); }
	else
		{ word.high |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 64 ); }

	return word;
}

char lastCharInKmer ( Kmer kmer )
{
	return ( char ) ( kmer.low & 0x3 );
}

char firstCharInKmer ( Kmer kmer )
{
	if ( 2 * ( overlaplen - 1 ) < 64 )
	{
		kmer.low >>= 2 * ( overlaplen - 1 );
		return kmer.low;// & 3;
	}
	else
	{
		kmer.high >>= 2 * ( overlaplen - 1 ) - 64;
		return kmer.high;// & 3;
	}
}

Kmer createFilter ( int overlaplen )
{
	Kmer word;
	word.high = word.low = 0;

	if ( 2 * overlaplen < 64 )
		{ word.low = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen ) ) - 1; }
	else
	{
		word.low = ~word.low;

		if ( 2 * overlaplen > 64 )
			{ word.high = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 64 ) ) - 1; }
	}

	return word;
}

Kmer KmerRightBitMove ( Kmer word, int dis )
{
	if ( dis < 64 )
	{
		ubyte8 mask = ( ( ( ubyte8 ) 1 ) << dis ) - 1;
		ubyte8 temp = ( word.high & mask ) << ( 64 - dis );
		word.high >>= dis;
		word.low >>= dis;
		word.low |= temp;
		return word;
	}

	word.high >>= ( dis - 64 );
	word.low = word.high;
	word.high = 0;
	return word;
}

void printKmerSeq ( FILE * fp, Kmer kmer )
{
	int i, bit1, bit2;
	char ch;
	char kmerSeq[64];
	bit2 = overlaplen > 32 ? 32 : overlaplen;
	bit1 = overlaplen > 32 ? overlaplen - 32 : 0;

	for ( i = bit1 - 1; i >= 0; i-- )
	{
		ch = kmer.high & 0x3;
		kmer.high >>= 2;
		kmerSeq[i] = ch;
	}

	for ( i = bit2 - 1; i >= 0; i-- )
	{
		ch = kmer.low & 0x3;
		kmer.low >>= 2;
		kmerSeq[i + bit1] = ch;
	}

	for ( i = 0; i < overlaplen; i++ )
		{ fprintf ( fp, "%c", int2base ( ( int ) kmerSeq[i] ) ); }
}

void print_kmer ( FILE * fp, Kmer kmer, char c )
{
	fprintf ( fp, "%llx %llx", kmer.high, kmer.low );
	fprintf ( fp, "%c", c );
}

static Kmer fastReverseComp ( Kmer seq, char seq_size )
{
	seq.low ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low = ( ( seq.low & 0x3333333333333333LLU ) << 2 ) | ( ( seq.low & 0xCCCCCCCCCCCCCCCCLLU ) >> 2 );
	seq.low = ( ( seq.low & 0x0F0F0F0F0F0F0F0FLLU ) << 4 ) | ( ( seq.low & 0xF0F0F0F0F0F0F0F0LLU ) >> 4 );
	seq.low = ( ( seq.low & 0x00FF00FF00FF00FFLLU ) << 8 ) | ( ( seq.low & 0xFF00FF00FF00FF00LLU ) >> 8 );
	seq.low = ( ( seq.low & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( seq.low & 0xFFFF0000FFFF0000LLU ) >> 16 );
	seq.low = ( ( seq.low & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( seq.low & 0xFFFFFFFF00000000LLU ) >> 32 );

	if ( seq_size < 32 )
	{
		seq.low >>= ( 64 - ( seq_size << 1 ) );
		return seq;
	}

	seq.high ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high = ( ( seq.high & 0x3333333333333333LLU ) << 2 ) | ( ( seq.high & 0xCCCCCCCCCCCCCCCCLLU ) >> 2 );
	seq.high = ( ( seq.high & 0x0F0F0F0F0F0F0F0FLLU ) << 4 ) | ( ( seq.high & 0xF0F0F0F0F0F0F0F0LLU ) >> 4 );
	seq.high = ( ( seq.high & 0x00FF00FF00FF00FFLLU ) << 8 ) | ( ( seq.high & 0xFF00FF00FF00FF00LLU ) >> 8 );
	seq.high = ( ( seq.high & 0x0000FFFF0000FFFFLLU ) << 16 ) | ( ( seq.high & 0xFFFF0000FFFF0000LLU ) >> 16 );
	seq.high = ( ( seq.high & 0x00000000FFFFFFFFLLU ) << 32 ) | ( ( seq.high & 0xFFFFFFFF00000000LLU ) >> 32 );
	ubyte8 temp = seq.high;
	seq.high = seq.low;
	seq.low = temp;
	seq = KmerRightBitMove ( seq, 128 - ( seq_size << 1 ) );
	return seq;
}

Kmer reverseComplementVerbose ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}

Kmer reverseComplement ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}
