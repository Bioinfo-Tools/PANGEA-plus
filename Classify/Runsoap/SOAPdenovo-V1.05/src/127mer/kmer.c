/*
 * 127mer/kmer.c
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

/*
U256b Kmer2int256(Kmer seq)
{
    U256b temp;
    temp.low = seq.high2;
    temp.low <<= 64;
    temp.low |= seq.low2;
    temp.high = seq.high1;
    temp.high <<= 64;
    temp.high |= seq.low1;

    //temp.low  = seq.high2 << 64 | seq.low2;
    //temp.high = seq.high1 << 64 | seq.low1;
    return temp;
} */


boolean KmerSmaller ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high1 != kmer2.high1 )  { return ( kmer1.high1 < kmer2.high1 ); }
	else
	{
		if ( kmer1.low1 != kmer2.low1 ) { return ( kmer1.low1 < kmer2.low1 ); }
		else
		{
			if ( kmer1.high2 != kmer2.high2 )
				{ return ( kmer1.high2 < kmer2.high2 ); }
			else  { return ( kmer1.low2 < kmer2.low2 ); }
		}
	}

	/*  if(kmer1.high1<kmer2.high1)
	        return 1;
	    else if(kmer1.high1==kmer2.high1){
	        if(kmer1.low1<kmer2.low1)
	            return 1;
	        else if(kmer1.low1==kmer2.low1){
	            if(kmer1.high2<kmer2.high2)
	                return 1;
	            else if(kmer1.high2==kmer2.high2){
	                if(kmer1.low2<kmer2.low2)
	                    return 1;
	                else return 0;
	            }else return 0;
	        }else return 0;
	    }else return 0;
	*/
}

boolean KmerLarger ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.high1 != kmer2.high1 )  { return ( kmer1.high1 > kmer2.high1 ); }
	else
	{
		if ( kmer1.low1 != kmer2.low1 ) { return ( kmer1.low1 > kmer2.low1 ); }
		else
		{
			if ( kmer1.high2 != kmer2.high2 )
				{ return ( kmer1.high2 > kmer2.high2 ); }
			else  { return ( kmer1.low2 > kmer2.low2 ); }
		}
	}

	/*  if(kmer1.high1>kmer2.high1)
	        return 1;
	    else if(kmer1.high1==kmer2.high1){
	        if(kmer1.low1>kmer2.low1)
	            return 1;
	        else if(kmer1.low1==kmer2.low1){
	            if(kmer1.high2>kmer2.high2)
	                return 1;
	            else if(kmer1.high2==kmer2.high2){
	                if(kmer1.low2>kmer2.low2)
	                    return 1;
	                else return 0;
	            }else return 0;
	        }else return 0;
	    }else return 0;
	*/
}

boolean KmerEqual ( Kmer kmer1, Kmer kmer2 )
{
	if ( kmer1.low2 != kmer2.low2 || kmer1.high2 != kmer2.high2 || kmer1.low1 != kmer2.low1 || kmer1.high1 != kmer2.high1 )
		{ return 0; }
	else { return 1; }
}
// kmer1 = kmer1 & kmer2
Kmer KmerAnd ( Kmer kmer1, Kmer kmer2 )
{
	kmer1.high1 &= kmer2.high1;
	kmer1.low1 &= kmer2.low1;
	kmer1.high2 &= kmer2.high2;
	kmer1.low2 &= kmer2.low2;
	return kmer1;
}
// kmer <<= 2
Kmer KmerLeftBitMoveBy2 ( Kmer word )
{
	word.high1 = ( word.high1 << 2 ) | ( word.low1 >> 62 );
	word.low1  = ( word.low1  << 2 ) | ( word.high2 >> 62 );
	word.high2 = ( word.high2 << 2 ) | ( word.low2 >> 62 );
	word.low2 <<= 2;
	return word;
}

// kmer >>= 2
Kmer KmerRightBitMoveBy2 ( Kmer word )
{
	word.low2  = ( word.low2  >> 2 ) | ( word.high2 & 0x3 ) << 62;
	word.high2 = ( word.high2 >> 2 ) | ( word.low1  & 0x3 ) << 62;
	word.low1  = ( word.low1  >> 2 ) | ( word.high1 & 0x3 ) << 62;
	word.high1 >>= 2;
	return word;
}

Kmer KmerPlus ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word.low2 |= ch;
	return word;
}

Kmer nextKmer ( Kmer prev, char ch )
{
	Kmer word = KmerLeftBitMoveBy2 ( prev );
	word = KmerAnd ( word, WORDFILTER );
	word.low2 |= ch;
	return word;
}

Kmer prevKmer ( Kmer next, char ch )
{
	Kmer word = KmerRightBitMoveBy2 ( next );

	switch ( overlaplen )
	{
		case 1  ... 32:
			word.low2 |= ( ( ( ubyte8 ) ch ) << 2 * ( overlaplen - 1 ) );
			break;
		case 33 ... 64:
			word.high2 |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 64 );
			break;
		case 65 ... 96 :
			word.low1 |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 128 );
			break;
		case 97 ... 128:
			word.high1 |= ( ( ubyte8 ) ch ) << ( 2 * ( overlaplen - 1 ) - 192 );
			break;
	}

	return word;
}

char lastCharInKmer ( Kmer kmer )
{
	return ( char ) ( kmer.low2 & 0x3 );
}

char firstCharInKmer ( Kmer kmer )
{
	switch ( overlaplen )
	{
		case 1  ... 32:
			kmer.low2 >>= 2 * ( overlaplen - 1 );
			return kmer.low2;// & 3;
		case 33 ... 64:
			kmer.high2 >>= 2 * ( overlaplen - 1 ) - 64;
			return kmer.high2;// & 3;
		case 65 ... 96 :
			kmer.low1 >>= 2 * ( overlaplen - 1 ) - 128;
			return kmer.low1;
		case 97 ... 128:
			kmer.high1 >>= 2 * ( overlaplen - 1 ) - 192;
			return kmer.high1;
	}
}

Kmer createFilter ( int overlaplen )
{
	Kmer word;
	word.high1 = word.low1 = word.high2 = word.low2 = 0;

	switch ( overlaplen )
	{
		case 1  ... 31:
			word.low2 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen ) ) - 1;
			break;
		case 32 ... 63:
			word.low2 = ~word.low2;
			word.high2 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 64 ) ) - 1;
			break;
		case 64 ... 95:
			word.high2 = word.low2 = ~word.low2;
			word.low1 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 128 ) ) - 1;
			break;
		case 96 ... 127:
			word.low1 = word.high2 = word.low2 = ~word.low2;
			word.high1 = ( ( ( ubyte8 ) 1 ) << ( 2 * overlaplen - 192 ) ) - 1;
			break;
	}

	return word;
}

Kmer KmerRightBitMove ( Kmer word, int dis )
{
	ubyte8 mask;

	switch ( dis )
	{
		case 1 ... 63:
			mask = ( ( ( ubyte8 ) 1 ) << dis ) - 1;
			word.low2  = ( word.low2  >> dis ) | ( word.high2 & mask ) << ( 64 - dis );
			word.high2 = ( word.high2 >> dis ) | ( word.low1  & mask ) << ( 64 - dis );
			word.low1  = ( word.low1  >> dis ) | ( word.high1 & mask ) << ( 64 - dis );
			word.high1 >>= dis;
			return word;
		case 64 ... 127:
			mask = ( ( ( ubyte8 ) 1 ) << ( dis - 64 ) ) - 1;
			word.low2  = word.high2 >> ( dis - 64 ) | ( word.low1 & mask ) << ( 128 - dis );
			word.high2 = word.low1  >> ( dis - 64 ) | ( word.high1 & mask ) << ( 128 - dis );
			word.low1  = word.high1 >> ( dis - 64 );
			word.high1 = 0;
			return word;
		case 128 ... 191:
			mask = ( ( ( ubyte8 ) 1 ) << ( dis - 128 ) ) - 1;
			word.low2  = word.low1  >> ( dis - 128 ) | ( word.high1 & mask ) << ( 192 - dis );
			word.high2 = word.high1 >> ( dis - 128 );
			word.high1 = word.low1 = 0;
			return word;
		case 192 ... 255:
			word.low2 = word.high1 >> ( dis - 192 );
			word.high1 = word.low1 = word.high2 = 0;
			return word;
	}

	/*  if(dis<64){
	        ubyte8 mask = (((ubyte8) 1) << dis) - 1;
	        ubyte8 temp1 = (word.high1&mask)<<(64-dis);
	        ubyte8 temp2 = (word.low1&mask)<<(64-dis);
	        ubyte8 temp3 = (word.high2&mask)<<(64-dis);
	        word.high1 >>= dis;
	        word.low1 >>= dis;
	        word.high2 >>= dis;
	        word.low2 >>= dis;
	        word.low1 |= temp1;
	        word.high2 |= temp2;
	        word.low2 |= temp3;
	        return word;
	    }
	    if(dis>=64 && dis<128){
	        ubyte8 mask = (((ubyte8) 1) << (dis-64)) - 1;
	        ubyte8 temp1 = (word.high1&mask)<<(128-dis);
	        ubyte8 temp2 = (word.low1&mask)<<(128-dis);
	        word.high1 >>= (dis-64);
	        word.low1 >>= (dis-64);
	        word.high2 >>= (dis-64);
	        word.high2 |= temp2;
	        word.low2 = word.high2;
	        word.low1 |= temp1;
	        word.high2 = word.low1;
	        word.low1 = word.high1;
	        word.high1 = 0;
	        return word;
	    }
	    if(dis>=128 && dis<192){
	        ubyte8 mask = (((ubyte8) 1) << (dis-128)) - 1;
	        ubyte8 temp1 = (word.high1&mask)<<(192-dis);
	        word.high1 >>= (dis-128);
	        word.low1 >>= (dis-128);
	        word.low1 |= temp1;
	        word.low2 = word.low1;
	        word.high2 = word.high1;
	        word.low1 =0;
	        word.high1 = 0;
	        return word;
	    }
	    if(dis>=192 && dis<256){
	        word.high1 >>= (dis-192);
	        word.low2 = word.high1;
	        word.high2 = 0;
	        word.low1 =0;
	        word.high1 = 0;
	        return word;
	    }*/
}

void printKmerSeq ( FILE * fp, Kmer kmer )
{
	int i, bit1, bit2, bit3, bit4;
	bit4 = bit3 = bit2 = bit1 = 0;
	char kmerSeq[128];

	switch ( overlaplen )
	{
		case 1  ... 31:
			bit4 = overlaplen;
			break;
		case 32 ... 63:
			bit4 = 32;
			bit3 = overlaplen - 32;
			break;
		case 64 ... 95:
			bit4 = bit3 = 32;
			bit2 = overlaplen - 64;
			break;
		case 96 ... 127:
			bit4 = bit3 = bit2 = 32;
			bit1 = overlaplen - 96;
			break;
	}

	for ( i = bit1 - 1; i >= 0; i-- )
	{
		kmerSeq[i] = kmer.high1 & 0x3;
		kmer.high1 >>= 2;
	}

	for ( i = bit2 - 1; i >= 0; i-- )
	{
		kmerSeq[i + bit1] = kmer.low1 & 0x3;
		kmer.low1 >>= 2;
	}

	for ( i = bit3 - 1; i >= 0; i-- )
	{
		kmerSeq[i + bit1 + bit2] = kmer.high2 & 0x3;
		kmer.high2 >>= 2;
	}

	for ( i = bit4 - 1; i >= 0; i-- )
	{
		kmerSeq[i + bit1 + bit2 + bit3] = kmer.low2 & 0x3;
		kmer.low2 >>= 2;
	}

	for ( i = 0; i < overlaplen; i++ )
		{ fprintf ( fp, "%c", int2base ( ( int ) kmerSeq[i] ) ); }
}

void print_kmer ( FILE * fp, Kmer kmer, char c )
{
	fprintf ( fp, "%llx %llx %llx %llx", kmer.high1, kmer.low1, kmer.high2, kmer.low2 );
	fprintf ( fp, "%c", c );
}

static const ubyte2 BitReverseTable [65536] =
{
# define R2(n)      n,     n + 1*16384,    n + 2*16384,    n + 3*16384
# define R4(n)   R2(n), R2(n + 1*4096), R2(n + 2*4096), R2(n + 3*4096)
# define R6(n)   R4(n), R4(n + 1*1024), R4(n + 2*1024), R4(n + 3*1024)
# define R8(n)   R6(n), R6(n + 1*256 ), R6(n + 2*256 ), R6(n + 3*256 )
# define R10(n)  R8(n), R8(n + 1*64  ), R8(n + 2*64  ), R8(n + 3*64  )
# define R12(n)  R10(n),R10(n + 1*16), R10(n + 2*16 ), R10(n + 3*16  )
# define R14(n)  R12(n),R12(n + 1*4 ), R12(n + 2*4  ), R12(n + 3*4  )
	R14 ( 0 ),  R14 ( 1 ),  R14 ( 2 ),   R14 ( 3 )
};



static Kmer fastReverseComp ( Kmer seq, char seq_size )
{
	seq.low2 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low2 = ( ( ubyte8 ) BitReverseTable[ seq.low2 & 0xffff] << 48 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low2 >> 16 ) & 0xffff] << 32 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low2 >> 32 ) & 0xffff] << 16 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low2 >> 48 ) & 0xffff] );

	if ( seq_size < 32 )
	{
		seq.low2 >>= ( 64 - ( seq_size << 1 ) );
		return seq;
	}

	seq.high2 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high2 = ( ( ubyte8 ) BitReverseTable[ seq.high2 & 0xffff] << 48 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high2 >> 16 ) & 0xffff] << 32 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high2 >> 32 ) & 0xffff] << 16 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high2 >> 48 ) & 0xffff] );

	if ( seq_size < 64 )
	{
		seq.high2 = seq.high2 ^ seq.low2;
		seq.low2  = seq.high2 ^ seq.low2;
		seq.high2 = seq.high2 ^ seq.low2;
		seq = KmerRightBitMove ( seq, 128 - ( seq_size << 1 ) );
		return seq;
	}

	seq.low1 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.low1 = ( ( ubyte8 ) BitReverseTable[ seq.low1 & 0xffff] << 48 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low1 >> 16 ) & 0xffff] << 32 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low1 >> 32 ) & 0xffff] << 16 ) |
	           ( ( ubyte8 ) BitReverseTable[ ( seq.low1 >> 48 ) & 0xffff] );

	if ( seq_size < 96 )
	{
		seq.low1 = seq.low1 ^ seq.low2;
		seq.low2 = seq.low1 ^ seq.low2;
		seq.low1 = seq.low1 ^ seq.low2;
		seq = KmerRightBitMove ( seq, 192 - ( seq_size << 1 ) );
		return seq;
	}

	seq.high1 ^= 0xAAAAAAAAAAAAAAAALLU;
	seq.high1 = ( ( ubyte8 ) BitReverseTable[ seq.high1 & 0xffff] << 48 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high1 >> 16 ) & 0xffff] << 32 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high1 >> 32 ) & 0xffff] << 16 ) |
	            ( ( ubyte8 ) BitReverseTable[ ( seq.high1 >> 48 ) & 0xffff] );
	seq.low1  = seq.low1 ^ seq.high2;
	seq.high2 = seq.low1 ^ seq.high2;
	seq.low1  = seq.low1 ^ seq.high2;
	seq.low2  = seq.low2 ^ seq.high1;
	seq.high1 = seq.low2 ^ seq.high1;
	seq.low2  = seq.low2 ^ seq.high1;
	seq = KmerRightBitMove ( seq, 256 - ( seq_size << 1 ) );
	return seq;
}
/*
seq.low2 ^= 0xAAAAAAAAAAAAAAAALLU;
seq.low2 = ((seq.low2 & 0x3333333333333333LLU)<< 2) | ((seq.low2 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
seq.low2 = ((seq.low2 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.low2 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
seq.low2 = ((seq.low2 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.low2 & 0xFF00FF00FF00FF00LLU)>> 8);
seq.low2 = ((seq.low2 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.low2 & 0xFFFF0000FFFF0000LLU)>>16);
seq.low2 = ((seq.low2 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.low2 & 0xFFFFFFFF00000000LLU)>>32);
if(seq_size<32){
    seq.low2 >>= (64 - (seq_size<<1));
    return seq;
}
seq.high2 ^= 0xAAAAAAAAAAAAAAAALLU;
seq.high2 = ((seq.high2 & 0x3333333333333333LLU)<< 2) | ((seq.high2 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
seq.high2 = ((seq.high2 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.high2 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
seq.high2 = ((seq.high2 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.high2 & 0xFF00FF00FF00FF00LLU)>> 8);
seq.high2 = ((seq.high2 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.high2 & 0xFFFF0000FFFF0000LLU)>>16);
seq.high2 = ((seq.high2 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.high2 & 0xFFFFFFFF00000000LLU)>>32);
ubyte8 temp;
if(seq_size<64){
    temp = seq.high2;
    seq.high2 = seq.low2;
    seq.low2 = temp;
    seq = KmerRightBitMove(seq,128-(seq_size<<1));
    return seq;
}
seq.low1 ^= 0xAAAAAAAAAAAAAAAALLU;
seq.low1 = ((seq.low1 & 0x3333333333333333LLU)<< 2) | ((seq.low1 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
seq.low1 = ((seq.low1 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.low1 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
seq.low1 = ((seq.low1 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.low1 & 0xFF00FF00FF00FF00LLU)>> 8);
seq.low1 = ((seq.low1 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.low1 & 0xFFFF0000FFFF0000LLU)>>16);
seq.low1 = ((seq.low1 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.low1 & 0xFFFFFFFF00000000LLU)>>32);
if(seq_size<96){
    temp = seq.low2;
    seq.low2 = seq.low1;
    seq.low1 = temp;
    seq = KmerRightBitMove(seq,192-(seq_size<<1));
    return seq;
}
seq.high1 ^= 0xAAAAAAAAAAAAAAAALLU;
seq.high1 = ((seq.high1 & 0x3333333333333333LLU)<< 2) | ((seq.high1 & 0xCCCCCCCCCCCCCCCCLLU)>> 2);
seq.high1 = ((seq.high1 & 0x0F0F0F0F0F0F0F0FLLU)<< 4) | ((seq.high1 & 0xF0F0F0F0F0F0F0F0LLU)>> 4);
seq.high1 = ((seq.high1 & 0x00FF00FF00FF00FFLLU)<< 8) | ((seq.high1 & 0xFF00FF00FF00FF00LLU)>> 8);
seq.high1 = ((seq.high1 & 0x0000FFFF0000FFFFLLU)<<16) | ((seq.high1 & 0xFFFF0000FFFF0000LLU)>>16);
seq.high1 = ((seq.high1 & 0x00000000FFFFFFFFLLU)<<32) | ((seq.high1 & 0xFFFFFFFF00000000LLU)>>32);
ubyte8 temp_t;
temp = seq.high2;
seq.high2 = seq.low1;
seq.low1 = temp;
temp_t = seq.high1;
seq.high1 = seq.low2;
seq.low2 = temp_t;
seq = KmerRightBitMove(seq,256-(seq_size<<1));
return seq;*/


Kmer reverseComplementVerbose ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}

Kmer reverseComplement ( Kmer word, int overlap )
{
	return fastReverseComp ( word, overlap );
}
