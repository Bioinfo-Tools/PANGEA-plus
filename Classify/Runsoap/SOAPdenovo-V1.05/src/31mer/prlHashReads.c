/*
 * 31mer/prlHashReads.c
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


//debugging variables
static long long * tips;
static long long * kmerCounter;

static long long ** kmerFreq;

//buffer related varibles for chop kmer
static int read_c;
static char ** rcSeq;
static char ** seqBuffer;
static int * lenBuffer;
static int * indexArray;

//buffer related varibles for splay tree
static int buffer_size = 10000000;
static volatile int kmer_c;
static Kmer * kmerBuffer, *hashBanBuffer;
static char * nextcBuffer, *prevcBuffer;

static void thread_mark ( KmerSet * set, unsigned char thrdID );
static void Mark1in1outNode ( unsigned char * thrdSignal );
static void thread_delow ( KmerSet * set, unsigned char thrdID );
static void deLowCov ( unsigned char * thrdSignal );

static void singleKmer ( int t, KmerSet * kset );
static void chopKmer4read ( int t, int threadID );

static void freqStat ( char * outfile );

static void threadRoutine ( void * para )
{
	PARAMETER * prm;
	int i;
	unsigned char id;
	prm = ( PARAMETER * ) para;
	id = prm->threadID;

	//printf("%dth thread with threadID %d, hash_table %p\n",id,prm.threadID,prm.hash_table);
	while ( 1 )
	{
		if ( * ( prm->selfSignal ) == 1 )
		{
			for ( i = 0; i < kmer_c; i++ )
			{
				//if((unsigned char)(magic_seq(hashBanBuffer[i])%thrd_num)!=id)
				//if((kmerBuffer[i]%thrd_num)!=id)
				if ( ( hashBanBuffer[i] % thrd_num ) != id )
					{ continue; }

				kmerCounter[id + 1]++;
				singleKmer ( i, KmerSets[id] );
			}

			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 2 )
		{
			for ( i = 0; i < read_c; i++ )
			{
				if ( i % thrd_num != id )
					{ continue; }

				chopKmer4read ( i, id + 1 );
			}

			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 3 )
		{
			* ( prm->selfSignal ) = 0;
			break;
		}
		else if ( * ( prm->selfSignal ) == 4 )
		{
			thread_mark ( KmerSets[id], id );
			* ( prm->selfSignal ) = 0;
		}
		else if ( * ( prm->selfSignal ) == 5 )
		{
			thread_delow ( KmerSets[id], id );
			* ( prm->selfSignal ) = 0;
		}

		usleep ( 1 );
	}
}

static void singleKmer ( int t, KmerSet * kset )
{
	kmer_t * pos;
	put_kmerset ( kset, kmerBuffer[t], prevcBuffer[t], nextcBuffer[t], &pos );
}

static void creatThrds ( pthread_t * threads, PARAMETER * paras )
{
	unsigned char i;
	int temp;

	for ( i = 0; i < thrd_num; i++ )
	{
		//printf("to create %dth thread\n",(*(char *)&(threadID[i])));
		if ( ( temp = pthread_create ( &threads[i], NULL, ( void * ) threadRoutine, & ( paras[i] ) ) ) != 0 )
		{
			printf ( "create threads failed\n" );
			exit ( 1 );
		}
	}

	printf ( "%d thread created\n", thrd_num );
}

static void thread_wait ( pthread_t * threads )
{
	int i;

	for ( i = 0; i < thrd_num; i++ )
		if ( threads[i] != 0 )
			{ pthread_join ( threads[i], NULL ); }
}

static void chopKmer4read ( int t, int threadID )
{
	char * src_seq = seqBuffer[t];
	char * bal_seq = rcSeq[threadID];
	int len_seq = lenBuffer[t];
	int j, bal_j;
	Kmer hash_ban, bal_hash_ban;
	Kmer word, bal_word;
	int index;
	char InvalidCh = 4;
	word = 0;

	for ( index = 0; index < overlaplen; index++ )
	{
		word <<= 2;
		word += src_seq[index];
	}

	reverseComplementSeq ( src_seq, len_seq, bal_seq );
	// complementary node
	bal_word = reverseComplement ( word, overlaplen );
	bal_j = len_seq - 0 - overlaplen; //  0;
	index = indexArray[t];

	if ( word < bal_word )
	{
		hash_ban = hash_kmer ( word );
		hashBanBuffer[index] = hash_ban;
		kmerBuffer[index] = word;
		prevcBuffer[index] = InvalidCh;
		nextcBuffer[index++] = src_seq[0 + overlaplen];
	}
	else
	{
		bal_hash_ban = hash_kmer ( bal_word );
		hashBanBuffer[index] = bal_hash_ban;
		kmerBuffer[index] = bal_word;
		prevcBuffer[index] = bal_seq[bal_j - 1];
		nextcBuffer[index++] = InvalidCh;
	}

	for ( j = 1; j <= len_seq - overlaplen; j ++ )
	{
		word = nextKmer ( word, src_seq[j - 1 + overlaplen] );
		bal_j = len_seq - j - overlaplen; //  j;
		bal_word = prevKmer ( bal_word, bal_seq[bal_j] );

		if ( word < bal_word )
		{
			hash_ban = hash_kmer ( word );
			hashBanBuffer[index] = hash_ban;
			kmerBuffer[index] = word;
			prevcBuffer[index] = src_seq[j - 1];

			if ( j < len_seq - overlaplen )
				{ nextcBuffer[index++] = src_seq[j + overlaplen]; }
			else
				{ nextcBuffer[index++] = InvalidCh; }

			//printf("%dth: %p with %p\n",kmer_c-1,word,hashBanBuffer[kmer_c-1]);
		}
		else
		{
			// complementary node
			bal_hash_ban = hash_kmer ( bal_word );
			hashBanBuffer[index] = bal_hash_ban;
			kmerBuffer[index] = bal_word;

			if ( bal_j > 0 )
				{ prevcBuffer[index] = bal_seq[bal_j - 1]; }
			else
				{ prevcBuffer[index] = InvalidCh; }

			nextcBuffer[index++] = bal_seq[bal_j + overlaplen];
			//printf("%dth: %p with %p\n",kmer_c-1,bal_word,hashBanBuffer[kmer_c-1]);
		}
	}
}

static void sendWorkSignal ( unsigned char SIG, unsigned char * thrdSignals )
{
	int t;

	for ( t = 0; t < thrd_num; t++ )
		{ thrdSignals[t + 1] = SIG; }

	while ( 1 )
	{
		usleep ( 10 );

		for ( t = 0; t < thrd_num; t++ )
			if ( thrdSignals[t + 1] )
				{ break; }

		if ( t == thrd_num )
			{ break; }
	}
}

boolean prlRead2HashTable ( char * libfile, char * outfile )
{
	long long i;
	char * next_name, name[256];
	FILE * fo;
	time_t start_t, stop_t;
	int maxReadNum;
	int libNo;
	pthread_t threads[thrd_num];
	unsigned char thrdSignal[thrd_num + 1];
	PARAMETER paras[thrd_num];
	boolean flag, pairs = 0;
	WORDFILTER = ( ( ( Kmer ) 1 ) << ( 2 * overlaplen ) ) - 1;
	maxReadLen = 0;
	maxNameLen = 256;
	scan_libInfo ( libfile );
	alloc_pe_mem ( num_libs );

	if ( !maxReadLen )
		{ maxReadLen = 100; }

	maxReadLen4all = maxReadLen;
	printf ( "In %s, %d libs, max seq len %d, max name len %d\n\n",
	         libfile, num_libs, maxReadLen, maxNameLen );
	next_name = ( char * ) ckalloc ( ( maxNameLen + 1 ) * sizeof ( char ) );
	kmerBuffer = ( Kmer * ) ckalloc ( buffer_size * sizeof ( Kmer ) );
	hashBanBuffer = ( Kmer * ) ckalloc ( buffer_size * sizeof ( Kmer ) );
	prevcBuffer = ( char * ) ckalloc ( buffer_size * sizeof ( char ) );
	nextcBuffer = ( char * ) ckalloc ( buffer_size * sizeof ( char ) );
	maxReadNum = buffer_size / ( maxReadLen - overlaplen + 1 );
	//printf("buffer size %d, max read len %d, max read num %d\n",buffer_size,maxReadLen,maxReadNum);
	seqBuffer = ( char ** ) ckalloc ( maxReadNum * sizeof ( char * ) );
	lenBuffer = ( int * ) ckalloc ( maxReadNum * sizeof ( int ) );
	indexArray = ( int * ) ckalloc ( maxReadNum * sizeof ( int ) );

	for ( i = 0; i < maxReadNum; i++ )
		{ seqBuffer[i] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) ); }

	rcSeq = ( char ** ) ckalloc ( ( thrd_num + 1 ) * sizeof ( char * ) );

	if ( 1 )
	{
		kmerCounter = ( long long * ) ckalloc ( ( thrd_num + 1 ) * sizeof ( long long ) );
		KmerSets = ( KmerSet ** ) ckalloc ( thrd_num * sizeof ( KmerSet * ) );
		ubyte8 init_size = 1024;
		ubyte8 k = 0;

		if ( initKmerSetSize )
		{
			init_size = ( ubyte8 ) ( ( double ) initKmerSetSize * 1024.0f * 1024.0f * 1024.0f / ( double ) thrd_num / 16 );

			do
			{
				++k;
			}
			while ( k * 0xFFFFFFLLU < init_size );
		}

		for ( i = 0; i < thrd_num; i++ )
		{
			KmerSets[i] = init_kmerset ( ( ( initKmerSetSize ) ? ( k * 0xFFFFFFLLU ) : ( init_size ) ), 0.77f );
			ubyte8 tmp = ( initKmerSetSize ) ? ( k * 0xFFFFFFLLU ) : ( init_size );
			thrdSignal[i + 1] = 0;
			paras[i].threadID = i;
			paras[i].mainSignal = &thrdSignal[0];
			paras[i].selfSignal = &thrdSignal[i + 1];
			kmerCounter[i + 1] = 0;
			rcSeq[i + 1] = ( char * ) ckalloc ( maxReadLen * sizeof ( char ) );
		}

		creatThrds ( threads, paras );
	}

	thrdSignal[0] = kmerCounter[0] = 0;
	time ( &start_t );
	kmer_c = n_solexa = read_c = i = libNo = readNumBack = gradsCounter = 0;

	while ( ( flag = read1seqInLib ( seqBuffer[read_c], next_name, & ( lenBuffer[read_c] ), &libNo, pairs, 1 ) ) != 0 )
	{
		if ( ( ++i ) % 100000000 == 0 )
			{ printf ( "--- %lldth reads\n", i ); }

		if ( lenBuffer[read_c] < 0 )
			{ printf ( "read len %d\n", lenBuffer[read_c] ); }

		if ( lenBuffer[read_c] < overlaplen + 1 )
			{ continue; }

		/*
		if(lenBuffer[read_c]>70)
		    lenBuffer[read_c] = 50;
		else if(lenBuffer[read_c]>40)
		    lenBuffer[read_c] = 40;
		*/
		indexArray[read_c] = kmer_c;
		kmer_c += lenBuffer[read_c] - overlaplen + 1;
		read_c++;

		if ( read_c == maxReadNum )
		{
			kmerCounter[0] += kmer_c;
			sendWorkSignal ( 2, thrdSignal );
			sendWorkSignal ( 1, thrdSignal );
			kmer_c = read_c = 0;
		}
	}

	if ( read_c )
	{
		kmerCounter[0] += kmer_c;
		sendWorkSignal ( 2, thrdSignal );
		sendWorkSignal ( 1, thrdSignal );
	}

	time ( &stop_t );
	printf ( "time spent on hash reads: %ds, %lld reads processed\n", ( int ) ( stop_t - start_t ), i );

	//record insert size info
	if ( pairs )
	{
		if ( gradsCounter )
			printf ( "%d pe insert size, the largest boundary is %lld\n\n",
			         gradsCounter, pes[gradsCounter - 1].PE_bound );
		else
			{ printf ( "no paired reads found\n" ); }

		sprintf ( name, "%s.peGrads", outfile );
		fo = ckopen ( name, "w" );
		fprintf ( fo, "grads&num: %d\t%lld\n", gradsCounter, n_solexa );

		for ( i = 0; i < gradsCounter; i++ )
			{ fprintf ( fo, "%d\t%lld\t%d\n", pes[i].insertS, pes[i].PE_bound, pes[i].rank ); }

		fclose ( fo );
	}

	free_pe_mem();
	free_libs();

	if ( 1 )
	{
		unsigned long long alloCounter = 0;
		unsigned long long allKmerCounter = 0;

		for ( i = 0; i < thrd_num; i++ )
		{
			alloCounter += count_kmerset ( ( KmerSets[i] ) );
			allKmerCounter += kmerCounter[i + 1];
			free ( ( void * ) rcSeq[i + 1] );
		}

		printf ( "%lli nodes allocated, %lli kmer in reads, %lli kmer processed\n"
		         , alloCounter, kmerCounter[0], allKmerCounter );
		fflush ( stdout );
	}

	free ( ( void * ) rcSeq );
	free ( ( void * ) kmerCounter );

	for ( i = 0; i < maxReadNum; i++ )
		{ free ( ( void * ) seqBuffer[i] ); }

	free ( ( void * ) seqBuffer );
	free ( ( void * ) lenBuffer );
	free ( ( void * ) indexArray );
	free ( ( void * ) kmerBuffer );
	free ( ( void * ) hashBanBuffer );
	free ( ( void * ) nextcBuffer );
	free ( ( void * ) prevcBuffer );
	free ( ( void * ) next_name );

	//printf("done hashing nodes\n");
	if ( deLowKmer )
	{
		time ( &start_t );
		deLowCov ( thrdSignal );
		time ( &stop_t );
		printf ( "time spent on delowcvgNode %ds\n", ( int ) ( stop_t - start_t ) );
	}

	time ( &start_t );
	Mark1in1outNode ( thrdSignal );
	freqStat ( outfile );
	time ( &stop_t );
	printf ( "time spent on marking linear nodes %ds\n", ( int ) ( stop_t - start_t ) );
	fflush ( stdout );
	sendWorkSignal ( 3, thrdSignal );
	thread_wait ( threads );
	/*
	    Kmer word = 0x21c3ca82c734c8d0;
	    Kmer hash_ban = hash_kmer(word);
	    int setPicker = hash_ban%thrd_num;
	    kmer_t *node;
	    boolean found = search_kmerset(KmerSets[setPicker], word, &node);
	    if(!found)
	        printf("kmer %llx not found,\n",word);
	    else{
	        printf("kmer %llx, linear %d\n",word,node->linear);
	        for(i=0;i<4;i++){
	            if(get_kmer_right_cov(*node,i)>0)
	                printf("right %d, kmer %llx\n",i,nextKmer(node->seq,i));
	            if(get_kmer_left_cov(*node,i)>0)
	                printf("left %d, kmer %llx\n",i,prevKmer(node->seq,i));
	        }

	    }
	*/
	return 1;
}

static void thread_delow ( KmerSet * set, unsigned char thrdID )
{
	int i, in_num, out_num, cvgSingle;
	int l_cvg, r_cvg;
	kmer_t * rs;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			in_num = out_num = l_cvg = r_cvg = 0;
			rs = set->array + set->iter_ptr;

			for ( i = 0; i < 4; i++ )
			{
				cvgSingle = get_kmer_left_cov ( *rs, i );

				if ( cvgSingle > 0 && cvgSingle <= deLowKmer )
					{ set_kmer_left_cov ( *rs, i, 0 ); }

				cvgSingle = get_kmer_right_cov ( *rs, i );

				if ( cvgSingle > 0 && cvgSingle <= deLowKmer )
					{ set_kmer_right_cov ( *rs, i, 0 ); }
			}

			if ( rs->l_links == 0 && rs->r_links == 0 )
			{
				rs->deleted = 1;
				tips[thrdID]++;
			}
		}

		set->iter_ptr ++;
	}

	//printf("%lld single nodes, %lld linear\n",counter,tips[thrdID]);
}

static void deLowCov ( unsigned char * thrdSignal )
{
	int i;
	long long counter = 0;
	tips = ( long long * ) ckalloc ( thrd_num * sizeof ( long long ) );

	for ( i = 0; i < thrd_num; i++ )
		{ tips[i] = 0; }

	sendWorkSignal ( 5, thrdSignal ); //mark linear nodes

	for ( i = 0; i < thrd_num; i++ )
		{ counter += tips[i]; }

	free ( ( void * ) tips );
	printf ( "%lld kmer removed\n", counter );
}

static void printKmer ( Kmer kmer )
{
	int i;
	char kmerSeq[32], ch;

	for ( i = overlaplen - 1; i >= 0; i-- )
	{
		ch = kmer & 3;
		kmer >>= 2;
		kmerSeq[i] = ch;
	}

	for ( i = 0; i < overlaplen; i++ )
		{ printf ( "%c", int2base ( ( int ) kmerSeq[i] ) ); }

	printf ( "\n" );
}

static void thread_mark ( KmerSet * set, unsigned char thrdID )
{
	int i, in_num, out_num, cvgSingle;
	int l_cvg, r_cvg;
	kmer_t * rs;
	long long counter = 0;
	set->iter_ptr = 0;

	while ( set->iter_ptr < set->size )
	{
		if ( !is_kmer_entity_null ( set->flags, set->iter_ptr ) )
		{
			in_num = out_num = l_cvg = r_cvg = 0;
			rs = set->array + set->iter_ptr;

			for ( i = 0; i < 4; i++ )
			{
				cvgSingle = get_kmer_left_cov ( *rs, i );

				if ( cvgSingle > 0 )
				{
					in_num++;
					l_cvg += cvgSingle;
				}

				cvgSingle = get_kmer_right_cov ( *rs, i );

				if ( cvgSingle > 0 )
				{
					out_num++;
					r_cvg += cvgSingle;
				}
			}

			if ( rs->single )
			{
				kmerFreq[thrdID][1]++;
				counter++;
			}
			else
				{ kmerFreq[thrdID][ ( l_cvg > r_cvg ? l_cvg : r_cvg )]++; }

			if ( in_num == 1 && out_num == 1 )
			{
				rs->linear = 1;
				tips[thrdID]++;
			}
		}

		set->iter_ptr ++;
	}

	//printf("%lld single nodes, %lld linear\n",counter,tips[thrdID]);
}

static void Mark1in1outNode ( unsigned char * thrdSignal )
{
	int i;
	long long counter = 0;
	kmerFreq = ( long long ** ) ckalloc ( thrd_num * sizeof ( long long * ) );
	tips = ( long long * ) ckalloc ( thrd_num * sizeof ( long long ) );

	for ( i = 0; i < thrd_num; i++ )
	{
		kmerFreq[i] = ( long long * ) ckalloc ( 257 * sizeof ( long long ) );
		memset ( kmerFreq[i], 0, 257 * sizeof ( long long ) );
		tips[i] = 0;
	}

	sendWorkSignal ( 4, thrdSignal ); //mark linear nodes

	for ( i = 0; i < thrd_num; i++ )
		{ counter += tips[i]; }

	free ( ( void * ) tips );
	printf ( "%lld linear nodes\n", counter );
}

static void freqStat ( char * outfile )
{
	FILE * fo;
	char name[256];
	int i, j;
	long long sum;
	sprintf ( name, "%s.kmerFreq", outfile );
	fo = ckopen ( name, "w" );

	for ( i = 1; i < 256; i++ )
	{
		sum = 0;

		for ( j = 0; j < thrd_num; j++ )
			{ sum += kmerFreq[j][i]; }

		fprintf ( fo, "%lld\n", sum );
	}

	for ( i = 0; i < thrd_num; i++ )
		{ free ( ( void * ) kmerFreq[i] ); }

	free ( ( void * ) kmerFreq );
	fclose ( fo );
}


