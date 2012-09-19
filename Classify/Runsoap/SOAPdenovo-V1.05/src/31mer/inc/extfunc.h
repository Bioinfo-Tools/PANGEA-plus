/*
 * 31mer/inc/extfunc.h
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

#include "check.h"
#include "extfunc2.h"
extern NODE ** seq2nodes_with_pair ( char * seqfile, char * outfile );
extern NODE ** prlSeq2nodes_with_pair ( char * seqfile, char * outfile );
extern void readseq1by1 ( char * src_seq, char * src_name, int * len_seq, FILE * fp, long long num_seq );
extern void readseqPbyP ( char * src_seq, char * src_name, int * insertS, int * len_seq, FILE * fp, long long num_seq );
extern void nodes2edges_with_pair ( NODE ** hash_table, EDGE_PT ** edge_list, char * outfile );
extern int findOrInsertOccurenceInNodeTree ( Kmer kmer, NODE ** T );
extern NODE * SplayNodeTree ( NODE * T, Kmer kmer );
extern Kmer reverseComplement ( Kmer word, int overlap );
extern Kmer hash_kmer ( Kmer kmer );
extern void link2next ( NODE * node, char ch );
extern unsigned char check_link2next ( NODE * node, char ch );
extern void unlink2next ( NODE * node, char ch );
extern void link2prev ( NODE * node, char ch );
extern unsigned char check_link2prev ( NODE * node, char ch );
extern void unlink2prev ( NODE * node, char ch );
extern int count_link2next ( NODE * node );
extern int count_link2prev ( NODE * node );
extern Kmer nextKmer ( Kmer prev, char ch );
extern Kmer prevKmer ( Kmer next, char ch );
extern long long readseqpar ( int * max_len, int * min_leg, int * max_name_len, FILE * fp );
extern void destroyNodeHash ( NODE ** hash_table );
extern void free_edge_list ( EDGE_PT * el );
extern void reverseComplementSeq ( char * seq, int len, char * bal_seq );
extern void free_node_list ( NODE_PT * np );
extern NODE * SplayNodeTree_FILTER ( NODE * T, Kmer kmer );
extern NODE * allocateNode_cvg ( Kmer kmer );
extern int findOrInsertOccurenceInNodeTree_cvg ( Kmer kmer, NODE ** T );
extern void free_edge_array ( EDGE * ed_array, int ed_num );
extern void free_lightctg_array ( LIGHTCTG * ed_array, int ed_num );
extern char getCharInTightString ( char * tightSeq, int pos );
extern void writeChar2tightSting ( char nt, char * tightSeq, int pos );
extern void short_reads_sum();
extern void read_one_sequence ( FILE * fp, long long * T, char ** X );
extern void output_edges ( preEDGE * ed_array, int ed_num, char * outfile );
extern void read2edge ( char * seqfile, NODE ** hash_table, char * outfile );
extern void loadVertex ( char * graphfile );
extern int kmer2vt ( Kmer kmer );
extern void loadEdge ( char * graphfile );
extern boolean loadPath ( char * graphfile );
extern READINTERVAL * allocateRV ( int readid, int edgeid );
extern void createRVmemo();
extern void dismissRV ( READINTERVAL * rv );
extern void destroyReadIntervMem();
extern void destroyConnectMem();
extern void u2uConcatenate();
extern void unlink2all ( NODE * node, NODE ** hash_table );
extern void cutTip ( NODE ** hash_table );
extern void output_contig ( EDGE * ed_array, unsigned int ed_num, char * outfile, int cut_len );
extern void printTightString ( char * tightSeq, int len );
extern int roughUniqueness ( unsigned int edgeno, char ignore_cvg, char * ignored );
extern void outputReadPos ( char * graphfile, int min_len );
extern NODE * reverseComplementNode ( NODE * node1, NODE ** hash_table );
extern void testSearch();
extern void print_kmer ( FILE * fp, Kmer kmer, char c );
extern void allpathConcatenate();
extern void output_updated_edges ( char * outfile );
extern void output_updated_vertex ( char * outfile );
extern void loadUpdatedEdges ( char * graphfile );
extern void loadUpdatedVertex ( char * graphfile );
extern void connectByPE ( char * infile );
extern void output_cntGVZ ( char * outfile );
extern void output_graph ( char * outfile );
extern void removeUnreliable ( NODE ** hash_talbe );
extern void testLinearC2C();
extern void output_contig_graph ( char * outfile );
extern void scaffolding ( unsigned int cut_len, char * outfile );
extern int cmp_int ( const void * a, const void * b );
extern CONNECT * allocateCN ( unsigned int contigId, int gap );
extern int recoverRep();
extern void loadPEgrads ( char * infile );
extern int putInsertS ( long long readid, int size, int * currGrads );
extern int getInsertS ( long long readid, int * readlen );
extern int connectByPE_grad ( FILE * fp, int peGrad, char * line );
extern void PEgradsScaf ( char * infile );
extern void reorderAnnotation ( char * infile, char * outfile );
extern int count_ends ( NODE ** hash_table );
extern void output_1edge ( preEDGE * edge, FILE * fp );
extern void prlRead2edge ( char * libfile, char * outfile );
extern int count_edges ( NODE ** hash_table );
extern int prlFindOrInsertOccurenceInNodeTree_cvg ( Kmer kmer, NODE ** T, MEM_MANAGER * node_mem_manager );
extern void prlDestroyNodeHash ( NODE ** hash_table );
extern void annotFileTrans ( char * infile, char * outfile );
extern void prlLoadPath ( char * graphfile );
extern void misCheck ( char * infile, char * outfile );
extern int uniqueLenSearch ( unsigned int * len_array, unsigned int * flag_array, int num, unsigned int target );
extern int cmp_vertex ( const void * a, const void * b );
extern void linkContig2Vts();
extern int bisearch ( VERTEX * vts, int num, Kmer target );
extern int connectByPE_gradPatch ( FILE * fp1, FILE * fp2, int peGrad, char * line1, char * line2 );
extern void scaftiging ( char * graphfile, int len_cut );
extern void gapFilling ( char * graphfile, int cut_len );
extern ARC * getArcBetween ( unsigned int from_ed, unsigned int to_ed );
extern void bubblePinch ( double simiCutoff, char * outfile, int M );
extern void linearConcatenate();
extern unsigned char setArcMulti ( unsigned int from_ed, unsigned int to_ed, unsigned char value );
extern ARC * allocateArc ( unsigned int edgeid );
extern void cutTipsInGraph ( int cutLen, boolean strict );
extern ARC * deleteArc ( ARC * arc_list, ARC * arc );
extern void compactEdgeArray();
extern void dismissArc ( ARC * arc );
extern void createArcMemo();
extern ARC * getArcBetween ( unsigned int from_ed, unsigned int to_ed );
extern ARC * allocateArc ( unsigned int edgeid );
extern void unlink2prevUncertain ( NODE * node, char ch, boolean smaller );
extern char firstCharInKmer ( Kmer kmer );
extern void writeChar2tightString ( char nt, char * tightSeq, int pos );
extern Kmer reverseComplementVerbose ( Kmer word, int overlap );
extern Kmer KmerPlus ( Kmer prev, char ch );
extern void output_heavyArcs ( char * outfile );
extern preARC * allocatePreArc ( unsigned int edgeid );
extern void destroyPreArcMem();
extern void traceAlongArc ( unsigned int destE, unsigned int currE, int max_steps, int min, int max, int index, int len, int * num_route );
extern void freeContig_array();
extern void output_scafSeq ( char * graphfile, int len_cut );
extern void putArcInHash ( unsigned int from_ed, unsigned int to_ed );
extern boolean DoesArcExist ( unsigned int from_ed, unsigned int to_ed );
extern void recordArcInHash();
extern void destroyArcHash();
extern void removeWeakEdges ( int lenCutoff, unsigned int multiCutoff );
extern void createArcLookupTable();
extern void deleteArcLookupTable();
extern void putArc2LookupTable ( unsigned int from_ed, ARC * arc );
extern void removeArcInLookupTable ( unsigned int from_ed, unsigned int to_ed );
extern ARC * arcCount ( unsigned int edgeid, unsigned int * num );
extern void mapFileTrans ( char * infile );
extern void solveReps();
extern void removeDeadArcs();
extern void destroyArcMem();
extern int count_link2prevB ( NODE * node );
extern int count_link2nextB ( NODE * node );
extern void getCntsInFile ( char * infile );
extern void scafByCntInfo ( char * infile );
extern CONNECT * add1Connect ( unsigned int e1, unsigned int e2, int gap, int weight, boolean inherit );
extern void getScaff ( char * infile );
extern void traceAlongMaskedCnt ( unsigned int destE, unsigned int currE, int max_steps, int min, int max, int index, int len, int * num_route );
extern void createPreArcMemManager();
extern boolean loadPathBin ( char * graphfile );
extern void analyzeTips ( NODE ** hash_table, char * graphfile );
extern void recordArcsInLookupTable();
extern FILE * multiFileRead1seq ( char * src_seq, char * src_name, int * len_seq, FILE * fp, FILE * freads );
extern void multiFileSeqpar ( FILE * fp );
extern long long multiFileParse ( int * max_leg, int * min_leg, int * max_name_leg, FILE * fp );
extern CONNECT * getCntBetween ( unsigned int from_ed, unsigned int to_ed );
extern void createCntMemManager();
extern void destroyConnectMem();
extern void createCntLookupTable();
extern void deleteCntLookupTable();
extern void putCnt2LookupTable ( unsigned int from_c, CONNECT * cnt );
extern int prlFindOrInsertOccurenceInEdonTree ( Kmer kmer, EDON ** T, MEM_MANAGER * node_mem_manager );
extern EDON * SplayEdonTree ( EDON * T, Kmer kmer );
extern void prlDestroyEdonHash ( EDON ** hash_table );
extern void prlRead2Ctg ( char * seqfile, char * outfile );
extern void prlLongRead2Ctg ( char * libfile, char * outfile );
extern boolean prlContig2nodes ( char * grapfile, int len_cut );
extern void scan_libInfo ( char * libfile );
extern int getMaxLongReadLen ( int num_libs );
extern void free_libs();
extern boolean read1seqInLib ( char * src_seq, char * src_name, int * len_seq,
                               int * libNo, boolean pair, unsigned char purpose );
extern NODE ** prlEdge2nodes ( char * grapfile );
extern void prlRead2graph ( char * libfile, NODE ** hash_table, char * outfile );
extern void save4laterSolve();
extern void solveRepsAfter();
extern void free_pe_mem();
extern void alloc_pe_mem ( int gradsCounter );
extern NODE * searchNodeTree ( NODE * T, Kmer kmer );
extern EDON * searchEdonTree ( EDON * T, Kmer kmer );
extern void prlDestroyPreArcMem();
extern preARC * prlAllocatePreArc ( unsigned int edgeid, MEM_MANAGER * manager );
extern boolean prlRead2HashTable ( char * libfile, char * outfile );
extern void free_allSets();
extern void removeSingleTips();
extern void removeMinorTips();
extern void kmer2edges ( char * outfile );
extern void output_vertex ( char * outfile );
extern boolean prlRead2HashTable ( char * libfile, char * outfile );
extern void Links2Scaf ( char * infile );
extern void PE2Links ( char * infile );
extern void basicContigInfo ( char * infile );
extern unsigned int getTwinCtg ( unsigned int ctg );
extern boolean isSmallerThanTwin ( unsigned int ctg );
extern boolean isLargerThanTwin ( unsigned int ctg );
extern boolean isSameAsTwin ( unsigned int ctg );
extern boolean loadMarkerBin ( char * graphfile );
extern void readsCloseGap ( char * graphfile );
extern void prlReadsCloseGap ( char * graphfile );
extern void locateReadOnScaf ( char * graphfile );
extern unsigned int getTwinEdge ( unsigned int edge );
extern boolean EdSmallerThanTwin ( unsigned int edge );
extern boolean EdLargerThanTwin ( unsigned int edge );
extern boolean EdSameAsTwin ( unsigned int edge );
extern void removeLowCovEdges ( int lenCutoff, unsigned short covCutoff );
extern int localGraph ( READNEARBY * rdArray, int num, CTGinSCAF * ctg1, CTGinSCAF * ctg2,
                        int origOverlap, Kmer * kmerCtg1, Kmer * kmerCtg2,
                        int overlap, DARRAY * gapSeqArray, char * seqCtg1, char * seqCtg2, char * seqGap );


