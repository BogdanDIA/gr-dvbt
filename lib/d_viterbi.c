/*
 * Copyright 1995 Phil Karn, KA9Q
 * Copyright 2008 Free Software Foundation, Inc.
 *
 * This file is part of GNU Radio
 *
 * GNU Radio is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3, or (at your option)
 * any later version.
 *
 * GNU Radio is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GNU Radio; see the file COPYING.  If not, write to
 * the Free Software Foundation, Inc., 51 Franklin Street,
 * Boston, MA 02110-1301, USA.
 */

/*
 * Viterbi decoder for K=7 rate=1/2 convolutional code
 * Some modifications from original Karn code by Matt Ettus
 */

#include "d_viterbi.h"
#include <xmmintrin.h>

/* The two generator polynomials for the NASA Standard K=7 code.
 * Since these polynomials are known to be optimal for this constraint
 * length there is not much point in changing them. But if you do, you
 * will have to regenerate the BUTTERFLY macro calls in viterbi()
 */
#define	POLYA	0x4f
#define	POLYB	0x6d

/* The basic Viterbi decoder operation, called a "butterfly"
 * operation because of the way it looks on a trellis diagram. Each
 * butterfly involves an Add-Compare-Select (ACS) operation on the two nodes
 * where the 0 and 1 paths from the current node merge at the next step of
 * the trellis.
 *
 * The code polynomials are assumed to have 1's on both ends. Given a
 * function encode_state() that returns the two symbols for a given
 * encoder state in the low two bits, such a code will have the following
 * identities for even 'n' < 64:
 *
 * 	encode_state(n) = encode_state(n+65)
 *	encode_state(n+1) = encode_state(n+64) = (3 ^ encode_state(n))
 *
 * Any convolutional code you would actually want to use will have
 * these properties, so these assumptions aren't too limiting.
 *
 * Doing this as a macro lets the compiler evaluate at compile time the
 * many expressions that depend on the loop index and encoder state and
 * emit them as immediate arguments.
 * This makes an enormous difference on register-starved machines such
 * as the Intel x86 family where evaluating these expressions at runtime
 * would spill over into memory.
 */

union branchtab27 { unsigned char c[32]; __m128i v[2];} Branchtab27_sse2[2];

unsigned char mmresult[64] __attribute__((aligned(16))) = {0};
unsigned char ppresult[64] __attribute__((aligned(16))) = {0};

unsigned char mr0[16] __attribute__((aligned(16))) = {0};
unsigned char mr1[16] __attribute__((aligned(16))) = {0};
unsigned char mr2[16] __attribute__((aligned(16))) = {0};
unsigned char mr3[16] __attribute__((aligned(16))) = {0};

int start_offset;
int count = 0;


#if 0
  if (count > OFFSET) {\
	printf("gen - m0: %i, m1: %i, m2: %i, m3: %i\n", m0, m1, m2, m3); \
	printf("next[%i].metric: %i, next[%i].path: %x\n", 2*i, next[2*i].metric, 2*i, next[2*i].path); \
	printf("next[%i].metric: %i, next[%i].path: %x\n", 2*i+1, next[2*i+1].metric, 2*i+1, next[2*i+1].path); \
  }
#endif


#define OFFSET (1526746-10)*8

#define	BUTTERFLY(i,sym) { \
	int m0,m1,m2,m3;\
	/* ACS for 0 branch */\
	m0 = state[i].metric + mets[sym];	/* 2*i */\
	m1 = state[i+32].metric + mets[3^sym];	/* 2*i + 64 */\
	if(m0 > m1){\
		next[2*i].metric = m0;\
		next[2*i].path = state[i].path << 1;\
	} else {\
		next[2*i].metric = m1;\
		next[2*i].path = (state[i+32].path << 1)|1;\
	}\
	/* ACS for 1 branch */\
	m2 = state[i].metric + mets[3^sym];	/* 2*i + 1 */\
	m3 = state[i+32].metric + mets[sym];	/* 2*i + 65 */\
	if(m2 > m3){\
		next[2*i+1].metric = m2;\
		next[2*i+1].path = state[i].path << 1;\
	} else {\
		next[2*i+1].metric = m3;\
		next[2*i+1].path = (state[i+32].path << 1)|1;\
	}\
}

extern unsigned char d_Partab[];	/* Parity lookup table */

/* Convolutionally encode data into binary symbols */
unsigned char
d_encode(unsigned char *symbols,
       unsigned char *data,
       unsigned int nbytes,
       unsigned char encstate)
{
  int i;

  while(nbytes-- != 0){
    for(i=7;i>=0;i--){
      encstate = (encstate << 1) | ((*data >> i) & 1);
      *symbols++ = d_Partab[encstate & POLYA];
      *symbols++ = d_Partab[encstate & POLYB];
    }
    data++;
  }

  return encstate;
}

/* Viterbi decoder */
int
d_viterbi(unsigned long *metric,	/* Final path metric (returned value) */
	unsigned char *data,	/* Decoded output data */
	unsigned char *symbols,	/* Raw deinterleaved input symbols */
	unsigned int nbits,	/* Number of output bits */
	int mettab[2][256]	/* Metric table, [sent sym][rx symbol] */
	){
  unsigned int bitcnt = 0;
  int mets[4];
  long bestmetric;
  int beststate,i;
  struct viterbi_state state0[64],state1[64],*state,*next;

  state = state0;
  next = state1;

  /* Initialize starting metrics to prefer 0 state */
  state[0].metric = 0;
  for(i=1;i<64;i++)
    state[i].metric = -999999;
  state[0].path = 0;

  for(bitcnt = 0;bitcnt < nbits;bitcnt++){
    /* Read input symbol pair and compute all possible branch
     * metrics
     */
    mets[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
    mets[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
    mets[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
    mets[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
    symbols += 2;

#if 0
    /* These macro calls were generated by genbut.c */
    BUTTERFLY(0,0);
    BUTTERFLY(1,1);
    BUTTERFLY(2,3);
    BUTTERFLY(3,2);
    BUTTERFLY(4,3);
    BUTTERFLY(5,2);
    BUTTERFLY(6,0);
    BUTTERFLY(7,1);
    BUTTERFLY(8,0);
    BUTTERFLY(9,1);
    BUTTERFLY(10,3);
    BUTTERFLY(11,2);
    BUTTERFLY(12,3);
    BUTTERFLY(13,2);
    BUTTERFLY(14,0);
    BUTTERFLY(15,1);
    BUTTERFLY(16,2);
    BUTTERFLY(17,3);
    BUTTERFLY(18,1);
    BUTTERFLY(19,0);
    BUTTERFLY(20,1);
    BUTTERFLY(21,0);
    BUTTERFLY(22,2);
    BUTTERFLY(23,3);
    BUTTERFLY(24,2);
    BUTTERFLY(25,3);
    BUTTERFLY(26,1);
    BUTTERFLY(27,0);
    BUTTERFLY(28,1);
    BUTTERFLY(29,0);
    BUTTERFLY(30,2);
    BUTTERFLY(31,3);
#else

  	  BUTTERFLY(0,0)
	  BUTTERFLY(1,2)
	  BUTTERFLY(2,3)
	  BUTTERFLY(3,1)
	  BUTTERFLY(4,3)
	  BUTTERFLY(5,1)
	  BUTTERFLY(6,0)
	  BUTTERFLY(7,2)
	  BUTTERFLY(8,0)
	  BUTTERFLY(9,2)
	  BUTTERFLY(10,3)
	  BUTTERFLY(11,1)
	  BUTTERFLY(12,3)
	  BUTTERFLY(13,1)
	  BUTTERFLY(14,0)
	  BUTTERFLY(15,2)
	  BUTTERFLY(16,1)
	  BUTTERFLY(17,3)
	  BUTTERFLY(18,2)
	  BUTTERFLY(19,0)
	  BUTTERFLY(20,2)
	  BUTTERFLY(21,0)
	  BUTTERFLY(22,1)
	  BUTTERFLY(23,3)
	  BUTTERFLY(24,1)
	  BUTTERFLY(25,3)
	  BUTTERFLY(26,2)
	  BUTTERFLY(27,0)
	  BUTTERFLY(28,2)
	  BUTTERFLY(29,0)
	  BUTTERFLY(30,1)
	  BUTTERFLY(31,3)
#endif


    /* Swap current and next states */
    if(bitcnt & 1){
      state = state0;
      next = state1;
    } else {
      state = state1;
      next = state0;
    }
    // ETTUS
    //if(bitcnt > nbits-7){
    /* In tail, poison non-zero nodes */
    //for(i=1;i<64;i += 2)
    //	state[i].metric = -9999999;
    //}
    /* Produce output every 8 bits once path memory is full */
    if((bitcnt % 8) == 5 && bitcnt > 32){
      /* Find current best path */
      bestmetric = state[0].metric;
      beststate = 0;
      for(i=1;i<64;i++){
	if(state[i].metric > bestmetric){
	  bestmetric = state[i].metric;
	  beststate = i;
	}
      }
#ifdef	notdef
      printf("metrics[%d] = %d state = %lx\n",beststate,
	     state[beststate].metric,state[beststate].path);
#endif
      *data++ = state[beststate].path >> 24;
    }

  }
  /* Output remaining bits from 0 state */
  // ETTUS  Find best state instead
  bestmetric = state[0].metric;
  beststate = 0;
  for(i=1;i<64;i++){
    if(state[i].metric > bestmetric){
      bestmetric = state[i].metric;
      beststate = i;
    }
  }
  if((i = bitcnt % 8) != 6)
    state[beststate].path <<= 6-i;

  *data++ = state[beststate].path >> 24;
  *data++ = state[beststate].path >> 16;
  *data++ = state[beststate].path >> 8;
  *data = state[beststate].path;
  //printf ("BS = %d\tBSM = %d\tM0 = %d\n",beststate,state[beststate].metric,state[0].metric);
  *metric = state[beststate].metric;
  return 0;
}

void
d_viterbi_chunks_init(struct viterbi_state* state) {
  // Initialize starting metrics to prefer 0 state
  int i;
  state[0].metric = 0;
  state[0].path = 0;
  for(i=1;i<64;i++)
  {
    state[i].metric = 0;//-999999;
    state[i].path = 0;
  }
}

void
d_viterbi_chunks_init_sse2(__m128i *mm0, __m128i *pp0) {
  // Initialize starting metrics to prefer 0 state
  int i;

  for (i = 0; i < 4; i++)
  {
    mm0[i] = _mm_setzero_si128();
    pp0[i] = _mm_setzero_si128();
  }

  int polys[2] = { POLYA, POLYB };
  for(i=0; i < 32; i++)
  {
    Branchtab27_sse2[0].c[i] = (polys[0] < 0) ^ d_Partab[(2*i) & abs(polys[0])] ? 1 : 0;
    Branchtab27_sse2[1].c[i] = (polys[1] < 0) ^ d_Partab[(2*i) & abs(polys[1])] ? 1 : 0;
  }
}


void
d_viterbi_butterfly8(unsigned char *symbols, int mettab[2][256], struct viterbi_state *state0, struct viterbi_state *state1)
{
  unsigned int bitcnt;
  int mets[4];

  struct viterbi_state *state, *next;
  state = state0;
  next = state1;
  // Operate on 16 symbols (8 bits) at a time
  for(bitcnt = 0;bitcnt < 8;bitcnt++){
    // Read input symbol pair and compute all possible branch metrics
    mets[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
    mets[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
    mets[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
    mets[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
    symbols += 2;

    // These macro calls were generated by genbut.c
    BUTTERFLY(0,0);BUTTERFLY(1,1);BUTTERFLY(2,3);BUTTERFLY(3,2);
    BUTTERFLY(4,3);BUTTERFLY(5,2);BUTTERFLY(6,0);BUTTERFLY(7,1);
    BUTTERFLY(8,0);BUTTERFLY(9,1);BUTTERFLY(10,3);BUTTERFLY(11,2);
    BUTTERFLY(12,3);BUTTERFLY(13,2);BUTTERFLY(14,0);BUTTERFLY(15,1);
    BUTTERFLY(16,2);BUTTERFLY(17,3);BUTTERFLY(18,1);BUTTERFLY(19,0);
    BUTTERFLY(20,1);BUTTERFLY(21,0);BUTTERFLY(22,2);BUTTERFLY(23,3);
    BUTTERFLY(24,2);BUTTERFLY(25,3);BUTTERFLY(26,1);BUTTERFLY(27,0);
    BUTTERFLY(28,1);BUTTERFLY(29,0);BUTTERFLY(30,2);BUTTERFLY(31,3);

    // Swap current and next states
    if(bitcnt & 1){
      state = state0;
      next = state1;
    } else {
      state = state1;
      next = state0;
    }
  }
}

void
d_viterbi_butterfly2(unsigned char *symbols, int mettab[2][256], struct viterbi_state *state0, struct viterbi_state *state1)
{
  int i, j;
  int mets[4];

  struct viterbi_state *state, *next;
  state = state0;
  next = state1;

  // Operate on 4 symbols (2 bits) at a time

  // Read input symbol pair and compute all possible branch metrics
  mets[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
  mets[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
  mets[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
  mets[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];

#if 0
  if (count > OFFSET)
  {
  for (i = 0; i < 4; i++)
    printf("mets[%i]: %i\n", i, mets[i]);
  }
#endif

  // These macro calls were generated by genbut.c
  BUTTERFLY(0,0)
  BUTTERFLY(1,2)
  BUTTERFLY(2,3)
  BUTTERFLY(3,1)
  BUTTERFLY(4,3)
  BUTTERFLY(5,1)
  BUTTERFLY(6,0)
  BUTTERFLY(7,2)
  BUTTERFLY(8,0)
  BUTTERFLY(9,2)
  BUTTERFLY(10,3)
  BUTTERFLY(11,1)
  BUTTERFLY(12,3)
  BUTTERFLY(13,1)
  BUTTERFLY(14,0)
  BUTTERFLY(15,2)
  BUTTERFLY(16,1)
  BUTTERFLY(17,3)
  BUTTERFLY(18,2)
  BUTTERFLY(19,0)
  BUTTERFLY(20,2)
  BUTTERFLY(21,0)
  BUTTERFLY(22,1)
  BUTTERFLY(23,3)
  BUTTERFLY(24,1)
  BUTTERFLY(25,3)
  BUTTERFLY(26,2)
  BUTTERFLY(27,0)
  BUTTERFLY(28,2)
  BUTTERFLY(29,0)
  BUTTERFLY(30,1)
  BUTTERFLY(31,3)

  state = state1;
  next = state0;

  // Read input symbol pair and compute all possible branch metrics
  mets[0] = mettab[0][symbols[2]] + mettab[0][symbols[3]];
  mets[1] = mettab[0][symbols[2]] + mettab[1][symbols[3]];
  mets[2] = mettab[1][symbols[2]] + mettab[0][symbols[3]];
  mets[3] = mettab[1][symbols[2]] + mettab[1][symbols[3]];

  // These macro calls were generated by genbut.c
  BUTTERFLY(0,0)
  BUTTERFLY(1,2)
  BUTTERFLY(2,3)
  BUTTERFLY(3,1)
  BUTTERFLY(4,3)
  BUTTERFLY(5,1)
  BUTTERFLY(6,0)
  BUTTERFLY(7,2)
  BUTTERFLY(8,0)
  BUTTERFLY(9,2)
  BUTTERFLY(10,3)
  BUTTERFLY(11,1)
  BUTTERFLY(12,3)
  BUTTERFLY(13,1)
  BUTTERFLY(14,0)
  BUTTERFLY(15,2)
  BUTTERFLY(16,1)
  BUTTERFLY(17,3)
  BUTTERFLY(18,2)
  BUTTERFLY(19,0)
  BUTTERFLY(20,2)
  BUTTERFLY(21,0)
  BUTTERFLY(22,1)
  BUTTERFLY(23,3)
  BUTTERFLY(24,1)
  BUTTERFLY(25,3)
  BUTTERFLY(26,2)
  BUTTERFLY(27,0)
  BUTTERFLY(28,2)
  BUTTERFLY(29,0)
  BUTTERFLY(30,1)
  BUTTERFLY(31,3)
}

void
d_viterbi_butterfly2_sse2(unsigned char *symbols, __m128i *mm0, __m128i *mm1, __m128i *pp0, __m128i *pp1)
{
  int i, j;

  __m128i *metric0, *metric1;
  __m128i *path0, *path1;

  metric0 = mm0;
  path0 = pp0;
  metric1 = mm1;
  path1 = pp1;

  // Operate on 4 symbols (2 bits) at a time

  __m128i m0, m1, m2, m3, decision0, decision1, survivor0, survivor1;
  __m128i metsv, metsvm;
  __m128i shift0, shift1;
  __m128i tmp0, tmp1;
  __m128i sym0v, sym1v;
  __m128i min;

  for (i = 0; i < 2; i++)
  {
    sym0v = _mm_set1_epi8(symbols[0]);
    sym1v = _mm_set1_epi8(symbols[1]);

    metsvm = _mm_add_epi8(_mm_xor_si128(Branchtab27_sse2[0].v[i],sym0v),_mm_xor_si128(Branchtab27_sse2[1].v[i],sym1v));
    metsv = _mm_sub_epi8(_mm_set1_epi8(2),metsvm);

    m0 = _mm_add_epi8(metric0[i], metsv);
    m1 = _mm_add_epi8(metric0[i+2], metsvm);
    m2 = _mm_add_epi8(metric0[i], metsvm);
    m3 = _mm_add_epi8(metric0[i+2], metsv);

#if 0
    _mm_store_si128((__m128i *) &mr0, m0);
    _mm_store_si128((__m128i *) &mr1, m1);
    _mm_store_si128((__m128i *) &mr2, m2);
    _mm_store_si128((__m128i *) &mr3, m3);

    if (count > OFFSET)
    {
    for (j = 0; j < 16; j++)
      printf("count: %i, m0: %i, m1: %i, m2: %i, m3: %i\n", count, mr0[j], mr1[j], mr2[j], mr3[j]);
    }
#endif

    decision0 = _mm_cmpgt_epi8(_mm_sub_epi8(m0,m1),_mm_setzero_si128());
    decision1 = _mm_cmpgt_epi8(_mm_sub_epi8(m2,m3),_mm_setzero_si128());
    survivor0 = _mm_or_si128(_mm_and_si128(decision0,m0),_mm_andnot_si128(decision0,m1));
    survivor1 = _mm_or_si128(_mm_and_si128(decision1,m2),_mm_andnot_si128(decision1,m3));

    shift0 = _mm_slli_epi16(path0[i], 1);
    shift1 = _mm_slli_epi16(path0[2+i], 1);
    shift1 = _mm_add_epi8(shift1, _mm_set1_epi8(1));

    metric1[2*i] = _mm_unpacklo_epi8(survivor0,survivor1);
    tmp0 = _mm_or_si128(_mm_and_si128(decision0,shift0),_mm_andnot_si128(decision0,shift1));

    metric1[2*i+1] = _mm_unpackhi_epi8(survivor0,survivor1);
    tmp1 = _mm_or_si128(_mm_and_si128(decision1,shift0),_mm_andnot_si128(decision1,shift1));

    path1[2*i] = _mm_unpacklo_epi8(tmp0, tmp1);
    path1[2*i+1] = _mm_unpackhi_epi8(tmp0, tmp1);


#if 0
    if (count > OFFSET)
    {
    _mm_store_si128((__m128i *) &mr0, metric1[2*i]);
    _mm_store_si128((__m128i *) &mr1, path1[2*i]);
    _mm_store_si128((__m128i *) &mr2, metric1[2*i+1]);
    _mm_store_si128((__m128i *) &mr3, path1[2*i+1]);
    for (j = 0; j < 16; j++)
      printf("count: %i, metric1: %x, path1: %x\n", count, mr0[j], mr1[j]);
    for (j = 0; j < 16; j++)
      printf("count: %i, metric1: %x, path1: %x\n", count, mr2[j], mr3[j]);
    }
#endif
  }

  metric0 = mm1;
  path0 = pp1;
  metric1 = mm0;
  path1 = pp0;

  count++;

  for (i = 0; i < 2; i++)
  {
    sym0v = _mm_set1_epi8(symbols[2]);
    sym1v = _mm_set1_epi8(symbols[3]);

    metsvm = _mm_add_epi8(_mm_xor_si128(Branchtab27_sse2[0].v[i],sym0v),_mm_xor_si128(Branchtab27_sse2[1].v[i],sym1v));
    metsv = _mm_sub_epi8(_mm_set1_epi8(2),metsvm);

    m0 = _mm_add_epi8(metric0[i], metsv);
    m1 = _mm_add_epi8(metric0[i+2], metsvm);
    m2 = _mm_add_epi8(metric0[i], metsvm);
    m3 = _mm_add_epi8(metric0[i+2], metsv);

#if 0
    _mm_store_si128((__m128i *) &mr0, m0);
    _mm_store_si128((__m128i *) &mr1, m1);
    _mm_store_si128((__m128i *) &mr2, m2);
    _mm_store_si128((__m128i *) &mr3, m3);

    if (count > OFFSET)
    {
    for (j = 0; j < 16; j++)
      printf("count: %i, m0: %i, m1: %i, m2: %i, m3: %i\n", count, mr0[j], mr1[j], mr2[j], mr3[j]);
    }
#endif

    decision0 = _mm_cmpgt_epi8(_mm_sub_epi8(m0,m1),_mm_setzero_si128());
    decision1 = _mm_cmpgt_epi8(_mm_sub_epi8(m2,m3),_mm_setzero_si128());
    survivor0 = _mm_or_si128(_mm_and_si128(decision0,m0),_mm_andnot_si128(decision0,m1));
    survivor1 = _mm_or_si128(_mm_and_si128(decision1,m2),_mm_andnot_si128(decision1,m3));

    shift0 = _mm_slli_epi16(path0[i], 1);
    shift1 = _mm_slli_epi16(path0[2+i], 1);
    shift1 = _mm_add_epi8(shift1, _mm_set1_epi8(1));

    metric1[2*i] = _mm_unpacklo_epi8(survivor0,survivor1);
    tmp0 = _mm_or_si128(_mm_and_si128(decision0,shift0),_mm_andnot_si128(decision0,shift1));

    metric1[2*i+1] = _mm_unpackhi_epi8(survivor0,survivor1);
    tmp1 = _mm_or_si128(_mm_and_si128(decision1,shift0),_mm_andnot_si128(decision1,shift1));

    path1[2*i] = _mm_unpacklo_epi8(tmp0, tmp1);
    path1[2*i+1] = _mm_unpackhi_epi8(tmp0, tmp1);

#if 0
    if (count > OFFSET)
    {
    _mm_store_si128((__m128i *) &mr0, metric1[2*i]);
    _mm_store_si128((__m128i *) &mr1, path1[2*i]);
    _mm_store_si128((__m128i *) &mr2, metric1[2*i+1]);
    _mm_store_si128((__m128i *) &mr3, path1[2*i+1]);
    for (j = 0; j < 16; j++)
      printf("count: %i, metric1: %x, path1: %x\n", count, mr0[j], mr1[j]);
    for (j = 0; j < 16; j++)
      printf("count: %i, metric1: %x, path1: %x\n", count, mr2[j], mr3[j]);
    }
#endif
  }

  count++;
}


unsigned char
d_viterbi_get_output(struct viterbi_state *state, unsigned char *outbuf) {
  // Produce output every 8 bits once path memory is full
  //  if((bitcnt % 8) == 5 && bitcnt > 32) {

  //  Find current best path
  unsigned int i, j, beststate = 0;
  int bestmetric;

  bestmetric = state[0].metric;
  beststate = 0;
  for(i=1;i<64;i++)
  {
    if(state[i].metric > bestmetric) {
      bestmetric = state[i].metric;
      beststate = i;
    }
#if 0
    if (count > OFFSET)
      printf("state[%i]: %i, path: %x\n", i, state[i].metric, state[i].path);
#endif
  }

#if 0
  if (count > OFFSET)
    printf("d_viterbi_get_output: %i, %x\n", beststate, state[beststate].path >> 24);
#endif

  *outbuf = state[beststate].path >> 24;

  return bestmetric;
}

unsigned char
d_viterbi_get_output_sse2(__m128i *mm0, __m128i *pp0, unsigned char *outbuf) {
  // Produce output every 8 bits once path memory is full
  //  if((bitcnt % 8) == 5 && bitcnt > 32) {

  //  Find current best path
  unsigned int i, j, beststate;
  int bestmetric, minmetric;

  static unsigned int shift_reg = 0;
  static unsigned int cc = 0;

  // TODO - find anither way to extract the value
  for (i = 0; i < 4; i++)
  {
    _mm_store_si128((__m128i *) &mmresult[i*16], mm0[i]);
    _mm_store_si128((__m128i *) &ppresult[i*16], pp0[i]);
  }

  beststate = 0;

  bestmetric = mmresult[beststate];
  minmetric = mmresult[beststate];

  for (i = 1; i < 64; i++)
  {
#if 0
    if (count > OFFSET)
      printf("metric[%i]: %i, path: %x\n", i, mmresult[i], ppresult[i]);
#endif

    if (mmresult[i] > bestmetric)
    {
      bestmetric = mmresult[i];
      beststate = i;
    }
    if (mmresult[i] < minmetric)
      minmetric = mmresult[i];
  }


#if 1
  // Offset the initial sequence with 3 bytes
  shift_reg = (shift_reg << 8) | ppresult[beststate];
  *outbuf = (shift_reg >> 24) & 0xff;
#else
  *outbuf = ppresult[beststate];
#endif

  // Zero out the path variable
  for (i = 0; i < 4; i++)
  {
    pp0[i] = _mm_setzero_si128();
    mm0[i] = _mm_sub_epi8(mm0[i], _mm_set1_epi8(minmetric));
  }

#if 0
  if (count > OFFSET)
    printf("d_viterbi_get_output_sse2: %i, beststate: %i, %x\n", cc, beststate, *outbuf);

  cc++;
#endif

  return bestmetric;
}




//printf ("BS = %d\tBSM = %d\tM0 = %d\n",beststate,state[beststate].metric,state[0].metric);
// In tail, poison non-zero nodes
//if(bits_out > packet_size-7)
//  for(i=1;i<64;i += 2)
//    state[i].metric = -9999999;

