/*
 * Copyright (C) 2018, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori.
 * 
 * This file is part of the DDalphaAMG solver library.
 * 
 * The DDalphaAMG solver library is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * The DDalphaAMG solver library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * 
 * You should have received a copy of the GNU General Public License
 * along with the DDalphaAMG solver library. If not, see http://www.gnu.org/licenses/.
 * 
 */

#ifndef K_VECTOR_H
#define K_VECTOR_H

#ifdef K_OPT
//#include "hpc_emmintrin.h"
#include "../../other_libs/K/hpc_emmintrin.h"

typedef struct { _fjsp_v2r4 u[6];  } ssu3colvec;
typedef struct { _fjsp_v2r4 u[9];  } ssu3gauge;
typedef struct { _fjsp_v2r4 y[12]; } ssu3ferm;

typedef struct { _fjsp_v2r4 y[6]; } ssu3ferm2;
typedef struct { _fjsp_v2r4 c[18]; } ssu3clov2;

typedef struct { _fjsp_v2r4 x[2]; } block2vec;
typedef struct { _fjsp_v2r4 x[4]; } block4vec;
typedef struct { _fjsp_v2r4 x[8]; } block4diag;
typedef struct { _fjsp_v2r4 x[16]; } block4mat;


#define load_block4vec( d, s ){	\
  d ## 0 =  _fjsp_stod_v2r8(s->x[0]);	\
  d ## 1 =  _fjsp_stod_v2r8(s->x[1]);	\
  d ## 2 =  _fjsp_stod_v2r8(s->x[2]);	\
  d ## 3 =  _fjsp_stod_v2r8(s->x[3]);	\
}

#define store_block4vec( s, d ){	\
    s->x[0] = _fjsp_dtos_v2r4(d ## 0 );	\
    s->x[1] = _fjsp_dtos_v2r4(d ## 1 );	\
    s->x[2] = _fjsp_dtos_v2r4(d ## 2 );	\
    s->x[3] = _fjsp_dtos_v2r4(d ## 3 );	\
  }

#define load_block2vec( d, s ){	\
  d ## 0 =  _fjsp_stod_v2r8(s->x[0]);	\
  d ## 1 =  _fjsp_stod_v2r8(s->x[1]);	\
}

#define store_block2vec( s, d ){	\
    s->x[0] = _fjsp_dtos_v2r4(d ## 0 );	\
    s->x[1] = _fjsp_dtos_v2r4(d ## 1 );	\
  }

#define load_block4diag( d, s ){	\
  d ## 0  =  _fjsp_stod_v2r8(s->x[0]);	\
  d ## 1  =  _fjsp_stod_v2r8(s->x[1]);	\
  d ## 01 =  _fjsp_stod_v2r8(s->x[2]);	\
  d ## 02 =  _fjsp_stod_v2r8(s->x[3]);	\
  d ## 03 =  _fjsp_stod_v2r8(s->x[4]);	\
  d ## 12 =  _fjsp_stod_v2r8(s->x[5]);	\
  d ## 13 =  _fjsp_stod_v2r8(s->x[6]);	\
  d ## 23 =  _fjsp_stod_v2r8(s->x[7]);	\
}

#define load_block4mat( d, s ){	\
  d ## 00 =  _fjsp_stod_v2r8(s->x[0]);	\
  d ## 01 =  _fjsp_stod_v2r8(s->x[1]);	\
  d ## 02 =  _fjsp_stod_v2r8(s->x[2]);	\
  d ## 03 =  _fjsp_stod_v2r8(s->x[3]);	\
  d ## 10 =  _fjsp_stod_v2r8(s->x[4]);	\
  d ## 11 =  _fjsp_stod_v2r8(s->x[5]);	\
  d ## 12 =  _fjsp_stod_v2r8(s->x[6]);	\
  d ## 13 =  _fjsp_stod_v2r8(s->x[7]);	\
  d ## 20 =  _fjsp_stod_v2r8(s->x[8]);	\
  d ## 21 =  _fjsp_stod_v2r8(s->x[9]);	\
  d ## 22 =  _fjsp_stod_v2r8(s->x[10]);	\
  d ## 23 =  _fjsp_stod_v2r8(s->x[11]);	\
  d ## 30 =  _fjsp_stod_v2r8(s->x[12]);	\
  d ## 31 =  _fjsp_stod_v2r8(s->x[13]);	\
  d ## 32 =  _fjsp_stod_v2r8(s->x[14]);	\
  d ## 33 =  _fjsp_stod_v2r8(s->x[15]);	\
}

#define load_block4mat_ABCD( A,B,C,D, s ){	\
  A ## 00 =  _fjsp_stod_v2r8(s->x[0]);	\
  A ## 01 =  _fjsp_stod_v2r8(s->x[1]);	\
  A ## 10 =  _fjsp_stod_v2r8(s->x[2]);	\
  A ## 11 =  _fjsp_stod_v2r8(s->x[3]);	\
  B ## 00 =  _fjsp_stod_v2r8(s->x[4]);	\
  B ## 01 =  _fjsp_stod_v2r8(s->x[5]);	\
  B ## 10 =  _fjsp_stod_v2r8(s->x[6]);	\
  B ## 11 =  _fjsp_stod_v2r8(s->x[7]);	\
  C ## 00 =  _fjsp_stod_v2r8(s->x[8]);	\
  C ## 01 =  _fjsp_stod_v2r8(s->x[9]);	\
  C ## 10 =  _fjsp_stod_v2r8(s->x[10]);	\
  C ## 11 =  _fjsp_stod_v2r8(s->x[11]);	\
  D ## 00 =  _fjsp_stod_v2r8(s->x[12]);	\
  D ## 01 =  _fjsp_stod_v2r8(s->x[13]);	\
  D ## 10 =  _fjsp_stod_v2r8(s->x[14]);	\
  D ## 11 =  _fjsp_stod_v2r8(s->x[15]);	\
}




#define load_gauge( d, s ){		\
  d ## 00 =  _fjsp_stod_v2r8(s->u[0]);	\
  d ## 01 =  _fjsp_stod_v2r8(s->u[1]);	\
  d ## 02 =  _fjsp_stod_v2r8(s->u[2]);	\
  d ## 10 =  _fjsp_stod_v2r8(s->u[3]);	\
  d ## 11 =  _fjsp_stod_v2r8(s->u[4]);	\
  d ## 12 =  _fjsp_stod_v2r8(s->u[5]);	\
  d ## 20 =  _fjsp_stod_v2r8(s->u[6]);	\
  d ## 21 =  _fjsp_stod_v2r8(s->u[7]);	\
  d ## 22 =  _fjsp_stod_v2r8(s->u[8]);	\
}


#define load_ferm( d, s ){			\
  d ## c0sp1 =  _fjsp_stod_v2r8(s->y[0]);	\
  d ## c1sp1 =  _fjsp_stod_v2r8(s->y[1]);	\
  d ## c2sp1 =  _fjsp_stod_v2r8(s->y[2]);	\
  d ## c0sp2 =  _fjsp_stod_v2r8(s->y[3]);	\
  d ## c1sp2 =  _fjsp_stod_v2r8(s->y[4]);	\
  d ## c2sp2 =  _fjsp_stod_v2r8(s->y[5]);	\
  d ## c0sp3 =  _fjsp_stod_v2r8(s->y[6]);	\
  d ## c1sp3 =  _fjsp_stod_v2r8(s->y[7]);	\
  d ## c2sp3 =  _fjsp_stod_v2r8(s->y[8]);	\
  d ## c0sp4 =  _fjsp_stod_v2r8(s->y[9]);	\
  d ## c1sp4 =  _fjsp_stod_v2r8(s->y[10]);	\
  d ## c2sp4 =  _fjsp_stod_v2r8(s->y[11]);	\
}

#define load_ferm2( d, s ){		\
  d ## 0 =  _fjsp_stod_v2r8(s->y[0]);	\
  d ## 1 =  _fjsp_stod_v2r8(s->y[1]);	\
  d ## 2 =  _fjsp_stod_v2r8(s->y[2]);	\
  d ## 3 =  _fjsp_stod_v2r8(s->y[3]);	\
  d ## 4 =  _fjsp_stod_v2r8(s->y[4]);	\
  d ## 5 =  _fjsp_stod_v2r8(s->y[5]);	\
}


#define store_ferm( s, d ){			\
    s->y[0] = _fjsp_dtos_v2r4(d ## c0sp1);	\
    s->y[1] = _fjsp_dtos_v2r4(d ## c1sp1);	\
    s->y[2] = _fjsp_dtos_v2r4(d ## c2sp1);	\
    s->y[3] = _fjsp_dtos_v2r4(d ## c0sp2);	\
    s->y[4] = _fjsp_dtos_v2r4(d ## c1sp2);	\
    s->y[5] = _fjsp_dtos_v2r4(d ## c2sp2);	\
    s->y[6] = _fjsp_dtos_v2r4(d ## c0sp3);	\
    s->y[7] = _fjsp_dtos_v2r4(d ## c1sp3);	\
    s->y[8] = _fjsp_dtos_v2r4(d ## c2sp3);	\
    s->y[9] = _fjsp_dtos_v2r4(d ## c0sp4);	\
    s->y[10] = _fjsp_dtos_v2r4(d ## c1sp4);	\
    s->y[11] = _fjsp_dtos_v2r4(d ## c2sp4);	\
}


#define store_ferm2( s, d ){		\
    s->y[0] = _fjsp_dtos_v2r4(d ## 0 );	\
    s->y[1] = _fjsp_dtos_v2r4(d ## 1 );	\
    s->y[2] = _fjsp_dtos_v2r4(d ## 2 );	\
    s->y[3] = _fjsp_dtos_v2r4(d ## 3 );	\
    s->y[4] = _fjsp_dtos_v2r4(d ## 4 );	\
    s->y[5] = _fjsp_dtos_v2r4(d ## 5 );	\
  }


#define load_clov2( d, s ){		\
  d ## 0 =  _fjsp_stod_v2r8(s->c[0]);	\
  d ## 1 =  _fjsp_stod_v2r8(s->c[1]);	\
  d ## 2 =  _fjsp_stod_v2r8(s->c[2]);	\
  d ## 01 =  _fjsp_stod_v2r8(s->c[3]);	\
  d ## 02 =  _fjsp_stod_v2r8(s->c[4]);	\
  d ## 03 =  _fjsp_stod_v2r8(s->c[5]);	\
  d ## 04 =  _fjsp_stod_v2r8(s->c[6]);	\
  d ## 05 =  _fjsp_stod_v2r8(s->c[7]);	\
  d ## 12 =  _fjsp_stod_v2r8(s->c[8]);	\
  d ## 13 =  _fjsp_stod_v2r8(s->c[9]);	\
  d ## 14 =  _fjsp_stod_v2r8(s->c[10]);	\
  d ## 15 =  _fjsp_stod_v2r8(s->c[11]);	\
  d ## 23 =  _fjsp_stod_v2r8(s->c[12]);	\
  d ## 24 =  _fjsp_stod_v2r8(s->c[13]);	\
  d ## 25 =  _fjsp_stod_v2r8(s->c[14]);	\
  d ## 34 =  _fjsp_stod_v2r8(s->c[15]);	\
  d ## 35 =  _fjsp_stod_v2r8(s->c[16]);	\
  d ## 45 =  _fjsp_stod_v2r8(s->c[17]); \
}


// c=a*b as real
#define rmult(c,a,b){	    \
  c=_fjsp_mul_v2r8(a,b);    \
}


// c=a*b
#define cmult(c,a,b){			    \
    c=_fjsp_madd_cp_v2r8(a,b,_fjsp_setzero_v2r8()); \
    c=_fjsp_nmsub_cp_neg_sr12_v2r8(a,b,c);	    \
}
// c+=a*b
#define accum_cmult(c,a,b){			\
    c=_fjsp_madd_cp_v2r8(a,b,c);		\
    c=_fjsp_nmsub_cp_neg_sr12_v2r8(a,b,c);	\
}

// c-=a*b
#define accum_n_cmult(c,a,b){		\
    c=_fjsp_nmsub_cp_v2r8(a,b,c);	\
    c=_fjsp_madd_cp_neg_sr12_v2r8(a,b,c);	\
}


// c=conj(a) *b
#define conj_cmult(c,a,b){		          \
  c=_fjsp_madd_cp_v2r8(a,b,_fjsp_setzero_v2r8()); \
  c=_fjsp_madd_cp_neg_sr12_v2r8(a,b,c);	\
}

// c+=conj(a) *b
#define accum_conj_cmult(c,a,b){	\
  c=_fjsp_madd_cp_v2r8(a,b,c);		\
  c=_fjsp_madd_cp_neg_sr12_v2r8(a,b,c);	\
}

// c-=conj(a) *b
#define accum_n_conj_cmult(c,a,b){	\
  c=_fjsp_nmsub_cp_v2r8(a,b,c);		\
  c=_fjsp_nmsub_cp_neg_sr12_v2r8(a,b,c);\
}


#define prepare_I\
 register _fjsp_v2r8 one_one= _fjsp_set_v2r8(1.0,1.0); { } \


// c=I
#define set_I(c) {		\
    c=_fjsp_set_v2r8(0.0,1.0);	\
  }
// c=-I
#define set_minusI(c) {		\
    c=_fjsp_set_v2r8(0.0,-1.0);	\
  }


// c+= I*a
#define accum_times_I(c,a) {			\
    c = _fjsp_nmsub_neg_sr2_v2r8(one_one, a, c);\
}

// c = I*a
#define times_I(a)  _fjsp_nmsub_neg_sr2_v2r8(one_one, a, _fjsp_setzero_v2r8())

// c+= -I*a
#define accum_times_minusI(c,a) {		\
    c = _fjsp_madd_neg_sr2_v2r8(one_one, a,  c);\
}

// c= -I*a
#define times_minusI(a)  _fjsp_madd_neg_sr2_v2r8(one_one, a, _fjsp_setzero_v2r8())


#endif
#endif
