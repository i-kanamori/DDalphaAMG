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

#ifndef DIRAC_K_H
#define DIRAC_K_H
#ifdef K_OPT

#include "K_vector.h"

void K_set_clover_double( double *out, complex_double *in );
void K_set_clover_float( float *out, complex_double *in );
void K_site_clover_double( double *eta, const double *phi, const double *clover );
void K_site_clover_float( float *eta, const float *phi, float *clover );

void K_site_clover_invert_double( double *clover_in, double *clover_out );
void K_site_clover_invert_float( float *clover_in, float *clover_out );
void block_oddeven_coupling_double( double *eta,  double *D, double *phi, 
				   int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor );
void block_oddeven_coupling_float( float *eta,  float *D, float *phi, 
				   int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor );


void block_n_oddeven_coupling_double( double *eta,  double *D, double *phi, 
				   int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor );
void block_n_oddeven_coupling_float( float *eta,  float *D, float *phi, 
				   int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor );

void K_block_d_float( complex_float *ceta, complex_float *cphi, int start, schwarz_float_struct *s);

void K_d_plus_clover_aggregate_float( vector_float eta1, vector_float eta2, vector_float phi0, schwarz_float_struct *s, level_struct *l, struct Thread *threading );

void K_d_neighbor_aggregate_float( vector_float eta1, vector_float eta2, vector_float phi, const int mu, schwarz_float_struct *s, level_struct *l, struct Thread *threading );

  
static inline void K_site_clover_2spin_float( float *eta, const float *phi, float *clover ) {

    register _fjsp_v2r8 in0,in1,in2,in3,in4,in5;
    register _fjsp_v2r8 out0,out1,out2,out3,out4,out5;

    register _fjsp_v2r8 clov0,clov1,clov2;
    register _fjsp_v2r8 clov01,clov02,clov03,clov04,clov05,
      clov12,clov13,clov14,clov15,
      clov23,clov24,clov25,
      clov34,clov35,
      clov45;
    register _fjsp_v2r8 diag0,diag1,diag2,diag3,diag4,diag5;

    // load input spinor
    register ssu3ferm2* in_ptr=(ssu3ferm2*)phi;
    load_ferm2(in,in_ptr);
  
    // load clover term
    register ssu3clov2* clov_ptr=(ssu3clov2*)clover;
    load_clov2(clov, clov_ptr);

    // extract diagnaols (=they are real)
    diag0 =_fjsp_unpacklo_v2r8(clov0,clov0);
    diag1 =_fjsp_unpackhi_v2r8(clov0,clov0);
    diag2 =_fjsp_unpacklo_v2r8(clov1,clov1);
    diag3 =_fjsp_unpackhi_v2r8(clov1,clov1);
    diag4 =_fjsp_unpacklo_v2r8(clov2,clov2);
    diag5 =_fjsp_unpackhi_v2r8(clov2,clov2);

    // 1st line
    rmult(out0,diag0,in0);
    accum_cmult(out0, clov01, in1);
    accum_cmult(out0, clov02, in2);
    accum_cmult(out0, clov03, in3);
    accum_cmult(out0, clov04, in4);
    accum_cmult(out0, clov05, in5);

    // 2nd line
    rmult(out1,diag1,in1);
    accum_conj_cmult(out1,clov01,in0);
    accum_cmult(out1, clov12, in2);
    accum_cmult(out1, clov13, in3);
    accum_cmult(out1, clov14, in4);
    accum_cmult(out1, clov15, in5);

    // 3rd line
    rmult(out2,diag2,in2);
    accum_conj_cmult(out2,clov02,in0);
    accum_conj_cmult(out2,clov12,in1);
    accum_cmult(out2, clov23, in3);
    accum_cmult(out2, clov24, in4);
    accum_cmult(out2, clov25, in5);

    // 4th line
    rmult(out3,diag3,in3);
    accum_conj_cmult(out3,clov03,in0);
    accum_conj_cmult(out3,clov13,in1);
    accum_conj_cmult(out3,clov23,in2);
    accum_cmult(out3, clov34, in4);
    accum_cmult(out3, clov35, in5);

    // 5th line
    rmult(out4,diag4,in4);
    accum_conj_cmult(out4,clov04,in0);
    accum_conj_cmult(out4,clov14,in1);
    accum_conj_cmult(out4,clov24,in2);
    accum_conj_cmult(out4,clov34,in3);
    accum_cmult(out4, clov45, in5);

    // 6th line
    rmult(out5,diag5,in5);
    accum_conj_cmult(out5,clov05,in0);
    accum_conj_cmult(out5,clov15,in1);
    accum_conj_cmult(out5,clov25,in2);
    accum_conj_cmult(out5,clov35,in3);
    accum_conj_cmult(out5,clov45,in4);

    register ssu3ferm2* out_ptr=(ssu3ferm2*)eta;
    store_ferm2(out_ptr, out);

  return;
}



#if 0
void prp_double( complex_double *prn[4], complex_double *phi, int start, int end );
void prp_float( complex_float *prn[4], complex_float *phi, int start, int end );
void prn_su3_double( complex_double *prp[4], complex_double *phi, operator_double_struct *op, int *neighbor, int start, int end );
void prn_su3_float( complex_float *prp[4], complex_float *phi, operator_float_struct *op, int *neighbor, int start, int end );
void pbn_double( complex_double *eta, complex_double *prp[4], int start, int end );
void pbn_float( complex_float *eta, complex_float *prp[4], int start, int end );
void su3_pbp_double( complex_double* eta, complex_double *prn[4], operator_double_struct *op, int *neighbor, int start, int end );
void su3_pbp_float( complex_float* eta, complex_float *prn[4], operator_float_struct *op, int *neighbor, int start, int end );

void block_oddeven_plus_coupling_double( double *eta, double *D, double *phi, int mu,
                                         int start, int end, int *ind, int *neighbor );
void block_oddeven_plus_coupling_float( float *eta, float *D, float *phi, int mu,
                                        int start, int end, int *ind, int *neighbor );
void block_oddeven_minus_coupling_double( double *eta, double *D, double *phi, int mu,
                                          int start, int end, int *ind, int *neighbor );
void block_oddeven_minus_coupling_float( float *eta, float *D, float *phi, int mu,
                                         int start, int end, int *ind, int *neighbor );
void block_oddeven_nplus_coupling_double( double *eta, double *D, double *phi, int mu,
                                         int start, int end, int *ind, int *neighbor );
void block_oddeven_nplus_coupling_float( float *eta, float *D, float *phi, int mu,
                                        int start, int end, int *ind, int *neighbor );
void block_oddeven_nminus_coupling_double( double *eta, double *D, double *phi, int mu,
                                          int start, int end, int *ind, int *neighbor );
void block_oddeven_nminus_coupling_float( float *eta, float *D, float *phi, int mu,
                                         int start, int end, int *ind, int *neighbor );
void boundary_minus_coupling_double( double *eta, double *D, double *phi, int mu,
                                              int start, int end, int *ind, int *neighbor );
void boundary_minus_coupling_float( float *eta, float *D, float *phi, int mu,
                                             int start, int end, int *ind, int *neighbor );
void boundary_plus_coupling_double( double *eta, double *D, double *phi, int mu,
                                             int start, int end, int *ind, int *neighbor );
void boundary_plus_coupling_float( float *eta, float *D, float *phi, int mu,
                                            int start, int end, int *ind, int *neighbor );
void boundary_nminus_coupling_double( double *eta, double *D, double *phi, int mu,
                                              int start, int end, int *ind, int *neighbor );
void boundary_nminus_coupling_float( float *eta, float *D, float *phi, int mu,
                                             int start, int end, int *ind, int *neighbor );
void boundary_nplus_coupling_double( double *eta, double *D, double *phi, int mu,
                                             int start, int end, int *ind, int *neighbor );
void boundary_nplus_coupling_float( float *eta, float *D, float *phi, int mu,
                                            int start, int end, int *ind, int *neighbor );

void sse_set_clover_double( double *out, complex_double *in );
void sse_set_clover_float( float *out, complex_double *in );
void sse_site_clover_double( double *eta, const double *phi, const double *clover );
void sse_site_clover_float( float *eta, const float *phi, float *clover );

void sse_site_clover_invert_double( double *clover_in, double *clover_out );
void sse_site_clover_invert_float( float *clover_in, float *clover_out );


static inline void sse_mvm_double_simd_length( const complex_double *eta, const complex_double *D, const complex_double *phi ) {}

static inline void sse_mvm_float_simd_length(
    const complex_float *eta, const complex_float *D, const complex_float *phi ) {
#ifdef SSE
  __m128 gauge_re;
  __m128 gauge_im;
  __m128 in_re[3];
  __m128 in_im[3];
  __m128 out_re[3];
  __m128 out_im[3];

  int elements = SIMD_LENGTH_float;

  // j runs over all right-hand sides, using vectorization
  for(int i=0; i<3; i++) {
    in_re[i] = _mm_load_ps((float *)(phi+i*elements));
    in_im[i] = _mm_load_ps((float *)(phi+i*elements)+elements);
  }

  for(int i=0; i<3; i++) {
    gauge_re = _mm_set1_ps(creal(D[0+3*i]));
    gauge_im = _mm_set1_ps(cimag(D[0+3*i]));
    cmul(gauge_re, gauge_im, in_re[0], in_im[0], out_re+i, out_im+i);
    gauge_re = _mm_set1_ps(creal(D[1+3*i]));
    gauge_im = _mm_set1_ps(cimag(D[1+3*i]));
    cfmadd(gauge_re, gauge_im, in_re[1], in_im[1], out_re+i, out_im+i);
    gauge_re = _mm_set1_ps(creal(D[2+3*i]));
    gauge_im = _mm_set1_ps(cimag(D[2+3*i]));
    cfmadd(gauge_re, gauge_im, in_re[2], in_im[2], out_re+i, out_im+i);
  }

  for(int i=0; i<3; i++) {
    _mm_store_ps((float *)(eta+i*elements),          out_re[i]);
    _mm_store_ps((float *)(eta+i*elements)+elements, out_im[i]);
  }
#endif
}


static inline void sse_mvm_double( const complex_double *eta, const complex_double *D,
                                   const complex_double *phi, int elements ) {}
// spinors are vectorized, gauge is same for all (use for multiple rhs)
static inline void sse_mvm_float( const complex_float *eta, const complex_float *D,
                                  const complex_float *phi, int elements ) {
#ifdef SSE
  __m128 gauge_re;
  __m128 gauge_im;
  __m128 in_re[3];
  __m128 in_im[3];
  __m128 out_re[3];
  __m128 out_im[3];

  // j runs over all right-hand sides, using vectorization
  for(int j=0; j<elements; j+=SIMD_LENGTH_float) {
    for(int i=0; i<3; i++) {
      in_re[i] = _mm_load_ps((float *)(phi+i*elements)+j);
      in_im[i] = _mm_load_ps((float *)(phi+i*elements)+j+elements);
    }

    for(int i=0; i<3; i++) {
      gauge_re = _mm_set1_ps(creal(D[0+3*i]));
      gauge_im = _mm_set1_ps(cimag(D[0+3*i]));
      cmul(gauge_re, gauge_im, in_re[0], in_im[0], out_re+i, out_im+i);
      gauge_re = _mm_set1_ps(creal(D[1+3*i]));
      gauge_im = _mm_set1_ps(cimag(D[1+3*i]));
      cfmadd(gauge_re, gauge_im, in_re[1], in_im[1], out_re+i, out_im+i);
      gauge_re = _mm_set1_ps(creal(D[2+3*i]));
      gauge_im = _mm_set1_ps(cimag(D[2+3*i]));
      cfmadd(gauge_re, gauge_im, in_re[2], in_im[2], out_re+i, out_im+i);
    }

    for(int i=0; i<3; i++) {
      _mm_store_ps((float *)(eta+i*elements)+j,          out_re[i]);
      _mm_store_ps((float *)(eta+i*elements)+j+elements, out_im[i]);
    }
  }
#endif
}


static inline void sse_mvmh_double( const complex_double *eta, const complex_double *D,
                                    const complex_double *phi, int elements ) {}
// spinors are vectorized, gauge is same for all (use for multiple rhs)
static inline void sse_mvmh_float( const complex_float *eta, const complex_float *D,
                                   const complex_float *phi, int elements ) {
#ifdef SSE
  __m128 gauge_re;
  __m128 gauge_im;
  __m128 in_re[3];
  __m128 in_im[3];
  __m128 out_re[3];
  __m128 out_im[3];


  // j runs over all right-hand sides, using vectorization
  for(int j=0; j<elements; j+=SIMD_LENGTH_float) {
    for(int i=0; i<3; i++) {
      in_re[i] = _mm_load_ps((float *)(phi+i*elements)+j);
      in_im[i] = _mm_load_ps((float *)(phi+i*elements)+j+elements);
    }

    for(int i=0; i<3; i++) {
      gauge_re = _mm_set1_ps(creal(D[0+i]));
      gauge_im = _mm_set1_ps(cimag(D[0+i]));
      cmul_conj(gauge_re, gauge_im, in_re[0], in_im[0], out_re+i, out_im+i);
      gauge_re = _mm_set1_ps(creal(D[3+i]));
      gauge_im = _mm_set1_ps(cimag(D[3+i]));
      cfmadd_conj(gauge_re, gauge_im, in_re[1], in_im[1], out_re+i, out_im+i);
      gauge_re = _mm_set1_ps(creal(D[6+i]));
      gauge_im = _mm_set1_ps(cimag(D[6+i]));
      cfmadd_conj(gauge_re, gauge_im, in_re[2], in_im[2], out_re+i, out_im+i);
    }

    for(int i=0; i<3; i++) {
      _mm_store_ps((float *)(eta+i*elements)+j,          out_re[i]);
      _mm_store_ps((float *)(eta+i*elements)+j+elements, out_im[i]);
    }
  }
#endif
}


static inline void sse_twospin_double( complex_double *out_spin0and1, complex_double *out_spin2and3, const complex_double *in, int elements, int mu, double sign ) {}
// mu is according to the enum for T,Z,Y,X defined in clifford.h
static inline void sse_twospin_float( complex_float *out_spin0and1, complex_float *out_spin2and3, const complex_float *in, int elements, int mu, double sign ) {

#ifdef SSE
  __m128 scale_re;
  __m128 scale_im;
  complex_float *out;
  // components 0-5  are subtracted from out_spin0_and1
  // components 6-11 are subtracted from out_spin2_and3

  // 6 complex components = 12
  for(int i=0; i<12*elements; i+=SIMD_LENGTH_float) {
    __m128 tmp = _mm_load_ps((float *)in            + i);
    __m128 out = _mm_load_ps((float *)out_spin0and1 + i);

    out = _mm_sub_ps(out, tmp);

    _mm_store_ps((float *)out_spin0and1 + i, out);
  }
  for(int i=12*elements; i<24*elements; i+=SIMD_LENGTH_float) {
    __m128 tmp = _mm_load_ps((float *)in            + i);
    __m128 out = _mm_load_ps((float *)out_spin2and3 + i);

    out = _mm_sub_ps(out, tmp);

    _mm_store_ps((float *)out_spin2and3 + i, out);
  }

  out = out_spin2and3;
  for(int spin=0; spin<4; spin++) {
    if(spin == 2)
      out = out_spin0and1;
    scale_re = _mm_set1_ps(sign*creal(gamma_val[mu][spin]));
    scale_im = _mm_set1_ps(sign*cimag(gamma_val[mu][spin]));
    // factors of 2 are for complex
    for(int j=0; j<3; j++) {
      for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
        __m128 in_re  = _mm_load_ps((float *)in  + i + (2*(3*gamma_co[mu][spin]+j)+0)*elements);
        __m128 in_im  = _mm_load_ps((float *)in  + i + (2*(3*gamma_co[mu][spin]+j)+1)*elements);
        __m128 out_re = _mm_load_ps((float *)out + i + (2*(3*spin+j)+0)*elements);
        __m128 out_im = _mm_load_ps((float *)out + i + (2*(3*spin+j)+1)*elements);

        cfmadd(scale_re, scale_im, in_re, in_im, &out_re, &out_im);

        _mm_store_ps((float *)out + i + (2*(3*spin+j)+0)*elements, out_re);
        _mm_store_ps((float *)out + i + (2*(3*spin+j)+1)*elements, out_im);
      }
    }
  }
#endif
}


static inline void sse_twospin2_p_double_simd_length( complex_double *out_spin0and1, complex_double *out_spin2and3, const complex_double *in, int mu ) {}
// mu is according to the enum for T,Z,Y,X defined in clifford.h
static inline void sse_twospin2_p_float_simd_length( complex_float *out_spin0and1, complex_float *out_spin2and3, const complex_float *in, int mu ) {
#ifdef SSE
  __m128 scale_re;
  __m128 scale_im;
//   __m128 out_re;
//   __m128 out_im;
  complex_float *out;
  int elements = SIMD_LENGTH_float;
  // components 0-5  are subtracted from out_spin0_and1
  // components 6-11 are subtracted from out_spin2_and3

  // 6 complex components = 12
  for(int i=0; i<12*elements; i+=SIMD_LENGTH_float) {
    __m128 tmp = _mm_load_ps((float *)in + i);
    _mm_store_ps((float *)out_spin0and1 + i, tmp);
  }
  for(int i=12*elements; i<24*elements; i+=SIMD_LENGTH_float) {
    __m128 tmp = _mm_load_ps((float *)in + i);
    _mm_store_ps((float *)out_spin2and3 + i, tmp);
  }

  out = out_spin2and3;
  for(int spin=0; spin<4; spin++) {
    if(spin == 2)
      out = out_spin0and1;
    scale_re = _mm_set1_ps(-creal(gamma_val[mu][spin]));
    scale_im = _mm_set1_ps(-cimag(gamma_val[mu][spin]));
    // factors of 2 are for complex
    for(int j=0; j<3; j++) {
      __m128 in_re  = _mm_load_ps((float *)in   + (2*(3*gamma_co[mu][spin]+j)+0)*elements);
      __m128 in_im  = _mm_load_ps((float *)in   + (2*(3*gamma_co[mu][spin]+j)+1)*elements);
      __m128 out_re = _mm_load_ps((float *)out  + (2*(3*spin+j)+0)*elements);
      __m128 out_im = _mm_load_ps((float *)out  + (2*(3*spin+j)+1)*elements);

      cmul(scale_re, scale_im, in_re, in_im, &out_re, &out_im);

      _mm_store_ps((float *)out  + (2*(3*spin+j)+0)*elements, out_re);
      _mm_store_ps((float *)out  + (2*(3*spin+j)+1)*elements, out_im);
    }
  }
#endif
}


static inline void sse_twospin2_p_double( complex_double *out_spin0and1, complex_double *out_spin2and3, const complex_double *in, int elements, int mu ) {}
// mu is according to the enum for T,Z,Y,X defined in clifford.h
static inline void sse_twospin2_p_float( complex_float *out_spin0and1, complex_float *out_spin2and3, const complex_float *in, int elements, int mu ) {

#ifdef SSE
  __m128 scale_re;
  __m128 scale_im;
//   __m128 out_re;
//   __m128 out_im;
  complex_float *out;
  // components 0-5  are subtracted from out_spin0_and1
  // components 6-11 are subtracted from out_spin2_and3

  // 6 complex components = 12
  for(int i=0; i<12*elements; i+=SIMD_LENGTH_float) {
    __m128 tmp = _mm_load_ps((float *)in + i);
    _mm_store_ps((float *)out_spin0and1 + i, tmp);
  }
  for(int i=12*elements; i<24*elements; i+=SIMD_LENGTH_float) {
    __m128 tmp = _mm_load_ps((float *)in + i);
    _mm_store_ps((float *)out_spin2and3 + i, tmp);
  }

  out = out_spin2and3;
  for(int spin=0; spin<4; spin++) {
    if(spin == 2)
      out = out_spin0and1;
    scale_re = _mm_set1_ps(-creal(gamma_val[mu][spin]));
    scale_im = _mm_set1_ps(-cimag(gamma_val[mu][spin]));
    // factors of 2 are for complex
    for(int j=0; j<3; j++) {
      for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
        __m128 in_re  = _mm_load_ps((float *)in  + i + (2*(3*gamma_co[mu][spin]+j)+0)*elements);
        __m128 in_im  = _mm_load_ps((float *)in  + i + (2*(3*gamma_co[mu][spin]+j)+1)*elements);
        __m128 out_re = _mm_load_ps((float *)out + i + (2*(3*spin+j)+0)*elements);
        __m128 out_im = _mm_load_ps((float *)out + i + (2*(3*spin+j)+1)*elements);

        cmul(scale_re, scale_im, in_re, in_im, &out_re, &out_im);

        _mm_store_ps((float *)out + i + (2*(3*spin+j)+0)*elements, out_re);
        _mm_store_ps((float *)out + i + (2*(3*spin+j)+1)*elements, out_im);
      }
    }
  }
#endif
}


static inline void sse_spin0and1_site_clover_double( const complex_double *eta, const complex_double *phi, const config_double clover, double shift, int elements ) {}

static inline void sse_spin0and1_site_clover_float( const complex_float *eta, const complex_float *phi, const config_float clover, double shift, int elements ) {
#ifdef SSE
  // offset computations 2*index+0/1 are for real and imaginary parts

  // diagonal
  if ( g.csw == 0.0 ) {
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
      for(int j=0; j<6; j++) {
        __m128 factor = _mm_set1_ps((float)shift);
        __m128 in_re  = _mm_load_ps((float *)phi + i + (2*j+0)*elements);
        __m128 in_im  = _mm_load_ps((float *)phi + i + (2*j+1)*elements);

        in_re = _mm_mul_ps( factor, in_re );
        in_im = _mm_mul_ps( factor, in_im );

        _mm_store_ps((float *)eta + i + (2*j+0)*elements, in_re);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, in_im);
      }
      __m128 zero = _mm_setzero_ps();
      for(int j=6; j<12; j++) {
        _mm_store_ps((float *)eta + i + (2*j+0)*elements, zero);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, zero);
      }
    }
  } else {
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
      for(int j=0; j<6; j++) {
        __m128 clover_re = _mm_set1_ps(creal(clover[j]));
        __m128 clover_im = _mm_set1_ps(cimag(clover[j]));
        __m128 in_re  = _mm_load_ps((float *)phi + i + (2*j+0)*elements);
        __m128 in_im  = _mm_load_ps((float *)phi + i + (2*j+1)*elements);
        __m128 out_re; __m128 out_im;
        
        cmul(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

        _mm_store_ps((float *)eta + i + (2*j+0)*elements, out_re);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, out_im);
      }
      __m128 zero = _mm_setzero_ps();
      for(int j=6; j<12; j++) {
        _mm_store_ps((float *)eta + i + (2*j+0)*elements, zero);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, zero);
      }
    }

    // spin 0 and 1
    __m128 clover_re;
    __m128 clover_im;
    __m128 in_re;
    __m128 in_im;
    __m128 out_re;
    __m128 out_im;
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {

      // io = index out, ii = index in, ic = index clover
      int ic = 12;
      for(int io=0; io<5; io++) {
        out_re = _mm_load_ps((float *)eta + i + (2*io+0)*elements);
        out_im = _mm_load_ps((float *)eta + i + (2*io+1)*elements);
        for(int ii=io+1; ii<=5; ii++) {
          clover_re = _mm_set1_ps(creal(clover[ic]));
          clover_im = _mm_set1_ps(cimag(clover[ic]));
          ic++;
          in_re  = _mm_load_ps((float *)phi + i + (2*ii+0)*elements);
          in_im  = _mm_load_ps((float *)phi + i + (2*ii+1)*elements);

          cfmadd(clover_re, clover_im, in_re, in_im, &out_re, &out_im);
        }
        _mm_store_ps((float *)eta + i + (2*io+0)*elements, out_re);
        _mm_store_ps((float *)eta + i + (2*io+1)*elements, out_im);
      }

      ic = 12;
      for(int ii=0; ii<5; ii++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*ii+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*ii+1)*elements);
        for(int io=ii+1; io<=5; io++) {
          clover_re = _mm_set1_ps(creal(clover[ic]));
          clover_im = _mm_set1_ps(cimag(clover[ic]));
          ic++;
          out_re = _mm_load_ps((float *)eta + i + (2*io+0)*elements);
          out_im = _mm_load_ps((float *)eta + i + (2*io+1)*elements);

          cfmadd_conj(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta + i + (2*io+0)*elements, out_re);
          _mm_store_ps((float *)eta + i + (2*io+1)*elements, out_im);
        }
      }
    }
  }
#endif
}


static inline void sse_spin2and3_site_clover_double( const complex_double *eta, const complex_double *phi, const config_double clover, double shift, int elements ) {}

static inline void sse_spin2and3_site_clover_float( const complex_float *eta, const complex_float *phi, const config_float clover, double shift, int elements ) {
#ifdef SSE
  // offset computations 2*index+0/1 are for real and imaginary parts

  // diagonal
  if ( g.csw == 0.0 ) {
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
      __m128 zero = _mm_setzero_ps();
      for(int j=0; j<6; j++) {
        _mm_store_ps((float *)eta + i + (2*j+0)*elements, zero);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, zero);
      }
      for(int j=6; j<12; j++) {
        __m128 factor = _mm_set1_ps((float)shift);
        __m128 in_re  = _mm_load_ps((float *)phi + i + (2*j+0)*elements);
        __m128 in_im  = _mm_load_ps((float *)phi + i + (2*j+1)*elements);

        in_re = _mm_mul_ps( factor, in_re );
        in_im = _mm_mul_ps( factor, in_im );

        _mm_store_ps((float *)eta + i + (2*j+0)*elements, in_re);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, in_im);
      }
    }
    
    
    
  } else {
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {
      __m128 zero = _mm_setzero_ps();
      for(int j=0; j<6; j++) {
        _mm_store_ps((float *)eta + i + (2*j+0)*elements, zero);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, zero);
      }
      for(int j=6; j<12; j++) {
        __m128 clover_re = _mm_set1_ps(creal(clover[j]));
        __m128 clover_im = _mm_set1_ps(cimag(clover[j]));
        __m128 in_re  = _mm_load_ps((float *)phi + i + (2*j+0)*elements);
        __m128 in_im  = _mm_load_ps((float *)phi + i + (2*j+1)*elements);
        __m128 out_re; __m128 out_im;

        cmul(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

        _mm_store_ps((float *)eta + i + (2*j+0)*elements, out_re);
        _mm_store_ps((float *)eta + i + (2*j+1)*elements, out_im);
      }
    }

    // spin 0 and 1
    __m128 clover_re;
    __m128 clover_im;
    __m128 in_re;
    __m128 in_im;
    __m128 out_re;
    __m128 out_im;
    for(int i=0; i<elements; i+=SIMD_LENGTH_float) {

      // io = index out, ii = index in, ic = index clover
      int ic = 27;
      for(int io=6; io<11; io++) {
        out_re = _mm_load_ps((float *)eta + i + (2*io+0)*elements);
        out_im = _mm_load_ps((float *)eta + i + (2*io+1)*elements);
        for(int ii=io+1; ii<=11; ii++) {
          clover_re = _mm_set1_ps(creal(clover[ic]));
          clover_im = _mm_set1_ps(cimag(clover[ic]));
          ic++;
          in_re  = _mm_load_ps((float *)phi + i + (2*ii+0)*elements);
          in_im  = _mm_load_ps((float *)phi + i + (2*ii+1)*elements);

          cfmadd(clover_re, clover_im, in_re, in_im, &out_re, &out_im);
        }
        _mm_store_ps((float *)eta + i + (2*io+0)*elements, out_re);
        _mm_store_ps((float *)eta + i + (2*io+1)*elements, out_im);
      }

      ic = 27;
      for(int ii=6; ii<11; ii++) {
        in_re  = _mm_load_ps((float *)phi + i + (2*ii+0)*elements);
        in_im  = _mm_load_ps((float *)phi + i + (2*ii+1)*elements);
        for(int io=ii+1; io<=11; io++) {
          clover_re = _mm_set1_ps(creal(clover[ic]));
          clover_im = _mm_set1_ps(cimag(clover[ic]));
          ic++;
          out_re = _mm_load_ps((float *)eta + i + (2*io+0)*elements);
          out_im = _mm_load_ps((float *)eta + i + (2*io+1)*elements);

          cfmadd_conj(clover_re, clover_im, in_re, in_im, &out_re, &out_im);

          _mm_store_ps((float *)eta + i + (2*io+0)*elements, out_re);
          _mm_store_ps((float *)eta + i + (2*io+1)*elements, out_im);
        }
      }
    }
  }
#endif
}

#endif //0
#endif
#endif // DIRAC_SSE_H
