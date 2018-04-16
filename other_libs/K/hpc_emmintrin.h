/*

Copyright (C) 2015, Ken-Ichi Ishikawa

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
*/

#ifndef HPC_EMMINTRIN_H
#define HPC_EMMINTRIN_H

#ifdef _FJ_HPC

#ifdef _FJ_HPC_EMU

//
// Emulate Fujitsu HPC-ACE intrinsics (not full)
//

typedef struct {
  double lo;
  double hi;
} _fjsp_v2r8;

typedef struct {
  float lo;
  float hi;
} _fjsp_v2r4;

static inline _fjsp_v2r8 _fjsp_set_v2r8(double a, double b){
  _fjsp_v2r8 out;
  out.lo = b;
  out.hi = a;
  return out;
}

static inline _fjsp_v2r8 _fjsp_setzero_v2r8(void){
  _fjsp_v2r8 out;
  out.lo = 0.0;
  out.hi = 0.0;
  return out;
}

static inline _fjsp_v2r8 _fjsp_load_v2r8(double const *p){
  _fjsp_v2r8 out;
  out.lo = *(p+0);
  out.hi = *(p+1);
  return out;
}

static inline _fjsp_v2r8 _fjsp_loadh_v2r8(_fjsp_v2r8 a,double const *p){
  _fjsp_v2r8 out;
  out.lo = a.lo;
  out.hi = *(p+0);
  return out;
}

static inline _fjsp_v2r8 _fjsp_loadl_v2r8(_fjsp_v2r8 a,double const *p){
  _fjsp_v2r8 out;
  out.lo = *(p+0);
  out.hi = a.hi;
  return out;
}

static inline void _fjsp_store_v2r8(double *p, _fjsp_v2r8 a){
  *(p+0) = a.lo;
  *(p+1) = a.hi;
}

static inline void _fjsp_storeh_v2r8(double *p, _fjsp_v2r8 a){
  *(p+0) = a.hi;
}

static inline void _fjsp_storel_v2r8(double *p, _fjsp_v2r8 a){
  *(p+0) = a.lo;
}

static inline _fjsp_v2r8 _fjsp_add_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b){
  _fjsp_v2r8 out;
  out.lo = a.lo + b.lo;
  out.hi = a.hi + b.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_sub_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b){
  _fjsp_v2r8 out;
  out.lo = a.lo - b.lo;
  out.hi = a.hi - b.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_mul_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b){
  _fjsp_v2r8 out;
  out.lo = a.lo * b.lo;
  out.hi = a.hi * b.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_unpackhi_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b){
  _fjsp_v2r8 out;
  out.lo = a.hi;
  out.hi = b.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_unpacklo_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b){
  _fjsp_v2r8 out;
  out.lo = a.lo;
  out.hi = b.lo;
  return out;
}

static inline _fjsp_v2r8 _fjsp_shuffle_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, unsigned int n){
  _fjsp_v2r8 out;
  if (0 == (n & 1)) {
    out.lo = a.lo;
  } else {
    out.lo = a.hi;
  }
  if (0 == (n & 2)) {
    out.hi = b.lo;
  } else {
    out.hi = b.hi;
  }
  return out;
}

static inline _fjsp_v2r8 _fjsp_neg_v2r8(_fjsp_v2r8 a){
  _fjsp_v2r8 out;
  out.lo = -a.lo;
  out.hi = -a.hi;
  return out;
}

static inline _fjsp_v2r4 _fjsp_dtos_v2r4(_fjsp_v2r8 a){
  _fjsp_v2r4 out;
  out.lo = (float )a.lo;
  out.hi = (float )a.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_stod_v2r8(_fjsp_v2r4 a){
  _fjsp_v2r8 out;
  out.lo = (double )a.lo;
  out.hi = (double )a.hi;
  return out;
}

//
// fused multiply and add
//
static inline _fjsp_v2r8 _fjsp_madd_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // a * b + c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.lo + c.lo;
  out.hi = a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // a * b - c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.lo - c.lo;
  out.hi = a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // -(a * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.lo - c.lo;
  out.hi = - a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // -(a * b - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.lo + c.lo;
  out.hi = - a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // a.lo * b + c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.lo + c.lo;
  out.hi = a.lo*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // a.lo * b - c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.lo - c.lo;
  out.hi = a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  - (a.lo * b + c)
  _fjsp_v2r8 out;
  out.hi = - a.lo*b.hi - c.hi;
  out.lo = - a.lo*b.lo - c.lo;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  - (a.lo * b - c)
  _fjsp_v2r8 out;
  out.hi = - a.lo*b.hi + c.hi;
  out.lo = - a.lo*b.lo + c.lo;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a * b + c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.lo + c.lo;
  out.hi = - a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a * b - c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.lo - c.lo;
  out.hi = - a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.lo - c.lo;
  out.hi =   a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a * b - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.lo + c.lo;
  out.hi =   a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.lo * b + c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.lo + c.lo;
  out.hi = - a.lo*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.lo * b - c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.lo - c.lo;
  out.hi = - a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a.lo * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.lo - c.lo;
  out.hi =   a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_neg_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a.lo * b - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.lo + c.lo;
  out.hi =   a.lo*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a^t * b + c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.lo + c.lo;
  out.hi = a.lo*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a^t * b - c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.lo - c.lo;
  out.hi = a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a^t * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo - c.lo;
  out.hi = - a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a^t * b - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo + c.lo;
  out.hi = - a.lo*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a.hi * b + c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.lo + c.lo;
  out.hi = a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a.hi * b - c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.lo - c.lo;
  out.hi = a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a.hi * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo - c.lo;
  out.hi = - a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a.hi * b - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo + c.lo;
  out.hi = - a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a^t * b + c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.lo + c.lo;
  out.hi = - a.lo*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a^t * b - c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.lo - c.lo;
  out.hi = - a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a^t * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo - c.lo;
  out.hi =   a.lo*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a^t * b - c)
  _fjsp_v2r8 out;
  out.hi = - a.lo*b.hi + c.hi;
  out.lo =   a.hi*b.lo + c.lo;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //   (-hi) a.hi * b + c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.lo + c.lo;
  out.hi = - a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.hi * b - c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.lo - c.lo;
  out.hi = - a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi)a.hi * b + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo - c.lo;
  out.hi =   a.hi*b.hi - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_neg_sr1_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi)a.hi * b - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.lo + c.lo;
  out.hi =   a.hi*b.hi + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a * b^t + c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.hi + c.lo;
  out.hi = a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a * b^t - c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.hi - c.lo;
  out.hi = a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi - c.lo;
  out.hi = - a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  // -(a * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi + c.lo;
  out.hi = - a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a.lo * b^t + c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.hi + c.lo;
  out.hi = a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a.lo * b^t - c
  _fjsp_v2r8 out;
  out.lo = a.lo*b.hi - c.lo;
  out.hi = a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a.lo * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi - c.lo;
  out.hi = - a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a.lo * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi + c.lo;
  out.hi = - a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a * b^t + c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.hi + c.lo;
  out.hi = - a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a * b^t - c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.hi - c.lo;
  out.hi = - a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi - c.lo;
  out.hi =   a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi + c.lo;
  out.hi =   a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.lo * b^t + c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.hi + c.lo;
  out.hi = - a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.lo * b^t - c
  _fjsp_v2r8 out;
  out.lo =   a.lo*b.hi - c.lo;
  out.hi = - a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a.lo * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi - c.lo;
  out.hi =   a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_neg_sr2_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a.lo * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.lo*b.hi + c.lo;
  out.hi =   a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a^t * b^t + c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.hi + c.lo;
  out.hi = a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a^t * b^t - c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.hi - c.lo;
  out.hi = a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a^t * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi - c.lo;
  out.hi = - a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a^t * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi + c.lo;
  out.hi = - a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a.hi * b^t + c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.hi + c.lo;
  out.hi = a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  a.hi * b^t - c
  _fjsp_v2r8 out;
  out.lo = a.hi*b.hi - c.lo;
  out.hi = a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a.hi * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi - c.lo;
  out.hi = - a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -(a.hi * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi + c.lo;
  out.hi = - a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a^t * b^t + c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.hi + c.lo;
  out.hi = - a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a^t * b^t - c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.hi - c.lo;
  out.hi = - a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a^t * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi - c.lo;
  out.hi =   a.lo*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a^t * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi + c.lo;
  out.hi =   a.lo*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_madd_cp_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.hi * b^t + c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.hi + c.lo;
  out.hi = - a.hi*b.lo + c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_msub_cp_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  (-hi) a.hi * b^t - c
  _fjsp_v2r8 out;
  out.lo =   a.hi*b.hi - c.lo;
  out.hi = - a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmadd_cp_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a.hi * b^t + c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi - c.lo;
  out.hi =   a.hi*b.lo - c.hi;
  return out;
}

static inline _fjsp_v2r8 _fjsp_nmsub_cp_neg_sr12_v2r8(_fjsp_v2r8 a, _fjsp_v2r8 b, _fjsp_v2r8 c){
  //  -((-hi) a.hi * b^t - c)
  _fjsp_v2r8 out;
  out.lo = - a.hi*b.hi + c.lo;
  out.hi =   a.hi*b.lo + c.hi;
  return out;
}

#else

//
// Fujitsu compiler HPC-ACE intrinsics
//
#include <emmintrin.h>

#endif

#endif

#endif
