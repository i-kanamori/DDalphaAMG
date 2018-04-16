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

/*

Hopping from Negative direction

requirements
  MU: 0...3, direction 
  
  BOUNDARY:  if defined, special tremtment for indexing
  NEGATIVE_COUPLING: if defined, accmulate with extra (-1)
  
 */

#ifdef BASIS2
#if defined(SPIN12)
#include "K_spinproj_basis2_spin12.h"
#elif defined(SPIN34)
#include "K_spinproj_basis2_spin34.h"
#else
#include "K_spinproj_basis2.h"
#endif
#else
#error wrong cliford basis: only BASIS2 is avilable on K
#endif

#ifdef ACCUM_PLUS
#if (MU==0)
#define prn prn_T
#define pbn plus_pbn_T
#endif
#if (MU==1)
 #define prn prn_Z
 #define pbn plus_pbn_Z
#endif
#if (MU==2)
 #define prn prn_Y
 #define pbn plus_pbn_Y
#endif
#if (MU==3)
 #define prn prn_X
 #define pbn plus_pbn_X
#endif

#else

#if (MU==0)
#define prn prn_T
#define pbn pbn_T
#endif
#if (MU==1)
 #define prn prn_Z
 #define pbn pbn_Z
#endif
#if (MU==2)
 #define prn prn_Y
 #define pbn pbn_Y
#endif
#if (MU==3)
 #define prn prn_X
 #define pbn pbn_X
#endif

#endif

{
#ifdef BOUNDARY  
  for ( int i=start; i<end; i+=2 ) {
#else
  for ( int i=start; i<end; i++ ) {
#endif

    // indexing: this site
    int k=ind[i];

    // indexing: next site
#ifdef BOUNDARY
    int k_next=ind[i+1];
#else
    int k_next=neighbor[4*k+MU];
#endif
  
    // point to the input spinor
    const float *pt = phi+24*k;

    // pointer to the gauge field
    const float *ptD = D + 72*k + 18*MU;
    
    // pointer to the output
    const float *pt_out = eta +24*k_next;

    // input spinor
    register _fjsp_v2r8 in_c0sp1, in_c0sp2, in_c0sp3, in_c0sp4;
    register _fjsp_v2r8 in_c1sp1, in_c1sp2, in_c1sp3, in_c1sp4;
    register _fjsp_v2r8 in_c2sp1, in_c2sp2, in_c2sp3, in_c2sp4;

    // output spinor
    register _fjsp_v2r8 out_c0sp1, out_c0sp2, out_c0sp3, out_c0sp4;
    register _fjsp_v2r8 out_c1sp1, out_c1sp2, out_c1sp3, out_c1sp4;
    register _fjsp_v2r8 out_c2sp1, out_c2sp2, out_c2sp3, out_c2sp4;

    // gauge field
    register _fjsp_v2r8 u00,u01,u02;
    register _fjsp_v2r8 u10,u11,u12;
    register _fjsp_v2r8 u20,u21,u22;

    // work vector
    register _fjsp_v2r8 c0sp1, c1sp1, c2sp1;
    register _fjsp_v2r8 c0sp2, c1sp2, c2sp2;

    // load input spinor
    register ssu3ferm* in_ptr=(ssu3ferm*)pt;
    load_ferm(in_,in_ptr);
      
    // load gauge
    register ssu3gauge* u_ptr=(ssu3gauge*)ptD;
    load_gauge(u, u_ptr);

    // load output spinor
    register ssu3ferm* out_ptr=(ssu3ferm*)pt_out;
    load_ferm(out_,out_ptr);
      
    // spin projectin
    prn(in_c0);
    prn(in_c1);
    prn(in_c2);

    // multiply gauge field: to sp1
    conj_cmult(      c0sp1,u00,in_c0sp1);
    accum_conj_cmult(c0sp1,u10,in_c1sp1);
    accum_conj_cmult(c0sp1,u20,in_c2sp1);
    conj_cmult(      c1sp1,u01,in_c0sp1);
    accum_conj_cmult(c1sp1,u11,in_c1sp1);
    accum_conj_cmult(c1sp1,u21,in_c2sp1);
    conj_cmult(      c2sp1,u02,in_c0sp1);
    accum_conj_cmult(c2sp1,u12,in_c1sp1);
    accum_conj_cmult(c2sp1,u22,in_c2sp1);

    // multiply gauge field: to sp2
    conj_cmult(      c0sp2,u00,in_c0sp2);
    accum_conj_cmult(c0sp2,u10,in_c1sp2);
    accum_conj_cmult(c0sp2,u20,in_c2sp2);
    conj_cmult(      c1sp2,u01,in_c0sp2);
    accum_conj_cmult(c1sp2,u11,in_c1sp2);
    accum_conj_cmult(c1sp2,u21,in_c2sp2);
    conj_cmult(      c2sp2,u02,in_c0sp2);
    accum_conj_cmult(c2sp2,u12,in_c1sp2);
    accum_conj_cmult(c2sp2,u22,in_c2sp2);

    // accumulate
    pbn(out_c0,c0);
    pbn(out_c1,c1);
    pbn(out_c2,c2);
    store_ferm(out_ptr,out_);
    
  }  
}

#undef prn
#undef pbn
#undef prp_T
#undef prp_Z
#undef prp_Y
#undef prp_X
#undef prn_T
#undef prn_Z
#undef prn_Y
#undef prn_X
