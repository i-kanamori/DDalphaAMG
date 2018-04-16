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

Hopping from Positive direction: for boundary

requirements
  MU: 0...3, direction 
  
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


#if (MU==0)
  #define prp prp_T
  #define pbp boundary_pbp_T
#endif
#if (MU==1)
  #define prp prp_Z
  #define pbp boundary_pbp_Z
#endif
#if (MU==2)
  #define prp prp_Y
  #define pbp boundary_pbp_Y
#endif
#if (MU==3)
  #define prp prp_X
  #define pbp boundary_pbp_X
#endif

{
  for ( int i=start; i<end; i++ ) {

    // indexing: this site
    int k=ind[i];

    // indexing: next site
    int k_next=neighbor[i];
  
    // point to the input spinor
    const float *pt = phi+24*k_next;

    // pointer to the gauge field
    const float *ptD = D + 72*k + 18*MU;
    
    // pointer to the output
    const float *pt_out = eta + 24*k;

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

    // output spinor
    register ssu3ferm* out_ptr=(ssu3ferm*)pt_out;

    // spin projectin
    prp(in_c0);
    prp(in_c1);
    prp(in_c2);

    // multiply gauge field: to sp1
    cmult(      c0sp1,u00,in_c0sp1);
    accum_cmult(c0sp1,u01,in_c1sp1);
    accum_cmult(c0sp1,u02,in_c2sp1);
    cmult(      c1sp1,u10,in_c0sp1);
    accum_cmult(c1sp1,u11,in_c1sp1);
    accum_cmult(c1sp1,u12,in_c2sp1);
    cmult(      c2sp1,u20,in_c0sp1);
    accum_cmult(c2sp1,u21,in_c1sp1);
    accum_cmult(c2sp1,u22,in_c2sp1);

    // multiply gauge field: to sp2
    cmult(      c0sp2,u00,in_c0sp2);
    accum_cmult(c0sp2,u01,in_c1sp2);
    accum_cmult(c0sp2,u02,in_c2sp2);
    cmult(      c1sp2,u10,in_c0sp2);
    accum_cmult(c1sp2,u11,in_c1sp2);
    accum_cmult(c1sp2,u12,in_c2sp2);
    cmult(      c2sp2,u20,in_c0sp2);
    accum_cmult(c2sp2,u21,in_c1sp2);
    accum_cmult(c2sp2,u22,in_c2sp2);

    // accumulate
    pbp(out_c0,c0);
    pbp(out_c1,c1);
    pbp(out_c2,c2);
    store_ferm(out_ptr,out_);
    
  }  
}

#undef prp
#undef pbp
#undef prp_T
#undef prp_Z
#undef prp_Y
#undef prp_X
#undef prn_T
#undef prn_Z
#undef prn_Y
#undef prn_X
