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

#ifdef K_OPT

///////////////////////////////////////////////////////
//
// T direction:
//
/* gamma_T =
 *  0  0  1  0
 *  0  0  0  1
 *  1  0  0  0
 *  0  1  0  0  
*/
///////////////////////////////////////////////////////

// 1 - gamma_T
// -sp3
// -sp4
#define prp_T(psi){			    \
  psi ## sp1 = _fjsp_neg_v2r8(psi ## sp3 ); \
  psi ## sp2 = _fjsp_neg_v2r8(psi ## sp4 ); \
} 

// -(1 - gamma_T)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3+sp1
// out4 = out4+sp2
#define pbp_T(out, psi){			       \
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_add_v2r8(out ## sp3, psi ## sp1 ); \
  out ## sp4 = _fjsp_add_v2r8(out ## sp4, psi ## sp2 ); \
}

// +(1 - gamma_T)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3-sp1
// out4 = out4-sp2
#define plus_pbp_T(out, psi){			       \
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_sub_v2r8(out ## sp3, psi ## sp1 ); \
  out ## sp4 = _fjsp_sub_v2r8(out ## sp4, psi ## sp2 ); \
}


// +(1 - gamma_T)
// out1 = sp1
// out2 = sp2
// out3 = -sp1
// out4 = -sp2
#define boundary_pbp_T(out, psi){		\
  out ## sp1 =  psi ## sp1;      		\
  out ## sp2 =  psi ## sp2;			\
  out ## sp3 = _fjsp_neg_v2r8(psi ## sp1 );	\
  out ## sp4 = _fjsp_neg_v2r8(psi ## sp2 );	\
}

// 1 + gamma_T
// sp3
// sp4
#define prn_T(psi){\
  psi ## sp1 = psi ## sp3; \
  psi ## sp2 = psi ## sp4; \
} 

// -(1 + gamma_T)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3-sp1
// out4 = out4-sp2
#define pbn_T(out, psi){			       \
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_sub_v2r8(out ## sp3, psi ## sp1 ); \
  out ## sp4 = _fjsp_sub_v2r8(out ## sp4, psi ## sp2 ); \
}

// +(1 + gamma_T)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3+sp1
// out4 = out4+sp2
#define plus_pbn_T(out, psi){			       \
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_add_v2r8(out ## sp3, psi ## sp1 ); \
  out ## sp4 = _fjsp_add_v2r8(out ## sp4, psi ## sp2 ); \
}



///////////////////////////////////////////////////////
//
// Z direction:
//
/* gamma_Z =
 *  0  0  I  0
 *  0  0  0 -I
 * -I  0  0  0
 *  0  I  0  0  
*/
///////////////////////////////////////////////////////

// 1 - gamma_Z
// - I sp3
// + I sp4
#define prp_Z(psi){\
  psi ## sp1 = times_minusI(psi ## sp3);	\
  psi ## sp2 = times_I(psi ## sp4);		\
} 

// -(1 - gamma_Z)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3- I sp1
// out4 = out4+ I sp2
#define pbp_Z(out, psi){				\
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 );	\
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_minusI(out ## sp3, psi## sp1 );		\
  accum_times_I(     out ## sp4, psi## sp2 );		\
}

// +(1 - gamma_Z)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3+ I sp1
// out4 = out4- I sp2
#define plus_pbp_Z(out, psi){				\
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 );	\
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_I(      out ## sp3, psi## sp1 );		\
  accum_times_minusI( out ## sp4, psi## sp2 );		\
}

// +(1 - gamma_Z)
// out1 = sp1
// out2 = sp2
// out3 = +I sp1
// out4 = -I sp2
#define boundary_pbp_Z(out, psi){	 \
  out ## sp1 = psi ## sp1;		 \
  out ## sp2 = psi ## sp2;		 \
  out ## sp3 = times_I(psi ## sp1);	 \
  out ## sp4 = times_minusI(psi ## sp2); \
}

// 1 + gamma_Z
// + I sp3
// - I sp4
#define prn_Z(psi){			\
  psi ## sp1 = times_I(psi ## sp3);	\
  psi ## sp2 = times_minusI(psi ## sp4);\
} 

// -(1 + gamma_Z)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3+ I sp1
// out4 = out4- I sp2
#define pbn_Z(out, psi){				\
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 );	\
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_I(     out ## sp3, psi ## sp1 );		\
  accum_times_minusI(out ## sp4, psi ## sp2 );		\
}

// +(1 + gamma_Z)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3- I sp1
// out4 = out4+ I sp2
#define plus_pbn_Z(out, psi){				\
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 );	\
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_minusI( out ## sp3, psi ## sp1 );		\
  accum_times_I(      out ## sp4, psi ## sp2 );		\
}



///////////////////////////////////////////////////////
//
// Y direction:
//
/* gamma_Y =
 *  0  0  0 -1
 *  0  0  1  0
 *  0  1  0  0
 * -1  0  0  0  
*/
///////////////////////////////////////////////////////
// 1 - gamma_Y
// +sp4
// -sp3
#define prp_Y(psi){				\
  psi ## sp1 =  psi ## sp4;			\
  psi ## sp2 =  _fjsp_neg_v2r8(psi ## sp3 );	\
}

// -(1 - gamma_Y)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3+sp2
// out4 = out4-sp1
#define pbp_Y(out, psi){				\
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_add_v2r8(out ## sp3, psi ## sp2 ); \
  out ## sp4 = _fjsp_sub_v2r8(out ## sp4, psi ## sp1 ); \
}

// +(1 - gamma_Y)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3-sp2
// out4 = out4+sp1
#define plus_pbp_Y(out, psi){				\
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_sub_v2r8(out ## sp3, psi ## sp2 ); \
  out ## sp4 = _fjsp_add_v2r8(out ## sp4, psi ## sp1 ); \
}

// +(1 - gamma_Y)
// out1 = sp1
// out2 = sp2
// out3 = -sp2
// out4 = +sp1
#define boundary_pbp_Y(out, psi){		\
  out ## sp1 = psi ## sp1 ;			\
  out ## sp2 = psi ## sp2 ;			\
  out ## sp3 = _fjsp_neg_v2r8(psi ## sp2 );	\
  out ## sp4 = psi ## sp1 ;			\
}

// 1 + gamma_Y
// -sp4
// +sp3
#define prn_Y(psi){				\
  psi ## sp1 = _fjsp_neg_v2r8(psi ## sp4 );	\
  psi ## sp2 = psi ## sp3;			\
}

// -(1 + gamma_Y)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3-sp2
// out4 = out4+sp1
#define pbn_Y(out, psi){				\
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_sub_v2r8(out ## sp3, psi ## sp2 ); \
  out ## sp4 = _fjsp_add_v2r8(out ## sp4, psi ## sp1 ); \
}

// +(1 + gamma_Y)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3+sp2
// out4 = out4-sp1
#define plus_pbn_Y(out, psi){				\
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 ); \
  out ## sp3 = _fjsp_add_v2r8(out ## sp3, psi ## sp2 ); \
  out ## sp4 = _fjsp_sub_v2r8(out ## sp4, psi ## sp1 ); \
}


///////////////////////////////////////////////////////
//
// X direction:
//
/* gamma_X =
 *  0  0  0  I
 *  0  0  I  0
 *  0 -I  0  0
 * -I  0  0  0  
*/
///////////////////////////////////////////////////////
// 1 - gamma_X
// - I sp4
// - I sp3
#define prp_X(psi){				\
  psi ## sp1 = times_minusI(psi ## sp4 );	\
  psi ## sp2 = times_minusI(psi ## sp3 );	\
}

// -(1 - gamma_X)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3- I sp2
// out4 = out4- I sp1
#define pbp_X(out, psi){				\
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_minusI(out ## sp3, psi ## sp2 ); \
  accum_times_minusI(out ## sp4, psi ## sp1 ); \
}

// +(1 - gamma_X)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3+ I sp2
// out4 = out4+ I sp1
#define plus_pbp_X(out, psi){				\
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_I( out ## sp3, psi ## sp2 );		\
  accum_times_I( out ## sp4, psi ## sp1 );		\
}

// +(1 - gamma_X)
// out1 = sp1
// out2 = sp2
// out3 =  I sp2
// out4 =  I sp1
#define boundary_pbp_X(out, psi){	\
  out ## sp1 = psi ## sp1 ;		\
  out ## sp2 = psi ## sp2 ;		\
  out ## sp3 = times_I(psi ## sp2);	\
  out ## sp4 = times_I(psi ## sp1);	\
}

// 1 + gamma_X
// + I sp4
// + I sp3
#define prn_X(psi){			\
  psi ## sp1 = times_I( psi ## sp4 );	\
  psi ## sp2 = times_I( psi ## sp3 );	\
}

// -(1 + gamma_X)
// out1 = out1-sp1
// out2 = out2-sp2
// out3 = out3+ I sp2
// out4 = out4+ I sp1
#define pbn_X(out, psi){				\
  out ## sp1 = _fjsp_sub_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_sub_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_I(out ## sp3, psi ## sp2 );		\
  accum_times_I(out ## sp4, psi ## sp1 );		\
}

// +(1 + gamma_X)
// out1 = out1+sp1
// out2 = out2+sp2
// out3 = out3- I sp2
// out4 = out4- I sp1
#define plus_pbn_X(out, psi){				\
  out ## sp1 = _fjsp_add_v2r8(out ## sp1, psi ## sp1 ); \
  out ## sp2 = _fjsp_add_v2r8(out ## sp2, psi ## sp2 );	\
  accum_times_minusI(out ## sp3, psi ## sp2 );	\
  accum_times_minusI(out ## sp4, psi ## sp1 );	\
}

#endif
