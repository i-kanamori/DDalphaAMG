/*
 * Copyright (C) 2018, Matthias Rottmann, Artur Strebel, Simon Heybrock, Simone Bacchio, Bjoern Leder, Issaku Kanamori, Ken-Ichi Ishikawa.
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

#include "main.h"

#ifdef OPTIMIZED_LINALG_float
complex_float local_xy_over_xx_float( vector_float phi, vector_float psi, int start, int end, level_struct *l  ) {
  
  complex_float numerator = 0.0; float denominator = 0.0;

  // TODO rewrite with intrinsics
  if ( l->depth == 0 ) {
    // (end-start) % 12 == 0
    //    for(int i=start; i<endl; i++){
    //      numerator += conj_float(phi[i])*psi[i];
    //      denominator += NORM_SQUARE_float(phi[i])
    //    }
    for (int i=start; i<end; i+=12 ){
      int ii=i;
      FOR12( numerator += conj_float(phi[ii])*psi[ii];
	     denominator += NORM_SQUARE_float(phi[ii]);
	    ; ii++; ) 
	}

  } else {
    // (end-start) % 2 == 0
    //    for(int i=start; i<endl; i++){
    //      numerator += conj_float(phi[i])*psi[i];
    //      denominator += NORM_SQUARE_float(phi[i])
    //    }
    for (int i=start; i<end; i+=2 ){
      int ii=i;
      FOR2( numerator += conj_float(phi[ii])*psi[ii];
	     denominator += NORM_SQUARE_float(phi[ii]);
	    ; ii++; ) 
	}
    
  }
  if ( abs_float(denominator) < EPS_float ) {
    return 0.0;
  }
  
  return numerator/denominator;
}
#endif

void gram_schmidt_on_aggregates_float_vectorized( complex_float *V_dummy, const int num_vect, level_struct *l, struct Thread *threading ) {
  
  PROF_float_START( _GRAM_SCHMIDT_ON_AGGREGATES, threading );
  SYNC_CORES(threading)
  SYNC_HYPERTHREADS(threading)
  int i, j, k, k1, k2, num_aggregates = l->s_float.num_aggregates,
      aggregate_size = l->inner_vector_size / num_aggregates, offset = l->num_lattice_site_var/2;
      
  complex_float alpha1, alpha2;
  vector_float v_pt1, v_pt2;
  float norm1, norm2;
  
  vector_float *V=l->is_float.interpolation;
    
  for ( j=threading->n_thread*threading->core+threading->thread; j<num_aggregates; j+=threading->n_thread*threading->n_core ) {
    for ( k1=0; k1<num_vect; k1++ ) {
      v_pt1 = V[k1] + j*aggregate_size;
      
      for ( k2=0; k2<k1; k2++ ) {
        v_pt2 = V[k2] + j*aggregate_size;
        alpha1 = 0; alpha2 = 0;
        // V[k1] -= <V[k2],V[k1]> V[k2] | 2*j-th and 2*j+1-st aggregate
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            alpha1 += conj_float(v_pt2[i]) * v_pt1[i];
          for ( k=0; k<offset; k++, i++ )
            alpha2 += conj_float(v_pt2[i]) * v_pt1[i];
        }
        for ( i=0; i<aggregate_size; ) {
          for ( k=0; k<offset; k++, i++ )
            v_pt1[i] -=  alpha1 * v_pt2[i];
          for ( k=0; k<offset; k++, i++ )
            v_pt1[i] -=  alpha2 * v_pt2[i];
        }
      }
      
      norm1 = 0; norm2 = 0;
      // V[k1] = V[k1]/norm(V[k1]) | 2*j-th and 2*j+1-st aggregate    
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          norm1 += NORM_SQUARE_float(v_pt1[i]);
        for ( k=0; k<offset; k++, i++ )
          norm2 += NORM_SQUARE_float(v_pt1[i]);
      }
      norm1 = 1/sqrt(norm1); norm2 = 1/sqrt(norm2);
      for ( i=0; i<aggregate_size; ) {
        for ( k=0; k<offset; k++, i++ )
          v_pt1[i] =  norm1 * creal_float(v_pt1[i]) + I*norm1* cimag_float(v_pt1[i]);
        for ( k=0; k<offset; k++, i++ )
          v_pt1[i] =  norm2 * creal_float(v_pt1[i]) + I*norm2* cimag_float(v_pt1[i]);
      }
    }
  }
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  PROF_float_STOP( _GRAM_SCHMIDT_ON_AGGREGATES, 1, threading );
}

#endif
