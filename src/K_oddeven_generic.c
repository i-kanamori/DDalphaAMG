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
 
#include "main.h"
#ifdef K_OPT

#include "K_dirac.h"

static inline void perform_fwd_bwd_subs_PRECISION( vector_PRECISION x, vector_PRECISION b, config_PRECISION L ) {

/*********************************************************************************
* Solves L*(L^H)*x = b for x, i.e., the clover coupling for a single lattice 
* site.
* - vector_PRECISION b: Right hand side.
* - vector_PRECISION x: Solution.
* - config_PRECISION L: Cholesky factor ( lower triangular matrix )
*********************************************************************************/
  
  register int i, j;
  int n;

  for ( n=0; n<2; n++ ) {
    // forward substitution with L
    for ( i=0; i<6; i++ ) {
      x[i] = b[i];
      for ( j=0; j<i; j++ ) {
        x[i] = x[i] - *L * x[j]; L++;
      }
      x[i] = x[i] / *L; L++;
    }
    L -= 21;
    // backward substitution with L^H
    for ( i=5; i>=0; i-- ) {
      for ( j=i+1; j<6; j++ ) {
        x[i] = x[i] - conj_PRECISION(L[(j*(j+1))/2 + i]) * x[j];
      }
      x[i] = x[i] / conj_PRECISION(L[(i*(i+1))/2 + i]);
    }
    x+=6;
    b+=6;
    L+=21;
  }
}


static inline void LLH_multiply_PRECISION( vector_PRECISION y, vector_PRECISION x, config_PRECISION L ) {

/*********************************************************************************
* Applies the clover coupling term to a vector, by multiplying L^H 
* and then L. 
* - vector_PRECISION x: Input vector.
* - vector_PRECISION y: Output vector.
* - config_PRECISION L: Cholesky factor ( lower triangular matrix )
*********************************************************************************/
  
  register int i, j;
  int n;
  complex_PRECISION z[6];
  
  for ( n=0; n<2; n++ ) {
    // z = L^H x
    for ( j=0; j<6; j++ ) { // columns
      for ( i=0; i<j; i++ ) { // rows
        z[i] += conj_PRECISION(*L)*x[j]; L++;
      }
      z[j] = conj_PRECISION(*L)*x[j]; L++;
    }
    L-=21;
    // y = L*z;
    for ( i=0; i<6; i++ ) { // rows
      y[i] = *L * z[0]; L++;
      for ( j=1; j<=i; j++ ) { // columns
        y[i] += *L * z[j]; L++;
      }
    }
    x+=6;
    y+=6;
  }
}


#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION

void hopping_term_PRECISION( vector_PRECISION eta, vector_PRECISION phi, operator_PRECISION_struct *op,
                             const int amount, level_struct *l, struct Thread *threading ) {
  
  int start_even, end_even, start_odd, end_odd;
  compute_core_start_end_custom(0, op->num_even_sites, &start_even, &end_even, l, threading, 1 );
  compute_core_start_end_custom(op->num_even_sites, op->num_even_sites+op->num_odd_sites, &start_odd, &end_odd, l, threading, 1 );

  int i, n = l->num_inner_lattice_sites, *neighbor = op->neighbor_table, *nb_pt,
      start=0, plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  complex_PRECISION pbuf[6];
  vector_PRECISION phi_pt, eta_pt, end_pt;
  config_PRECISION D_pt;
  
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = start_odd, n = end_odd;
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    start = start_even, n = end_even;
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  // project in negative directions
  for ( i=6*start, phi_pt=phi+12*start; i<6*n; i+=6, phi_pt+=12 ) {
    prp_T_PRECISION( op->prnT+i, phi_pt );
    prp_Z_PRECISION( op->prnZ+i, phi_pt );
    prp_Y_PRECISION( op->prnY+i, phi_pt );
    prp_X_PRECISION( op->prnX+i, phi_pt );
  }
  // start communication in negative direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  // project plus dir and multiply with U dagger
  for ( phi_pt=phi+12*start, end_pt=phi+12*n, D_pt = op->D+36*start, nb_pt=neighbor+4*start; phi_pt<end_pt; phi_pt+=12 ) {
    // T dir
    i = 6*(*nb_pt); nb_pt++;
    prn_T_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpT+i, D_pt, pbuf );
    mvmh_PRECISION( op->prpT+i+3, D_pt, pbuf+3 ); D_pt += 9;
    // Z dir
    i = 6*(*nb_pt); nb_pt++;
    prn_Z_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpZ+i, D_pt, pbuf );
    mvmh_PRECISION( op->prpZ+i+3, D_pt, pbuf+3 ); D_pt += 9;
    // Y dir
    i = 6*(*nb_pt); nb_pt++;
    prn_Y_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpY+i, D_pt, pbuf );
    mvmh_PRECISION( op->prpY+i+3, D_pt, pbuf+3 ); D_pt += 9;
    // X dir
    i = 6*(*nb_pt); nb_pt++;
    prn_X_PRECISION( pbuf, phi_pt );
    mvmh_PRECISION( op->prpX+i, D_pt, pbuf );
    mvmh_PRECISION( op->prpX+i+3, D_pt, pbuf+3 ); D_pt += 9;
  }
  if ( amount == _EVEN_SITES ) {
    start = start_even, n = end_even;
  } else if ( amount == _ODD_SITES ) {
    start = start_odd, n = end_odd;
  }  
  // start communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_sendrecv_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
  ghost_sendrecv_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
  // wait for communication in negative direction
  ghost_wait_PRECISION( op->prnT, T, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnZ, Z, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnY, Y, -1, &(op->c), minus_dir_param, l );
  ghost_wait_PRECISION( op->prnX, X, -1, &(op->c), minus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  // multiply with U and lift up minus dir
  for ( eta_pt=eta+12*start, end_pt=eta+12*n, D_pt = op->D+36*start, nb_pt=neighbor+4*start; eta_pt<end_pt; eta_pt+=12 ) {
    // T dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnT+i );
    mvm_PRECISION( pbuf+3, D_pt, op->prnT+i+3 );
    pbp_su3_T_PRECISION( pbuf, eta_pt ); D_pt += 9;
    // Z dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnZ+i );
    mvm_PRECISION( pbuf+3, D_pt, op->prnZ+i+3 );
    pbp_su3_Z_PRECISION( pbuf, eta_pt ); D_pt += 9;
    // Y dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnY+i );
    mvm_PRECISION( pbuf+3, D_pt, op->prnY+i+3 );
    pbp_su3_Y_PRECISION( pbuf, eta_pt ); D_pt += 9;
    // X dir
    i = 6*(*nb_pt); nb_pt++;
    mvm_PRECISION( pbuf, D_pt, op->prnX+i );
    mvm_PRECISION( pbuf+3, D_pt, op->prnX+i+3 );
    pbp_su3_X_PRECISION( pbuf, eta_pt ); D_pt += 9;
  }
  // wait for communication in positive direction
  START_LOCKED_MASTER(threading)
  ghost_wait_PRECISION( op->prpT, T, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpZ, Z, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpY, Y, +1, &(op->c), plus_dir_param, l );
  ghost_wait_PRECISION( op->prpX, X, +1, &(op->c), plus_dir_param, l );
  END_LOCKED_MASTER(threading) 
  // lift up plus dir
  for ( i=6*start, eta_pt=eta+12*start; i<6*n; i+=6, eta_pt+=12 ) {
    pbn_su3_T_PRECISION( op->prpT+i, eta_pt );
    pbn_su3_Z_PRECISION( op->prpZ+i, eta_pt );
    pbn_su3_Y_PRECISION( op->prpY+i, eta_pt );
    pbn_su3_X_PRECISION( op->prpX+i, eta_pt );
  }

  SYNC_CORES(threading)
}

#endif

// ---- block odd even ---------------------------------------------------

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void selfcoupling_cholesky_decomposition_PRECISION( const config_PRECISION , config_double );
#endif

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void schwarz_PRECISION_oddeven_setup( operator_PRECISION_struct *op, level_struct *l ) {
  
  PRECISION *clover_pt = op->clover_vectorized, *oe_clover_pt = op->oe_clover_vectorized;

  int mu, i, d0, c0, b0, a0, d1, c1, b1, a1, t, z, y, x, agg_split[4], block_split[4], block_size[4];

  
  if ( g.csw ) {
    for ( mu=0; mu<4; mu++ ) {
      agg_split[mu] = l->local_lattice[mu]/l->coarsening[mu];
      block_split[mu] = l->coarsening[mu]/l->block_lattice[mu];
      block_size[mu] = l->block_lattice[mu];
    }
    
    for ( d0=0; d0<agg_split[T]; d0++ )
      for ( c0=0; c0<agg_split[Z]; c0++ )
        for ( b0=0; b0<agg_split[Y]; b0++ )
          for ( a0=0; a0<agg_split[X]; a0++ )
            
            for ( d1=d0*block_split[T]; d1<(d0+1)*block_split[T]; d1++ )
              for ( c1=c0*block_split[Z]; c1<(c0+1)*block_split[Z]; c1++ )
                for ( b1=b0*block_split[Y]; b1<(b0+1)*block_split[Y]; b1++ )
                  for ( a1=a0*block_split[X]; a1<(a0+1)*block_split[X]; a1++ ) {
                    
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 0 ) {                      
                              for ( i=0; i<72; i++ )
                                oe_clover_pt[i] = clover_pt[i];
                              clover_pt += 72;
                              oe_clover_pt += 72;

                            }
                          }
                    for ( t=d1*block_size[T]; t<(d1+1)*block_size[T]; t++ )
                      for ( z=c1*block_size[Z]; z<(c1+1)*block_size[Z]; z++ )
                        for ( y=b1*block_size[Y]; y<(b1+1)*block_size[Y]; y++ )
                          for ( x=a1*block_size[X]; x<(a1+1)*block_size[X]; x++ ) {
                            if (((t-d1*block_size[T])+(z-c1*block_size[Z])+
                                (y-b1*block_size[Y])+(x-a1*block_size[X]))%2 == 1 ) {
                              K_site_clover_invert_PRECISION( clover_pt, oe_clover_pt);
                              clover_pt += 72;
                              oe_clover_pt += 72;
                            }
                          }
                  }
  } else {
    vector_PRECISION_copy( op->oe_clover, op->clover, 0, l->inner_vector_size, l );
  }
  op->shift = 4+l->dirac_shift;
}
#endif


#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void block_diag_ee_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)  

#pragma omp master
  {
    PROF_PRECISION_START( _TUNE2 );
  }

  int i, n1 = s->num_block_even_sites;
  PRECISION *clover_vectorized = s->op.oe_clover_vectorized + (start/12)*72;  
  config_PRECISION clover = (g.csw==0.0)?s->op.oe_clover+start:s->op.oe_clover+(start/12)*42;
  vector_PRECISION lphi = phi+start, leta = eta+start;


  // diagonal blocks applied to the even sites of a block
  if ( g.csw ) {
    for ( i=0; i<n1; i++ ) {
      K_site_clover_PRECISION( (PRECISION*)leta, (const PRECISION*)lphi, clover_vectorized );
      leta+=12; lphi+=12; clover_vectorized+=72;
    }
  } else {
    for ( i=0; i<12*n1; i++ )
      leta[i] = lphi[i]*clover[i];
  }

#pragma omp master
  {
    PROF_PRECISION_STOP( _TUNE2, 1 );
  }


  END_UNTHREADED_FUNCTION(threading)
}

#endif

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION

void block_diag_oo_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
  config_PRECISION clover = (g.csw==0.0)?s->op.oe_clover+start:s->op.oe_clover+(start/12)*42;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  // diagonal blocks applied to the odd sites of a block
  if ( g.csw ) {
    error0("block_diag_oo_PRECISION is not available when using K_OPT\n");

    leta += n1*12; lphi += n1*12; clover += n1*42;
    for ( i=0; i<n2; i++ ) {
      LLH_multiply_PRECISION( leta, lphi, clover );
      leta+=12; lphi+=12; clover+=42;
    }
  } else {
    leta += n1*12; lphi += n1*12; clover += n1*12;
    for ( i=0; i<12*n2; i++ )
      leta[i] = lphi[i]*clover[i];
  }

  END_UNTHREADED_FUNCTION(threading)
}

#endif

#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
void block_diag_oo_inv_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)

#pragma omp master
  {
    PROF_PRECISION_START( _TUNE1 );
  }

  
  PRECISION *clover_vectorized = s->op.oe_clover_vectorized + (start/12)*72;
  int i, n1 = s->num_block_even_sites, n2 = s->num_block_odd_sites;
  config_PRECISION clover = (g.csw==0.0)?s->op.oe_clover+start:s->op.oe_clover+(start/12)*42;
  vector_PRECISION lphi = phi+start, leta = eta+start;
  // inverted diagonal blocks applied to the odd sites of a block
  if ( g.csw ) {
    leta += n1*12; lphi += n1*12;
    clover_vectorized += n1*72; // use the inverse
    for ( i=0; i<n2; i++ ) {
      K_site_clover_PRECISION( (PRECISION*)leta, (const PRECISION*)lphi, clover_vectorized );
      leta+=12; lphi+=12; clover_vectorized+=72;
    }
  } else {
    leta += n1*12; lphi += n1*12; clover += n1*12;
    for ( i=0; i<12*n2; i++ )
      leta[i] = lphi[i]/clover[i];
  }

#pragma omp master
  {
    PROF_PRECISION_STOP( _TUNE1,1 );
  }

  END_UNTHREADED_FUNCTION(threading)
}

#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION

void block_hopping_term_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {

#pragma omp master
    {
      PROF_PRECISION_START( _TUNE7 );
    }
  
  START_UNTHREADED_FUNCTION(threading)

  int *length_even = s->dir_length_even, *length_odd = s->dir_length_odd,
    **index = s->oe_index, *neighbor = s->op.neighbor_table;
  config_PRECISION D = s->op.D + (start/12)*36;

  int a1[4], a2[4];
  int n1[4], n2[4];
  for ( int mu=0; mu<4; mu++ ) {
    if ( amount == _EVEN_SITES ) {
      a1[mu] = 0; n1[mu] = length_even[mu];
      a2[mu] = n1[mu]; n2[mu] = a2[mu] + length_odd[mu];
    } else if ( amount == _ODD_SITES ) {
      a1[mu] = length_even[mu]; n1[mu] = a1[mu] + length_odd[mu];
      a2[mu] = 0; n2[mu] = a1[mu];
    } else {
      a1[mu] = 0; n1[mu] = length_even[mu]+length_odd[mu];
      a2[mu] = 0; n2[mu] = n1[mu];
    }
  }

  block_oddeven_coupling_PRECISION( (PRECISION*)(eta+start), 
				    (PRECISION*)D,
				    (PRECISION*)(phi+start),
				    a1,a2,n1,n2,index, neighbor);

  END_UNTHREADED_FUNCTION(threading)

#pragma omp master
    {
      PROF_PRECISION_STOP( _TUNE7,1 );
    }

}

#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
void block_n_hopping_term_PRECISION( vector_PRECISION eta, vector_PRECISION phi,
    int start, int amount, schwarz_PRECISION_struct *s, level_struct *l, struct Thread *threading ) {
  
  START_UNTHREADED_FUNCTION(threading)

  int *length_even = s->dir_length_even, *length_odd = s->dir_length_odd,
      **index = s->oe_index, *neighbor = s->op.neighbor_table;
  config_PRECISION D = s->op.D + (start/12)*36;

  int a1[4], a2[4];
  int n1[4], n2[4];
  for ( int mu=0; mu<4; mu++ ) {
    if ( amount == _EVEN_SITES ) {
      a1[mu] = 0; n1[mu] = length_even[mu];
      a2[mu] = n1[mu]; n2[mu] = a2[mu] + length_odd[mu];
    } else if ( amount == _ODD_SITES ) {
      a1[mu] = length_even[mu]; n1[mu] = a1[mu] + length_odd[mu];
      a2[mu] = 0; n2[mu] = a1[mu];
    } else {
      a1[mu] = 0; n1[mu] = length_even[mu]+length_odd[mu];
      a2[mu] = 0; n2[mu] = n1[mu];
    }
  }

  block_n_oddeven_coupling_PRECISION( (PRECISION*)(eta+start), 
				    (PRECISION*)D,
				    (PRECISION*)(phi+start),
				    a1,a2,n1,n2,index, neighbor);
  
  END_UNTHREADED_FUNCTION(threading)
}

#endif

#endif // K_OPT

