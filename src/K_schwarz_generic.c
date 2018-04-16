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

#include "main.h"

#ifdef K_OPT

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
// from schwarz_genereic.c
void block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi, int k,
                                  schwarz_PRECISION_struct *s, level_struct *l ) {

  // k: number of current block
  int i, mu, index, neighbor_index, *bbl = s->block_boundary_length;
  complex_PRECISION buf1[12], *buf2=buf1+6;
  config_PRECISION D_pt, D = s->op.D;
  vector_PRECISION phi_pt, eta_pt;
  
  mu=T;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_T_PRECISION( buf1, phi_pt );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_T_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_T_PRECISION( buf1, phi_pt );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_T_PRECISION( buf2, eta_pt );
  }
  
  mu=Z;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_Z_PRECISION( buf1, phi_pt );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_Z_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_Z_PRECISION( buf1, phi_pt );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_Z_PRECISION( buf2, eta_pt );
  }
  
  mu=Y;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_Y_PRECISION( buf1, phi_pt );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_Y_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_Y_PRECISION( buf1, phi_pt );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_Y_PRECISION( buf2, eta_pt );
  }
  
  mu=X;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_X_PRECISION( buf1, phi_pt );
    mvm_PRECISION( buf2, D_pt, buf1 );
    mvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_X_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_X_PRECISION( buf1, phi_pt );
    mvmh_PRECISION( buf2, D_pt, buf1 );
    mvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_X_PRECISION( buf2, eta_pt );
  }  
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
// from schwarz_genereic.c
void n_block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi, int k,
                                    schwarz_PRECISION_struct *s, level_struct *l ) {

  // k: number of current block
  int i, mu, index, neighbor_index, *bbl = s->block_boundary_length;
  complex_PRECISION buf1[12], *buf2=buf1+6;
  config_PRECISION D_pt, D = s->op.D;
  vector_PRECISION phi_pt, eta_pt;
  
  mu=T;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_T_PRECISION( buf1, phi_pt );
    nmvm_PRECISION( buf2, D_pt, buf1 );
    nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_T_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_T_PRECISION( buf1, phi_pt );
    nmvmh_PRECISION( buf2, D_pt, buf1 );
    nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_T_PRECISION( buf2, eta_pt );
  }
  
  mu=Z;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_Z_PRECISION( buf1, phi_pt );
    nmvm_PRECISION( buf2, D_pt, buf1 );
    nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_Z_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_Z_PRECISION( buf1, phi_pt );
    nmvmh_PRECISION( buf2, D_pt, buf1 );
    nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_Z_PRECISION( buf2, eta_pt );
  }
  
  mu=Y;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_Y_PRECISION( buf1, phi_pt );
    nmvm_PRECISION( buf2, D_pt, buf1 );
    nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_Y_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_Y_PRECISION( buf1, phi_pt );
    nmvmh_PRECISION( buf2, D_pt, buf1 );
    nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_Y_PRECISION( buf2, eta_pt );
  }
  
  mu=X;
  // plus mu direction
  for ( i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prp_X_PRECISION( buf1, phi_pt );
    nmvm_PRECISION( buf2, D_pt, buf1 );
    nmvm_PRECISION( buf2+3, D_pt, buf1+3 );
    pbp_su3_X_PRECISION( buf2, eta_pt );
  }
  // minus mu direction
  for ( i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
    index = s->block[k].bt[i];
    neighbor_index = s->block[k].bt[i+1];
    D_pt = D + 36*neighbor_index + 9*mu;
    phi_pt = phi + 12*neighbor_index;
    eta_pt = eta + 12*index;
    prn_X_PRECISION( buf1, phi_pt );
    nmvmh_PRECISION( buf2, D_pt, buf1 );
    nmvmh_PRECISION( buf2+3, D_pt, buf1+3 );
    pbn_su3_X_PRECISION( buf2, eta_pt );
  }  
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
// from schwarz_generic.c
void coarse_block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi,
                                         int k, schwarz_PRECISION_struct *s, level_struct *l ) {

  // k: number of current block
  int *bbl = s->block_boundary_length, n = l->num_lattice_site_var;
  config_PRECISION D = s->op.D;
  int link_size = SQUARE(l->num_lattice_site_var), site_size=4*link_size;
  
  for ( int mu=0; mu<4; mu++ ) {
    // plus mu direction
    for ( int i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      config_PRECISION D_pt = D + site_size*index + link_size*mu;
      coarse_hopp_PRECISION( eta_pt, phi_pt, D_pt, l );
    }
    // minus mu direction
    for ( int i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      config_PRECISION D_pt = D + site_size*neighbor_index + link_size*mu;
      coarse_daggered_hopp_PRECISION( eta_pt, phi_pt, D_pt, l );
    }
  }
}
#endif

#ifdef OPTIMIZED_NEIGHBOR_COUPLING_PRECISION
// from schwarz_generic.c
void n_coarse_block_PRECISION_boundary_op( vector_PRECISION eta, vector_PRECISION phi,
                                           int k, schwarz_PRECISION_struct *s, level_struct *l ) {

  // k: number of current block
  int *bbl = s->block_boundary_length, n = l->num_lattice_site_var;
  int link_size = SQUARE(l->num_lattice_site_var), site_size=4*link_size;
  config_PRECISION D = s->op.D;
  
  for ( int mu=0; mu<4; mu++ ) {
    // plus mu direction
    for ( int i=bbl[2*mu]; i<bbl[2*mu+1]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      config_PRECISION D_pt = D + site_size*index + link_size*mu;
      coarse_n_hopp_PRECISION( eta_pt, phi_pt, D_pt, l );
    }
    // minus mu direction
    for ( int i=bbl[2*mu+1]; i<bbl[2*mu+2]; i+=2 ) {
      int index = s->block[k].bt[i];
      int neighbor_index = s->block[k].bt[i+1];
      vector_PRECISION phi_pt = phi + n*neighbor_index;
      vector_PRECISION eta_pt = eta + n*index;
      config_PRECISION D_pt = D + site_size*neighbor_index + link_size*mu;
      coarse_n_daggered_hopp_PRECISION( eta_pt, phi_pt, D_pt, l );
    }
  }
}
#endif

#if defined(OPTIMIZED_NEIGHBOR_COUPLING_PRECISION) || defined(OPTIMIZED_SELF_COUPLING_PRECISION)
// from schwarz_generic.c

void schwarz_PRECISION_setup( schwarz_PRECISION_struct *s, operator_double_struct *op_in, level_struct *l ) {

/*********************************************************************************  
* Copies the Dirac operator and the clover term from op_in into the Schwarz 
* struct (this function is depth 0 only).
* - operator_double_struct *op_in: Input operator.                                  
*********************************************************************************/

  int i, index, n = l->num_inner_lattice_sites, *tt = s->op.translation_table;
  config_PRECISION D_out_pt, clover_out_pt;
  config_double D_in_pt = op_in->D, clover_in_pt = op_in->clover;
  s->op.shift = op_in->shift;
  
  // copy as it is
  for ( i=0; i<n; i++ ) {
    index = tt[i];
    D_out_pt = s->op.D + 36*index;
    FOR36( *D_out_pt = (complex_PRECISION) *D_in_pt; D_out_pt++; D_in_pt++; )
  }
  

  if ( g.csw != 0 ) {
    for ( i=0; i<n; i++ ) {
      index = tt[i];
      clover_out_pt = s->op.clover + 42*index;
#ifdef OPTIMIZED_SELF_COUPLING_PRECISION
      PRECISION *clover_out_vectorized_pt = s->op.clover_vectorized + 72*index;
      K_set_clover_PRECISION( clover_out_vectorized_pt, clover_in_pt );
#endif
      FOR42( *clover_out_pt = (complex_PRECISION) *clover_in_pt; clover_out_pt++; clover_in_pt++; )
    }
  } else {
    for ( i=0; i<n; i++ ) {
      index = tt[i];
      clover_out_pt = s->op.clover + 12*index;
      FOR12( *clover_out_pt = (complex_PRECISION) *clover_in_pt; clover_out_pt++; clover_in_pt++; )
    }
  }
  
  if ( g.odd_even )
    schwarz_PRECISION_oddeven_setup( &(s->op), l );
  
  schwarz_PRECISION_boundary_update( s, l );
}
#endif

#endif // K_OPT
