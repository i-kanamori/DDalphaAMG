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

void copy_coarse_operator_clover_to_vectorized_layout_double(config_double clover,
							     OPERATOR_TYPE_double *clover_vectorized,
							     int num_aggregates, int num_eig_vect){}




void copy_coarse_operator_clover_block4x4(const complex_float *clover,
					  float *clover_vectorized,
					  int num_aggregates, int num_eig_vect)
{
  // in: clover on a site wtth a default packing out: clover on a site
  // with an optimzied packing
  //
  //  in = (A, D, B)
  //     with  [ A B ]
  //           [ C D ],  C=-B^dag, A=A^dag, D=D^dag
  //       A: upper triangle, column order
  //       D: upper triangle, column order
  //       B: column order
  //  column order: B_{ij} = B[i+n*j]
  // 
  //  out (D1, D2,...., U1,U2,...)
  //  
  //  [ D1 U1 U2 ... U? ]
  //  [    D2 U? ... U? ]  = [ A      B]
  //  [ ...             ]    [ B^dag -D]
  //  [              D? ]
  //
  //  D? and U? are 4x4 matrix
  //
  
  int n = num_eig_vect;
  int N=2*n;
  //assert(N % 4 ==0);

  static const int block_size=4; 
  int block_num=N/block_size;


  // offset between blocks in clover
  int offset_to_D = (n*n+n)/2; // upper triangle of A including diagonal
  int offset_to_B = 2*offset_to_D; // B comes after A and D
 
  float *work=NULL;
  const int work_size=2*N*N;
 
  MALLOC( work, float, work_size );

  for ( int a=0; a<num_aggregates; a++ ) {
    const complex_float *in=clover+(offset_to_B + n*n)*a;
    const float *in_real = (float*)in;
    float *out=clover_vectorized+N*N*a;
   
    // once reconstruct the whole matrix
    int ii=0; // index to read

    // A
    for(int j=0; j<n; j++){
#pragma loop noalias
      for(int i=0; i<j+1; i++){
	int ij=i*N+j;
	int ji=j*N+i;
	work[2*ij] = in_real[ii];
	work[2*ji] = in_real[ii];
	ii++;
	work[2*ij+1] = in_real[ii];
	work[2*ji+1] = -in_real[ii];
	ii++;
      }
    }
    
   // D as -D
   for(int j=n; j<N; j++){
#pragma loop noalias
     for(int i=n; i<j+1; i++){
       int ij=i*N+j;
       int ji=j*N+i;
       work[2*ij] = -in_real[ii];
       work[2*ji] = -in_real[ii];
       ii++;
       work[2*ij+1] = -in_real[ii];
       work[2*ji+1] = +in_real[ii];
       ii++;
     }
   }
   
  
   // B
   for(int j=n; j<N; j++){
#pragma loop noalias
     for(int i=0; i<n; i++){
       int ij=i*N+j;
       int ji=j*N+i;
       work[2*ij] = in_real[ii];
       work[2*ji] = in_real[ii];
       ii++;
       work[2*ij+1] = in_real[ii];
       work[2*ji+1] = -in_real[ii];
       ii++;
     }
   }

   // now the whole matrix is ready
   full_to_block4x4_Hermite(out, (const float*)work, block_num);

  }// a-loop
  
  FREE( work, float, work_size );
}


void copy_coarse_operator_clover_to_vectorized_layout_float(const complex_float *clover,
							    float *clover_vectorized,
							    int num_aggregates, int num_eig_vect)
{
  copy_coarse_operator_clover_block4x4(clover, clover_vectorized, num_aggregates,num_eig_vect);
}



void K_coarse_sc_inverse(const float *clover_in, float *clover_out, int n, float *work){

  //  in= [ A B ]
  //      [ C D ],  C=B^dag, A=A^dag, D=D^dag
  //format:
  //  ([daig, real part only], [offdiag A], [offdiag D], B ) 
  //
  //  for 
  //   in^{-1} = [ A' B' ]
  //             [ C' D' ] with C'=B'^dag
  //
  //  out = [A'^dag C'^dag ]
  //        [B'^dag D'^dag ]
  //

  int N=2*n;
  static const int block_size=4;
  int block_num=N/block_size;

  // working buffer
  complex_float *A    =(complex_float*)(work + 0    );
  complex_float *A_inv=(complex_float*)(work + 2*N*N);
  complex_float *b    =(complex_float*)(work + 4*N*N);

  block4x4_to_full_Hermite((float*)A, clover_in, block_num);
  // full matrix is ready

  // LU decomp in A
  for ( int k=0; k<N-1; k++ ) {
    for ( int i=k+1; i<N; i++ ) {
      // alpha = A_ik/A_kk
      complex_float alpha = A[i+k*N]/A[k+k*N];
      A[i+k*N] = alpha;
      for (int j=k+1; j<N; j++ ) {
	// A_ij = A_ij - alpha * A_kj
	A[i+j*N] -= alpha* A[k+j*N];
      }
    }   
  }
    
  complex_float *x;
  
  for (int k=0; k<N; k++ ) {
    b[k] = 0;
  }
    
  for (int k=0; k<N; k++ ) {
    x = A_inv+k*N;
    b[k] = 1;
    if ( k>0 )
      b[k-1] = 0;
    
    for (int i=0; i<N; i++ ) {
      x[i] = b[i];
      for (int j=0; j<i; j++ ) {
	// x_i = x_i - A_ij + x_j
	x[i] = x[i] - A[i+j*N]*x[j];
      }
    } // i
    
    for (int i=N-1; i>=0; i-- ) {
      for (int j=i+1; j<N; j++ ) {
	// x_i = x_i - A_ij * x_j
	x[i] = x[i] - A[i+j*N]*x[j];
      }
      // x_i = x_i / A_ii
      x[i] = x[i]/A[i+i*N];
    } // i
  } // k

  // A_inv is ready!

  //  [A' B']  to [ A' -B']
  //  [C' D']     [-C'  D']
  for(int i=0; i<n; i++){
    for(int j=n; j<N; j++){
      int ij=i*N+j;
      int ji=j*N+i;
      A_inv[ij] = - A_inv[ij];
      A_inv[ji] = - A_inv[ji];
    }
  }
  full_to_block4x4_Hermite(clover_out, (const float*)A_inv, block_num);

}


void  copy_coarse_operator_to_vectorized_layout_float(const complex_float *matrix_in,
						      float *matrix_out,
						      int n_per_core,
						      int nvec){
  // assumes nvec is even
  // input
  // [ A B ]
  // [ C D ]  each in column order
  // 
  //  A = [ a00 a01 a02..]
  //      [ a10 a22...   ]
  //      [...           ]  with 2x2 a?? 
  // out =
  // [ (a00 b00 c00 d00) (a01 b01 c01 d01)...] each in row order
  //   

  static const int blocksize=4;
  int bs2=blocksize/2;
  int numblock=2*nvec/blocksize;
  int n=nvec;
  int N=2*n;

  for ( int a=0; a<4*n_per_core; a++ ) { // 4 is for 4 directions
    const complex_float *in=matrix_in+N*N*a;
    complex_float *out=(complex_float*)(matrix_out+2*N*N*a); // 2 is for complex
    int offset2=0;

    // A
    for(int j=0; j<n; j++){
      for(int i=0; i<n; i++){
	int ij=i+ n*j;
	int ib=i/2;
	int i0=i % 2;
	int jb=j/2;
	int j0=j % 2;
	int iijj=i0*bs2+j0;
	int offset= (ib*numblock+jb)*blocksize*blocksize;
	out[offset+iijj]=in[ij];
      }}
    in+=n*n;
    offset2=2*bs2*bs2;
    // C
    for(int j=0; j<n; j++){
      for(int i=0; i<n; i++){
	int ij=i+ n*j;
	int ib=i/2;
	int i0=i % 2;
	int jb=j/2;
	int j0=j % 2;
	int iijj=i0*bs2+j0 + offset2;
	int offset= (ib*numblock+jb)*blocksize*blocksize;
	out[offset+iijj]=in[ij];
      }}
    in+=n*n;
    offset2=bs2*bs2;
    // B
    for(int j=0; j<n; j++){
      for(int i=0; i<n; i++){
	int ij=i+ n*j;
	int ib=i/2;
	int i0=i % 2;
	int jb=j/2;
	int j0=j % 2;
	int iijj=i0*bs2+j0 + offset2;
	int offset= (ib*numblock+jb)*blocksize*blocksize;
	out[offset+iijj]=in[ij];
      }}
    in+=n*n;
    offset2=3*bs2*bs2;
    // D
    for(int j=0; j<n; j++){
      for(int i=0; i<n; i++){
	int ij=i+ n*j;
	int ib=i/2;
	int i0=i % 2;
	int jb=j/2;
	int j0=j % 2;
	int iijj=i0*bs2+j0 + offset2;
	int offset= (ib*numblock+jb)*blocksize*blocksize;
	out[offset+iijj]=in[ij];
      }}
  }//a    
}


#ifdef VECTORIZE_COARSE_OPERATOR_float
void coarse_block_operator_float( vector_float eta, vector_float phi, int start, schwarz_float_struct *s, level_struct *l, struct Thread *threading ) {

  START_UNTHREADED_FUNCTION(threading)

  int n = s->num_block_sites, *length = s->dir_length, **index = s->index,
      *ind, *neighbor = s->op.neighbor_table, m = l->num_lattice_site_var;
  vector_float lphi = phi+start, leta = eta+start;
  int hopp_size = 4 * SQUARE( l->num_lattice_site_var );
  config_float D_pt, D = s->op.D + (start/m)*hopp_size;
  
  // site-wise self coupling
  int clov_size = ( (l->num_lattice_site_var*(l->num_lattice_site_var+1))/2 );
  config_float clover = s->op.clover + (start/m)*clov_size;
  coarse_self_couplings_float( leta, lphi, clover, n*m, l );
  
  // inner block couplings
  for ( int mu=0; mu<4; mu++ ) {
    ind = index[mu]; // mu direction
    for ( int i=0; i<length[mu]; i++ ) {
      int k = ind[i]; int j = neighbor[5*k+mu+1];
      D_pt = D + hopp_size*k + (hopp_size/4)*mu;
      coarse_hopp_float( leta+m*k, lphi+m*j, D_pt, l );
      coarse_daggered_hopp_float( leta+m*j, lphi+m*k, D_pt, l );
    }
  }

  END_UNTHREADED_FUNCTION(threading)
}
#endif


#ifdef VECTORIZE_COARSE_OPERATOR_float
void apply_coarse_operator_float( vector_float eta, vector_float phi, operator_float_struct *op,
                                      level_struct *l, struct Thread *threading ) {
  
  PROF_float_START( _SC, threading );
  START_LOCKED_MASTER(threading)
  coarse_self_couplings_float( eta, phi, op->clover, l->inner_vector_size, l );
  END_LOCKED_MASTER(threading)
  PROF_float_STOP( _SC, 1, threading );
  PROF_float_START( _NC, threading );
  coarse_hopping_term_float( eta, phi, op, _FULL_SYSTEM, l, threading );
  PROF_float_STOP( _NC, 1, threading );
}
#endif

#ifdef VECTORIZE_COARSE_OPERATOR_float
void coarse_operator_float_set_couplings( operator_float_struct *op, level_struct *l, struct Thread *threading ) { 

  int n = l->num_inner_lattice_sites;
  int sc_size = (l->num_lattice_site_var/2)*(l->num_lattice_site_var+1);
  int nc_size = SQUARE(l->num_lattice_site_var);
  int n_site_var=l->num_lattice_site_var;
  int nvecs=n_site_var/2;

  int n1, n2;
  if ( l->depth > 0 ) {
    n1 = l->num_lattice_sites;
    n2 = 2*l->num_lattice_sites-l->num_inner_lattice_sites;
  } else {
    n1 = l->num_inner_lattice_sites;
    n2 = l->num_inner_lattice_sites;
  }
    
  START_LOCKED_MASTER(threading)
  if( op->D_vectorized == NULL ) {
    MALLOC_HUGEPAGES( op->D_vectorized, OPERATOR_TYPE_float, 2*4*n_site_var*n_site_var*n2, 64 );
    MALLOC_HUGEPAGES( op->clover_vectorized, OPERATOR_TYPE_float, n_site_var*n_site_var*n, 64 );
  }
  END_LOCKED_MASTER(threading)

  int start, end;
  compute_core_start_end_custom(0, n, &start, &end, l, threading, 1);
  int n_per_core = end-start;
  int offset_v = n_site_var*n_site_var;

  // self coupling: store
  copy_coarse_operator_clover_to_vectorized_layout_float(
      op->clover + start*sc_size,
      op->clover_vectorized + start*offset_v,
      n_per_core, nvecs);
  SYNC_CORES(threading);

  // neighborhood coupling: store
  copy_coarse_operator_to_vectorized_layout_float(
      op->D + 4*start*nc_size,              // D is complex_float*
      op->D_vectorized + 2*4*start*nc_size, // 2 is for complex: D_vectorized is float* 
      n_per_core, nvecs);

  SYNC_CORES(threading)
  
  // vectorize negative boundary
  if ( l->depth > 0 ) {
    compute_core_start_end_custom(n1, n2, &start, &end, l, threading, 1);
    n_per_core = end-start;

  // neighborhood coupling: store
  copy_coarse_operator_to_vectorized_layout_float(
      op->D + 4*start*nc_size,              // D is complex_float*
      op->D_vectorized + 2*4*start*nc_size, // 2 is for complex: D_vectorized is float* 
      n_per_core, nvecs);
    SYNC_CORES(threading)
  }


}
#endif


#ifdef VECTORIZE_COARSE_OPERATOR_double
void coarse_oddeven_setup_double_set_couplings( operator_PRECISION_struct *in, int reorder, level_struct *l, struct Thread *threading ) {
  printf0("!!!!! double precision is not ready: %s in %s\n",__func__, __FILE__);
  abort();
}
#endif


#ifdef VECTORIZE_COARSE_OPERATOR_float
void coarse_oddeven_setup_float_set_couplings( operator_float_struct *in, 
						 int reorder, level_struct *l, 
						 struct Thread *threading ) {

  int i, j;
  int n=l->num_inner_lattice_sites;
  int sc_size = (l->num_lattice_site_var/2)*(l->num_lattice_site_var+1);
  int nc_size = SQUARE(l->num_lattice_site_var);
  int t, z, y, x;
  operator_float_struct *op = &(l->oe_op_float);
  config_float sc_in = in->clover;
  config_float nc_in = in->D;
  config_float Aee = NULL, Aoo = NULL;
  int *le = l->local_lattice;
  int oe_offset = op->oe_offset;

  if(reorder) {
    printf0("hoge: reorder is set, %s in %s\n",__func__,__FILE__);
    fflush(0);
    abort();
  }

  Aee = op->clover;
  Aoo = op->clover + op->num_even_sites*sc_size;

  // self coupling  
  if ( reorder ) {
  START_LOCKED_MASTER(threading)
    int k=0, index, *it = in->index_table, *dt = in->table_dim;
    j=0;
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              for ( i=0; i<sc_size; i++, j++ )
                Aoo[j] = sc_in[ sc_size*index+i ];
            } else {
              for ( i=0; i<sc_size; i++, k++ )
                Aee[k] = sc_in[ sc_size*index+i ];
            }
          }
  END_LOCKED_MASTER(threading)
  } else {
    j = op->num_even_sites*sc_size;
    int start=0;
    int end=0;
    compute_core_start_end_custom(0, j, &start, &end, l, threading, 1);

    for ( i=start; i<end; i++ )
      Aee[i] = sc_in[i]; // even sites
    
    compute_core_start_end_custom(j, n*sc_size, &start, &end, l, threading, 1);
    for ( i=start; i<end; i++ )
      Aee[i] = sc_in[i]; // even sites

    SYNC_CORES(threading);
  }
  
  // neighbor couplings
  if ( reorder ) {
  START_LOCKED_MASTER(threading)

    int k=0, index, *it = in->index_table, *dt = in->table_dim, site_size=4*nc_size;
    config_float oAe=op->D, eAo=(op->D)+site_size*op->num_even_sites;
    j=0;
    for ( t=0; t<le[T]; t++ )
      for ( z=0; z<le[Z]; z++ )
        for ( y=0; y<le[Y]; y++ )
          for ( x=0; x<le[X]; x++ ) {
            index = site_index( t, z, y, x, dt, it );
            if ( (t+z+y+x+oe_offset)%2 == 1 ) {
              for ( i=0; i<site_size; i++, j++ ) {
                eAo[j] = nc_in[ site_size*index+i ];
              }
            } else {
              for ( i=0; i<site_size; i++, k++ ) {
                oAe[k] = nc_in[ site_size*index+i ];
              }
            }
          }
  END_LOCKED_MASTER(threading)          
  } else {

    j = n*4*nc_size;
    int start=0;
    int end=0;
    compute_core_start_end_custom(0, j, &start, &end, l, threading, 1);
    for ( i=start; i<end; i++ )
      op->D[i] = nc_in[i];

    SYNC_CORES(threading);

  }


  // prepare vectorized data
  int start=0;
  int end=0;
  compute_core_start_end_custom(0, n, &start, &end, l, threading, 1);
  int n_per_core = end-start;

  int n_site_var=l->num_lattice_site_var;
  int nvecs=n_site_var/2;
  int offset_v = n_site_var*n_site_var;

 
  // self coupling: store
  copy_coarse_operator_clover_to_vectorized_layout_float(
      op->clover + start*sc_size,
      op->clover_vectorized + start*offset_v,
      n_per_core, nvecs);
  SYNC_CORES(threading);

  // neighborhood coupling: store
  copy_coarse_operator_to_vectorized_layout_float(
      op->D + 4*start*nc_size,              // D is complex_float*
      op->D_vectorized + 2*4*start*nc_size, // 2 is for complex: D_vectorized is float* 
      n_per_core, nvecs);

  // self couplint:  make inverse for odd sites
  float *work=NULL;
  int work_size=2*(2*n_site_var*n_site_var + n_site_var);
  MALLOC( work, float, work_size );
  compute_core_start_end_custom(op->num_even_sites, n, &start, &end, l, threading, 1);
  for(int a=start; a<end; a++) {
    K_coarse_sc_inverse(
        op->clover_vectorized + a*offset_v,
	op->clover_vectorized + a*offset_v, nvecs, work);
  }
  FREE( work, float, work_size );  

  SYNC_CORES(threading);

}
#endif


#ifdef VECTORIZE_COARSE_OPERATOR_float

void coarse_n_hopping_term_float_vectorized( vector_float out, vector_float in, 
					       operator_float_struct *op, 
					       const int amount,
					       level_struct *l,
					       struct Thread *threading ) {
  //
  int mu, i, index;
  int num_site_var=l->num_lattice_site_var;
  int num_4link_var=4*l->num_lattice_site_var*l->num_lattice_site_var;
  int num_link_var=l->num_lattice_site_var*l->num_lattice_site_var;
  int start=0;
  int num_lattice_sites=l->num_inner_lattice_sites;
  int plus_dir_param=_FULL_SYSTEM;
  int minus_dir_param=_FULL_SYSTEM;
  const complex_float* in_pt;
  complex_float *out_pt;
  const float *vD_pt;
  int core_start;
  int core_end;


  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_float( out, 0, l, threading );

  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }

  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_float( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites;
    num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0;
    num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);


  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 0*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    K_coarse_n_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }

  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 1*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    K_coarse_n_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 2*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    K_coarse_n_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  SYNC_CORES(threading)
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 3*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    K_coarse_n_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_float( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_float( in, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];

    K_coarse_n_hopp_float( out_pt, in_pt, vD_pt, l );
    
    vD_pt += 2*num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    K_coarse_n_hopp_float( out_pt, in_pt, vD_pt, l );
    
    vD_pt += 2*num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    K_coarse_n_hopp_float( out_pt, in_pt, vD_pt, l );
    
    vD_pt += 2*num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    K_coarse_n_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_float( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)


}
#endif


#ifdef VECTORIZE_COARSE_OPERATOR_float

void coarse_hopping_term_float_vectorized( vector_float out, vector_float in, 
					   operator_float_struct *op, 
					   const int amount,
					   level_struct *l,
					   struct Thread *threading ) {
  //
  START_NO_HYPERTHREADS(threading);
  
  int mu, i, index;
  int num_site_var=l->num_lattice_site_var;
  int num_4link_var=4*l->num_lattice_site_var*l->num_lattice_site_var;
  int num_link_var=l->num_lattice_site_var*l->num_lattice_site_var;
  int start=0;
  int num_lattice_sites=l->num_inner_lattice_sites;
  int plus_dir_param=_FULL_SYSTEM, minus_dir_param=_FULL_SYSTEM;
  const complex_float* in_pt;
  complex_float* out_pt;
  float* vD_pt;

  int core_start;
  int core_end;
  
  // assumptions (1) self coupling has already been performed
  //          OR (2) "out" is initialized with zeros
  set_boundary_float( out, 0, l, threading );
  
  if ( amount == _EVEN_SITES ) {
    minus_dir_param = _ODD_SITES;
    plus_dir_param = _EVEN_SITES;
  } else if ( amount == _ODD_SITES ) {
    minus_dir_param = _EVEN_SITES;
    plus_dir_param = _ODD_SITES;
  }
  
  START_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in -mu direction
      ghost_sendrecv_float( in, mu, -1, &(op->c), minus_dir_param, l );
    }
  }
  END_MASTER(threading)
  SYNC_CORES(threading)
  
  if ( amount == _EVEN_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  } else if ( amount == _ODD_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
  
  // compute U_mu^dagger coupling
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 0*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+T];
    K_coarse_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  SYNC_CORES(threading);
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 1*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Z];
    K_coarse_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  SYNC_CORES(threading);
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 2*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+Y];
    K_coarse_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  SYNC_CORES(threading);
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    in_pt = in + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*(num_4link_var*op->neighbor_table[index] + 3*num_link_var);
    index++;
    out_pt = out + num_site_var*op->neighbor_table[index+X];
    K_coarse_daggered_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // communicate in +mu direction
      ghost_sendrecv_float( out, mu, +1, &(op->c), plus_dir_param, l );    
    }
    for ( mu=0; mu<4; mu++ ) {
      // wait for -mu direction
      ghost_wait_float( in, mu, -1, &(op->c), minus_dir_param, l );    
    }
  }
  END_LOCKED_MASTER(threading)
      
  if ( amount == _EVEN_SITES ) {
    start = 0; num_lattice_sites = op->num_even_sites;
  } else if ( amount == _ODD_SITES ) {
    start = op->num_even_sites, num_lattice_sites = op->num_odd_sites;
  }
  compute_core_start_end_custom(start, start+num_lattice_sites, &core_start, &core_end, l, threading, 1);
 
  // compute U_mu couplings
  for ( i=core_start; i<core_end; i++ ) {
    index = 5*i;
    out_pt = out + num_site_var*op->neighbor_table[index];
    vD_pt = op->D_vectorized + 2*num_4link_var*op->neighbor_table[index];
    index++;
    in_pt = in + num_site_var*op->neighbor_table[index+T];
    K_coarse_hopp_float( out_pt, in_pt, vD_pt, l );
  
    vD_pt += 2*num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Z];
    K_coarse_hopp_float( out_pt, in_pt, vD_pt, l );
    
    vD_pt += 2*num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+Y];
    K_coarse_hopp_float( out_pt, in_pt, vD_pt, l );
    
    vD_pt += 2*num_link_var;
    in_pt = in + num_site_var*op->neighbor_table[index+X];
    K_coarse_hopp_float( out_pt, in_pt, vD_pt, l );
  }
  
  START_LOCKED_MASTER(threading)
  if ( op->c.comm ) {
    for ( mu=0; mu<4; mu++ ) {
      // wait for +mu direction
      ghost_wait_float( out, mu, +1, &(op->c), plus_dir_param, l );
    }
  }
  END_LOCKED_MASTER(threading)

  END_NO_HYPERTHREADS(threading)

}
#endif

void K_coarse_aggregate_self_couplings_float( vector_float eta1, vector_float eta2, vector_float phi,
						  schwarz_float_struct *s, level_struct *l,
						  struct Thread *threading) {
  
  int i, mu, index1, index2, length;
  int *index_dir;
  int *neighbor = s->op.neighbor_table;
  int n = l->num_lattice_site_var;
  int Dls = n*n;
  int Dss = 4*n*n;
  vector_float eta1_pt, eta2_pt, phi_pt;
  config_float D_pt;
  config_float D = s->op.D;

  START_LOCKED_MASTER(threading)

  vector_float_define( eta1, 0, 0, l->vector_size, l );
  vector_float_define( eta2, 0, 0, l->vector_size, l );  
  coarse_spinwise_self_couplings_float( eta1, eta2, phi, s->op.clover, l->inner_vector_size, l );
  
  for ( mu=0; mu<4; mu++ ) { // direction mu
    length = l->is_float.agg_length[mu];
    index_dir = l->is_float.agg_index[mu];
    
    for ( i=0; i<length; i++ ) {
      index1 = index_dir[i];
      index2 = neighbor[5*index1+mu+1];
      D_pt = D + Dss*index1 + Dls*mu;
      phi_pt = phi + n*index2;
      eta1_pt = eta1 + n*index1;
      eta2_pt = eta2 + n*index1;
      coarse_spinwise_n_hopp_float( eta1_pt, eta2_pt, phi_pt, D_pt, l );
      phi_pt = phi + n*index1;
      eta1_pt = eta1 + n*index2;
      eta2_pt = eta2 + n*index2;
      coarse_spinwise_n_daggered_hopp_float( eta1_pt, eta2_pt, phi_pt, D_pt, l );
    }
  }
  END_LOCKED_MASTER(threading)

}


void K_set_coarse_self_coupling_float( vector_float spin_0_1, vector_float spin_2_3, 
				       const vector_float *V, const int n, level_struct *l,
				       struct Thread *threading,
				       int start, int end) {
  
  int i, j, k, m, k1, k2;
  int num_aggregates = l->is_float.num_agg;
  int num_eig_vect = l->next_level->num_lattice_site_var/2;
  int aggregate_size = l->inner_vector_size / num_aggregates;
  int offset = l->num_lattice_site_var/2;
  int clover_site_size = (l->next_level->num_lattice_site_var*(l->next_level->num_lattice_site_var+1))/2;
  vector_float spin_0_1_pt, spin_2_3_pt, interpolation_data;
  config_float clover_pt;
  config_float clover = l->next_level->op_float.clover;  

  
  // U(x) = [ A B      , A=A*, D=D*, C = -B*
  //          C D ]
  // storage order: upper triangle of A, upper triangle of D, B, columnwise
  // diagonal coupling
  
  
  for ( j=start; j<end; j++ ) {
    spin_0_1_pt = spin_0_1 + j*aggregate_size;
    spin_2_3_pt = spin_2_3 + j*aggregate_size;
    clover_pt = clover + j*clover_site_size;
      
    for ( k=0; k<=n; k++ ) {
      k1 = (n*(n+1))/2+k; k2 = (n*(n+1))/2+k+(num_eig_vect*(num_eig_vect+1))/2;

      interpolation_data = V[k] + j*aggregate_size;
      
      for ( i=0; i<aggregate_size; ) {
        // A
        for ( m=0; m<offset; m++, i++ ){
          clover_pt[ k1 ] += conj_float( interpolation_data[i] ) * spin_0_1_pt[i];
	}
        // D
        for ( m=0; m<offset; m++, i++ ){
          clover_pt[ k2 ] += conj_float( interpolation_data[i] ) * spin_2_3_pt[i];
	}
      }
    }

    for ( k=0; k<num_eig_vect; k++ ) {
      k1 = num_eig_vect*(num_eig_vect+1+n) + k;
      interpolation_data = V[k] + j*aggregate_size;
        
      for ( i=0; i<aggregate_size; ) {
        // B
        for ( m=0; m<offset; m++, i++ ){
          clover_pt[ k1 ] += conj_float( interpolation_data[i] ) * spin_2_3_pt[i];
	}
        i += offset;
      }
    }
  }
}


void K_coarse_aggregate_neighbor_couplings_float( vector_float eta1, vector_float eta2, vector_float phi,
						  const int mu, schwarz_float_struct *s, level_struct *l,
						  struct Thread *threading){
  
  int i, index1, index2;
  int length = l->is_float.agg_boundary_length[mu];
  int *index_dir = l->is_float.agg_boundary_index[mu];
  int *neighbor = l->is_float.agg_boundary_neighbor[mu];
  int n = l->num_lattice_site_var;
  int Dls = n*n;
  int Dss = 4*n*n;
  vector_float eta1_pt, eta2_pt, phi_pt;
  config_float D_pt;
  config_float D = s->op.D;
  
  START_LOCKED_MASTER(threading)
  vector_float_define( eta1, 0, 0, l->vector_size, l );
  vector_float_define( eta2, 0, 0, l->vector_size, l ); 

  // requires the positive boundaries of phi to be communicated befor
  for ( i=0; i<length; i++ ) {
    index1 = index_dir[i];
    index2 = neighbor[i];
    D_pt = D + Dss*index1 + Dls*mu;
    phi_pt = phi + n*index2;
    eta1_pt = eta1 + n*index1;
    eta2_pt = eta2 + n*index1;
    coarse_spinwise_hopp_float( eta1_pt, eta2_pt, phi_pt, D_pt, l );
  }
  END_LOCKED_MASTER(threading)

}

void K_set_coarse_neighbor_coupling_float( vector_float spin_0_1, vector_float spin_2_3, 
					   vector_float *V, const int mu, const int n, level_struct *l, struct Thread *threading, int start, int end ) {
  
  int i, i1, j, k, m;
  int num_aggregates = l->is_float.num_agg;
  int num_eig_vect = l->next_level->num_lattice_site_var/2;
  int offset = l->num_lattice_site_var/2;
  int nlsv = l->num_lattice_site_var;
  int D_link_size = num_eig_vect*num_eig_vect*4;
  int  *index_dir = l->is_float.agg_boundary_index[mu];
  int aggregate_boundary_sites = l->is_float.agg_boundary_length[mu]/num_aggregates;
      
  vector_float spin_0_1_pt, spin_2_3_pt;
  vector_float interpolation_data;
  config_float D_pt;
  config_float D = l->next_level->op_float.D;
  
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // storage order: A, C, B, D, each column wise
  //  for ( j=0; j<num_aggregates; j++ ) {
  for ( j=start; j<end; j++ ) {
    D_pt = D+(j*4+mu)*D_link_size;
      
    for ( k=0; k<num_eig_vect; k++ ) {
      int k1 = n*num_eig_vect + k;
      int k2 = (n+num_eig_vect)*num_eig_vect + k;
      int k3 = (n+2*num_eig_vect)*num_eig_vect + k;
      int k4 = (n+3*num_eig_vect)*num_eig_vect + k;

      for ( i=0; i<aggregate_boundary_sites; i++ ) {
	i1=j*aggregate_boundary_sites+i;
        spin_0_1_pt = spin_0_1 + nlsv*index_dir[i1];
        interpolation_data = V[k] + nlsv*index_dir[i1];
        // A
        for ( m=0; m<offset; m++ ){
          D_pt[ k1 ] += conj_float( interpolation_data[m] ) * spin_0_1_pt[m];
	}
        // C
        for ( ; m<2*offset; m++ ){
          D_pt[ k2 ] += conj_float( interpolation_data[m] ) * spin_0_1_pt[m];
	}

        spin_2_3_pt = spin_2_3 + nlsv*index_dir[i1];
        // B
        for ( m=0; m<offset; m++ )
          D_pt[ k3 ] += conj_float( interpolation_data[m] ) * spin_2_3_pt[m];
        // D
        for ( ; m<2*offset; m++ )
          D_pt[ k4 ] += conj_float( interpolation_data[m] ) * spin_2_3_pt[m];
      }
    }
  }
}


#ifdef SETUP_OPTIMIZED_float
void coarse_operator_float_setup_vectorized(vector_float *V, level_struct *l, struct Thread *threading ) {
  
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

  double t0, t1;
#pragma omp master
  {
    t0 = MPI_Wtime();
  }

  vector_float buffer1 = l->vbuf_float[4], buffer2 = l->vbuf_float[5];
  
  int mu;
  int n = l->num_eig_vect;
  int i, j;
  int num_aggregates = l->is_float.num_agg;
  int clover_site_size = (l->next_level->num_lattice_site_var*(l->next_level->num_lattice_site_var+1))/2;
  int D_link_size = 4*l->num_eig_vect*l->num_eig_vect*4;  // size of links in all 4 directions

  START_LOCKED_MASTER(threading)
  operator_float_define( &(l->next_level->op_float), l->next_level );
  END_LOCKED_MASTER(threading)
  SYNC_HYPERTHREADS(threading)

  int start, end;
  compute_core_start_end(0, num_aggregates, &start, &end, l, threading);
  
  // each thread loops overs its aggregates and then over internal d.o.f.
  for(int a=start; a<end; a++){
    for ( j=0; j<D_link_size; j++ ){
      l->next_level->op_float.D[a*D_link_size + j] = _COMPLEX_float_ZERO;
    }
    for ( j=0; j<clover_site_size; j++ ){
      l->next_level->op_float.clover[a*clover_site_size +j] = _COMPLEX_float_ZERO;
    }
  }

  // Dirc operators
  //  depth==0: fine dirac operator
  //  depth> 0: coarse dirac operator

  // self coupling
  void (*aggregate_self_coupling)() = (l->depth==0)?K_d_plus_clover_aggregate_float:K_coarse_aggregate_self_couplings_float;

  // hopping
  void (*aggregate_neighbor_coupling)() = (l->depth==0)?K_d_neighbor_aggregate_float:K_coarse_aggregate_neighbor_couplings_float;

  SYNC_CORES(threading);
  // for all test vectors V[i]:
  for ( i=0; i<n; i++ ) {
    START_LOCKED_MASTER(threading)
    for ( mu=0; mu<4; mu++ ) {
      // update ghost cells of V[i]
      negative_sendrecv_float( V[i], mu, &(l->s_float.op.c), l );
    }
    END_LOCKED_MASTER(threading)

    // apply self coupling of block-and-2spin-restricted dirac operator for each aggregate
    aggregate_self_coupling( buffer1, buffer2, V[i], &(l->s_float), l, threading);
    SYNC_CORES(threading);

    K_set_coarse_self_coupling_float( buffer1, buffer2, V, i, l, threading, start, end);
    SYNC_CORES(threading);

    for ( mu=0; mu<4; mu++ ) {
      START_LOCKED_MASTER(threading)
      // finish updating ghostcells of V[i]
      negative_wait_float( mu, &(l->s_float.op.c), l );      
      // apply 2spin-restricted dirac operator for direction mu for all aggregates
      END_LOCKED_MASTER(threading)
      aggregate_neighbor_coupling( buffer1, buffer2, V[i], mu, &(l->s_float), l, threading );      
      SYNC_CORES(threading);
      K_set_coarse_neighbor_coupling_float( buffer1, buffer2, V, mu, i, l, threading, start, end );
      SYNC_CORES(threading);
    }
  }


#pragma omp master
  {
  t1 = MPI_Wtime();
  if ( g.print > 0 ) printf0("depth: %d, time spent for setting up next coarser operator: %lf seconds\n", l->depth, t1-t0 );
  }

}
#endif

#endif
