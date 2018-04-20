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

#include "K_vector.h"
#include "K_dirac.h"
#include "dirac_float.h" // to be removed

void K_set_clover_double( double *out, complex_double *in ) { 

}

void K_set_clover_float( float *out, complex_double *in ) {
  // in:   clover on a site wtth a default packing
  // out:  clover on a site with an optimzied packing
  //
  // format for out
  //   0--35: upper 6x6 block, Hermite Matrix 
  //          diangonal and upper right half
  //     0-- 5: diagnal, real
  //     6--15:  H_{12} -- H_{16}
  //    16--23:  H_{23} -- H_{26}
  //    24--29:  H_{34} -- H_{36}
  //    30--33:  H_{45} -- H_{46}
  //    34--35:  H_{56}
  //
  //   36--71: lower 6x6 block, Hermite Matrix 
  //          diangonal and upper right half
  //     (upper block + 36)


  {// upper block
    int offset_out=0;
    int offset_in=0;
#pragma loop noalias
    for(int i=0; i<6; i++){
      out[offset_out+i]=(float)creal(in[offset_in+i]);
    }
    offset_out=6;
    offset_in=12;
    const double *in_real=(double*)(in+offset_in);
#pragma loop noalias
    for(int i=0; i<30; i++){
      out[offset_out+i]=(float)in_real[i];
    }
  }
  {  // lower block
    int offset_out=36;
    int offset_in=6;
#pragma loop noalias
    for(int i=0; i<6; i++){
      out[offset_out+i]=(float)creal(in[offset_in+i]);
    }
    offset_out=42;
    offset_in=27;
    const double *in_real=(double*)(in+offset_in);
#pragma loop noalias
    for(int i=0; i<30; i++){
      out[offset_out+i]=(float)in_real[i];
    }
  }
}


void K_clover_double( vector_double eta, vector_double phi, operator_double_struct *op, int start, int end,
                        level_struct *l, struct Thread *threading ) {
  
  clover_double( eta+start, phi+start, op->clover+((start/12)*42), end-start, l, threading );
  
}


void K_clover_float( vector_float eta, vector_float phi, operator_float_struct *op, int start, int end,
                       level_struct *l, struct Thread *threading ) {
  float *clov = op->clover_vectorized+start*72;
  for ( int i=start; i<end; i+=12 ) {
    K_site_clover_float( (float*)(eta+i), (float*)(phi+i), clov+6*i );
  }
}


void K_site_clover_double( double *eta, const double *phi, const double *clover ) { 

}

void K_site_clover_float_test( float *eta, const float *phi, float *clover, complex_float *clover_org ) {

  
  for(int block=0; block<2; block++){ 
    // block=0: upper 6x6
    // block=1: lower  6x6
    const complex_float *in=(complex_float*)phi + 6*block;
    complex_float *out=(complex_float*)eta + 6*block;
    const float *diag=clover+36*block;
    const float *clov=diag+6;
    complex_float c[6][6];

    c[0][0]=diag[0];
    c[0][1]=clov[0]+I*clov[1];
    c[0][2]=clov[2]+I*clov[3];
    c[0][3]=clov[4]+I*clov[5];
    c[0][4]=clov[6]+I*clov[7];
    c[0][5]=clov[8]+I*clov[9];

    c[1][0]=conj(c[0][1]);
    c[1][1]=diag[1];
    c[1][2]=clov[10]+I*clov[11];
    c[1][3]=clov[12]+I*clov[13];
    c[1][4]=clov[14]+I*clov[15];
    c[1][5]=clov[16]+I*clov[17];

    c[2][0]=conj(c[0][2]);
    c[2][1]=conj(c[1][2]);
    c[2][2]=diag[2];
    c[2][3]=clov[18]+I*clov[19];
    c[2][4]=clov[20]+I*clov[21];
    c[2][5]=clov[22]+I*clov[23];

    c[3][0]=conj(c[0][3]);
    c[3][1]=conj(c[1][3]);
    c[3][2]=conj(c[2][3]);
    c[3][3]=diag[3];
    c[3][4]=clov[24]+I*clov[25];
    c[3][5]=clov[26]+I*clov[27];

    c[4][0]=conj(c[0][4]);
    c[4][1]=conj(c[1][4]);
    c[4][2]=conj(c[2][4]);
    c[4][3]=conj(c[3][4]);
    c[4][4]=diag[4];
    c[4][5]=clov[28]+I*clov[29];

    c[5][0]=conj(c[0][5]);
    c[5][1]=conj(c[1][5]);
    c[5][2]=conj(c[2][5]);
    c[5][3]=conj(c[3][5]);
    c[5][4]=conj(c[4][5]);
    c[5][5]=diag[5];

    for(int i=0; i<6; i++){
      out[i]=0.0;
      for(int j=0; j<6; j++){
	out[i]+=c[i][j]*in[j];
      }
    }
  } // block

  /*
  complex_float eta_org[12];
  site_clover_float(eta_org, (complex_float*)phi, clover_org);


  double diff2=0;
  float *eta_test = (float*)eta_org;
  for(int i=0; i<24; i++){
    diff2+=(eta_test[i]-eta[i])*(eta_test[i]-eta[i]);
    printf("i=%d, %f %f\n",i, eta_test[i], eta[i]);
  }

  if(diff2>1e-6){
    printf("diff2=%f, %s\n",diff2, __func__);
    fflush(0);
    exit(1);
  }
  */
}

void K_site_clover_float( float *eta, const float *phi, float *clover ) {

  // upper 6x6
  K_site_clover_2spin_float(eta, phi, clover);

  // lower 6x6
  phi+=12;
  eta+=12;
  clover+=36;
  K_site_clover_2spin_float(eta, phi, clover);
  return;
}


void K_site_clover_invert_double( double *clover_in, double *clover_out ) { }

void K_site_clover_invert_float( float *clover_in, float *clover_out ) {


  // in:  compact format with 72 real numbers
  // out: compact format with 72 real numbers

  static const int N=6;
  complex_float  A[N*N];
  complex_float  A_inv[N*N];
  float *Atmp=(float*)A;


  for(int block=0; block<2; block++){
    const float *AA=clover_in + 36*block;
    float *AA_inv=clover_out + 36*block;

    for(int i=0; i<N; i++){
      int ij=i+N*i;
      Atmp[2*ij  ] = AA[i];
      Atmp[2*ij+1] = 0.0;
      //      A[ij] = AA[i];
    }
    int ii=N;
    for(int i=0; i<N; i++){
      for(int j=i+1; j<N; j++){
	int ij=i+N*j;
	int ji=j+N*i;
	//	A[ij]=AA[ii]+I*AA[ii+1];
	//	A[ji]=AA[ii]-I*AA[ii+1];
	//	ii+=2;
	Atmp[2*ij]   = AA[ii];
	Atmp[2*ji]   = AA[ii];
	ii++;
	Atmp[2*ij+1] = AA[ii];
	Atmp[2*ji+1] = -AA[ii];
	ii++;
      }
    }


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
    
    complex_float b[N];
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

    for(int i=0; i<N; i++){
      AA_inv[i]=creal(A_inv[N*i+i]);
    }
    ii=N;
    for(int i=0; i<N; i++){
      for(int j=i+1; j<N; j++){
	int ij= i+N*j;
	AA_inv[ii] = crealf(A_inv[ij]);
	ii++;
	AA_inv[ii] = cimagf(A_inv[ij]);
	ii++;
      }
    }

  } // block
}

void block_oddeven_coupling_double( double *eta, double *D, double *phi, 
				    int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor ) { }

void block_oddeven_coupling_float( float *eta,  float *D, float *phi, 
				    int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor ) {

  prepare_I; // for accum_times_I, accum_times_minusI

  {  // T direction
#define MU 0
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
#include "K_dirac_su3local_n.h"
#undef MU
  }

  {  //  Z direction
#define MU 1
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
#include "K_dirac_su3local_n.h"
#undef MU
  }

  {  //  Y direction
#define MU 2
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
    #include "K_dirac_su3local_n.h"
#undef MU
  }

  {  //  X direction
#define MU 3
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
#include "K_dirac_su3local_n.h"
#undef MU
  }
 }

void block_n_oddeven_coupling_double( double *eta, double *D, double *phi, 
				    int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor ) { }

void block_n_oddeven_coupling_float( float *eta,  float *D, float *phi, 
				    int *start1, int *start2, int *end1, int *end2, int **index, int *neighbor ) {

  prepare_I; // for accum_times_I, accum_times_minusI
#define ACCUM_PLUS
  {  // T direction
#define MU 0
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
#include "K_dirac_su3local_n.h"
#undef MU
  }

  {  //  Z direction
#define MU 1
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
#include "K_dirac_su3local_n.h"
#undef MU
  }

  {  //  Y direction
#define MU 2
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
    #include "K_dirac_su3local_n.h"
#undef MU
  }

  {  //  X direction
#define MU 3
    int *ind=index[MU];
    int start=start1[MU];
    int end=end1[MU];
#include "K_dirac_su3local_p.h"
    start=start2[MU];
    end=end2[MU];
#include "K_dirac_su3local_n.h"
#undef MU
  }
#undef ACCUM_PLUS
 }


void K_block_d_double( vector_double *eta, vector_double *phi, int start, schwarz_double_struct *s){};

void K_block_d_float( complex_float *ceta, complex_float *cphi, int start, schwarz_float_struct *s){

  int *length = s->dir_length;
  int **index = s->index;
  int *neighbor = s->op.neighbor_table;

  float *eta = (float*)ceta;
  const float *phi = (const float*)cphi;
  float* D = (float*)(s->op.D + (start/12)*36);

  prepare_I; // for accum_times_I, accum_times_minusI

  {  // T direction
#define MU 0
    int *ind=index[MU];
    int start=0;
    int end=length[MU];
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef MU
   }
  {  // Z direction
#define MU 1
    int *ind=index[MU];
    int start=0;
    int end=length[MU];
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef MU
   }
  {  // Y direction
#define MU 2
    int *ind=index[MU];
    int start=0;
    int end=length[MU];
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef MU
   }
  {  // X direction
#define MU 3
    int *ind=index[MU];
    int start=0;
    int end=length[MU];
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef MU
   }
};


void K_d_plus_clover_aggregate_float( vector_float eta1, vector_float eta2, vector_float phi0, schwarz_float_struct *s, level_struct *l, struct Thread *threading ) {

  int i, length;
  int *ind;
  int *neighbor = s->op.neighbor_table;
  float *D = (float*)(s->op.D);

  prepare_I; // for accum_times_I, accum_times_minusI

  int start=0;
  int end=0;
  compute_core_start_end_custom(0, l->num_inner_lattice_sites, &start, &end, l, threading, 1);

  // add clover term/shift
  if ( g.csw == 0.0 ) {
    //    spinwise_float_skalarmultiply( eta1, eta2, phi, s->op.shift, 0, l->inner_vector_size, l );
    error0("csw==0 is not supported in %s\n",__func__);
  } else {
    float *pt_clover = s->op.clover_vectorized + 72*start;
    float *pt_eta1 = (float*)(eta1+12*start);
    float *pt_eta2 = (float*)(eta2+12*start);
    float *pt_phi= (float*)(phi0+12*start);
    for(i=start; i<end; i++) {
      K_site_clover_2spin_float(pt_eta1, pt_phi, pt_clover);
      for(int ii=0; ii<12; ii++){
	pt_eta2[ii]=0.0;
      }
      pt_eta1+=12;
      pt_eta2+=12;
      pt_phi+=12;
      pt_clover+=36;

      K_site_clover_2spin_float(pt_eta2, pt_phi, pt_clover);
      for(int ii=0; ii<12; ii++){
	pt_eta1[ii]=0.0;
      }
      pt_eta1+=12;
      pt_eta2+=12;
      pt_phi+=12;
      pt_clover+=36;
    }
  }
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

  // T dir
  {
#define MU 0
    length = l->is_float.agg_length[MU];
    compute_core_start_end(0, length, &start, &end, l, threading);
    ind = l->is_float.agg_index[MU];
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN34
#undef MU
  }
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)
  // Z dir
  {
#define MU 1
    length = l->is_float.agg_length[MU];
    compute_core_start_end(0, length, &start, &end, l, threading);
    ind = l->is_float.agg_index[MU];
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN34
#undef MU
  }
  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

  // Y dir
  {
#define MU 2
    length = l->is_float.agg_length[MU];
    compute_core_start_end(0, length, &start, &end, l, threading);
    ind = l->is_float.agg_index[MU];
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN34
#undef MU
  }

  // X dir
  {
#define MU 3
    length = l->is_float.agg_length[MU];
    compute_core_start_end(0, length, &start, &end, l, threading);
    ind = l->is_float.agg_index[MU];
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_p.h"
#include "K_dirac_su3local_n.h"
#undef SPIN34
#undef MU
  }

  SYNC_HYPERTHREADS(threading)
  SYNC_CORES(threading)

}

void K_d_neighbor_aggregate_float( vector_float eta1, vector_float eta2, vector_float phi0, const int mu, schwarz_float_struct *s, level_struct *l, struct Thread *threading ) {
  
  float *D = (float*)(s->op.D);
  int start=0;
  int end=0;
  int length = l->is_float.agg_boundary_length[mu];
  int *ind = l->is_float.agg_boundary_index[mu];
  int *neighbor = l->is_float.agg_boundary_neighbor[mu];

  compute_core_start_end(0, length, &start, &end, l, threading);

  prepare_I;

  // requires the positive boundaries of phi to be communicated befor
  if ( mu == T ) {
    // T dir
#define MU 0
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_boundary.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_boundary.h"
#undef SPIN34
#undef MU

  } else if ( mu == Z ) {
    // Z dir
#define MU 1
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_boundary.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_boundary.h"
#undef SPIN34
#undef MU
  } else if ( mu == Y ) {
    // Y dir
#define MU 2
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_boundary.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_boundary.h"
#undef SPIN34
#undef MU
  } else if ( mu == X ) {
    // X dir
#define MU 3
    float *eta=(float*)eta1;
    float *phi=(float*)phi0;
#define SPIN12
#include "K_dirac_su3local_boundary.h"
#undef SPIN12
    eta=(float*)eta2;
    phi=(float*)phi0;
#define SPIN34
#include "K_dirac_su3local_boundary.h"
#undef SPIN34
#undef MU
  }
}



#endif // K_OPT
