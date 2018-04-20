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


#include "K_vector.h"

void copy_coarse_operator_clover_to_vectorized_layout_double(const config_double clover,
							     OPERATOR_TYPE_double *clover_vectorized,
							     int num_aggregates, int num_eig_vect);

void copy_coarse_operator_clover_to_vectorized_layout_float(const complex_float *clover,
							    float *clover_vectorized,
							    int num_aggregates, int num_eig_vect);

void K_coarse_sc_inverse(const float *clover_in, float *clover_out, int n, float *Atmp);


void coarse_self_couplings_double_vectorized( vector_double eta, vector_double phi, 
					      OPERATOR_TYPE_double *clover, int start, int end, level_struct *l );

void coarse_operator_float_setup_vectorized(vector_float *V, level_struct *l, struct Thread *threading );

#ifdef VECTORIZE_COARSE_OPERATOR_double
void coarse_oddeven_setup_double_set_couplings( operator_PRECISION_struct *in, int reorder, level_struct *l, struct Thread *threading );
#endif


#ifdef VECTORIZE_COARSE_OPERATOR_float
void coarse_oddeven_setup_float_set_couplings( operator_float_struct *in, 
					       int reorder, level_struct *l, 
					       struct Thread *threading );
#endif



static inline void block4x4_to_full_Hermite(float *work, const float *in, int block_num){

  static const int block_size=4;
  int N=block_size*block_num;

  int k=0;  // index to the input

  // diagonal block matrices
  for(int ib=0; ib<block_num; ib++){
    // diangonal real
#pragma loop noalias
    for(int i=0; i<block_size; i++){
      int ii=ib*block_size+i;
      int ij=ii*N+ii;
      work[2*ij]  =in[k];
      work[2*ij+1]=0.0;   // real
      k++;
    }
    // offdiangoal
    for(int i=0; i<block_size; i++){
#pragma loop noalias
      for(int j=i+1; j<block_size; j++){
	int ii=ib*block_size+i;
	int jj=ib*block_size+j;
	int ij=ii*N + jj;
	int ji=jj*N + ii;
	work[2*ij]=in[k];
	work[2*ji]=in[k];
	k++;
	work[2*ij+1]=in[k];
	work[2*ji+1]=-in[k];  // Hermite
	k++;
      }
    }
  } // ib: end of diagonal blocks


  for(int ib=0; ib<block_num; ib++){
    for(int jb=ib+1; jb<block_num; jb++){
      for(int i=0; i<block_size; i++){
#pragma loop noalias
	for(int j=0; j<block_size; j++){
	  int ii=ib*block_size+i;
	  int jj=jb*block_size+j;
	  int ij=ii*N + jj;
	  int ji=jj*N + ii;
	  work[2*ij]=in[k];
	  work[2*ji]=in[k];
	  k++;
	  work[2*ij+1]=in[k];
	  work[2*ji+1]=-in[k]; // Hermite
	  k++;
	} 
      } // i,j
    }}// ib,jb
}


static inline void full_to_block4x4_Hermite(float *out, const float *work, int block_num){

  static const int block_size=4;
  int N=block_size*block_num;

  int k=0;  // index to the output

  // diagonal block matrices
  for(int ib=0; ib<block_num; ib++){
    // diangonal real
#pragma loop noalias
    for(int i=0; i<block_size; i++){
      int ii=ib*block_size+i;
      int ij=ii*N+ii;
      out[k]=work[2*ij];
      k++;
    }
    // offdiangoal
    for(int i=0; i<block_size; i++){
#pragma loop noalias
      for(int j=i+1; j<block_size; j++){
	int ii=ib*block_size+i;
	int jj=ib*block_size+j;
	int ij=ii*N + jj;
	out[k]=work[2*ij];
	k++;
	out[k]=work[2*ij+1];
	k++;
      }
    }
  } // ib: end of diagonal blocks

  for(int ib=0; ib<block_num; ib++){
    for(int jb=ib+1; jb<block_num; jb++){
      for(int i=0; i<block_size; i++){
#pragma loop noalias
	for(int j=0; j<block_size; j++){
	  int ii=ib*block_size+i;
	  int jj=jb*block_size+j;
	  int ij=ii*N+jj;
	  out[k]=work[2*ij];
	  k++;
	  out[k]=work[2*ij+1];
	  k++;
	} 
      } // i,j
    }}// ib,jb

}



static inline void K_sc_site_block4x4_reference(float *out, float *clover, const float *in, const int N){

  static const int block_size=4;
  int block_num=N/block_size;
  const int size_diag   =  block_size*block_size;
  const int size_offdiag=2*block_size*block_size;

  const complex_float *cin = (const complex_float*)in;
  complex_float *cout = (complex_float*)out;

  // diagonal blocks
  const float *A=clover;
  for(int ib=0; ib<block_num; ib++){
    // diangonal real
#pragma loop noalias
    for(int i=0; i<block_size; i++){
      int ii=ib*block_size+i;
      float Aii=A[i];
      cout[ii] =Aii*cin[ii];
    }
    // offdiangonal
    const complex_float *AA=(const complex_float*)(A+block_size);
    int k=0;
    for(int i=0; i<block_size; i++){
#pragma loop noalias
      for(int j=i+1; j<block_size; j++){
	int ii=ib*block_size+i;
	int jj=ib*block_size+j;
	complex_float Aij=AA[k];
	cout[ii] += Aij*cin[ii];
	cout[jj] += conj(Aij)*cin[jj];
	k++;
      }
    }
    A+=size_diag;
  }

  const complex_float *AA=(const complex_float*)A;
  // off-diagonal blocks
  for(int ib=0; ib<block_num; ib++){
    for(int jb=ib+1; jb<block_num; jb++){
      int k=0;
      for(int i=0; i<block_size; i++){
#pragma loop noalias
	for(int j=0; j<block_size; j++){
	  int ii=ib*block_size+i;
	  int jj=jb*block_size+j;

	  complex_float Aij=AA[k];
	  cout[ii] += Aij*cin[jj];
	  cout[jj] += conj(Aij)*cin[ii];
	  k++;

	} 
      } // i,j
      AA+=size_offdiag/2;
    }}// ib,jb

  // multiply gamma 5
  for(int i=N/2; i<N; i++){
    cout[i]=-cout[i];
  }

}



static inline void K_sc_site_block4x4(float *out, float *clover, const float *in, const int N){

  static const int block_size=4;
  int block_num=N/block_size;
  const int size_diag   =  block_size*block_size;
  const int size_offdiag=2*block_size*block_size;

  // diagonal blocks: Hermite
  const float *A=clover;
  for(int ib=0; ib<block_num; ib++){
    register _fjsp_v2r8 in0,in1,in2,in3;
    register _fjsp_v2r8 out0,out1,out2,out3;
    register _fjsp_v2r8 A0,A1;
    register _fjsp_v2r8 A01,A02,A03,A12,A13,A23;
    register _fjsp_v2r8 diag0,diag1,diag2,diag3;

    int offset=2*ib*block_size; // 2 is for complex

    // load input vector
    register block4vec *in_ptr = (block4vec*)(in+offset);
    load_block4vec(in, in_ptr);

    // load matrix
    register block4diag *A_ptr = (block4diag*)(A);
    load_block4diag(A, A_ptr);

    // extract diagnal real
    diag0 =_fjsp_unpacklo_v2r8(A0,A0);
    diag1 =_fjsp_unpackhi_v2r8(A0,A0);
    diag2 =_fjsp_unpacklo_v2r8(A1,A1);
    diag3 =_fjsp_unpackhi_v2r8(A1,A1);

    // 1st line
    rmult(out0,diag0,in0);
    accum_cmult(out0, A01, in1);
    accum_cmult(out0, A02, in2);
    accum_cmult(out0, A03, in3);

    // 2nd line
    rmult(out1,diag1,in1);
    accum_conj_cmult(out1,A01,in0);
    accum_cmult(out1, A12, in2);
    accum_cmult(out1, A13, in3);

    // 3rd line
    rmult(out2,diag2,in2);
    accum_conj_cmult(out2,A02,in0);
    accum_conj_cmult(out2,A12,in1);
    accum_cmult(out2, A23, in3);

    // 4th line
    rmult(out3,diag3,in3);
    accum_conj_cmult(out3,A03,in0);
    accum_conj_cmult(out3,A13,in1);
    accum_conj_cmult(out3,A23,in2);

    // store output vector
    register block4vec* out_ptr = (block4vec*)(out+offset);
    store_block4vec(out_ptr, out);

    A+=size_diag;
  }


  // off-diagonal blocks
  for(int ib=0; ib<block_num; ib++){
    int ii=2*ib*block_size; // 2 is for complex

    // input vector "in_i"
    register _fjsp_v2r8 iin0,iin1,iin2,iin3;
    register block4vec *iin_ptr = (block4vec*)(in+ii);
    load_block4vec(iin, iin_ptr);

    // output vector "out_i"
    register _fjsp_v2r8 iout0,iout1,iout2,iout3;
    register block4vec *iout_ptr = (block4vec*)(out+ii);
    load_block4vec(iout, iout_ptr);


    for(int jb=ib+1; jb<block_num; jb++){

      int jj=2*jb*block_size; // 2 is for complex

      // input vector "in_j"
      register _fjsp_v2r8 jin0,jin1,jin2,jin3;
      register block4vec *jin_ptr = (block4vec*)(in+jj);
      load_block4vec(jin, jin_ptr);

      // output vector "out_j"
      register _fjsp_v2r8 jout0,jout1,jout2,jout3;
      register block4vec *jout_ptr = (block4vec*)(out+jj);
      load_block4vec(jout, jout_ptr);

      // matrix element Aij
      register _fjsp_v2r8 A00,A01,A02,A03;
      register _fjsp_v2r8 A10,A11,A12,A13;
      register _fjsp_v2r8 A20,A21,A22,A23;
      register _fjsp_v2r8 A30,A31,A32,A33;

      register block4mat *A_ptr = (block4mat*)(A);
      load_block4mat(A, A_ptr);

      // A00
      accum_cmult(     iout0, A00, jin0);
      accum_conj_cmult(jout0, A00, iin0);

      // A01
      accum_cmult(     iout0, A01, jin1);
      accum_conj_cmult(jout1, A01, iin0);

      // A02
      accum_cmult(     iout0, A02, jin2);
      accum_conj_cmult(jout2, A02, iin0);

      // A03
      accum_cmult(     iout0, A03, jin3);
      accum_conj_cmult(jout3, A03, iin0);

      // A10
      accum_cmult(     iout1, A10, jin0);
      accum_conj_cmult(jout0, A10, iin1);

      // A11
      accum_cmult(     iout1, A11, jin1);
      accum_conj_cmult(jout1, A11, iin1);

      // A12
      accum_cmult(     iout1, A12, jin2);
      accum_conj_cmult(jout2, A12, iin1);

      // A13
      accum_cmult(     iout1, A13, jin3);
      accum_conj_cmult(jout3, A13, iin1);

      // A20
      accum_cmult(     iout2, A20, jin0);
      accum_conj_cmult(jout0, A20, iin2);

      // A21
      accum_cmult(     iout2, A21, jin1);
      accum_conj_cmult(jout1, A21, iin2);

      // A22
      accum_cmult(     iout2, A22, jin2);
      accum_conj_cmult(jout2, A22, iin2);

      // A23
      accum_cmult(     iout2, A23, jin3);
      accum_conj_cmult(jout3, A23, iin2);

      // A30
      accum_cmult(     iout3, A30, jin0);
      accum_conj_cmult(jout0, A30, iin3);

      // A31
      accum_cmult(     iout3, A31, jin1);
      accum_conj_cmult(jout1, A31, iin3);

      // A32
      accum_cmult(     iout3, A32, jin2);
      accum_conj_cmult(jout2, A32, iin3);

      // A33
      accum_cmult(     iout3, A33, jin3);
      accum_conj_cmult(jout3, A33, iin3);

      store_block4vec(jout_ptr, jout);

      A+=size_offdiag;
    }
    store_block4vec(iout_ptr, iout);

  }// ib,jb

  // multiply gamma 5
  for(int i=N; i<2*N; i++){
    out[i]=-out[i];
  }

}


static inline void coarse_self_couplings_float_vectorized( vector_float eta, vector_float phi, 
                                                                 OPERATOR_TYPE_float *clover, int start, int end, level_struct *l ) {
#ifdef VECTORIZE_COARSE_OPERATOR_float
  int site_size = l->num_lattice_site_var;
  for(int i=start; i<end; i++) {
    //    K_sc_site_block4x4_reference((float *)(eta+i*site_size),
    //    		       clover+i*site_size*site_size, 
    //    		       (const float *)(phi+i*site_size), site_size);

    K_sc_site_block4x4((float *)(eta+i*site_size),
       		       clover+i*site_size*site_size, 
		       (const float *)(phi+i*site_size), site_size);


  }
#endif
}


static inline void K_coarse_n_daggered_hopp_float_reference( complex_float *eta, const complex_float *phi,
						   const float *D, level_struct *l ) {
    
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // note: minus sign of D = self_coupling - hopping_term is added here
  
  //  A = [ a00 a01 a02..]
  //      [ a10 a22...   ]
  //      [...           ]  with 2x2 a?? 
  // stored as
  // [ (a00 b00 c00 d00) (a01 b01 c01 d01)...] each in row order
  //   

  static const int blocksize=4;
  int bs2=blocksize/2;
  int N=l->num_lattice_site_var;
  int n=N/2;
  int numblock=N/blocksize;

  const complex_float *A=(const complex_float*)D;
  for(int ib=0; ib<numblock; ib++){
    const complex_float* in1_ii= phi    + ib*2;
    const complex_float* in2_ii=(phi+n) + ib*2;

    for(int jb=0; jb<numblock; jb++){

      complex_float* out1_jj= eta    + jb*2;
      complex_float* out2_jj=(eta+n) + jb*2;
      int offset=0;

      // A* : out_1j += (Aji)^* in_1i 
      for(int i=0; i<2; i++){
	for(int j=0; j<2; j++){
	  int ij=2*i+j;
	  out1_jj[j]+=conj(A[ij])*in1_ii[i];
	}}

      // -B* : out_2j += (-Aji)^* in_1i
      offset=bs2*bs2;
      for(int i=0; i<2; i++){
	for(int j=0; j<2; j++){
	  int ij=2*i+j + offset;
	  out2_jj[j]-=conj(A[ij])*in1_ii[i];
	}}

      // -C* : out_1j += (-Aji)^* in_2i
       offset=2*bs2*bs2;
      for(int i=0; i<2; i++){
	for(int j=0; j<2; j++){
	  int ij=2*i+j+offset;
	  out1_jj[j]-=conj(A[ij])*in2_ii[i];
	}}

      // D* : out_2j += (Aji)^* in_2i 
      offset=3*bs2*bs2;
      for(int i=0; i<2; i++){
	for(int j=0; j<2; j++){
	  int ij=bs2*i+j + offset;
	  out2_jj[j]+=conj(A[ij])*in2_ii[i];
	}}

      A+=blocksize*blocksize;
    }
  }
}




static inline void K_coarse_n_daggered_hopp_float( complex_float *eta, const complex_float *phi,
						   const float *D, level_struct *l ) {
    
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // note: minus sign of D = self_coupling - hopping_term is added here
  
  //  A = [ a00 a01 a02..]
  //      [ a10 a22...   ]
  //      [...           ]  with 2x2 a?? 
  // stored as
  // [ (a00 b00 c00 d00) (a01 b01 c01 d01)...] each in row order
  //   

  static const int blocksize=4;
  int N=l->num_lattice_site_var;
  int n=N/2;
  int numblock=N/blocksize;
  //  int bs2=blocksize/2;

  const complex_float *A=(const complex_float*)D;
  for(int ib=0; ib<numblock; ib++){
    // input vector
    //    const complex_float* in1_ii= phi    + ib*2;
    //    const complex_float* in2_ii=(phi+n) + ib*2;
    register _fjsp_v2r8 in1_ii0, in1_ii1;
    register _fjsp_v2r8 in2_ii0, in2_ii1;
    register block2vec *in1_ptr = (block2vec*)(phi    + ib*2);
    register block2vec *in2_ptr = (block2vec*)((phi+n)+ ib*2);
    load_block2vec(in1_ii, in1_ptr);
    load_block2vec(in2_ii, in2_ptr);

    for(int jb=0; jb<numblock; jb++){
      // output vector
      //      complex_float* out1_jj= eta    + jb*2;
      //      complex_float* out2_jj=(eta+n) + jb*2;
      register _fjsp_v2r8 out1_jj0, out1_jj1;
      register _fjsp_v2r8 out2_jj0, out2_jj1;
      register block2vec *out1_ptr = (block2vec*)(eta    + jb*2);
      register block2vec *out2_ptr = (block2vec*)((eta+n)+ jb*2);
      load_block2vec(out1_jj, out1_ptr);
      load_block2vec(out2_jj, out2_ptr);

      // matrix
      register _fjsp_v2r8 A00, A01, A10, A11;
      register _fjsp_v2r8 B00, B01, B10, B11;
      register _fjsp_v2r8 C00, C01, C10, C11;
      register _fjsp_v2r8 D00, D01, D10, D11;
      register block4mat *mat_ptr=(block4mat*)A;
      load_block4mat_ABCD(A,B,C,D, mat_ptr);
      
      // A* : out_1j += (Aji)^* in_1i 
      accum_conj_cmult(out1_jj0, A00, in1_ii0);
      accum_conj_cmult(out1_jj0, A10, in1_ii1);
      accum_conj_cmult(out1_jj1, A01, in1_ii0);
      accum_conj_cmult(out1_jj1, A11, in1_ii1);
      //int offset=0;
      //      for(int i=0; i<2; i++){
      //      	for(int j=0; j<2; j++){
      //      	  int ij=bs2*i+j;
      //      	  out1_jj[j]+=conj(A[ij])*in1_ii[i];
      //      	}}

      // -B* : out_2j += (-Aji)^* in_1i 
      accum_n_conj_cmult(out2_jj0, B00, in1_ii0);
      accum_n_conj_cmult(out2_jj0, B10, in1_ii1);
      accum_n_conj_cmult(out2_jj1, B01, in1_ii0);
      accum_n_conj_cmult(out2_jj1, B11, in1_ii1);
      //      offset=bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out2_jj[j]-=conj(A[ij])*in1_ii[i];
      //	}}

      // -C* : out_1j += (-Aji)^* in_2i 
      accum_n_conj_cmult(out1_jj0, C00, in2_ii0);
      accum_n_conj_cmult(out1_jj0, C10, in2_ii1);
      accum_n_conj_cmult(out1_jj1, C01, in2_ii0);
      accum_n_conj_cmult(out1_jj1, C11, in2_ii1);
      //      offset=2*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out1_jj[j]-=conj(A[ij])*in2_ii[i];
      //	}}
      
      // D* : out_2j += (Aji)^* in_2i 
      accum_conj_cmult(out2_jj0, D00, in2_ii0);
      accum_conj_cmult(out2_jj0, D10, in2_ii1);
      accum_conj_cmult(out2_jj1, D01, in2_ii0);
      accum_conj_cmult(out2_jj1, D11, in2_ii1);
      //      offset=3*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //      	for(int j=0; j<2; j++){
      //      	  int ij=bs2*i+j+offset;
      //      	  out2_jj[j]+=conj(A[ij])*in2_ii[i];
      //      	}}

      // store output
      store_block2vec(out1_ptr, out1_jj);
      store_block2vec(out2_ptr, out2_jj);

      A+=blocksize*blocksize;
    }
  }
}

static inline void K_coarse_daggered_hopp_float( complex_float *eta, const complex_float *phi,
						  const float *D, level_struct *l ) {
    
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // note: minus sign of D = self_coupling - hopping_term is added here
  
  //  A = [ a00 a01 a02..]
  //      [ a10 a22...   ]
  //      [...           ]  with 2x2 a?? 
  // stored as
  // [ (a00 b00 c00 d00) (a01 b01 c01 d01)...] each in row order
  //   

  static const int blocksize=4;
  int N=l->num_lattice_site_var;
  int n=N/2;
  int numblock=N/blocksize;
  //  int bs2=blocksize/2;

  const complex_float *A=(const complex_float*)D;
  for(int ib=0; ib<numblock; ib++){
    // input  ector
    //    const complex_float* in1_ii= phi    + ib*2;
    //    const complex_float* in2_ii=(phi+n) + ib*2;
    register _fjsp_v2r8 in1_ii0, in1_ii1;
    register _fjsp_v2r8 in2_ii0, in2_ii1;
    register block2vec *in1_ptr = (block2vec*)(phi    + ib*2);
    register block2vec *in2_ptr = (block2vec*)((phi+n)+ ib*2);
    load_block2vec(in1_ii, in1_ptr);
    load_block2vec(in2_ii, in2_ptr);

    for(int jb=0; jb<numblock; jb++){
      // output vector
      //      complex_float* out1_jj= eta    + jb*2;
      //      complex_float* out2_jj=(eta+n) + jb*2;
      //      int offset=0;
      register _fjsp_v2r8 out1_jj0, out1_jj1;
      register _fjsp_v2r8 out2_jj0, out2_jj1;
      register block2vec *out1_ptr = (block2vec*)(eta    + jb*2);
      register block2vec *out2_ptr = (block2vec*)((eta+n)+ jb*2);
      load_block2vec(out1_jj, out1_ptr);
      load_block2vec(out2_jj, out2_ptr);

      // matrix
      register _fjsp_v2r8 A00, A01, A10, A11;
      register _fjsp_v2r8 B00, B01, B10, B11;
      register _fjsp_v2r8 C00, C01, C10, C11;
      register _fjsp_v2r8 D00, D01, D10, D11;
      register block4mat *mat_ptr=(block4mat*)A;
      load_block4mat_ABCD(A,B,C,D, mat_ptr);

      // A* : out_1j -= (Aji)^* in_1i 
      accum_n_conj_cmult(out1_jj0, A00, in1_ii0);
      accum_n_conj_cmult(out1_jj0, A10, in1_ii1);
      accum_n_conj_cmult(out1_jj1, A01, in1_ii0);
      accum_n_conj_cmult(out1_jj1, A11, in1_ii1);
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j;
      //	  out1_jj[j]-=conj(A[ij])*in1_ii[i];
      //	}}

      // -B* : out_2j -= (-Aji)^* in_1i 
      accum_conj_cmult(out2_jj0, B00, in1_ii0);
      accum_conj_cmult(out2_jj0, B10, in1_ii1);
      accum_conj_cmult(out2_jj1, B01, in1_ii0);
      accum_conj_cmult(out2_jj1, B11, in1_ii1);
      //      offset=bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j + offset;
      //	  out2_jj[j]+=conj(A[ij])*in1_ii[i];
      //	}}

      // -C* : out_1j -= (-Aji)^* in_2i 
      accum_conj_cmult(out1_jj0, C00, in2_ii0);
      accum_conj_cmult(out1_jj0, C10, in2_ii1);
      accum_conj_cmult(out1_jj1, C01, in2_ii0);
      accum_conj_cmult(out1_jj1, C11, in2_ii1);
      //      offset=2*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j + offset;
      //	  out1_jj[j]+=conj(A[ij])*in2_ii[i];
      //	}}
      
      // D* : out_2j -= (Aji)^* in_2i 
      accum_n_conj_cmult(out2_jj0, D00, in2_ii0);
      accum_n_conj_cmult(out2_jj0, D10, in2_ii1);
      accum_n_conj_cmult(out2_jj1, D01, in2_ii0);
      accum_n_conj_cmult(out2_jj1, D11, in2_ii1);
      //      offset=3*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j + offset;
      //	  out2_jj[j]-=conj(A[ij])*in2_ii[i];
      //	}}
	
      // store output
      store_block2vec(out1_ptr, out1_jj);
      store_block2vec(out2_ptr, out2_jj);

      A+=blocksize*blocksize;
    }
  }
}


static inline void K_coarse_n_hopp_float( complex_float *eta, const complex_float *phi,
                                          const float *D, level_struct *l ) {
    
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // note: minus sign of D = self_coupling - hopping_term is added here
  
  //  A = [ a00 a01 a02..]
  //      [ a10 a22...   ]
  //      [...           ]  with 2x2 a?? 
  // stored as 
  // [ (a00 b00 c00 d00) (a01 b01 c01 d01)...] each in row order
  //   

  static const int blocksize=4;
  //  int bs2=blocksize/2;
  int N=l->num_lattice_site_var;
  int n=N/2;
  int numblock=N/blocksize;

  const complex_float *A=(const complex_float*)D;

  for(int ib=0; ib<numblock; ib++){
    // output vector
    //    complex_float* out1_ii= eta    + ib*2;
    //    complex_float* out2_ii=(eta+n) + ib*2;
    register _fjsp_v2r8 out1_ii0, out1_ii1;
    register _fjsp_v2r8 out2_ii0, out2_ii1;
    register block2vec *out1_ptr = (block2vec*)(eta    + ib*2);
    register block2vec *out2_ptr = (block2vec*)((eta+n)+ ib*2);
    load_block2vec(out1_ii, out1_ptr);
    load_block2vec(out2_ii, out2_ptr);

    for(int jb=0; jb<numblock; jb++){
      // input vector
      //      const complex_float* in1_jj= phi   + jb*2;
      //      const complex_float* in2_jj=(phi+n)+ jb*2;
      register _fjsp_v2r8 in1_jj0, in1_jj1;
      register _fjsp_v2r8 in2_jj0, in2_jj1;
      register block2vec *in1_ptr = (block2vec*)(phi    + jb*2);
      register block2vec *in2_ptr = (block2vec*)((phi+n)+ jb*2);
      load_block2vec(in1_jj, in1_ptr);
      load_block2vec(in2_jj, in2_ptr);

      // matrix
      register _fjsp_v2r8 A00, A01, A10, A11;
      register _fjsp_v2r8 B00, B01, B10, B11;
      register _fjsp_v2r8 C00, C01, C10, C11;
      register _fjsp_v2r8 D00, D01, D10, D11;
      register block4mat *mat_ptr=(block4mat*)A;
      load_block4mat_ABCD(A,B,C,D, mat_ptr);

      // A : out_1i += (Aij) in_1j
      accum_cmult(out1_ii0, A00, in1_jj0);
      accum_cmult(out1_ii0, A01, in1_jj1);
      accum_cmult(out1_ii1, A10, in1_jj0);
      accum_cmult(out1_ii1, A11, in1_jj1);
      //      int offset=0;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j;
      //	  out1_ii[i]+=A[ij]*in1_jj[j];
      //	}}

      // B : out_1i += (Aij) in_2j 
      accum_cmult(out1_ii0, B00, in2_jj0);
      accum_cmult(out1_ii0, B01, in2_jj1);
      accum_cmult(out1_ii1, B10, in2_jj0);
      accum_cmult(out1_ii1, B11, in2_jj1);
      //      offset=bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out1_ii[i]+=A[ij]*in2_jj[j];
      //	}}

      // C : out_2i += (Aij) in_1j 
      accum_cmult(out2_ii0, C00, in1_jj0);
      accum_cmult(out2_ii0, C01, in1_jj1);
      accum_cmult(out2_ii1, C10, in1_jj0);
      accum_cmult(out2_ii1, C11, in1_jj1);
      //      offset=2*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out2_ii[i]+=A[ij]*in1_jj[j];
      //	}}
      
      // D : out_2j += (Aij) in_2i 
      accum_cmult(out2_ii0, D00, in2_jj0);
      accum_cmult(out2_ii0, D01, in2_jj1);
      accum_cmult(out2_ii1, D10, in2_jj0);
      accum_cmult(out2_ii1, D11, in2_jj1);
      //      offset=3*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out2_ii[i]+=A[ij]*in2_jj[j];
      //	}}
	
      A+=blocksize*blocksize;
    }
    // store output
    store_block2vec(out1_ptr, out1_ii);
    store_block2vec(out2_ptr, out2_ii);
  }
}

static inline void K_coarse_hopp_float( complex_float *eta, const complex_float *phi,
                                        const float *D, level_struct *l ) {
    
  // U_mu(x) = [ A B      , U_-mu(x+muhat) = [ A* -C*
  //             C D ]                        -B*  D* ]
  // note: minus sign of D = self_coupling - hopping_term is added here
  
  //  A = [ a00 a01 a02..]
  //      [ a10 a22...   ]
  //      [...           ]  with 2x2 a?? 
  // stored as 
  // [ (a00 b00 c00 d00) (a01 b01 c01 d01)...] each in row order
  //   

  static const int blocksize=4;
  int N=l->num_lattice_site_var;
  int n=N/2;
  int numblock=N/blocksize;
  //  int bs2=blocksize/2;

  const complex_float *A=(const complex_float*)D;

  for(int ib=0; ib<numblock; ib++){
    // output vector
    //    complex_float* out1_ii= eta    + ib*2;
    //    complex_float* out2_ii=(eta+n) + ib*2;
    register _fjsp_v2r8 out1_ii0, out1_ii1;
    register _fjsp_v2r8 out2_ii0, out2_ii1;
    register block2vec *out1_ptr = (block2vec*)(eta    + ib*2);
    register block2vec *out2_ptr = (block2vec*)((eta+n)+ ib*2);
    load_block2vec(out1_ii, out1_ptr);
    load_block2vec(out2_ii, out2_ptr);

    for(int jb=0; jb<numblock; jb++){
      // input vector
      //      const complex_float* in1_jj= phi   + jb*2;
      //      const complex_float* in2_jj=(phi+n)+ jb*2;
      register _fjsp_v2r8 in1_jj0, in1_jj1;
      register _fjsp_v2r8 in2_jj0, in2_jj1;
      register block2vec *in1_ptr = (block2vec*)(phi    + jb*2);
      register block2vec *in2_ptr = (block2vec*)((phi+n)+ jb*2);
      load_block2vec(in1_jj, in1_ptr);
      load_block2vec(in2_jj, in2_ptr);

      // matrix
      register _fjsp_v2r8 A00, A01, A10, A11;
      register _fjsp_v2r8 B00, B01, B10, B11;
      register _fjsp_v2r8 C00, C01, C10, C11;
      register _fjsp_v2r8 D00, D01, D10, D11;
      register block4mat *mat_ptr=(block4mat*)A;
      load_block4mat_ABCD(A,B,C,D, mat_ptr);

      // A : out_1i -= (Aij) in_1j
      accum_n_cmult(out1_ii0, A00, in1_jj0);
      accum_n_cmult(out1_ii0, A01, in1_jj1);
      accum_n_cmult(out1_ii1, A10, in1_jj0);
      accum_n_cmult(out1_ii1, A11, in1_jj1);
      //      int offset=0;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j;
      //	  out1_ii[i]-=A[ij]*in1_jj[j];
      //	}}

      // B : out_1i -= (Aij) in_2j 
      accum_n_cmult(out1_ii0, B00, in2_jj0);
      accum_n_cmult(out1_ii0, B01, in2_jj1);
      accum_n_cmult(out1_ii1, B10, in2_jj0);
      accum_n_cmult(out1_ii1, B11, in2_jj1);
      //      offset=bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out1_ii[i]-=A[ij]*in2_jj[j];
      //	}}

      // C : out_2i -= (Aij) in_1j 
      accum_n_cmult(out2_ii0, C00, in1_jj0);
      accum_n_cmult(out2_ii0, C01, in1_jj1);
      accum_n_cmult(out2_ii1, C10, in1_jj0);
      accum_n_cmult(out2_ii1, C11, in1_jj1);
      //      offset=2*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out2_ii[i]-=A[ij]*in1_jj[j];
      //	}}
      
      // D : out_2j -= (Aij) in_2i 
      accum_n_cmult(out2_ii0, D00, in2_jj0);
      accum_n_cmult(out2_ii0, D01, in2_jj1);
      accum_n_cmult(out2_ii1, D10, in2_jj0);
      accum_n_cmult(out2_ii1, D11, in2_jj1);
      //      offset=3*bs2*bs2;
      //      for(int i=0; i<2; i++){
      //	for(int j=0; j<2; j++){
      //	  int ij=bs2*i+j+offset;
      //	  out2_ii[i]-=A[ij]*in2_jj[j];
      //	}}
	
      A+=blocksize*blocksize;
    }
    // store output
    store_block2vec(out1_ptr, out1_ii);
    store_block2vec(out2_ptr, out2_ii);
  }
}

#endif  // K_OPT
