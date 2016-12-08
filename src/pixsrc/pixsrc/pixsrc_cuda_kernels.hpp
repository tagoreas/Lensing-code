#ifndef PS_SIT
#ifdef INTEGER_TYPE_LL
typedef long long PS_SIT;
typedef unsigned long long PS_unsignedSIT;
#endif
#ifdef INTEGER_TYPE_L
typedef long PS_SIT;
typedef unsigned long PS_unsignedSIT;
#endif
#ifdef INTEGER_TYPE_I
typedef int PS_SIT;
typedef unsigned int PS_unsignedSIT;
#endif
#endif

#ifndef PS_FPT
#ifdef SINGLE_PRECISION
typedef float PS_FPT;
#endif
#ifdef DOUBLE_PRECISION
typedef double PS_FPT;
#endif
#endif



#ifndef PIXSRC_CUDA_KERNELS_HPP_
#define PIXSRC_CUDA_KERNELS_HPP_

// for CUDA
#include <driver_types.h>


void   psc_fill                   ( PS_FPT*, PS_FPT, int, cudaStream_t* );
void   psc_cpy_mat_col_to_vec     ( int, int*, int*, PS_FPT*, PS_FPT* );
void   psc_undo_cpy_mat_col_to_vec( int, int*, int*, PS_FPT* );
void   psc_cpy_mat_col_to_vec     ( int, int, int*, int*, PS_FPT*, PS_FPT* );
PS_FPT psc_sum_vector             ( PS_FPT*, int );
void   psc_tri_eq_solve_cusp      ( int*, int*, PS_FPT*, PS_FPT*, PS_FPT*, int, int, void* );
void   psc_get_preconditioner_cusp( int*, int*, PS_FPT*, int, int, void** );
void   psc_destr_preconditioner_cusp( void* );
void   psc_tr_inv_a_b_matrixfree  ( int*, int*, PS_FPT*, int*, int*, PS_FPT*,
                                    PS_FPT*, PS_FPT*, int, void*, void*, int, void*         );

#endif
