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



#ifdef __USE_PIXSRC_CUDA__

#ifndef PIXSRC_CUDA_HPP_
#define PIXSRC_CUDA_HPP_

#define CUDA pixsrc_cuda::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_vector.hpp"
#include "pixsrc_matrix.hpp"

#include <driver_types.h>

// flag that indicates in which direction data is to be copied
// either from CPU to GPU or the other way
typedef enum
{
    ccd_cpu2gpu = 0,
    ccd_gpu2cpu = 1

} cuda_cpy_direction;


class pixsrc_cuda
{
public:

    static void detectgpus     ( commoninputdata*                                  );
    static void create_handle  ( commoninputdata*, int                             );
    static void mult           ( pixsrc_matrix*, pixsrc_vector*, pixsrc_vector*,
                                 bool                                              );
    static void mult           ( pixsrc_matrix*, pixsrc_matrix*, pixsrc_matrix*,
                                 bool, bool                                        );
    static void plus           ( pixsrc_matrix*, pixsrc_matrix*, pixsrc_matrix*,
                                 bool, bool, PS_FPT, PS_FPT                        );
    static void plus           ( pixsrc_vector*, pixsrc_vector*, pixsrc_vector*,
                                 PS_FPT, PS_FPT                                    );
    static void minus          ( pixsrc_vector*, pixsrc_vector*, pixsrc_vector*,
                                 PS_FPT, PS_FPT                                    );
    static PS_FPT innerproduct ( pixsrc_vector*, pixsrc_vector*                    );
    static void set_device     ( int                                               );
    static void benchmark      ( commoninputdata*                                  );
    static void ps_start_clock ( void**                                            );
    static void ps_stop_clock  ( void*, float*                                     );
    static void creatercform   ( pixsrc_matrix*                                    );
    static void compile        ( pixsrc_matrix*                                    );
    static void ps_cuda_free   ( int*                                              );
    static void ps_cuda_free   ( PS_FPT*                                           );
    static void fill_vec       ( pixsrc_vector*, PS_FPT                            );
    static PS_FPT ps_cuda_max  ( pixsrc_vector*                                    );
    static void dissolve_scalar( pixsrc_vector*                                    );
    static void dissolve_scalar( pixsrc_matrix*, pm_format                         );
    static void tr_inv_a_b  ( pixsrc_matrix*, pixsrc_matrix*, PS_FPT*           );



    static void sync_vector   ( pixsrc_vector*, cuda_cpy_direction                 );
    static void sync_matrix   ( pixsrc_matrix*, pm_format, cuda_cpy_direction      );
    static void ps_cuda_malloc ( commoninputdata*, int**,    int, int              );
    static void ps_cuda_malloc ( commoninputdata*, PS_FPT**, int, int              );

    static void tri_eq_solve  ( pixsrc_matrix*, pixsrc_vector*, pixsrc_vector*, bool );

    static void wait_for_stream ( MATRIX*                                          );
    static void wait_for_stream ( VECTOR*                                          );
    static void create_stream ( MATRIX*                                            );
    static void create_stream ( VECTOR*                                            );
    static void destroy_stream( MATRIX*                                            );
    static void destroy_stream( VECTOR*                                            );
    static void destroy_cusp_prec( MATRIX*                                         );

    // currently unused
    /*
     * static void ps_wrap_cusp   ( pixsrc_matrix*                                 );
     * static void ps_wrap_cusp   ( pixsrc_vector*                                 );
     * static void destroy_cusp_wrapper   ( pixsrc_matrix*                         );
     * static void destroy_cusp_wrapper   ( pixsrc_vector*                         );
     * static void inc_chol_factor_cusparse( pixsrc_matrix*                        );
     * static void tri_eq_solve   ( pixsrc_matrix*, pixsrc_vector*, pixsrc_vector* );
     * static void create_sym_form( pixsrc_matrix*                                 );
     * static void ps_dest_csrsvi ( pixsrc_matrix*                                 );
     * static void ps_set_descr   ( pixsrc_matrix*                                 );
     * static void ps_dest_descr  ( pixsrc_matrix*                                 );
     * static void copy_csc_c2g_h ( pixsrc_matrix*                                 );
     * static void copy_csr_c2g_s ( pixsrc_matrix*                                 );
     * static void inc_chol_factor( pixsrc_matrix*                                 );
     * static void cusparse_tr_inva_b  ( pixsrc_matrix*, pixsrc_matrix*, PS_FPT*   );
     * static void analyze_upp_low( pixsrc_matrix*                                 );
     * static void csr2coo        ( pixsrc_matrix*                                 );
     * static void *tr_invA_B_kernel( void*                                        );
     */

private:

    static void ps_cuda_memcpy ( commoninputdata*, int*, int*, int, int, cudaStream_t*, int       );
    static void ps_cuda_memcpy ( commoninputdata*, PS_FPT*, PS_FPT*, int, int, cudaStream_t*, int );
    static void errorcheck     ( commoninputdata*, void*, int, const char*         );
    static void postprocess    ( pixsrc_matrix*, pixsrc_matrix*, pixsrc_matrix*    );
    static void *tr_inv_a_b_kernel( void*                                     );
};

#endif

#endif
