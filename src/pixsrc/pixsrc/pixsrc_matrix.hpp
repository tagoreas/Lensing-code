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



/////////////////////////////////////////////////////////////////////
// NOTE: a lot of the variables here are not currently being used. //
//       See pixsrc_matrix.cpp for details of what the below do    //
//       and which are actually used.                              //
/////////////////////////////////////////////////////////////////////



#ifndef PIXSRC_MATRIX_HPP_
#define PIXSRC_MATRIX_HPP_

#define MATRIX pixsrc_matrix

#include "pixsrc_inputdata.hpp"
#include "pixsrc_vector.hpp"
#include <pthread.h>

#include <string>
using std::string;

// flag for where the matrix data is up to date
typedef enum
{
    pms_neither = -1,
    pms_both    = 0,
    pms_cpu     = 1,
    pms_gpu     = 2

} pm_update_status;

// flag for the compression format of the matrix
typedef enum
{
    pmf_csr = 0,
    pmf_csc = 1,
    pmf_dense = 2

} pm_format;


// matrix class
class pixsrc_matrix
{

    // let cuda access private variables directly
#ifdef __USE_PIXSRC_CUDA__
    friend class pixsrc_cuda;
#endif


public:

    pixsrc_matrix          ( commoninputdata*, inputdata*,
                             PS_SIT, PS_SIT, PS_SIT, PS_SIT, inputdata*);
    ~pixsrc_matrix         (                                           );
    void pm_init           ( commoninputdata*, inputdata*,
                             PS_SIT, PS_SIT, PS_SIT, PS_SIT        );

    bool are_we_using_cuda;

    void   init_cpu_vecs   (                                           );
    double trace_invA_B    ( pixsrc_matrix*, PS_SIT, PS_SIT*                 );
    void   noise_invA      ( pixsrc_vector*, PS_SIT, PS_SIT*                 );
    void   noise_invA_fullmat      ( double*, PS_SIT, PS_SIT*                 );
    void   set             ( PS_SIT, PS_SIT, double                          );
    void   set             ( PS_SIT, double                               );
    double get             ( PS_SIT, PS_SIT                            );
    PS_SIT get_nnz         (                                           );
    void   mult            ( double                                    );
    void   mult            ( pixsrc_vector*, pixsrc_vector*,
                             bool, PS_SIT                                 );
    void   submultadd      (pixsrc_matrix *b, pixsrc_matrix *c, bool tr1, bool tr2,
                            PS_SIT am, PS_SIT an, PS_SIT bm, PS_SIT bn,
                            PS_SIT cm, PS_SIT cn, PS_SIT cstartm, PS_SIT cstartn, PS_SIT numprocesses);
    void   inv_dense_mult  ( pixsrc_vector*, pixsrc_vector*, PS_SIT       );
    void   inv_matrix_mult ( pixsrc_matrix*, pixsrc_matrix*, PS_SIT);
    void   mult            ( pixsrc_matrix*, pixsrc_matrix*,
                             bool, bool, PS_SIT                           );
    void   plus            ( pixsrc_matrix*, pixsrc_matrix*,
                             bool, bool, double, double, PS_SIT           );
    double logdet          (                                           );
    void   linequationsolve( pixsrc_vector*, pixsrc_vector*, bool      );
    void   update_gpu      ( pm_format                                 );
    void   update_cpu      ( pm_format                                 );

    void   atasbc          ( double*, PS_SIT                              );
    void   print_me        ( PS_SIT, string, string, bool, string, bool     );

    void   dissolve_scalar (                                           );
    void   remove_duplicates();
    void   inplace_inverse ();
    void   inplace_lu_fac ();
    void   inplace_transpose ();
    void   remove_last_row_col ();
    void   inplace_leq (pixsrc_matrix*);


    /*
     * void   copycolumn      ( PS_SIT, pixsrc_vector*                       );
     * void   copymatrix      ( pixsrc_matrix*                            );
     */

//private:

    bool is_dense;
    bool got_det;
    bool got_lu_fac;
    bool got_inv;
    int/*PS_SIT*/ *ipiv_dense;
    PS_FPT *mat_dense;
    PS_FPT *inv_dense;
    PS_FPT *lu_fac_dense;
    PS_FPT det_dense;
    pthread_mutex_t dense_lu_mutex;
    pthread_mutex_t dense_inv_mutex;
    void cleardense ();
    void resetdense (PS_SIT, PS_SIT);
    void dense_matrix_stuff (MATRIX*, MATRIX*, VECTOR*, VECTOR*, double*, PS_SIT);

    ////////////////////////////////////////////////////
    /////////////////// MATRIX TYPE ////////////////////
    ////////////////////////////////////////////////////
    //// 0:  general                                ////
    //// 1:  symmetric                              ////
    //// 2:  Hermitian  (upper)                     ////
    //// 3:  triangular (upper, non-unity diagonal) ////
    //// 4:  triangular (upper, unity diagonal)     ////
    //// 10: diagonal                               ////
    //// 11: proportional to identity               ////
    ////////////////////////////////////////////////////
    PS_SIT matrix_type;

    commoninputdata *cdata_;
    inputdata *data_;

    pm_update_status up2date_csr;
    pm_update_status up2date_csc;

    PS_FPT  scalar;              // scalar to multiply matrix by;
    PS_SIT  ncol, nrow, nnz;     // # of {columns, rows, nonzeroes}
    bool    rcform, ludone;
    bool    hasbeencompiled;
    PS_SIT *apc, *apr;           // ap[1+1]-ap[i] = num of entries in col[i]/row[i]
    PS_SIT *aic, *air;           // row/column indices
    PS_FPT *axc, *axr;           // matrix entries
    void   *symbolic, *numeric;  // used by umfpack;
    PS_SIT vecsize;             // size of below
    PS_SIT *tiv, *tjv;           // triplets - row, column
    PS_FPT *tkv;                 // triplets

    bool is_sym_good;
    PS_SIT *air_sym, *apr_sym;
    PS_FPT *axr_sym;

    PS_SIT *icf_ptr_cpu;
    PS_SIT *icf_ind_cpu;
    PS_FPT *icf_val_cpu;


#ifdef __USE_PIXSRC_CUDA__

    // cuda temporary vectors (for working)
    PS_SIT num_ctv;
    PS_FPT **ctv;

    // CUSPARSE CSR and CSC format pointers
    PS_SIT nnz_chol;
    PS_SIT    *coo_col;
    PS_FPT *coo_val;
    PS_SIT    *csx;
    PS_SIT    *coo_row_tr;
    PS_FPT *coo_val_tr;
    PS_SIT    *csx_tr;

    // for CHOLMOD
    void *chol_fac;
    void *chol_com;
    void *chol_mat;

    /*
    // These vectors are for storage of symmetric matrices
    PS_SIT *coo_col_sym;
    PS_SIT *csx_sym;
    PS_FPT *coo_val_sym;
    */

    // incomplete Cholesky factorization. The following PS_FPT array
    // corresponds to coo_col and csx. The voids correspond to the
    // CUSPARSE analysis of the upper and lower triangular matrices
    PS_SIT *icf_ptr;
    PS_SIT *icf_ind;
    PS_FPT *icf_val;
    bool    icf_done;
    void *icf_upper_sva;
    void *icf_lower_sva;

    // CUSPARSE CSRSM analysis info pointer
    //void *csrsvi;
    //bool  csrsvi_done;

    // CUSPARSE matrix descriptor
    //void *descr;

    // CUDA mutexes
    pthread_mutex_t pm_mutex;
    //pthread_mutex_t csrsvi_mutex;

    // CUDA related functions
    void clearcudavec();
    void cuda_reset  ();

    // CUSP object
    void *cusp_prec;

    // CUDA stream
    void *stream;

#endif

    pthread_attr_t  attr;
    pthread_mutex_t compilemutex;
    pthread_mutex_t dissolvemutex;
    pthread_mutex_t numsymbmutex;
    pthread_mutex_t rcformmutex;

    pixsrc_matrix(const pixsrc_matrix &m);
    pixsrc_matrix operator=(const pixsrc_matrix &ps);

    void resize(PS_SIT**,PS_SIT**,PS_FPT**,PS_SIT,PS_SIT*);
    static void *trace_invA_B_body( void* );
    static void *noise_invA_body( void* );
    static void *noise_invA_body_fullmat( void* );
    static void *multmatrixmatrix ( void* );
    static void *multmatrixvector ( void* );

    void compile();
    void creatercform();
    void createnumsymb();
    void clearrows();
    void clearcols();
    void clearvoids();
    void clearvecs();
    void clearall();

    void   selective_copy  ( PS_SIT*, PS_SIT*, PS_FPT*,
                             PS_SIT*, PS_SIT*, PS_FPT*, PS_SIT                  );
};

#endif
