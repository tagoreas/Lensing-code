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



///
/// The __USE_PIXSRC_CUDA__ variable passed on to CXX
/// decides whether this source code will be compiled
///

#include "pixsrc_cuda.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_printer.hpp"
#include <cstring>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>

#include "pixsrc_cuda_kernels.hpp"

#include <pthread.h>

PS_FPT pixsrc_cuda::innerproduct( VECTOR *a, VECTOR *b )
{
    // compute inner product between two vectors

    // update GPU data
    a->update_gpu();
    b->update_gpu();

    // set GPU device and CUDA handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );

    commoninputdata *cdata_ = a->cdata_;
    cublasStatus_t stat;
    cublasHandle_t *handle_cublas  = (cublasHandle_t*)cdata_->gpudevices[dev].handle_cublas;

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cublas_lock );

    // going to use VECTOR a's stream for calculations
    // so we have to wait for it to finish afterwards too.
    cudaStream_t *stream = (cudaStream_t*)a->stream;
    cublasSetStream( *handle_cublas, *stream );

    PS_FPT res;

    // compute inner product
    if( a == b )
    {

#ifdef SINGLE_PRECISION
        stat = cublasSnrm2( *handle_cublas, a->size, a->cuda_vec, 1, &res );
#else
        stat = cublasDnrm2( *handle_cublas, a->size, a->cuda_vec, 1, &res );
#endif

        // square because the norm was computed above
        res *= res;
    }
    else
    {

#ifdef SINGLE_PRECISION
        stat = cublasSdot ( *handle_cublas, a->size, a->cuda_vec, 1,
                            b->cuda_vec, 1, &res );
#else
        stat = cublasDdot ( *handle_cublas, a->size, a->cuda_vec, 1,
                            b->cuda_vec, 1, &res );
#endif

    }

    // error check
    CUDA errorcheck( cdata_, (void*)&stat, 2, "v-v-p 1" );

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cublas_lock );

    // update (update will wait for GPU to finish, so no explicit call is needed)
    a->update_gpu();

    return res * a->scalar * b->scalar;
}

void pixsrc_cuda::plus( VECTOR *a, VECTOR *b, VECTOR *c, PS_FPT alpha, PS_FPT beta )
{
    // add two vectors

    // update GPU data
    a->update_gpu();
    b->update_gpu();

    // set GPU device and handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;
    cublasStatus_t stat;
    cublasHandle_t *handle_cublas  = (cublasHandle_t*)cdata_->gpudevices[dev].handle_cublas;

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cublas_lock );

    // result will be stores in a if c is NULL
    VECTOR *res = (c) ? c : a;

    // use res's stream for calculations
    cudaStream_t *stream = (cudaStream_t*)res->stream;
    cublasSetStream( *handle_cublas, *stream );

    if (c)
    {

        // copy a into c
#ifdef SINGLE_PRECISION
        stat = cublasScopy( *handle_cublas, a->size,
                            a->cuda_vec, 1, c->cuda_vec, 1 );
#else
        stat = cublasDcopy( *handle_cublas, a->size,
                            a->cuda_vec, 1, c->cuda_vec, 1 );
#endif

        // error check
        CUDA errorcheck( cdata_, (void*)&stat, 2, "v-v-p copy" );
    }

    alpha *= a->scalar;
    beta  *= b->scalar;
    res->scalar = 1;

    if( alpha != 1.0 )
    {

        // scale res by alpha before sum is performed
#ifdef SINGLE_PRECISION
        stat = cublasSscal( *handle_cublas, a->size,
                            &alpha, res->cuda_vec, 1 );
#else
        stat = cublasDscal( *handle_cublas, a->size,
                            &alpha, res->cuda_vec, 1 );
#endif

        // error check
        CUDA errorcheck( cdata_, (void*)&stat, 2, "v-scale" );
    }

    // add b to res
#ifdef SINGLE_PRECISION
    stat = cublasSaxpy( *handle_cublas, a->size, &beta,
                        b->cuda_vec, 1, res->cuda_vec, 1 );
#else
    stat = cublasDaxpy( *handle_cublas, a->size, &beta,
                        b->cuda_vec, 1, res->cuda_vec, 1 );
#endif

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cublas_lock );

    // error check
    CUDA errorcheck( cdata_, (void*)&stat, 2, "v-v-p 2" );

    // update
    res->up2date = pvs_gpu;
}

void pixsrc_cuda::minus( VECTOR *a, VECTOR *b, VECTOR *c, PS_FPT alpha, PS_FPT beta )
{
    // wrapper
    CUDA plus( a, b, c, alpha, -beta );
}

void pixsrc_cuda::mult( MATRIX *a, MATRIX *b, MATRIX *c, bool tr1, bool tr2 )
{
    // multiply two matrices together

    // set GPU device
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;

    if( tr2 )
        PRINTER printerror( "cuda", "cuda error: can't transpose second matrix"
                            " in matrix-matrix multiplication yet",
                            cdata_->print2screenmutex );

    // I handle transposing the matrix myself, since CSR is the transpose of CSC.
    // So, figure out format of matrix that we need
    pm_format format;
    format = ( tr1 ) ? pmf_csc : pmf_csr;
    a->update_gpu( format );
    format = ( tr2 ) ? pmf_csc : pmf_csr;
    b->update_gpu( format );

    cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

    cusparseHandle_t   *handle   = (cusparseHandle_t  *)cdata_->gpudevices[dev].handle;

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cusparse_lock );

    // use c's stream for calculations
    cudaStream_t *stream = (cudaStream_t*)c->stream;
    cusparseSetStream( *handle, *stream );

    // get matrix desciptors
    cusparseMatDescr_t *descr_a  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;
    cusparseMatDescr_t *descr_b  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;
    cusparseMatDescr_t *descr_c  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;

    // these definitions os trans_? are overridden below.
    cusparseOperation_t trans_a = ( tr1 ) ?
        CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
    cusparseOperation_t trans_b = ( tr2 ) ?
        CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;

    // I transpose the matrix myself instead of letting
    // CUDA do it it internally every time I perform
    // a calculation. I think doing it myself and
    // saving the result of the transposing step
    // is faster for multiple operations.
    trans_a = trans_b = CUSPARSE_OPERATION_NON_TRANSPOSE;

    // reset matrix to put answer into
    c->cuda_reset();

    // rows and columns of the matrices
    int m = c->nrow;
    int n = c->ncol;
    int k = ( tr1 ) ? a->nrow : a->ncol;


    // actual multiplication
    if( !tr1 && !tr2 )
    {
        if( !(c->csx) )
        {
            // allocate memory for C's CSR row pointer on GPU
            // and figure out number of nonzero entries
            CUDA ps_cuda_malloc( cdata_, &c->csx, m+1, dev );

            status = cusparseXcsrgemmNnz( *handle, trans_a, trans_b,
                                          m, n, k,
                                          *descr_a, a->nnz, a->csx, a->coo_col,
                                          *descr_b, b->nnz, b->csx, b->coo_col,
                                          *descr_c,         c->csx, &c->nnz     );

            CUDA errorcheck( cdata_, (void*)&status, 1, "first pass" );

        }

        // allocate memory for C's tjv and tkv if already haven't
        if( !c->coo_col )
            CUDA ps_cuda_malloc( cdata_, &c->coo_col, c->nnz, dev );

        if( !c->coo_val )
            CUDA ps_cuda_malloc( cdata_, &c->coo_val, c->nnz, dev );

        // perform multiplication

#ifdef SINGLE_PRECISION
        status = cusparseScsrgemm( *handle, trans_a, trans_b,
                                   m, n, k,
                                   *descr_a, a->nnz,
                                   a->coo_val, a->csx, a->coo_col,
                                   *descr_b, b->nnz,
                                   b->coo_val, b->csx, b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col  );
#else
        status = cusparseDcsrgemm( *handle, trans_a, trans_b,
                                   m, n, k,
                                   *descr_a, a->nnz,
                                   a->coo_val, a->csx, a->coo_col,
                                   *descr_b, b->nnz,
                                   b->coo_val, b->csx, b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col  );
#endif

        CUDA errorcheck( cdata_, (void*)&status, 1, "matrix-matrix multiplication" );

    }
    else if( tr1 && !tr2 )
    {
        if( !(c->csx) )
        {
            // allocate memory for C's CSR row pointer on GPU
            // and figure out number of nonzero entries
            CUDA ps_cuda_malloc( cdata_, &c->csx, m+1, dev );

            status = cusparseXcsrgemmNnz( *handle, trans_a, trans_b,
                                          m, n, k,
                                          *descr_a, a->nnz, a->csx_tr, a->coo_row_tr,
                                          *descr_b, b->nnz, b->csx,    b->coo_col,
                                          *descr_c,         c->csx,   &c->nnz         );

            CUDA errorcheck( cdata_, (void*)&status, 1, "first pass" );

        }

        // allocate memory for C's tjv and tkv if already haven't
        if( !c->coo_col )
            CUDA ps_cuda_malloc( cdata_, &c->coo_col, c->nnz, dev );

        if( !c->coo_val )
            CUDA ps_cuda_malloc( cdata_, &c->coo_val, c->nnz, dev );

        // perform multiplication

#ifdef SINGLE_PRECISION
        status = cusparseScsrgemm( *handle, trans_a, trans_b,
                                   m, n, k,
                                   *descr_a, a->nnz,
                                   a->coo_val_tr, a->csx_tr, a->coo_row_tr,
                                   *descr_b, b->nnz,
                                   b->coo_val,    b->csx,    b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col        );
#else
        status = cusparseDcsrgemm( *handle, trans_a, trans_b,
                                   m, n, k,
                                   *descr_a, a->nnz,
                                   a->coo_val_tr, a->csx_tr, a->coo_row_tr,
                                   *descr_b, b->nnz,
                                   b->coo_val,    b->csx,    b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col        );
#endif

        CUDA errorcheck( cdata_, (void*)&status, 1, "matrix-matrix multiplication" );

    }

    pthread_mutex_unlock( cdata_->gpudevices[dev].cusparse_lock );

    // do post-processing of matrix c
    c->scalar = a->scalar * b->scalar;
    CUDA postprocess( a, b, c );
}

void pixsrc_cuda::mult( MATRIX *a, VECTOR *b, VECTOR *c, bool tr )
{
    // multiply matrix and vector

    // set GPU device and handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;

    cusparseStatus_t status   = CUSPARSE_STATUS_SUCCESS;

    // I handle transposing the matrix myself, since CSR is the transpose of CSC.
    // So, figure out format of matrix that we need
    pm_format format = ( tr ) ? pmf_csc : pmf_csr;
    a->update_gpu( format );
    // update GPU data
    b->update_gpu();

    cusparseHandle_t   *handle   = (cusparseHandle_t  *)cdata_->gpudevices[dev].handle;

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cusparse_lock );

    // use c's stream for calculations
    cudaStream_t *stream = (cudaStream_t*)c->stream;
    cusparseSetStream( *handle, *stream );

    // matrix descriptor
    cusparseMatDescr_t *descr_a  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;

    // I transpose the matrix myself instead of letting
    // CUDA do it it internally every time I perform
    // a calculation. I think doing it myself and
    // saving the result of the transposing step
    // is faster for multiple operations.
    cusparseOperation_t trans_a = CUSPARSE_OPERATION_NON_TRANSPOSE;

    // reset vector c to put answer into
    c->cuda_reset();

    // rows and columns of matrix op(a)
    // "op" refers to the transposition or conjugation operation
    int m = c->size;
    int n = b->size;


    // actual multiplication

    PS_FPT dummy_zero = 0;
    PS_FPT dummy_one  = 1;

    if( !tr )
    {

#ifdef SINGLE_PRECISION
        status = cusparseScsrmv( *handle, trans_a,
                                 m, n, a->nnz, &dummy_one,
                                 *descr_a,
                                 a->coo_val, a->csx, a->coo_col,
                                 b->cuda_vec,
                                 &dummy_zero, c->cuda_vec        );
#else
        status = cusparseDcsrmv( *handle, trans_a,
                                 m, n, a->nnz, &dummy_one,
                                 *descr_a,
                                 a->coo_val, a->csx, a->coo_col,
                                 b->cuda_vec,
                                 &dummy_zero, c->cuda_vec        );
#endif

        CUDA errorcheck( cdata_, (void*)&status, 1, "M-V mult 0" );
    }
    else if( tr )
    {

#ifdef SINGLE_PRECISION
        status = cusparseScsrmv( *handle, trans_a,
                                 m, n, a->nnz, &dummy_one,
                                 *descr_a,
                                 a->coo_val_tr, a->csx_tr, a->coo_row_tr,
                                 b->cuda_vec,
                                 &dummy_zero, c->cuda_vec                 );
#else
        status = cusparseDcsrmv( *handle, trans_a,
                                 m, n, a->nnz, &dummy_one,
                                 *descr_a,
                                 a->coo_val_tr, a->csx_tr, a->coo_row_tr,
                                 b->cuda_vec,
                                 &dummy_zero, c->cuda_vec                 );
#endif

        CUDA errorcheck( cdata_, (void*)&status, 1, "M-V mult 1" );
    }

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cusparse_lock );

    // update
    c->scalar = a->scalar * b->scalar;
    c->up2date = pvs_gpu;
}

void pixsrc_cuda::plus( MATRIX *a, MATRIX *b, MATRIX *c, bool tr1, bool tr2,
                        PS_FPT scalar1, PS_FPT scalar2                       )
{
    // add two matrices

    // I handle transposing the matrix myself, since CSR is the transpose of CSC.
    // So, figure out format of matrix that we need
    pm_format format;
    format = ( tr1 ) ? pmf_csc : pmf_csr;
    a->update_gpu( format );
    format = ( tr2 ) ? pmf_csc : pmf_csr;
    b->update_gpu( format );

    // set GPU device and handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;

    if( tr2 )
        PRINTER printerror( "cuda", "cuda error: can't transpose second matrix"
                            " in matrix-matrix addition yet",
                            cdata_->print2screenmutex );

    cusparseStatus_t status = CUSPARSE_STATUS_SUCCESS;

    cusparseHandle_t   *handle   = (cusparseHandle_t  *)cdata_->gpudevices[dev].handle;

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cusparse_lock );

    // use c's stream for calculations
    cudaStream_t *stream = (cudaStream_t*)c->stream;
    cusparseSetStream( *handle, *stream );

    // matrix descriptors
    cusparseMatDescr_t *descr_a  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;
    cusparseMatDescr_t *descr_b  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;
    cusparseMatDescr_t *descr_c  = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;

    //cusparseOperation_t trans_a = ( tr1 ) ?
    //    CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;
    //cusparseOperation_t trans_b = ( tr2 ) ?
    //    CUSPARSE_OPERATION_TRANSPOSE : CUSPARSE_OPERATION_NON_TRANSPOSE;

    // I transpose the matrix myself instead of letting
    // CUDA do it it internally every time I perform
    // a calculation. I think doing it myself and
    // saving the result of the transposing step
    // is faster for multiple operations.
    // Turns out, CUSPARSE doesn't do transposed matrices
    // on its own.
    //trans_a = trans_b = CUSPARSE_OPERATION_NON_TRANSPOSE;

    // reset matrix to put answer into
    c->cuda_reset();

    // rows and columns of the matrices
    int m = c->nrow;
    int n = c->ncol;

    // factors multiplying MATRIX a & b
    PS_FPT alpha = scalar1 * a->scalar;
    PS_FPT beta  = scalar2 * b->scalar;

    // actual multiplication

    if( !tr1 && !tr2 )
    {
        if( !c->csx )
        {
            // allocate memory for C's CSR row pointer on GPU
            // and figure out number of nonzero entries
            CUDA ps_cuda_malloc( cdata_, &c->csx, m+1, dev );

            status = cusparseXcsrgeamNnz( *handle, m, n,
                                          *descr_a, a->nnz, a->csx, a->coo_col,
                                          *descr_b, b->nnz, b->csx, b->coo_col,
                                          *descr_c,         c->csx, &c->nnz     );

            CUDA errorcheck( cdata_, (void*)&status, 1, "first pass" );

        }

        // allocate memory for C's tjv and tkv if already haven't
        if( !c->coo_col )
            CUDA ps_cuda_malloc( cdata_, &c->coo_col, c->nnz, dev );

        if( !c->coo_val )
            CUDA ps_cuda_malloc( cdata_, &c->coo_val, c->nnz, dev );

        // perform multiplication

#ifdef SINGLE_PRECISION
        status = cusparseScsrgeam( *handle, m, n,
                                   &alpha, *descr_a, a->nnz,
                                   a->coo_val, a->csx, a->coo_col,
                                   &beta,  *descr_b, b->nnz,
                                   b->coo_val, b->csx, b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col  );
#else
        status = cusparseDcsrgeam( *handle, m, n,
                                   &alpha, *descr_a, a->nnz,
                                   a->coo_val, a->csx, a->coo_col,
                                   &beta,  *descr_b, b->nnz,
                                   b->coo_val, b->csx, b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col  );
#endif

        CUDA errorcheck( cdata_, (void*)&status, 1, "matrix-matrix multiplication" );

    }
    else if( tr1 && !tr2 )
    {
        if( !(c->csx) )
        {
            // allocate memory for C's CSR row pointer on GPU
            // and figure out number of nonzero entries
            CUDA ps_cuda_malloc( cdata_, &c->csx, m+1, dev );

            status = cusparseXcsrgeamNnz( *handle, m, n,
                                          *descr_a, a->nnz, a->csx_tr, a->coo_row_tr,
                                          *descr_b, b->nnz, b->csx,    b->coo_col,
                                          *descr_c,         c->csx,   &c->nnz         );

            CUDA errorcheck( cdata_, (void*)&status, 1, "first pass" );

        }

        // allocate memory for C's tjv and tkv if already haven't
        if( !c->coo_col )
            CUDA ps_cuda_malloc( cdata_, &c->coo_col, c->nnz, dev );

        if( !c->coo_val )
            CUDA ps_cuda_malloc( cdata_, &c->coo_val, c->nnz, dev );

        // perform multiplication

#ifdef SINGLE_PRECISION
        status = cusparseScsrgeam( *handle, m, n,
                                   &alpha, *descr_a, a->nnz,
                                   a->coo_val_tr, a->csx_tr, a->coo_row_tr,
                                   &beta,  *descr_b, b->nnz,
                                   b->coo_val,    b->csx,    b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col           );
#else
        status = cusparseDcsrgeam( *handle, m, n,
                                   &alpha, *descr_a, a->nnz,
                                   a->coo_val_tr, a->csx_tr, a->coo_row_tr,
                                   &beta,  *descr_b, b->nnz,
                                   b->coo_val,    b->csx,    b->coo_col,
                                   *descr_c,
                                   c->coo_val, c->csx, c->coo_col           );
#endif

        CUDA errorcheck( cdata_, (void*)&status, 1, "matrix-matrix multiplication" );

    }

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cusparse_lock );

    CUDA postprocess (a, b, c);

}

PS_FPT pixsrc_cuda::ps_cuda_max( VECTOR *a )
{
    // get maximum value of vector

    // update GPU data
    a->update_gpu();

    // set GPU device and handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;
    cublasHandle_t *handle_cublas  = (cublasHandle_t*)cdata_->gpudevices[dev].handle_cublas;
    cublasStatus_t stat;

    // set mutex lock
    pthread_mutex_lock( cdata_->gpudevices[dev].cublas_lock );

    // use a's stream
    cudaStream_t *stream = (cudaStream_t*)a->stream;
    cublasSetStream( *handle_cublas, *stream );

    // find index of maximum of vector
    int ind;
    PS_FPT res;
    if( a->scalar > 0 )
    {

#ifdef SINGLE_PRECISION
        stat = cublasIsamax( *handle_cublas, a->size, a->cuda_vec, 1, &ind );
#else
        stat = cublasIdamax( *handle_cublas, a->size, a->cuda_vec, 1, &ind );
#endif

    }
    else
    {

#ifdef SINGLE_PRECISION
        stat = cublasIsamin( *handle_cublas, a->size, a->cuda_vec, 1, &ind );
#else
        stat = cublasIdamin( *handle_cublas, a->size, a->cuda_vec, 1, &ind );
#endif

    }

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cublas_lock );

    CUDA errorcheck( cdata_, (void*)&stat, 0, "vec max" );

    // copy maximum number into res
    CUDA ps_cuda_memcpy( cdata_, &res, a->cuda_vec + ind, 1, 1, stream, dev );

    // not sure why this is being called here, but I think it is to make the code
    // wait until the GPU is finished with all calculations
    a->update_gpu();

    return a->scalar * res;
}

void pixsrc_cuda::tri_eq_solve( MATRIX *a, VECTOR *b, VECTOR *c, bool zeroit )
{
    // solve linear equation

    // Note: I don't think I'm using this anymore.
    // I can't remember why I decided not to use cusp.


    zeroit = 1;

    /*
      b->set_scalar( 1.0 );
      CUDA fill_vec( b, b->max() );
    */

    if( zeroit )
        b->zeromeout();
    else
        b->update_gpu();

    a->update_gpu( pmf_csr );
    c->update_gpu();

    a->dissolve_scalar();
    c->dissolve_scalar();


    /*
      void *clock;
      float time;
      CUDA ps_start_clock( &clock );
    */

    // get cusp preconditoner
    if( !a->cusp_prec )
    {
        psc_get_preconditioner_cusp( a->csx, a->coo_col, a->coo_val,
                                     a->nrow, a->nnz, &a->cusp_prec  );
    }

    a->update_gpu( pmf_csr );

    psc_tri_eq_solve_cusp( a->csx, a->coo_col, a->coo_val,
                           b->cuda_vec, c->cuda_vec, a->nrow, a->nnz, a->cusp_prec );

    /*
      CUDA ps_stop_clock( clock, &time );
      std::cout << "cusp took " << time << " ms" << std::endl;
    */

    cudaDeviceSynchronize();

    b->set_status( pvs_gpu );

    /*
      CUDA ps_wrap_cusp( a );
      CUDA ps_wrap_cusp( b );
      CUDA ps_wrap_cusp( c );

      typedef typename cusp::array1d_view< thrust::device_ptr<PS_SIT>    > DeviceIndexArrayView;
      typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

      typedef cusp::csr_matrix_view<DeviceIndexArrayView,
      DeviceIndexArrayView,
      DeviceValueArrayView> DeviceView;

      DeviceValueArrayView *res = (DeviceValueArrayView*)b->cusp_wrapper;
      DeviceValueArrayView *rhs = (DeviceValueArrayView*)c->cusp_wrapper;
      DeviceView           *mat = (DeviceView*)          a->cusp_csr;
    */
}
/*
  void resolve_dep( MATRIX *a, int pos, vector<PS_SIT> *dep )
  {
  a->update_cpu( pmf_csr );
  for( int g=a->apr[pos]; g<a->apr[pos+1]; ++g )
  {
  dep->push_back(a->air[g]);
  }

  for( unsigned int j=0; j<dep->size(); ++j )
  {
  int r = (*dep)[j];
  for( int g=a->apr[r]; g<a->apr[r+1]; ++g )
  {
  bool addit = 1;
  for( unsigned int jj=0; jj<dep->size(); ++jj )
  {
  if(a->air[g] == (*dep)[jj])
  {
  addit = 0;
  break;
  }
  }
  if( addit )
  {
  dep->push_back(a->air[g]);
  }
  }
  }
  }
*/

PS_FPT *bleh = 0;
void pixsrc_cuda::tr_inv_a_b( MATRIX *a, MATRIX *b, PS_FPT *tr )
{
    // this function calculates the trace of the following:
    // tr[ a^(-1)*b  ]
    // this is the cycle function

    // Note: I don't use this anymore.
    // I can't remember why I stopped.

    a->update_gpu( pmf_csr );
    b->update_gpu( pmf_csc );

    commoninputdata  *cdata_ = a->cdata_;
    int dev = a->data_->gpu2use[0];

    int l_side = a->nrow;



/*
  vector< vector<PS_SIT> > dep(l_side);
  for(int c=0; c<l_side; ++c)
  {
  resolve_dep( a, c, &dep[c] );
  }

  for( int j=0; j<l_side; ++j )
  std::cout << dep[j].size() << " ";
  std::cout << std::endl;
  exit(1);
*/
    int l_side2 = l_side * l_side;

    PS_FPT *res, *rhs;
    if( !bleh )
    {
        CUDA ps_cuda_malloc( cdata_, &res, l_side2, dev );
        bleh = res;
    }
    else
    {
        res = bleh;
    }
    CUDA ps_cuda_malloc( cdata_, &rhs, l_side2, dev );

    if( !a->cusp_prec )
    {
        psc_get_preconditioner_cusp( a->csx, a->coo_col, a->coo_val,
                                     a->nrow, a->nnz, &a->cusp_prec  );
    }

    psc_tr_inv_a_b_matrixfree( a->csx, a->coo_col, a->coo_val,
                               b->csx, b->coo_col, b->coo_val,
                               res, rhs, l_side,
                               cdata_->gpudevices[dev].handle,
                               cdata_->gpudevices[dev].descr, a->nnz, a->cusp_prec );

    CUDA ps_cuda_memcpy( cdata_, tr, res, 1, 1, NULL, dev );
    cudaDeviceSynchronize();

//    CUDA ps_cuda_free( res );
    CUDA ps_cuda_free( rhs );










/*
  int l_side = a->nrow;
  int l_side2 = l_side * l_side;

  PS_FPT *res, *rhs;
  CUDA ps_cuda_malloc( cdata_, &res, l_side2 );
  CUDA ps_cuda_malloc( cdata_, &rhs, l_side2 );

  psc_tr_inv_a_b_matrixfree( a->csx, a->coo_col, a->coo_val,
  b->csx, b->coo_col, b->coo_val,
  res, rhs, l_side,
  cdata_->gpudevices[dev].handle,
  cdata_->gpudevices[dev].descr, a->nnz );

  CUDA ps_cuda_memcpy( cdata_, tr, res, 1, 1 );

  CUDA ps_cuda_free( res );
  CUDA ps_cuda_free( rhs );
*/
}

///////////////////////////
// CURRENTLY UNUSED CODE //
///////////////////////////
// it's not all useless, fyi


// code that uses self-written conjugate gradient solver
// it currently requires CHOLMOD support, which requires
// double precision
// this could be circumvented using CUSPARSE's incomplete
// Cholesky factorization, but that doesn't seem to be
// working correctly.
// I use CUSP instead. Plus, CUSP's linear solver allows
// matrix-free methods which is great for calculating
// tr[A^-1 C]
/*
  void pixsrc_cuda::tri_eq_solve( MATRIX *a, VECTOR *b, VECTOR *c )
  {
  b->set_scalar( 1.0 );

  // transposing MATRIX a is not supported yet
  // c is right-hand-side vector
  // b is solution vector

  // few things we'll need
  int dev                           = a->data_->gpu2use[0];
  CUDA set_device( dev );
  commoninputdata  *cdata_          = a->cdata_;
  cusparseHandle_t *handle          = (cusparseHandle_t*)  cdata_->gpudevices[dev].handle;
  cublasHandle_t *handle_cublas     = (cublasHandle_t*)cdata_->gpudevices[dev].handle_cublas;
  cusparseMatDescr_t *descr_gen     = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;
  cusparseOperation_t trans         = CUSPARSE_OPERATION_TRANSPOSE;
  cusparseOperation_t nontrans      = CUSPARSE_OPERATION_NON_TRANSPOSE;
  //cusparseSolveAnalysisInfo_t *info = (cusparseSolveAnalysisInfo_t*)a->csrsvi;
  cusparseStatus_t status;
  cusparseMatDescr_t *descrT = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr_tri_upper;

  int l_side = a->nrow;

  // get working space on gpu
  for( int j=0; j<5; ++j )
  if( !a->ctv[j] )
  CUDA ps_cuda_malloc( cdata_, &a->ctv[j], l_side );

  PS_FPT *temp_vec_r, *temp_vec_t, *temp_vec_z, *temp_vec_q, *temp_vec_p;
  temp_vec_r = a->ctv[0];
  temp_vec_t = a->ctv[1];
  temp_vec_z = a->ctv[2];
  temp_vec_q = a->ctv[3];
  temp_vec_p = a->ctv[4];

  // make sure matrix a has CSR format ready and
  // vectors have been prepped
  a->update_gpu( pmf_csr );
  c->update_gpu();

  //CUDA sync_vector( cdata_, b, 1 );
  //CUDA sync_vector( cdata_, c, 1 );

  // perform incomplete Cholesky factorization
  CUDA inc_chol_factor( a );

  CUDA analyze_upp_low( a );
  cusparseSolveAnalysisInfo_t *infoU = (cusparseSolveAnalysisInfo_t*)a->icf_upper_sva;
  cusparseSolveAnalysisInfo_t *infoL = (cusparseSolveAnalysisInfo_t*)a->icf_lower_sva;

  // allocate memory for solution vector and temporar vectors
  int max_iter = 2000;
  PS_FPT alpha, beta, nrmr0, nrmr, rho, rhop, temp_num;
  PS_FPT conv_tol = 1e-6;
  PS_FPT dummy_zero  =  0.0;
  PS_FPT dummy_m_one = -1.0;
  PS_FPT dummy_p_one =  1.0;








  // compute initial guess

  // solve equation with lower triangle
  alpha = c->scalar / a->scalar;

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, c->cuda_vec, temp_vec_r          );
  #else
  status = cusparseDcsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, c->cuda_vec, temp_vec_r          );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system L" );

  // solve equation with upper triangle

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_r, b->cuda_vec            );
  #else
  status = cusparseDcsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_r, b->cuda_vec            );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system U" );

  // compute initial residual from initial guess

  #ifdef SINGLE_PRECISION
  cusparseScsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  b->cuda_vec, &dummy_zero, temp_vec_r                    );

  cublasSscal( *handle_cublas, l_side, &dummy_m_one, temp_vec_r,  1                );
  cublasSaxpy( *handle_cublas, l_side,  &c->scalar,  c->cuda_vec, 1, temp_vec_r, 1 );

  cublasSnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr0 );
  #else
  cusparseDcsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  b->cuda_vec, &dummy_zero, temp_vec_r                    );

  cublasDscal( *handle_cublas, l_side, &dummy_m_one, temp_vec_r,  1                );
  cublasDaxpy( *handle_cublas, l_side,  &c->scalar,  c->cuda_vec, 1, temp_vec_r, 1 );

  cublasDnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr0 );
  #endif


  nrmr0 *= nrmr0;







  // start iterative solving

  int evals = 0;
  for( int iter=0; iter<max_iter; ++iter )
  {
  ++evals;
  // solve linear systems

  // solve equation with lower triangle
  alpha = 1.0 / a->scalar;

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, temp_vec_r, temp_vec_t            );
  #else
  status = cusparseDcsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, temp_vec_r, temp_vec_t            );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system L" );

  // solve equation with upper triangle

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_t, temp_vec_z            );
  #else
  status = cusparseDcsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_t, temp_vec_z            );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system U" );


  // compute weights
  rhop = rho;

  #ifdef SINGLE_PRECISION
  cublasSdot( *handle_cublas, l_side, temp_vec_r, 1, temp_vec_z, 1, &rho );
  #else
  cublasDdot( *handle_cublas, l_side, temp_vec_r, 1, temp_vec_z, 1, &rho );
  #endif

  // get next correction
  if( !iter )
  {

  #ifdef SINGLE_PRECISION
  cublasScopy( *handle_cublas, l_side, temp_vec_z, 1, temp_vec_p, 1 );
  #else
  cublasDcopy( *handle_cublas, l_side, temp_vec_z, 1, temp_vec_p, 1 );
  #endif

  }
  else
  {
  beta = rho / rhop;

  #ifdef SINGLE_PRECISION
  cublasSaxpy( *handle_cublas, l_side, &beta, temp_vec_p, 1, temp_vec_z, 1 );
  cublasScopy( *handle_cublas, l_side,        temp_vec_z, 1, temp_vec_p, 1 );
  #else
  cublasDaxpy( *handle_cublas, l_side, &beta, temp_vec_p, 1, temp_vec_z, 1 );
  cublasDcopy( *handle_cublas, l_side,        temp_vec_z, 1, temp_vec_p, 1 );
  #endif

  }

  #ifdef SINGLE_PRECISION
  cusparseScsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  temp_vec_p, &dummy_zero, temp_vec_q                  );
  #else
  cusparseDcsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  temp_vec_p, &dummy_zero, temp_vec_q                  );
  #endif

  #ifdef SINGLE_PRECISION
  // get next best solutions
  cublasSdot( *handle_cublas, l_side, temp_vec_p, 1, temp_vec_q, 1, &temp_num );
  alpha = rho / temp_num;
  cublasSaxpy( *handle_cublas, l_side, &alpha, temp_vec_p, 1, b->cuda_vec, 1 );
  alpha *= -1.0;
  cublasSaxpy( *handle_cublas, l_side, &alpha, temp_vec_q, 1, temp_vec_r,  1 );

  // check for convergence
  cublasSnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr );

  #else
  // get next best solutions
  cublasDdot( *handle_cublas, l_side, temp_vec_p, 1, temp_vec_q, 1, &temp_num );
  alpha = rho / temp_num;
  cublasDaxpy( *handle_cublas, l_side, &alpha, temp_vec_p, 1, b->cuda_vec, 1 );
  alpha *= -1.0;
  cublasDaxpy( *handle_cublas, l_side, &alpha, temp_vec_q, 1, temp_vec_r,  1 );

  // check for convergence
  cublasDnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr );
  #endif

  if( nrmr*nrmr / nrmr0 < conv_tol )
  {
  break;
  }
  }

  b->set_status( pvs_gpu );
  }
*/


// code that uses self-written conjugate gradient solver
// it currently requires CHOLMOD support, which requires
// double precision
// this could be circumvented using CUSPARSE's incomplete
// Cholesky factorization, but that doesn't seem to be
// working correctly.
// I use CUSP instead. Plus, CUSP's linear solver allows
// matrix-free methods which is great for calculating
// tr[A^-1 C]
/*
  typedef struct
  {
  PS_FPT **tempvecs;
  int lower, upper;
  VECTOR *res, *rhs;
  PS_FPT *traces;
  MATRIX *a, *b;
  cudaStream_t *streamID;

  } tr_invA_B_kernel_struct;


  void *pixsrc_cuda::tr_invA_B_kernel( void *args)
  {
  tr_invA_B_kernel_struct     *st        = (tr_invA_B_kernel_struct*)args;
  PS_FPT                      **tempvecs = st->tempvecs;
  int                          lower     = st->lower;
  int                          upper     = st->upper;
  PS_FPT                      *res       = st->res->cuda_vec;
  PS_FPT                      *rhs       = st->rhs->cuda_vec;
  PS_FPT                      *traces    = st->traces;
  MATRIX                      *a         = st->a;
  MATRIX                      *c         = st->b;
  cudaStream_t                *streamID  = st->streamID;

  int                          dev       = a->data_->gpu2use[0];
  commoninputdata             *cdata_    = a->cdata_;
  cusparseMatDescr_t          *descr_gen   = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr;
  cusparseMatDescr_t *descrT = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr_tri_upper;
  cusparseSolveAnalysisInfo_t *infoU = (cusparseSolveAnalysisInfo_t*)a->icf_upper_sva;
  cusparseSolveAnalysisInfo_t *infoL = (cusparseSolveAnalysisInfo_t*)a->icf_lower_sva;

  cusparseOperation_t          trans     = CUSPARSE_OPERATION_TRANSPOSE;
  cusparseOperation_t          nontrans  = CUSPARSE_OPERATION_NON_TRANSPOSE;

  cusparseStatus_t status;
  cublasStatus_t stat;

  // get handle
  cusparseHandle_t *handle;
  MEMORY ps_malloc( &handle, 1 );
  status = cusparseCreate( handle );
  CUDA errorcheck( cdata_, (void*)&status, 1, "create handle"     );

  // get cublas handle
  cublasHandle_t *handle_cublas;
  MEMORY ps_malloc( &handle_cublas, 1 );
  stat = cublasCreate( handle_cublas );
  CUDA errorcheck( cdata_, (void*)&stat, 2, "cublas handle" );

  CUDA set_device( dev );

  status = cusparseSetStream( *handle,        *streamID );
  CUDA errorcheck( cdata_, (void*)&status, 1, "cusparse stream set" );

  stat = cublasSetStream  ( *handle_cublas, *streamID );
  CUDA errorcheck( cdata_, (void*)&status, 2, "cublas stream set" );

  PS_FPT *temp_vec_r = tempvecs[0];
  PS_FPT *temp_vec_t = tempvecs[1];
  PS_FPT *temp_vec_z = tempvecs[2];
  PS_FPT *temp_vec_q = tempvecs[3];
  PS_FPT *temp_vec_p = tempvecs[4];

  int max_iter = 2000;
  int l_side = a->nrow;
  PS_FPT conv_tol = 1e-6;
  PS_FPT dummy_zero  =  0.0;
  PS_FPT dummy_m_one = -1.0;
  PS_FPT dummy_p_one =  1.0;
  PS_FPT alpha, beta, nrmr0, nrmr, rho, rhop, temp_num;

  for( int col=lower; col<upper; ++col )
  {
  psc_cpy_mat_col_to_vec ( col, c->csx_tr, c->coo_row_tr, c->coo_val_tr, rhs );


  // compute initial guess

  // solve equation with lower triangle
  alpha = c->scalar / a->scalar;

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, rhs, temp_vec_r             );
  #else
  status = cusparseDcsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, rhs, temp_vec_r             );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system L" );

  // solve equation with upper triangle

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_r, res                  );
  #else
  status = cusparseDcsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_r, res                  );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system U" );

  // compute initial residual from initial guess

  #ifdef SINGLE_PRECISION
  cusparseScsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  res, &dummy_zero, temp_vec_r                            );


  cublasSscal( *handle_cublas, l_side, &dummy_m_one, temp_vec_r,  1                );
  cublasSaxpy( *handle_cublas, l_side,  &c->scalar,  rhs        , 1, temp_vec_r, 1 );

  cublasSnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr0 );
  #else
  cusparseDcsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  res, &dummy_zero, temp_vec_r                            );


  cublasDscal( *handle_cublas, l_side, &dummy_m_one, temp_vec_r,  1                );
  cublasDaxpy( *handle_cublas, l_side,  &c->scalar,  rhs        , 1, temp_vec_r, 1 );

  cublasDnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr0 );
  #endif

  nrmr0 *= nrmr0;

  for( int iter=0; iter<max_iter; ++iter )
  {
  // solve linear systems

  // solve equation with lower triangle
  alpha = 1.0 / a->scalar;

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, temp_vec_r, temp_vec_t            );
  #else
  status = cusparseDcsrsv_solve( *handle, trans,
  l_side, &alpha, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoL, temp_vec_r, temp_vec_t            );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system L" );

  // solve equation with upper triangle

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_t, temp_vec_z            );
  #else
  status = cusparseDcsrsv_solve( *handle, nontrans,
  l_side, &dummy_p_one, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *infoU, temp_vec_t, temp_vec_z            );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "solve triangular system U" );


  // compute weights
  rhop = rho;

  #ifdef SINGLE_PRECISION
  cublasSdot( *handle_cublas, l_side, temp_vec_r, 1, temp_vec_z, 1, &rho );
  #else
  cublasDdot( *handle_cublas, l_side, temp_vec_r, 1, temp_vec_z, 1, &rho );
  #endif

  // get next correction
  if( !iter )
  {
  #ifdef SINGLE_PRECISION
  cublasScopy( *handle_cublas, l_side, temp_vec_z, 1, temp_vec_p, 1 );
  #else
  cublasDcopy( *handle_cublas, l_side, temp_vec_z, 1, temp_vec_p, 1 );
  #endif
  }
  else
  {
  beta = rho / rhop;

  #ifdef SINGLE_PRECISION
  cublasSaxpy( *handle_cublas, l_side, &beta, temp_vec_p, 1, temp_vec_z, 1 );
  cublasScopy( *handle_cublas, l_side,        temp_vec_z, 1, temp_vec_p, 1 );
  #else
  cublasDaxpy( *handle_cublas, l_side, &beta, temp_vec_p, 1, temp_vec_z, 1 );
  cublasDcopy( *handle_cublas, l_side,        temp_vec_z, 1, temp_vec_p, 1 );
  #endif
  }

  #ifdef SINGLE_PRECISION
  cusparseScsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  temp_vec_p, &dummy_zero, temp_vec_q                  );
  #else
  cusparseDcsrmv( *handle, nontrans,
  l_side, l_side, a->nnz, &a->scalar, *descr_gen,
  a->coo_val, a->csx, a->coo_col,
  temp_vec_p, &dummy_zero, temp_vec_q                  );
  #endif

  #ifdef SINGLE_PRECISION
  // get next best solutions

  cublasSdot( *handle_cublas, l_side, temp_vec_p, 1, temp_vec_q, 1, &temp_num );
  alpha = rho / temp_num;
  cublasSaxpy( *handle_cublas, l_side, &alpha, temp_vec_p, 1, res,          1 );
  alpha *= -1.0;
  cublasSaxpy( *handle_cublas, l_side, &alpha, temp_vec_q, 1, temp_vec_r,   1 );

  // check for convergence

  cublasSnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr );
  #else
  // get next best solutions
  cublasDdot( *handle_cublas, l_side, temp_vec_p, 1, temp_vec_q, 1, &temp_num );
  alpha = rho / temp_num;
  cublasDaxpy( *handle_cublas, l_side, &alpha, temp_vec_p, 1, res,          1 );
  alpha *= -1.0;
  cublasDaxpy( *handle_cublas, l_side, &alpha, temp_vec_q, 1, temp_vec_r,   1 );

  // check for convergence
  cublasDnrm2( *handle_cublas, l_side, temp_vec_r, 1, &nrmr );
  #endif

  if( nrmr*nrmr / nrmr0 < conv_tol )
  {
  break;
  }
  }

  CUDA ps_cuda_memcpy( cdata_, traces+col, res+col, 1, 2 );

  psc_undo_cpy_mat_col_to_vec ( col, c->csx_tr, c->coo_row_tr, rhs );
  }

  cudaDeviceSynchronize();

  cublasDestroy( *handle_cublas );
  cusparseDestroy( *handle );
  MEMORY ps_free( handle_cublas );
  MEMORY ps_free( handle );

  return NULL;
  }

  void pixsrc_cuda::cusparse_tr_inva_b( MATRIX *a, MATRIX *b, PS_FPT *tr )
  {
  // this function calculates the trace of the following:
  // tr[ a^(-1)*b  ]
  // this is the cycle function

  a->update_gpu( pmf_csr );
  b->update_gpu( pmf_csc );

  // perform incomplete Cholesky factorization
  CUDA inc_chol_factor( a );

  CUDA analyze_upp_low( a );

  // few things we'll need
  int dev                    = a->data_->gpu2use[0];
  CUDA set_device( dev );
  commoninputdata  *cdata_   = a->cdata_;
  inputdata *data_           = a->data_;
  int num_kernels = 1;//cdata_->gpudevices[dev].num_kernels;

  cublasHandle_t *handle_cublas = (cublasHandle_t*)(cdata_->gpudevices[dev].handle_cublas);
  cusparseHandle_t *handle = (cusparseHandle_t*)(cdata_->gpudevices[dev].handle);

  // error recorders
  //cusparseStatus_t status;
  cudaError_t      cudaStat;

  int l_side = a->nrow;
  PS_FPT alpha = b->scalar / a->scalar;
  int interval = b->ncol / num_kernels;

  // 1st vector to hold the right-hand side
  //  rhs has an extra element thats set to zero
  // 2nd to hold the result
  // 3rd to hold diagonal elements of res
  // 4th is a dummy filled with zeros
  VECTOR **rhs_, **res_, **rhs_dummy, **res_dummy;
  PS_FPT *inv_a_b_diag;

  // allocate and fill dummy
  // diagonal elements of trace will be placed in first nrow elements

  // allocate gpu vector memory
  MEMORY ps_malloc( &rhs_dummy, num_kernels, 1 );
  MEMORY ps_malloc( &res_dummy, num_kernels, 1 );
  MEMORY ps_malloc( &rhs_,      num_kernels    );
  MEMORY ps_malloc( &res_,      num_kernels    );

  for( int k=0; k<num_kernels; ++k )
  {
  rhs_[k] = new ( rhs_dummy[k] ) VECTOR( cdata_, data_, l_side );
  res_[k] = new ( res_dummy[k] ) VECTOR( cdata_, data_, l_side );
  }

  CUDA ps_cuda_malloc( cdata_, &inv_a_b_diag, l_side );

  // copy zeros to vectors
  for( int k=0; k<num_kernels; ++k )
  {
  psc_fill( res_[k]->cuda_vec, 0.0, l_side );
  psc_fill( rhs_[k]->cuda_vec, 0.0, l_side );
  }


  // initializing streams to launch simultaneous kernels
  cudaStream_t *streamID;
  MEMORY ps_malloc( &streamID, num_kernels );
  for( int k=0; k<num_kernels; ++k )
  cudaStreamCreate( &streamID[k] );

  // getting temporary vecs

  PS_FPT **tempvecs;
  MEMORY ps_malloc( &tempvecs, num_kernels*5 );

  for( int j=0; j<5; ++j )
  {
  if( !a->ctv[j] )
  CUDA ps_cuda_malloc( cdata_, &a->ctv[j], l_side );
  tempvecs[j] = a->ctv[j];
  }

  for(int k=1; k<num_kernels; ++k)
  {
  for( int p=0; p<5; ++p )
  CUDA ps_cuda_malloc( cdata_, &tempvecs[k*5+p], l_side );
  }


  // setup structs for pthreads
  tr_invA_B_kernel_struct *sts;
  MEMORY ps_malloc( &sts, num_kernels );
  for( int k=0; k<num_kernels; ++k )
  {
  sts[k].tempvecs  = &tempvecs[k*5];
  sts[k].res       = res_[k];
  sts[k].rhs       = rhs_[k];
  sts[k].traces    = inv_a_b_diag;
  sts[k].a         = a;
  sts[k].b         = b;
  sts[k].streamID  = &streamID[k];
  sts[k].lower     = k   *interval;
  sts[k].upper     = (k+1)*interval;

  if( k == num_kernels-1 )
  sts[k].upper = b->ncol;
  }

  // create pthreads
  pthread_t *threads;
  MEMORY ps_malloc( &threads, num_kernels);

  // launch pthreads
  for( int k=0; k<num_kernels; ++k )
  {
  pthread_create( &threads[k], cdata_->attrjoinable,
  CUDA tr_invA_B_kernel, &sts[k] );
  }


  // wait for threads to finish
  for( int k=0; k<num_kernels; ++k )
  pthread_join( threads[k], NULL );

  *tr = psc_sum_vector ( inv_a_b_diag, l_side );


  // memory cleanup

  for(int k=1; k<num_kernels; ++k)
  {
  for( int p=0; p<5; ++p )
  CUDA ps_cuda_free( tempvecs[k*5+p] );
  }
  MEMORY ps_free( tempvecs );

  for( int k=0; k<num_kernels; ++k )
  {
  res_[k]->~VECTOR();
  rhs_[k]->~VECTOR();
  }
  MEMORY ps_free( rhs_dummy, num_kernels );
  MEMORY ps_free( res_dummy, num_kernels );
  MEMORY ps_free( rhs_ );
  MEMORY ps_free( res_ );

  CUDA   ps_cuda_free( inv_a_b_diag );

  for( int k=0; k<num_kernels; ++k )
  cudaStreamDestroy( streamID[k] );
  MEMORY ps_free( streamID );

  MEMORY ps_free( threads );
  MEMORY ps_free( sts     );
  }
*/

/*
  typedef struct
  {
  int lower, upper;
  VECTOR *res, *rhs;
  PS_FPT *traces;
  MATRIX *a, *b;
  cudaStream_t *streamID;
  void *prec;

  } tr_invA_B_kernel_struct_cusp;

  void *pixsrc_cuda::tr_invA_B_kernel_cusp( void *args)
  {
  tr_invA_B_kernel_struct_cusp     *st        = (tr_invA_B_kernel_struct_cusp*)args;
  int                          lower     = st->lower;
  int                          upper     = st->upper;
  PS_FPT                      *res       = st->res->cuda_vec;
  PS_FPT                      *rhs       = st->rhs->cuda_vec;
  PS_FPT                      *traces    = st->traces;
  MATRIX                      *a         = st->a;
  MATRIX                      *c         = st->b;
  cudaStream_t                *streamID  = st->streamID;
  void                        *prec      = st->prec;

  int                          dev       = a->data_->gpu2use[0];
  commoninputdata             *cdata_    = a->cdata_;

  CUDA set_device( dev );

  for( int col=lower; col<upper; ++col )
  {
  std::cout << col << std::endl;

  psc_cpy_mat_col_to_vec ( col, c->csx_tr, c->coo_row_tr, c->coo_val_tr, rhs );

  psc_tri_eq_solve_cusp( a->csx, a->coo_col, a->coo_val,
  res, rhs, a->nrow, a->nnz, prec );

  /-*
  psc_tr_inv_a_b_cusp( a->csx, a->coo_col, a->coo_val,
  res, rhs, a->nrow, a->nnz, streamID );
  *-/

  CUDA ps_cuda_memcpy( cdata_, traces+col, res+col, 1, 2 );

  psc_undo_cpy_mat_col_to_vec ( col, c->csx_tr, c->coo_row_tr, rhs );
  }

  cudaDeviceSynchronize();

  return NULL;
  }

  void pixsrc_cuda::cusp_tr_inva_b_cusp( MATRIX *a, MATRIX *b, PS_FPT *tr )
  {
  // this function calculates the trace of the following:
  // tr[ a^(-1)*b  ]
  // this is the cycle function

  a->update_gpu( pmf_csr );
  b->update_gpu( pmf_csc );

  // few things we'll need
  int dev                    = a->data_->gpu2use[0];
  CUDA set_device( dev );
  commoninputdata  *cdata_   = a->cdata_;
  inputdata *data_           = a->data_;
  int num_kernels = 1;//cdata_->gpudevices[dev].num_kernels;

  cublasHandle_t *handle_cublas = (cublasHandle_t*)(cdata_->gpudevices[dev].handle_cublas);
  cusparseHandle_t *handle = (cusparseHandle_t*)(cdata_->gpudevices[dev].handle);

  // error recorders
  //cusparseStatus_t status;
  cudaError_t      cudaStat;

  int l_side = a->nrow;
  PS_FPT alpha = b->scalar / a->scalar;
  int interval = b->ncol / num_kernels;

  // 1st vector to hold the right-hand side
  //  rhs has an extra element thats set to zero
  // 2nd to hold the result
  // 3rd to hold diagonal elements of res
  // 4th is a dummy filled with zeros
  VECTOR **rhs_, **res_, **rhs_dummy, **res_dummy;
  PS_FPT *inv_a_b_diag;

  // allocate and fill dummy
  // diagonal elements of trace will be placed in first nrow elements

  // allocate gpu vector memory
  MEMORY ps_malloc( &rhs_dummy, num_kernels, 1 );
  MEMORY ps_malloc( &res_dummy, num_kernels, 1 );
  MEMORY ps_malloc( &rhs_,      num_kernels    );
  MEMORY ps_malloc( &res_,      num_kernels    );

  for( int k=0; k<num_kernels; ++k )
  {
  rhs_[k] = new ( rhs_dummy[k] ) VECTOR( cdata_, data_, l_side );
  res_[k] = new ( res_dummy[k] ) VECTOR( cdata_, data_, l_side );
  }

  CUDA ps_cuda_malloc( cdata_, &inv_a_b_diag, l_side );

  // copy zeros to vectors
  for( int k=0; k<num_kernels; ++k )
  {
  psc_fill( res_[k]->cuda_vec, 0.0, l_side );
  psc_fill( rhs_[k]->cuda_vec, 0.0, l_side );
  }


  // initializing streams to launch simultaneous kernels
  cudaStream_t *streamID;
  MEMORY ps_malloc( &streamID, num_kernels );
  for( int k=0; k<num_kernels; ++k )
  cudaStreamCreate( &streamID[k] );

  // get preconditioner
  void *prec;
  psc_get_preconditioner_cusp(a->csx, a->coo_col, a->coo_val,
  a->nrow, a->nnz, &prec );

  // setup structs for pthreads
  tr_invA_B_kernel_struct_cusp *sts;
  MEMORY ps_malloc( &sts, num_kernels );
  for( int k=0; k<num_kernels; ++k )
  {
  sts[k].res       = res_[k];
  sts[k].rhs       = rhs_[k];
  sts[k].traces    = inv_a_b_diag;
  sts[k].a         = a;
  sts[k].b         = b;
  sts[k].streamID  = &streamID[k];
  sts[k].lower     = k   *interval;
  sts[k].upper     = (k+1)*interval;
  sts[k].prec      = prec;

  if( k == num_kernels-1 )
  sts[k].upper = b->ncol;
  }

  // create pthreads
  pthread_t *threads;
  MEMORY ps_malloc( &threads, num_kernels);

  // launch pthreads
  for( int k=0; k<num_kernels; ++k )
  {
  pthread_create( &threads[k], cdata_->attrjoinable,
  CUDA tr_invA_B_kernel_cusp, &sts[k] );
  }


  // wait for threads to finish
  for( int k=0; k<num_kernels; ++k )
  pthread_join( threads[k], NULL );

  *tr = psc_sum_vector ( inv_a_b_diag, l_side );


  // memory cleanup

  for( int k=0; k<num_kernels; ++k )
  {
  res_[k]->~VECTOR();
  rhs_[k]->~VECTOR();
  }
  MEMORY ps_free( rhs_dummy, num_kernels );
  MEMORY ps_free( res_dummy, num_kernels );
  MEMORY ps_free( rhs_ );
  MEMORY ps_free( res_ );

  CUDA   ps_cuda_free( inv_a_b_diag );

  for( int k=0; k<num_kernels; ++k )
  cudaStreamDestroy( streamID[k] );
  MEMORY ps_free( streamID );

  MEMORY ps_free( threads );
  MEMORY ps_free( sts     );
  }
*/
