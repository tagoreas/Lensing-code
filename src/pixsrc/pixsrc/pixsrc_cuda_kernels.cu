// -*- C++ -*-
// the previous line is for emacs .. so that it formats
// this file according to the c++ format style
// specified in ~/.emacs

// this is CUDA code that needs to be compiled by NVCC
// the resulting ibject file can be linked with g++
// I think most of this code is unused currently.

#ifndef PS_FPT
#ifdef SINGLE_PRECISION
#define PS_FPT float
#endif
#ifdef DOUBLE_PRECISION
#define PS_FPT double
#endif
#endif

#include <cusparse_v2.h>
#include <cuda_runtime.h>

#include <thrust/device_ptr.h>
#include <thrust/fill.h>

#include <cusp/csr_matrix.h>
#include <cusp/krylov/cg.h>
#include <cusp/precond/ainv.h>
#include <cusp/precond/aggregation/smoothed_aggregation.h>
#include <cusp/krylov/gmres.h>
#include <cusp/print.h>


// fill vector with constant value
__global__ void psc_fill_kernel( PS_FPT *vec, PS_FPT val )
{
    int i = blockIdx.x;

    vec[i] = val;

    return;
}

__global__ void psc_cpy_mat_col_to_vec_kernel( int col,
                                               int *apc, int *aic, PS_FPT *axc, PS_FPT *vec )
{
    for( int g=apc[col]; g<apc[col+1]; ++g )
        vec[aic[g]] = axc[g];

    return;
}

__global__ void psc_undo_cpy_mat_col_to_vec_kernel( int col,
                                                    int *apc, int *aic, PS_FPT *vec )
{
    for( int g=apc[col]; g<apc[col+1]; ++g )
        vec[aic[g]] = 0;

    return;
}

__global__ void psc_cpy_mat_col_to_vec_kernel( int col, int nrow,
                                               int *apc, int *aic, PS_FPT *axc, PS_FPT *vec_ )
{
    int i = blockIdx.x;

    col = i;
    PS_FPT *vec = vec_ + i*nrow;

    for( int g=apc[col]; g<apc[col+1]; ++g )
        vec[aic[g]] = axc[g];

    return;
}

__global__ void psc_cpy_mat_diag_to_vec_kernel( int col, int nrow,
                                                int *apc, int *aic, PS_FPT *axc, PS_FPT *vec_ )
{
    int i = blockIdx.x;

    col = i;
    PS_FPT *vec = vec_;

    for( int g=apc[col]; g<apc[col+1]; ++g )
        if( col == aic[g])
            vec[aic[g]] = axc[g];

    return;
}

__global__ void sum_vec_for_trace( PS_FPT *vec, int l_side )
{
    for( int g=1; g<l_side; ++g )
        vec[0] += vec[l_side*g+g];

    return;
}









template <typename Monitor>
void report_status(Monitor& monitor)
{
    if (monitor.converged())
    {
        std::cout << "Solver converged to " << monitor.tolerance() << " tolerance";
        std::cout << " after " << monitor.iteration_count() << " iterations";
        std::cout << " (" << monitor.residual_norm() << " final residual)" << std::endl;
    }
    else
    {
        std::cout << "Solver reached iteration limit " << monitor.iteration_limit() << " before converging";
        std::cout << " to " << monitor.tolerance() << " tolerance ";
        std::cout << " (" << monitor.residual_norm() << " final residual)" << std::endl;
    }
}









void psc_fill( PS_FPT *vec, PS_FPT val, int size, cudaStream_t *str )
{
    psc_fill_kernel<<<size,1,0,*str>>>( vec, val );
/*
  thrust::device_ptr<PS_FPT> dev_ptr( vec );
  thrust::fill( dev_ptr, dev_ptr + (size_t)size, val );
*/
}

void psc_cpy_mat_col_to_vec( int col, int *apc, int *aic, PS_FPT *axc, PS_FPT *vec )
{
    psc_cpy_mat_col_to_vec_kernel <<< 1,1 >>> ( col, apc, aic, axc, vec );

    //for( int g=apc[col]; g<apc[col+1]; ++g )
    //     vec[aic[g]] = axc[g];
}

void psc_undo_cpy_mat_col_to_vec( int col, int *apc, int *aic, PS_FPT *vec )
{
    psc_undo_cpy_mat_col_to_vec_kernel <<< 1,1 >>> ( col, apc, aic, vec );

    //for( int g=apc[col]; g<apc[col+1]; ++g )
    //     vec[aic[g]] = 0;
}

PS_FPT psc_sum_vector( PS_FPT *vec, int size )
{
    thrust::device_ptr<PS_FPT> dev_ptr( vec );
    return thrust::reduce( dev_ptr, dev_ptr + (size_t)size,
                           (PS_FPT)0, thrust::plus<PS_FPT>() );
}

void psc_tri_eq_solve_cusp( int *p, int *i, PS_FPT *x,
                            PS_FPT *res_, PS_FPT *rhs_, int l_side, int nnz, void *prec )
{
    prec = 0;

    // typedef's
    typedef typename cusp::array1d_view< thrust::device_ptr<int>    > DeviceIndexArrayView;
    typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

    typedef cusp::csr_matrix_view<DeviceIndexArrayView,
        DeviceIndexArrayView,
        DeviceValueArrayView> DeviceView;

    // wrapping matrix

    thrust::device_ptr<int>    wrapped_device_Ap( p );
    thrust::device_ptr<int>    wrapped_device_Aj( i );
    thrust::device_ptr<PS_FPT> wrapped_device_Ax( x );

    DeviceIndexArrayView row_offsets   ( wrapped_device_Ap, wrapped_device_Ap + l_side+1 );
    DeviceIndexArrayView column_indices( wrapped_device_Aj, wrapped_device_Aj + nnz      );
    DeviceValueArrayView values        ( wrapped_device_Ax, wrapped_device_Ax + nnz      );

    DeviceView mat( l_side, l_side, nnz, row_offsets, column_indices, values );

    //cusp::print(mat);

    // wrapping vectors

    thrust::device_ptr<PS_FPT> wrapped_device_rhs( rhs_ );
    thrust::device_ptr<PS_FPT> wrapped_device_res( res_ );

    DeviceValueArrayView rhs( wrapped_device_rhs, wrapped_device_rhs + l_side );
    DeviceValueArrayView res( wrapped_device_res, wrapped_device_res + l_side );

    // set stopping criteria
    //  iteration_limit    = 1000
    //  relative_tolerance = 1e-3
    cusp::default_monitor<PS_FPT> monitor(rhs, 1000, 1e-6);

    // cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory> M(mat, .1);
    // cusp::krylov::cg(mat, res, rhs, monitor ,M );
    // cusp::krylov::cg(mat, res, rhs, monitor /*,M*/ );

    if( prec )
    {
        cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory> *preccer =
            (cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory>*) prec;
        cusp::krylov::cg(mat, res, rhs, monitor , *preccer );
    }
    else
    {
        cusp::krylov::cg(mat, res, rhs, monitor /*,M*/ );
    }

    //report_status(monitor);
}

void psc_destr_preconditioner_cusp( void *prec )
{
    typedef typename cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory> cusp_prec;
    cusp_prec *preccer = (cusp_prec*)prec;

    delete preccer;
    prec = 0;
}

void psc_get_preconditioner_cusp( int *p, int *i, PS_FPT *x,
                                  int l_side, int nnz, void **prec )
{


    // typedef's

    typedef typename cusp::array1d_view< thrust::device_ptr<int>    > DeviceIndexArrayView;
    typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

    typedef cusp::csr_matrix_view<DeviceIndexArrayView,
        DeviceIndexArrayView,
        DeviceValueArrayView> DeviceView;

    typedef typename cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory> cusp_prec;

    // wrapping matrix

    thrust::device_ptr<int>    wrapped_device_Ap( p );
    thrust::device_ptr<int>    wrapped_device_Aj( i );
    thrust::device_ptr<PS_FPT> wrapped_device_Ax( x );

    DeviceIndexArrayView row_offsets   ( wrapped_device_Ap, wrapped_device_Ap + l_side+1 );
    DeviceIndexArrayView column_indices( wrapped_device_Aj, wrapped_device_Aj + nnz      );
    DeviceValueArrayView values        ( wrapped_device_Ax, wrapped_device_Ax + nnz      );

    DeviceView mat( l_side, l_side, nnz, row_offsets, column_indices, values );

    cusp_prec *M = new cusp_prec(mat, 0.1);

    *prec = (void*)M;
}








float tot=0;
class matrixfree_op : public cusp::linear_operator<PS_FPT,cusp::device_memory>
{
public:
    typedef cusp::linear_operator<PS_FPT,cusp::device_memory> super;

    int N;
    cusparseHandle_t   *handle;
    cusparseMatDescr_t *descr_a;
    int *ap, *ai;
    int a_nnz;
    PS_FPT *ax;

    // constructor
    matrixfree_op( int N, void *handle_, void *descr_, int *ap_, int *ai_, PS_FPT *ax_, int a_nnz_ )
        : super(N*N,N*N),
          N(N),
          handle ((cusparseHandle_t*  )handle_),
          descr_a((cusparseMatDescr_t*)descr_ ),
          ap(ap_), ai(ai_), ax(ax_), a_nnz(a_nnz_) {}

    template <typename VectorType1,
              typename VectorType2>
    void operator()(const VectorType1& x, VectorType2& y) const
        {
            // obtain a raw pointer to device memory
            const PS_FPT* x_ptr = thrust::raw_pointer_cast(&x[0]);
            PS_FPT* y_ptr       = thrust::raw_pointer_cast(&y[0]);

            cusparseOperation_t trans_a = CUSPARSE_OPERATION_NON_TRANSPOSE;
            PS_FPT dummy_one  = 1;
            PS_FPT dummy_zero = 0;

            cudaEvent_t evt[2];
            cudaEventCreate( &evt[0] );
            cudaEventCreate( &evt[1] );
            cudaEventRecord( evt[0], 0 );

            for( int j=0; j<N; ++j )
            {
                //cusparseSetStream( *handle, stream[str] );

#ifdef SINGLE_PRECISION
                cusparseScsrmv( *handle, trans_a,
                                N, N, a_nnz, &dummy_one,
                                *descr_a,
                                ax, ap, ai,
                                x_ptr+j*N,
                                &dummy_zero, y_ptr+j*N        );
#else
                cusparseDcsrmv( *handle, trans_a,
                                N, N, a_nnz, &dummy_one,
                                *descr_a,
                                ax, ap, ai,
                                x_ptr+j*N,
                                &dummy_zero, y_ptr+j*N        );
#endif

            }

            cudaEventRecord( evt[1], 0 );
            cudaEventSynchronize( evt[1] );
            float elapsedTime;
            cudaEventElapsedTime( &elapsedTime, evt[0], evt[1] );
            cudaEventDestroy( evt[0] );
            cudaEventDestroy( evt[1] );
            tot += elapsedTime;
//          std::cout << "took " << elapsedTime << " ms" << std::endl;

        }
};

void psc_tr_inv_a_b_matrixfree( int *ap, int *ai, PS_FPT *ax,
                                int *cp, int *ci, PS_FPT *cx,
                                PS_FPT *res, PS_FPT *rhs, int l_side,
                                void *handle, void *descr, int a_nnz, void *prec  )
{
    cudaEvent_t evt[2];
    cudaEventCreate( &evt[0] );
    cudaEventCreate( &evt[1] );
    cudaEventRecord( evt[0], 0 );

    int l_side2 = l_side * l_side;

    // typedef's
    typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

    // zero out arrays
    thrust::device_ptr<PS_FPT> res_wrapped( res );
    thrust::device_ptr<PS_FPT> rhs_wrapped( rhs );
    thrust::fill( res_wrapped, res_wrapped + l_side, 1 );
    thrust::fill( rhs_wrapped, rhs_wrapped + l_side, 0 );

    // fill in rhs vector
    psc_cpy_mat_col_to_vec_kernel<<<l_side,1>>>( -1, l_side, cp, ci, cx, rhs );

    // wrapping vectors
    DeviceValueArrayView res_wrapped2( res_wrapped, res_wrapped + l_side2 );
    DeviceValueArrayView rhs_wrapped2( rhs_wrapped, rhs_wrapped + l_side2 );

    // linear operator
    matrixfree_op op( l_side, handle, descr, ap, ai, ax, a_nnz );

    cusp::default_monitor<PS_FPT> monitor(rhs_wrapped2, 1000, 1e-3);

    if( prec )
    {
        cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory> *preccer =
            (cusp::precond::scaled_bridson_ainv<PS_FPT, cusp::device_memory>*) prec;
        cusp::krylov::cg(op, res_wrapped2, rhs_wrapped2, monitor , *preccer );
    }
    else
    {
        cusp::krylov::cg(op, res_wrapped2, rhs_wrapped2, monitor);
    }

    sum_vec_for_trace<<<1,1>>>( res, l_side );

    cudaEventRecord( evt[1], 0 );
    cudaEventSynchronize( evt[1] );
    float elapsedTime;
    cudaEventElapsedTime( &elapsedTime, evt[0], evt[1] );
    cudaEventDestroy( evt[0] );
    cudaEventDestroy( evt[1] );
    std::cout << "total took " << elapsedTime << " ms " << tot << std::endl;
    //exit(1);
}





/*
  class matrixfree_op : public cusp::linear_operator<PS_FPT,cusp::device_memory>
  {
  public:
  typedef cusp::linear_operator<PS_FPT,cusp::device_memory> super;

  int N;
  cusparseHandle_t   *handle;
  cusparseMatDescr_t *descr_a;
  int *ap, *ai;
  int a_nnz;
  PS_FPT *ax;

  // constructor
  matrixfree_op( int N, void *handle_, void *descr_, int *ap_, int *ai_, PS_FPT *ax_, int a_nnz_ )
  : super(N*N,N*N),
  N(N),
  handle ((cusparseHandle_t*  )handle_),
  descr_a((cusparseMatDescr_t*)descr_ ),
  ap(ap_), ai(ai_), ax(ax_), a_nnz(a_nnz_) {}

  template <typename VectorType1,
  typename VectorType2>
  void operator()(const VectorType1& x, VectorType2& y) const
  {
  // obtain a raw pointer to device memory
  const PS_FPT* x_ptr = thrust::raw_pointer_cast(&x[0]);
  PS_FPT* y_ptr       = thrust::raw_pointer_cast(&y[0]);

  cusparseOperation_t trans_a = CUSPARSE_OPERATION_NON_TRANSPOSE;
  PS_FPT dummy_one  = 1;
  PS_FPT dummy_zero = 0;

  for( int j=0; j<N; ++j )
  {
  //cusparseSetStream( *handle, stream[str] );

  #ifdef SINGLE_PRECISION
  cusparseScsrmv( *handle, trans_a,
  N, N, a_nnz, &dummy_one,
  *descr_a,
  ax, ap, ai,
  x_ptr+j*N,
  &dummy_zero, y_ptr+j*N        );
  #else
  cusparseDcsrmv( *handle, trans_a,
  N, N, a_nnz, &dummy_one,
  *descr_a,
  ax, ap, ai,
  x_ptr+j*N,
  &dummy_zero, y_ptr+j*N        );
  #endif

  }

  }
  };

  void psc_tr_inv_a_b_matrixfree( int *ap, int *ai, PS_FPT *ax,
  int *cp, int *ci, PS_FPT *cx,
  PS_FPT *res, PS_FPT *rhs, int l_side,
  void *handle, void *descr, int a_nnz  )
  {
  int l_side2 = l_side * l_side;

  // typedef's
  typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

  // zero out arrays
  psc_fill( res, 0, l_side2 );
  psc_fill( rhs, 0, l_side2 );

  // fill in rhs vector
  psc_cpy_mat_col_to_vec_kernel<<<l_side,1>>>( -1, l_side, cp, ci, cx, rhs );

  // wrapping vectors
  thrust::device_ptr<PS_FPT> res_wrapped( res );
  thrust::device_ptr<PS_FPT> rhs_wrapped( rhs );
  DeviceValueArrayView res_wrapped2( res_wrapped, res_wrapped + l_side2 );
  DeviceValueArrayView rhs_wrapped2( rhs_wrapped, rhs_wrapped + l_side2 );

  // linear operator
  matrixfree_op op( l_side, handle, descr, ap, ai, ax, a_nnz );

  cusp::verbose_monitor<PS_FPT> monitor(rhs_wrapped2, 1000, 1e-3);

  cusp::krylov::cg(op, res_wrapped2, rhs_wrapped2, monitor);

  sum_vec_for_trace<<<1,1>>>( res, l_side );
  }

*/
