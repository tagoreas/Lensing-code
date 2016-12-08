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
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_printer.hpp"

#include "pixsrc_cuda_kernels.hpp"

#include <cstring>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cusparse_v2.h>
#include <cublas_v2.h>

//#include <cholmod.h>

#include <pthread.h>
#include <cmath>

// disabled currently
//#include <cusp/csr_matrix.h>

void pixsrc_cuda::benchmark( commoninputdata *cdata_ )
{
    // this will perform some basic matrix calculations and time them
    // benchmarking multiple GPU's at one time fails .. so just do them one at a time

    inputdata *data_;
    MEMORY ps_malloc( &data_, cdata_->numgpudevices );
    for( int dev=0; dev<cdata_->numgpudevices; ++dev )
    {
        MEMORY ps_malloc( &(data_[dev].gpu2use), 1 );
        data_[dev].numgpu2use = 1;
        data_[dev].gpu2use[0] = dev;

        if( !cdata_->gpudevices[dev].ignoreme )
            CUDA create_handle( cdata_, dev );
    }

    int m_a = 5000;
    int n_a = 5000;
    int m_b = 5000;
    int n_b = 5000;

    /*
      int m_a = 8000;
      int n_a = 7000;
      int m_b = 7000;
      int n_b = 10000;
    */

    char tr1 = 1;
    char tr2 = 0;

    int m_c = (!tr1) ? m_a : n_a;
    int n_c = (!tr2) ? n_b : m_b;

    MATRIX *a_c, **a_g, *b_c, **b_g, *c, *a_cdummy, *a_gdummy, *b_cdummy, *b_gdummy, *cdummy;

    MEMORY ps_malloc( &a_cdummy,            1          );
    MEMORY ps_malloc( &b_cdummy,            1          );
    MEMORY ps_malloc( &a_gdummy, cdata_->numgpudevices );
    MEMORY ps_malloc( &b_gdummy, cdata_->numgpudevices );
    MEMORY ps_malloc( &cdummy,              1          );

    MEMORY ps_malloc( &a_g,      cdata_->numgpudevices );
    MEMORY ps_malloc( &b_g,      cdata_->numgpudevices );


    a_c = new (a_cdummy) MATRIX( cdata_, NULL, m_a, n_a, m_a*11, 0, NULL );
    b_c = new (b_cdummy) MATRIX( cdata_, NULL, m_b, n_b, m_b*11, 0, NULL );


    for( int dev=0; dev<cdata_->numgpudevices; ++dev )
    {
        a_g[dev] = new (&(a_gdummy[dev])) MATRIX( cdata_, &(data_[dev]),
                                                  m_a, n_a, m_a*11, 0, NULL );
        b_g[dev] = new (&(b_gdummy[dev])) MATRIX( cdata_, &(data_[dev]),
                                                  m_b, n_b, m_b*11, 0, NULL );
    }

    for( int j=0; j<m_a; ++j )
    {
        for( int jj=-5; jj<=5; ++jj )
        {
            int jjj = j+jj;
            if( jjj>=0 && jjj<n_a )
            {
                a_c->set( j, jjj, jj+6 );
                for( int dev=0; dev<cdata_->numgpudevices; ++dev )
                    a_g[dev]->set( j, jjj, jj+6 );
            }
        }
    }
    for( int j=0; j<m_b; ++j )
    {
        for( int jj=-5; jj<=5; ++jj )
        {
            int jjj = j+jj;
            if( jjj>=0 && jjj<n_b )
            {
                b_c->set( j, jjj, jj+8 );
                for( int dev=0; dev<cdata_->numgpudevices; ++dev )
                    b_g[dev]->set( j, jjj, jj+8 );
            }
        }
    }

    for( int loop=0; loop<2; ++loop )
    {
        if( !loop )
            PRINTER print2screen( "pixsrc", "priming cpu and gpu(s) ..\n",
                                  cdata_->print2screenmutex   );

        for( int dev=-1; dev<cdata_->numgpudevices; ++dev )
        {

#ifdef SINGLE_PRECISION
            if( dev == -1 )
                continue;
#endif

            if( dev!=-1 && cdata_->gpudevices[dev].ignoreme )
                continue;

            if( dev != -1 )
            {
                CUDA set_device( dev );
                c = new (cdummy) MATRIX( cdata_, &(data_[dev]), m_c, n_c, 0, 0, NULL );
            }
            else
            {
                c = new (cdummy) MATRIX( cdata_, NULL,          m_c, n_c, 0, 0, NULL );
            }

            if( dev == -1 )
            {
                if( loop )
                    PRINTER print2screen( "pixsrc", "starting matrix multiplication on cpu",
                                          cdata_->print2screenmutex                         );
                cudaEvent_t start,stop;
                cudaEventCreate( &start    );
                cudaEventCreate( &stop     );
                cudaEventRecord(  start, 0 );

                a_c->mult( b_c, c, tr1, tr2, 16 );
                a_c->compile();
                a_c->creatercform();

                cudaEventRecord(stop, 0);
                cudaEventSynchronize(stop);
                float elapsedTime;
                cudaEventElapsedTime(&elapsedTime, start, stop);
                cudaEventDestroy(start);
                cudaEventDestroy(stop);

                if( loop )
                    PRINTER print2screen( "pixsrc", "done with matrix multiplication on cpu: " +
                                          OPERA tostring(elapsedTime) + " ms elapsed",
                                          cdata_->print2screenmutex                             );
            }
            else
            {
                if( loop )
                    PRINTER print2screen( "pixsrc", "starting matrix multiplication on " +
                                          string(cdata_->gpudevices[dev].name),
                                          cdata_->print2screenmutex                        );

                cudaEvent_t start,stop;
                cudaEventCreate( &start    );
                cudaEventCreate( &stop     );
                cudaEventRecord(  start, 0 );

                a_g[dev]->mult( b_g[dev], c, tr1, tr2, 0 );

                //CUDA mult( a_g[dev], b_g[dev], c, tr1, tr2 );
                //a_g[dev]->compile();
                //a_g[dev]->creatercform();

                cudaEventRecord(stop, 0);
                cudaEventSynchronize(stop);
                float elapsedTime;
                cudaEventElapsedTime(&elapsedTime, start, stop);
                cudaEventDestroy(start);
                cudaEventDestroy(stop);

                if( loop )
                    PRINTER print2screen( "pixsrc", "done with matrix multiplication on " +
                                          string(cdata_->gpudevices[dev].name) + ": " +
                                          OPERA tostring(elapsedTime) + " ms elapsed",
                                          cdata_->print2screenmutex                           );
            }

            if( ( 0 && loop ) )
            {
                std::cout << std::endl;
                string filler;

                if(dev==-1)
                {
                    for( int nr=0; nr<a_c->nrow; ++nr )
                    {
                        for( int nc=0; nc<a_c->ncol; ++nc )
                            std::cout << a_c->get(nr,nc) << " ";
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;

                    for( int nr=0; nr<b_c->nrow; ++nr )
                    {
                        for( int nc=0; nc<b_c->ncol; ++nc )
                            std::cout << b_c->get(nr,nc) << " ";
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;

                    for( int nr=0; nr<c->nrow; ++nr )
                    {
                        for( int nc=0; nc<c->ncol; ++nc )
                            std::cout << c->get(nr,nc) << " ";
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
                else
                {
                    for( int nr=0; nr<c->nrow; ++nr )
                    {
                        for( int nc=0; nc<c->ncol; ++nc )
                            std::cout << c->get(nr,nc) << " ";
                        std::cout << std::endl;
                    }
                    std::cout << std::endl;
                }
            }

            c->~MATRIX();
        }
    }

    for( int dev=0; dev<cdata_->numgpudevices; ++dev )
    {
        a_g[dev]->~MATRIX();
        b_g[dev]->~MATRIX();
    }
    a_c->~MATRIX();
    b_c->~MATRIX();

    MEMORY ps_free( a_cdummy );
    MEMORY ps_free( b_cdummy );
    MEMORY ps_free( a_gdummy );
    MEMORY ps_free( b_gdummy );
    MEMORY ps_free(   cdummy );

    for( int dev=0; dev<cdata_->numgpudevices; ++dev )
    {
        MEMORY ps_free( data_[dev].gpu2use );
    }
    MEMORY ps_free( data_  );
}

void pixsrc_cuda::detectgpus( commoninputdata *cdata_ )
{
    if (cdata_->gpudevices)
        return;

    // read in GPU devices
    cudaError_t error_id = cudaGetDeviceCount( &cdata_->numgpudevices );
    CUDA errorcheck( cdata_, (void*)&error_id, 0, "device count (possibly no gpu's found)" );

    // must have at least one GPU
    if( !cdata_->numgpudevices )
    {
        PRINTER printerror("pixsrc", "no gpu devices found", cdata_->print2screenmutex);
    }

    // read in GPU names
    MEMORY ps_malloc( &(cdata_->gpudevices), cdata_->numgpudevices );
    for ( int dev=0; dev<cdata_->numgpudevices; ++dev )
    {
        //cusparseStatus_t  status;

        // set device
        CUDA set_device( dev );

        // get properties
        cudaDeviceProp deviceProp;
        cudaGetDeviceProperties(&deviceProp, dev);

        // get device name
        MEMORY ps_malloc( &(cdata_->gpudevices[dev].name), OPERA sizestring(deviceProp.name) +1 );
        strcpy( cdata_->gpudevices[dev].name, deviceProp.name );

        cdata_->gpudevices[dev].busid = deviceProp.pciBusID;
        cdata_->gpudevices[dev].locid = deviceProp.pciDeviceID;

        // set device to be possibly used
        cdata_->gpudevices[dev].ignoreme = 0;

        // setting number of sumultaneously kernels to maximum
        // allowed by CUDA version 5.0
        cdata_->gpudevices[dev].num_kernels = 16;

        // convert to lowercase
        char *p = cdata_->gpudevices[dev].name;
        for( ; *p; ++p )
            *p = tolower(*p);
        // remove whitespace
        OPERA trim( cdata_->gpudevices[dev].name, NULL );
        p = cdata_->gpudevices[dev].name;
        for( ; *p; ++p )
            if( *p == ' ' )
                *p = '_';
    }
}

void pixsrc_cuda::sync_vector( VECTOR *a, cuda_cpy_direction cpy )
{
    // sync data between CPU and GPU
    // cpy determines the direction of the copying

    // lock mutex
    pthread_mutex_lock( &(a->pv_mutex) );

    // if both CPU and GPU are up to date, return
    if( a->up2date == pvs_both )
    {
        pthread_mutex_unlock( &(a->pv_mutex) );
        return;
    }

    // set GPU device to use
    commoninputdata *cdata_ = a->cdata_;
    inputdata *data_ = a->data_;
    int dev = data_->gpu2use[0];

    // create CPU vector is it doens't exist
    bool makenew = 0;
    if( !a->vec )
    {
        MEMORY ps_malloc( &a->vec, a->size );
        makenew = 1;
    }

    // sanity check
    if( makenew && ( a->up2date == pvs_cpu || a->up2date == pvs_both ) )
    {
        PRINTER printerror( data_->print2screenname,
                            "pixsrc_vector up2date",
                            cdata_->print2screenmutex );
    }

    // if we want to copy memory to gpu
    if( cpy == ccd_cpu2gpu )
    {
        CUDA ps_cuda_memcpy( cdata_, a->cuda_vec, a->vec, a->size, 0, (cudaStream_t*)a->stream, dev );
    }
    // if we want to copy memory to cpu
    else if( cpy == ccd_gpu2cpu )
    {
        CUDA ps_cuda_memcpy( cdata_, a->vec, a->cuda_vec, a->size, 1, (cudaStream_t*)a->stream, dev );
    }

    // update and unlock mutex
    a->up2date = pvs_both;
    pthread_mutex_unlock( &(a->pv_mutex) );

}

void pixsrc_cuda::sync_matrix( MATRIX *a, pm_format format, cuda_cpy_direction cpy )
{
    // sync data between CPU and GPU
    // cpy determines the direction of the copying

    // lock mutex
    pthread_mutex_lock( &a->pm_mutex );

    // if both CPU and GPU are up to date, return
    if( format == pmf_csr && a->up2date_csr == pms_both )
    {
        pthread_mutex_unlock( &a->pm_mutex );
        return;
    }
    if( format == pmf_csc && a->up2date_csc == pms_both )
    {
        pthread_mutex_unlock( &a->pm_mutex );
        return;
    }

    // set GPU device
    commoninputdata *cdata_ = a->cdata_;
    inputdata *data_ = a->data_;
    int dev = data_->gpu2use[0];

    bool make_new_cpu = 0;
    bool make_new_gpu = 0;

    // CSR format (row compressed)
    if( format == pmf_csr )
    {
        // create CSR format
        if( !a->rcform )
            a->creatercform();

        // create vectors on CPU if they don't exist
        if( !a->air )
        {
            MEMORY ps_malloc( &a->air, a->nnz );
            make_new_cpu = 1;
        }
        if( !a->apr )
        {
            MEMORY ps_malloc( &a->apr, a->nrow+1 );
            make_new_cpu = 1;
        }
        if( !a->axr )
        {
            MEMORY ps_malloc( &a->axr, a->nnz );
            make_new_cpu = 1;
        }

        // create vectors on GPU if they don't exist
        if( !a->coo_col )
        {
            CUDA ps_cuda_malloc( cdata_, &a->coo_col, a->nnz, dev );
            make_new_gpu = 1;
        }
        if( !a->coo_val )
        {
            CUDA ps_cuda_malloc( cdata_, &a->coo_val, a->nnz, dev );
            make_new_gpu = 1;
        }
        if( !a->csx )
        {
            CUDA ps_cuda_malloc( cdata_, &a->csx, a->nrow+1, dev );
            make_new_gpu = 1;
        }

        // sanity check
        if( make_new_gpu && ( a->up2date_csr == pms_gpu || a->up2date_csr == pms_both ) )
        {
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix up2date 1",
                                cdata_->print2screenmutex );
        }
        if( make_new_cpu && ( a->up2date_csr == pms_cpu || a->up2date_csr == pms_both ) )
        {
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix up2date 2",
                                cdata_->print2screenmutex );
        }

        // if we want to copy memory to gpu
        if( cpy == ccd_cpu2gpu )
        {
            CUDA ps_cuda_memcpy( cdata_, a->coo_col, a->air, a->nnz,    0, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->coo_val, a->axr, a->nnz,    0, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->csx,     a->apr, a->nrow+1, 0, (cudaStream_t*)a->stream, dev );
        }

        // if we want to copy memory to cpu
        else if( cpy == ccd_gpu2cpu )
        {
            CUDA ps_cuda_memcpy( cdata_, a->air, a->coo_col, a->nnz,    1, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->axr, a->coo_val, a->nnz,    1, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->apr, a->csx,     a->nrow+1, 1, (cudaStream_t*)a->stream, dev );
        }

        // update
        a->up2date_csr = pms_both;

    }
    // CSC format (column compressed)
    else
    {
        // create CSC format
        if( !a->hasbeencompiled )
            a->compile();

        // create vectors on CPU if they don't exist
        if( !a->aic )
        {
            MEMORY ps_malloc( &a->aic, a->nnz );
            make_new_cpu = 1;
        }
        if( !a->apc )
        {
            MEMORY ps_malloc( &a->apc, a->ncol+1 );
            make_new_cpu = 1;
        }
        if( !a->axc )
        {
            MEMORY ps_malloc( &a->axc, a->nnz );
            make_new_cpu = 1;
        }

        // create vectors on GPU if they don't exist
        if( !a->coo_row_tr )
        {
            CUDA ps_cuda_malloc( cdata_, &a->coo_row_tr, a->nnz, dev );
            make_new_gpu = 1;
        }
        if( !a->coo_val_tr )
        {
            CUDA ps_cuda_malloc( cdata_, &a->coo_val_tr, a->nnz, dev );
            make_new_gpu = 1;
        }
        if( !a->csx_tr )
        {
            CUDA ps_cuda_malloc( cdata_, &a->csx_tr, a->ncol+1, dev );
            make_new_gpu = 1;
        }

        // sanity check
        if( make_new_gpu && ( a->up2date_csc == pms_gpu || a->up2date_csc == pms_both ) )
        {
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix up2date 3",
                                cdata_->print2screenmutex );
        }
        if( make_new_cpu && ( a->up2date_csc == pms_cpu || a->up2date_csc == pms_both ) )
        {
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix up2date 4",
                                cdata_->print2screenmutex );
        }

        // if we want to copy memory to gpu
        if( cpy == ccd_cpu2gpu )
        {
            CUDA ps_cuda_memcpy( cdata_, a->coo_row_tr, a->aic, a->nnz,    0, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->coo_val_tr, a->axc, a->nnz,    0, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->csx_tr,     a->apc, a->ncol+1, 0, (cudaStream_t*)a->stream, dev );
        }

        // if we want to copy memory to cpu
        else if( cpy == ccd_gpu2cpu )
        {
            CUDA ps_cuda_memcpy( cdata_, a->aic, a->coo_row_tr, a->nnz,    1, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->axc, a->coo_val_tr, a->nnz,    1, (cudaStream_t*)a->stream, dev );
            CUDA ps_cuda_memcpy( cdata_, a->apc, a->csx_tr,     a->ncol+1, 1, (cudaStream_t*)a->stream, dev );
        }

        // update
        a->up2date_csc = pms_both;

    }

    // unlock mutex
    pthread_mutex_unlock( &a->pm_mutex );

}

void pixsrc_cuda::postprocess( MATRIX *a, MATRIX *b, MATRIX *c )
{
    // this function does post-processing on matrices
    // created by CUDA mult and CUDA plus

    // diagonal matrices are processed by cpu and not gpu

    // commoninputdata *cdata_ = c->cdata_;

    // the product of two diagonal matrices is another diagonal
    if( (a->matrix_type == 10 /*|| a->matrix_type == 11*/) &&
        (b->matrix_type == 10 /*|| b->matrix_type == 11*/) )
        c->matrix_type = 10;

    if( c->matrix_type == 10)
    {
        // for a diagonal matrix, we can easily create CSC format from CSR format

        c->coo_row_tr = c->coo_col;
        c->coo_val_tr = c->coo_val;
        c->csx_tr     = c->csx;

        c->hasbeencompiled = 1;
        c->up2date_csc = pms_gpu;
    }

    /*
      if( c->matrix_type == 10 /-*|| c->matrix_type == 11*-/ )
      {
      if( c->tkv && c->vecsize < c->nnz )
      {
      MEMORY ps_free( c->tkv );
      c->tkv = 0;
      }

      if( !c->tkv )
      {
      MEMORY ps_malloc( &c->tkv, c->nnz );
      }

      CUDA ps_cuda_memcpy( cdata_, c->tkv, c->coo_val, c->nnz, 1 );

      return;
      }
    */

    // processing non-diagonal matrices now

    /*
      if( !c->air )
      MEMORY ps_malloc( &c->air, c->nnz    );
      if( !c->axr )
      MEMORY ps_malloc( &c->axr, c->nnz    );
      if( !c->apr )
      MEMORY ps_malloc( &c->apr, c->nrow+1 );

      CUDA ps_cuda_memcpy( cdata_, c->air, c->coo_col, c->nnz,    1 );
      CUDA ps_cuda_memcpy( cdata_, c->apr, c->csx,     c->nrow+1, 1 );
      CUDA ps_cuda_memcpy( cdata_, c->axr, c->coo_val, c->nnz,    1 );
    */

    c->rcform        = 1;
    c->icf_done      = 0;
    //c->csrsvi_done = 0;
    c->up2date_csr = pms_gpu;
}

void pixsrc_cuda::compile( MATRIX *a )
{
    // this function is only called from pixsrc_matrix so
    // there's no need to put a lock on a mutex because the
    // mutex should already be locked
    // also, creatercform() is assumed to have already been
    // called from pixsrc_matrix before calling this function

    // CUDA create CSR
    // this creates a CSC format of the matrix.

    // if CSC already exists, return
    if( a->hasbeencompiled )
        return;

    // update the GPU
    a->update_gpu( pmf_csr );

    // set GPU and handles
    int dev                    = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata  *cdata_   = a->cdata_;
    cusparseHandle_t *handle   = (cusparseHandle_t*)(cdata_->gpudevices[dev].handle);

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cusparse_lock );

    cusparseStatus_t status;

    // allocate memory on gpu
    if( !a->coo_row_tr )
        CUDA ps_cuda_malloc( cdata_, &a->coo_row_tr, a->nnz, dev );

    if( !a->coo_val_tr )
        CUDA ps_cuda_malloc( cdata_, &a->coo_val_tr, a->nnz, dev );

    if( !a->csx_tr     )
        CUDA ps_cuda_malloc( cdata_, &a->csx_tr, a->ncol + 1, dev );

    cudaStream_t *stream = (cudaStream_t*)a->stream;
    cusparseSetStream( *handle, *stream );

    // perform format conversion
#ifdef SINGLE_PRECISION
    status = cusparseScsr2csc( *handle, a->nrow, a->ncol, a->nnz,
                               a->coo_val, a->csx, a->coo_col,
                               a->coo_val_tr, a->coo_row_tr, a->csx_tr,
                               CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO );
#else
    status = cusparseDcsr2csc( *handle, a->nrow, a->ncol, a->nnz,
                               a->coo_val, a->csx, a->coo_col,
                               a->coo_val_tr, a->coo_row_tr, a->csx_tr,
                               CUSPARSE_ACTION_NUMERIC, CUSPARSE_INDEX_BASE_ZERO );
#endif

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cusparse_lock );

    // error check
    CUDA errorcheck( cdata_, (void*)&status, 1, "csr2csc" );

    // update
    a->up2date_csc = pms_gpu;
    a->hasbeencompiled = 1;
}

void pixsrc_cuda::creatercform( MATRIX *a )
{
    // this function is only called from pixsrc_matrix so
    // there's no need to put a lock on a mutex because the
    // mutex should already be locked

    // if we're using CUDA, then the only way to create csr format
    // is to process the t?v COO format vectors

    // if CSR exists, return
    if( a->rcform )
        return;

    // set GPU device and handles
    int dev                    = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata  *cdata_   = a->cdata_;
    cusparseHandle_t *handle   = (cusparseHandle_t*)(cdata_->gpudevices[dev].handle);

    cusparseStatus_t status;

    // allocate memory on gpu
    if( !a->coo_col )
        CUDA ps_cuda_malloc( cdata_, &a->coo_col, a->nnz, dev );

    if( !a->coo_val )
        CUDA ps_cuda_malloc( cdata_, &a->coo_val, a->nnz, dev );

    if( !a->csx     )
        CUDA ps_cuda_malloc( cdata_, &a->csx, a->nrow + 1, dev );

    int *temp;
    CUDA ps_cuda_malloc( cdata_, &temp, a->nnz, dev );

    // copy data to gpu memory
    CUDA ps_cuda_memcpy( cdata_, a->coo_col, a->tjv, a->nnz, 0, (cudaStream_t*)a->stream, dev );
    CUDA ps_cuda_memcpy( cdata_, temp,       a->tiv, a->nnz, 0, (cudaStream_t*)a->stream, dev );
    CUDA ps_cuda_memcpy( cdata_, a->coo_val, a->tkv, a->nnz, 0, (cudaStream_t*)a->stream, dev );

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cusparse_lock );

    cudaStream_t *stream = (cudaStream_t*)a->stream;
    cusparseSetStream( *handle, *stream );

    // perform format conversion
    status = cusparseXcoo2csr( *handle, temp, a->nnz, a->nrow,
                               a->csx, CUSPARSE_INDEX_BASE_ZERO );

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cusparse_lock );

    // error check
    CUDA errorcheck( cdata_, (void*)&status, 1, "coo2csr" );

    CUDA ps_cuda_free( temp );

    // update
    a->up2date_csr = pms_gpu;
    a->rcform = 1;
}

void pixsrc_cuda::create_handle( commoninputdata *cdata_, int dev )
{
    // This function creates handles to CUDA libraries
    // It also creates mutex locks, so multiple threads don't cause problems.
    // It also creates descriptors for the matrices.

    cusparseStatus_t  status;

    CUDA set_device( dev );

    // get handle
    cusparseHandle_t *handle;
    MEMORY ps_malloc( &handle, 1 );
    status = cusparseCreate( handle );
    CUDA errorcheck( cdata_, (void*)&status, 1, "create handle"     );
    cdata_->gpudevices[dev].handle = (void*)handle;


    // get cublas handle
    cublasHandle_t *cublashandle;
    MEMORY ps_malloc( &cublashandle, 1 );
    cublasStatus_t stat = cublasCreate( cublashandle );
    CUDA errorcheck( cdata_, (void*)&stat, 2, "cublas handle" );
    cdata_->gpudevices[dev].handle_cublas = (void*)cublashandle;

    // get handle locks
    MEMORY ps_malloc( &cdata_->gpudevices[dev].cusparse_lock, 1 );
    MEMORY ps_malloc( &cdata_->gpudevices[dev].cublas_lock,   1 );
    pthread_mutex_init( cdata_->gpudevices[dev].cusparse_lock, NULL );
    pthread_mutex_init( cdata_->gpudevices[dev].cublas_lock,   NULL );

    // create descriptors
    cusparseMatDescr_t *descr;
    MEMORY ps_malloc( &descr, 1 );
    status = cusparseCreateMatDescr( descr );
    CUDA errorcheck( cdata_, (void*)&status, 1, "create descriptor" );
    cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_GENERAL   );
    cusparseSetMatIndexBase( *descr, CUSPARSE_INDEX_BASE_ZERO       );
    cdata_->gpudevices[dev].descr = (void*)descr;

    cusparseMatDescr_t *descr_sym;
    MEMORY ps_malloc( &descr_sym, 1 );
    status = cusparseCreateMatDescr( descr_sym );
    CUDA errorcheck( cdata_, (void*)&status, 1, "create descriptor"     );
    cusparseSetMatType     ( *descr_sym, CUSPARSE_MATRIX_TYPE_SYMMETRIC );
    cusparseSetMatIndexBase( *descr_sym, CUSPARSE_INDEX_BASE_ZERO       );
    cusparseSetMatFillMode ( *descr_sym, CUSPARSE_FILL_MODE_UPPER       );
    cdata_->gpudevices[dev].descr_sym = (void*)descr_sym;

    cusparseMatDescr_t *descr_tri_upper;
    MEMORY ps_malloc( &descr_tri_upper, 1 );
    status = cusparseCreateMatDescr( descr_tri_upper );
    CUDA errorcheck( cdata_, (void*)&status, 1, "create descriptor"     );
    cusparseSetMatType     ( *descr_tri_upper, CUSPARSE_MATRIX_TYPE_TRIANGULAR );
    cusparseSetMatIndexBase( *descr_tri_upper, CUSPARSE_INDEX_BASE_ZERO        );
    cusparseSetMatFillMode ( *descr_tri_upper, CUSPARSE_FILL_MODE_UPPER        );
    cdata_->gpudevices[dev].descr_tri_upper = (void*)descr_tri_upper;

}

void pixsrc_cuda::set_device( int dev )
{
    // set GPU device
    cudaSetDevice( dev );
}

void pixsrc_cuda::errorcheck( commoninputdata *cdata_, void *var_, int which, const char *errstr )
{
    // var    = the status code
    // which  = the CUDA library used to generate status code
    // errstr = an error message to print

    // check for error
    bool err = 0;
    string code;
    if( which == 0 )
    {
        cudaError_t cudaStat = *((cudaError_t*)var_);

        if( cudaStat != cudaSuccess )
        {
            err = 1;
            code = OPERA tostring( cudaStat );
        }
    }
    else if( which == 1 )
    {
        cusparseStatus_t status = *((cusparseStatus_t*)var_);

        if( status != CUSPARSE_STATUS_SUCCESS )
        {
            err = 1;
            code = OPERA tostring( status );
        }
    }
    else if( which == 2 )
    {
        cublasStatus_t status = *((cublasStatus_t*)var_);

        if( status != CUBLAS_STATUS_SUCCESS )
        {
            err = 1;
            code = OPERA tostring( status );
        }
    }

    // print error
    if( err )
    {
        PRINTER printerror( "pixsrc",
                            "cuda error: " + string(errstr) + ": code " + code,
                            cdata_->print2screenmutex                           );
    }

}

// memory allocation and de-allocation code

void pixsrc_cuda::ps_cuda_free( int *vec )
{
    if( vec )
        cudaFree( vec );
    vec = 0;
}

void pixsrc_cuda::ps_cuda_free( PS_FPT *vec )
{
    if( vec )
        cudaFree( vec );
    vec = 0;
}

void pixsrc_cuda::ps_cuda_malloc( commoninputdata *cdata_, int **vec, int size, int dev )
{
    CUDA set_device( dev );

    cudaError_t cudaStat =
        cudaMalloc( (void**)&*vec, size * sizeof(PS_SIT) );

    CUDA errorcheck( cdata_, (void*)&cudaStat, 0, "device malloc int" );
}

void pixsrc_cuda::ps_cuda_malloc( commoninputdata *cdata_, PS_FPT **vec, int size, int dev )
{
    CUDA set_device( dev );

    cudaError_t cudaStat =
        cudaMalloc( (void**)&*vec, size * sizeof(PS_FPT) );

    CUDA errorcheck( cdata_, (void*)&cudaStat, 0, "device malloc PS_FPT" );
}

// Note: the following memcpy codes could be combined into a template function
void pixsrc_cuda::ps_cuda_memcpy( commoninputdata *cdata_, int *dest, int *src, int size, int dir, cudaStream_t *str, int dev )
{
    // copy memory between CPU and GPU

    // set GPU device
    CUDA set_device( dev );

    // figure out the direction in which to copy data
    cudaMemcpyKind cuda_dir = (cudaMemcpyKind)0;
    switch( dir )
    {
    case 0: cuda_dir = cudaMemcpyHostToDevice;
        break;
    case 1: cuda_dir = cudaMemcpyDeviceToHost;
        break;
    case 2: cuda_dir = cudaMemcpyDeviceToDevice;
        break;
    default:
        PRINTER printerror( "cuda", "memcpy error", cdata_->print2screenmutex );
        break;
    }

    cudaError_t cudaStat;

    // copy
    if( str )
        cudaStat = cudaMemcpyAsync( dest, src, (size_t)( size * sizeof(PS_SIT) ), cuda_dir, *str );
    else
        cudaStat = cudaMemcpy( dest, src, (size_t)( size * sizeof(PS_SIT) ), cuda_dir );

    // error check
    CUDA errorcheck( cdata_, (void*)&cudaStat, 0, "memcpy int" );
}

void pixsrc_cuda::ps_cuda_memcpy( commoninputdata *cdata_, PS_FPT *dest, PS_FPT *src, int size, int dir, cudaStream_t *str, int dev )
{
    // see the "int" version of this function for comments

    CUDA set_device( dev );

    cudaMemcpyKind cuda_dir = (cudaMemcpyKind)0;
    switch( dir )
    {
    case 0: cuda_dir = cudaMemcpyHostToDevice;
        break;
    case 1: cuda_dir = cudaMemcpyDeviceToHost;
        break;
    case 2: cuda_dir = cudaMemcpyDeviceToDevice;
        break;
    default:
        PRINTER printerror( "cuda", "memcpy error", cdata_->print2screenmutex );
        break;
    }

    cudaError_t cudaStat;

    if( str )
        cudaStat = cudaMemcpyAsync( dest, src, (size_t)( size * sizeof(PS_FPT) ), cuda_dir, *str );
    else
        cudaStat = cudaMemcpy( dest, src, (size_t)( size * sizeof(PS_FPT) ), cuda_dir );

    CUDA errorcheck( cdata_, (void*)&cudaStat, 0, "memcpy PS_FPT" );
}

// code for timing GPU and/or CPU calculations

void pixsrc_cuda::ps_start_clock( void **ptr_ )
{
    cudaEvent_t **ptr = (cudaEvent_t**)ptr_;

    MEMORY ps_malloc( &(*ptr), 2 );

    cudaEventCreate( &((*ptr)[0]) );
    cudaEventCreate( &((*ptr)[1]) );

    cudaEventRecord(  (*ptr)[0], 0 );
}

void pixsrc_cuda::ps_stop_clock( void *ptr_, float *ans )
{
    cudaEvent_t *ptr = (cudaEvent_t*)ptr_;

    cudaEventRecord( ptr[1], 0 );
    cudaEventSynchronize( ptr[1] );

    float elapsedTime;

    cudaEventElapsedTime( &elapsedTime, ptr[0], ptr[1] );
    cudaEventDestroy( ptr[0] );
    cudaEventDestroy( ptr[1] );

    MEMORY ps_free( ptr );

    *ans = elapsedTime;
}

void pixsrc_cuda::fill_vec( VECTOR *a, PS_FPT val )
{
    // fill CUDA vector with a constant

    //commoninputdata *cdata_ = a->cdata_;

    psc_fill( a->cuda_vec, val, a->size, (cudaStream_t*)a->stream );
    a->scalar = 1.0;
    a->up2date = pvs_gpu;
}

void pixsrc_cuda::dissolve_scalar( VECTOR *a )
{
    // multiply each element of vector by scalar and set scalar to 1

    // if already dissolved, return
    if( a->scalar == 1.0 )
        return;

    // set GPU device and handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;
    cublasHandle_t *handle_cublas  = (cublasHandle_t*)cdata_->gpudevices[dev].handle_cublas;
    cublasStatus_t stat;

    // update GPU data
    a->update_gpu();

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cublas_lock );

    // multiply vector
    cudaStream_t *stream = (cudaStream_t*)a->stream;
    cublasSetStream( *handle_cublas, *stream );

#ifdef SINGLE_PRECISION
    stat = cublasSscal( *handle_cublas, a->size, &a->scalar, a->cuda_vec, 1 );
#else
    stat = cublasDscal( *handle_cublas, a->size, &a->scalar, a->cuda_vec, 1 );
#endif

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cublas_lock );

    // error check
    CUDA errorcheck( cdata_, (void*)&stat, 2, "dissolve scalar" );

    // update
    a->scalar = 1.0;
    a->up2date = pvs_gpu;
}

void pixsrc_cuda::dissolve_scalar( MATRIX *a, pm_format format )
{
    // multiply each element of matrix by scalar and set scalar to 1
    // no resetting of up2date_<format> is done here
    // because this function is only (as of now) being called from
    // MATRIX::dissolve_vector, which DOES NOT change pm_status

    // if already dissolved, return
    if( a->scalar == 1.0 )
        return;

    // set GPU device and handles
    int dev = a->data_->gpu2use[0];
    CUDA set_device( dev );
    commoninputdata *cdata_ = a->cdata_;
    cublasHandle_t *handle_cublas  = (cublasHandle_t*)cdata_->gpudevices[dev].handle_cublas;
    cublasStatus_t stat;

    // update GPU data
    a->update_gpu( format );

    // lock mutex
    pthread_mutex_lock( cdata_->gpudevices[dev].cublas_lock );

    // multiply matrix values for CSR and CSC formats
    cudaStream_t *stream = (cudaStream_t*)a->stream;
    cublasSetStream( *handle_cublas, *stream );

#ifdef SINGLE_PRECISION
    if( a->hasbeencompiled )
        stat = cublasSscal( *handle_cublas, a->nnz, &a->scalar, a->coo_val_tr, 1 );
    if( a->rcform )
        stat = cublasSscal( *handle_cublas, a->nnz, &a->scalar, a->coo_val,    1 );
#else
    if( a->hasbeencompiled )
        stat = cublasDscal( *handle_cublas, a->nnz, &a->scalar, a->coo_val_tr, 1 );
    if( a->rcform )
        stat = cublasDscal( *handle_cublas, a->nnz, &a->scalar, a->coo_val,    1 );
#endif

    // unlock mutex
    pthread_mutex_unlock( cdata_->gpudevices[dev].cublas_lock );

    // error check
    CUDA errorcheck( cdata_, (void*)&stat, 2, "dissolve scalar" );

    a->scalar = 1.0;
}

void pixsrc_cuda::wait_for_stream( MATRIX *a )
{
    // wait for CUDA stream to finish calculations
    cudaStream_t *str = (cudaStream_t*)a->stream;
    cudaStreamSynchronize( *str );
}

void pixsrc_cuda::wait_for_stream( VECTOR *a )
{
    // wait for CUDA stream to finish calculations
    cudaStream_t *str = (cudaStream_t*)a->stream;
    cudaStreamSynchronize( *str );
}

void pixsrc_cuda::create_stream( MATRIX *a )
{
    // create CUDA stream
    cudaStream_t *str;
    MEMORY ps_malloc( &str, 1 );
    cudaStreamCreate( str );
    a->stream = (void*)str;
}

void pixsrc_cuda::create_stream( VECTOR *a )
{
    // create CUDA stream
    cudaStream_t *str;
    MEMORY ps_malloc( &str, 1 );
    cudaStreamCreate( str );
    a->stream = (void*)str;
}

void pixsrc_cuda::destroy_stream( MATRIX *a )
{
    // destroy CUDA stream

    if( !a->stream )
        return;

    cudaStream_t *str = (cudaStream_t*)a->stream;
    cudaStreamDestroy( *str );
    MEMORY ps_free( str );
}

void pixsrc_cuda::destroy_stream( VECTOR *a )
{
    // destroy CUDA stream
    if( !a->stream )
        return;

    cudaStream_t *str = (cudaStream_t*)a->stream;
    cudaStreamDestroy( *str );
    MEMORY ps_free( str );
}


// CUSP preconditioning is disabled for now
void pixsrc_cuda::destroy_cusp_prec( MATRIX *a )
{
//    if( a->cusp_prec )
    //       psc_destr_preconditioner_cusp( a->cusp_prec );
}





///////////////////////////
// CURRENTLY UNUSED CODE //
///////////////////////////

// it's not all useless, fyi

/*
  void pixsrc_cuda::ps_dest_csrsvi( MATRIX *a )
  {
  void *info_ = a->csrsvi;
  commoninputdata *cdata_ = a->cdata_;

  if( !info_ )
  return;

  cusparseStatus_t status;
  cusparseSolveAnalysisInfo_t *info = (cusparseSolveAnalysisInfo_t*)info_;

  status = cusparseDestroySolveAnalysisInfo( *info );

  CUDA errorcheck( cdata_, (void*)&status, 1, "destroy csrsv info" );

  MEMORY ps_free( info );
  }
*/
/*
  void pixsrc_cuda::ps_set_descr( MATRIX *a )
  {
  commoninputdata *cdata_ = a->cdata_;
  int matrix_type = a->matrix_type;

  cusparseMatDescr_t *descr;
  MEMORY ps_malloc( &descr, 1 );

  cusparseStatus_t status = cusparseCreateMatDescr( descr );

  CUDA errorcheck( cdata_, (void*)&status, 1, "create descriptor" );

  cusparseSetMatIndexBase( *descr, CUSPARSE_INDEX_BASE_ZERO     );

  switch( matrix_type )
  {

  case 0:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_GENERAL    );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_NON_UNIT     );
  break;

  case 1:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_SYMMETRIC  );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_NON_UNIT     );
  cusparseSetMatFillMode ( *descr, CUSPARSE_FILL_MODE_UPPER        );
  break;

  case 2:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_HERMITIAN  );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_NON_UNIT     );
  cusparseSetMatFillMode ( *descr, CUSPARSE_FILL_MODE_UPPER        );
  break;

  case 3:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_TRIANGULAR );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_NON_UNIT     );
  cusparseSetMatFillMode ( *descr, CUSPARSE_FILL_MODE_UPPER        );
  break;

  case 4:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_TRIANGULAR );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_UNIT         );
  cusparseSetMatFillMode ( *descr, CUSPARSE_FILL_MODE_UPPER        );
  break;

  case 10:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_TRIANGULAR );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_NON_UNIT     );
  cusparseSetMatFillMode ( *descr, CUSPARSE_FILL_MODE_UPPER        );
  break;

  case 11:
  cusparseSetMatType     ( *descr, CUSPARSE_MATRIX_TYPE_TRIANGULAR );
  cusparseSetMatDiagType ( *descr, CUSPARSE_DIAG_TYPE_UNIT         );
  cusparseSetMatFillMode ( *descr, CUSPARSE_FILL_MODE_UPPER        );
  break;

  default:
  break;

  }

  a->descr = (void*)descr;

  }
*/
/*
  void pixsrc_cuda::ps_dest_descr( MATRIX *a )
  {
  commoninputdata    *cdata_ = a->cdata_;
  cusparseMatDescr_t *descr  = (cusparseMatDescr_t*)a->descr;

  cusparseStatus_t status = cusparseDestroyMatDescr( *descr );

  CUDA errorcheck( cdata_, (void*)&status, 1, "create descriptor" );

  MEMORY ps_free( descr );
  }
*/

/*
  void pixsrc_cuda::copy_csc_c2g_h( MATRIX *a )
  {
  commoninputdata *cdata_ = a->cdata_;

  // allocate memory

  if( !a->coo_row_tr )
  CUDA ps_cuda_malloc( cdata_, &a->coo_row_tr, a->nnz      );

  if( !a->coo_val_tr )
  CUDA ps_cuda_malloc( cdata_, &a->coo_val_tr, a->nnz      );

  if( !a->csx_tr )
  CUDA ps_cuda_malloc( cdata_, &a->csx_tr,     a->ncol + 1 );

  // copy data

  CUDA ps_cuda_memcpy( cdata_, a->coo_row_tr, a->aic, a->nnz,      0 );
  CUDA ps_cuda_memcpy( cdata_, a->coo_val_tr, a->axc, a->nnz,      0 );
  CUDA ps_cuda_memcpy( cdata_, a->csx_tr,     a->apc, a->ncol + 1, 0 );

  a->up2date = 0;
  }
*/

/*
  void pixsrc_cuda::copy_csr_c2g_s( MATRIX *a )
  {
  // get rid of  memory if already allocated

  if( a->coo_col )
  CUDA ps_cuda_free( a->coo_col );

  if( a->coo_val )
  CUDA ps_cuda_free( a->coo_val );

  if( a->csx )
  CUDA ps_cuda_free( a->csx     );

  // create soft links

  a->coo_col = a->coo_row_tr;
  a->coo_val = a->coo_val_tr;
  a->csx     = a->csx_tr;

  a->up2date_csr = pms_both;
  }
*/
/*
  void pixsrc_cuda::create_sym_form( MATRIX *a )
  {
  if( a->is_sym_good )
  return;

  // commoninputdata *cdata_ = a->cdata_;

  int l_side = a->nrow;
  int nnz_sym = ( a->nnz - l_side ) / 2 + l_side;

  a->update_cpu( pmf_csr );

  // Converting form matrix_type_general to matrix_type_symmetric
  // create structure for symmetric matrix

  bool copy_structure = 0;

  // cpu check
  if( !a->apr_sym )
  {
  MEMORY ps_malloc( &a->apr_sym, l_side+1 );
  copy_structure = 1;
  }
  if( !a->air_sym )
  {
  MEMORY ps_malloc( &a->air_sym, nnz_sym  );
  copy_structure = 1;
  }
  if( !a->axr_sym )
  {
  MEMORY ps_malloc( &a->axr_sym, nnz_sym  );
  }

  /-*
  // gpu check
  if( !a->coo_col_sym )
  {
  CUDA ps_cuda_malloc( cdata_, &a->coo_col_sym, nnz_sym  );
  copy_structure = 1;
  }
  if( !a->coo_val_sym )
  {
  CUDA ps_cuda_malloc( cdata_, &a->coo_val_sym, nnz_sym  );
  }
  if( !a->csx_sym )
  {
  CUDA ps_cuda_malloc( cdata_, &a->csx_sym,     l_side+1 );
  copy_structure = 1;
  }
  *-/

  // if copying data to gpu
  if( copy_structure )
  {
  // create symmetric matrix on cpu
  a->apr_sym[0] = 0;
  int c, index = 0;
  for( int r=0; r<l_side; ++r )
  {
  a->apr_sym[r+1] = a->apr_sym[r];
  for( int i=a->apr[r]; i<a->apr[r+1]; ++i )
  {
  c = a->air[i];
  // if in upper triangle
  if( c >= r )
  {
  a->axr_sym[index] = a->axr[i];
  a->air_sym[index] = a->air[i];
  ++a->apr_sym[r+1];
  ++index;
  }
  }
  }

  // copy symmetric matrix structure to gpu
  //CUDA ps_cuda_memcpy( cdata_, a->coo_col_sym, a->air_sym,     nnz_sym,  0 );
  //CUDA ps_cuda_memcpy( cdata_, a->csx_sym,     a->apr_sym,     l_side+1, 0 );

  }
  else
  {
  // create symmetric matrix on cpu
  int c, index = 0;
  for( int r=0; r<l_side; ++r )
  {
  for( int i=a->apr[r]; i<a->apr[r+1]; ++i )
  {
  c = a->air[i];
  // if in upper triangle
  if( c >= r )
  {
  a->axr_sym[index++] = a->axr[i];
  }
  }
  }
  }

  // copy symmetric data values and delete cpu memory
  //CUDA ps_cuda_memcpy( cdata_, a->coo_val_sym, a->axr_sym, nnz_sym,  0 );

  a->is_sym_good = 1;
  }

  PS_FPT get( int *air, int *apr, PS_FPT *axr, int row, int col )
  {
  for(int g=apr[row]; g<apr[row+1]; g++)
  {
  if(air[g]==col)
  {
  return axr[g];
  }
  else if(air[g]>col)
  {
  return 0;
  }
  }

  return 0;
  }
*/

/*
  void pixsrc_cuda::inc_chol_factor_cusparse( MATRIX *a )
  {
  if( a->icf_done )
  return;

  CUDA create_sym_form( a );

  // few things we'll need
  int dev                           = a->data_->gpu2use[0];
  CUDA set_device( dev );
  commoninputdata  *cdata_          = a->cdata_;
  cusparseHandle_t *handle          = (cusparseHandle_t*)  cdata_->gpudevices[dev].handle;
  cusparseMatDescr_t *descr_sym     = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr_sym;
  cusparseOperation_t nontrans      = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseSolveAnalysisInfo_t *info = (cusparseSolveAnalysisInfo_t*)a->csrsvi;
  cusparseStatus_t status;


  int l_side = a->nrow;
  int nnz_sym = ( a->nnz - l_side ) / 2 + l_side;

  // perform analysis on MATRIX a if it hasn't been done yet
  pthread_mutex_lock( &a->csrsvi_mutex );

  if( !a->csrsvi_done && a->csrsvi )
  {
  status = cusparseDestroySolveAnalysisInfo( *info );

  CUDA errorcheck( cdata_, (void*)&status, 1, "destroy csrsv info" );
  }

  if( !a->csrsvi_done )
  {
  // create csrsv info object
  cusparseSolveAnalysisInfo_t *info_temp;

  if( !a->csrsvi )
  MEMORY ps_malloc( &info_temp, 1 );
  else
  info_temp = (cusparseSolveAnalysisInfo_t*)a->csrsvi;

  status = cusparseCreateSolveAnalysisInfo( info_temp );
  info = info_temp;
  a->csrsvi = (void*)info_temp;

  CUDA errorcheck( cdata_, (void*)&status, 1, "create csrsv info" );

  // perform analysis
  status = cusparseDcsrsv_analysis( *handle, nontrans,
  l_side, nnz_sym, *descr_sym,
  a->coo_val_sym, a->csx_sym, a->coo_col_sym,
  *info                                       );

  CUDA errorcheck( cdata_, (void*)&status, 1, "analyze csrsv info" );

  cudaDeviceSynchronize();

  a->csrsvi_done = 1;
  }

  pthread_mutex_unlock( &a->csrsvi_mutex );



  // allocate memory for incomplete Cholesky factrization if we haven't yet
  if( !a->icf )
  CUDA ps_cuda_malloc( cdata_, &a->icf, nnz_sym );

  // perform incomplete Cholesky factorization if we haven't yet
  if( !a->icf_done )
  {


  // This would need to be done if we were using cusparse to do
  // the Cholesky decomposition but I just can't get it fully working

  CUDA ps_cuda_memcpy( cdata_, a->icf, a->coo_val_sym, nnz_sym, 2 );

  status = cusparseDcsric0( *handle, nontrans,
  l_side, *descr_sym,
  a->icf, a->csx_sym, a->coo_col_sym,
  *info                                 );

  CUDA errorcheck( cdata_, (void*)&status, 1, "inc Cholesky decomposition" );

  a->nnz_chol = nnz_sym;

  a->icf_done = 1;
  }
  }
*/
/*
  void pixsrc_cuda::analyze_upp_low( MATRIX *a )
  {
  // for analyzigin a Cholesky factorization

  if( a->icf_upper_sva && a->icf_lower_sva )
  return;

  int dev                           = a->data_->gpu2use[0];
  CUDA set_device( dev );
  commoninputdata  *cdata_          = a->cdata_;
  cusparseHandle_t *handle          = (cusparseHandle_t*)  cdata_->gpudevices[dev].handle;
  cusparseMatDescr_t *descrT = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr_tri_upper;
  cusparseOperation_t nontrans      = CUSPARSE_OPERATION_NON_TRANSPOSE;
  cusparseOperation_t trans         = CUSPARSE_OPERATION_TRANSPOSE;
  cusparseStatus_t status;

  int l_side = a->nrow;


  if( !a->icf_upper_sva )
  {
  cusparseSolveAnalysisInfo_t *sv;
  MEMORY ps_malloc( &sv, 1 );

  status = cusparseCreateSolveAnalysisInfo( sv );

  CUDA errorcheck( cdata_, (void*)&status, 1, "analyze csrsv info ic-U - create" );

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_analysis( *handle, nontrans,
  l_side, a->nnz_chol, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *sv                                 );
  #else
  status = cusparseDcsrsv_analysis( *handle, nontrans,
  l_side, a->nnz_chol, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *sv                                 );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "analyze csrsv info ic-U" );

  a->icf_upper_sva = (void*)sv;
  }

  if( !a->icf_lower_sva )
  {
  cusparseSolveAnalysisInfo_t *sv;
  MEMORY ps_malloc( &sv, 1 );

  status = cusparseCreateSolveAnalysisInfo( sv );

  CUDA errorcheck( cdata_, (void*)&status, 1, "analyze csrsv info ic-L - create" );

  #ifdef SINGLE_PRECISION
  status = cusparseScsrsv_analysis( *handle, trans,
  l_side, a->nnz_chol, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *sv                                 );
  #else
  status = cusparseDcsrsv_analysis( *handle, trans,
  l_side, a->nnz_chol, *descrT,
  a->icf_val, a->icf_ptr, a->icf_ind,
  *sv                                 );
  #endif

  CUDA errorcheck( cdata_, (void*)&status, 1, "analyze csrsv info ic-L" );

  a->icf_lower_sva = (void*)sv;
  }
  }
*/
/*
  void pixsrc_cuda::inc_chol_factor( MATRIX *a )
  {

  commoninputdata  *cdata_          = a->cdata_;
  inputdata        *data_           = a->data_;

  #ifdef SINGLE_PRECISION

  PRINTER printerror( data_->print2screenname,
  "double precision required for CHOLMOD support",
  cdata_->print2screenmutex                        );

  #else

  bool usefull = 0;

  //inc_chol_factor_cusparse( a );

  if( a->icf_done )
  return;

  CUDA create_sym_form( a );

  // few things we'll need
//    int dev                           = a->data_->gpu2use[0];
//    CUDA set_device( dev );
//    cusparseHandle_t *handle          = (cusparseHandle_t*)  cdata_->gpudevices[dev].handle;
//    cusparseMatDescr_t *descr_sym     = (cusparseMatDescr_t*)cdata_->gpudevices[dev].descr_sym;
//    cusparseOperation_t nontrans      = CUSPARSE_OPERATION_NON_TRANSPOSE;
//    cusparseSolveAnalysisInfo_t *info = (cusparseSolveAnalysisInfo_t*)a->csrsvi;
//    cusparseStatus_t status;


int l_side = a->nrow;
int nnz_sym = ( a->nnz - l_side ) / 2 + l_side;

// resulting factorization
int fac_nnz, *fac_air, *fac_apr;
PS_FPT *fac_axr;


// CHOLDMOD section
{
bool first_time = ( a->chol_com ) ? 0 : 1;

cholmod_sparse *mat;
cholmod_factor *fac;
cholmod_common *com;

if( first_time )
{
MEMORY ps_malloc( &com, 1 );
a->chol_com = (void*)com;

cholmod_start( com );
#ifdef SINGLE_PRECISION
com->dtype      = CHOLMOD_SINGLE;     // not double

#else

com->dtype      = CHOLMOD_DOUBLE;     // not float
#endif

com->itype      = CHOLMOD_INT;        // not long
com->final_ll   = 1;                  // Cholesky and not LDL'
com->supernodal = CHOLMOD_SIMPLICIAL; // not supernodal


mat = cholmod_allocate_sparse( l_side, l_side, nnz_sym,
1, 1, -1, CHOLMOD_REAL, com );
MEMORY ps_free( (int*)mat->p );
a->chol_mat = (void*)mat;
}
else
{
com = (cholmod_common*)a->chol_com;
mat = (cholmod_sparse*)a->chol_mat;
}

mat->p = (void*)a->apr_sym;
mat->i = (void*)a->air_sym;
mat->x = (void*)a->axr_sym;

if( first_time )
{
fac = cholmod_analyze( mat, com );
a->chol_fac = (void*)fac;
}
else
{
fac = (cholmod_factor*)a->chol_fac;
}

cholmod_factorize( mat, fac, com );

fac_nnz = (PS_SIT)    fac->nzmax;
fac_air = (int*)   fac->i;
fac_apr = (int*)   fac->p;
fac_axr = (PS_FPT*)fac->x;

/-*
// not making sure this works
if( usefull )
{
a->nnz_chol = fac_nnz;
CUDA ps_cuda_malloc( cdata_, &a->coo_col_sym, fac_nnz );
CUDA ps_cuda_malloc( cdata_, &a->icf,         fac_nnz );
CUDA ps_cuda_memcpy( cdata_,  a->coo_col_sym, fac_air, fac_nnz,  0, dev );
CUDA ps_cuda_memcpy( cdata_,  a->icf,         fac_axr, fac_nnz,  0, dev );
CUDA ps_cuda_memcpy( cdata_,  a->csx_sym,     fac_apr, l_side+1, 0, dev );
}
*-/

/-*
fac->i = 0;
fac->x = 0;
fac->p = 0;
mat->i = 0;
mat->x = 0;
mat->p = 0;

cholmod_free_sparse( &mat, com );
cholmod_free_factor( &fac, com );

cholmod_finish     (       &com );
*-/
}

if( !usefull )
{
if( !a->icf_ptr_cpu )
{
MEMORY ps_malloc( &a->icf_ptr_cpu, l_side+1 );
}
if( !a->icf_ind_cpu )
{
MEMORY ps_malloc( &a->icf_ind_cpu, nnz_sym );
}
if( !a->icf_val_cpu )
{
MEMORY ps_malloc( &a->icf_val_cpu, nnz_sym );
}
std::copy( a->apr_sym, a->apr_sym + l_side+1, a->icf_ptr_cpu );
std::copy( a->air_sym, a->air_sym + nnz_sym, a->icf_ind_cpu );

a->selective_copy( a->icf_ind_cpu, a->icf_ptr_cpu, a->icf_val_cpu,
fac_air,    fac_apr,    fac_axr,    l_side  );
a->nnz_chol = a->icf_ptr_cpu[l_side];


// allocate memory for incomplete Cholesky factrization if we haven't yet
if( !a->icf_ptr )
CUDA ps_cuda_malloc( cdata_, &a->icf_ptr, l_side+1 );
if( !a->icf_ind )
CUDA ps_cuda_malloc( cdata_, &a->icf_ind, nnz_sym );
if( !a->icf_val )
CUDA ps_cuda_malloc( cdata_, &a->icf_val, nnz_sym );

// copy factorization to gpu
CUDA ps_cuda_memcpy( cdata_, a->icf_val, a->icf_val_cpu, a->nnz_chol,  0, dev );
CUDA ps_cuda_memcpy( cdata_, a->icf_ind, a->icf_ind_cpu, a->nnz_chol, 0, dev );
CUDA ps_cuda_memcpy( cdata_, a->icf_ptr, a->icf_ptr_cpu, l_side+1,  0, dev );
}
else
{
a->nnz_chol = fac_nnz;

if( !a->icf_ptr )
CUDA ps_cuda_malloc( cdata_, &a->icf_ptr, l_side+1 );
if( !a->icf_ind )
CUDA ps_cuda_malloc( cdata_, &a->icf_ind, fac_nnz );
if( !a->icf_val )
CUDA ps_cuda_malloc( cdata_, &a->icf_val, fac_nnz );

// copy factorization to gpu
CUDA ps_cuda_memcpy( cdata_, a->icf_val, fac_axr, fac_nnz,  0, dev );
CUDA ps_cuda_memcpy( cdata_, a->icf_ind, fac_air, fac_nnz, 0, dev );
CUDA ps_cuda_memcpy( cdata_, a->icf_ptr, fac_apr, l_side+1,  0, dev );
}

a->icf_done = 1;























/-*

// This would need to be done if we were using cusparse to do
// the Cholesky decomposition but I just can't get it fully working

// perform analysis on MATRIX a if it hasn't been done yet
pthread_mutex_lock( &a->csrsvi_mutex );

if( !a->csrsvi_done && a->csrsvi )
{
status = cusparseDestroySolveAnalysisInfo( *info );

CUDA errorcheck( cdata_, (void*)&status, 1, "destroy csrsv info" );
}

if( !a->csrsvi_done )
{
// create csrsv info object
cusparseSolveAnalysisInfo_t *info_temp;

if( !a->csrsvi )
MEMORY ps_malloc( &info_temp, 1 );
else
info_temp = (cusparseSolveAnalysisInfo_t*)a->csrsvi;

status = cusparseCreateSolveAnalysisInfo( info_temp );
info = info_temp;
a->csrsvi = (void*)info_temp;

CUDA errorcheck( cdata_, (void*)&status, 1, "create csrsv info" );

// perform analysis
status = cusparseDcsrsv_analysis( *handle, nontrans,
l_side, nnz_sym, *descr_sym,
a->coo_val_sym, a->csx_sym, a->coo_col_sym,
*info                                       );

CUDA errorcheck( cdata_, (void*)&status, 1, "analyze csrsv info" );

cudaDeviceSynchronize();

a->csrsvi_done = 1;
}

pthread_mutex_unlock( &a->csrsvi_mutex );



// allocate memory for incomplete Cholesky factrization if we haven't yet
if( !a->icf )
CUDA ps_cuda_malloc( cdata_, &a->icf, nnz_sym );

// perform incomplete Cholesky factorization if we haven't yet
if( !a->icf_done )
{


// This would need to be done if we were using cusparse to do
// the Cholesky decomposition but I just can't get it fully working

CUDA ps_cuda_memcpy( cdata_, a->icf, a->coo_val_sym, nnz_sym, 2, dev );

status = cusparseDcsric0( *handle, nontrans,
l_side, *descr_sym,
a->icf, a->csx_sym, a->coo_col_sym,
*info                                 );

CUDA errorcheck( cdata_, (void*)&status, 1, "inc Cholesky decomposition" );

//    }return;{








int *of_od_air, *of_od_apr;
PS_FPT *of_od_axr;
int nnz_of_od = 0;
{
// analysis phase

bool *nz;
int nz_n = ( l_side*l_side - l_side ) / 2 + l_side;
MEMORY ps_malloc( &nz, nz_n );
std::fill( nz, nz + nz_n, 0 );

// flag element as used
for( int r=0; r<l_side; ++r )
for( int i=a->apr_sym[r]; i<a->apr_sym[r+1]; ++i )
{
nz[ r*l_side + a->air_sym[i] - (r*(r+1))/2 ] = 1;
}

// go through algorithm and see which elements will be
// changed from zero to nonzero

// iterate over all rows
for( int r=0; r<l_side; ++r )
{
// iterate over rows under this row
for( int r2=r+1; r2<l_side; ++r2 )
{
// iterate over columns
for( int c2=r2; c2<l_side; ++c2 )
{
// if both are non-zero then we have a possibly new nonzero
if( nz[ r *l_side + r2 - (r *(r +1))/2 ] &&
nz[ r *l_side + c2 - (r *(r +1))/2 ]   )
nz[ r2*l_side + c2 - (r2*(r2+1))/2 ] = 1;
}
}
}

// find total number of nonzeros off-diagonal and off first row
nnz_of_od = 0;
for( int j=0; j< nz_n; ++j )
if( nz[j] )
++nnz_of_od;

// allocate memory for these to be found nonzeros in CSR format
MEMORY ps_malloc( &of_od_apr, l_side+1  );
MEMORY ps_malloc( &of_od_air, nnz_of_od );
MEMORY ps_malloc( &of_od_axr, nnz_of_od );
int r = 0;
int c = 0;
int index = 0;
of_od_apr[0] = of_od_apr[1] = 0;
for( int j=0; j<nz_n; ++j )
{
if( nz[j] )
{
of_od_air[index] = c;
of_od_axr[index] = get(a->air_sym,a->apr_sym,a->axr_sym,r,c);
++of_od_apr[r+1];
++index;
}

++c;
if( c == l_side )
{
++r;
of_od_apr[r+1] = of_od_apr[r];
c = r;
}
}

MEMORY ps_free( nz );
}

int c1, c2;
PS_FPT diag_entry, val1, val2;
for( int r=0; r<l_side; ++r )
{
std::cout << r << std::endl;

// change element on diagonal in this row
diag_entry = std::sqrt( get(of_od_air,of_od_apr,of_od_axr,r,r) );
of_od_axr[of_od_apr[r]] = diag_entry;

// change every non-diagonal element on this row
for(int i=of_od_apr[r]+1; i<of_od_apr[r+1]; ++i )
{
of_od_axr[i] = of_od_axr[i] / diag_entry;
}

// iterate over rows under this row
for( int r2=r+1; r2<l_side; ++r2 )
{
// iterate over columns
for( int i2=of_od_apr[r2]; i2<of_od_apr[r2+1]; ++i2 )
{
// find product of two entries in original row
// that this entry will be reduced by

// two column indices to match
c1 = r2;
c2 = of_od_air[i2];

val1 = get(of_od_air,of_od_apr,of_od_axr,r,c1);
val2 = get(of_od_air,of_od_apr,of_od_axr,r,c2);

of_od_axr[i2] -= val1*val2;
}
}
}

for( int r=0; r<l_side; ++r )
for( int i=a->apr_sym[r]; i<a->apr_sym[r+1]; ++i )
a->axr_sym[i] = get(of_od_air,of_od_apr,of_od_axr,r,a->air_sym[i]);

CUDA ps_cuda_memcpy( cdata_, a->icf, a->axr_sym, nnz_sym, 0, dev );

a->icf_done = 1;

MEMORY ps_free( of_od_air );
MEMORY ps_free( of_od_apr );
MEMORY ps_free( of_od_axr );
}
*-/

#endif

}
*/
 /*
   void pixsrc_cuda::ps_wrap_cusp( MATRIX *a )
   {
   if( a->cusp_csr )
   return;

   a->update_gpu( pmf_csr );

   thrust::device_ptr<PS_SIT>    wrapped_device_Ap( a->csx     );
   thrust::device_ptr<PS_SIT>    wrapped_device_Aj( a->coo_col );
   thrust::device_ptr<PS_FPT> wrapped_device_Ax( a->coo_val );

   typedef typename cusp::array1d_view< thrust::device_ptr<PS_SIT>    > DeviceIndexArrayView;
   typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

   DeviceIndexArrayView row_offsets   ( wrapped_device_Ap, wrapped_device_Ap + a->nrow+1 );
   DeviceIndexArrayView column_indices( wrapped_device_Aj, wrapped_device_Aj + a->nnz    );
   DeviceValueArrayView values        ( wrapped_device_Ax, wrapped_device_Ax + a->nnz    );

   typedef cusp::csr_matrix_view<DeviceIndexArrayView,
   DeviceIndexArrayView,
   DeviceValueArrayView> DeviceView;

   DeviceView *dummy, *a_cusp;
   MEMORY ps_malloc( &dummy, 1 );
   a->cusp_csr_ptr = (void*)dummy;

   a_cusp = new (dummy) DeviceView( a->nrow, a->ncol, a->nnz,
   row_offsets, column_indices, values );
   a->cusp_csr = (void*)a_cusp;
   }

   void pixsrc_cuda::ps_wrap_cusp( VECTOR *a )
   {
   if( a->cusp_wrapper )
   return;

   a->update_gpu();

   thrust::device_ptr<PS_FPT> wrapped_device_vec( a->cuda_vec );

   typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

   DeviceValueArrayView *dummy, *a_cusp;
   MEMORY ps_malloc( &dummy, 1 );
   a->cusp_wrapper_ptr = (void*)dummy;

   a_cusp = new (dummy) DeviceValueArrayView( wrapped_device_vec, wrapped_device_vec + a->size );
   a->cusp_wrapper = (void*)a_cusp;
   }

   void pixsrc_cuda::destroy_cusp_wrapper( MATRIX *a )
   {
   if( !a->cusp_csr )
   return;

   typedef typename cusp::array1d_view< thrust::device_ptr<PS_SIT>    > DeviceIndexArrayView;
   typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

   typedef cusp::csr_matrix_view<DeviceIndexArrayView,
   DeviceIndexArrayView,
   DeviceValueArrayView> DeviceView;

   DeviceView *dv  = (DeviceView*)a->cusp_csr;
   DeviceView *dvp = (DeviceView*)a->cusp_csr_ptr;
   dv->~DeviceView();
   MEMORY ps_free( dvp );
   }

   void pixsrc_cuda::destroy_cusp_wrapper( VECTOR *a )
   {
   if( !a->cusp_wrapper )
   return;

   typedef typename cusp::array1d_view< thrust::device_ptr<PS_FPT> > DeviceValueArrayView;

   DeviceValueArrayView *dv  = (DeviceValueArrayView*)a->cusp_wrapper;
   DeviceValueArrayView *dvp = (DeviceValueArrayView*)a->cusp_wrapper_ptr;
   dv->~DeviceValueArrayView();
   MEMORY ps_free( dvp );
   }
 */
