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



#include "pixsrc_matrix.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_printer.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_cuda.hpp"
#include <umfpack.h>
#include <cmath>
#include <algorithm>
#include <cstdlib> // for exit
#include <iostream>
#include <cstring>
#include <iomanip>

// this function is passed to std::transform
PS_SIT op_decrease (PS_SIT i) {return --i;}

// constructor
pixsrc_matrix::pixsrc_matrix( commoninputdata *cdata__, inputdata *data__,
                              PS_SIT nrow_, PS_SIT ncol_, PS_SIT vecsize_, PS_SIT type_, inputdata *dummy_)
{
    // if fitting in visibility space (intferometric data), then use dense matrices
    if (data__ && data__->is_uvdata)
    {
        vecsize_ = -1;
        type_ = -1;
    }

    // matrix type is really not used right now
    if( type_ == 11 )
        PRINTER printerror( data__->print2screenname,
                            "identity matrix not supported yet",
                            cdata__->print2screenmutex           );

    pm_init( cdata__, data__, nrow_, ncol_, vecsize_, type_ );
}

void pixsrc_matrix::pm_init( commoninputdata *cdata__, inputdata *data__,
                             PS_SIT nrow_, PS_SIT ncol_, PS_SIT vecsize_, PS_SIT type_ )
{
    cdata_ = cdata__;
    data_  = data__;

    ncol = ncol_;
    nrow = nrow_;

    vecsize = vecsize_<0 ? 0 : (PS_SIT)vecsize_;
//    vecsize = 0 ? -1 : vecsize_;

    matrix_type = type_;

    // dense matrix stuff added for use with shapelets
    is_dense = -1==vecsize_ || (data_ && data_->use_shapelets) ? 1 : 0;
    lu_fac_dense = (PS_FPT*)0;
    mat_dense    = (PS_FPT*)0;
    inv_dense    = (PS_FPT*)0;
    ipiv_dense   = 0;
    got_det = got_lu_fac = got_inv = 0;
    if (is_dense)
    {
        pthread_mutex_init (&dense_lu_mutex, NULL);
        pthread_mutex_init (&dense_inv_mutex, NULL);
        nnz = nrow*ncol;
        MEMORY ps_malloc (&mat_dense, (PS_SIT)nnz);
        std::fill (mat_dense, mat_dense+nnz, (PS_FPT)0);
    }
    else
    {
        nnz             = 0;
    }


    scalar          = (PS_FPT)1;
    rcform          = 0;
    ludone          = 0;
    hasbeencompiled = 0;

    apc      = 0;
    apr      = 0;
    aic      = 0;
    air      = 0;
    axc      = (PS_FPT*)0;
    axr      = (PS_FPT*)0;
    icf_ptr_cpu = 0;
    icf_ind_cpu = 0;
    icf_val_cpu = (PS_FPT*)0;
    symbolic = 0;
    numeric  = 0;

    is_sym_good = 0;
    air_sym  = 0;
    apr_sym  = 0;
    axr_sym  = (PS_FPT*)0;

    are_we_using_cuda = ( data_ && data_->numgpu2use ) ? 1 : 0;

#ifdef __USE_PIXSRC_CUDA__

    num_ctv = 5;
    MEMORY ps_malloc( &ctv, (PS_SIT)num_ctv );
    for( PS_SIT j=0; j<num_ctv; ++j )
        ctv[j] = (PS_FPT*)0;

    /*
      csx_sym     = 0;
      coo_col_sym = 0;
      coo_val_sym = 0;
    */

    chol_fac = 0;
    chol_com = 0;
    chol_mat = 0;

    csx         = 0;
    coo_col     = 0;
    coo_val     = (PS_FPT*)0;
    csx_tr      = 0;
    coo_row_tr  = 0;
    coo_val_tr  = (PS_FPT*)0;

    icf_ind  = 0;
    icf_ptr  = 0;
    icf_val  = (PS_FPT*)0;
    icf_done = 0;
    icf_lower_sva = 0;
    icf_upper_sva = 0;

    /*
      csrsvi      = 0;
      csrsvi_done = 0;

      CUDA ps_set_descr( this );
    */

    // cusp preconditioner
    cusp_prec = 0;

    if( are_we_using_cuda )
    {
        CUDA create_stream( this );
    }

#endif

    tiv = 0;
    tjv = 0;
    tkv = (PS_FPT*)0;

    if( matrix_type == 10 /*|| matrix_type == 11*/ )
        vecsize = nnz = ncol_;

    if( !are_we_using_cuda )
    {
        if( vecsize )
        {
            if( !(matrix_type == 10 /*|| matrix_type == 11*/) )
            {
                MEMORY ps_malloc( &tiv, (PS_SIT)vecsize );
                MEMORY ps_malloc( &tjv, (PS_SIT)vecsize );
            }
            MEMORY ps_malloc( &tkv, (PS_SIT)vecsize );
        }

        if( matrix_type == 10 )
            std::fill( tkv, tkv + nnz, (PS_FPT)0 );

        /*
          if( matrix_type == 11 )
          std::fill( tkv, tkv + nnz, 1 );
        */

    }

    up2date_csr = up2date_csc = pms_neither;

    pthread_attr_init          ( &attr                                  );
    pthread_attr_setdetachstate( &attr,         PTHREAD_CREATE_JOINABLE );
    pthread_mutex_init         ( &compilemutex, NULL                    );
    pthread_mutex_init         ( &dissolvemutex, NULL                    );
    pthread_mutex_init         ( &numsymbmutex, NULL                    );
    pthread_mutex_init         ( &rcformmutex,  NULL                    );


#ifdef __USE_PIXSRC_CUDA__
    pthread_mutex_init         ( &pm_mutex,     NULL                    );
    /*
     * pthread_mutex_init         ( &csrsvi_mutex, NULL                    );
     */
#endif


}

// destructor
pixsrc_matrix::~pixsrc_matrix()
{
    clearall();
}

void pixsrc_matrix::init_cpu_vecs()
{
    // t?v vectors hold matrix coordinates and values for sparse matrix
    // after "compiling," they are converted to a compressed format

    if (is_dense)
        return;

    if( vecsize>0 )
    {
        if( !(matrix_type == 10 /*|| matrix_type == 11*/) )
        {
            if( !tiv )
            {
                MEMORY ps_malloc( &tiv, (PS_SIT)vecsize );
            }
            if( !tjv )
            {
                MEMORY ps_malloc( &tjv, (PS_SIT)vecsize );
            }
        }
        if( !tkv )
        {
            MEMORY ps_malloc( &tkv, (PS_SIT)vecsize );

            if( matrix_type == 10 )
                std::fill( tkv, tkv + nnz, (PS_FPT)0 );
        }
    }
}

void pixsrc_matrix::set( PS_SIT r, double val )
{
    // no mutex locking is performed.
    // This assumes only one thread tries to write to a matrix at a time

    // this function should only be used by diagonal matrices
    // but not those proportional to identity

    // create t?v vectors
    this->init_cpu_vecs();

    // cannot set matrix value after thas been compiled
    if( (hasbeencompiled || rcform) && !is_dense )
    {
        PRINTER printerror( data_->print2screenname,
                            "cannot change matrix after it has been compiled (used)",
                            cdata_->print2screenmutex                                 );
    }

    if (is_dense)
        mat_dense[r*nrow+r] = (PS_FPT)val / scalar;
    else
        tkv[r] = (PS_FPT)val / scalar;
}

void pixsrc_matrix::set(PS_SIT y, PS_SIT x, double val)
{
    // no mutex locking is performed.
    // This assumes only one thread tries to write to a matrix at a time

    // this function is NOT meant for diagonal matrices

    // create t?v vectors
    this->init_cpu_vecs();

    // cannot set matrix value after thas been compiled
    if( (hasbeencompiled || rcform) && !is_dense )
    {
        PRINTER printerror( data_->print2screenname,
                            "cannot change matrix after it has been compiled (used)",
                            cdata_->print2screenmutex                                 );
    }

    // check if we need to make t?v vectors larger
    if( vecsize == nnz )
    {
        this->resize( &tiv, &tjv, &tkv, nnz, &vecsize );
    }

    if (is_dense)
    {
        mat_dense[x*nrow+y] = (PS_FPT)val / scalar;
    }
    else
    {
        tiv[nnz] = y;
        tjv[nnz] = x;
        tkv[nnz] = (PS_FPT)val / scalar;
        ++nnz;
    }
}

void pixsrc_matrix::resize( PS_SIT **ptr1, PS_SIT **ptr2, PS_FPT **ptr3, PS_SIT oldsize, PS_SIT *newsize )
{
    // make t?v vectors larger, so that more entries can be stored

    if (is_dense)
        return;

    if( !*newsize )
        ++*newsize;
    *newsize*=2;

    PS_SIT    *dum1,*dum2;
    PS_FPT *dum3;
    MEMORY ps_malloc( &dum1, (PS_SIT)(*newsize) );
    MEMORY ps_malloc( &dum2, (PS_SIT)(*newsize) );
    MEMORY ps_malloc( &dum3, (PS_SIT)(*newsize) );
    std::copy( *ptr1, *ptr1+oldsize, dum1 );
    std::copy( *ptr2, *ptr2+oldsize, dum2 );
    std::copy( *ptr3, *ptr3+oldsize, dum3 );
    MEMORY ps_free( *ptr1 );
    MEMORY ps_free( *ptr2 );
    MEMORY ps_free( *ptr3 );
    *ptr1 = dum1;
    *ptr2 = dum2;
    *ptr3 = dum3;
}

void pixsrc_matrix::update_gpu( pm_format format )
{
    // make sure the matrix on the GPU is up to date
    // and has the right compression (row- versus column- compressed)

    if (is_dense)
    {
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: CUDA not supported yet",
                            cdata_->print2screenmutex);
        return;
    }


#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU to be free
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

    switch( format )
    {

        // if we want compressed rows
    case pmf_csr:
    {
        // compress
        if( !rcform )
            creatercform();

        if( up2date_csr == pms_both || up2date_csr == pms_gpu )
        {
            // if there's nothing to do, return
            return;
        }
        else if( up2date_csr == pms_cpu )
        {

#ifdef __USE_PIXSRC_CUDA__
            // copy to GPU
            CUDA sync_matrix( this, format, ccd_cpu2gpu );
#endif

        }
        else
        {
            // no data has been stored in matrix
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix: up2date gpu csr",
                                cdata_->print2screenmutex );
        }

        break;
    }

    // if we want compressed columns
    case pmf_csc:
    {
        // compress
        if( !hasbeencompiled )
            compile();

        if( up2date_csc == pms_both || up2date_csc == pms_gpu )
        {
            // if there's nothing to do, return
            return;
        }
        else if( up2date_csc == pms_cpu )
        {

#ifdef __USE_PIXSRC_CUDA__
            // copy to GPU
            CUDA sync_matrix( this, format, ccd_cpu2gpu );
#endif

        }
        else
        {
            // no data has been stored in matrix
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix: up2date gpu csc",
                                cdata_->print2screenmutex );
        }

        break;
    }

    default:
    {
        // format not recognized
        PRINTER printerror( data_->print2screenname,
                            "pixsrc_matrix: up2date gpu",
                            cdata_->print2screenmutex );
        break;
    }

    }

#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU to finish
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

}

void pixsrc_matrix::update_cpu( pm_format format )
{
    if (is_dense)
        return;


#ifdef __USE_PIXSRC_CUDA__
    if( are_we_using_cuda )
        // wait for GPU to be free
        CUDA wait_for_stream( this );
#endif

    switch( format )
    {

        // if we want compressed rows
    case pmf_csr:
    {
        // compress
        if( !rcform )
            creatercform();

        if( up2date_csr == pms_both || up2date_csr == pms_cpu )
        {
            // if there's nothing to do, return
            return;
        }
        else if( up2date_csr == pms_gpu )
        {

#ifdef __USE_PIXSRC_CUDA__
            // copy to GPU
            CUDA sync_matrix( this, format, ccd_gpu2cpu );
#endif

        }
        else
        {
            // no data has been stored in matrix
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix: up2date cpu csr",
                                cdata_->print2screenmutex );
        }

        break;
    }

    // if we want compressed columns
    case pmf_csc:
    {
        // compress
        if( !hasbeencompiled )
            compile();

        if( up2date_csc == pms_both || up2date_csc == pms_cpu )
        {
            // if there's nothing to do, return
            return;
        }
        else if( up2date_csc == pms_gpu )
        {

#ifdef __USE_PIXSRC_CUDA__
            // copy to GPU
            CUDA sync_matrix( this, format, ccd_gpu2cpu );
#endif

        }
        else
        {
            // no data has been stored in matrix
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_matrix: up2date cpu csc",
                                cdata_->print2screenmutex );
        }

        break;
    }

    default:
    {
        // format not recognized
        PRINTER printerror( data_->print2screenname,
                            "pixsrc_matrix: up2date cpu",
                            cdata_->print2screenmutex );
        break;
    }

    }

#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU to finish
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

}

void pixsrc_matrix::dissolve_scalar()
{
    // set scalar to 1
    // This multiplies all elements in matrix by scalar first.

    // scalar is already 1
    if( scalar == (PS_FPT)1 )
        return;

    // lock out any other threads trying to edit matrix
    pthread_mutex_lock (&dissolvemutex);
    if (scalar==(PS_FPT)1)
    {
        pthread_mutex_unlock (&dissolvemutex);
        return;
    }

    // dense matrix code (for shapelets)
    if (is_dense)
    {
        std::transform (mat_dense, mat_dense+nnz, mat_dense,
                        std::bind1st(std::multiplies<PS_FPT>(),scalar));
        scalar = (PS_FPT)1;
        pthread_mutex_unlock (&dissolvemutex);
        return;
    }

    // if matrix has been column compressed
    if( hasbeencompiled )
    {
        // find out whether the CPU or GPU has the matrix up to date
        // and dissolve scalar on that device

        if( up2date_csc == pms_cpu )
        {
            std::transform (axc, axc+nnz, axc,
                            std::bind1st(std::multiplies<PS_FPT>(),scalar));
        }
        else if( up2date_csc == pms_gpu )
        {

#ifdef __USE_PIXSRC_CUDA
            CUDA dissolve_scalar( this, pmf_csc );
#endif

        }
        else if( up2date_csc == pms_both )
        {
            std::transform (axc, axc+nnz, axc,
                            std::bind1st(std::multiplies<PS_FPT>(),scalar));

#ifdef __USE_PIXSRC_CUDA
            CUDA dissolve_scalar( this, pmf_csc );
#endif

        }
    }

    // if matrix has been row compressed
    if( rcform )
    {
        // find out whether the CPU or GPU has the matrix up to date
        // and dissolve scalar on that device

        if( up2date_csr == pms_cpu )
        {
            std::transform (axr, axr+nnz, axr,
                            std::bind1st(std::multiplies<PS_FPT>(),scalar));
        }
        else if( up2date_csr == pms_gpu )
        {

#ifdef __USE_PIXSRC_CUDA
            CUDA dissolve_scalar( this, pmf_csr );
#endif

        }
        else if( up2date_csr == pms_both )
        {
            std::transform (axr, axr+nnz, axr,
                            std::bind1st(std::multiplies<PS_FPT>(),scalar));

#ifdef __USE_PIXSRC_CUDA
            CUDA dissolve_scalar( this, pmf_csr );
#endif

        }
    }

    // set scalar to 1 and unlock mutex lock
    scalar = (PS_FPT)1;
    pthread_mutex_unlock (&dissolvemutex);
}

void pixsrc_matrix::compile()
{
    // column compress matrix

    if (is_dense)
        return;

    // cdata_ is NULL unless this function is called from pixsrc_cuda

    pthread_mutex_lock( &compilemutex );

    // if there's nothing to do, return
    if( hasbeencompiled )
    {
        pthread_mutex_unlock( &compilemutex );
        return;
    }

    // if matrix is diagonal
    if( matrix_type == 10 /*|| matrix_type == 11*/ )
    {
        // first up, we can easily create CSR and CSC formats
        // for a diagonal matrix

        // allocate compressed format vectors
        if( !apc )
            MEMORY ps_malloc( &apc, (PS_SIT)ncol+1 );

        if( !aic )
            MEMORY ps_malloc( &aic, (PS_SIT)nnz    );

        if(  axc && axc != tkv )
            MEMORY ps_free  (  axc         );

        // matrix values copy over
        axc = tkv;
        tkv = (PS_FPT*)0;

        // simple compression for diagonal matrix
        PS_SIT *ptr;
        ptr = apc;
        for( PS_SIT c=0; c<ncol+1; ++c )
        {
            *ptr = c;
            ++ptr;
        }
        ptr = aic;
        for( PS_SIT r=0; r<nrow; ++r )
        {
            *ptr = r;
            ++ptr;
        }

        up2date_csc = pms_cpu;

#ifdef __USE_PIXSRC_CUDA__

        // copy compressed matrix to the GPU
        if( are_we_using_cuda )
            CUDA sync_matrix( this, pmf_csc, ccd_cpu2gpu );

#endif

        // creatercform now (since it's simple)
        // for diagonal matrices, CSR and CSC forms are identical
        pthread_mutex_lock(&rcformmutex);
        if( !rcform )
        {
            apr = apc;
            air = aic;
            axr = axc;

            up2date_csr = pms_cpu;

#ifdef __USE_PIXSRC_CUDA__

            if( are_we_using_cuda )
            {
                // get rid of  memory if already allocated
                if( coo_col )
                    CUDA ps_cuda_free( coo_col );

                if( coo_val )
                    CUDA ps_cuda_free( coo_val );

                if( csx )
                    CUDA ps_cuda_free( csx     );

                // create soft links
                coo_col = coo_row_tr;
                coo_val = coo_val_tr;
                csx     = csx_tr;

                // update
                up2date_csr = pms_both;
            }

#endif

            rcform = 1;
        }
        pthread_mutex_unlock(&rcformmutex);

    }

#ifdef __USE_PIXSRC_CUDA__

    else if( are_we_using_cuda )
    {
        // second up, we let the gpu do the conversion
        // must create CSR format first before we can make CSC format

        CUDA compile( this );
    }

#endif

    else
    {
        // third up, we let umfpack do the conversion

        // allocate vectors
        if( !apc )
            MEMORY ps_malloc( &apc, (PS_SIT)ncol+1 );
        if( !aic )
            MEMORY ps_malloc( &aic, (PS_SIT)nnz    );
        if( !axc )
            MEMORY ps_malloc( &axc, (PS_SIT)nnz    );

#ifndef DOUBLE_PRECISION
        // UMFPACK only works with double precision
        PRINTER printerror( data_->print2screenname,
                            "double precision required for UMFPACK: coo2csc",
                            cdata_->print2screenmutex                );

#else

        // call UMFPACK
#ifndef INTEGER_TYPE_I
        umfpack_dl_triplet_to_col( nrow, ncol, nnz,
                                   tiv, tjv, tkv,
                                   apc, aic, axc, (PS_SIT*)NULL );
#else
        umfpack_di_triplet_to_col( nrow, ncol, nnz,
                                   tiv, tjv, tkv,
                                   apc, aic, axc, (PS_SIT*)NULL );
#endif

#endif

        // update
        up2date_csc = pms_cpu;
        hasbeencompiled = 1;

    }

    // don't need t?v vectors anymore
    clearvecs();
    // unlock mutex
    pthread_mutex_unlock( &compilemutex );
}

double pixsrc_matrix::get( PS_SIT row, PS_SIT col )
{
    // return value of matrix at row, col

    // dense stuff (for shapelets)
    if (is_dense)
        return scalar*mat_dense[col*nrow+row];

    // compile if we haven't yet
    if( !rcform && !hasbeencompiled )
    {
        // prefer row compression for CUDA
        if( are_we_using_cuda )
            this->creatercform();
        // prefer column compression for UMFPACK
        else
            this->compile();
    }

    // update
    pm_format format = pmf_csr;

    // figure out the best way to find entry in matrix
    // this depends on compression format (column versus row) and on
    // whether the CPU or GPU is up to date
    if( rcform && ( up2date_csr == pms_both || up2date_csr == pms_cpu ) )
    {
        format = pmf_csr;
    }
    else if( hasbeencompiled && ( up2date_csc == pms_both || up2date_csc == pms_cpu ) )
    {
        format = pmf_csc;
    }
    else if( rcform && up2date_csr == pms_gpu )
    {
        format = pmf_csr;
        // update CPU to get entry
        this->update_cpu( format );
    }
    else if( hasbeencompiled && up2date_csc == pms_gpu )
    {
        format = pmf_csc;
        // update CPU to get entry
        this->update_cpu( format );
    }
    else
    {
        PRINTER printerror( data_->print2screenname,
                            " pixsrc_matrix get",
                            cdata_->print2screenmutex );
    }

    // now that the CPU is up to date and we have figured out the compression format,
    // search the matrix and get the entry
    if( format == pmf_csr )
    {
        for(PS_SIT g=apr[row]; g<apr[row+1]; g++)
        {
            if(air[g]==col)
            {
                return (double)( axr[g]*scalar );
            }
            else if(air[g]>col)
                return (double)0.0;
        }
        return (double)0.0;
    }
    else
    {
        for(PS_SIT g=apc[col]; g<apc[col+1]; g++)
        {
            if(aic[g]==row)
            {
                return (double)( axc[g]*scalar );
            }
            else if(aic[g]>row)
                return (double)0.0;
        }
        return (double)0.0;
    }
}

PS_SIT pixsrc_matrix::get_nnz()
{
    // number of non-zero entries
    return nnz;
}

void pixsrc_matrix::remove_last_row_col()
{
    if (!is_dense)
        PRINTER printerror( data_->print2screenname,
                            "removing last row and col for sparse matrix attempted",
                            cdata_->print2screenmutex                );

    // chop off last row
    PS_SIT copysize = (nrow-1)*sizeof(mat_dense[0]);
    for (PS_SIT n=0; n<ncol; ++n)
        std::memmove (mat_dense+n*(nrow-1), mat_dense+n*nrow, copysize);

    nrow -= 1;
    ncol -= 1;
    nnz = nrow*ncol;
}

void pixsrc_matrix::inplace_transpose()
{
    // transpose matrix
    if (!is_dense || ncol!=nrow)
        PRINTER printerror( data_->print2screenname,
                            "in place transposition of non-square or sparse matrix attempted",
                            cdata_->print2screenmutex                );

    for (PS_SIT m=0; m<nrow-1; ++m)
        for(PS_SIT n=m+1; n<ncol; ++n)
            OPERA swap (&mat_dense[m*nrow+n], &mat_dense[n*nrow+m]);
}

void pixsrc_matrix::mult( double val )
{
    // this could be accounted for, but I don't need it right now
    if (is_dense && got_lu_fac)
        PRINTER printerror( data_->print2screenname,
                            "can't multiply matrix after LU factorization",
                            cdata_->print2screenmutex                );

    scalar *= (PS_FPT)val;
}

void pixsrc_matrix::creatercform()
{
    // create row-compressed matrix

    // dense stuff (for shapelets)
    if (is_dense)
        return;

    // cdata_ is NULL unless this function is called from pixsrc_cuda
    // or from compile() if it was called from within pixsrc_cuda

    // for a diagonal matrix, CSR and CSC formatting is done all at once
    if( matrix_type == 10 /*|| matrix_type == 11*/ )
    {
        // column compressed function will handle this
        compile();
        return;
    }

    // lock mutex
    pthread_mutex_lock(&rcformmutex);

    // if there's nothing to do, return
    if( rcform )
    {
        pthread_mutex_unlock(&rcformmutex);
        return;
    }

#ifdef __USE_PIXSRC_CUDA__

    if( are_we_using_cuda )
    {
        // CUDA will handle this
        CUDA creatercform( this );
    }

    else

#endif

    {
        // update matrix on CPU
        update_cpu( pmf_csc );

        // allocate vectors
        if( !apr )
            MEMORY ps_malloc( &apr, (PS_SIT)nrow+1 );
        if( !air )
            MEMORY ps_malloc( &air, (PS_SIT)nnz    );
        if ( !axr )
            MEMORY ps_malloc( &axr, (PS_SIT)nnz    );

#ifndef DOUBLE_PRECISION
        PRINTER printerror( data_->print2screenname,
                            "double precision required for UMFPACK: csc2csr",
                            cdata_->print2screenmutex                );

#else

        // let UMFPACK do this
#ifndef INTEGER_TYPE_I
        umfpack_dl_transpose( nrow, ncol, apc, aic, axc,
                              (PS_SIT*)NULL, (PS_SIT*)NULL,
                              apr, air, axr              );
#else
        umfpack_di_transpose( nrow, ncol, apc, aic, axc,
                              (PS_SIT*)NULL, (PS_SIT*)NULL,
                              apr, air, axr              );
#endif

#endif

        // update
        up2date_csr = pms_cpu;
        rcform = 1;
    }

    // unlock mutex
    pthread_mutex_unlock(&rcformmutex);
}

// structure for multhreading matrix-matrix multiplication
struct multmatrixmatrixstruct
{
    bool tr1, tr2;
    PS_SIT lowercol, uppercol;
    PS_SIT    *tiv;
    PS_SIT    *tjv;
    PS_FPT *tkv;
    PS_SIT nnz, vecsize;
    pixsrc_matrix *a, *b;
};

void* pixsrc_matrix::multmatrixmatrix( void *args )
{
    // this function is passed to pthreads for multithreading a
    // matrix-matrix multiplication

    multmatrixmatrixstruct *mmmstruct =  (multmatrixmatrixstruct*)args;
    MATRIX *a = mmmstruct->a;
    MATRIX *b = mmmstruct->b;

    if( !mmmstruct->tr1 && !mmmstruct->tr2 )
    {
        for(PS_SIT ccb=mmmstruct->lowercol; ccb<mmmstruct->uppercol; ccb++)
        {
            for(PS_SIT rr=0; rr<a->nrow; rr++)
            {
                PS_SIT startat = b->apc[ccb];
                PS_FPT val = 0;
                PS_SIT index1=a->apr[rr];
                while(1)
                {
                loop3ff:
                    if(!(index1<(PS_SIT)a->apr[rr+1]))
                        break;
                    for(PS_SIT index2=startat; index2<(PS_SIT)b->apc[ccb+1]; index2++)
                    {
                        if(a->air[index1] == b->aic[index2])
                        {
                            val += a->axr[index1]*b->axc[index2];
                            startat=index2+1;
                            index1++;
                            goto loop3ff;
                        }
                        if(b->aic[index2] > a->air[index1])
                        {
                            startat=index2-1;
                            if(startat < (PS_SIT)b->apc[ccb])
                                startat = b->apc[ccb];
                            index1++;
                            goto loop3ff;
                        }
                    }
                    break;
                }
                if( val )
                {
                    if(mmmstruct->vecsize==mmmstruct->nnz)
                    {
                        // no particular reason for calling a's method .. just need an object
                        a->resize( &(mmmstruct->tiv),&(mmmstruct->tjv),&(mmmstruct->tkv),
                                   mmmstruct->nnz,&(mmmstruct->vecsize) );
                    }

                    mmmstruct->tiv[mmmstruct->nnz] = rr;
                    mmmstruct->tjv[mmmstruct->nnz] = ccb;
                    mmmstruct->tkv[mmmstruct->nnz] = val;
                    mmmstruct->nnz++;
                }
            }
        }
    }
    else if (mmmstruct->tr1 && !mmmstruct->tr2)
    {
        for(PS_SIT ccb=mmmstruct->lowercol; ccb<mmmstruct->uppercol; ccb++)
        {
            for(PS_SIT rr=0; rr<a->ncol; rr++)
            {
                PS_SIT startat = b->apc[ccb];
                PS_FPT val = 0;
                PS_SIT index1=a->apc[rr];
                while(1)
                {
                loop3tf:
                    if(!(index1<(PS_SIT)a->apc[rr+1]))
                        break;
                    for(PS_SIT index2=startat; index2<(PS_SIT)b->apc[ccb+1]; index2++)
                    {
                        if(a->aic[index1] == b->aic[index2])
                        {
                            val += a->axc[index1]*b->axc[index2];
                            startat=index2+1;
                            index1++;
                            goto loop3tf;
                        }
                        else if(b->aic[index2] > a->aic[index1])
                        {
                            startat=index2-1;
                            if(startat < (PS_SIT)b->apc[ccb])
                                startat = b->apc[ccb];
                            index1++;
                            goto loop3tf;
                        }
                    }
                    break;
                }
                if( val )
                {
                    if(mmmstruct->vecsize==mmmstruct->nnz)
                    {
                        // no particular reason for calling a's method .. just need an object
                        a->resize( &(mmmstruct->tiv),&(mmmstruct->tjv),&(mmmstruct->tkv),
                                   mmmstruct->nnz,&(mmmstruct->vecsize) );
                    }

                    mmmstruct->tiv[mmmstruct->nnz] = rr;
                    mmmstruct->tjv[mmmstruct->nnz] = ccb;
                    mmmstruct->tkv[mmmstruct->nnz] = val;
                    mmmstruct->nnz++;
                }
            }
        }
    }
    else
    {} // do not need this yet

    return NULL;
}

void pixsrc_matrix::submultadd (pixsrc_matrix *b, pixsrc_matrix *c, bool tr1, bool tr2,
                                PS_SIT am, PS_SIT an, PS_SIT bm, PS_SIT bn,
                                PS_SIT cm, PS_SIT cn, PS_SIT cstartm, PS_SIT cstartn, PS_SIT numprocesses)
{
    //

    if (!is_dense)
        PRINTER printerror (data_->print2screenname,
                            "cannot do submatrix multiplication for sparse matrix",
                            cdata_->print2screenmutex);

    char t = 'T';
    char n = 'N';
    double one = 1.0;
    double *cptr = c->mat_dense + cstartn*c->nrow + cstartm;
    PS_SIT tcnrow = (PS_SIT)c->nrow;

    // let BLAS do the multiplication
    EXTERNAL ps_dgemm_ (tr1 ? &t : &n, tr2 ? &t : &n,
                        &cm, &cn, tr1 ? &am : &an,
                        &one, mat_dense, &am,
                        b->mat_dense, &bm,
                        &one, cptr, &tcnrow);
    c->scalar = scalar*b->scalar;
}

void pixsrc_matrix::mult( pixsrc_matrix *b, pixsrc_matrix *c,
                          bool tr1, bool tr2, PS_SIT numprocesses )
{
    // matrix-matrix multiplication

    // dense stuff (for shapelets)
    if (is_dense)
    {
        char t = 'T';
        char n = 'N';
        double one = 1.0;
        double zer = 0.0;
        PS_SIT tbnrow = (PS_SIT)b->nrow;
        PS_SIT tcnrow = (PS_SIT)c->nrow;
        PS_SIT tcncol = (PS_SIT)c->ncol;
        PS_SIT tnrow = (PS_SIT)nrow;
        PS_SIT tncol = (PS_SIT)ncol;

        // let BLAS do the multiplication
        EXTERNAL ps_dgemm_ (tr1 ? &t : &n, tr2 ? &t : &n,
                            &tcnrow, &tcncol, tr1 ? &tnrow : &tncol,
                            &one, mat_dense, &tnrow,
                            b->mat_dense, &tbnrow,
                            &zer, c->mat_dense, &tcnrow);
        c->scalar = scalar*b->scalar;
        return;
    }



#ifdef __USE_PIXSRC_CUDA__

    if( are_we_using_cuda )
    {
        // let CUDA do the multiplication
        CUDA mult( this, b, c, tr1, tr2 );

        return;
    }

#endif

    // update data on the CPU
    if( !tr1 && !tr2 )
    {
        this->update_cpu( pmf_csr );
        b   ->update_cpu( pmf_csc );
    }
    else if( tr1 && !tr2 )
    {
        this->update_cpu( pmf_csc );
        b   ->update_cpu( pmf_csc );
    }
    c->init_cpu_vecs();
    c->scalar = (double)1.0;

    // could just have the multilication threaded every time but this saves
    // the operations of copying 3 possibly large arrays in the case that
    // there is only one processor.
    if(numprocesses==1)
    {
        if(!tr1 && !tr2)
        {
            for(PS_SIT ccb=0; ccb<b->ncol; ccb++)
            {
                for(PS_SIT rr=0; rr<nrow; rr++)
                {
                    PS_SIT startat = b->apc[ccb];
                    PS_FPT val = 0;
                    PS_SIT index1=apr[rr];
                    while(1)
                    {
                    loop3ff:
                        if(!(index1<(PS_SIT)apr[rr+1]))
                            break;
                        for(PS_SIT index2=startat; index2<(PS_SIT)b->apc[ccb+1]; index2++)
                        {
                            if(air[index1] == b->aic[index2])
                            {
                                val += axr[index1]*b->axc[index2];
                                startat=index2+1;
                                index1++;
                                goto loop3ff;
                            }
                            if(b->aic[index2] > air[index1])
                            {
                                startat=index2-1;
                                if(startat < (PS_SIT)b->apc[ccb])
                                    startat = b->apc[ccb];
                                index1++;
                                goto loop3ff;
                            }
                        }
                        break;
                    }
                    if( val )
                        c->set(rr,ccb,(double)val);
                }
            }
        }
        else if(tr1 && !tr2)
        {
            for(PS_SIT ccb=0; ccb<b->ncol; ccb++)
            {
                for(PS_SIT rr=0; rr<ncol; rr++)
                {
                    PS_SIT startat = b->apc[ccb];
                    PS_FPT val = 0;
                    PS_SIT index1=apc[rr];
                    while(1)
                    {
                    loop3tf:
                        if(!(index1<(PS_SIT)apc[rr+1]))
                            break;
                        for(PS_SIT index2=startat; index2<(PS_SIT)b->apc[ccb+1]; index2++)
                        {
                            if(aic[index1] == b->aic[index2])
                            {
                                val += axc[index1]*b->axc[index2];
                                startat=index2+1;
                                index1++;
                                goto loop3tf;
                            }
                            else if(b->aic[index2] > aic[index1])
                            {
                                startat=index2-1;
                                if(startat < (PS_SIT)b->apc[ccb])
                                    startat = b->apc[ccb];
                                index1++;
                                goto loop3tf;
                            }
                        }
                        break;
                    }
                    if( val )
                        c->set(rr,ccb,(double)val);
                }
            }
        }
        else
        {} // not need yet
    }
    else // if numprocesses != 1
    {
        // multi-threaded code

        pthread_t             threadshere[numprocesses];
        multmatrixmatrixstruct mmmstructs[numprocesses];

        PS_SIT interval  = b->ncol/numprocesses;
        PS_SIT startsize = c->vecsize/numprocesses;
        for(PS_SIT proc=0; proc<numprocesses; proc++)
        {
            mmmstructs[proc].tr1 = tr1;
            mmmstructs[proc].tr2 = tr2;
            mmmstructs[proc].nnz    =0;
            mmmstructs[proc].vecsize=startsize;
            MEMORY ps_malloc( &(mmmstructs[proc].tiv), startsize );
            MEMORY ps_malloc( &(mmmstructs[proc].tjv), startsize );
            MEMORY ps_malloc( &(mmmstructs[proc].tkv), startsize );
            mmmstructs[proc].a=this;
            mmmstructs[proc].b=b;
            mmmstructs[proc].lowercol= proc   *interval;
            mmmstructs[proc].uppercol=(proc+1)*interval;
            if(proc==numprocesses-1)
                mmmstructs[proc].uppercol = b->ncol;

            pthread_create(&threadshere[proc] ,&attr ,multmatrixmatrix,&mmmstructs[proc]);
        }
        for(PS_SIT proc=0; proc<numprocesses; proc++)
            pthread_join(threadshere[proc],NULL);

        // copy results from threads into matrix c
        PS_SIT totalsize = 0;
        for(PS_SIT proc=0; proc<numprocesses; proc++)
            totalsize += mmmstructs[proc].nnz;

        if( totalsize > c->vecsize )
        {
            c->vecsize = c->nnz = totalsize;
            MEMORY ps_free  (   c->tiv             );
            MEMORY ps_free  (   c->tjv             );
            MEMORY ps_free  (   c->tkv             );
            MEMORY ps_malloc( &(c->tiv), (PS_SIT)totalsize );
            MEMORY ps_malloc( &(c->tjv), (PS_SIT)totalsize );
            MEMORY ps_malloc( &(c->tkv), (PS_SIT)totalsize );
        }

        PS_SIT start4copy=0;
        c->nnz = totalsize;
        for(PS_SIT proc=0; proc<numprocesses; proc++)
        {
            std::copy(mmmstructs[proc].tiv,
                      mmmstructs[proc].tiv + mmmstructs[proc].nnz,
                      c->tiv + start4copy);
            std::copy(mmmstructs[proc].tjv,
                      mmmstructs[proc].tjv + mmmstructs[proc].nnz,
                      c->tjv + start4copy);
            std::copy(mmmstructs[proc].tkv,
                      mmmstructs[proc].tkv + mmmstructs[proc].nnz,
                      c->tkv + start4copy);

            MEMORY ps_free( mmmstructs[proc].tiv );
            MEMORY ps_free( mmmstructs[proc].tjv );
            MEMORY ps_free( mmmstructs[proc].tkv );

            start4copy+=mmmstructs[proc].nnz;
        }
    }

    c->mult( (double)( scalar * b->scalar ) );
    // no need to play with up2date_* because they're
    // initialized to -1 (no CSC or CSR)
}

// structure for multhreading matrix-vector multiplication
struct multmatrixvectorstruct
{
    bool tr;
    PS_SIT lowercol, uppercol;
    PS_FPT *subvec;
    pixsrc_matrix *a;
    pixsrc_vector *b;
    pixsrc_vector *c;
};

void* pixsrc_matrix::multmatrixvector(void *args)
{
    // this function is passed to pthreads for multithreading a
    // matrix-vector multiplication

    multmatrixvectorstruct *mmvstruct =  (multmatrixvectorstruct*)args;
    MATRIX *a = mmvstruct->a;
    VECTOR *b = mmvstruct->b;
    VECTOR *c = mmvstruct->c;

    if( !mmvstruct->tr )
    {
        PS_FPT *subvec = mmvstruct->subvec;

        for( PS_SIT cc=mmvstruct->lowercol; cc<mmvstruct->uppercol; ++cc )
            for(PS_SIT index=a->apc[cc]; index<a->apc[cc+1]; index++)
            {
                subvec[a->aic[index]] = subvec[a->aic[index]] +
                    a->axc[index] * (PS_FPT)b->get(cc);
            }
    }
    else
    {
        PS_FPT val;
        for( PS_SIT cc=mmvstruct->lowercol; cc<mmvstruct->uppercol; ++cc )
        {
            val=0;
            for(PS_SIT index=a->apc[cc]; index<a->apc[cc+1]; index++)
                val += a->axc[index]*(PS_FPT)b->get(a->aic[index]);
            c->set(cc, (double)val);
        }
    }

    return NULL;
}

void pixsrc_matrix::inv_dense_mult (pixsrc_vector *b, pixsrc_vector *c, PS_SIT numprocesses)
{
    if (!is_dense)
        PRINTER printerror( data_->print2screenname,
                            "sparse inverse multiplication not supported",
                            cdata_->print2screenmutex           );

    char tr = 'N';
    double one = 1.0;
    double zer = 0.0;
    PS_SIT onei = 1;
    PS_SIT tnrow = (PS_SIT)nrow;
    PS_SIT tncol = (PS_SIT)ncol;

    // let BLAS do the multiplication
    EXTERNAL ps_dgemv_ (&tr, &tnrow, &tncol, &one, inv_dense, &tnrow,
                        b->vec, &onei, &zer, c->vec, &onei);
    c->scalar = scalar*b->scalar;
}

void pixsrc_matrix::inv_matrix_mult (pixsrc_matrix *b, pixsrc_matrix *c, PS_SIT numprocesses)
{
    if (!is_dense)
        PRINTER printerror( data_->print2screenname,
                            "sparse inverse matrix multiplication not supported",
                            cdata_->print2screenmutex           );

    // update CPU
    //this->update_cpu();
    //b->update_cpu();

    // dense stuff (for shapelets)
    dense_matrix_stuff (b,c,0,0,0,5);
}

void pixsrc_matrix::mult( pixsrc_vector *b, pixsrc_vector *c,
                          bool tr1, PS_SIT numprocesses)
{
    // vector-vector multiplication

    // dense stuff (for shapelets)
    if (is_dense)
    {
        char tr = tr1 ? 'T' : 'N';
        double one = 1.0;
        double zer = 0.0;
        PS_SIT onei = 1;
        PS_SIT tnrow = (PS_SIT)nrow;
        PS_SIT tncol = (PS_SIT)ncol;

        // let BLAS do the multiplication
        EXTERNAL ps_dgemv_ (&tr, &tnrow, &tncol, &one, mat_dense, &tnrow,
                            b->vec, &onei, &zer, c->vec, &onei);
        c->scalar = scalar*b->scalar;
        return;
    }

#ifdef __USE_PIXSRC_CUDA__

    // if using CUDA
    if( are_we_using_cuda )
    {
        // let CUDA do the multiplication
        CUDA mult( this, b, c, tr1 );

        return;
    }

#endif

    // update data on the CPU
    this->update_cpu( pmf_csc );
    b->update_cpu();
    c->init_cpu_vec();
    c->set_scalar( (double)1.0 );

    // ERROR .. can't get multithreading to work .. so its disabled
    // if only one thread
    if( 1 || numprocesses == 1 )
    {
        if(!tr1)
        {
            // zeromeout is only called here because other methods
            // will end up zeroing c out eventually.
            c->zeromeout();
            for(PS_SIT cc=0; cc<ncol; cc++)
                for(PS_SIT index=apc[cc]; index<(PS_SIT)apc[cc+1]; index++)
                    c->set(aic[index],
                           (double)((PS_FPT)c->get(aic[index]) + axc[index]*(PS_FPT)b->get(cc)));
        }
        else
        {
            for(PS_SIT cc=0; cc<ncol; cc++)
            {
                PS_FPT val=0;
                for(PS_SIT index=apc[cc]; index<(PS_SIT)apc[cc+1]; index++)
                    val += axc[index]*(PS_FPT)b->get(aic[index]);
                c->set(cc, (double)val);
            }
        }
    }
    /*
    // else if more than one thread
    else
    {
    pthread_t             threadshere[numprocesses];
    multmatrixvectorstruct mmvstructs[numprocesses];

    PS_SIT interval  = ncol/numprocesses;
    for(PS_SIT proc=0; proc<numprocesses; proc++)
    {
    mmvstructs[proc].tr = tr1;
    if( !tr1 )
    {
    MEMORY ps_malloc( &(mmvstructs[proc].subvec), (PS_SIT)(c->get_size()) );
    std::fill( mmvstructs[proc].subvec, mmvstructs[proc].subvec+c->get_size(), 0 );
    }
    mmvstructs[proc].a=this;
    mmvstructs[proc].b=b;
    mmvstructs[proc].c=c;
    mmvstructs[proc].lowercol= proc   *interval;
    mmvstructs[proc].uppercol=(proc+1)*interval;
    if(proc==numprocesses-1)
    mmvstructs[proc].uppercol = ncol;

    pthread_create(&threadshere[proc] ,&attr ,multmatrixvector, &mmvstructs[proc]);
    }
    for(PS_SIT proc=0; proc<numprocesses; proc++)
    pthread_join(threadshere[proc],NULL);

    // copy results from threads into vector c
    if( !tr1 )
    {
    PS_FPT sum;
    for( PS_SIT v=0; v<c->get_size(); ++v )
    {
    sum = 0;
    for( PS_SIT proc=0; proc<numprocesses; ++proc )
    sum += mmvstructs[proc].subvec[v];

    c->vec[v] = sum;
    }

    for(PS_SIT proc=0; proc<numprocesses; proc++)
    MEMORY ps_free( mmvstructs[proc].subvec );
    }
    }
    */

    // only multiply by matrix a's scalar because VECTOR::get() includes scalar
    c->set_scalar( scalar );
    c->set_status( pvs_cpu);
}

void pixsrc_matrix::plus(pixsrc_matrix *b, pixsrc_matrix *c, bool tr1, bool tr2, double scalar1_, double scalar2_, PS_SIT numprocesses )
{
    // add two matrix together

    if (tr1 || tr2 || scalar1_!=1.0)
        PRINTER printerror( data_->print2screenname,
                            "matrix addition: transpose not supported or\n"
                            "asked for unsupported scalar multiplier",
                            cdata_->print2screenmutex           );

    PS_FPT scalar1 = (PS_FPT)scalar1_;
    PS_FPT scalar2 = (PS_FPT)scalar2_;

    // dense stuff (for shapelets)
    if (is_dense)
    {
        // remove scalar from this matrix
        dissolve_scalar();
        // copy matrix b into matrix c
        std::copy (b->mat_dense, b->mat_dense+nnz, c->mat_dense);
        // dissolve scalar for matrix c
        std::transform (c->mat_dense, c->mat_dense+nnz, c->mat_dense,
                        std::bind1st(std::multiplies<PS_FPT>(),b->scalar*scalar2));
        // add this matrix to matrix c
        std::transform(c->mat_dense, c->mat_dense+nnz, mat_dense,
                       c->mat_dense, std::plus<double>());
        c->scalar = (PS_FPT)1;
        return;
    }

#ifdef __USE_PIXSRC_CUDA__

    if( are_we_using_cuda )
    {
        // let CUDA handle this
        CUDA plus( this, b, c, tr1, tr2, scalar1, scalar2 );
        return;
    }

#endif

    // update CPU
    this->update_cpu( pmf_csc );
    b   ->update_cpu( pmf_csc );
    c->init_cpu_vecs();
    c->scalar = (double)1.0;

    // no need yet for transpose stuff

    for(PS_SIT cc = 0; cc<ncol; cc++)
    {
        PS_SIT index1 = apc[cc];
        PS_SIT index2 = b->apc[cc];
        PS_SIT upper1 = apc[cc+1];
        PS_SIT upper2 = b->apc[cc+1];
        bool done1 = false;
        bool done2 = false;

        for(;;)
        {
            if(index1==upper1 && index2==upper2)
            {
                done1 = done2 = true;
                break;
            }
            else if(index1==upper1)
            {
                done1 = true;
                break;
            }
            else if(index2==upper2)
            {
                done2 = true;
                break;
            }

            if(aic[index1] == b->aic[index2])
            {
                c->set(aic[index1],cc,
                       (double)(scalar1*scalar*axc[index1]+scalar2*b->scalar*b->axc[index2]));
                index1++;
                index2++;
            }
            else if(aic[index1] < b->aic[index2])
            {
                while(aic[index1] < b->aic[index2])
                {
                    if(index1==upper1)
                        break;
                    c->set(aic[index1],cc,(double)(scalar1*scalar*axc[index1]));
                    index1++;
                }
            }
            else
            {
                while(b->aic[index2] < aic[index1])
                {
                    if(index2==upper2)
                        break;
                    c->set(b->aic[index2],cc,(double)(scalar2*b->scalar*b->axc[index2]));
                    index2++;
                }
            }
        }

        if(!done1)
            for(; index1<upper1; index1++)
                c->set(aic[index1],cc,(double)(scalar1*scalar*axc[index1]));
        else if(!done2)
            for(; index2<upper2; index2++)
                c->set(b->aic[index2],cc,(double)(scalar2*b->scalar*b->axc[index2]));
    }
}

void pixsrc_matrix::createnumsymb()
{
    // create factorizations of the matrix

    // dense stuff (for shapelets)
    if (is_dense)
    {
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: 1",
                            cdata_->print2screenmutex);
        return;
    }

    // lock mutex
    pthread_mutex_lock(&numsymbmutex);

    // if there's nothing to do, return
    if( ludone )
    {
        pthread_mutex_unlock(&numsymbmutex);
        return;
    }

    // update CPU (UMFPACK needs column compression)
    this->update_cpu( pmf_csc );

#ifndef DOUBLE_PRECISION
    PRINTER printerror( data_->print2screenname,
                        "double precision required for UMFPACK: analysis",
                        cdata_->print2screenmutex                );

#else

    // call UMFPACK
#ifndef INTEGER_TYPE_I
    umfpack_dl_symbolic( nrow, ncol, apc, aic, axc,
                         &symbolic, (PS_FPT*)NULL, (PS_FPT*)NULL );
    umfpack_dl_numeric ( apc,  aic,  axc,
                         symbolic, &numeric,
                         (PS_FPT*)NULL, (PS_FPT*)NULL );
#else
    umfpack_di_symbolic( nrow, ncol, apc, aic, axc,
                         &symbolic, (PS_FPT*)NULL, (PS_FPT*)NULL );
    umfpack_di_numeric ( apc,  aic,  axc,
                         symbolic, &numeric,
                         (PS_FPT*)NULL, (PS_FPT*)NULL );
#endif

#endif

    // update and unlock mutex
    ludone = 1;
    pthread_mutex_unlock( &numsymbmutex );
}

double pixsrc_matrix::logdet()
{
    // get logarithm of determinant

    // dense stuff (for shapelets)
    if (is_dense)
    {
        dense_matrix_stuff (0,0,0,0,0,1);
        return det_dense;
    }

/*
  if( matrix_type == 11 )
  {
  return nrow * log( fabs( scalar ) );
  }


  else */if( matrix_type == 10 )
    {
        // compute logarithm of diagonal matrix

        // update cpu
        this->update_cpu( pmf_csc );

        PS_FPT val = (PS_FPT)0;
        for( PS_SIT r=0; r<nrow; ++r )
            val += log( fabs( tkv[r] ) );
        return (double)( nrow * log( fabs( scalar ) ) + val );
    }

    else
    {
        // factorize matrix
        if(!ludone)
            createnumsymb();

        // let UMFPACK get determinant
        PS_FPT m[1];
        PS_FPT e[1];
        e[0] = (PS_FPT)(-100);

#ifndef DOUBLE_PRECISION
        PRINTER printerror( data_->print2screenname,
                            "double precision required for UMFPACK: determinant",
                            cdata_->print2screenmutex                );

#else

        // call UMFPACK, which returns mantissa and exponent
#ifndef INTEGER_TYPE_I
        umfpack_dl_get_determinant (m, e, numeric, NULL);
#else
        umfpack_di_get_determinant (m, e, numeric, NULL);
#endif

#endif

        // compute natural logarithm
        return (double)( nrow*log(fabs(scalar)) + log(fabs(m[0])) +
                         (PS_FPT)2.302585092994046*e[0]             );
    }
}

void pixsrc_matrix::inplace_inverse ()
{
    if (!is_dense)
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: sparse in place inverse attempted",
                            cdata_->print2screenmutex);

    // get inverse of matrix in place
    int *ipiv;
    PS_SIT info, lwork=-1, tnrow=(PS_SIT)nrow, tncol=(PS_SIT)ncol;
    double wkopt, *work;
    MEMORY ps_malloc (&ipiv, (PS_SIT)nrow+1);
    EXTERNAL ps_dgetrf_ (&tnrow, &tncol, mat_dense, &tnrow, ipiv, &info);
    EXTERNAL ps_dgetri_ (&tnrow, mat_dense, &tnrow, ipiv, &wkopt, &lwork, &info);
    lwork = (PS_SIT)wkopt;
    MEMORY ps_malloc (&work, (PS_SIT)lwork);
    EXTERNAL ps_dgetri_ (&tnrow, mat_dense, &tnrow, ipiv, work, &lwork, &info);
    MEMORY ps_free (work);
    MEMORY ps_free (ipiv);
}

void pixsrc_matrix::inplace_lu_fac ()
{
    if (!is_dense)
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: sparse in place LU factrization attempted",
                            cdata_->print2screenmutex);

    this->dense_matrix_stuff (0, 0, 0, 0, 0, -1);
}

void pixsrc_matrix::inplace_leq (MATRIX *b)
{
    if (!is_dense)
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: sparse in place leq attempted",
                            cdata_->print2screenmutex);

    this->dense_matrix_stuff (b, 0, 0, 0, 0, 6);
}

void pixsrc_matrix::dense_matrix_stuff (MATRIX *b1, MATRIX *c1, pixsrc_vector *x, pixsrc_vector *b, double *ans, PS_SIT option)
{
    // dense matrix stuff (for shapelets)
    // option determiners which operation will be performed.

    if (!is_dense)
        return;

    // update cpu and remove scalar
    update_cpu (pmf_dense);
    dissolve_scalar();
    PS_SIT tnrow = (PS_SIT)nrow, tncol = (PS_SIT)ncol;

    // do in-place LU factorization
    if (!got_lu_fac && -1==option)
    {
        pthread_mutex_lock (&dense_lu_mutex);
        if (!got_lu_fac)
        {
            PS_SIT info;
            if (!ipiv_dense)
                MEMORY ps_malloc (&ipiv_dense, (PS_SIT)nrow+1);
            if (lu_fac_dense)
                MEMORY ps_free (lu_fac_dense);
            lu_fac_dense = mat_dense;
            mat_dense = NULL;
            EXTERNAL ps_dgetrf_ (&tnrow, &tncol, lu_fac_dense, &tnrow, ipiv_dense, &info);
            got_lu_fac = 1;
        }
        pthread_mutex_unlock (&dense_lu_mutex);
    }

    // do LU factorization
    if (!got_lu_fac && option)
    {
        pthread_mutex_lock (&dense_lu_mutex);
        if (!got_lu_fac)
        {
            PS_SIT info;
            if (!ipiv_dense)
                MEMORY ps_malloc (&ipiv_dense, (PS_SIT)nrow+1);
            if (!lu_fac_dense)
                MEMORY ps_malloc (&lu_fac_dense, (PS_SIT)nnz);
            std::copy (mat_dense, mat_dense+nnz, lu_fac_dense);
            EXTERNAL ps_dgetrf_ (&tnrow, &tncol, lu_fac_dense, &tnrow, ipiv_dense, &info);
            got_lu_fac = 1;
        }
        pthread_mutex_unlock (&dense_lu_mutex);
    }

    // get determinant
    if (!got_det && option==1)
    {
        det_dense = (PS_FPT)0;
        PS_FPT* ptr = lu_fac_dense;
        for (PS_SIT s=0; s<nrow; ++s)
        {
            det_dense += log (fabs (*ptr));
            ptr += nrow+1;
        }
        got_det = 1;
    }

    // get inverse
    else if (!got_inv && option==2)
    {
        pthread_mutex_lock (&dense_inv_mutex);
        if (!got_inv)
        {
            if (!inv_dense)
                MEMORY ps_malloc (&inv_dense, (PS_SIT)nnz);
            std::copy (lu_fac_dense, lu_fac_dense+nnz, inv_dense);
            PS_SIT info;
            PS_SIT lwork=-1;     // "workspace"-related variable
            double wkopt;     // workspace neeed
            double *work=0;   // actual workspace vector
            EXTERNAL ps_dgetri_ (&tnrow, inv_dense, &tnrow, ipiv_dense, &wkopt, &lwork, &info);
            // allocate workspace
            lwork = (PS_SIT)wkopt;
            MEMORY ps_malloc (&work, (PS_SIT)lwork);
            // perform matrix inversion
            EXTERNAL ps_dgetri_ (&tnrow, inv_dense, &tnrow, ipiv_dense, work, &lwork, &info);
            MEMORY ps_free (work);
            got_inv = 1;
        }
        pthread_mutex_unlock (&dense_inv_mutex);
    }

    // solve linear equation
    else if (option==3)
    {
        char tr = 'N';
        PS_SIT rhs = 1, info;
        b->update_cpu();
        x->update_cpu();
        b->dissolve_scalar();
        std::copy (b->vec, b->vec+b->size, x->vec);
        EXTERNAL ps_dgetrs_ (&tr, &tnrow, &rhs, lu_fac_dense,
                             &tnrow, ipiv_dense, x->vec,
                             &tnrow, &info);
        x->scalar = (PS_FPT)1;
    }

    // solve system of linear equations and get trace
    else if (option==4)
    {
        // solve equations
        char tr = 'N';
        PS_SIT rhs = b1->ncol, info;
        b1->update_cpu (pmf_dense);
        double *product;
        MEMORY ps_malloc (&product, (PS_SIT)(b1->nnz));
        std::copy (b1->mat_dense, b1->mat_dense+b1->nnz, product);
        EXTERNAL ps_dgetrs_ (&tr, &tnrow, &rhs, lu_fac_dense,
                             &tnrow, ipiv_dense, product,
                             &tnrow, &info);
        // get trace
        *ans = 0;
        double *ptr = product;
        for (PS_SIT s=0; s<b1->ncol; ++s)
        {
            *ans += *ptr;
            ptr += b1->nrow+1;
        }
        *ans *= b1->scalar;
        MEMORY ps_free (product);
    }

    // solve system of linear equations
    else if (option==5)
    {
        // solve equations
        char tr = 'N';
        PS_SIT rhs = (PS_SIT)b1->ncol, info;
        b1->update_cpu (pmf_dense);
        std::copy (b1->mat_dense, b1->mat_dense+b1->nnz, c1->mat_dense);
        EXTERNAL ps_dgetrs_ (&tr, &tnrow, &rhs, lu_fac_dense,
                             &tnrow, ipiv_dense, c1->mat_dense,
                             &tnrow, &info);
        c1->mult (b1->scalar);
    }

    // solve system of linear equations
    else if (option==6)
    {
        // solve equations
        char tr = 'N';
        PS_SIT rhs = (PS_SIT)b1->ncol, info;
        b1->update_cpu (pmf_dense);
        EXTERNAL ps_dgetrs_ (&tr, &tnrow, &rhs, lu_fac_dense,
                             &tnrow, ipiv_dense, b1->mat_dense,
                             &tnrow, &info);
    }
}

void pixsrc_matrix::linequationsolve(pixsrc_vector *x, pixsrc_vector *b, bool zeroit )
{
    // solve linear matrix equation

    // dense matrix stuff (for shapelets)
    if (is_dense)
    {
        dense_matrix_stuff (0,0,x,b,0,3);
        return;
    }

    zeroit = 1;

    ///////////////////////////////////////
    // DISABLING CUDA FOR THIS OPERATION //
    ///////////////////////////////////////
    //
    // I think there were some performance issues with standard CUDA library.
    // cusp library may not have been fully working.
    // Cholesky decomposition for manually solving equation fails with CUDA for large matrices,
    // which may possibly be a shortage of GPU memory issue.

#ifdef __USE_PIXSRC_CUDA__

    if( 0 )
        if( are_we_using_cuda )
        {
            CUDA tri_eq_solve( this, x, b, zeroit );

            return;
        }

#endif


    // factorize
    if( !ludone )
        createnumsymb();

    // prepare vectors
    double *tempx;
    MEMORY ps_malloc( &tempx, (PS_SIT)ncol );

    b->dissolve_scalar();

    const double *tempb = b->get_vec_ptr();

    /*
      void *clock;
      PS_FPT time;
      CUDA ps_start_clock( &clock );
    */

#ifndef DOUBLE_PRECISION
    PRINTER printerror( data_->print2screenname,
                        "double precision required for UMFPACK: solve",
                        cdata_->print2screenmutex                );

#else

    // call UMFPACK
#ifndef INTEGER_TYPE_I
    umfpack_dl_solve ( UMFPACK_A, apc, aic, axc,
                       tempx, tempb, numeric, (PS_FPT*)NULL, (PS_FPT*)NULL );
#else
    umfpack_di_solve ( UMFPACK_A, apc, aic, axc,
                       tempx, tempb, numeric, (PS_FPT*)NULL, (PS_FPT*)NULL );
#endif

#endif

    /*
      CUDA ps_stop_clock( clock, &time );
      std::cout << "umfpack took " << time << " ms" << std::endl;
    */

    // copy solution into x
    for(PS_SIT g=0; g<ncol; g++)
        x->set( g, tempx[g] );
    x->mult( (double)( (PS_FPT)1 / scalar ) );

    MEMORY ps_free( tempx );

}

// structure for inverting matrix
struct noise_struct
{
    MATRIX *a;
    PS_SIT lower, upper;
    VECTOR *noise;
};

void* pixsrc_matrix::noise_invA_body( void *args )
{
    // code passed to pthreads for multithreading, for inverting matrix

    noise_struct *ts = (noise_struct*)args;

    VECTOR *dummy1, *dummy2, *res, *col;
    MEMORY ps_malloc( &dummy1, 1 );
    MEMORY ps_malloc( &dummy2, 1 );

    res = new (dummy1) VECTOR( ts->a->cdata_, ts->a->data_, ts->a->ncol );
    col = new (dummy2) VECTOR( ts->a->cdata_, ts->a->data_, ts->a->nrow );

    if( res->are_we_using_cuda )
    {
        // have to initilalize vecs here because ideally, this would
        // be done on the gpu and not on the cpu
        res->init_cpu_vec();
        col->init_cpu_vec();
    }

    for(PS_SIT g=ts->lower; g<ts->upper; g++)
    {
        std::fill( col->vec, col->vec + col->get_size(), (PS_FPT)0 );
        col->vec[g] = (PS_FPT)1;;

        ts->a->linequationsolve(res,col,1);
        ts->noise->set(g,std::sqrt(res->get(g)));
    }

    res->~VECTOR();
    col->~VECTOR();

    MEMORY ps_free( dummy1 );
    MEMORY ps_free( dummy2 );

    return NULL;
}

void pixsrc_matrix::noise_invA( VECTOR *noise, PS_SIT numthreads, PS_SIT *err  )
{
    // error check
    if( err && *err )
    {
        return;
    }

    // update CPU
    noise->update_cpu();

    // dense stuff (for shapelets)
    if (is_dense)
    {
        dense_matrix_stuff (0,0,0,0,0,2);
        PS_FPT *ptr = noise->vec;
        for (PS_SIT s=0; s<noise->size; ++s)
        {
            *ptr = (PS_FPT)inv_dense[s*nrow+s];
            ++ptr;
        }
        return;
    }

    // set up and launch multiple threads
    PS_SIT interval = ncol/numthreads;
    noise_struct noise_structs[numthreads];
    pthread_t threadshere[numthreads];

    for(PS_SIT proc=0; proc<numthreads; proc++)
    {
        noise_structs[proc].a  = this;
        noise_structs[proc].lower= proc   *interval;
        noise_structs[proc].upper=(proc+1)*interval;
        noise_structs[proc].noise=noise;
        if(proc==numthreads-1)
            noise_structs[proc].upper = ncol;

        pthread_create(&threadshere[proc], &attr, noise_invA_body, &noise_structs[proc]);
    }

    for(PS_SIT proc=0; proc<numthreads; proc++)
        pthread_join(threadshere[proc],NULL);
}

// structure for calculating trace of inverse(A)*B
struct tracestruct
{
    MATRIX *a, *b;
    PS_SIT lower, upper;
    PS_FPT tr;
};

void* pixsrc_matrix::trace_invA_B_body( void *args )
{
    // code passed to pthreads for multithreading, for calculating trace of inverse(A)*B

    tracestruct *ts = (tracestruct*)args;

    VECTOR *dummy1, *dummy2, *res, *col;
    MEMORY ps_malloc( &dummy1, 1 );
    MEMORY ps_malloc( &dummy2, 1 );

    res = new (dummy1) VECTOR( ts->a->cdata_, ts->a->data_, ts->a->ncol );
    col = new (dummy2) VECTOR( ts->a->cdata_, ts->a->data_, ts->b->nrow );

    if( res->are_we_using_cuda )
    {
        // have to initilalize vecs here because ideally, this would
        // be done on the gpu and not on the cpu
        res->init_cpu_vec();
        col->init_cpu_vec();
    }

    col->set_scalar( ts->b->scalar );

    for(PS_SIT g=ts->lower; g<ts->upper; g++)
    {
        /*
          void *clock_ptr;
          PS_FPT elapsed_time;
          CUDA ps_start_clock( &clock_ptr );
        */

        // copy column of MATRIX *b into VECTOR *col
        std::fill( col->vec, col->vec + col->get_size(), (PS_FPT)0 );
        for(PS_SIT g2=ts->b->apc[g]; g2<ts->b->apc[g+1]; g2++)
            col->vec[ts->b->aic[g2]] = ts->b->axc[g2];

        ts->a->linequationsolve(res,col,1);
        ts->tr += (PS_FPT)res->get(g);

        /*
          CUDA ps_stop_clock( clock_ptr, &elapsed_time );
          string et = "chicca took " + OPERA tostring(elapsed_time) + " ms to execute\n";
          PRINTER print2screen("", et, ts->b->cdata_->print2screenmutex);
        */
    }

    res->~VECTOR();
    col->~VECTOR();

    MEMORY ps_free( dummy1 );
    MEMORY ps_free( dummy2 );

    return NULL;
}

double pixsrc_matrix::trace_invA_B( MATRIX *b, PS_SIT numthreads, PS_SIT *err  )
{
    // calculate trace of inverse(A)*B

    // error check
    if( err && *err )
    {
        return (double)0.0;
    }

    /*
      void *clock_ptr;
      PS_FPT elapsed_time;
      CUDA ps_start_clock( &clock_ptr );
    */

    // dense stuff (for shapelets)
    if (is_dense)
    {
        double tr;
        dense_matrix_stuff (b,0,0,0,&tr,4);
        return tr;
    }


    PS_FPT tr = (PS_FPT)0;

    ///////////////////////////////////////
    // DISABLING CUDA FOR THIS OPERATION //
    ///////////////////////////////////////

#ifdef __USE_PIXSRC_CUDA__

    if( 0 )
        if( are_we_using_cuda )
        {
            CUDA tr_inv_a_b( this, b, &tr );

            //CUDA tr_invA_B_setup( this, b      );
            //CUDA tr_invA_B_cycle( this, b, &tr );

            /*
              CUDA ps_stop_clock( clock_ptr, &elapsed_time );
              string et = "total time " + OPERA tostring(elapsed_time) + " ms to execute\n";
              PRINTER print2screen("", et, cdata_->print2screenmutex);
            */

            return (double)tr;
        }

#endif


    // set up and launch multiple threads
    b->update_cpu( pmf_csc );

    PS_SIT interval = ncol/numthreads;
    tracestruct tracestructs[numthreads];
    pthread_t threadshere[numthreads];

    for(PS_SIT proc=0; proc<numthreads; proc++)
    {
        tracestructs[proc].tr = 0;
        tracestructs[proc].a  = this;
        tracestructs[proc].b  = b;
        tracestructs[proc].lower= proc   *interval;
        tracestructs[proc].upper=(proc+1)*interval;
        if(proc==numthreads-1)
            tracestructs[proc].upper = ncol;

        pthread_create(&threadshere[proc], &attr, trace_invA_B_body, &tracestructs[proc]);
    }

    for(PS_SIT proc=0; proc<numthreads; proc++)
        pthread_join(threadshere[proc],NULL);
    for(PS_SIT proc=0; proc<numthreads; proc++)
        tr += tracestructs[proc].tr;

    /*
      CUDA ps_stop_clock( clock_ptr, &elapsed_time );
      string et = "total time " + OPERA tostring(elapsed_time) + " ms to execute\n";
      PRINTER print2screen("", et, cdata_->print2screenmutex);
    */

    return (double)tr;
}

void pixsrc_matrix::atasbc( double *adj, PS_SIT num )
{
    // this function adds matrix entries to a matrix thate already contain
    // entries for that matrix element.
    // This is only used in pixsrc_reg_greengauss1

    if (is_dense)
    {
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: 2",
                            cdata_->print2screenmutex);
        return;
    }

    // atasbc = add_to_already_set_before_compiling

    // this function is for matrices that have not yet had compile()
    // or creatercform() called on them.
    // it adjusts existing matrix entries by adding user-given values

    PS_SIT curr_row_start = 0;
    PS_SIT curr_row_val;
    PS_SIT curr_row_end = -1;

    while( curr_row_end != nnz )
    {
        curr_row_val = tiv[curr_row_start];

        for (PS_SIT n=curr_row_start; n<nnz; ++n )
        {
            if( (PS_SIT)tiv[n] != curr_row_val )
            {
                curr_row_end = n;
                break;
            }
            if( n == nnz-1 )
                curr_row_end = nnz;
        }

        for( PS_SIT g=0; g<num; ++g )
        {
            if( (PS_SIT)adj[g*3] != curr_row_val )
                continue;

            for( PS_SIT n=curr_row_start; n<curr_row_end; ++n )
            {
                if( (PS_SIT)adj[g*3+1] == tjv[n] )
                {
                    tkv[n] += (PS_FPT)adj[g*3+2];
                    adj[g*3] = -1;
                    break;
                }
            }
        }

        curr_row_start = curr_row_end;
    }

    // checking if any entries in adj didnt get used
    for( PS_SIT g=0; g<num; ++g )
    {
        if( (PS_SIT)adj[g*3] != -1 )
        {
            PRINTER printerror( data_->print2screenname,
                                "unused entries in atasbc: "+OPERA tostring(adj[g*3]),
                                cdata_->print2screenmutex   );
        }
    }
}

void pixsrc_matrix::print_me( PS_SIT tracker, string basename, string imagename, bool num,
                              string filename, bool nonstandardname )
{
    // print contents of matrix -- mainly for debugging.

    if (is_dense)
    {
        PRINTER printerror (data_->print2screenname,
                            "dense matrix error: 3",
                            cdata_->print2screenmutex);
        return;
    }
    update_cpu( pmf_csc );

    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    std::ofstream file;
    file.open(filename.c_str());
    for(PS_SIT j=0; j<ncol; j++)
    {
        for(PS_SIT k=apc[j]; k<apc[j+1]; k++)
        {
            file << std::setprecision(15) << aic[k] << "\t" << j << "\t" <<
                axc[k]*scalar << "\n";
        }
    }
}

/////////////////////////////////////
// memory cleanup code starts here //
/////////////////////////////////////

void pixsrc_matrix::clearcols()
{
    MEMORY ps_free( apc );
    MEMORY ps_free( aic );
    MEMORY ps_free( axc );

    apc = (PS_SIT*   )NULL;
    aic = (PS_SIT*   )NULL;
    axc = (PS_FPT*)NULL;

    hasbeencompiled = 0;
}

void pixsrc_matrix::clearrows()
{
    MEMORY ps_free( apr );
    MEMORY ps_free( air );
    MEMORY ps_free( axr );

    apr = (PS_SIT*)NULL;
    air = (PS_SIT*)NULL;
    axr = (PS_FPT*)NULL;

    rcform = 0;

    MEMORY ps_free( apr_sym );
    MEMORY ps_free( air_sym );
    MEMORY ps_free( axr_sym );

    apr_sym = (PS_SIT*)NULL;
    air_sym = (PS_SIT*)NULL;
    axr_sym = (PS_FPT*)NULL;

}

void pixsrc_matrix::clearvoids()
{
    if(ludone)
    {

#ifndef DOUBLE_PRECISION
        PRINTER printerror( data_->print2screenname,
                            "double precision required for UMFPACK: free analysis",
                            cdata_->print2screenmutex                );
#endif

#ifndef INTEGER_TYPE_I
        umfpack_dl_free_symbolic(&symbolic);
        umfpack_dl_free_numeric (&numeric );
#else
        umfpack_di_free_symbolic(&symbolic);
        umfpack_di_free_numeric (&numeric );
#endif
    }

    symbolic = (void*)NULL;
    numeric =  (void*)NULL;
}

void pixsrc_matrix::resetdense(PS_SIT nrow_, PS_SIT ncol_)
{
    // This function removes anything the CPU has done, so that
    // I don't have to re-create a matrix do do multiple operations on it.

    if (!is_dense)
        return;

    if (nrow_!=nrow || ncol_!=ncol)
    {
        nrow = nrow_;
        ncol = ncol_;
        nnz = nrow*ncol;
        MEMORY ps_free (mat_dense);
        MEMORY ps_free (inv_dense);
        MEMORY ps_free (lu_fac_dense);
        mat_dense = (PS_FPT*)0;
        inv_dense = (PS_FPT*)0;
        lu_fac_dense = (PS_FPT*)0;
        MEMORY ps_malloc (&mat_dense, (PS_SIT)nnz);
    }

    got_det = got_lu_fac = got_inv = 0;
    std::fill (mat_dense, mat_dense+nnz, (PS_FPT)0);
}

void pixsrc_matrix::cleardense()
{
    MEMORY ps_free (mat_dense);
    MEMORY ps_free (inv_dense);
    MEMORY ps_free (ipiv_dense);
    MEMORY ps_free (lu_fac_dense);

    mat_dense    = (PS_FPT*)0;
    inv_dense    = (PS_FPT*)0;
    lu_fac_dense = (PS_FPT*)0;
    ipiv_dense   = (int*/*PS_SIT**/)0;
}

void pixsrc_matrix::clearvecs()
{
    //if(hasbeencompiled)
    //    return;

    MEMORY ps_free( tiv );
    MEMORY ps_free( tjv );
    MEMORY ps_free( tkv );

    tiv = 0;
    tjv = 0;
    tkv = (PS_FPT*)0;

    vecsize = 0;
}

void pixsrc_matrix::clearall()
{
    // don't want to double-free so eliminating
    // multiple pointers to same memory

    PS_SIT i_num = 10;
    PS_SIT d_num = 5;

    PS_SIT **i_ptr[i_num];
    PS_FPT **d_ptr[d_num];

    i_ptr[0] = &apc;
    i_ptr[1] = &apr;
    i_ptr[2] = &aic;
    i_ptr[3] = &air;
    i_ptr[4] = &tiv;
    i_ptr[5] = &tjv;
    i_ptr[6] = &apr_sym;
    i_ptr[7] = &air_sym;
    i_ptr[8] = &icf_ptr_cpu;
    i_ptr[9] = &icf_ind_cpu;

    d_ptr[0] = &axc;
    d_ptr[1] = &axr;
    d_ptr[2] = &tkv;
    d_ptr[3] = &axr_sym;
    d_ptr[4] = &icf_val_cpu;

    // remove degenerate int pointers
    for( PS_SIT i1=0; i1<i_num; ++i1 )
    {
        if( !*i_ptr[i1] )
            continue;

        for( PS_SIT i2=i1+1; i2<i_num; ++i2 )
            if( *i_ptr[i2] && *i_ptr[i1] == *i_ptr[i2] )
                *i_ptr[i2] = 0;
    }

    // remove degenerate PS_FPT pointers
    for( PS_SIT d1=0; d1<d_num; ++d1 )
    {
        if( !*d_ptr[d1] )
            continue;

        for( PS_SIT d2=d1+1; d2<d_num; ++d2 )
            if( *d_ptr[d2] && *d_ptr[d1] == *d_ptr[d2] )
                *d_ptr[d2] = 0;
    }

    MEMORY ps_free( icf_ptr_cpu );
    MEMORY ps_free( icf_ind_cpu );
    MEMORY ps_free( icf_val_cpu );
    icf_ptr_cpu = 0;
    icf_ind_cpu = 0;
    icf_val_cpu = (PS_FPT*)0;

    // clear em all!!
    clearvecs();
    clearcols();
    clearrows();
    clearvoids();
    cleardense();


#ifdef __USE_PIXSRC_CUDA__

    clearcudavec();

    /*
      CUDA ps_dest_csrsvi( this );
      CUDA ps_dest_descr ( this );
    */

    if( are_we_using_cuda )
    {
        CUDA destroy_stream( this );
        CUDA destroy_cusp_prec( this );
    }

#endif

}

#ifdef __USE_PIXSRC_CUDA__

void pixsrc_matrix::clearcudavec()
{
    // don't want to PS_FPT-free so eliminating
    // multiple pointers to same memory

    PS_SIT i_num = 6;
    PS_SIT d_num = 3 + num_ctv;

    PS_SIT    **i_ptr[i_num];
    PS_FPT **d_ptr[d_num];

    i_ptr[0] = &csx;
    i_ptr[1] = &csx_tr;
    i_ptr[2] = &coo_col;
    i_ptr[3] = &coo_row_tr;
    i_ptr[4] = &icf_ptr;
    i_ptr[5] = &icf_ind;
    //i_ptr[4] = &csx_sym;
    //i_ptr[5] = &coo_col_sym;

    d_ptr[0] = &coo_val;
    d_ptr[1] = &coo_val_tr;
    d_ptr[2] = &icf_val;
    //d_ptr[3] = &coo_val_sym;
    for( PS_SIT j=3; j<3+num_ctv; ++j)
        d_ptr[j] = &ctv[j-3];

    // remove degenerate int pointers
    for( PS_SIT i1=0; i1<i_num; ++i1 )
    {
        if( !*i_ptr[i1] )
            continue;

        for( PS_SIT i2=i1+1; i2<i_num; ++i2 )
            if( *i_ptr[i2] && *i_ptr[i1] == *i_ptr[i2] )
                *i_ptr[i2] = 0;
    }

    // remove degenerate PS_FPT pointers
    for( PS_SIT d1=0; d1<d_num; ++d1 )
    {
        if( !*d_ptr[d1] )
            continue;

        for( PS_SIT d2=d1+1; d2<d_num; ++d2 )
            if( *d_ptr[d2] && *d_ptr[d1] == *d_ptr[d2] )
                *d_ptr[d2] = 0;
    }

    CUDA ps_cuda_free( csx         );
    CUDA ps_cuda_free( coo_col     );
    CUDA ps_cuda_free( coo_val     );
    CUDA ps_cuda_free( csx_tr      );
    CUDA ps_cuda_free( coo_row_tr  );
    CUDA ps_cuda_free( coo_val_tr  );
    CUDA ps_cuda_free( icf_ind     );
    CUDA ps_cuda_free( icf_ptr     );
    CUDA ps_cuda_free( icf_val     );
    /*
      CUDA ps_cuda_free( coo_col_sym );
      CUDA ps_cuda_free( coo_val_sym );
      CUDA ps_cuda_free( csx_sym     );
    */
    for( PS_SIT j=0; j<num_ctv; ++j )
        CUDA ps_cuda_free( ctv[j] );
    MEMORY ps_free( ctv );


    csx         = 0;
    coo_col     = 0;
    coo_val     = (PS_FPT*)0;
    csx_tr      = 0;
    coo_row_tr  = 0;
    coo_val_tr  = (PS_FPT*)0;
    icf_ptr     = 0;
    icf_ind     = 0;
    icf_val     = (PS_FPT*)0;
    /*
      coo_col_sym = 0;
      coo_val_sym = 0;
      csx_sym     = 0;
    */

}

void pixsrc_matrix::cuda_reset()
{
    // This function removes anythis CUDA has done, so that
    // I don't have to re-create a matrix do do multiple operations on it.

    // it is assumed that csx and csx_tr variables
    // will remain unchanged
    // it is also assumed that the NNZ is    /* will eventually be */
    // the same as it is now, so that re-allocation
    // of all vectors isn't necessary

    // I don't know the structure of the following,
    // so I destroy them
    clearvoids();
    symbolic = (void*)NULL;
    numeric =  (void*)NULL;

    //nnz             = 0;
    scalar          = (PS_FPT)1;
    rcform          = 0;
    ludone          = 0;
    hasbeencompiled = 0;
    icf_done = 0;


    /*
      clearall( this );

      scalar = 1;
      nnz = 0;
      rcform = false;
      ludone = false;
      hasbeencompiled = false;

      apc = (PS_SIT*)NULL;
      apr = (PS_SIT*)NULL;
      aic = (PS_SIT*)NULL;
      air = (PS_SIT*)NULL;
      axc = (PS_FPT*)NULL;
      axr = (PS_FPT*)NULL;
      symbolic = (void*)NULL;
      numeric = (void*)NULL;

      csx        = 0;
      coo_col    = 0;
      coo_val    = 0;
      csx_tr     = 0;
      coo_row_tr = 0;
      coo_val_tr = 0;

      csrsvi = 0;

      tiv = 0;
      tjv = 0;
      tkv = 0;
      vecsize = 0;
    */
}

#endif
/*
  void pixsrc_matrix::reset()
  {

  #ifdef __USE_PIXSRC_CUDA__

  cuda_reset();
  return;

  #endif

  clearvoids();
  symbolic = (void*)NULL;
  numeric =  (void*)NULL;

  }
*/

/*
  void pixsrc_matrix::selective_copy( PS_SIT *d_i, PS_SIT *d_p, PS_FPT *d_x,
  PS_SIT *s_i, PS_SIT *s_p, PS_FPT *s_x, PS_SIT size )
  {
  // This function is really used anymore.

  if (is_dense)
  {
  PRINTER printerror (data_->print2screenname,
  "dense matrix error: 4",
  cdata_->print2screenmutex);
  return;
  }

  PS_SIT s_ind;

  // loop through dest rows
  for( PS_SIT r=0; r<size; ++r )
  {
  // set source column index
  s_ind = s_p[r];

  // loop through dest columns
  for( PS_SIT d_ind=d_p[r]; d_ind<d_p[r+1]; ++d_ind )
  {
  // loop through source columns
  for( ; ; ++s_ind )
  {
  // if source and dest column index match
  if( d_i[d_ind] == s_i[s_ind] )
  {
  // copy this element
  d_x[d_ind] = s_x[s_ind];
  ++s_ind;
  break;
  }
  // if it can't possibly find it
  else if( s_i[s_ind] > d_i[d_ind] )
  {
  d_x[d_ind] = (PS_FPT)0;
  break;
  }
  }
  }
  }
  //return;
  PS_SIT row=0;
  for( PS_SIT j=0; j<d_p[size]; ++j )
  {
  //
  //  for( PS_SIT g=1; g<size+1; ++g )
  //  if( d_p[g] > j )
  //  {
  //  row = g-1;
  //  break;
  //  }
  //

  //
  //if( j == d_p[row+1] )
  //    ++row;
  //

  if( d_x[j] == (PS_FPT)0 )
  {
  std::transform( d_p + row+1, d_p + size+1, d_p + row+1, op_decrease );
  std::memmove( d_i + j, d_i + j+1, (d_p[size]-j)*sizeof( PS_SIT  ));
  std::memmove( d_x + j, d_x + j+1, (d_p[size]-j)*sizeof(PS_FPT));
  --j;
  }
  }
  }
*/

/*
  void pixsrc_matrix::remove_duplicates()
  {
  // This function has been replaced by atasbc

  if (is_dense)
  return;

  for( PS_SIT j=0; j<nnz-1; ++j )
  {
  if( tiv[j] == tiv[j+1] && tjv[j] == tjv[j+1] )
  {
  tkv[j] += tkv[j+1];
  std::memmove( tiv + j+1, tiv + j+2, (nnz-j-2)*sizeof(  PS_SIT ));
  std::memmove( tjv + j+1, tjv + j+2, (nnz-j-2)*sizeof(  PS_SIT ));
  std::memmove( tkv + j+1, tkv + j+2, (nnz-j-2)*sizeof(PS_FPT));
  --nnz;
  --j;
  }
  }
  }
*/

/*
  void pixsrc_matrix::copycolumn(PS_SIT col, pixsrc_vector *v)
  {
  // this function assumes the matrix it is called on has had
  // update_cpu( 1 ) called on it.
  // it also assumes the scalar for v has been equated to this scalar.

  std::fill( v->vec, v->vec + v->size, 0 );

  for(PS_SIT g=apc[col]; g<apc[col+1]; g++)
  v->vec[aic[g]] = axc[g];
  }
*/

/*
  void pixsrc_matrix::copymatrix( MATRIX *dest )
  {
  if(!hasbeencompiled)
  compile();

  for(PS_SIT cc=0; cc<ncol; cc++)
  for(PS_SIT index=apc[cc]; index<apc[cc+1]; index++)
  dest->set( aic[index], cc, axc[index] );
  dest->scalar = scalar;
  }
*/

/*
  PS_FPT pixsrc_matrix::get( PS_SIT r )
  {
  // only diagonal matrices should use this function
  // matrices proportional to identity should not

  if( !hasbeencompiled )
  compile();

  return scalar*axc[r];
  }
*/
