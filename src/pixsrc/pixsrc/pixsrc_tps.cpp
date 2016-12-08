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


// This class creates a thin plate spline

#include "pixsrc_tps.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer.hpp"
#include "pixsrc_wcs.hpp"
#include <gsl/gsl_sf_erf.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_bessel.h>


pixsrc_tps::pixsrc_tps(commoninputdata *cdata__, inputdata *data__, const double *x, const double *y, const double *z, PS_SIT xstride, PS_SIT ystride, PS_SIT zstride, PS_SIT numpts_, ps_rbf_type type, double rbf_parm, PS_SIT integrateme, PS_SIT addlinearterms, PS_SIT transmode)
{
    is_exact = 1;
    numpts = numpts_;
    this->tps_init (cdata__, data__, type, rbf_parm, integrateme, addlinearterms, transmode);
    this->init_exact (x, y, z, xstride, ystride, zstride, numpts_, addlinearterms);
}

pixsrc_tps::pixsrc_tps(commoninputdata *cdata__, inputdata *data__, const double *x, const double *y, const double *z, PS_SIT xstride, PS_SIT ystride, PS_SIT zstride, PS_SIT numtotalpts_, bool *is_control, ps_rbf_type type, double rbf_parm, PS_SIT integrateme, PS_SIT addlinearterms, PS_SIT transmode)
{
    // this constructor is for creating an approximate TPS
    // only a subset of the points will be used as TPS basis functions
    // but all the points will be used to constrain the TPS weights
    // if optpts is true, then the z-values that minimize the chi^2 (with 0 regularization)
    // will be solved for
    // currently , it's only known to work when called from pixsrc_init_uvplane -- i.e. when z=NULL
    // though itmay work if z!=NULL

    // count how many points are control points
    numtotalpts = numtotalpts_;
    numpts = 0;
    if (!is_control)
        numpts = numtotalpts;
    else
        for (PS_SIT i=0; i<numtotalpts; ++i)
            if (is_control[i])
                ++numpts;
    if (numpts==numtotalpts)
    {
        is_exact = 1;
        this->tps_init (cdata__, data__, type, rbf_parm, integrateme, addlinearterms, transmode);
        this->init_exact (x, y, z, xstride, ystride, zstride, numpts, addlinearterms);
        return;
    }

    is_exact = 0;
    this->tps_init (cdata__, data__, type, rbf_parm, integrateme, addlinearterms, transmode);

    MEMORY ps_malloc (&tpsmatptr, 1);
    tpsmat = new (tpsmatptr) MATRIX (cdata_, data_, dim, dim, -1, -1, data_);
    double *tps_mat = tpsmat->mat_dense;
    MEMORY ps_malloc (&tmatptr, 1);
    tmat = new (tmatptr) MATRIX (cdata_, data_, dim, numtotalpts, -1, -1, data_);
    double *t_mat = tmat->mat_dense;

    MEMORY ps_malloc (&xcoords, (PS_SIT)numpts);
    MEMORY ps_malloc (&ycoords, (PS_SIT)numpts);

    // get characteristic length scale of control points
    // and fill in K matrix
    meandist2 = 0;
    PS_SIT indexA, indexB, indexN=-1, index;
    double r2;
    for (PS_SIT m=0; m<numtotalpts; ++m)
    {
        // track column positions
        if (is_control && !is_control[m])
            continue;
        ++indexN;

        indexA = indexB = 0;
        for (PS_SIT n=0; n<numtotalpts; ++n)
        {
            // track row positions
            if (!is_control || is_control[n])
                index = indexA++;
            else
                index = numpts + indexB++;

            r2 = OPERA distance2 (x[m*xstride], y[m*ystride], x[n*xstride], y[n*ystride]);
            double limits[4] = {x[m*xstride]-x[n*xstride]-0.5,
                                x[m*xstride]-x[n*xstride]+0.5,
                                y[m*ystride]-y[n*ystride]-0.5,
                                y[m*ystride]-y[n*ystride]+0.5};
            if (integrate_rbf)
                tmat->set (indexN, index, (double)this->kernel_integral (limits));
            else
                tmat->set (indexN, index, (double)this->kernel (r2));

            if (m!=n && (!is_control || is_control[n]))
                meandist2 += std::sqrt (r2);
        }
    }
    meandist2 /= (numpts*numpts-numpts);
    meandist2 *= meandist2;

    // fill in P matrix into R matrix
    if(-1!=uselinearterms)
        for (PS_SIT m=0; m<numtotalpts; ++m)
            tmat->set (numpts, m, 1);
    index=0;
    // get control points
    for (PS_SIT m=0; m<numtotalpts; ++m)
        if (!is_control || is_control[m])
        {
            if (1==uselinearterms)
            {
                tmat->set (numpts+1, index, (double)x[m*xstride]);
                tmat->set (numpts+2, index, (double)y[m*ystride]);
            }
            xcoords[index]             = (double)x[m*xstride];
            ycoords[index]             = (double)y[m*ystride];
            ++index;
        }
    // now get the other points
    for (PS_SIT m=0; m<numtotalpts; ++m)
        if (is_control && !is_control[m])
        {
            if (1==uselinearterms)
            {
                tmat->set (numpts+1, index, (double)x[m*xstride]);
                tmat->set (numpts+2, index, (double)y[m*ystride]);
            }
            ++index;
        }

    // constants for BLAS routines
    char t     = 'T';
    char n     = 'N';
    PS_SIT onei   = 1;
    double one      = 1.0;
    double zer      = 0.0;
    PS_SIT tnumtotalpts = (PS_SIT)numtotalpts;
    PS_SIT tdim = (PS_SIT)dim;
    // compute R^t.R = TPS-MATRIX
    my_pstpsgemm (&n, &t, &tdim, &tdim, &tnumtotalpts,
                  &one, t_mat, &tdim,
                  t_mat, &tdim,
                  &zer, tps_mat, &tdim);

    if (z)
    {
        double *z_vec;
        MEMORY ps_malloc (&a_mat,   (PS_SIT)numpts*numpts);
        MEMORY ps_malloc (&b_vec,   (PS_SIT)dim);
        MEMORY ps_malloc (&z_vec,   (PS_SIT)numtotalpts);

        // copy top square part of K matrix into A matrix
        for (PS_SIT m=0; m<numpts; ++m)
            for (PS_SIT n=0; n<numpts; ++n)
                a_mat[m*numpts+n] = tmat->get(m,n);

        // store z-values of all points (control pts first, then others)
        index=0;
        for (PS_SIT i=0; i<numtotalpts; ++i)
            if (!is_control || is_control[i])
                z_vec[index++] = (double)z[i*zstride];
        for (PS_SIT i=0; i<numtotalpts; ++i)
            if (is_control && !is_control[i])
                z_vec[index++] = (double)z[i*zstride];

        // calculate b vector = R^t.v
        my_pstpsgemm (&n, &n, &tdim, &onei, &tnumtotalpts,
                      &one, t_mat, &tdim,
                      z_vec, &tnumtotalpts,
                      &zer, b_vec, &tdim);

        // cleanup
        MEMORY ps_free (z_vec);

        // allocate memory that will be needed
        MEMORY ps_malloc (&tps_mat_fac, (PS_SIT)dim*dim);
        MEMORY ps_malloc (&x_vec,       (PS_SIT)dim);
        MEMORY ps_malloc (&ipiv,        (PS_SIT)(dim+1));
    }

    if (!z)
    {
        PS_SIT numextraRBF;
        if (-1==uselinearterms)
            numextraRBF = 0;
        else if (0==uselinearterms)
            numextraRBF = 1;
        else
            numextraRBF = 3;
        // get LU decomposition
        tpsmat->inplace_lu_fac ();
        // cut off non-control points
        tmat->ncol = numpts+numextraRBF;
        tmat->nnz = tmat->nrow*tmat->ncol;
        MEMORY ps_realloc (&tmat->mat_dense, tmat->nnz);
        t_mat = tmat->mat_dense;
        if (0==uselinearterms)
        {
            std::fill (t_mat+numpts*dim, t_mat+numpts*dim+numpts, 1);
            t_mat[dim*dim-1] = 0;
        }
        // compute (k^t k)^-1 k^t
        tpsmat->inplace_leq (tmat);
        // uv plane code assumes that tpsmat is the correct matrix
        MATRIX *dummy=tpsmat, *dummyptr=tpsmatptr;
        tpsmat    = tmat;
        tpsmatptr = tmatptr;
        tmat    = dummy;
        tmatptr = dummyptr;


        // some more processing
        // get inverse of TPS matrix and multiply by TPS R-matrix -- (R^t.R)^-1.R^t
        if (0==transmode)
        {
            if (!data_->uv_padzeros)
                this->tpsmat->inplace_inverse ();
            // need to chop off one column to throw away constant term in RBF
            this->tpsmat->ncol -= numextraRBF;
            this->tpsmat->nnz   = this->tpsmat->ncol*this->tpsmat->nrow;
        }
        else if (!data_->uv_padzeros)
        {
            this->tpsmat->inplace_lu_fac ();
        }


        MEMORY ps_malloc (&x_rad_offset, (PS_SIT)numpts);
        MEMORY ps_malloc (&y_rad_offset, (PS_SIT)numpts);
        pthread_mutex_lock (cdata_->wcsmutex);
        for (PS_SIT b=0; b<numpts; ++b)
        {
            // get arcsecond offsets
            double ra0, dec0, pos[2];
            HEADER getimgwcscoord  (data_->wcs, data_->imgy, xcoords[b], ycoords[b], &ra0, &dec0);
            HEADER getwcslfromwcss (ra0, dec0, data_->pra, data_->pdec,
                                    data_->r1, data_->r2, pos);
            // arcseconds to radians
            pos[0] *= CONSTANT arc2rad;
            pos[1] *= CONSTANT arc2rad;
            // account for R.A. increasing leftwards
            pos[0] *= -1;
            // I don't know why I have to do the following negation, but it works
            // It could fail if pointing center/coordinate system are not at (0,0)
            pos[1] *= -1;
            x_rad_offset[b] = pos[0];
            y_rad_offset[b] = pos[1];
        }
        pthread_mutex_unlock (cdata_->wcsmutex);
    }

    // cleanup
    tmat->~MATRIX();
    MEMORY ps_free (tmatptr);
    tmat    = NULL;
    tmatptr = NULL;
}

void pixsrc_tps::init_exact (const double *x, const double *y, const double *z, PS_SIT xstride, PS_SIT ystride, PS_SIT zstride, PS_SIT numpts_, PS_SIT addlinearterms)
{
    MEMORY ps_malloc (&tpsmatptr, 1);
    tpsmat = new (tpsmatptr) MATRIX (cdata_, data_, dim, dim, -1, -1, data_);
    double *tps_mat = tpsmat->mat_dense;

    // get characteristic length scale of control points
    // and fill in TPS matrix
    meandist2 = 0;
    PS_SIT index = 0;
    double r2;
    for (PS_SIT m=0; m<numpts; ++m)
    {
        for (PS_SIT n=0; n<numpts; ++n)
        {
            r2 = OPERA distance2 (x[m*xstride], y[m*ystride], x[n*xstride], y[n*ystride]);
            double limits[4] = {x[m*xstride]-x[n*xstride]-0.5,
                                x[m*xstride]-x[n*xstride]+0.5,
                                y[m*ystride]-y[n*ystride]-0.5,
                                y[m*ystride]-y[n*ystride]+0.5};
            //std::cout << m << " " << n << " " << index << " " << dim*dim << std::endl;
            if (n<=m)
            {
                if (integrate_rbf)
                    tps_mat[index] = (double)this->kernel_integral (limits);
                else
                    tps_mat[index] = (double)this->kernel (r2);
            }
            //if (n==0)
            //  std::cout << limits[0] << " " << limits[1] << " " << limits[2] << " "
            //<< limits[3] << " " << tps_mat[index] << std::endl;
            ++index;

            if (m!=n)
                meandist2 += std::sqrt (r2);
        }

        if (0==uselinearterms)
        {
            tps_mat[index++] = (double)1;
        }
        else if (1==uselinearterms)
        {
            tps_mat[index++] = (double)1;
            tps_mat[index++] = (double)x[m*xstride];
            tps_mat[index++] = (double)y[m*ystride];
        }
    }

    // copy other half of top left block
    for (PS_SIT m=0; m<numpts; ++m)
        for (PS_SIT n=0; n<numpts; ++n)
            if (n>m)
                tps_mat[m*dim+n] = tps_mat[n*dim+m];

    if (-1==uselinearterms)
    {}
    else if (0==uselinearterms)
    {
        std::fill (tps_mat+index, tps_mat+index+numpts, (double)1);
        index += numpts;
        std::fill (tps_mat+index, tps_mat+index+1, (double)0);
        index += 1;
    }
    else
    {
        std::fill (tps_mat+index, tps_mat+index+numpts, (double)1);
        index += numpts;
        std::fill (tps_mat+index, tps_mat+index+3, (double)0);
        index += 3;
        for (PS_SIT n=0; n<numpts; ++n)
        {
            tps_mat[index++] = (double)x[n*xstride];
        }
        std::fill (tps_mat+index, tps_mat+index+3, (double)0);
        index += 3;
        for (PS_SIT n=0; n<numpts; ++n)
        {
            tps_mat[index++] = (double)y[n*ystride];
        }
        std::fill (tps_mat+index, tps_mat+index+3, (double)0);
        index +=3;
    }
    meandist2 /= (numpts*numpts-numpts);
    meandist2 *= meandist2;

    if (z)
    {
        MEMORY ps_malloc (&b_vec, (PS_SIT)dim);
        for (PS_SIT m=0; m<numpts; ++m)
            b_vec[m] = (double)z[m*zstride];

        if (-1==uselinearterms)
        {}
        else if (0==uselinearterms)
            b_vec[dim-1] = 0;
        else
            b_vec[dim-1] = b_vec[dim-2] = b_vec[dim-3] = (double)0;

        MEMORY ps_malloc (&tps_mat_fac, (PS_SIT)dim*dim);
        MEMORY ps_malloc (&x_vec,       (PS_SIT)dim);
        MEMORY ps_malloc (&ipiv,        (PS_SIT)(dim+1));
    }

    if (uselinearterms<=0)
    {
        MEMORY ps_malloc (&xcoords, (PS_SIT)numpts);
        MEMORY ps_malloc (&ycoords, (PS_SIT)numpts);
        for (PS_SIT b=0; b<numpts; ++b)
        {
            xcoords[b] = x[b*xstride];
            ycoords[b] = y[b*ystride];
        }
    }

    if (!z)
    {
        PS_SIT numextraRBF;
        if (-1==uselinearterms)
            numextraRBF = 0;
        else if (0==uselinearterms)
            numextraRBF = 1;
        else
            numextraRBF = 3;
        // some more processing
        // get inverse of TPS matrix and multiply by TPS R-matrix -- (R^t.R)^-1.R^t
        if (0==transmode)
        {
            if (!data_->uv_padzeros)
                this->tpsmat->inplace_inverse ();
            // need to chop off one column to throw away constant term in RBF
            this->tpsmat->ncol -= numextraRBF;
            this->tpsmat->nnz   = this->tpsmat->ncol*this->tpsmat->nrow;
        }
        else if (!data_->uv_padzeros)
        {
            this->tpsmat->inplace_lu_fac ();
        }

        MEMORY ps_malloc (&x_rad_offset, (PS_SIT)numpts);
        MEMORY ps_malloc (&y_rad_offset, (PS_SIT)numpts);
        for (PS_SIT b=0; b<numpts; ++b)
        {
            // get arcsecond offsets
            // I should be putting a mutex lock around this b/c it's not thread-safe,
            // but it's not necessary yet
            double ra0, dec0, pos[2];
            HEADER getimgwcscoord  (data_->wcs, data_->imgy, x[b*xstride], y[b*ystride], &ra0, &dec0);
            HEADER getwcslfromwcss (ra0, dec0, data_->pra, data_->pdec,
                                    data_->r1, data_->r2, pos);
            // arcseconds to radians
            pos[0] *= CONSTANT arc2rad;
            pos[1] *= CONSTANT arc2rad;
            // account for R.A. increasing leftwards
            pos[0] *= -1;
            // I don't know why I have to do the following negation, but it works
            // It could fail if pointing center/coordinate system are not at (0,0)
            pos[1] *= -1;
            x_rad_offset[b] = pos[0];
            y_rad_offset[b] = pos[1];
        }
    }

    /*
      for (PS_SIT m=0; m<dim; ++m)
      for (PS_SIT n=0; n<dim; ++n)
      std::cout << "MAT " << m << " " << n << " " << tps_mat[m*dim+n] << std::endl;
    */
}

void pixsrc_tps::tps_init(commoninputdata *cdata__, inputdata *data__, ps_rbf_type type, double rbf_parm, PS_SIT integrateme, PS_SIT addlinearterms, PS_SIT transmode_)
{
    transmode = transmode_;
    reg         = 0;
    meandist2   = 0;
    ipiv        = NULL;
    tps_mat_fac = NULL;
    x_vec       = NULL;
    b_vec       = NULL;
    xcoords     = NULL;
    ycoords     = NULL;
    x_rad_offset= NULL;
    y_rad_offset= NULL;
    a_mat       = NULL;
    this->set_precision (tps_mat_fac);
    tpsmatptr = NULL;
    tpsmat = NULL;
    tmatptr = NULL;
    tmat = NULL;

    cdata_ = cdata__;
    data_ = data__;
    psrbftype = type;
    integrate_rbf = integrateme;
    uselinearterms = addlinearterms;
    gsl_workspace_tps_x = (void*)gsl_integration_workspace_alloc (10000);
    gsl_workspace_tps_y = (void*)gsl_integration_workspace_alloc (10000);

    // Fourier transform stuff
    gauss_sigma    = rbf_parm;
    ft_gauss_sigma = rbf_parm *data_->pix2arc *CONSTANT arc2rad;
    dxdyfac        = data_->pix2arc *CONSTANT arc2rad;
    dxdyfac        = 1.0 /(dxdyfac*dxdyfac);
    if (ps_rbf_gaussian==psrbftype)
        fac1 = CONSTANT twopi *ft_gauss_sigma*ft_gauss_sigma;
    else if (ps_rbf_linear_weighted==psrbftype)
        fac1 = std::sqrt (2.) *std::pow(CONSTANT pi,1.5)
            *ft_gauss_sigma*ft_gauss_sigma*ft_gauss_sigma;

    if (-1==uselinearterms)
        dim = numpts;
    else if (0==uselinearterms)
        dim = numpts+1;
    else
        dim = numpts+3;
}
/*
  void pixsrc_tps::set_precision (float *dummy)
  {
  // set single precision BLAS functions
  my_pstpsgemm  = EXTERNAL ps_sgemm_;
  my_pstpssytrf = EXTERNAL ps_ssytrf_;
  my_pstpssytrs = EXTERNAL ps_ssytrs_;
  }
*/
void pixsrc_tps::set_precision (double *dummy)
{
    // set double precision BLAS functions
    my_pstpsgemm  = EXTERNAL ps_dgemm_;
    my_pstpssytrf = EXTERNAL ps_dsytrf_;
    my_pstpssytrs = EXTERNAL ps_dsytrs_;
}

void pixsrc_tps::set_regularization (double reg_)
{
    reg = reg_;
}

PS_SIT pixsrc_tps::get_num_ctrl_pts ()
{
    return numpts;
}

void pixsrc_tps::get_tps_weights ()
{
    double *tps_mat = tpsmat->mat_dense;
    // copy tps matrix and set regularization strengths
    std::copy (tps_mat, tps_mat+dim*dim, tps_mat_fac);
    // copy vector into which solutions will be stored
    std::copy (b_vec, b_vec+dim, x_vec);

    // add regularization terms
    if (is_exact)
    {
        if (reg)
            for (PS_SIT m=0; m<numpts; ++m)
                tps_mat_fac[m*dim+m] += (double)(meandist2*reg);
    }
    else
    {
        // should use BLAS (LDA set properly) to do this addition instead
        if (reg)
            for (PS_SIT m=0; m<numpts; ++m)
                for (PS_SIT n=0; n<numpts; ++n)
                    tps_mat_fac[m*dim+n] += (double)(meandist2*reg)*a_mat[m*numpts+n];
    }

    PS_SIT info, lwork=-1, nrhs=1;
    PS_SIT LDA=dim, LDB=dim;
    char UPLO = 'U';
    double *work, wkopt;
    PS_SIT tdim = (PS_SIT)dim;

    // estimate work space needed
    my_pstpssytrf (&UPLO, &tdim, tps_mat_fac, &LDA,
                   ipiv, &wkopt, &lwork, &info);

    lwork = (PS_SIT)wkopt;
    MEMORY ps_malloc (&work, (PS_SIT)lwork);

    // do LU factorization of tps matrix
    my_pstpssytrf (&UPLO, &tdim, tps_mat_fac, &LDA,
                   ipiv, work, &lwork, &info);

    MEMORY ps_free (work);

    // solve for tps weights
    my_pstpssytrs (&UPLO, &tdim, &nrhs, tps_mat_fac,
                   &LDA, ipiv, x_vec, &LDB, &info);

    //std::cout << "0th " << x_vec[dim-1] << std::endl;

    //double *tmp = (double*) malloc (dim*sizeof(double));
    //char n='N';double a=1; PS_SIT onei=1; double moi=-1;
    //EXTERNAL ps_dgemv_ (&n,&tdim,&tdim,&a,tps_mat,&tdim,x_vec,&onei,&moi,b_vec,&onei);
    //for (PS_SIT r=0; r<tdim; ++r)
    //std::cout << "res: " << b_vec[r] << std::endl;
    //exit(0);
}

void pixsrc_tps::cleanup ()
{
    // if LU factorization is complete, then you can delete large arrays
    // but we need to save TPS control point positions
    double *tps_mat = tpsmat->mat_dense;
    if (!xcoords && uselinearterms>0)
    {
        MEMORY ps_malloc (&xcoords, (PS_SIT)numpts);
        std::copy (tps_mat+dim*(numpts+1), tps_mat+dim*(numpts+1)+numpts, xcoords);
    }
    if (!ycoords && uselinearterms>0)
    {
        MEMORY ps_malloc (&ycoords, (PS_SIT)numpts);
        std::copy (tps_mat+dim*(numpts+2), tps_mat+dim*(numpts+2)+numpts, ycoords);
    }

    if (tpsmat)
    {
        tpsmat-> ~MATRIX ();
        tpsmat  = NULL;
        MEMORY ps_free (tpsmatptr);
    }
    if (tmat)
    {
        tmat-> ~MATRIX ();
        tmat  = NULL;
        MEMORY ps_free (tmatptr);
    }

    MEMORY ps_free (ipiv);
    MEMORY ps_free (b_vec);
    MEMORY ps_free (tps_mat_fac);
    MEMORY ps_free (a_mat);
    ipiv        = NULL;
    b_vec       = NULL;
    tps_mat_fac = NULL;
    a_mat       = NULL;
}

double pixsrc_tps::interpolate (double x, double y)
{
    double r2, zval=0;
    double *gridx, *gridy;
    double *tps_mat = tpsmat->mat_dense;
    // get pointer to control point coordinates
    if (tps_mat && is_exact && uselinearterms>0)
    {
        gridx = &tps_mat[dim*(numpts+1)];
        gridy = &tps_mat[dim*(numpts+2)];
    }
    else
    {
        gridx = xcoords;
        gridy = ycoords;
    }

    if (-1==uselinearterms)
    {}
    else if (0==uselinearterms)
        zval += x_vec[dim-1];
    else
        zval += x_vec[dim-3] + x_vec[dim-2]*x + x_vec[dim-1]*y;

    for (PS_SIT v=0; v<numpts; ++v)
    {
        r2 = OPERA distance2 (x, y, gridx[v], gridy[v]);
        zval += x_vec[v] *this->kernel (r2);
    }

    return zval;
}

void pixsrc_tps::precompute_fourier (double *pos, PS_SIT ndp, double *fac2)
{
    if (ps_rbf_gaussian==psrbftype)
    {
        for (PS_SIT uv=0; uv<ndp; ++uv)
        {
            double uv2 = pos[uv*2]*pos[uv*2] + pos[uv*2+1]*pos[uv*2+1];
            fac2[uv] = std::exp (-CONSTANT pi *fac1 *uv2);
        }
    }
    else if (ps_rbf_linear_weighted==psrbftype)
    {
        for (PS_SIT uv=0; uv<ndp; ++uv)
        {
            double uv2 = pos[uv*2]*pos[uv*2] + pos[uv*2+1]*pos[uv*2+1];
            double tmpfac = CONSTANT pi*CONSTANT pi*ft_gauss_sigma*ft_gauss_sigma *uv2;
            fac2[uv] = std::exp (-tmpfac) *(
                (1.0-2.0*tmpfac)*gsl_sf_bessel_I0(tmpfac) +
                2.0*tmpfac*gsl_sf_bessel_I1(tmpfac)
                );
        }
    }
}

void pixsrc_tps::get_fourier_weights (double *pos, MATRIX *weights, PS_SIT ndp)
{
    // get weights for interpolation in Fourier space
    double wts[2];
    for (PS_SIT d=0; d<ndp; ++d)
    {
        for (PS_SIT v=0; v<numpts; ++v)
        {
            this->fourier_kernel (x_rad_offset[v], y_rad_offset[v],
                                  pos[d*2+0], pos[d*2+1], wts, NULL);
            weights->set (v, d*2+0, wts[0]);
            weights->set (v, d*2+1, wts[1]);
        }
        if (-1!=uselinearterms)
        {
            weights->set (numpts, d*2+0, 0);
            weights->set (numpts, d*2+1, 0);
        }
    }
}

void pixsrc_tps::get_fourier_weights_col (PS_SIT col, double *pos, double *weights, PS_SIT ndp, double *fac2, PS_SIT start, PS_SIT end)
{
    // get weights for interpolation in Fourier space
    double wts[2];
    for (PS_SIT d=start; d<end; ++d)
    {
        if (col!=numpts)
        {
            this->fourier_kernel (x_rad_offset[col], y_rad_offset[col],
                                  pos[d*2+0], pos[d*2+1], wts, &fac2[d]);
            weights[(d-start)*2+0] = wts[0];
            weights[(d-start)*2+1] = wts[1];
        }
        else
        {
            weights[(d-start)*2+0] = weights[(d-start)*2+1] = 0;
        }
    }
}

void pixsrc_tps::interpolate (double *pos, PS_SIT *r4r, VECTOR *img)
{
    PS_SIT ptrack=-1;
    for (PS_SIT r=0; r<img->get_size(); ++r)
    {
        ++ptrack;
        // find position of first good pixel
        while (-1==r4r[ptrack])
            ++ptrack;
        img->set (r, this->interpolate(pos[ptrack*2], pos[ptrack*2+1]));
    }
}

void pixsrc_tps::interpolate (double **pos, PS_SIT *r4r, VECTOR *img, PS_SIT subsample)
{
    PS_SIT ptrack=-1;
    for (PS_SIT r=0; r<img->get_size(); ++r)
    {
        ++ptrack;
        // find position of first good pixel
        while (-1==r4r[ptrack])
            ++ptrack;

        // subsample pixel
        double val=0;
        for (PS_SIT s=0; s<subsample*subsample; ++s)
            val += this->interpolate(pos[ptrack][s*2], pos[ptrack][s*2+1]);
        img->set (r, val);
    }
    img->mult (1.0/(subsample*subsample));
}

double pixsrc_tps::kernel (double r2)
{
    if (ps_rbf_gaussian==psrbftype)
        return std::exp(-0.5*r2/(gauss_sigma*gauss_sigma));
    else if (ps_rbf_tps_weighted==psrbftype)
        return OPERA equalszero(r2) ? 0 : r2 *0.5 *std::log(r2) *std::exp(-0.5*r2/(gauss_sigma*gauss_sigma));
    else if (ps_rbf_linear_weighted==psrbftype)
        return std::sqrt(r2) *std::exp(-0.5*r2/(gauss_sigma*gauss_sigma));
    else if (ps_rbf_tps==psrbftype)
        return OPERA equalszero(r2) ? 0 : r2 *0.5 *std::log(r2);

    // gug: throw error here
    return 0;
    //return std::sqrt(r2)*std::exp(-0.5*r2/(gauss_sigma*gauss_sigma));
}

void pixsrc_tps::fourier_kernel (double x, double y, double u, double v, double *wts, double *fac2_)
{
    if (ps_rbf_tps==psrbftype)
        PRINTER printerror (data_->print2screenname,
                            "cannot Fourier transform TPS -- try weighted TPS.",
                            cdata_->print2screenmutex);

    if (ps_rbf_gaussian==psrbftype)
    {
        double uv2  = u*u + v*v;
        double fac2;
        if (fac2_)
            fac2 = *fac2_;
        else
            fac2 = std::exp (-CONSTANT pi *fac1 *uv2);
        // get real and imaginary parts
        double angle = CONSTANT twopi *(u*x + v*y);
        wts[0] =  fac1 *fac2 *dxdyfac *std::cos (angle);
        wts[1] = -fac1 *fac2 *dxdyfac *std::sin (angle);
    }
    else if (ps_rbf_linear_weighted==psrbftype)
    {
        double uv2  = u*u + v*v;
        double fac2;
        if (fac2_)
        {
            fac2 = *fac2_;
        }
        else
        {
            double tmpfac = CONSTANT pi*CONSTANT pi*ft_gauss_sigma*ft_gauss_sigma *uv2;
            fac2 = std::exp (-tmpfac) *(
                (1.0-2.0*tmpfac)*gsl_sf_bessel_I0(tmpfac) +
                2.0*tmpfac*gsl_sf_bessel_I1(tmpfac)
                );
        }
        // get real and imaginary parts
        double angle = CONSTANT twopi *(u*x + v*y);
        wts[0] =  fac1 *fac2 *dxdyfac *std::cos (angle);
        wts[1] = -fac1 *fac2 *dxdyfac *std::sin (angle);
    }
}

struct ps_rbf_integration_struct
{
    pixsrc_tps *thisclass;
    double x;
    double *lim;
};
double pixsrc_tps::ps_rbf_integrate_y (double y, void *parms)
{
    ps_rbf_integration_struct *pris = (ps_rbf_integration_struct*)parms;
    double r2 = pris->x*pris->x + y*y;
    return pris->thisclass->kernel (r2);
}

double pixsrc_tps::ps_rbf_integrate_x (double x, void *parms)
{
    ps_rbf_integration_struct *pris = (ps_rbf_integration_struct*)parms;
    pris->x = x;
    double integral, interr;
    gsl_integration_workspace *wly = (gsl_integration_workspace*)pris->thisclass->gsl_workspace_tps_y;
    gsl_function F1;
    F1.function = &pixsrc_tps::ps_rbf_integrate_y;
    F1.params = parms;
    gsl_integration_qags (&F1, pris->lim[2], pris->lim[3], 0, 1.0e-3, 10000, wly, &integral, &interr);
    return integral;
}

PS_SIT pixsrc_tps::integrate_quadrant_query (PS_SIT quadrant, double *lim)
{
    // move this quadrant to first quadrant
    switch (quadrant)
    {
    case 0:
        break;
    case 1:
        lim[0] *= -1;
        lim[1] *= -1;
        OPERA swap (&lim[0], &lim[1]);
        break;
    case 2:
        lim[0] *= -1;
        lim[1] *= -1;
        OPERA swap (&lim[0], &lim[1]);
        lim[2] *= -1;
        lim[3] *= -1;
        OPERA swap (&lim[2], &lim[3]);
        break;
    case 3:
        lim[2] *= -1;
        lim[3] *= -1;
        OPERA swap (&lim[2], &lim[3]);
        break;
    default:
        break;
    }
    // cut off parts of square not in first quadrant
    lim[0] = std::max (0.0, lim[0]);
    lim[2] = std::max (0.0, lim[2]);
    // if no overlap of square with this quadrant
    if (lim[1]<=lim[0] || lim[3]<=lim[2])
        return 0;
    else
        return 1;
}

double tps_quadrant_integral (double *lim)
{
    double integral;
    double pi  = CONSTANT pi;
    double x1  = lim[0];
    double x2  = lim[1];
    double y1  = lim[2];
    double y2  = lim[3];
    double x12 = x1*x1;
    double x22 = x2*x2;
    double y12 = y1*y1;
    double y22 = y2*y2;
    //double y23 = y22*y2;
    double x14 = x12*x12;
    double x24 = x22*x22;
    double y14 = y12*y12;
    double y24 = y22*y22;
    if (OPERA equalszero (x1) && OPERA equalszero (y1))
        integral = 1./18.*(3.*pi*x24-10.*x2*y2*(x22+y22)
                           +6.*(-x24+y24)*std::atan(x2/y2)
                           +6.*x2*y2*(x22+y22)*std::log(x22+y22));
    else if (OPERA equalszero (x1))
        integral = 1./9.*(3.*(x24-y14)*std::atan(x2/y1)
                          +3.*(-x24+y24)*std::atan(x2/y2)
                          +x2*(5.*(y1-y2)*(x22+y12+y1*y2+y22)
                               -3.*y1*(x22+y12)*std::log(x22+y12)
                               +3.*y2*(x22+y22)*std::log(x22+y22)));
    else if (OPERA equalszero (y1))
        integral = 1./18.*(3.*pi*(-x14+x24)
                           +10.*(x1-x2)*y2*(x12+x1*x2+x22+y22)
                           +6.*(x14-y24)*std::atan(x1/y2)
                           +6.*(-x24+y24)*std::atan(x2/y2)
                           -6.*x1*y2*(x12+y22)*std::log(x12+y22)
                           +6.*x2*y2*(x22+y22)*std::log(x22+y22));
    else
        integral = 1./9.*(-(x12+y12)*(3.*(x1-y1)*(x1+y1)
                                      *std::atan(x1/y1)+x1*y1*(5.-3.*std::log(x12+y12)))
                          +(x22+y12)*(3.*(x2-y1)*(x2+y1)
                                      *std::atan(x2/y1)+x2*y1*(5.-3.*std::log(x22+y12)))
                          +(x12+y22)*(3.*(x1-y2)*(x1+y2)
                                      *std::atan(x1/y2)+x1*y2*(5.-3.*std::log(x12+y22)))
                          -(x22+y22)*(3.*(x2-y2)*(x2+y2)
                                      *std::atan(x2/y2)+x2*y2*(5.-3.*std::log(x22+y22))));

    return 0.5*integral;
}

double linear_quadrant_integral (double *lim)
{
    double integral;
    //double pi  = CONSTANT pi;
    double x1  = lim[0];
    double x2  = lim[1];
    double y1  = lim[2];
    double y2  = lim[3];
    double x12 = x1*x1;
    double x13 = x12*x1;
    double x22 = x2*x2;
    double x23 = x22*x2;
    double y12 = y1*y1;
    double y13 = y12*y1;
    double y22 = y2*y2;
    double y23 = y22*y2;
    //double x14 = x12*x12;
    //double x24 = x22*x22;
    //double y14 = y12*y12;
    //double y24 = y22*y22;

    if (OPERA equalszero (x1) && OPERA equalszero (y1))
        integral = 1./6.*(2.*x2*y2*std::sqrt(x22+y22)
                          +y23*std::log((x2+std::sqrt(x22+y22))/y2)
                          +x23*std::log((y2+std::sqrt(x22+y22))/x2));
    else
        integral = 1./6.*(2.*x1*y1*std::sqrt(x12+y12)
                          -2.*x2*y1*std::sqrt(x22+y12)
                          -2.*x1*y2*std::sqrt(x12+y22)
                          +2.*x2*y2*std::sqrt(x22+y22)
                          +y13*std::log((x1+std::sqrt(x12+y12))
                                        /(x2+std::sqrt(x22+y12)))
                          +x13*std::log((y1+std::sqrt(x12+y12))
                                        /(y2+std::sqrt(x12+y22)))
                          +y23*std::log((x2+std::sqrt(x22+y22))
                                        /(x1+std::sqrt(x12+y22)))
                          +x23*std::log((y2+std::sqrt(x22+y22))
                                        /(y1+std::sqrt(x22+y12))));

    return integral;
}

double pixsrc_tps::kernel_integral (double *lim)
{
    if (ps_rbf_gaussian==psrbftype)
    {
        double erffac = 1.0 / (gauss_sigma *std::sqrt (2.0));
        return CONSTANT piby2 *gauss_sigma*gauss_sigma
            *(gsl_sf_erf(erffac*lim[0])-gsl_sf_erf(erffac*lim[1]))
            *(gsl_sf_erf(erffac*lim[2])-gsl_sf_erf(erffac*lim[3]));
    }
    else if (ps_rbf_tps_weighted==psrbftype)
    {
        double minx = std::min (std::abs(lim[0]),std::abs(lim[1]));
        double miny = std::min (std::abs(lim[2]),std::abs(lim[3]));
        double r2min = (minx*minx+miny*miny) /(gauss_sigma*gauss_sigma);
        if (r2min>5.0*5.0)
        {
            double ctrx, ctry, integral=0, limcopy[4];
            for (PS_SIT q=0; q<4; ++q)
            {
                std::copy (lim, lim+4, limcopy);
                if (integrate_quadrant_query(q, limcopy))
                {
                    ctrx = 0.5 *(limcopy[0]+limcopy[1]);
                    ctry = 0.5 *(limcopy[2]+limcopy[3]);
                    integral += tps_quadrant_integral (limcopy)
                        *std::exp(-0.5*(ctrx*ctrx+ctry*ctry)/(gauss_sigma*gauss_sigma));
                }
            }
            return integral;
        }

        ps_rbf_integration_struct pris;
        pris.thisclass = this;
        pris.lim = lim;
        double integral, interr;
        gsl_integration_workspace *wlx = (gsl_integration_workspace*)gsl_workspace_tps_x;
        gsl_function F1;
        F1.function = &pixsrc_tps::ps_rbf_integrate_x;
        F1.params = (void*)&pris;
        gsl_integration_qags (&F1, lim[0], lim[1], 0, 1.0e-3, 10000, wlx, &integral, &interr);
        return integral;
    }
    else if (ps_rbf_linear_weighted==psrbftype)
    {
        double minx = std::min (std::abs(lim[0]),std::abs(lim[1]));
        double miny = std::min (std::abs(lim[2]),std::abs(lim[3]));
        double r2min = (minx*minx+miny*miny) /(gauss_sigma*gauss_sigma);
        if (r2min>5.0*5.0)
        {
            double ctrx, ctry, integral=0, limcopy[4];
            for (PS_SIT q=0; q<4; ++q)
            {
                std::copy (lim, lim+4, limcopy);
                if (integrate_quadrant_query(q, limcopy))
                {
                    ctrx = 0.5 *(limcopy[0]+limcopy[1]);
                    ctry = 0.5 *(limcopy[2]+limcopy[3]);
                    integral += linear_quadrant_integral (limcopy)
                        *std::exp(-0.5*(ctrx*ctrx+ctry*ctry)/(gauss_sigma*gauss_sigma));
                }
            }
            return integral;
        }

        ps_rbf_integration_struct pris;
        pris.thisclass = this;
        pris.lim = lim;
        double integral, interr;
        gsl_integration_workspace *wlx = (gsl_integration_workspace*)gsl_workspace_tps_x;
        gsl_function F1;
        F1.function = &pixsrc_tps::ps_rbf_integrate_x;
        F1.params = (void*)&pris;
        gsl_integration_qags (&F1, lim[0], lim[1], 0, 1.0e-3, 10000, wlx, &integral, &interr);
        return integral;
    }
    else if (ps_rbf_tps==psrbftype)
    {
        double integral=0, limcopy[4];
        for (PS_SIT q=0; q<4; ++q)
        {
            std::copy (lim, lim+4, limcopy);
            if (integrate_quadrant_query(q, limcopy))
                integral += tps_quadrant_integral (limcopy);
        }
        return integral;
    }
    return 0;
}

pixsrc_tps::~pixsrc_tps ()
{
    MEMORY ps_free (ipiv);
    MEMORY ps_free (b_vec);
    MEMORY ps_free (x_vec);
    MEMORY ps_free (xcoords);
    MEMORY ps_free (ycoords);
    MEMORY ps_free (x_rad_offset);
    MEMORY ps_free (y_rad_offset);
    MEMORY ps_free (tps_mat_fac);
    MEMORY ps_free (a_mat);
    if (tpsmat)
    {
        tpsmat->~MATRIX();
        MEMORY ps_free (tpsmatptr);
    }
    if (tmat)
    {
        tmat->~MATRIX();
        MEMORY ps_free (tmatptr);
    }

    gsl_integration_workspace_free ((gsl_integration_workspace*)gsl_workspace_tps_x);
    gsl_integration_workspace_free ((gsl_integration_workspace*)gsl_workspace_tps_y);
}

