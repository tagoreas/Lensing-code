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



#include "pixsrc_common_adaptive.hpp"
#include "pixsrc_shapelets_operations.hpp"
#include "pixsrc_common.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_regularization.hpp"
#include "pixsrc_shapelets.hpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_d2matrix.cpp"
#include <unistd.h> // for usleep
#include <cstring>
#include <cmath>
#include <algorithm>

#include "pixsrc_clipper.hpp"
using namespace ClipperLib;

void pixsrc_common_adaptive::createcbody( inputdata *data_, commoninputdata *cdata_, lensvar *vars_ )
{
    if (vars_->fatalerror)
    {
        return;
    }

    // if no regularization is being done
    if (!data_->reg)
    {
        // this is so the program doesn't hang forever elsewhere
        pthread_create(&(vars_->pthreads[1]), cdata_->attrdetached,
                       STATISTIC computedetc, vars_->dav);
        return;
    }

    // waiting for triangulation and unused pixel removal
    PS_SIT dim = 1;
    PS_SIT waitlist[1] = {7};
    OPERA pthreadswait(vars_,dim,waitlist);

    // waiting for shapelet parameters to be optimized
    if (data_->use_shapelets)
    {
        // sleep for 10 ms
        while( !vars_->shapeletopt )
            usleep( CONSTANT longwaitms );
    }

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "computing regularization",
                             cdata_->print2screenmutex);

    PS_SIT initsize = !data_->numgpu2use ? vars_->lonc : 0;
    MEMORY ps_malloc( &(vars_->c1ptr), 1 );
    vars_->c1 = new (vars_->c1ptr)
        MATRIX (cdata_, data_, vars_->lonc, vars_->lonc, initsize, 1, data_);

    // waiting for source to be optimized
    if( data_->usersetsrc )
    {
        while( !vars_->srcopt )
        {
            // sleep for 10 ms
            usleep( CONSTANT longwaitms );
        }

        // set all pixel fluxes below a minimum value to that value
        double min = std::sqrt(vars_->variance) / 10;
        for( PS_SIT x=0; x<vars_->mps->get_size(); ++x )
            if( vars_->mps->get(x) < min )
                vars_->mps->set( x, min );
        // normalize
        vars_->mps->mult( 1.0 / vars_->mps->max() );
        // remove image plane analytic source
        if (vars_->srcexactvec)
        {
            vars_->srcexactvec->~VECTOR ();
            MEMORY ps_free (vars_->srcexactvecptr);
            vars_->srcexactvecptr = NULL;
        }
    }

    if( data_->usersetsrc && data_->regaccuracy==0 )
    {
        REGULAR getsersicregmin( data_, cdata_, vars_ );
    }
    else if( data_->usersetsrc && data_->regaccuracy==1 )
    {
        REGULAR getsersicregmax( data_, cdata_, vars_ );
    }
    else if (data_->use_shapelets)
    {
        REGULAR regshapelets (data_, cdata_, vars_);
    }
    else if( data_->regorder==0 )
    {
        REGULAR getder0( data_, cdata_, vars_ );
    }
    else
    {
        if( data_->regaccuracy==0 )
        {
            REGULAR getder12( data_, cdata_, vars_ );
        }
        else if( data_->regaccuracy==1 )
        {
            if( data_->regorder==2 )
                REGULAR greengaussgetder2( data_, cdata_, vars_ );
            else if( data_->regorder==1 )
                REGULAR greengaussgetder1( data_, cdata_, vars_ );
        }
    }

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "done computing regularization",
                             cdata_->print2screenmutex);

    pthread_create(&(vars_->pthreads[1]), cdata_->attrdetached,
                   STATISTIC computedetc, vars_->dav);

    if (data_->debug)
    {
        if(vars_->h1)
            PRINTER print( vars_->tracker, cdata_->basename,
                           data_->name, true, "h1.MATRIX",
                           vars_->lonc, vars_->lonc, vars_->h1, 0);

        if(vars_->h1x)
            PRINTER print( vars_->tracker, cdata_->basename,
                           data_->name, true, "h1x.MATRIX",
                           vars_->lonc, vars_->lonc, vars_->h1x, 0);

        if(vars_->h1y)
            PRINTER print( vars_->tracker, cdata_->basename,
                           data_->name, true, "h1y.MATRIX",
                           vars_->lonc, vars_->lonc, vars_->h1y, 0);

        PRINTER print( vars_->tracker, cdata_->basename,
                       data_->name, true, "c1.MATRIX",
                       vars_->lonc, vars_->lonc, vars_->c1, 0);
    }
}

void* pixsrc_common_adaptive::createc(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[0] = 1;

    COMMONADAPTIVE createcbody( data_, cdata_, vars_ );

    vars_->pthreadstracker[0] = 2;

    return NULL;
}

void pixsrc_common_adaptive::printsourcebody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    VECTOR *noise_dummy=0, *noise=0;

    if (data_->noisemap)
        if (data_->lambdaguess &&
            (!data_->usersetsrc || data_->reg) &&
            data_->printvec)
        {
            if( data_->verbose!=1 )
                PRINTER print2screen (data_->print2screenname,
                                      "computing (optional) source plane noise; could take a while",
                                      cdata_->print2screenmutex);

            // compute noise
            MEMORY ps_malloc( &noise_dummy, 1 );
            noise = new (noise_dummy) VECTOR( cdata_, data_, vars_->lonc );
            vars_->a1->noise_invA( noise, cdata_->numthreads, NULL );

            if (0)
            {

                PS_SIT no_chuck=1;
                double **a_chuck;
                if (!no_chuck)
                {
                    // load in Chuck's matrices for testing method
                    char *fname;
                    char **vecin;
                    PS_SIT vecinsize;
                    char ext[]="temp1_4";
                    const char *listcc[2];
                    listcc[0] = CONSTANT dir_in;
                    listcc[1] = ext;
                    OPERA concatenate( listcc, 2, &fname );
                    OPERA readfile( fname, &vecin, &vecinsize, NULL);
                    MEMORY ps_free( fname );
                    MEMORY ps_malloc( &a_chuck, vars_->a1->nrow, vars_->a1->ncol);
                    MATRIX achuck(cdata_, data_, vars_->lonc,
                                  vars_->lonc, vars_->lonc, 1,  data_);
                    for (PS_SIT v=0; v<vecinsize; ++v)
                    {
                        char *tok = strtok(vecin[v]," ");
                        const char *zzz = tok;
                        PS_SIT fir = OPERA convert_string <PS_SIT> (zzz);
                        tok = strtok(NULL," ");
                        zzz = tok;
                        PS_SIT sec = OPERA convert_string <PS_SIT> (zzz);
                        tok = strtok(NULL," ");
                        zzz = tok;
                        double thi = OPERA convert_string <double> (zzz);
                        a_chuck[fir][sec]=thi;
                        if(thi)
                            achuck.set (fir,sec,thi);
                    }
                    PRINTER print ( vars_->tracker, cdata_->basename, data_->name, true,
                                    "a1_chuck.MATRIX", vars_->lonc, vars_->lonc, &achuck, 0);
                }





                // compute source-plane resolution

                if( data_->verbose!=1 )
                {
                    PRINTER print2screen( data_->print2screenname,
                                          "doing dense matrix stuff to get source plane resolution",
                                          cdata_->print2screenmutex);
                }

                // perform SVD --- A = U D V^T

                char jobu  = 'O'; // overwrite matrix A
                char jobvt = 'N'; // do not compute
                PS_SIT m = vars_->a1->nrow;
                PS_SIT n = vars_->a1->ncol;
                PS_SIT lda  = m;     // leading dimension of matrix A
                PS_SIT ldu  = 1;     // not computing; not needed
                PS_SIT ldvt = 1;     // not computing; not needed
                PS_SIT info=0;         // status upon return (if XERBLA deosn't crash)
                PS_SIT lwork;        // "workspace"-related variable
                double wkopt;     // workspace neeed
                double *work;     // actual workspace vector
                double *u  =  0;  // left matrix in factorization (null, not needed)
                double *vt =  0;  // right matrix in factorization (null; not needed)
                double *a  =  0;  // matrix to be factored
                double *s  =  0;  // singular values

                // store vars_->a1 into dense format
                MEMORY ps_malloc( &a, m*n );
                MEMORY ps_malloc( &s, m   );
                for (PS_SIT i=0; i<vars_->a1->ncol; ++i)
                    for (PS_SIT j=0; j<vars_->a1->nrow; ++j)
                        if (no_chuck)
                            a[i*m+j] = vars_->a1->get(j,i);
                        else
                            a[i*m+j] = a_chuck[j][i];

                // calculate workspace needed
                lwork = -1;
                EXTERNAL ps_dgesvd_( &jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt,
                                     &ldvt, &wkopt, &lwork, &info );

                // allocate workspace
                lwork = (PS_SIT)wkopt;
                MEMORY ps_malloc( &work, lwork );

                // actually perform SVD
                EXTERNAL ps_dgesvd_( &jobu, &jobvt, &m, &n, a, &lda, s, u, &ldu, vt,
                                     &ldvt, work,   &lwork, &info );

                // check for success
                if (info)
                {
                    PRINTER printerror( data_->print2screenname,
                                        "SVD failed",
                                        cdata_->print2screenmutex );
                }

                // calculate inverse of left matrix --- u_inv = u^-1

                int    *ipiv;  // pivot indices (for efficient LU factorization)
                double *u_inv; // inverse of left matrix in SVD

                // allocate memory and store copy of left matrix u for inversion
                MEMORY ps_malloc (&ipiv,  m  );
                MEMORY ps_malloc (&u_inv, m*n);
                std::copy (a, a+m*n, u_inv);

                // perform LU factorization
                EXTERNAL ps_dgetrf_ (&m, &n, u_inv, &lda, ipiv, &info);

                // check for success
                if (info)
                {
                    PRINTER printerror( data_->print2screenname,
                                        "dense LU factorization failed",
                                        cdata_->print2screenmutex );
                }

                // calculate workspace needed
                lwork = -1;
                EXTERNAL ps_dgetri_ (&n, u_inv, &lda, ipiv, &wkopt, &lwork, &info);

                // allocate workspace
                lwork = (PS_SIT)wkopt;
                MEMORY ps_free( work );
                work = 0;
                MEMORY ps_malloc( &work, lwork );

                // perform matrix inversion
                EXTERNAL ps_dgetri_ (&n, u_inv, &lda, ipiv, work, &lwork, &info);

                // check for success
                if (info)
                {
                    PRINTER printerror( data_->print2screenname,
                                        "dense matrix inversion failed",
                                        cdata_->print2screenmutex );
                }

                // perform matrix matrix multiplications --- R = u D^-1/2 u_inv

                // multiplying diagonal matrix by inverted matrix
                double fac;
                for (PS_SIT i=0; i<m; ++i)
                    if (s[i])
                    {
                        fac = 1.0 / std::sqrt(s[i]);
                        for (PS_SIT j=0; j<m; ++j)
                            u_inv[j*m+i] = u_inv[j*m+i]*fac;
                    }

                double *r;
                MEMORY ps_malloc (&r, m*n);
                char trans = 'N';
                double alpha = 1.0;
                double beta  = 0.0;

                // perform matrix matrix multiplication
                EXTERNAL ps_dgemm_ ( &trans, &trans, &m, &m, &m,
                                     &alpha, a, &m, u_inv, &m, &beta, r, &m);

                if( data_->verbose!=1 )
                {
                    PRINTER print2screen( data_->print2screenname,
                                          "done doing dense matrix stuff to get source plane resolution",
                                          cdata_->print2screenmutex);
                }

                // print out result
                double **vec;
                MEMORY ps_malloc( &vec, m*n, 3 );
                for (PS_SIT i=0; i<n; ++i)
                    for (PS_SIT j=0; j<m; ++j)
                    {
                        vec[i*m+j][0] = i;
                        vec[i*m+j][1] = j;
                        vec[i*m+j][2] = r[i*m+j];
                    }
                PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                         "sourcebeam.VECTOR", vec, m*n, 3, 0, data_->precision, NULL);

                // also print out a1^-1 as well..
                // a1^-1 = r*r, so we do the multiplication
                // we will re-use "u_inv" variable here
                EXTERNAL ps_dgemm_ ( &trans, &trans, &m, &m, &m,
                                     &alpha, r, &m, r, &m, &beta, u_inv, &m);
                for (PS_SIT i=0; i<n; ++i)
                    for (PS_SIT j=0; j<m; ++j)
                        vec[i*m+j][2] = u_inv[i*m+j];
                PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                         "inv_a1.VECTOR", vec, m*n, 3, 0, data_->precision, NULL);

                // also print out a1^-1*a as well as a sanity check
                // we will re-use "r" variable here
                // first, restore a1 into variable "a"
                for (PS_SIT i=0; i<vars_->a1->ncol; ++i)
                    for (PS_SIT j=0; j<vars_->a1->nrow; ++j)
                        if (no_chuck)
                            a[i*m+j] = vars_->a1->get(j,i);
                        else
                            a[i*m+j] = a_chuck[j][i];
                EXTERNAL ps_dgemm_ ( &trans, &trans, &m, &m, &m,
                                     &alpha, u_inv, &m, a, &m, &beta, r, &m);
                for (PS_SIT i=0; i<n; ++i)
                    for (PS_SIT j=0; j<m; ++j)
                        vec[i*m+j][2] = r[i*m+j];
                PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                         "inv_a1_a1.VECTOR", vec, m*n, 3, 0, data_->precision, NULL);


                // memory cleanup
                MEMORY ps_free( vec, m*n );
                MEMORY ps_free( work     );
                MEMORY ps_free(  a       );
                MEMORY ps_free(  s       );
                MEMORY ps_free(  r       );
                MEMORY ps_free( ipiv     );
                MEMORY ps_free( u_inv    );
                if (!no_chuck)
                    MEMORY ps_free (a_chuck, vars_->a1->nrow);
            }
        }

    if( data_->printvec==1 || data_->printvec==3 )
    {
        vars_->mps->update_cpu();

        double **mps;
        MEMORY ps_malloc( &mps, vars_->lonc, 3 );

        for(PS_SIT x=0; x<vars_->lonc; x++)
        {
            if (data_->use_shapelets)
            {
                if (data_->use_shapelets==1)
                {
                    mps[x][0] = x/vars_->num_shapelets2;
                    mps[x][1] = x%vars_->num_shapelets2;
                }
                else if (data_->use_shapelets==2)
                {
                    PRINTER printerror(data_->print2screenname,
                                       "polar shapelets disabled",
                                       cdata_->print2screenmutex);
                    /*
                      PS_SIT mps_shapelet_ind[2] = {0,0};
                      mps[x][0] = mps_shapelet_ind[0];
                      mps[x][1] = mps_shapelet_ind[1];

                      // go to next radial mode if ready
                      if (mps_shapelet_ind[0]==-mps_shapelet_ind[1])
                      {
                      ++mps_shapelet_ind[0];
                      mps_shapelet_ind[1] = (PS_SIT)mps_shapelet_ind[0]%2;
                      }
                      // otherwise increment angular mode if zero
                      else if (mps_shapelet_ind[1]==0)
                      {
                      mps_shapelet_ind[1]+=2;
                      }
                      // otherwise increment angular mode if negative
                      else if (mps_shapelet_ind[1]<0)
                      {
                      mps_shapelet_ind[1] = -mps_shapelet_ind[1]+2;
                      }
                      // otherwise get imaginary angular mode if positive
                      else
                      {
                      mps_shapelet_ind[1] *= -1;
                      }
                    */
                }
            }
            else
            {
                mps[x][0]=vars_->triout->pointlist[x*2  ];
                mps[x][1]=vars_->triout->pointlist[x*2+1];
            }
            mps[x][2]=vars_->mps->get(x);
        }
        PRINTER print <double> ( vars_->tracker, cdata_->basename,
                                 data_->name, true,"mps.VECTOR",mps,vars_->lonc,3, 0, data_->precision, NULL);

        if (data_->noisemap)
            if( data_->lambdaguess && ( !data_->usersetsrc || data_->reg ) )
            {
                for(PS_SIT x=0; x<vars_->lonc; x++)
                    mps[x][2]=noise->get(x);
                PRINTER print <double> ( vars_->tracker, cdata_->basename,
                                         data_->name, true,"noise.VECTOR",mps,vars_->lonc,3, 0, data_->precision, NULL);
                for(PS_SIT x=0; x<vars_->lonc; x++)
                    mps[x][2] = vars_->mps->get(x) / noise->get(x);
                PRINTER print <double> ( vars_->tracker, cdata_->basename,
                                         data_->name, true,"s2n.VECTOR",mps,vars_->lonc,3, 0, data_->precision, NULL);
            }
        MEMORY ps_free( mps, vars_->lonc );
    }

    if( data_->printvec==2 || data_->printvec==3 )
    {
        actualprintsource(data_, cdata_, vars_, noise);
    }

    if (data_->noisemap)
        if( data_->lambdaguess && ( !data_->usersetsrc || data_->reg ) && data_->printvec )
        {
            noise->~VECTOR();
            MEMORY ps_free( noise_dummy );
        }
}

void* pixsrc_common_adaptive::printsource(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[9] = 1;

    COMMONADAPTIVE printsourcebody( data_, cdata_, vars_ );

    vars_->pthreadstracker[9] = 2;

    return NULL;
}

void pixsrc_common_adaptive::getmagnificationbody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->maguncertainty ) {

        //double mag_errs[7] = {-1,-1,-1,-1,-1,-1,-1};

        if (data_->fullsrccov)
            getmagfromtriangulation_with_err(data_, cdata_, vars_);
        else
            getmagfromtriangulation_with_err_diagonals(data_, cdata_, vars_);

        //PRINTER writeoutstream <double> (mag_errs,7,data_->mag_errs->stream,data_->mag_errs->lock, data_->precision, NULL);
    }

    if(data_->magparams || data_->penaltyquery[1])
    {
        if( data_->magparams || data_->penaltyquery[1] ) {

            getmagfromtriangulation(data_, cdata_, vars_);

            if(data_->verbose!=1)
                PRINTER print2screen(data_->print2screenname,
                                     "magnification #" + OPERA tostring(vars_->tracker) +
                                     " = " + OPERA tostring(vars_->magger),
                                     cdata_->print2screenmutex);

            if(data_->magparams)
            {
                double writeout[2];
                writeout[0] = vars_->lambda1;
                writeout[1] = vars_->magger;
                PRINTER writeoutstream <double> (writeout,2,data_->mags->stream,data_->mags->lock, data_->precision, NULL);
            }

            COMMON computemagpenalty(data_, cdata_, vars_);
        }
    }
}

void* pixsrc_common_adaptive::getmagnification(void* args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[10] = 1;

    COMMONADAPTIVE getmagnificationbody( data_, cdata_, vars_ );

    vars_->pthreadstracker[10] = 2;

    return NULL;
}

void pixsrc_common_adaptive::getmagfromtriangulation(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    double srcflux=0,imgflux=0;

    if (data_->use_shapelets)
    {
        // get square roots of binomial coefficients
        double *binom;
        PS_SIT max_num = std::max (vars_->num_shapelets1,vars_->num_shapelets2);
        MEMORY ps_malloc (&binom, (max_num+1)/2);
        SHAPELETSOPERA get_binomials (data_, cdata_, vars_, binom, max_num);

        // integrate shapelets to get total flux (pixel^2 units)
        for (PS_SIT s1=0; s1<vars_->num_shapelets1; s1+=2)
        {
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; s2+=2)
            {
                PS_SIT sindex = s1*vars_->num_shapelets2 + s2;
                srcflux += std::pow (2.0,0.5*(2-s1-s2))
                    * binom[s1/2] * binom[s2/2]
                    * vars_->mps->get (sindex);
            }
        }
        srcflux *= CONSTANT sqrtpi * vars_->shapelet_scale;

        MEMORY ps_free (binom);
    }
    else
    {
        double coords[6];
        for(PS_SIT t=0; t<vars_->triout->numberoftriangles; t++)
        {
            for(PS_SIT j=0; j<3; j++)
            {
                coords[j*2]   = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2  ];
                coords[j*2+1] = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2+1];
            }

            srcflux += GEOM areatriangle(coords) *
                (
                    vars_->mps->get(vars_->triout->trianglelist[t*3  ])+
                    vars_->mps->get(vars_->triout->trianglelist[t*3+1])+
                    vars_->mps->get(vars_->triout->trianglelist[t*3+2])
                    ) / 3.0;
        }
    }

    //for(PS_SIT r=0; r<vars_->lonr; r++)
    //imgflux += vars_->data->get(r);

    // get the flux that's already been got

    VECTOR *dummy;
    MEMORY ps_malloc( &dummy, 1 );
    VECTOR *lensedmpsnobo = new (dummy) VECTOR(cdata_, data_, vars_->lonr);
    COMMON lensgalaxy( data_, cdata_, vars_, vars_->lensingoperatornobo,
                       vars_->mps, lensedmpsnobo);

    for(PS_SIT r=0; r<vars_->lonr; r++)
        imgflux += lensedmpsnobo->get(r);

    lensedmpsnobo->~VECTOR();
    MEMORY ps_free( dummy );

    if (data_->fullmag)
    {
        // get the flux that hasn't been got for whatever reason
        // NOTE: CURRENTLY NOT TAKING INTERPOLATION ERRORS INTO ACCOUNT
        //          FOR DATA PIXELS NOT CONTAINED WITHIN LENSING OPERATOR.
        if (data_->use_shapelets)
        {
            // create work-space
            double **workspace;
            MEMORY ps_malloc (&workspace, 2, vars_->numberofshapelets);

            // if we haven't' gotten the flux yet, get it
            for (PS_SIT r=0; r<data_->ndp; ++r)
                if (-1==vars_->r4r[r])
                    imgflux += SHAPELETSOPERA get_flux_one_pixel (data_, cdata_, vars_, r, workspace);

            MEMORY ps_free (workspace, 2);
        }
        else
        {
            double pos[2][3];
            double weights3[3];
            char *gotitalready;
            MEMORY ps_malloc( &(gotitalready), data_->ndp );
            std::fill(gotitalready,gotitalready+data_->ndp,0);

            for(PS_SIT tri = 0; tri < vars_->triout->numberoftriangles; tri++)
            {
                // NOTE: rewrite: use the newer, faster ray-trace function in pixsrc::geometry
                for(PS_SIT f=0; f<3; f++)
                {
                    pos[0][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2];
                    pos[1][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2+1];
                }
                for(PS_SIT r=0; r<data_->ndp; r++)
                {
                    if( vars_->r4r[r]==-1 && !gotitalready[r] &&
                        GEOM isintri(pos,vars_->newloc[r*2],vars_->newloc[r*2+1]))
                    {
                        OPERA planarinterpolation3pts(vars_->newloc[r*2],vars_->newloc[r*2+1],
                                                      pos,weights3);
                        for(PS_SIT f=0; f<3; f++)
                            imgflux += weights3[f]*vars_->mps->get(vars_->triout->trianglelist[tri*3+f]);
                        gotitalready[r] = 1;
                    }
                }
            }

            MEMORY ps_free( gotitalready );
        }
    }

    vars_->magger = imgflux/srcflux;
}

void pixsrc_common_adaptive::getmagfromtriangulation_with_err(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    double srcflux=0;

    int numsamples = data_->magsamples;
    PS_SIT nummagmasks = *(std::max_element(data_->magmasks, data_->magmasks + data_->ndp)) + 1;
    if (0 == nummagmasks)
      return;

    double *allmags = new double[numsamples * (nummagmasks+1)];
    double imgflux[nummagmasks+1];

    // compute noise
    //VECTOR *noise_dummy=0,  *noise=0;
    double *noise = new double[vars_->lonc * vars_->lonc];
    std::fill(noise, noise+vars_->lonc*vars_->lonc, 0);
    VECTOR *newmps_dummy=0, *newmps=0; 
    //MEMORY ps_malloc( &noise_dummy, 1 );
    MEMORY ps_malloc( &newmps_dummy, 1 );
    //noise  = new (noise_dummy)  VECTOR( cdata_, data_, vars_->lonc );
    newmps = new (newmps_dummy) VECTOR( cdata_, data_, vars_->lonc );
    vars_->a1->noise_invA_fullmat( noise, cdata_->numthreads, NULL );

    for (int magno = 0; magno < numsamples; ++magno) {

      std::copy( vars_->mps->vec, vars_->mps->vec + vars_->lonc, newmps->vec );
      double *randvals = new double[vars_->lonc];
      for (int jj=0; jj<vars_->lonc; ++jj) {
	  randvals[jj] = OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      }
      for (int jj=0; jj<vars_->lonc; ++jj) {
	double newval = 0;
	for (int jjj=0; jjj<=jj; ++jjj)
	  newval += noise[jj*vars_->lonc +jjj] * randvals[jjj];
	newmps->set(jj, vars_->mps->get(jj) + newval);
      }
      delete [] randvals;


    if (data_->use_shapelets)
    {
        // get square roots of binomial coefficients
        double *binom;
        PS_SIT max_num = std::max (vars_->num_shapelets1,vars_->num_shapelets2);
        MEMORY ps_malloc (&binom, (max_num+1)/2);
        SHAPELETSOPERA get_binomials (data_, cdata_, vars_, binom, max_num);

        // integrate shapelets to get total flux (pixel^2 units)
        for (PS_SIT s1=0; s1<vars_->num_shapelets1; s1+=2)
        {
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; s2+=2)
            {
                PS_SIT sindex = s1*vars_->num_shapelets2 + s2;
                srcflux += std::pow (2.0,0.5*(2-s1-s2))
                    * binom[s1/2] * binom[s2/2]
                    * newmps->get (sindex);
            }
        }
        srcflux *= CONSTANT sqrtpi * vars_->shapelet_scale;

        MEMORY ps_free (binom);
    }
    else
    {
        double coords[6];
        for(PS_SIT t=0; t<vars_->triout->numberoftriangles; t++)
        {
            for(PS_SIT j=0; j<3; j++)
            {
                coords[j*2]   = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2  ];
                coords[j*2+1] = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2+1];
            }

            srcflux += GEOM areatriangle(coords) *
                (
                    newmps->get(vars_->triout->trianglelist[t*3  ])+
                    newmps->get(vars_->triout->trianglelist[t*3+1])+
                    newmps->get(vars_->triout->trianglelist[t*3+2])
                    ) / 3.0;
        }
    }

    // get the flux that's already been got
    VECTOR *dummy;
    MEMORY ps_malloc( &dummy, 1 );
    VECTOR *lensedmpsnobo = new (dummy) VECTOR(cdata_, data_, vars_->lonr);
    COMMON lensgalaxy( data_, cdata_, vars_, vars_->lensingoperatornobo,
                       newmps, lensedmpsnobo);

    for(PS_SIT r=0; r<vars_->lonr; r++) {
      PS_SIT magmaskind = data_->magmasks[vars_->r4rback[r]];
      if (magmaskind != -1) {
        imgflux[magmaskind]  += lensedmpsnobo->get(r);
	imgflux[nummagmasks] += lensedmpsnobo->get(r);
      }
    }

    lensedmpsnobo->~VECTOR();
    MEMORY ps_free( dummy );

    if (data_->fullmag)
    {
        // get the flux that hasn't been got for whatever reason
        // NOTE: CURRENTLY NOT TAKING INTERPOLATION ERRORS INTO ACCOUNT
        //          FOR DATA PIXELS NOT CONTAINED WITHIN LENSING OPERATOR.
        if (data_->use_shapelets)
        {
            // create work-space
            double **workspace;
            MEMORY ps_malloc (&workspace, 2, vars_->numberofshapelets);

            // if we haven't' gotten the flux yet, get it
            for (PS_SIT r=0; r<data_->ndp; ++r)
	      if (-1==vars_->r4r[r]) {
		PS_SIT magmaskind = data_->magmasks[r];
		if (magmaskind != -1) {
		  double thenewflux = SHAPELETSOPERA get_flux_one_pixel (data_, cdata_, vars_, r, workspace, newmps);
		  imgflux[magmaskind]  += thenewflux;
		  imgflux[nummagmasks] += thenewflux;
		}
	      }

            MEMORY ps_free (workspace, 2);
        }
        else
        {
            double pos[2][3];
            double weights3[3];
            char *gotitalready;
            MEMORY ps_malloc( &(gotitalready), data_->ndp );
            std::fill(gotitalready,gotitalready+data_->ndp,0);

            for(PS_SIT tri = 0; tri < vars_->triout->numberoftriangles; tri++)
            {
                // NOTE: rewrite: use the newer, faster ray-trace function in pixsrc::geometry
                for(PS_SIT f=0; f<3; f++)
                {
                    pos[0][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2];
                    pos[1][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2+1];
                }
                for(PS_SIT r=0; r<data_->ndp; r++)
                {
                    if( vars_->r4r[r]==-1 && !gotitalready[r] &&
                        GEOM isintri(pos,vars_->newloc[r*2],vars_->newloc[r*2+1]))
                    {
                        OPERA planarinterpolation3pts(vars_->newloc[r*2],vars_->newloc[r*2+1],
                                                      pos,weights3);
                        for(PS_SIT f=0; f<3; f++) {
			  PS_SIT magmaskind = data_->magmasks[r];
			  if (magmaskind != -1) {
                            imgflux[magmaskind]  += weights3[f]*newmps->get(vars_->triout->trianglelist[tri*3+f]);
                            imgflux[nummagmasks] += weights3[f]*newmps->get(vars_->triout->trianglelist[tri*3+f]);
			  }
			}
                        gotitalready[r] = 1;
                    }
                }
            }

            MEMORY ps_free( gotitalready );
        }
    }

    for (int jj=0; jj<nummagmasks+1; ++jj) {
      allmags[numsamples*jj  + magno] = imgflux[jj] / srcflux;
    }


    if (0) {

        int ind1 = std::floor(0.025 * magno);
        int ind2 = std::floor(0.050 * magno);
        int ind3 = std::floor(0.160 * magno);
        int ind4 = std::floor(0.500 * magno);
        int ind5 = std::floor(0.840 * magno);
        int ind6 = std::floor(0.950 * magno);
        int ind7 = std::floor(0.975 * magno);

        for (int jj=0; jj<nummagmasks+1; ++jj) {
          
          std::sort(allmags + jj*numsamples, allmags + jj*numsamples + magno + 1);

          PRINTER print2screen(data_->print2screenname,
    			   "mags " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind1)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind2)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind3)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind4)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind5)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind6)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind7)),
    			   cdata_->print2screenmutex);
        }
    }


    
    }

    int ind0  = std::floor(0.000 * numsamples);
    int ind1  = std::floor(0.025 * numsamples);
    int ind2  = std::floor(0.050 * numsamples);
    int ind3  = std::floor(0.160 * numsamples);
    int ind4  = std::floor(0.300 * numsamples);
    int ind5  = std::floor(0.400 * numsamples);
    int ind6  = std::floor(0.500 * numsamples);
    int ind7  = std::floor(0.600 * numsamples);
    int ind8  = std::floor(0.700 * numsamples);
    int ind9  = std::floor(0.840 * numsamples);
    int ind10 = std::floor(0.950 * numsamples);
    int ind11 = std::floor(0.975 * numsamples);
    int ind12 = std::floor(0.99999999 * numsamples);

    for (int jj=0; jj<nummagmasks+1; ++jj) {
      
      std::sort(allmags + jj*numsamples, allmags + (jj+1)*numsamples);

      double writeout[13];
      writeout[0] = *(allmags + jj*numsamples + ind0);
      writeout[1] = *(allmags + jj*numsamples + ind1);
      writeout[2] = *(allmags + jj*numsamples + ind2);
      writeout[3] = *(allmags + jj*numsamples + ind3);
      writeout[4] = *(allmags + jj*numsamples + ind4);
      writeout[5] = *(allmags + jj*numsamples + ind5);
      writeout[6] = *(allmags + jj*numsamples + ind6);
      writeout[7] = *(allmags + jj*numsamples + ind7);
      writeout[8] = *(allmags + jj*numsamples + ind8);
      writeout[9] = *(allmags + jj*numsamples + ind9);
      writeout[10] = *(allmags + jj*numsamples + ind10);
      writeout[11] = *(allmags + jj*numsamples + ind11);
      writeout[12] = *(allmags + jj*numsamples + ind12);

      PRINTER writeoutstream <double> (writeout,13,data_->mag_errs->stream,data_->mag_errs->lock, data_->precision, NULL);
    }
    PRINTER writeoutstream ("",data_->mag_errs->stream,data_->mag_errs->lock, data_->precision, NULL);

    //noise->~VECTOR();
    //MEMORY ps_free( noise_dummy );
    newmps->~VECTOR();
    MEMORY ps_free( newmps_dummy );

    delete [] allmags;
    delete [] noise;
}

void pixsrc_common_adaptive::getmagfromtriangulation_with_err_diagonals(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    double srcflux=0;

    int numsamples = data_->magsamples;
    PS_SIT nummagmasks = *(std::max_element(data_->magmasks, data_->magmasks + data_->ndp)) + 1;
    if (0 == nummagmasks)
      return;

    double *allmags = new double[numsamples * (nummagmasks+1)];
    double imgflux[nummagmasks+1];

    // compute noise
    VECTOR *noise_dummy=0,  *noise=0;
    VECTOR *newmps_dummy=0, *newmps=0; 
    MEMORY ps_malloc( &noise_dummy, 1 );
    MEMORY ps_malloc( &newmps_dummy, 1 );
    noise  = new (noise_dummy)  VECTOR( cdata_, data_, vars_->lonc );
    newmps = new (newmps_dummy) VECTOR( cdata_, data_, vars_->lonc );
    vars_->a1->noise_invA( noise, cdata_->numthreads, NULL );

    for (int magno = 0; magno < numsamples; ++magno) {

      std::copy( vars_->mps->vec, vars_->mps->vec + vars_->lonc, newmps->vec );
      for (int jj=0; jj<vars_->lonc; ++jj) {
	double newval = -1;
	//while (newval <0)
	  newval = vars_->mps->get(jj) + noise->get(jj) * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
	  newmps->set(jj, newval);
      }


    if (data_->use_shapelets)
    {
        // get square roots of binomial coefficients
        double *binom;
        PS_SIT max_num = std::max (vars_->num_shapelets1,vars_->num_shapelets2);
        MEMORY ps_malloc (&binom, (max_num+1)/2);
        SHAPELETSOPERA get_binomials (data_, cdata_, vars_, binom, max_num);

        // integrate shapelets to get total flux (pixel^2 units)
        for (PS_SIT s1=0; s1<vars_->num_shapelets1; s1+=2)
        {
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; s2+=2)
            {
                PS_SIT sindex = s1*vars_->num_shapelets2 + s2;
                srcflux += std::pow (2.0,0.5*(2-s1-s2))
                    * binom[s1/2] * binom[s2/2]
                    * newmps->get (sindex);
            }
        }
        srcflux *= CONSTANT sqrtpi * vars_->shapelet_scale;

        MEMORY ps_free (binom);
    }
    else
    {
        double coords[6];
        for(PS_SIT t=0; t<vars_->triout->numberoftriangles; t++)
        {
            for(PS_SIT j=0; j<3; j++)
            {
                coords[j*2]   = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2  ];
                coords[j*2+1] = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2+1];
            }

            srcflux += GEOM areatriangle(coords) *
                (
                    newmps->get(vars_->triout->trianglelist[t*3  ])+
                    newmps->get(vars_->triout->trianglelist[t*3+1])+
                    newmps->get(vars_->triout->trianglelist[t*3+2])
                    ) / 3.0;
        }
    }

    // get the flux that's already been got
    VECTOR *dummy;
    MEMORY ps_malloc( &dummy, 1 );
    VECTOR *lensedmpsnobo = new (dummy) VECTOR(cdata_, data_, vars_->lonr);
    COMMON lensgalaxy( data_, cdata_, vars_, vars_->lensingoperatornobo,
                       newmps, lensedmpsnobo);

    for(PS_SIT r=0; r<vars_->lonr; r++) {
      PS_SIT magmaskind = data_->magmasks[vars_->r4rback[r]];
      if (magmaskind != -1) {
        imgflux[magmaskind]  += lensedmpsnobo->get(r);
	imgflux[nummagmasks] += lensedmpsnobo->get(r);
      }
    }

    lensedmpsnobo->~VECTOR();
    MEMORY ps_free( dummy );

    if (data_->fullmag)
    {
        // get the flux that hasn't been got for whatever reason
        // NOTE: CURRENTLY NOT TAKING INTERPOLATION ERRORS INTO ACCOUNT
        //          FOR DATA PIXELS NOT CONTAINED WITHIN LENSING OPERATOR.
        if (data_->use_shapelets)
        {
            // create work-space
            double **workspace;
            MEMORY ps_malloc (&workspace, 2, vars_->numberofshapelets);

            // if we haven't' gotten the flux yet, get it
            for (PS_SIT r=0; r<data_->ndp; ++r)
	      if (-1==vars_->r4r[r]) {
		PS_SIT magmaskind = data_->magmasks[r];
		if (magmaskind != -1) {
		  double thenewflux = SHAPELETSOPERA get_flux_one_pixel (data_, cdata_, vars_, r, workspace, newmps);
		  imgflux[magmaskind]  += thenewflux;
		  imgflux[nummagmasks] += thenewflux;
		}
	      }

            MEMORY ps_free (workspace, 2);
        }
        else
        {
            double pos[2][3];
            double weights3[3];
            char *gotitalready;
            MEMORY ps_malloc( &(gotitalready), data_->ndp );
            std::fill(gotitalready,gotitalready+data_->ndp,0);

            for(PS_SIT tri = 0; tri < vars_->triout->numberoftriangles; tri++)
            {
                // NOTE: rewrite: use the newer, faster ray-trace function in pixsrc::geometry
                for(PS_SIT f=0; f<3; f++)
                {
                    pos[0][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2];
                    pos[1][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2+1];
                }
                for(PS_SIT r=0; r<data_->ndp; r++)
                {
                    if( vars_->r4r[r]==-1 && !gotitalready[r] &&
                        GEOM isintri(pos,vars_->newloc[r*2],vars_->newloc[r*2+1]))
                    {
                        OPERA planarinterpolation3pts(vars_->newloc[r*2],vars_->newloc[r*2+1],
                                                      pos,weights3);
                        for(PS_SIT f=0; f<3; f++) {
			  PS_SIT magmaskind = data_->magmasks[r];
			  if (magmaskind != -1) {
                            imgflux[magmaskind]  += weights3[f]*newmps->get(vars_->triout->trianglelist[tri*3+f]);
                            imgflux[nummagmasks] += weights3[f]*newmps->get(vars_->triout->trianglelist[tri*3+f]);
			  }
			}
                        gotitalready[r] = 1;
                    }
                }
            }

            MEMORY ps_free( gotitalready );
        }
    }

    for (int jj=0; jj<nummagmasks+1; ++jj) {
      allmags[numsamples*jj  + magno] = imgflux[jj] / srcflux;
    }

    if (0) {

        int ind1 = std::floor(0.025 * magno);
        int ind2 = std::floor(0.050 * magno);
        int ind3 = std::floor(0.160 * magno);
        int ind4 = std::floor(0.500 * magno);
        int ind5 = std::floor(0.840 * magno);
        int ind6 = std::floor(0.950 * magno);
        int ind7 = std::floor(0.975 * magno);

        for (int jj=0; jj<nummagmasks+1; ++jj) {
          
          std::sort(allmags + jj*numsamples, allmags + jj*numsamples + magno + 1);

          PRINTER print2screen(data_->print2screenname,
    			   "mags " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind1)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind2)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind3)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind4)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind5)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind6)) + " " + 
    			   OPERA tostring(*(allmags + jj*numsamples + ind7)),
    			   cdata_->print2screenmutex);
        }
    }

    
    }

    int ind0  = std::floor(0.000 * numsamples);
    int ind1  = std::floor(0.025 * numsamples);
    int ind2  = std::floor(0.050 * numsamples);
    int ind3  = std::floor(0.160 * numsamples);
    int ind4  = std::floor(0.300 * numsamples);
    int ind5  = std::floor(0.400 * numsamples);
    int ind6  = std::floor(0.500 * numsamples);
    int ind7  = std::floor(0.600 * numsamples);
    int ind8  = std::floor(0.700 * numsamples);
    int ind9  = std::floor(0.840 * numsamples);
    int ind10 = std::floor(0.950 * numsamples);
    int ind11 = std::floor(0.975 * numsamples);
    int ind12 = std::floor(0.99999999 * numsamples);

    for (int jj=0; jj<nummagmasks+1; ++jj) {
      
      std::sort(allmags + jj*numsamples, allmags + (jj+1)*numsamples);

      double writeout[13];
      writeout[0] = *(allmags + jj*numsamples + ind0);
      writeout[1] = *(allmags + jj*numsamples + ind1);
      writeout[2] = *(allmags + jj*numsamples + ind2);
      writeout[3] = *(allmags + jj*numsamples + ind3);
      writeout[4] = *(allmags + jj*numsamples + ind4);
      writeout[5] = *(allmags + jj*numsamples + ind5);
      writeout[6] = *(allmags + jj*numsamples + ind6);
      writeout[7] = *(allmags + jj*numsamples + ind7);
      writeout[8] = *(allmags + jj*numsamples + ind8);
      writeout[9] = *(allmags + jj*numsamples + ind9);
      writeout[10] = *(allmags + jj*numsamples + ind10);
      writeout[11] = *(allmags + jj*numsamples + ind11);
      writeout[12] = *(allmags + jj*numsamples + ind12);

      PRINTER writeoutstream <double> (writeout,13,data_->mag_errs->stream,data_->mag_errs->lock, data_->precision, NULL);
    }
    PRINTER writeoutstream ("",data_->mag_errs->stream,data_->mag_errs->lock, data_->precision, NULL);

    noise->~VECTOR();
    MEMORY ps_free( noise_dummy );
    newmps->~VECTOR();
    MEMORY ps_free( newmps_dummy );

    delete [] allmags;
}

void pixsrc_common_adaptive::actualprintsource (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, VECTOR *noise)
{
    // this function print the source, 1-sigma errors, and a s/n map.

    if (vars_->fatalerror)
    {
        return;
    }

    vars_->mps->update_cpu();

    //PS_SIT div = 10;
    //if (data_->use_shapelets)
    //    div=4;
    //else if(data_->gridtype == 2)
    //div = std::min( (PS_SIT)std::ceil(
    //OPERA round( std::pow(2.0,(double)vars_->grid->size() ) ) /
    //(vars_->zeroethgridsize) ), 10 );

    PS_SIT div = data_->srcpixelscale0;
    if (data_->srcpixelscale)
        div = OPERA round (1.0/(data_->srcpixelscale[0]*data_->arc2pix));
    div = div<1 ? 1 : div;

    double minx2, maxx2, miny2, maxy2;
    OPERA assign_p_infinity( &minx2 );
    OPERA assign_p_infinity( &miny2 );
    OPERA assign_n_infinity( &maxx2 );
    OPERA assign_n_infinity( &maxy2 );

    for(PS_SIT cc=0; cc<vars_->triout->numberofpoints; cc++)
    {
        if(vars_->triout->pointlist[cc*2] < minx2)
            minx2 = vars_->triout->pointlist[cc*2];
        if(vars_->triout->pointlist[cc*2+1] < miny2)
            miny2 = vars_->triout->pointlist[cc*2+1];
        if(vars_->triout->pointlist[cc*2] > maxx2)
            maxx2 = vars_->triout->pointlist[cc*2];
        if(vars_->triout->pointlist[cc*2+1] > maxy2)
            maxy2 = vars_->triout->pointlist[cc*2+1];
    }

    PS_SIT startx = (PS_SIT)std::floor(minx2);
    PS_SIT starty = (PS_SIT)std::floor(miny2);
    PS_SIT endx   = (PS_SIT)std::ceil (maxx2);
    PS_SIT endy   = (PS_SIT)std::ceil (maxy2);

    // NOTE GUG
    // PRINT warning here instead
    /*
      if(startx<0)
      startx=0;
      if(starty<0)
      starty=0;
      if(endx>data_->imgx-1)
      endx=data_->imgx-1;
      if(endy>data_->imgy-1)
      endy=data_->imgy-1;
    */

    double *im;
    double *im_n=0, *im_s2n=0;
    bool passone = 1;
    bool failed;
    do
    {
        if(!passone)
        {
            PRINTER printwarning( data_->print2screenname,
                                  "not enough memory to use full resolution in \"*.fits\".."
                                  " decreasing resolution of \"*.fits\" by ~50%",
                                  cdata_->print2screenmutex                                     );

            div = std::max ((PS_SIT)1, (PS_SIT)std::floor( 0.707107 * div ));
        }
        else
        {
            passone = 0;
        }
        vars_->srcx = (endx-startx)*div+1;
        vars_->srcy = (endy-starty)*div+1;
        vars_->numsrcpoints= vars_->srcx*vars_->srcy;
        im = (double*)malloc(sizeof(double)*vars_->numsrcpoints);
        if (data_->use_shapelets)
        {
            im_n   = (double*)malloc(sizeof(double)*vars_->numsrcpoints);
            im_s2n = (double*)malloc(sizeof(double)*vars_->numsrcpoints);
        }
        failed = im && (!data_->use_shapelets || (im_n && im_s2n));
    } while (!failed);

    std::fill( im, im + vars_->numsrcpoints, 0 );
    if (data_->use_shapelets)
    {
        std::fill( im_n,   im_n   + vars_->numsrcpoints, 0 );
        std::fill( im_s2n, im_s2n + vars_->numsrcpoints, 0 );
    }

    // create working array for shapelets
    double *hermvals_vec = NULL;
    if (data_->use_shapelets)
        MEMORY ps_malloc (&hermvals_vec, vars_->numberofshapelets);

    // get convex hull of source pixels
    PS_SIT numvert;
    double *convexhull;
    GEOM getconvexhull (NULL, vars_->triout, &numvert, &convexhull, NULL);

    for(PS_SIT loop=0; loop<3; ++loop)
    {
        if (loop && !data_->noisemap)
            break;
        if( loop && !(data_->lambdaguess && ( !data_->usersetsrc || data_->reg )) )
            continue;
        if (loop && data_->use_shapelets)
            break;

        PS_SIT triind, triseed=0, xind, yind, indexsh;
        double xx,yy;
        for (xind=0; xind<vars_->srcx; xind+=1)
        {
            for (yind=0; yind<vars_->srcy; yind+=1)
            {
                xx = startx + (double)xind/div;
                yy = starty + (double)yind/div;

                // skip is point is not in convex hull
                if (!GEOM isinpoly(xx, yy, convexhull, numvert))
                    continue;

                indexsh = (vars_->srcy-yind-1)*vars_->srcx+xind;
                triind = GEOM search_triangle (vars_->triout, &triseed, xx, yy);

                if (-1==triind)
                    continue;

                if (!data_->use_shapelets)
                {
                    double pos[2][3];
                    double val[3];
                    for(PS_SIT f=0; f<3; f++)
                    {
                        pos[0][f] = vars_->triout->pointlist[vars_->triout->trianglelist[triind*3+f]*2  ];
                        pos[1][f] = vars_->triout->pointlist[vars_->triout->trianglelist[triind*3+f]*2+1];

                        switch (loop)
                        {
                        case 0: val[f] = vars_->mps->get(vars_->triout->trianglelist[triind*3+f]);
                            break;
                        case 1: val[f] = noise->get(vars_->triout->trianglelist[triind*3+f]);
                            break;
                        case 2: val[f] = vars_->mps->get(vars_->triout->trianglelist[triind*3+f]) /
                                noise->get(vars_->triout->trianglelist[triind*3+f]);
                            break;
                        default: break;
                        }
                    }
                    OPERA planarvalueinterpolation (pos, val, xx, yy, &im[indexsh]);
                }
                else
                {
                    if (data_->use_shapelets==1)
                    {
                        im[indexsh] =
                            SHAPELETSOPERA integrate_square (data_, cdata_, vars_, hermvals_vec,
                                                             vars_->num_shapelets1, vars_->num_shapelets2,
                                                             vars_->shapelet_ctr[0], vars_->shapelet_ctr[1],
                                                             vars_->shapelet_scale,
                                                             xx-1.0/(div*2.0), xx+1.0/(div*2.0),
                                                             yy-1.0/(div*2.0), yy+1.0/(div*2.0),
                                                             vars_->mps, &im_n[indexsh]);

                        // convert from flux to surface brightness
                        im[indexsh]   *= div*div;
                        im_n[indexsh] *= div*div;

                        if (data_->noisemap && im[indexsh] && im_n[indexsh])
                            im_s2n[indexsh] = im[indexsh] / im_n[indexsh];
                    }
                    else if (data_->use_shapelets==2)
                    {
                        PRINTER printerror(data_->print2screenname,
                                           "polar shapelets disabled",
                                           cdata_->print2screenmutex);
                    }
                }
            }
        }

        vars_->reduction = 1.0 / div;
        double wcscoord[2];
        pthread_mutex_lock( cdata_->wcsmutex );
        HEADER getimgwcscoord( data_->wcs, data_->imgy, startx, starty, &wcscoord[0], &wcscoord[1] );
        pthread_mutex_unlock( cdata_->wcsmutex );
        HEADER setsrcwcs( data_, vars_, 0, 0, wcscoord[0], wcscoord[1] );
        string outname="";
        switch(loop)
        {
        case 0: outname = "mps.fits";
            break;
        case 1: outname = "noise.fits";
            break;
        case 2: outname = "s2n.fits";
            break;
        default: break;
        }
        PRINTER printfitssrcplane( vars_->tracker, cdata_->basename, data_->name,
                                   data_->print2screenname, 1,outname,vars_->srcx,
                                   vars_->srcy,im,vars_->wcsinfo, cdata_->print2screenmutex, 1, 0, cdata_->fitsiomutex);
        if (data_->use_shapelets && data_->noisemap)
        {
            outname = "noise.fits";
            PRINTER printfitssrcplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,outname,vars_->srcx,
                                       vars_->srcy,im_n,vars_->wcsinfo, cdata_->print2screenmutex, 1, 0, cdata_->fitsiomutex);
            outname = "s2n.fits";
            PRINTER printfitssrcplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,outname,vars_->srcx,
                                       vars_->srcy,im_s2n,vars_->wcsinfo, cdata_->print2screenmutex, 1, 0, cdata_->fitsiomutex);
        }

    }

    MEMORY ps_free (convexhull);

    if (hermvals_vec)
        MEMORY ps_free (hermvals_vec);

    free(im);
    if (data_->use_shapelets)
    {
        free (im_n);
        free (im_s2n);
    }
}

// This struct compares angles for sorting in COMMONADAPTIVE createlo_s2i
struct angle_compare
{
    angle_compare (PS_SIT im, double *pointlist, D2MATRIX <double> *images, PS_SIT *tracker)
        {
            this->im = im;
            this->pointlist = pointlist;
            this->images = images;
            this->tracker = tracker;
        }
    bool operator () (PS_SIT i, PS_SIT j)
        {
            // i or j is im, then nudge the position of im
            // towards the center of all surrpoints
            double pos1[2] = {0,0}, pos2[2] = {0,0}, ctr[2] = {0,0};
            if (i==im || j==im)
            {
                for (PS_SIT ind=0; ind<images->size (tracker[im]); ++ind)
                {
                    ctr[0] += images->get (tracker[im], ind*2);
                    ctr[1] += images->get (tracker[im], ind*2+1);
                }
                ctr[0] /= images->size (tracker[im]);
                ctr[1] /= images->size (tracker[im]);
                ctr[0] -= pointlist[im*2];
                ctr[1] -= pointlist[im*2+1];
                ctr[0] *= -0.0001;
                ctr[1] *= -0.0001;
            }

            if (i==im)
            {
                pos1[0] = ctr[0];
                pos1[1] = ctr[1];
            }
            else
            {
                pos1[0] = pointlist[i*2]-pointlist[im*2];
                pos1[1] = pointlist[i*2+1]-pointlist[im*2+1];
            }
            if (j==im)
            {
                pos2[0] = ctr[0];
                pos2[1] = ctr[1];
            }
            else
            {
                pos2[0] = pointlist[j*2]-pointlist[im*2];
                pos2[1] = pointlist[j*2+1]-pointlist[im*2+1];
            }

            return atan2 (pos1[1], pos1[0]) < atan2 (pos2[1], pos2[0]);
        }

    PS_SIT im;
    double *pointlist;
    D2MATRIX <double> *images;
    PS_SIT *tracker;
};

void pixsrc_common_adaptive::createlo_s2i (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    // This function constructs the lensing operator by ray-tracing
    // pixels from the source plane to the image plane.
    // The fractional overlap of a source pixel's image with
    // an image pixel determines the weight for that particular
    // source and image pixel.

    // CUDA requires row-then-column ordering in lensing operator matrix.

    // this will hold the initial tabulation of lensed images
    D2MATRIX <double> *dummy1;
    MEMORY ps_malloc (&dummy1, 1);
    D2MATRIX <double> *images = new (dummy1) D2MATRIX <double> (vars_->lonc);

    // This will temporarily store lensing operator.
    D2MATRIX <double> **dummy3;
    D2MATRIX <double> **weights;
    MEMORY ps_malloc (&dummy3, vars_->lonr, 1);
    MEMORY ps_malloc (&weights, vars_->lonr);
    for (PS_SIT r=0; r<vars_->lonr; ++r)
        weights[r] = new (dummy3[r]) D2MATRIX <double> (vars_->lonc);

    // this will hold final list of images and triangulation
    struct triangulateio *triinloc, *trioutloc;
    MEMORY ps_malloc(            &triinloc , 1 );
    MEMORY ps_malloc(            &trioutloc, 1 );
    MEMORY triangulatestructinit( triinloc     );
    MEMORY triangulatestructinit( trioutloc    );

    // this is a reusable array for image finding
    PS_SIT numimg = 0;
    double **imgarr = NULL;

    // begin image finding
    PS_SIT total_img = 0;
    pthread_mutex_lock( cdata_->potdefmagmutex );
    pthread_mutex_lock( cdata_->wcsmutex       );
    for (PS_SIT s=0; s<vars_->lonc; ++s)
    {
        COMMON raytrace_s2i (data_, cdata_,
                             vars_->triout->pointlist[s*2],
                             vars_->triout->pointlist[s*2+1],
                             &numimg, &imgarr);
        total_img += numimg;

        // fortran-like indexing for numerical recipes
        for (PS_SIT im=1; im<=numimg; ++im)
        {
            images->pushback (s, imgarr[im][1]);
            images->pushback (s, imgarr[im][2]);
        }

        // clear up imgarr array
        if (numimg>0)
        {
            free (imgarr[1]);
            free (imgarr);
        }
    }
    pthread_mutex_unlock( cdata_->potdefmagmutex );
    pthread_mutex_unlock( cdata_->wcsmutex       );

    // This will keep track of which source and image number
    // a point in triangulation belongs to
    PS_SIT *tracker;
    MEMORY ps_malloc (&tracker, total_img);

    // fill in images into array for triangulation
    triinloc->numberofpointattributes = 0;
    triinloc->numberofpoints = total_img;
    MEMORY ps_malloc( &(triinloc->pointlist), total_img*2 );

    PS_SIT im_ind = 0;
    for (PS_SIT s=0; s<vars_->lonc; ++s)
    {
        for (PS_SIT im=0; im<images->size(s); im+=2)
        {
            triinloc->pointlist[im_ind*2]   = images->get (s, im);
            triinloc->pointlist[im_ind*2+1] = images->get (s, im+1);
            tracker[im_ind] = s;
            ++im_ind;
        }
    }

    ps_tri_triangulate ((char*)CONSTANT triswitchesnominangle,
                        triinloc, trioutloc, (struct triangulateio*)NULL);

    // for each point p, find all other points that share a triangle with p
    D2MATRIX <PS_SIT> *dummy2;
    MEMORY ps_malloc( &dummy2, 1 );
    D2MATRIX <PS_SIT> *surrpoints = new (dummy2) D2MATRIX<PS_SIT> (total_img);
    PS_SIT thisc, thisc1, thisc2;
    char addit;
    for( PS_SIT t=0; t<trioutloc->numberoftriangles; ++t )
    {
        for( PS_SIT v=0; v<3; ++v )
        {
            thisc = trioutloc->trianglelist[t*3+v];
            thisc1 = v+1;
            thisc2 = v+2;

            if(thisc1>2)
                thisc1 -= 3;
            if(thisc2>2)
                thisc2 -= 3;

            thisc1 = trioutloc->trianglelist[t*3+thisc1];
            thisc2 = trioutloc->trianglelist[t*3+thisc2];

            addit = 1;
            for( PS_SIT cc=0; cc<surrpoints->size(thisc); ++cc )
                if( surrpoints->get(thisc,cc)==thisc1 )
                {
                    addit = 0;
                    break;
                }
            if (addit)
                surrpoints->pushback( thisc, thisc1 );

            addit = 1;
            for( PS_SIT cc=0; cc<surrpoints->size(thisc); ++cc )
                if( surrpoints->get(thisc,cc)==thisc2 )
                {
                    addit = 0;
                    break;
                }
            if (addit)
                surrpoints->pushback( thisc, thisc2 );
        }
    }

    // For points on edge, add them to their own surrpoint list
    // if they haven't been added already.
    PS_SIT edge_i;
    for(PS_SIT j=0; j<trioutloc->numberofedges; j++)
    {
        if(trioutloc->edgemarkerlist[j]==1)
        {
            edge_i = trioutloc->edgelist[j*2];
            if (surrpoints->get (edge_i, surrpoints->size (edge_i)-1) != edge_i)
                surrpoints->pushback(edge_i, edge_i);

            edge_i = trioutloc->edgelist[j*2+1];
            if (surrpoints->get (edge_i, surrpoints->size (edge_i)-1) != edge_i)
                surrpoints->pushback(edge_i, edge_i);
        }
    }

    // For each point, sort the surrounding points into CCW order.
    for (PS_SIT im=0; im<total_img; ++im)
    {
        PS_SIT *this_splist = surrpoints->get_pointer (im);
        std::sort (this_splist, this_splist + surrpoints->size (im),
                   angle_compare (im, triinloc->pointlist, images, tracker));
    }

    // Find fractional overlap of image triangles with image plane pixels
    double minx, miny, maxx, maxy;
    for (PS_SIT im=0; im<total_img; ++im)
    {
        PS_SIT this_numsp = surrpoints->size (im);
        PS_SIT *this_splist = surrpoints->get_pointer (im);
        // create stencil around point im
        double *poly = NULL;
        MEMORY ps_malloc (&poly, this_numsp*2*2);
        PS_SIT ind = 0, pt_prev;

        double pos[2], ctr[2] = {0,0};
        for (PS_SIT pt=0; pt<this_numsp; ++pt)
        {
            pt_prev = pt==0 ? this_numsp-1 : pt-1;

            // add centroid of this triangle
            pos[0] = 1.0/3.0 * (triinloc->pointlist[this_splist[pt]*2] +
                                triinloc->pointlist[this_splist[pt_prev]*2] +
                                triinloc->pointlist[im*2]);
            pos[1] = 1.0/3.0 * (triinloc->pointlist[this_splist[pt]*2+1] +
                                triinloc->pointlist[this_splist[pt_prev]*2+1] +
                                triinloc->pointlist[im*2+1]);
            poly[ind*2]   = pos[0];
            poly[ind*2+1] = pos[1];
            ++ind;

            // add midpoint along this line segment
            pos[0] = 0.5 * (triinloc->pointlist[this_splist[pt]*2] +
                            triinloc->pointlist[im*2]);
            pos[1] = 0.5 * (triinloc->pointlist[this_splist[pt]*2+1] +
                            triinloc->pointlist[im*2+1]);
            if (im==this_splist[pt])
            {
                ctr[0] = ctr[1] = 0;
                for (PS_SIT pt=0; pt<this_numsp; ++pt)
                {
                    ctr[0] += triinloc->pointlist[this_splist[pt]*2];
                    ctr[0] += triinloc->pointlist[this_splist[pt]*2+1];
                }
                ctr[0] /= this_numsp;
                ctr[1] /= this_numsp;
                ctr[0] -= triinloc->pointlist[im*2];
                ctr[1] -= triinloc->pointlist[im*2+1];
                ctr[0] *= 0.0001;
                ctr[1] *= 0.0001;
                pos[0] -= ctr[0];
                pos[1] -= ctr[1];
            }
            poly[ind*2]   = pos[0];
            poly[ind*2+1] = pos[1];
            ++ind;
        }

        // compute smallest box of image pixels that covers images
        OPERA assign_p_infinity (&minx);
        OPERA assign_p_infinity (&miny);
        OPERA assign_n_infinity (&maxx);
        OPERA assign_n_infinity (&maxy);
        for (PS_SIT pt=0; pt<this_numsp*2; ++pt)
        {
            if (poly[pt*2]  <minx) minx = poly[pt*2];
            if (poly[pt*2+1]<miny) miny = poly[pt*2+1];
            if (poly[pt*2]  >maxx) maxx = poly[pt*2];
            if (poly[pt*2+1]>maxy) maxy = poly[pt*2+1];
        }
        minx = std::max (std::floor(minx)-1, 0.0);
        miny = std::max (std::floor(miny)-1, 0.0);
        maxx = std::min (std::ceil(maxx)+1, data_->imgx-1.0);
        maxy = std::min (std::ceil(maxy)+1, data_->imgy-1.0);

        // loop over box of image pixels and compute areas
        PS_SIT num_vert = 0;
        double *clipped = NULL;

        double sumarea = 0;
        for (PS_SIT x=(PS_SIT)minx; x<=(PS_SIT)maxx; ++x)
        {
            for (PS_SIT y=(PS_SIT)miny; y<=(PS_SIT)maxy; ++y)
            {
                PS_SIT r = x*data_->imgy+y;
                if (vars_->r4r[r]==-1)
                    continue;

                double clip_box[8] = {x-0.5, y-0.5, x+0.5, y-0.5,
                                      x+0.5, y+0.5, x-0.5, y+0.5};

                // The Clipper stuff is code from software named
                // clipper.
                Paths subj(1), clip(1), solution;
                for (PS_SIT p=0; p<4; ++p)
                    clip[0] << IntPoint((PS_SIT)(clip_box[p*2]*10000),
                                        (PS_SIT)(clip_box[p*2+1]*10000));
                for (PS_SIT p=0; p<this_numsp*2; ++p)
                    subj[0] << IntPoint((PS_SIT)(poly[p*2]*10000),
                                        (PS_SIT)(poly[p*2+1]*10000));
                Clipper c;
                c.AddPaths(subj, ptSubject, true);
                c.AddPaths(clip, ptClip, true);
                c.Execute(ctIntersection, solution, pftNonZero, pftNonZero);
                if (solution.size()==0)
                    continue;
                num_vert = solution[0].size();
                if (num_vert<3)
                    continue;
                MEMORY ps_malloc (&clipped, num_vert*2);
                for (PS_SIT p=0; p<num_vert; ++p)
                {
                    clipped[p*2]   = solution[0][p].X/10000.0;
                    clipped[p*2+1] = solution[0][p].Y/10000.0;
                }

                double area = GEOM areapoly (clipped, num_vert);
                sumarea += area;

                if (area>0)
                    weights[vars_->r4r[r]]->pushback (tracker[im], area);

                MEMORY ps_free (clipped);
                clipped = NULL;
            }
        }
        MEMORY ps_free (poly);
    }

    // set weights in actual lensing matrix
    for (PS_SIT r=0; r<vars_->lonr; ++r)
    {
        for (PS_SIT s=0; s<vars_->lonc; ++s)
        {
            double sum = 0;
            for (PS_SIT w=0; w<weights[r]->size (s); ++w)
                sum += weights[r]->get (s,w);
            if (sum>0.0001)
            {
                vars_->lensingoperatornobo->set (r,s,sum);
            }
        }
    }

    // memory cleanup
    MEMORY triangulatestructdestruct( triinloc  );
    MEMORY triangulatestructdestruct( trioutloc );
    MEMORY ps_free(                   triinloc  );
    MEMORY ps_free(                   trioutloc );
    MEMORY ps_free (tracker);

    images->~D2MATRIX();
    MEMORY ps_free (dummy1);

    surrpoints->~D2MATRIX();
    MEMORY ps_free (dummy2);

    for (PS_SIT r=0; r<vars_->lonr; ++r)
        weights[r]->~D2MATRIX();
    MEMORY ps_free (weights);
    MEMORY ps_free (dummy3, vars_->lonr);
}

void pixsrc_common_adaptive::createlo_i2s (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    double pos[2][3]; // used in planar interpolation

    double weights3[3]; // used in planar interpolation

    if( data_->subsampling == 1 )
    {
        // this holds the order in which column indices should be
        // set in the lensing operator.
        // CUDA requires row-then-column ordering.
        PS_SIT setorder[3];

        for( PS_SIT r=0; r<vars_->lonr; ++r )
        {
            if( vars_->pointer[vars_->r4rback[r]] != -1 )
            {
                for( PS_SIT f=0; f<3; ++f )
                {
                    pos[0][f] = vars_->triout->pointlist[
                        vars_->triout->trianglelist[
                            vars_->pointer[vars_->r4rback[r]]*3+f]*2  ];
                    pos[1][f] = vars_->triout->pointlist[
                        vars_->triout->trianglelist[
                            vars_->pointer[vars_->r4rback[r]]*3+f]*2+1];
                }

                OPERA planarinterpolation3pts( vars_->newloc[ vars_->r4rback[r] * 2    ],
                                               vars_->newloc[ vars_->r4rback[r] * 2 +1 ],
                                               pos, weights3                             );

                // recording column positions in setorder
                PS_SIT setpos = 0;
                for(PS_SIT f=0; f<3; f++)
                {
                    // keeping at least 99.99% of the flux
                    if( std::abs(weights3[f])>0.0001 )
                        setorder[setpos++] = vars_->triout->trianglelist[
                            vars_->pointer[ vars_->r4rback[r] ]*3 + f ];
                }

                // sorting column indices
                std::sort( setorder, setorder + setpos );

                // setting lensing operator
                for( PS_SIT sp=0; sp<setpos; ++sp )
                {
                    for( PS_SIT f=0; f<3; ++f )
                        if( setorder[sp] == vars_->triout->trianglelist[
                                vars_->pointer[ vars_->r4rback[r] ]*3 + f ] )
                        {
                            vars_->lensingoperatornobo->set(
                                r,vars_->triout->trianglelist[
                                    vars_->pointer[ vars_->r4rback[r] ]*3 + f ],
                                weights3[f] );
                            break;
                        }
                }
            }
        }
    }
    else
    {
        double ss2 = 1.0/(data_->subsampling*data_->subsampling);
        double *val, *valcopy;
        PS_SIT *setorder;
        MEMORY ps_malloc (&val,      vars_->lonc); // stores lensing operator weights
        MEMORY ps_malloc (&valcopy , vars_->lonc); // stores lensing operator weights
        MEMORY ps_malloc (&setorder, vars_->lonc);

        for(PS_SIT r=0; r<vars_->lonr; r++)
        {
            std::fill(val,val+vars_->lonc,0);
            for(PS_SIT x=0; x<data_->subsampling; x++)
            {
                for(PS_SIT y=0; y<data_->subsampling; y++)
                {
                    if(vars_->sspointer[vars_->r4rback[r]][x*data_->subsampling+y]!=-1)
                    {
                        for(PS_SIT f=0; f<3; f++)
                        {
                            pos[0][f] =
                                vars_->triout->pointlist[ vars_->triout->trianglelist[
                                    vars_->sspointer[vars_->r4rback[r]][x*data_->subsampling+y]*3+f ]*2   ];
                            pos[1][f] =
                                vars_->triout->pointlist[ vars_->triout->trianglelist[
                                    vars_->sspointer[vars_->r4rback[r]][x*data_->subsampling+y]*3+f ]*2+1 ];
                        }
                        OPERA planarinterpolation3pts(
                            vars_->newloc_sslo[ vars_->r4rback[r] ][(x*data_->subsampling+y)*2],
                            vars_->newloc_sslo[ vars_->r4rback[r] ][(x*data_->subsampling+y)*2+1], pos,weights3 );
                        for(PS_SIT f=0; f<3; f++)
                            val[vars_->triout->trianglelist[
                                    vars_->sspointer[ vars_->r4rback[r] ][x*data_->subsampling+y]*3+f]] +=
                                weights3[f];
                    }
                }
            }

            std::transform (val, val+vars_->lonc, val,
                            std::bind1st (std::multiplies<PS_FPT>(),ss2));
            std::copy (val, val+vars_->lonc, valcopy);
            std::sort (valcopy, valcopy+vars_->lonc);

            double runningsum = 0, totalsum = 0;
            double *it = valcopy + vars_->lonc-1;
            for (PS_SIT c=0; c<vars_->lonc; ++c)
            {
                if (!*it)
                    break;
                totalsum += *it;
                --it;
            }
            if (!totalsum)
                continue;

            PS_SIT minpos;
            it = valcopy+vars_->lonc-1;
            for( minpos=vars_->lonc-1; minpos>=0; --minpos )
            {
                // keeping at least 99.99% of the flux
                if( runningsum > 0.9999*totalsum )
                    break;
                runningsum += *it;
                --it;
            }
            if (-1==minpos)
                minpos = 0;
            double min = valcopy[minpos];

            it = val;
            for (PS_SIT c=0; c<vars_->lonc; ++c)
            {
                if (*it>=min && !OPERA equalszero(*it))
                    vars_->lensingoperatornobo->set (r, c, *it);
                ++it;
            }
        }

        MEMORY ps_free(   val    );
        MEMORY ps_free( valcopy  );
        MEMORY ps_free( setorder );
    }

    //vars_->lensingoperatornobo->creatercform();
}

void pixsrc_common_adaptive::createlo(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    PS_SIT dim = 1;
    PS_SIT waitlist[1] = {7};
    OPERA pthreadswait( vars_, dim, waitlist );

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "creating lensing operator",
                             cdata_->print2screenmutex);

    // initial lensing matrix
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonr*3 : vars_->lonr*3;

    MEMORY ps_malloc( &(vars_->lensingoperatornoboptr), 1 );

    vars_->lensingoperatornobo = new (vars_->lensingoperatornoboptr) MATRIX(
        cdata_, data_, vars_->lonr, vars_->lonc, initsize, 0, data_ );


    // fill in lensing operator based on user settings
    // shapelets expansion
    if (data_->use_shapelets)
    {
        SHAPELETSOPERA createlo_shapelets (data_, cdata_, vars_);
    }
    else if (0==data_->raydirection)
    {
        // based on the warren & dye method
        COMMONADAPTIVE createlo_s2i (data_, cdata_, vars_);
    }
    else
    {
        // based on koopmans method
        COMMONADAPTIVE createlo_i2s (data_, cdata_, vars_);
    }


    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "done creating lensing operator",
                             cdata_->print2screenmutex);

    // convolve image plane with PSF
    COMMON blurit(data_, cdata_, vars_);
}
