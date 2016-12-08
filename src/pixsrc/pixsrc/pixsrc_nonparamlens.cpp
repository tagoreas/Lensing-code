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



#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_external.hpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_geometry.hpp"
#include <cmath>

double get_evidence (inputdata*, commoninputdata*);

double pixsrc_nonparamlens::npl_penalty (const gsl_vector *v, void *params)
{
    // setup data
    gsl_mdm_params *params_ = (gsl_mdm_params*)params;
    inputdata        *data_ = params_-> data;
    commoninputdata *cdata_ = params_->cdata;
    npl_struct         *npl = data_->data_->npl_st;
    PS_SIT             numvary = npl->numvary;

    // copy potential perturbations for pixsrc
    for (PS_SIT i=0; i<numvary; ++i)
        npl->b_vec[i] = gsl_vector_get (v, npl->index[i]);

    // get TPS control point weights
    NONPARAMLENS get_tps_weights (data_, cdata_, gsl_vector_get (v, npl->index[numvary]));

    // do source reconstruction stuff
    double evi = -2 * get_evidence (data_, cdata_);
    //std::cout << "npl penalty = " << evi << std::endl;
    return evi;
}

double pixsrc_nonparamlens::gsl_src_penalty (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    double pen = 0;
    VECTOR *dummy;
    MEMORY ps_malloc (&dummy, 1);
    VECTOR *temp = new (dummy) VECTOR (cdata_, data_, vars_->lonc);

    vars_->c1->mult (vars_->mps, temp, 0, cdata_->numthreads);
    pen = vars_->mps->innerproduct (temp);

    temp->~VECTOR();
    MEMORY ps_free (dummy);

    return pen;
}

double pixsrc_nonparamlens::gsl_lens_penalty (inputdata *data_, commoninputdata *cdata_)
{
    // only compute penalty once (using first image)
    if (data_ != data_->data_)
        return 0;

    // npl_struct *npl = data_->data_->npl_st;

    double pen = 0;
    for (PS_SIT r=0; r<data_->ndp; ++r)
        if (data_->imagemasks[r]==2)
            pen += NONPARAMLENS get_grad_kappa (data_, cdata_,
                                                data_->oldloc_arc[r*2],
                                                data_->oldloc_arc[r*2+1]);

    // also add penalty for regularizing the TPS
    // pen += gsl_vector_get (npl->x, npl->index[npl->numvary]) * pen;

    return pen;
}

void pixsrc_nonparamlens::minimize (inputdata *data_, commoninputdata *cdata_)
{
    npl_struct *npl = data_->data_->npl_st;

    struct gsl_mdm mdm_;
    struct gsl_mdm *mdm = &mdm_;

    bool printit = 0;
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
        printit = printit || npl->printvec0[g];
    
    if( cdata_->npl_stepsize)
    {
        // set initial step sizes
        mdm->ss = gsl_vector_alloc (npl->numvary+1);
        gsl_vector_set_all (mdm->ss, cdata_->npl_stepsize);
        // this is a special stepsize for tps regularization
        gsl_vector_set (mdm->ss, npl->index[npl->numvary], cdata_->npl_reg_parms[1]);
        //gsl_vector_set (mdm->ss, npl->index[npl->numvary], 0);

        // more init
        mdm->minex_func.n = npl->numvary+1;
        mdm->minex_func.f = NONPARAMLENS npl_penalty;
        mdm->minex_func.params = npl->params;

        mdm->s = gsl_multimin_fminimizer_alloc (mdm->T, npl->numvary+1);

        gsl_multimin_fminimizer_set( mdm->s, &mdm->minex_func, npl->x, mdm->ss );

        PS_SIT print_temp = 1;
        do
        {
            ++mdm->iter;

            mdm->status = gsl_multimin_fminimizer_iterate (mdm->s);

            if (mdm->status)
                break;

            mdm->size = gsl_multimin_fminimizer_size (mdm->s);
            mdm->status = gsl_multimin_test_size (mdm->size, cdata_->npl_ftolsize);

            // if minimizer converged
            if (mdm->status == GSL_SUCCESS )
            {
                printf ("converged after %d evals to minimum at\n",
                        (int)cdata_->numcall);

                // copy best model to x
                gsl_vector *best = gsl_multimin_fminimizer_x (mdm->s);
                for (PS_SIT nv=0; nv<npl->numvary+1; ++nv)
                    gsl_vector_set (npl->x, npl->index[nv],
                                    gsl_vector_get (best, npl->index[nv]));
            }

            // periodically print out best model
            if (cdata_->numcall/100>=print_temp)
            {
                ++print_temp;
                gsl_vector *best = gsl_multimin_fminimizer_x (mdm->s);
                double **print_out;
                PS_SIT dimfile = 3;
                MEMORY ps_malloc( &print_out, npl->numvary+1, dimfile );
                for( PS_SIT r=0; r<npl->numvary; ++r )
                {
                    print_out[r][0] = npl->tps_grid[r*2];
                    print_out[r][1] = npl->tps_grid[r*2+1];
                    print_out[r][2] = gsl_vector_get (best, npl->index[r]);
                }
                print_out[npl->numvary][0] = print_out[npl->numvary][1] = 0;
                print_out[npl->numvary][2] = gsl_vector_get (best, npl->index[npl->numvary]);
                PRINTER print <double> (0, cdata_->basename,
                                        data_->name, true, "npl_nonopt_tps.dat",
                                        print_out, npl->numvary+1, dimfile, 0,
                                        data_->precision, NULL);
                MEMORY ps_free (print_out, npl->numvary+1);
            }

            printf ("%5d f() = %7.3f size = %.3f\n",
                    (int)mdm->iter, mdm->s->fval, mdm->size );
        }
        while (mdm->status == GSL_CONTINUE);

        if (mdm->status != GSL_SUCCESS)
        {
            PRINTER printerror ("pixsrc",
                                "non-parametric lens model minimization failed",
                                cdata_->print2screenmutex);
        }

        // cleanup
        gsl_vector_free (mdm->ss);
        gsl_multimin_fminimizer_free (mdm->s);
    }
    else if (!printit)
    {
        NONPARAMLENS npl_penalty (npl->x, npl->params);
    }
}

void pixsrc_nonparamlens::print_model (inputdata *data_, commoninputdata *cdata_)
{
    npl_struct *npl = data_->data_->npl_st;

    bool printit = 0;
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
        printit = printit || npl->printvec0[g];

    if (printit)
    {
        for (PS_SIT g=0; g<cdata_->numimages; ++g)
        {
            data_[g].printvec = npl->printvec0[g];
            data_[g].doreconstruction =
                !data_[g].noevinochi ||
                data_[g].penaltyquery[0] ||
                data_[g].penaltyquery[1] ||
                data_[g].penaltyquery[5] ||
                data_[g].magparams ||
                data_[g].printvec ? 1 : 0;
        }
        NONPARAMLENS npl_penalty (npl->x, npl->params);
    }
    else
    {
        // copy potential perturbations for pixsrc
        for (PS_SIT i=0; i<npl->numvary; ++i)
            npl->b_vec[i] = gsl_vector_get (npl->x, npl->index[i]);

        // do perturbation stuff
        NONPARAMLENS get_tps_weights (data_, cdata_,
                                      gsl_vector_get (npl->x, npl->index[npl->numvary]));
    }

    if (1 || cdata_->npl_stepsize)
    {
        // print out model
        PS_SIT dimfile = 11;
        double **print_out, pot, alpha[2], kappa, grad_kappa;
        MEMORY ps_malloc( &print_out, data_->ndp, dimfile );
        for( PS_SIT r=0; r<data_->ndp; ++r )
        {
            pot = NONPARAMLENS get_pot (data_, cdata_,
                                        data_->oldloc_arc[r*2],
                                        data_->oldloc_arc[r*2+1]);
            NONPARAMLENS get_alpha (data_, cdata_,
                                    data_->oldloc_arc[r*2],
                                    data_->oldloc_arc[r*2+1],
                                    &alpha[0], &alpha[1]);
            kappa = NONPARAMLENS get_kappa (data_, cdata_,
                                            data_->oldloc_arc[r*2],
                                            data_->oldloc_arc[r*2+1]);
            grad_kappa = NONPARAMLENS get_grad_kappa (data_, cdata_,
                                                      data_->oldloc_arc[r*2],
                                                      data_->oldloc_arc[r*2+1]);
            print_out[r][0] = data_->oldloc[r*2  ];
            print_out[r][1] = data_->oldloc[r*2+1];
            print_out[r][2] = data_->oldloc_wcs[r*2  ];
            print_out[r][3] = data_->oldloc_wcs[r*2+1];
            print_out[r][4] = data_->oldloc_arc[r*2  ];
            print_out[r][5] = data_->oldloc_arc[r*2+1];
            print_out[r][6] = pot;
            print_out[r][7] = kappa;
            print_out[r][8] = grad_kappa;
            print_out[r][9] = alpha[0];
            print_out[r][10] = alpha[1];
        }
        PRINTER print <double> (0, cdata_->basename,
                                data_->name, true, "npl.dat",
                                print_out, data_->ndp, dimfile, 0, data_->precision, NULL);
        MEMORY ps_free (print_out, data_->ndp);

        dimfile = 3;
        MEMORY ps_malloc( &print_out, npl->numvary+1, dimfile );
        for( PS_SIT r=0; r<npl->numvary; ++r )
        {
            print_out[r][0] = npl->tps_grid[r*2];
            print_out[r][1] = npl->tps_grid[r*2+1];
            print_out[r][2] = gsl_vector_get (npl->x, npl->index[r]);
        }
        print_out[npl->numvary][0] = print_out[npl->numvary][1] = 0;
        print_out[npl->numvary][2] = gsl_vector_get (npl->x, npl->index[npl->numvary]);
        PRINTER print <double> (0, cdata_->basename,
                                data_->name, true, "npl_tps.dat",
                                print_out, npl->numvary+1, dimfile, 0, data_->precision, NULL);
        MEMORY ps_free (print_out, npl->numvary+1);
    }

    // cleanup
    gsl_vector_free (npl->x);
    MEMORY ps_free (npl->params);
}

void pixsrc_nonparamlens::init_npl (inputdata *data_, commoninputdata *cdata_)
{
    MEMORY ps_malloc (&data_->data_->npl_st, 1);
    npl_struct *npl = data_->data_->npl_st;
    MEMORY ps_malloc (&npl->printvec0, cdata_->numimages);

    // turn off reference perturbations (initially)
    npl->dim_ref = npl->numvary_ref = 0;

    // temporarily store user-set parameters and change default values
    // so we don't waste time doing these during minimization
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
    {
        npl->printvec0[g] = data_[g].printvec;
        data_[g].printvec = 0;
        data_[g].doreconstruction =
            !data_[g].noevinochi ||
            data_[g].penaltyquery[0] ||
            data_[g].penaltyquery[1] ||
            data_[g].penaltyquery[5] ||
            data_[g].magparams ||
            data_[g].printvec ? 1 : 0;
    }

    // create struct to hold parameters for GSL function minimizer
    MEMORY ps_malloc (&npl->params, 1);
    npl->params->data  = data_;
    npl->params->cdata = cdata_;
    npl->numvary = 0;

    // if exclusively using tps grid from a previous run, read tps file
    char **gridvec=0;
    char ***vecsplit=0;
    PS_SIT vecsplitnum=0;
    double spacing_ra=0, spacing_dec=0, start_ra=0, start_dec=0;
    if (cdata_->npl==2)
    {
        // read file
        char *fname;
        const char *listcc[2];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = cdata_->npl_filename[0];
        OPERA concatenate (listcc, 2, &fname);
        OPERA readfile (fname, &gridvec, &npl->numvary, cdata_->print2screenmutex);
        MEMORY ps_free (fname);

        // parse lines of file
        MEMORY ps_malloc (&vecsplit, npl->numvary);
        for (PS_SIT v=0; v<npl->numvary; ++v)
            OPERA split (gridvec[v], " ", &vecsplit[v], &vecsplitnum);
        MEMORY ps_free (gridvec, npl->numvary);
    }
    else
    {
        if (cdata_->npl_num_pts[0]<2 || cdata_->npl_num_pts[1]<2)
            PRINTER printerror ("pixsrc",
                                "cannot have potential perturbation grid with "
                                "fewer than 2 points in R.A. or Dec.",
                                cdata_->print2screenmutex);

        // derived parameters
        spacing_ra  = cdata_->npl_size[0] / (cdata_->npl_num_pts[0]-1);
        spacing_dec = cdata_->npl_size[1] / (cdata_->npl_num_pts[1]-1);
        start_ra  = cdata_->npl_ctr[0] - cdata_->npl_size[0]/2.0;
        start_dec = cdata_->npl_ctr[1] - cdata_->npl_size[1]/2.0;

        // calculate number of points in tps grid that overlap with image grid
        double pos0[2], pos1[2];
        pthread_mutex_lock (cdata_->wcsmutex);
        for (PS_SIT m=0; m<cdata_->npl_num_pts[0]*cdata_->npl_num_pts[1]; ++m)
        {
            pos0[0] = start_ra  + (m/cdata_->npl_num_pts[1])*spacing_ra;
            pos0[1] = start_dec + (m%cdata_->npl_num_pts[1])*spacing_dec;
            HEADER getwcssfromwcsl (pos0[0], pos0[1], data_->pra, data_->pdec,
                                    data_->r1, data_->r2, pos1);
            HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, pos1[0], pos1[1],
                                   &pos0[0], &pos0[1]);
            PS_SIT r = OPERA round (pos0[0])*data_->imgy + OPERA round (pos0[1]);

            if (pos0[0]>0             && pos0[1]>0 &&
                pos0[0]<data_->imgx-1 && pos0[1]<data_->imgy-1 &&
                data_->imagemasks[r]==2)
            {
                ++npl->numvary;
            }
        }
        pthread_mutex_unlock (cdata_->wcsmutex);
    }

    if (npl->numvary<3)
        PRINTER printerror ("pixsrc",
                            "cannot have fewer than 3 TPS control points.",
                            cdata_->print2screenmutex);
    else
        PRINTER print2screen ("pixsrc",
                              "found " + OPERA tostring (npl->numvary) +
                              " TPS control points.",
                              cdata_->print2screenmutex);

    // init
    npl->trans = 'N';
    npl->dim = npl->numvary+3;
    npl->nrhs = 1;
    npl->LDA = npl->dim;
    npl->LDB = npl->dim;
    MEMORY ps_malloc (&npl->tps_grid,     npl->numvary*2);
    MEMORY ps_malloc (&npl->tps_mat,      npl->dim*npl->dim);
    MEMORY ps_malloc (&npl->tps_mat_fac,  npl->dim*npl->dim);
    MEMORY ps_malloc (&npl->x_vec,        npl->dim);
    MEMORY ps_malloc (&npl->b_vec,        npl->dim);
    MEMORY ps_malloc (&npl->ipiv,         npl->dim);
    MEMORY ps_malloc (&npl->index,        npl->numvary+1);
    for (PS_SIT ind=0; ind<npl->numvary+1; ++ind)
	npl->index[ind] = ind;
    // keep regularization strength at last index
    OPERA shuffle (npl->index, npl->numvary/*+1*/);

    // this vector holds initial guess.
    npl->x = gsl_vector_alloc (npl->numvary+1);

    // partially fill in tps rhs vector
    npl->b_vec[npl->dim-1] = npl->b_vec[npl->dim-2] = npl->b_vec[npl->dim-3] = 0;

    // calculate tps grid positions
    if (cdata_->npl==2)
    {
        double minx, maxx, miny, maxy;
        OPERA assign_p_infinity (&minx);
        OPERA assign_p_infinity (&miny);
        OPERA assign_n_infinity (&maxx);
        OPERA assign_n_infinity (&maxy);

        for (PS_SIT v=0; v<npl->numvary; ++v)
        {
            npl->tps_grid[v*2]   = OPERA convert_string <double> (vecsplit[v][0]);
            npl->tps_grid[v*2+1] = OPERA convert_string <double> (vecsplit[v][1]);
	    if (3==vecsplitnum)
		gsl_vector_set (npl->x, npl->index[v], OPERA convert_string <double> (vecsplit[v][2]));

            if (npl->tps_grid[v*2]<minx)
                minx = npl->tps_grid[v*2];
            if (npl->tps_grid[v*2]>maxx)
                maxx = npl->tps_grid[v*2];
            if (npl->tps_grid[v*2+1]<miny)
                miny = npl->tps_grid[v*2+1];
            if (npl->tps_grid[v*2+1]>maxy)
                maxy = npl->tps_grid[v*2+1];
        }
	if (3==vecsplitnum)
	    gsl_vector_set (npl->x, npl->index[npl->numvary], -15);
        npl->pert_region[0] = minx;
        npl->pert_region[1] = miny;
        npl->pert_region[2] = minx;
        npl->pert_region[3] = maxy;
        npl->pert_region[4] = maxx;
        npl->pert_region[5] = maxy;
        npl->pert_region[6] = maxx;
        npl->pert_region[7] = miny;
        MEMORY ps_free (vecsplit, npl->numvary, vecsplitnum);
    }
    else
    {
        pthread_mutex_lock (cdata_->wcsmutex);
        PS_SIT tpsindex = 0;
        double pos0[2], pos1[2];
        for (PS_SIT m=0; m<cdata_->npl_num_pts[0]*cdata_->npl_num_pts[1]; ++m)
        {
            pos0[0] = start_ra  + (m/cdata_->npl_num_pts[1])*spacing_ra;
            pos0[1] = start_dec + (m%cdata_->npl_num_pts[1])*spacing_dec;
            HEADER getwcssfromwcsl (pos0[0], pos0[1], data_->pra, data_->pdec,
                                    data_->r1, data_->r2, pos1);
            HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, pos1[0], pos1[1],
                                   &pos0[0], &pos0[1]);
            PS_SIT r = OPERA round (pos0[0])*data_->imgy + OPERA round (pos0[1]);

            if (m==0)
            {
                npl->pert_region[0] = pos0[0];
                npl->pert_region[1] = pos0[1];
            }
            else if (m==cdata_->npl_num_pts[1]-1)
            {
                npl->pert_region[2] = pos0[0];
                npl->pert_region[3] = pos0[1];
            }
            else if (m==cdata_->npl_num_pts[0]*cdata_->npl_num_pts[1]-1)
            {
                npl->pert_region[4] = pos0[0];
                npl->pert_region[5] = pos0[1];
            }
            else if (m==cdata_->npl_num_pts[0]*cdata_->npl_num_pts[1] -
                     cdata_->npl_num_pts[1])
            {
                npl->pert_region[6] = pos0[0];
                npl->pert_region[7] = pos0[1];
            }

            if (pos0[0]>0             && pos0[1]>0 &&
                pos0[0]<data_->imgx-1 && pos0[1]<data_->imgy-1 &&
                data_->imagemasks[r]==2)
            {
                npl->tps_grid[tpsindex*2]   = start_ra  + (m/cdata_->npl_num_pts[1])*spacing_ra;
                npl->tps_grid[tpsindex*2+1] = start_dec + (m%cdata_->npl_num_pts[1])*spacing_dec;
                ++tpsindex;
            }
        }
        pthread_mutex_unlock (cdata_->wcsmutex);
    }

    // fill in tps matrices
    // have to fill in column-major order
    // Since it's symmetric, row-major (needed for CUDA) works.
    npl->meandist2 = 0;
    PS_SIT index = 0;
    double pos_m[2], pos_n[2], r2;
    for (PS_SIT m=0; m<npl->numvary; ++m)
    {
        pos_m[0] = npl->tps_grid[m*2];
        pos_m[1] = npl->tps_grid[m*2+1];

        for (PS_SIT n=0; n<npl->numvary; ++n)
        {
            pos_n[0] = npl->tps_grid[n*2];
            pos_n[1] = npl->tps_grid[n*2+1];

            r2 = OPERA distance2 (pos_m[0], pos_m[1], pos_n[0], pos_n[1]);
            npl->tps_mat[index++] = OPERA equalszero(r2) ? 0 : r2 * 0.5*std::log10 (r2);

            if (m!=n)
                npl->meandist2 += std::sqrt (r2);
        }

        npl->tps_mat[index++] = 1;
        npl->tps_mat[index++] = pos_m[0];
        npl->tps_mat[index++] = pos_m[1];
    }
    std::fill (npl->tps_mat+index, npl->tps_mat+index+npl->numvary, 1);
    index += npl->numvary;
    std::fill (npl->tps_mat+index, npl->tps_mat+index+3, 0);
    index += 3;
    for (PS_SIT n=0; n<npl->numvary; ++n)
    {
        pos_n[0] = npl->tps_grid[n*2];
        npl->tps_mat[index++] = pos_n[0];
    }
    std::fill (npl->tps_mat+index, npl->tps_mat+index+3, 0);
    index += 3;
    for (PS_SIT n=0; n<npl->numvary; ++n)
    {
        pos_n[1] = npl->tps_grid[n*2+1];
        npl->tps_mat[index++] = pos_n[1];
    }
    std::fill (npl->tps_mat+index, npl->tps_mat+index+3, 0);
    npl->meandist2 /= (npl->numvary*npl->numvary-npl->numvary);
    npl->meandist2 *= npl->meandist2;

//    else if (cdata_->npl==1)
    {
        MEMORY ps_malloc (&npl->dim_ref,      cdata_->npl_num_ref_pot);
        MEMORY ps_malloc (&npl->numvary_ref,  cdata_->npl_num_ref_pot);
        MEMORY ps_malloc (&npl->x_vec_ref,    cdata_->npl_num_ref_pot);
        MEMORY ps_malloc (&npl->tps_grid_ref, cdata_->npl_num_ref_pot);

        // read in reference potential perturbations
        for (PS_SIT rp=0; rp<cdata_->npl_num_ref_pot; ++rp)
        {
            // read file
            char *fname, **gridvec2, ***vecsplit2;
            PS_SIT numvary0, info0;
            const char *listcc[2];
            listcc[0] = CONSTANT dir_in;
            listcc[1] = cdata_->npl==1 ? cdata_->npl_filename[rp] :
		cdata_->npl_filename[rp+1];
            OPERA concatenate (listcc, 2, &fname);
            OPERA readfile (fname, &gridvec2, &numvary0, cdata_->print2screenmutex);
            MEMORY ps_free (fname);

            // parse lines of file
            PS_SIT dummy;
            MEMORY ps_malloc (&vecsplit2, numvary0);
            for (PS_SIT v=0; v<numvary0; ++v)
                OPERA split (gridvec2[v], "\t ", &vecsplit2[v], &dummy);
            MEMORY ps_free (gridvec2, numvary0);

            double **tps0;
            MEMORY ps_malloc (&tps0, numvary0, 3);
            for (PS_SIT v=0; v<numvary0; ++v)
            {
                tps0[v][0] = OPERA convert_string <double> (vecsplit2[v][0]);
                tps0[v][1] = OPERA convert_string <double> (vecsplit2[v][1]);
                tps0[v][2] = OPERA convert_string <double> (vecsplit2[v][2]);
            }
            MEMORY ps_free (vecsplit2, numvary0, 3);
            --numvary0;

            // fill in tps matrices
            // have to fill in column-major order
            // Since it's symmetric, row-major (needed for CUDA) works.
            double *tps_mat0, *x_vec0;
            int/*PS_SIT*/ *ipiv0;
            PS_SIT dim0 = numvary0+3;
            MEMORY ps_malloc (&tps_mat0,      dim0*dim0);
            MEMORY ps_malloc (&x_vec0,        dim0);
            MEMORY ps_malloc (&ipiv0,         dim0);

            double meandist2 = 0;
            index = 0;
            for (PS_SIT m=0; m<numvary0; ++m)
            {
                for (PS_SIT n=0; n<numvary0; ++n)
                {
                    r2 = OPERA distance2 (tps0[m][0], tps0[m][1], tps0[n][0], tps0[n][1]);
                    tps_mat0[index++] = OPERA equalszero(r2) ? 0 : r2 * 0.5*std::log10 (r2);

                    if (m!=n)
                        meandist2 += std::sqrt (r2);
                }

                tps_mat0[index++] = 1;
                tps_mat0[index++] = tps0[m][0];
                tps_mat0[index++] = tps0[m][1];
            }
            std::fill (tps_mat0+index, tps_mat0+index+numvary0, 1);
            index += numvary0;
            std::fill (tps_mat0+index, tps_mat0+index+3, 0);
            index += 3;
            for (PS_SIT n=0; n<numvary0; ++n)
            {
                tps_mat0[index++] = tps0[n][0];
            }
            std::fill (tps_mat0+index, tps_mat0+index+3, 0);
            index += 3;
            for (PS_SIT n=0; n<numvary0; ++n)
            {
                tps_mat0[index++] = tps0[n][1];
            }
            std::fill (tps_mat0+index, tps_mat0+index+3, 0);
            meandist2 /= (numvary0*numvary0-numvary0);
            meandist2 *= meandist2;

            for (PS_SIT m=0; m<numvary0; ++m)
            {
                tps_mat0[m*dim0+m] = meandist2* std::pow (10,tps0[numvary0][2]);
                x_vec0[m] = tps0[m][2];
            }
            x_vec0[dim0-1] = x_vec0[dim0-2] = x_vec0[dim0-3] = 0;

            // do LU factorization of tps matrix
            EXTERNAL ps_dgetrf_ (&dim0, &dim0, tps_mat0,
                                 &dim0, ipiv0, &info0);

            // solve for tps weights
            EXTERNAL ps_dgetrs_ (&npl->trans, &dim0, &npl->nrhs, tps_mat0,
                                 &dim0, ipiv0, x_vec0,
                                 &dim0, &info0);

            // save tps weights and grid
            MEMORY ps_malloc (&npl->tps_grid_ref[rp], numvary0*2);
            for (PS_SIT v=0; v<numvary0; ++v)
            {
                npl->tps_grid_ref[rp][v*2]   = tps0[v][0];
                npl->tps_grid_ref[rp][v*2+1] = tps0[v][1];
            }
            npl->numvary_ref[rp] = numvary0;
            npl->dim_ref[rp]     = dim0;
            npl->x_vec_ref[rp]   = x_vec0;

            MEMORY ps_free (tps0);
            MEMORY ps_free (tps_mat0);
            MEMORY ps_free (ipiv0);

            /*
            // determine initial control points
            double val, x, y;
            for (PS_SIT k=0; k<npl->numvary; ++k)
            {
            x = npl->tps_grid[k*2];
            y = npl->tps_grid[k*2+1];
            val = x_vec0[dim0-3] +
            x_vec0[dim0-2]*x + x_vec0[dim0-1]*y;

            for (PS_SIT v=0; v<numvary0; ++v)
            {
            r2 = OPERA distance2 (x, y,
            tps0[v][0],
            tps0[v][1]);

            if (!OPERA equalszero (r2))
            val += x_vec0[v]* r2* std::log10(r2)/2.0;
            }
            gsl_vector_set (npl->x, k, npl->index[val]);
            }
            gsl_vector_set (npl->x, npl->numvary, npl->index[0]);

            MEMORY ps_free (tps0);
            MEMORY ps_free (tps_mat0);
            MEMORY ps_free (x_vec0);
            MEMORY ps_free (ipiv0);
            */
        }
    }

    // set initial gsl vector to zero except possibly regularization strength
    if (3!=vecsplitnum)
    {
	gsl_vector_set_all (npl->x, 0);
	gsl_vector_set     (npl->x, npl->index[npl->numvary], cdata_->npl_reg_parms[0]);
    }
}

void pixsrc_nonparamlens::get_tps_weights (inputdata *data_, commoninputdata *cdata_, double exp)
{
    npl_struct *npl = data_->data_->npl_st;

    // figure out regularization strength
    npl->reg = std::pow (10.0, exp);

    // copy tps matrix and set regularization strengths
    std::copy (npl->tps_mat, npl->tps_mat + npl->dim*npl->dim,
               npl->tps_mat_fac);
    for (PS_SIT m=0; m<npl->numvary; ++m)
        npl->tps_mat_fac[m*npl->dim+m] = npl->meandist2*npl->reg;

    // copy vector into which solutions will be stored
    std::copy (npl->b_vec, npl->b_vec + npl->dim, npl->x_vec);

    // do LU factorization of tps matrix
    EXTERNAL ps_dgetrf_ (&npl->dim, &npl->dim, npl->tps_mat_fac,
                         &npl->LDA, npl->ipiv, &npl->info);

    // solve for tps weights
    EXTERNAL ps_dgetrs_ (&npl->trans, &npl->dim, &npl->nrhs, npl->tps_mat_fac,
                         &npl->LDA, npl->ipiv, npl->x_vec,
                         &npl->LDB, &npl->info);
}

// I should check to see if rx/ry/r2 are equal to zero,
// but this is highly unlikely.
void pixsrc_nonparamlens::get_alpha (inputdata *data_, commoninputdata *cdata_,
                                     double x, double y, double *dx, double *dy)
{
    npl_struct *npl = data_->data_->npl_st;
    double rx, ry, r2, fac;
    *dx = *dy = 0;

    // add deflections from reference perturbations
    for (PS_SIT rp=0; rp<cdata_->npl_num_ref_pot; ++rp)
    {
        *dx += npl->x_vec_ref[rp][npl->dim_ref[rp]-2];
        *dy += npl->x_vec_ref[rp][npl->dim_ref[rp]-1];
        for (PS_SIT v=0; v<npl->numvary_ref[rp]; ++v)
        {
            rx = x - npl->tps_grid_ref[rp][v*2];
            ry = y - npl->tps_grid_ref[rp][v*2+1];
            r2 = rx*rx + ry*ry;
            fac = (1.0 + (OPERA equalszero(r2) ? 0 : std::log (r2))) / std::log (10.0);

            *dx += npl->x_vec_ref[rp][v] * rx * fac;
            *dy += npl->x_vec_ref[rp][v] * ry * fac;
        }
    }

    // deflections
    *dx += npl->x_vec[npl->dim-2];
    *dy += npl->x_vec[npl->dim-1];
    for (PS_SIT v=0; v<npl->numvary; ++v)
    {
        rx = x - npl->tps_grid[v*2];
        ry = y - npl->tps_grid[v*2+1];
        r2 = rx*rx + ry*ry;
        fac = (1.0 + (OPERA equalszero(r2) ? 0 : std::log (r2))) / std::log (10.0);

        *dx += npl->x_vec[v] * rx * fac;
        *dy += npl->x_vec[v] * ry * fac;
    }
}

double pixsrc_nonparamlens::get_pot (inputdata *data_, commoninputdata *cdata_,
                                     double x, double y)
{
    npl_struct *npl = data_->data_->npl_st;
    double r2, pot=0;


    // add potential from reference perturbations
    for (PS_SIT rp=0; rp<cdata_->npl_num_ref_pot; ++rp)
    {
        pot += npl->x_vec_ref[rp][npl->dim_ref[rp]-3] +
            npl->x_vec_ref[rp][npl->dim_ref[rp]-2]*x +
            npl->x_vec_ref[rp][npl->dim_ref[rp]-1]*y;

        for (PS_SIT v=0; v<npl->numvary_ref[rp]; ++v)
        {
            r2 = OPERA distance2 (x, y,
                                  npl->tps_grid_ref[rp][v*2],
                                  npl->tps_grid_ref[rp][v*2+1]);

            pot += npl->x_vec_ref[rp][v]* r2* (OPERA equalszero(r2) ? 0 : std::log10(r2))/2.0;
        }
    }

    // potential
    pot += npl->x_vec[npl->dim-3] +
        npl->x_vec[npl->dim-2]*x + npl->x_vec[npl->dim-1]*y;

    for (PS_SIT v=0; v<npl->numvary; ++v)
    {
        r2 = OPERA distance2 (x, y,
                              npl->tps_grid[v*2],
                              npl->tps_grid[v*2+1]);

        pot += npl->x_vec[v]* r2* (OPERA equalszero(r2) ? 0 : std::log10(r2))/2.0;
    }

    return pot;
}

double pixsrc_nonparamlens::get_kappa (inputdata *data_, commoninputdata *cdata_,
                                       double x, double y)
{
    npl_struct *npl = data_->data_->npl_st;
    double r2, fac, kappa=0;

    // add convergence from reference perturbations
    for (PS_SIT rp=0; rp<cdata_->npl_num_ref_pot; ++rp)
    {
        for (PS_SIT v=0; v<npl->numvary_ref[rp]; ++v)
        {
            r2 = OPERA distance2 (x, y,
                                  npl->tps_grid_ref[rp][v*2],
                                  npl->tps_grid_ref[rp][v*2+1]);
            fac = 2.0*(2.0 + (OPERA equalszero(r2) ? 0 : std::log(r2))) / std::log(10.0);

            kappa += npl->x_vec_ref[rp][v] * fac;
        }
    }

    // convergence
    for (PS_SIT v=0; v<npl->numvary; ++v)
    {
        r2 = OPERA distance2 (x, y,
                              npl->tps_grid[v*2],
                              npl->tps_grid[v*2+1]);
        fac = 2.0*(2.0 + (OPERA equalszero(r2) ? 0 : std::log(r2))) / std::log(10.0);

        kappa += npl->x_vec[v] * fac;
    }

    return kappa;
}

double pixsrc_nonparamlens::get_grad_kappa (inputdata *data_, commoninputdata *cdata_,
                                            double x, double y)
{
    npl_struct *npl = data_->data_->npl_st;
    double r2, fac, grad_kappa=0;

    // add gradient from reference perturbations
    for (PS_SIT rp=0; rp<cdata_->npl_num_ref_pot; ++rp)
    {
        for (PS_SIT v=0; v<npl->numvary_ref[rp]; ++v)
        {
            r2 = OPERA distance2 (x, y,
                                  npl->tps_grid_ref[rp][v*2],
                                  npl->tps_grid_ref[rp][v*2+1]);
            fac = OPERA equalszero(r2) ? 0 : 4.0 * std::sqrt(1/r2) / std::log (10.0);

            grad_kappa += std::abs (npl->x_vec_ref[rp][v]) * fac;
        }
    }

    // magnitude of gradient of convergence
    for (PS_SIT v=0; v<npl->numvary; ++v)
    {
        r2 = OPERA distance2 (x, y,
                              npl->tps_grid[v*2],
                              npl->tps_grid[v*2+1]);
        fac = OPERA equalszero(r2) ? 0 : 4.0 * std::sqrt(1/r2) / std::log (10.0);

        grad_kappa += std::abs (npl->x_vec[v]) * fac;
    }

    return grad_kappa;
}
