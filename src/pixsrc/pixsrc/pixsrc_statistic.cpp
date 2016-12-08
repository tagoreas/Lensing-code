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



#include "pixsrc_statistic.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_tps.hpp"
#include "pixsrc_cuda.hpp"

#include <cmath>
#include <iomanip>

bool trieqon2 = 0;

void* pixsrc_statistic::computees(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[4] = 1;

    STATISTIC computeesbody( data_, cdata_, vars_ );

    vars_->pthreadstracker[4] = 2;

    return NULL;
}

void* pixsrc_statistic::computeed(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[3] = 1;

    STATISTIC computeedbody( data_, cdata_, vars_ );

    vars_->pthreadstracker[3] = 2;

    return NULL;
}

void* pixsrc_statistic::computedeta(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[2] = 1;

    STATISTIC computedetabody( data_, cdata_, vars_ );

    vars_->pthreadstracker[2] = 2;

    return NULL;
}

void* pixsrc_statistic::computedetc(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[1] = 1;

    STATISTIC computedetcbody( data_, cdata_, vars_ );

    vars_->pthreadstracker[1] = 2;

    return NULL;
}

void* pixsrc_statistic::evidencecalculator(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[11] = 1;

    STATISTIC evidencecalculatorbody( data_, cdata_, vars_ );

    vars_->pthreadstracker[11] = 2;

    return NULL;
}

void pixsrc_statistic::computeinterperrors(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (!data_->interperr || !vars_->tpsfit)
    {
        return;
    }

    if(data_->verbose!=1)
        PRINTER print2screen(data_->print2screenname,
                             "integrating over TPS for interpolation errors",
                             cdata_->print2screenmutex);

    if (!vars_->interperrptr)
    {
        MEMORY ps_malloc (&vars_->interperrptr, 1);
        vars_->interperr = new (vars_->interperrptr) VECTOR (cdata_, data_, vars_->lonr);
    }
    if (!vars_->ieanalyticptr)
    {
        MEMORY ps_malloc (&vars_->ieanalyticptr, 1);
        vars_->ieanalytic = new (vars_->ieanalyticptr) VECTOR (cdata_, data_, vars_->lonr);
    }

    // get "true" lens model
    if (1==data_->interperr)
        ((TPS*)(vars_->tpsfit))->interpolate (vars_->newloc,      vars_->r4r, vars_->ieanalytic);
    else
        ((TPS*)(vars_->tpsfit))->interpolate (vars_->newloc_ssas, vars_->r4r, vars_->ieanalytic, data_->interperr);

    // get interpolation errors
    VECTOR *dummy;
    MEMORY ps_malloc( &dummy, 1 );
    VECTOR *lensedmpsnobo = new (dummy) VECTOR(cdata_, data_, vars_->lonr);
    vars_->lensingoperatornobo->mult (vars_->mps, lensedmpsnobo, 0, cdata_->numthreads);
    vars_->ieanalytic->minus (lensedmpsnobo, NULL, 1, 1);
    lensedmpsnobo->~VECTOR();
    MEMORY ps_free (dummy);

    // blur interpolation errors (or don't)
    if (!data_->nopsf)
    {
        if (data_->lowmem)
            PRINTER printerror (data_->print2screenname,
                                "cannot use interpolation errors with PSF and low memory on",
                                cdata_->print2screenmutex);

        vars_->interperr->zeromeout ();
        vars_->blurringoperator->mult (vars_->ieanalytic, vars_->interperr, 0, cdata_->numthreads);
    }
    else
    {
        VECTOR *temp2, *temp2ptr;
        temp2                = vars_->interperr;
        temp2ptr             = vars_->interperrptr;
        vars_->interperr     = vars_->ieanalytic;
        vars_->interperrptr  = vars_->ieanalyticptr;
        vars_->ieanalytic    = temp2;
        vars_->ieanalyticptr = temp2ptr;
    }

    if(data_->verbose!=1)
        PRINTER print2screen(data_->print2screenname,
                             "done integrating over TPS for interpolation errors",
                             cdata_->print2screenmutex);
}

void pixsrc_statistic::actualcomputeed (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    // this has been given its own function so that source parameter optimization
    // can call this routine without the rest of computeed(..)

    if(data_->noevinochi)
    {
        vars_->ed=1;
    }
    else
    {
        if (data_->is_uvdata && 1==data_->transmode)
        {
            vars_->truelonr = vars_->lonr;
            vars_->dofreduction = 0;

            VECTOR uvtmpvec (cdata_, data_, vars_->lonc);
            vars_->b1->mult (vars_->mps, &uvtmpvec, 0, cdata_->numthreads);

            vars_->ed = uvtmpvec.innerproduct (vars_->mps) +
                data_->uv_transform_scalar - 2.0 *vars_->mps->innerproduct (vars_->lod);
            return;
        }

        // if using interpolation errors in computing statistic
        if (data_->interperr && !data_->irscheme && vars_->tpsfit)
            STATISTIC computeinterperrors (data_, cdata_, vars_);

        vars_->ed = 0;

        if (data_->is_uvdata && data_->uv_newvariance)
        {
            vars_->truelonr = vars_->lonr;
            vars_->dofreduction = 0;

            // if doing uv-plane reconstruction with varying noise
            for (PS_SIT s=0; s<data_->uv_ndp*2; ++s)
            {
                double tmp = vars_->residual4stat->get(s);
                vars_->ed += tmp*tmp / data_->uv_newvariance[s/2];
            }
        }
        else if( data_->chi2mask )
        {
            vars_->truelonr = vars_->lonr;
            double temp, temperr;
            for(PS_SIT r=0; r<vars_->lonr; r++)
            {
                if(data_->chi2mask[vars_->r4rback[r]])
                {
                    temp = vars_->residual4stat->get(r);
                    if (data_->interperr && !data_->irscheme && vars_->tpsfit)
                    {
                        temperr = vars_->interperr->get(r);
                        vars_->ed += temp*temp / (temperr*temperr + vars_->datavariance);
                    }
                    else
                    {
                        vars_->ed += temp*temp;
                    }
                }
                else
                {
                    vars_->truelonr--;
                }
            }
        }
        else
        {
            if (data_->interperr && !data_->irscheme && vars_->tpsfit)
            {
                double temp, temperr;
                for(PS_SIT r=0; r<vars_->lonr; r++)
                {
                    temp    = vars_->residual4stat ->get(r);
                    temperr = vars_->interperr->get(r);
                    vars_->ed += temp*temp / (temperr*temperr + vars_->datavariance );
                }
            }
            else
            {
                vars_->ed = vars_->residual4stat->innerproduct(vars_->residual4stat);
            }

            vars_->truelonr = vars_->lonr;
        }

        if(data_->onlychi)
        {
            // if regularizing and have already optimized source (if doing so)
            if( data_->reg && !( data_->usersetsrc && !vars_->srcopt )  )
            {
                vars_->dofreduction = 0;
                /*
                  vars_->dofreduction = vars_->lonc -
                  vars_->lambda1*vars_->a1->
                  trace_invA_B( vars_->c1, cdata_->numthreads, &vars_->fatalerror );
                */
            }
            else
            {
                vars_->dofreduction = 0;
            }
        }

        // imposing penalties for sister images
        if(data_->penaltyquery[5])
        {
            //vars_->ed += vars_->residual4statout->innerproduct(vars_->residual4statout);
            //vars_->truelonr += vars_->lonrout;
        }

        // if the variance hasn't already been factored in
        if (!(data_->interperr && !data_->irscheme && vars_->tpsfit) &&
            !(data_->is_uvdata && data_->uv_newvariance))
            vars_->ed /= vars_->datavariance;
    }
}

void pixsrc_statistic::evidencecalculatorbody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        OPERA assign_n_infinity( &data_->evidence );
        return;
    }

    // wait for majority of threads to complete
    PS_SIT dim = 7;
    PS_SIT waitlist[7] = {1,2,3,4,5,6,10};
    OPERA pthreadswait(vars_,dim,waitlist);

    if(data_->verbose!=1)
        PRINTER print2screen(data_->print2screenname,
                             "calculating statistic",
                             cdata_->print2screenmutex);

    double detCdTerm=0,detATerm=0,detCTerm=0,EdTerm=0,EsTerm=0,lambdaTerm=0,PiTerm=0,evidence;

    if (vars_->fatalerror)
    {
        OPERA assign_n_infinity( &data_->evidence );
        return;
    }
    else if(vars_->terminatelensing)
    {
        evidence=0;
    }
    else if(0&&cdata_->npl && cdata_->npl_stepsize)
    {
        vars_->lambda_lens = data_->lambda_lens<0 ? -data_->lambda_lens : data_->lambda_lens;

        double pen[3] = {0,0,0};
        pen[0] = vars_->ed;
        if( vars_->lambda1 )
            pen[1] = 1e2* NONPARAMLENS gsl_src_penalty (data_, cdata_, vars_);
        if( vars_->lambda_lens )
            pen[2] = vars_->lambda_lens * NONPARAMLENS gsl_lens_penalty (data_, cdata_);

        evidence += -0.5 * (pen[0] + pen[1] + pen[2]);

//        std::cout << pen[0] << " " << pen[1] << " " << pen[2] << std::endl;
    }
    else if(data_->noevinochi)
    {
        evidence = 0;
        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "no chi^2 or evidence requested",
                                 cdata_->print2screenmutex);
    }
    else if(data_->onlychi)
    {
        // the (-1/2) factor is to compensate for how Chuck's code turns evidence into chi2
        // the dof reduction term is there b/c regularization reduces the effective dof.
        //evidence = -0.5*vars_->ed/(vars_->truelonr-vars_->dofreduction); //reduced chi^2
        evidence = -0.5*vars_->ed; //traditional chi^2
        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "chi^2 = " + OPERA tostring(-2.0*evidence),
                                 cdata_->print2screenmutex);
    }
    else
    {
        // calculating detCd term
        if (data_->is_uvdata && data_->uv_newvariance)
        {
            // if doing uv-plane reconstruction with varying noise
            detCdTerm = data_->uv_newvariance[0];
            //for (PS_SIT s=0; s<data_->uv_ndp*2; ++s)
            //    detCdTerm += std::log (data_->uv_newvariance[s/2]);
            detCdTerm *= -0.5;
        }
        else if( data_->interperr && !data_->irscheme )
        {
            if( data_->chi2mask )
            {
                for( PS_SIT r=0; r<vars_->lonr; ++r )
                {
                    if( data_->chi2mask[vars_->r4rback[r]] )
                        detCdTerm += std::log( vars_->interperr->get(r) *
                                               vars_->interperr->get(r) + vars_->datavariance );
                }
                detCdTerm *= -0.5;
            }
            else
            {
                for( PS_SIT r=0; r<vars_->lonr; ++r )
                {
                    detCdTerm += std::log( vars_->interperr->get(r) *
                                           vars_->interperr->get(r) + vars_->datavariance );
                }
                detCdTerm *= -0.5;
            }
        }
        else
        {
            detCdTerm = -0.5*vars_->truelonr*std::log(vars_->datavariance);
        }

        detATerm = -vars_->logdeta/2.0;
        detCTerm =  vars_->logdetc/2.0;
        EdTerm = -vars_->ed/2.0;
        EsTerm = -vars_->lambda1*vars_->es/2.0;
        lambdaTerm = vars_->lonc/2.0*std::log(vars_->lambda1);
        PiTerm = -vars_->truelonr/2.0*std::log(2*CONSTANT pi);
        evidence = detATerm + detCTerm + detCdTerm + EsTerm + EdTerm + lambdaTerm + PiTerm;

        //std::cout << detATerm << " " << detCTerm << " " << detCdTerm << " " << EdTerm << " " << EsTerm << " " << lambdaTerm << " " << PiTerm << std::endl;

        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "-2 ln(evidence) = " + OPERA tostring(-2.0*evidence),
                                 cdata_->print2screenmutex);
    }

    if(cdata_->npl && cdata_->npl_stepsize)
    {
        vars_->lambda_lens = data_->lambda_lens<0 ? -data_->lambda_lens : data_->lambda_lens;

        double pen[3] = {0,0,0};
        //pen[0] = vars_->ed;
        //if( vars_->lambda1 )
        //  pen[1] = 1e2* NONPARAMLENS gsl_src_penalty (data_, cdata_, vars_);
        if( vars_->lambda_lens )
            pen[2] = vars_->lambda_lens * NONPARAMLENS gsl_lens_penalty (data_, cdata_);

        evidence += -0.5 * (pen[0] + pen[1] + pen[2]);

        //std::cout << pen[0] << " " << pen[1] << " " << pen[2] << std::endl;
    }


    double pensum=0;
    for(PS_SIT j=0; j<data_->extlengths[5]; j++)
        pensum += vars_->penalties[j];

    evidence += -0.5 * (pensum);

    //std::cout << -2.0*evidence << std::endl;

    data_->evidence=evidence;

    if(data_->traceparams[2]>0)
    {
        PS_SIT dim = 9+data_->extlengths[5];
        double writeout[dim];
        writeout[0] = vars_->lambda1;
        writeout[1] = evidence;
        writeout[2] = detCdTerm;
        writeout[3] = detATerm;
        writeout[4] = detCTerm;
        writeout[5] = EdTerm;
        writeout[6] = EsTerm;
        writeout[7] = lambdaTerm;
        writeout[8] = PiTerm;
        for(PS_SIT j=0; j<data_->extlengths[5]; j++)
            writeout[9+j] = vars_->penalties[j];


        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "trace #" + OPERA tostring(vars_->tracker) +
                                 ": " + OPERA tostring(writeout[0]) + " => " +
                                 OPERA tostring(writeout[1]),
                                 cdata_->print2screenmutex);

        PRINTER writeoutstream(writeout,dim,data_->traces->stream,data_->traces->lock, data_->precision, NULL);
    }
}

void pixsrc_statistic::computedetabody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->noevinochi || data_->onlychi || !data_->reg )
    {
        vars_->logdeta=1;
    }
    else
    {
        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "calculating log(det(A))",
                                 cdata_->print2screenmutex);
        vars_->logdeta = vars_->a1->logdet();
    }
}

void pixsrc_statistic::computedetcbody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->noevinochi || data_->onlychi || !data_->reg )
    {
        vars_->logdetc=1;
    }
    else
    {
        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "calculating log(det(C))",
                                 cdata_->print2screenmutex);
        //      vars_->logdetc = 2*vars_->h1->logdet();
        vars_->logdetc = vars_->c1->logdet();
    }
}

void pixsrc_statistic::computeesbody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->noevinochi || data_->onlychi || !data_->reg )
    {
        vars_->es=1;
    }
    else if( data_->regorder==0 && data_->gridtype==0 )
    {
        vars_->es = vars_->mps->innerproduct(vars_->mps);
    }
    else
    {
        VECTOR *dummy;
        MEMORY ps_malloc( &dummy, 1 );
        VECTOR *temp = new (dummy) VECTOR( cdata_, data_, vars_->lonc);

        vars_->c1->mult( vars_->mps, temp, 0, cdata_->numthreads );
        vars_->es = vars_->mps->innerproduct( temp );

        temp->~VECTOR();
        MEMORY ps_free( dummy );

        if( vars_->es < 0 && !OPERA equalszero(vars_->es) )
        {
            if(data_->verbose!=1)
                PRINTER print2screen(data_->print2screenname,
                                     "fatal error: Es = " + OPERA tostring(vars_->es),
                                     cdata_->print2screenmutex);
            ++vars_->fatalerror;
            OPERA assign_p_infinity( &vars_->es );
        }
    }
}

void pixsrc_statistic::computeedbody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_ )
{
    if (vars_->fatalerror)
    {
        return;
    }

    actualcomputeed (data_, cdata_, vars_);
}
