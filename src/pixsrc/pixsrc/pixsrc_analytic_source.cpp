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



#include "pixsrc_analytic_source.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_common.hpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_geometry.hpp"
#include <cmath>

typedef struct
{
    dataandvars *dav;
    PS_SIT *paramindex;
} amoebastruct;

void pixsrc_analytic_source::getsource( inputdata* data_, commoninputdata* cdata_, lensvar* vars_ )
{
    if( !data_->usersetsrc )
        return;


    double *bestsourceoverall;
    MEMORY ps_malloc( &bestsourceoverall, data_->srcnvary );
    double bestchioverall;
    OPERA assign_p_infinity (&bestchioverall);
    PS_SIT restart = data_->srcrestart;



    // get best guess for position of source (only handles 1 source currently)
    // changes all source centers to this value
    double centroidx = 0, centroidy = 0;
    double ellip = 0, ellpa = 0;
    if( data_->guess_src_position )
    {
        //
        // this method uses brightness-weighted image pixel positions
        //
        double weight, f_total = 0;
        for( PS_SIT r=0; r<data_->ndp; ++r )
            if( vars_->r4r[r] != -1 )
            {
                weight = std::abs(data_->data[r]);
                centroidx += vars_->newloc[r*2  ] * weight;
                centroidy += vars_->newloc[r*2+1] * weight;
                f_total += weight;
            }
        centroidx /= f_total;
        centroidy /= f_total;

	f_total = 0;
	double qxx=0, qxy=0, qyy=0;
	for( PS_SIT r=0; r<data_->ndp; ++r )
	    if( vars_->r4r[r] != -1 )
	    {
		weight = std::abs(data_->data[r]);
		qxx += weight *(vars_->newloc[r*2+0]-centroidx)*(vars_->newloc[r*2+0]-centroidx);
		qxy += weight *(vars_->newloc[r*2+0]-centroidx)*(vars_->newloc[r*2+1]-centroidy);
		qyy += weight *(vars_->newloc[r*2+1]-centroidy)*(vars_->newloc[r*2+1]-centroidy);
		f_total += weight;
	    }
	qxx /= f_total;
	qxy /= f_total;
	qyy /= f_total;
	if (qxx*qyy-qxy*qxy<0)
	    qxy = std::sqrt(qxx*qyy);
	double denom = qxx+qyy+2.0*std::sqrt(qxx*qyy-qxy*qxy + CONSTANT smallnumber);
	double eps1 = (qxx-qyy)/denom;
	double eps2 = 2.0*qxy/denom;
	double tmp = std::sqrt(eps1*eps1+eps2*eps2);
	double q = (1.0-tmp)/(1.0+tmp);
	ellip = 1.0-q;
	ellpa = std::atan2 (eps2,eps1)*0.5;

        // convert centroid from pixel coordinates to WCS coordinates
        double ra0, dec0, pos[2];
        pthread_mutex_lock (cdata_->wcsmutex);
        HEADER getimgwcscoord( data_->wcs, data_->imgy, centroidx, centroidy, &ra0, &dec0 );
        HEADER getwcslfromwcss( ra0, dec0, data_->pra, data_->pdec, data_->r1, data_->r2, pos );
        pthread_mutex_unlock (cdata_->wcsmutex);
        centroidx = pos[0];
        centroidy = pos[1];

	ellpa = (ellpa - CONSTANT piby2 + data_->rollangle) /CONSTANT deg2rad;

    	if (ellip<0)
	    ellip = 0;
	while (ellpa<-90)
	    ellpa += 180;
	while (ellpa>90)
	    ellpa -= 180;
}


    while( restart )
    {
        if( data_->srcnvary )
        {
            gsl_vector *x  = gsl_vector_alloc (data_->srcnvary);
            gsl_vector *ss = gsl_vector_alloc (data_->srcnvary);
            gsl_vector_set_all (x, 0);
            gsl_vector_set_all (ss, 0);

            PS_SIT *paramindex;
            MEMORY ps_malloc (&paramindex, data_->srcnvary);

            PS_SIT tracker, index = 0;

            for( PS_SIT src=0; src<data_->usersetsrc[0]; ++src )
            {
                // switching profile type
                switch( (PS_SIT)data_->usersetsrc[1+src*17] )
                {
                    // if profile type = "sersic"
                case 0:
                {
                    tracker = 1;
                    for( PS_SIT pp=src*17+1+9; pp<src*17+1+17; ++pp )
                    {
                        if( data_->usersetsrc[pp] )
                        {
                            paramindex[index] = tracker;
                            gsl_vector_set (x, index, data_->usersetsrc[pp-8]);


                            // auto-position the centroid of source (only 1 source currently)
                            // changes all source centers to this value
                            if( data_->guess_src_position )
                            {
                                if( tracker==2 )
                                    gsl_vector_set (x, index, centroidx);
                                if( tracker==3 )
                                    gsl_vector_set (x, index, centroidy);
                                if( tracker==4 )
                                    gsl_vector_set (x, index, ellip);
                                if( tracker==5 )
                                    gsl_vector_set (x, index, ellpa);
                            }

                            if (tracker==1 || tracker==2 || tracker==3 || tracker==6)
                                gsl_vector_set (ss, index,
                                                data_->srcstepsizes ? data_->srcstepsizes[index] : 1.0);
                            else if( tracker==4 )
                                gsl_vector_set (ss, index,
                                                data_->srcstepsizes ? data_->srcstepsizes[index] : 0.25);
                            else if( tracker==5 )
                                gsl_vector_set (ss, index,
                                                data_->srcstepsizes ? data_->srcstepsizes[index] : 45);
                            else if( tracker==8 )
                                gsl_vector_set (ss, index,
                                                data_->srcstepsizes ? data_->srcstepsizes[index] : 3.0);
                            ++index;
                        }
                        ++tracker;
                    }
                    break;
                }
                // if profile type = "none"
                case -1:
                {
                    tracker = 1;
                    for( PS_SIT pp=src*17+1+9; pp<src*17+1+17; ++pp )
                    {
                        if( data_->usersetsrc[pp] )
                        {
                            paramindex[index] = tracker;
                            gsl_vector_set (x, index, data_->usersetsrc[pp-8]);
                            gsl_vector_set (ss, index, 0);
                            ++index;
                        }
                        ++tracker;
                    }
                    break;
                }
                default:
                {
                    // if profile type = "vector image"
                    if ((PS_SIT)data_->usersetsrc[1+src*17]>1000)
                    {
                        tracker = 1;
                        for( PS_SIT pp=src*17+1+9; pp<src*17+1+17; ++pp )
                        {
                            if( data_->usersetsrc[pp] )
                            {
                                paramindex[index] = tracker;
                                gsl_vector_set (x, index, data_->usersetsrc[pp-8]);

                                if (tracker==1 || tracker==2 || tracker==3 || tracker==4 || tracker==5 || tracker==6 || tracker==7)
                                    gsl_vector_set (ss, index,
                                                    data_->srcstepsizes ? data_->srcstepsizes[index] : 1.0);
                                ++index;
                            }
                            ++tracker;
                        }
                    }

                    break;
                }
                }
            }

            amoebastruct as;
            as.dav = vars_->dav;
            as.paramindex = paramindex;

            // perturb startup model if doing multiple minimizations
            if (restart != data_->srcrestart)
                for (PS_SIT nv=0; nv<data_->srcnvary; ++nv)
                    if (0!=gsl_vector_get(ss,nv))
                    {
                        double myran;
                        do myran = 0.1 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
                        while (std::abs(myran)>0.15);

                        gsl_vector_set (x, nv, gsl_vector_get(x,nv)
                                        + gsl_vector_get(ss,nv)*myran);
                    }

            if(data_->verbose!=1)
            {
                PRINTER print2screen(data_->print2screenname,
                                     "optimizing source parameters # " + 
				     OPERA tostring(data_->srcrestart-restart+1) + " of " +
				     OPERA tostring(data_->srcrestart),
                                     cdata_->print2screenmutex);
            }

            double bestchi;
            OPERA multidimmin (data_, cdata_, data_->srcnvary, data_->srcftol,
                               "source parameters", x, ss, &as, lenssourceandreturnchi2,
                               &bestchi, 5000, 25);

            // if this is the first run through or we did better
            if( restart == data_->srcrestart  ||
                (restart != data_->srcrestart && bestchi < bestchioverall) )
            {
                for (PS_SIT nv=0; nv<data_->srcnvary; ++nv)
                    bestsourceoverall[nv] = gsl_vector_get (x, nv);
                bestchioverall = bestchi;
            }

            //vars_->residual ->zeromeout();
            //vars_->lensedmps->zeromeout();

            MEMORY ps_free (paramindex);
            gsl_vector_free (x);
            gsl_vector_free (ss);
        }

        // if on last run through while loop
        if( restart == 1 )
        {
            ANALYTICSRC setusersource( data_, cdata_, vars_, bestsourceoverall,
                                       NULL, data_->srcbounds, data_->numsrcbounds,
                                       data_->usersetsrc, vars_->mps );

            // save best model
            ANALYTICSRC expandsrc( data_, cdata_, vars_, bestsourceoverall,
                                   &vars_->bestanalyticsrc, data_->usersetsrc );

            // print source to file
            if (data_->printvec)
            {
                // get source parameters with vary flags
                double *srclocflags = NULL;
                vector <string> outvec;
                string thisline;
                ANALYTICSRC expandsrcwithflags (data_, cdata_, vars_, bestsourceoverall,
                                                &srclocflags, data_->usersetsrc);
                // get source type and parameters
                for (PS_SIT j=0; j<srclocflags[0]; ++j)
                {
                    // first get source type
                    thisline = "";
                    switch ((PS_SIT)srclocflags[1+j*17])
                    {
                    case 0:
                        thisline += "sersic ";
                        break;
                    case -1:
                        thisline += "none ";
                        break;
                    default:
                        if (srclocflags[1+j*17]>1000)
                            thisline += "vector ";
                        else
                            PRINTER printerror(data_->print2screenname,
                                               "unrecognized source type",
                                               cdata_->print2screenmutex);
                        break;
                    }
                    // now get parameters
                    for (PS_SIT p=0; p<8; ++p)
                        thisline += OPERA tostring (srclocflags[1+j*17+1+p]) + " ";
                    outvec.push_back (thisline);
                }
                // get source vary flags
                for (PS_SIT j=0; j<srclocflags[0]; ++j)
                {
                    thisline = " ";
                    for (PS_SIT p=0; p<8; ++p)
                        thisline += OPERA tostring ((PS_SIT)srclocflags[1+j*17+1+8+p]) + " ";;
                    outvec.push_back (thisline);
                }
                // print and free memory
                PRINTER print <string> (vars_->tracker, cdata_->basename, data_->name, true,
                                        "srcopt.dat",outvec, 0, data_->precision, NULL);
                MEMORY ps_free (srclocflags);
            }

            if( data_->srcrestart>1 && data_->verbose!=1 && data_->srcnvary )
            {
                string bs = "";
                for(PS_SIT j=0; j<data_->srcnvary; ++j)
                {
                    bs += OPERA tostring(bestsourceoverall[j]);
                    bs += " ";
                }
                bs += ": " + OPERA tostring(bestchioverall);

                PRINTER print2screen(data_->print2screenname,
                                     "best source overall: " + bs,
                                     cdata_->print2screenmutex);
            }
        }

        --restart;
    }

    vars_->srcopt = 1;

    MEMORY ps_free( bestsourceoverall );
}

double pixsrc_analytic_source::lenssourceandreturnchi2(const gsl_vector *gslv, void *args)
{
    amoebastruct *as = (amoebastruct*)args;
    dataandvars *dav = as->dav;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;
    PS_SIT *paramindex = as->paramindex;

    // copy gsl vector
    double *source;
    MEMORY ps_malloc (&source, data_->srcnvary);
    for (PS_SIT nv=0; nv<data_->srcnvary; ++nv)
        source[nv] = gsl_vector_get (gslv, nv);

    double *srclocflags = 0;
    ANALYTICSRC expandsrcwithflags( data_, cdata_, vars_, source, &srclocflags, data_->usersetsrc );
    double penalty = 0;
    double emax = 0.999;
    PS_SIT srctype=-2, sum_vary;
    for( PS_SIT i=0; i<data_->srcnvary; ++i )
    {
        sum_vary=-1;
        for (PS_SIT j=0; j<srclocflags[0]; ++j)
            for (PS_SIT k=1+j*17+1+8; k<1+j*17+1+8+8; ++k)
                if (srclocflags[k])
                {
                    ++sum_vary;
                    if (sum_vary == i)
                    {
                        srctype = srclocflags[1+j*17];
                        j = srclocflags[0];
                        break;
                    }
                }

        switch (srctype)
        {
            // sersic case
        case 0:
        {
            if( paramindex[i]==1 && source[i]<0 )
                penalty += 1.0e100 * ( 1.0 + source[i]*source[i] );
            if( paramindex[i]==4 && std::abs(source[i])<0 )
                penalty += 1.0e100 * ( 1.0 + source[i]*source[i] );
            if( paramindex[i]==4 && std::abs(source[i])>emax )
                penalty += 1.0e100 * ( 1.0 + (source[i]-emax)*(source[i]-emax) );
            if( paramindex[i]==6 && source[i]<0 )
                penalty += 1.0e100 * ( 1.0 + source[i]*source[i] );
            if( paramindex[i]==8 && source[i]<0 )
                penalty += 1.0e100 * ( 1.0 + source[i]*source[i] );
            break;
        }
        // none case
        case -1:
        {
            break;
        }
        default:
        {
            // vector image case
            if (srctype>1000)
            {
                if( paramindex[i]==1 && source[i]<0 )
                    penalty += 1.0e100 * ( 1.0 + source[i]*source[i] );
            }
            else
            {
                PRINTER printerror(data_->print2screenname,
                                   "unrecognized source type",
                                   cdata_->print2screenmutex);
            }
            break;
        }
        }
    }
    MEMORY ps_free (srclocflags);
    if (penalty)
        return penalty;

    // this sets the source flux on the source-plane grid and checks user bounds
    ANALYTICSRC setusersource (data_, cdata_, vars_, source, &penalty,
                               data_->srcbounds, data_->numsrcbounds,
                               data_->usersetsrc, vars_->mps);

    if( penalty )
        return penalty;

    double *sourceloc = 0;
    ANALYTICSRC expandsrc( data_, cdata_, vars_, source, &sourceloc, data_->usersetsrc );


    // this lenses the source that has been placed on the source grid
    COMMON lensgalaxy (data_, cdata_, vars_, vars_->lensingoperator,
                       vars_->mps, vars_->lensedmps);

    // check for uv-plane data
    // gug: need to do error checking for uv transmode
    bool is_uv = data_->is_uvdata;
    vars_->residual4stat = is_uv ? vars_->uv_residual : vars_->residual;
    vars_->data->minus(vars_->lensedmps, vars_->residual, 1, 1);

    STATISTIC actualcomputeed (data_, cdata_, vars_);

    MEMORY ps_free (sourceloc);
    MEMORY ps_free (source);

    return vars_->ed;
}

void pixsrc_analytic_source::setusersource( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *source, double *penalty, double **bounds, PS_SIT numbounds, double *src0, VECTOR *fluxhere )
{
    double *sourceloc = 0;
    ANALYTICSRC expandsrcwithflags( data_, cdata_, vars_, source, &sourceloc, src0 );

    // test bounds
    if( bounds && penalty )
    {
        double sum;
        double basepenalty = 1.0e100;
        PS_SIT totmt = sourceloc[0]*8;

        for( PS_SIT n=0; n<numbounds; ++n )
        {
            sum = 0;
            for( PS_SIT j=0; j<totmt; ++j )
            {
                if( bounds[n][j] )
                {
                    sum += bounds[n][j] * sourceloc[1+(j/8)*17+1+j%8];
                }
            }

            if( sum<bounds[n][totmt  ] )
                *penalty += basepenalty * ( 1.0 +
                                            (sum-bounds[n][totmt  ])*
                                            (sum-bounds[n][totmt  ]) );
            if( sum>bounds[n][totmt+1] )
                *penalty += basepenalty * ( 1.0 +
                                            (sum-bounds[n][totmt+1])*
                                            (sum-bounds[n][totmt+1]) );
        }

        if( *penalty )
        {
            MEMORY ps_free( sourceloc );
            return;
        }
    }

    ANALYTICSRC srcflux( data_, cdata_, vars_, sourceloc, fluxhere, 1,
                         vars_->triout->pointlist);

    MEMORY ps_free( sourceloc );
}

void pixsrc_analytic_source::expandsrc( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *src, double **srcreturn, double *src0 )
{
    PS_SIT length = src0[0]*9+1;
    /*
      if( *srcreturn )
      MEMORY ps_free( *srcreturn );
    */
    if( !*srcreturn )
        MEMORY ps_malloc( &(*srcreturn), length );

    // expanding the abridged source parameters that amoeba uses for computeed
    PS_SIT trackerloc = 0;
    (*srcreturn)[0] = src0[0];
    for( PS_SIT j=0; j<src0[0]; ++j )
    {
        // copy user initial params
        std::copy( src0 + 1+j*17, src0 + 1+j*17+9,
                   *srcreturn + 1+j*9 );

        // copy newly varied params
        for( PS_SIT g=0; g<8; ++g )
        {
            if( src0[10+j*17+g] )
            {
                (*srcreturn)[1+j*9+1+g] = src[trackerloc++];
            }
        }
    }
}

void pixsrc_analytic_source::expandsrcwithflags( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *src, double **srcreturn, double *src0 )
{
    PS_SIT length = src0[0]*17+1;
    /*
      if( *srcreturn )
      MEMORY ps_free( *srcreturn );
    */
    if( !*srcreturn )
        MEMORY ps_malloc( &(*srcreturn), length );

    // expanding the abridged source parameters that amoeba uses
    std::copy( src0, src0+length, *srcreturn );

    // putting in newly varied parameters
    PS_SIT index = 0;
    for( PS_SIT s=0; s<src0[0]; ++s )
        for( PS_SIT p=s*17+1+9; p<s*17+1+17; ++p )
            if( src0[p] )
                (*srcreturn)[p-8] = src[index++];
}

// multithreading code .. decided not worth it
/*
  struct srcfluxstruct
  {
  inputdata *data_;
  PS_SIT srcindex;
  double xpos, ypos, val;
  double a, b, c, d;
  };

  void* pixsrc_analytic_source::srcfluxthread( void* )
  {
  srcfluxstruct sfs = (srcfluxstruct*)void;

  sfs->val = ANALYTICSRC sersicflux( sfs->a, xpos, ypos, sfs->b, sfs->c, sfs->d );
  }
*/

void pixsrc_analytic_source::srcflux( inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                      double *sourceloc_, VECTOR *source, PS_SIT sourcetype, double *pos)
{
    // this function is for setting flux onto source grid
    // sourceloc is the source
    // source is where flux will be stored
    // sourcetype indicates whether sourceloc contains vary flags or not
    // pos holds positions where flux is to be calculated (can be NULL if subsampling)

    // GUG: write multithreaded code here. get all subsampled positions (or not subsamples) in one array and split it up across cpus/gpus.

    source->zeromeout ();

    double x=0,y=0,val,qaxis;
    double cos1,sin1,xrot,yrot;

    PS_SIT numsrcs = sourceloc_[0];
    double *sourceloc;
    MEMORY ps_malloc( &sourceloc, 1 + numsrcs*9 + sourcetype*numsrcs*8  );
    std::copy(       sourceloc_,
                     sourceloc_ + 1 + numsrcs*9 + sourcetype*numsrcs*8,
                     sourceloc                                          );

    double *srcdata;

    // convert sources to pixel coordinates
    pthread_mutex_lock( cdata_->wcsmutex );
    for( PS_SIT src=0; src<numsrcs; ++src )
    {
        // current source
        srcdata = &(sourceloc[1+src*9+sourcetype*src*8]);

        switch ((PS_SIT)srcdata[0])
        {
            // if sersic
        case 0:
        {
            // converting center position from arcseconds to pixels
            double poss[2];
            HEADER getwcssfromwcsl( srcdata[2], srcdata[3], data_->pra, data_->pdec,
                                    data_->r1, data_->r2, poss                       );
            HEADER getimgpixcoord( data_->wcs, data_->imgy, cdata_, poss[0], poss[1],
                                   &srcdata[2], &srcdata[3]  );

            // convert scale radius from arcseconds to pixels
            srcdata[6] *= data_->arc2pix;
            break;
        }
        default:
        {
            // if vector image
            if (srcdata[0]>1000)
            {
                PS_SIT vecindex = srcdata[0]-1001;
                // convert shift from arcsec to pixel with a little hack

                // first find image position of first pixel (unshifted)
                double pos1[2] = {data_->vector_src_pos[vecindex][0],
                                  data_->vector_src_pos[vecindex][1]};
                double pos_orig[2], pos_shift[2], pos2[2];
                HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, pos1[0], pos1[1],
                                       &pos_orig[0], &pos_orig[1]);
                // now apply srcdata shift at position of first pixel and get its image position
                HEADER getwcssfromwcsl (srcdata[2], srcdata[3], pos1[0], pos1[1], 0, 0, pos2);
                HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, pos2[0], pos2[1],
                                       &pos_shift[0], &pos_shift[1]);
                // subtract image positions to get shifts in x and y (pixel coordinates)
                srcdata[2] = pos_shift[0] - pos_orig[0];
                srcdata[3] = pos_shift[1] - pos_orig[1];
            }

            break;
        }
        }
    }
    pthread_mutex_unlock( cdata_->wcsmutex );

    // vector sync occurs here

    for( PS_SIT src=0; src<numsrcs; ++src )
    {
        // this source
        srcdata = &(sourceloc[1+src*9+sourcetype*src*8]);

        switch ((PS_SIT)srcdata[0])
        {
            // if sersic
        case 0:
        {
            // if calculating lensed surface brightness exactly, skip putting
            // source on source grid (sersic only)
            if (data_->srcexact)
                break;

            // convert angle to one measured from +x axis
	    // handle negative ellipticity
	    double angle;
	    if (srcdata[4]<0)
	    {
		qaxis = 1.0+srcdata[4];
		angle = (srcdata[5]+90.0)*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2;
	    }
	    else
	    {
		qaxis = 1.0-srcdata[4];
		angle = srcdata[5]*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2;
	    }
            cos1 = std::cos (angle);
            sin1 = std::sin (angle);

            for(PS_SIT s=0; s<source->get_size(); ++s)
            {
                x = pos[s*2  ] - srcdata[2];
                y = pos[s*2+1] - srcdata[3];
                y *= -1;
                xrot = x*cos1 + y*sin1;
                yrot = -x*sin1 + y*cos1;

                val = ANALYTICSRC sersicflux (srcdata[1], xrot, yrot,
                                              qaxis, srcdata[6], srcdata[8]);
                source->set( s, source->get(s) + val );
            }
            break;
        }
        default:
        {
            // if vector image
            if (srcdata[0]>1000)
            {
                PS_SIT vecindex = srcdata[0]-1001;
                PS_SIT numpoints = data_->vector_src_numpix[vecindex];

                struct triangulateio *trioutloc = data_->vector_src_triout[vecindex];
                double *vertices_orig = trioutloc->pointlist;
                double *vertices_transformed;
                MEMORY ps_malloc( &(vertices_transformed), numpoints*2 );

                // rotate image coordinates to align with WCS,
                // apply linear transformation, and then rotate back
                double tmpmat[4];
                double cosr = std::cos (-data_->rollangle);
                double sinr = std::sin (-data_->rollangle);
                double roll1[4] = {cosr, -sinr,  sinr, cosr};
                double roll2[4] = {cosr,  sinr, -sinr, cosr};
                // extra minus signs are for x->R.A. refelection and y-axis relfection,
                // because pixsrc measures y increasing downwards
                for (PS_SIT m=0; m<4; ++m)
                {
                    roll1[m] *= -1;
                    roll2[m] *= -1;
                }
                double *rot1 = roll1, *rot2 = roll2;
                tmpmat[0]  = srcdata[4]*rot1[0] + srcdata[5]*rot1[2];
                tmpmat[1]  = srcdata[4]*rot1[1] + srcdata[5]*rot1[3];
                tmpmat[2]  = srcdata[6]*rot1[0] + srcdata[7]*rot1[2];
                tmpmat[3]  = srcdata[6]*rot1[1] + srcdata[7]*rot1[3];
                srcdata[4] = rot2[0]*tmpmat[0]  + rot2[1]*tmpmat[2];
                srcdata[5] = rot2[0]*tmpmat[1]  + rot2[1]*tmpmat[3];
                srcdata[6] = rot2[2]*tmpmat[0]  + rot2[3]*tmpmat[2];
                srcdata[7] = rot2[2]*tmpmat[1]  + rot2[3]*tmpmat[3];

                // rotate, stretch, and/or skew vector source pixel positions
                // and shift image back to original position on the sky
                for (PS_SIT s=0; s<numpoints; ++s)
                {
                    vertices_transformed[s*2]   =
                        srcdata[4]*vertices_orig[s*2] + srcdata[5]*vertices_orig[s*2+1] +
                        data_->rotation_axis[vecindex*2  ];
                    vertices_transformed[s*2+1] =
                        srcdata[6]*vertices_orig[s*2] + srcdata[7]*vertices_orig[s*2+1] +
                        data_->rotation_axis[vecindex*2+1];
                }
                // temporarily move transformed pixels into struct
                // original pixel positions are stored in vertices_orig
                trioutloc->pointlist = vertices_transformed;

                PS_SIT triseed = 0;
                for(PS_SIT s=0; s<source->get_size(); ++s)
                {
                    x = pos[s*2  ] - srcdata[2];
                    y = pos[s*2+1] - srcdata[3];

                    val = ANALYTICSRC vectorflux (srcdata[1], x, y, trioutloc,
                                                  data_->vector_src_flux[vecindex], &triseed);
                    source->set( s, source->get(s) + val );
                }

                // restore original pixel positions to triangle struct and destroy transformed pixels
                trioutloc->pointlist = vertices_orig;
                MEMORY ps_free (vertices_transformed);
            }

            break;
        }
        }
    }

    MEMORY ps_free( sourceloc );

    ANALYTICSRC imgflux (data_, cdata_, vars_, sourceloc_, sourcetype, vars_->newloc,
                         vars_->r4r, data_->subsampling, vars_->newloc_sslo);
}

void pixsrc_analytic_source::imgflux (inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                      double *sourceloc_, PS_SIT sourcetype,
                                      double *pos, PS_SIT *posflags, PS_SIT subsample, double **posss )
{
    // this flux is for setting flux in image plane (exact, vs lensing through matrix)
    // sourceloc is the source
    // source is where flux will be stored
    // sourcetype indicates whether sourceloc contains vary flags or not
    // pos holds positions where flux is to be calculated (can be NULL if subsampling)
    // posflags indicate whether to skip a certain position in source (NULL means use all)
    // subsample indicates whether to integrate flux over region or get flux at center
    // posss holds (source plane) positions of subsampled pixels (can be NULL if not subsampling)

    // GUG: write multithreaded code here. get all subsampled positions (or not subsamples) in one array and split it up across cpus/gpus.

    if (!data_->srcexact)
        return;

    if (!vars_->srcexactvec)
    {
        MEMORY ps_malloc (&vars_->srcexactvecptr, 1);
        vars_->srcexactvec = new (vars_->srcexactvecptr) VECTOR (cdata_, data_, vars_->lonr);
    }
    else
    {
        vars_->srcexactvec->zeromeout ();
    }
    VECTOR *source = vars_->srcexactvec;

    double ss2 = subsample*subsample;
    double x=0,y=0,xx,yy,val,qaxis;
    double cos1,sin1,xrot,yrot;

    PS_SIT numsrcs = sourceloc_[0];
    double *sourceloc;
    MEMORY ps_malloc( &sourceloc, 1 + numsrcs*9 + sourcetype*numsrcs*8  );
    std::copy(       sourceloc_,
                     sourceloc_ + 1 + numsrcs*9 + sourcetype*numsrcs*8,
                     sourceloc                                          );

    double *srcdata;

    // convert sources to pixel coordinates
    pthread_mutex_lock( cdata_->wcsmutex );
    for( PS_SIT src=0; src<numsrcs; ++src )
    {
        // current source
        srcdata = &(sourceloc[1+src*9+sourcetype*src*8]);

        switch ((PS_SIT)srcdata[0])
        {
            // if sersic
        case 0:
        {
            // converting center position from arcseconds to pixels
            double poss[2];
            HEADER getwcssfromwcsl( srcdata[2], srcdata[3], data_->pra, data_->pdec,
                                    data_->r1, data_->r2, poss                       );
            HEADER getimgpixcoord( data_->wcs, data_->imgy, cdata_, poss[0], poss[1],
                                   &srcdata[2], &srcdata[3]  );

            // convert scale radius from arcseconds to pixels
            srcdata[6] *= data_->arc2pix;
            break;
        }
        default: break;
        }
    }
    pthread_mutex_unlock( cdata_->wcsmutex );

    PS_SIT postracker;

    // vector sync occurs here

    for( PS_SIT src=0; src<numsrcs; ++src )
    {
        // this source
        srcdata = &(sourceloc[1+src*9+sourcetype*src*8]);

        switch ((PS_SIT)srcdata[0])
        {
            // if sersic
        case 0:
        {
            // convert angle to one measure from +x axis
	    if(srcdata[4]<0)
	    {
            qaxis = 1.0+srcdata[4];
            cos1 = std::cos((srcdata[5]+90.0)*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2);
	    sin1 = std::sin((srcdata[5]+90.0)*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2);
	    }
	    else
	    {
            qaxis = 1.0-srcdata[4];
            cos1 = std::cos(srcdata[5]*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2);
            sin1 = std::sin(srcdata[5]*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2);
	    }

            postracker = 0;
            for(PS_SIT s=0; s<source->get_size(); ++s)
            {
                if( !posflags )
                {
                    x = pos[s*2  ];
                    y = pos[s*2+1];
                }
                else
                {
                    // find next location to get flux at
                    while( posflags[postracker] == -1 )
                        ++postracker;

                    if( subsample == 1 )
                    {
                        x = pos[postracker*2  ];
                        y = pos[postracker*2+1];
                    }
                }

                if( subsample == 1 )
                {
                    x -= srcdata[2];
                    y -= srcdata[3];
                    y *= -1;

                    xrot = x*cos1 - y*sin1;
                    yrot = x*sin1 + y*cos1;

                    val = ANALYTICSRC sersicflux (srcdata[1], xrot, yrot,
                                                  qaxis, srcdata[6], srcdata[8]);
                    source->set( s, source->get(s) + val );
                }
                else
                {
                    val = 0;

                    // was gonna use this for multithreading
                    //if( cdata_->numthreads == 1 )
                    {
                        for( PS_SIT ss=0; ss<ss2; ++ss )
                        {
                            xx = posss[postracker][ss*2  ] - srcdata[2];
                            yy = posss[postracker][ss*2+1] - srcdata[3];

                            xrot = xx*cos1 - yy*sin1;
                            yrot = xx*sin1 + yy*cos1;

                            val += ANALYTICSRC sersicflux (srcdata[1], xrot, yrot,
                                                           qaxis, srcdata[6], srcdata[8]);
                        }
                    }
                    // was gonna use this for multithreading
                    /*
                      else
                      {

                      }
                    */

                    source->set( s, source->get(s) + val );
                }

                if( posflags )
                    ++postracker;
            }
            break;
        }
        default: break;
        }
    }

    if (subsample>1)
        source->set_scalar( 1.0 / ss2 );

    MEMORY ps_free( sourceloc );
}

double pixsrc_analytic_source::sersicflux( double I0, double xrot, double yrot,
                                           double q, double rs, double k )
{
    return I0
        * std::exp(
            -std::pow(
                ( xrot*xrot+yrot*yrot/(q*q) ) / ( rs*rs ) ,
                1.0/(2*k)
                )
            );
}

double pixsrc_analytic_source::vectorflux (double I0, double x, double y,
                                           struct triangulateio *trioutloc, double *flux, PS_SIT *triseed)
{
    // find where pixel x,y lies in the triangulation and return flux
    PS_SIT index = GEOM search_triangle( trioutloc, triseed, x, y);

    if (index==-1)
        return 0;

    // find vertices of triangle
    double pos[2][3], vals[3], result;
    for (PS_SIT f=0; f<3; ++f)
    {
        pos[0][f] = trioutloc->pointlist[trioutloc->trianglelist[index*3+f]*2  ];
        pos[1][f] = trioutloc->pointlist[trioutloc->trianglelist[index*3+f]*2+1];
        vals[f]   = flux[trioutloc->trianglelist[index*3+f]];
    }
    OPERA planarvalueinterpolation (pos, vals, x, y, &result);

    return result*I0;
}
