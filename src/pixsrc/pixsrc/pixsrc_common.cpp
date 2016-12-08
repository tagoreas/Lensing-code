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


#include "pixsrc_cuda.hpp"

#include "pixsrc_common.hpp"
#include "pixsrc_common_adaptive.hpp"
#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_external.hpp"
#include "pixsrc_tps.hpp"
#include "pixsrc_init.hpp"
#include "pixsrc_analytic_source.hpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_printer_templates.cpp"
#include <cmath>
#include <cstring>
#include <algorithm>
#include <exception>

bool trieqon = 0;

void pixsrc_common::computesizepenaltybody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror || !data_->penaltyquery[0] )
    {
        return;
    }

    double coords[6];
    double norm = 0;
    double sizepenalty = 0;
    for(PS_SIT t=0; t<vars_->triout->numberoftriangles; t++)
    {
        double flux = 0;
        for( PS_SIT v=0; v<3; ++v )
        {
            flux += std::abs(vars_->mps->get(vars_->triout->trianglelist[t*3+v]));
        }
        flux /= 3.0;
        norm += flux;

        if( flux )
        {
            for(PS_SIT j=0; j<3; j++)
            {
                coords[j*2  ] = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2  ];
                coords[j*2+1] = vars_->triout->pointlist[vars_->triout->trianglelist[t*3+j]*2+1];
            }
        }
        sizepenalty += GEOM areatriangle(coords)*flux;
    }

    sizepenalty /= norm;

    penaltypenalizer( data_, cdata_, vars_, 0, 0, sizepenalty, "brightness-weighted size" );

    if( data_->verbose != 1 && !vars_->terminatelensing )
        PRINTER print2screen(data_->print2screenname,
                             "brightness-weighted size penalty = " +
                             OPERA tostring(vars_->penalties[0]),
                             cdata_->print2screenmutex);
}

void* pixsrc_common::computesizepenalty(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[5] = 1;

    COMMON computesizepenaltybody( data_, cdata_, vars_ );

    vars_->pthreadstracker[5] = 2;

    return NULL;
}
void pixsrc_common::computemagpenalty(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if( vars_-> fatalerror || !data_->penaltyquery[1] )
    {
        return;
    }

    penaltypenalizer( data_, cdata_, vars_, 0, 1, vars_->magger, "magnification" );

    if( data_->verbose!=1 && !vars_->terminatelensing )
    {
        if( data_->penaltymatrix[0][1][0]==2 )
        {
            PRINTER print2screen(data_->print2screenname,
                                 "model passed magnification test",
                                 cdata_->print2screenmutex);
        }
        else
        {
            PRINTER print2screen(data_->print2screenname,
                                 "magnification penalty = " +
                                 OPERA tostring(vars_->penalties[1]),
                                 cdata_->print2screenmutex);
        }
    }
}

void pixsrc_common::computepenaltybody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    // i've incorporated this into computeed
}

void* pixsrc_common::computepenalty(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    commoninputdata *cdata_ = dav->cdata;
    inputdata *data_ = dav->data;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[6] = 1;

    COMMON computepenaltybody( data_, cdata_, vars_ );

    vars_->pthreadstracker[6] = 2;

    return NULL;
}

void pixsrc_common::startlensing(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    MEMORY initialize  (data_, cdata_, vars_);
    MEMORY ps_malloc( &(vars_->dav), 1 );
    vars_->dav->data  = data_;
    vars_->dav->cdata = cdata_;
    vars_->dav->vars  = vars_;

    vars_->lambda1 = 0;
    vars_->xoffset = vars_->yoffset = 0.0;
    vars_->srcopt = 0;
    vars_->shapeletopt = 0;
    vars_->fatalerror = vars_->nonfatalerror = 0;
    vars_->terminatelensing=0;
    std::fill(vars_->penalties,vars_->penalties+data_->extlengths[5],0);
    vars_->endedearly=0;
    vars_->modelrejected=0;

    std::fill(vars_->pthreadstracker,vars_->pthreadstracker+vars_->numpthreads,0);
    for(PS_SIT j=0; j<vars_->numpthreads; ++j)
        pthread_mutex_init(&(vars_->pthreadslock[j]),NULL);

    MEMORY ps_malloc( &(vars_->newloc ), data_->ndp*2 );
    MEMORY ps_malloc( &(vars_->r4r    ), data_->ndp   );
    MEMORY ps_malloc( &(vars_->r4rback), data_->ndp   );
    std::fill( vars_->r4r, vars_->r4r+data_->ndp, -1  );

    vars_->variance = data_->myvariance0;
    vars_->datavariance = data_->myvariance0;

    vars_->need2retriangulate = 0;

    if( data_->gridtype == 2 )
    {
        vars_->grid        = new vector< vector< vector< vector<bool> > > >();
        vars_->gridpointer = new vector< vector< vector< vector<PS_SIT > > > >();
    }

    if(data_->debug)
        PRINTER print <PS_SIT> ( vars_->tracker, cdata_->basename, data_->name,
                                 true,"oldloc.dat",data_->oldloc,data_->ndp*2, 0,
                                 data_->precision, NULL);

    if(!data_->useall)
    {
        // apply masks
        for(PS_SIT r=0; r<data_->ndp; r++)
            if(data_->imagemasks[r]==2)
                vars_->r4r[r]=0;

        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "image masks applied",
                                 cdata_->print2screenmutex);
    }
    else
    {
        std::fill( vars_->r4r, vars_->r4r+data_->ndp, 0 );
    }

    if(data_->gridtype==2 || data_->use_shapelets)
    {
        MEMORY ps_malloc( &(vars_->magnification), data_->ndp );
    }

    if(data_->verbose==3)
        PRINTER print2screen(data_->print2screenname,
                             "ray tracing data pixels",
                             cdata_->print2screenmutex);

    double magxx, magxy, magyx, magyy, pot;

    // only ray-trace for image pixels that are needed for testformmimages
    if( data_->penaltyquery[2] ||
        data_->penaltyquery[3] || data_->penaltyquery[4] )
    {
        pthread_mutex_lock (cdata_->potdefmagmutex);
        pthread_mutex_lock (cdata_->wcsmutex);
        for (PS_SIT r = 0; r < data_->ndp; r++)
        {
            bool doit = 0;
            for (PS_SIT src=0; src<data_->extlengths[11]; ++src)
                if (data_->mmimages[src][r])
                {
                    doit = 1;
                    break;
                }
            if (doit)
            {

                COMMON raytrace (data_, cdata_,
                                 data_->oldloc[r*2], data_->oldloc[r*2+1],
                                 &(vars_->newloc[r*2]), &(vars_->newloc[r*2+1]),
                                 &magxx, &magxy, &magyx, &magyy, &pot, r);

                if (data_->gridtype==2 || data_->use_shapelets)
                    vars_->magnification[r]=std::fabs(1/((1-magxx)*(1-magyy)-magxy*magyx));
            }
        }
        // ray-trace border of mmimage masks
        for (PS_SIT src=0; src<data_->extlengths[11]; ++src)
            for(PS_SIT img=0; img<data_->num_mmimages[src]; ++img)
                for(PS_SIT vert=0; vert<(PS_SIT)data_->mmborder[src][img][0]; ++vert)
                    COMMON raytrace (data_, cdata_,
                                     data_->mmborder[src][img][vert*2+0+1],
                                     data_->mmborder[src][img][vert*2+1+1],
                                     &(data_->mmborder_defl[src][img][vert*2+0+1]),
                                     &(data_->mmborder_defl[src][img][vert*2+1+1]),
                                     &magxx, &magxy, &magyx, &magyy, &pot, -1);
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );
    }

    COMMON testformismatchedimages(data_, cdata_, vars_);

    if (vars_->fatalerror || vars_->terminatelensing ||
        vars_->endedearly || !data_->doreconstruction   )
    {
        return;
    }

    // ray-trace rest of data pixels
    if( data_->penaltyquery[2] ||
        data_->penaltyquery[3] || data_->penaltyquery[4] )
    {
        pthread_mutex_lock( cdata_->potdefmagmutex );
        pthread_mutex_lock( cdata_->wcsmutex       );
        for(PS_SIT r = 0; r < data_->ndp; r++)
        {
            bool doit = 0;
            for( PS_SIT src=0; src<data_->extlengths[11]; ++src )
                if( data_->mmimages[src][r] )
                {
                    doit = 1;
                    break;
                }
            if (!doit &&
                (data_->magparams ||
                 data_->findsisterimages ||
                 !OPERA is_exactly_p_inf(data_->srcinputcircle[0][2]) ||
                 data_->imagemasks[r]==2)
                )
            {

                COMMON raytrace (data_, cdata_,
                                 data_->oldloc[r*2], data_->oldloc[r*2+1],
                                 &(vars_->newloc[r*2]), &(vars_->newloc[r*2+1]),
                                 &magxx, &magxy, &magyx, &magyy, &pot, r);

                if(data_->gridtype==2 || data_->use_shapelets)
                    vars_->magnification[r]=std::fabs(1/((1-magxx)*(1-magyy)-magxy*magyx));
            }
        }
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );
    }
    // if no mmimages, then ray-trace all data pixels
    else
    {
        pthread_mutex_lock( cdata_->potdefmagmutex );
        pthread_mutex_lock( cdata_->wcsmutex       );
        for(PS_SIT r = 0; r < data_->ndp; r++)
        {
            if (data_->magparams ||
                data_->findsisterimages ||
                !OPERA is_exactly_p_inf(data_->srcinputcircle[0][2]) ||
                data_->imagemasks[r]==2)
            {

                COMMON raytrace (data_, cdata_,
                                 data_->oldloc[r*2], data_->oldloc[r*2+1],
                                 &(vars_->newloc[r*2]), &(vars_->newloc[r*2+1]),
                                 &magxx, &magxy, &magyx, &magyy, &pot, r);

                if(data_->gridtype==2 || data_->use_shapelets)
                    vars_->magnification[r]=std::fabs(1/((1-magxx)*(1-magyy)-magxy*magyx));
            }
        }
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );
    }

    if(data_->debug)
    {
        PRINTER print <double> ( vars_->tracker, cdata_->basename,
                                 data_->name, true,"newloc.dat",
                                 vars_->newloc,data_->ndp*2, 0, data_->precision, NULL);
    }

    if(data_->verbose==3)
    {
        PRINTER print2screen(data_->print2screenname,
                             "done ray tracing data pixels",
                             cdata_->print2screenmutex);
    }

    OPERA assign_p_infinity( &vars_->minx );
    OPERA assign_p_infinity( &vars_->miny );
    OPERA assign_n_infinity( &vars_->maxx );
    OPERA assign_n_infinity( &vars_->maxy );

    for(PS_SIT m = 0; m<data_->ndp; m++)
        if (data_->magparams ||
            data_->findsisterimages ||
            !OPERA is_exactly_p_inf(data_->srcinputcircle[0][2]) ||
            data_->imagemasks[m]==2)
            //if( 1 || vars_->r4r[m]!=-1 )
        {
            if(vars_->newloc[m*2  ]>vars_->maxx)
                vars_->maxx=vars_->newloc[m*2  ];
            if(vars_->newloc[m*2+1]>vars_->maxy)
                vars_->maxy=vars_->newloc[m*2+1];
            if(vars_->newloc[m*2  ]<vars_->minx)
                vars_->minx=vars_->newloc[m*2  ];
            if(vars_->newloc[m*2+1]<vars_->miny)
                vars_->miny=vars_->newloc[m*2+1];
        }

    vars_->maxx += CONSTANT smallnumber;
    vars_->maxy += CONSTANT smallnumber;
    vars_->minx -= CONSTANT smallnumber;
    vars_->miny -= CONSTANT smallnumber;

    vars_->rangex=vars_->maxx-vars_->minx;
    vars_->rangey=vars_->maxy-vars_->miny;

    vars_->ctrx = (vars_->maxx+vars_->minx)/2;
    vars_->ctry = (vars_->maxy+vars_->miny)/2;

    COMMON subsample_ie( vars_->dav );
}

void pixsrc_common::raytrace (inputdata *data_, commoninputdata *cdata_,
                              double x, double y, double *xx, double *yy,
                              double *magxx, double *magxy,
                              double *magyx, double *magyy, double *pot, PS_SIT r)
{
    // **************************
    // ******** WARNING *********
    // **************************
    // must put mutex locks on POTDEFMAG and WCS before
    // calling this function

    double pos[2];
    double defx, defy;
    double ra0,dec0;
    double npl_alpha[2];

    if( !strcmp( cdata_->coordsysorig, "PIXEL" ) )
    {
        HEADER getlfroms( x, y, data_->px, data_->py, data_->r1, data_->r2, pos );
    }
    else if (r==-1)
    {
        HEADER getimgwcscoord (data_->wcs, data_->imgy, x, y, &ra0, &dec0);
        HEADER getwcslfromwcss (ra0, dec0, data_->pra, data_->pdec,
                                data_->r1, data_->r2, pos);
    }
    else
    {
        ra0    = data_->oldloc_wcs[r*2];
        dec0   = data_->oldloc_wcs[r*2+1];
        pos[0] = data_->oldloc_arc[r*2];
        pos[1] = data_->oldloc_arc[r*2+1];
    }

    // rotate and translate coordinate
    double oldpos[2] = {pos[0]-data_->rottrans[2],pos[1]-data_->rottrans[3]};
    pos[0] =  oldpos[0]*data_->rottrans[0] - oldpos[1]*data_->rottrans[1];
    pos[1] =  oldpos[0]*data_->rottrans[1] + oldpos[1]*data_->rottrans[0];
    pos[0] += data_->rottrans[4];
    pos[1] += data_->rottrans[5];

#ifdef PS_HAVE_GRAVLENS
    EXTERNAL ps_potdefmag (pos[0], pos[1], pot, &defx, &defy,
                           magxx, magyy, magxy, magyx);
#endif
#ifdef PS_HAVE_TRIAXIAL
    EXTERNAL ps_potdefmag (pos[0], pos[1], pot, &defx, &defy,
                           magxx, magyy, magxy, magyx, cdata_->tlmparms, cdata_->tlmtime, cdata_->tlmenvirogals);
#endif

    if( !strcmp( cdata_->coordsysorig, "PIXEL" ) )
    {
        *xx = x - defx;
        *yy = y + defy;
    }
    else
    {
        // if using non-parametric lens potential perturbations
        if( cdata_->npl )
        {
            // Note: I could calculate potential here as well,
            // but it's not needed right now.

            NONPARAMLENS get_alpha (data_, cdata_, pos[0], pos[1],
                                    &npl_alpha[0], &npl_alpha[1]);
            defx += npl_alpha[0];
            defy += npl_alpha[1];
        }

        HEADER getwcssfromwcsl( -defx, -defy, ra0, dec0, 0, 0, pos );
        HEADER getimgpixcoord( data_->wcs, data_->imgy, cdata_, pos[0], pos[1], xx, yy );

        *pot *= data_->arc2pix*data_->arc2pix;
    }
}

void pixsrc_common::raytrace_s2i (inputdata *data_, commoninputdata *cdata_,
                                  double u, double v, PS_SIT *numimg, double ***imgarr)
{
    // **************************
    // ******** WARNING *********
    // **************************
    // must put locks on POTDEFMAG and WCS before
    // calling this function

    if( cdata_->npl )
    {
        PRINTER printerror (data_->print2screenname,
                            "cannot compute images for source position "
                            "with non-parametric lens potential "
                            "perturbations.. yet",
                            cdata_->print2screenmutex);
    }

    double pos[2];
    double ra0,dec0;

    if( !strcmp( cdata_->coordsysorig, "PIXEL" ) )
    {
        HEADER getlfroms( u, v, data_->px, data_->py, data_->r1, data_->r2, pos );
    }
    else
    {
        HEADER getimgwcscoord( data_->wcs, data_->imgy, u, v, &ra0, &dec0 );
        HEADER getwcslfromwcss( ra0, dec0, data_->pra, data_->pdec,
                                data_->r1, data_->r2, pos );
    }

    EXTERNAL ps_find_img (pos[0], pos[1], numimg, imgarr);

    for (PS_SIT n=1; n<=*numimg; ++n)
    {
        if( !strcmp( cdata_->coordsysorig, "PIXEL" ) )
        {}
        else
        {
            HEADER getwcssfromwcsl ((*imgarr)[n][1], (*imgarr)[n][2],
                                    data_->pra, data_->pdec,
                                    data_->r1, data_->r2, pos);
            HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, pos[0], pos[1],
                                   &((*imgarr)[n][1]), &((*imgarr)[n][2]));
        }
    }
}

void pixsrc_common::subsample_ie( void *args )
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    // subsample pixel that have been masked as data
    if( data_->interperr > 1 )
    {

        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "oversampling pixels for interpolation errors",
                                 cdata_->print2screenmutex);

        double x, y, xxx, yyy;
        PS_SIT ss, ss2 = data_->interperr*data_->interperr;
        double one_2ss = 1.0 / ( 2 * data_->interperr );

        double magxx, magxy, magyx, magyy, pot;

        MEMORY ps_malloc( &vars_->newloc_ssas, data_->ndp );

        pthread_mutex_lock( cdata_->potdefmagmutex );
        pthread_mutex_lock( cdata_->wcsmutex       );
        for( PS_SIT r=0; r<data_->ndp; ++r )
        {
            if( data_->imagemasks[r] == 2 )
            {
                x = data_->oldloc[r*2  ];
                y = data_->oldloc[r*2+1];

                MEMORY ps_malloc( &(vars_->newloc_ssas[r]), ss2*2 );

                for( PS_SIT xx=0; xx<data_->interperr; ++xx )
                {
                    xxx = one_2ss * ( 2*data_->interperr*x - data_->interperr + 1 + 2*xx );

                    for( PS_SIT yy=0; yy<data_->interperr; ++yy )
                    {
                        yyy = one_2ss * ( 2*data_->interperr*y + data_->interperr - 1 - 2*yy );
                        ss = xx*data_->interperr + yy;

                        COMMON raytrace (data_, cdata_, xxx, yyy,
                                         &(vars_->newloc_ssas[r][ss*2  ]),
                                         &(vars_->newloc_ssas[r][ss*2+1]),
                                         &magxx, &magxy, &magyx, &magyy, &pot, -1);
                    }
                }
            }
            else
            {
                vars_->newloc_ssas[r] = 0;
            }
        }
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );

        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "done oversampling pixels for interpolation errors",
                                 cdata_->print2screenmutex);

    }
}

void pixsrc_common::resubsample_ie( void *args )
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    // subsample pixels that weren't masked as data but
    // weren't masked as bad pixels either
    // and were found to be sister pixels
    if( data_->interperr > 1 )
    {
        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "re-oversampling sister pixels for interpolation errors",
                                 cdata_->print2screenmutex);

        PS_SIT x, y, xxx, yyy, ss;
        PS_SIT ss2 = data_->interperr*data_->interperr;
        double one_2ss = 1.0 / ( 2 * data_->interperr );

        double magxx, magxy, magyx, magyy, pot;

        pthread_mutex_lock( cdata_->potdefmagmutex );
        pthread_mutex_lock( cdata_->wcsmutex       );
        for( PS_SIT r=0; r<data_->ndp; ++r )
        {
            if( data_->imagemasks[r] == 1 && vars_->r4r[r] != -1 )
            {
                x = data_->oldloc[r*2  ];
                y = data_->oldloc[r*2+1];

                MEMORY ps_malloc( &(vars_->newloc_ssas[r]), ss2*2 );

                for( PS_SIT xx=0; xx<data_->interperr; ++xx )
                {
                    xxx = one_2ss * ( 2*data_->interperr*x - data_->interperr + 1 + 2*xx );

                    for( PS_SIT yy=0; yy<data_->interperr; ++yy )
                    {
                        yyy = one_2ss * ( 2*data_->interperr*y + data_->interperr - 1 - 2*yy );
                        ss = xx*data_->interperr + yy;

                        COMMON raytrace (data_, cdata_, xxx, yyy,
                                         &(vars_->newloc_ssas[r][ss*2  ]),
                                         &(vars_->newloc_ssas[r][ss*2+1]),
                                         &magxx, &magxy, &magyx, &magyy, &pot, -1);
                    }
                }
            }
        }
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );

        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "done re-oversampling sister pixels for interpolation errors",
                                 cdata_->print2screenmutex);
    }
}

void pixsrc_common::penaltypenalizer( inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                      PS_SIT imgind, PS_SIT penaltyindex, double val, string testname )
{
    if (vars_->fatalerror)
    {
        return;
    }

    char fallsinrejectzone = (val<data_->penaltymatrix[imgind][penaltyindex][3] ||
                              val>data_->penaltymatrix[imgind][penaltyindex][4]) ? 1 : 0;

    if( data_->penaltymatrix[imgind][penaltyindex][0]==1 ||
        (data_->penaltymatrix[imgind][penaltyindex][0]==3 && !fallsinrejectzone) )
    {
        vars_->penalties[penaltyindex] += (val-data_->penaltymatrix[imgind][penaltyindex][1])*
            (val-data_->penaltymatrix[imgind][penaltyindex][1]) /
            ( data_->penaltymatrix[imgind][penaltyindex][2] * data_->penaltymatrix[imgind][penaltyindex][2] );
    }
    if( (data_->penaltymatrix[imgind][penaltyindex][0]==2 ||
         data_->penaltymatrix[imgind][penaltyindex][0]==3) && fallsinrejectzone )
    {
        double diff = std::min( std::fabs(val-data_->penaltymatrix[imgind][penaltyindex][3]),
                                std::fabs(val-data_->penaltymatrix[imgind][penaltyindex][4]) );
        vars_->penalties[penaltyindex] += (1+diff*diff)
            /data_->penaltymatrix[imgind][penaltyindex][5]/data_->penaltymatrix[imgind][penaltyindex][5];
        vars_->terminatelensing = 1;
        vars_->modelrejected++;
        if(data_->verbose != 1)
            PRINTER print2screen(data_->print2screenname,
                                 "model rejected by " + testname + " test",
                                 cdata_->print2screenmutex);
    }
}
void pixsrc_common::testformismatchedimages(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{



    // NOTE: add input parameter for controlling whether nested mmimages are given frac of 1 or <1
    // perhaps call it: nestedmmimages



    if (vars_->fatalerror || vars_->terminatelensing ||
        (!data_->penaltyquery[2] && !data_->penaltyquery[3]
         && !data_->penaltyquery[4]) )
    {
        return;
    }

    for(PS_SIT src=0; src<data_->extlengths[11]; ++src)
    {
        PS_SIT makesepconvexhulls = ( data_->penaltymatrix[src][2][0] ||
                                      data_->penaltymatrix[src][3][0] ) ? 1 : 0;

        PS_SIT numimgs = data_->num_mmimages[src];

        double **convexhulls = 0;
        PS_SIT numverts[numimgs];
        double areas[numimgs];

        if(makesepconvexhulls)
        {
            MEMORY ps_malloc( &(convexhulls), numimgs );

            for(PS_SIT img=0; img<numimgs; img++)
            {
                PS_SIT numpoints=0;
                for(PS_SIT r=0; r<data_->ndp; r++)
                    if(data_->mmimages[src][r]==img+1)
                        ++numpoints;

                // also include border pixels
                numpoints += (PS_SIT)data_->mmborder[src][img][0];

                struct triangulateio *triinloc, *trioutloc;
                MEMORY ps_malloc(            &triinloc , 1 );
                MEMORY ps_malloc(            &trioutloc, 1 );
                MEMORY triangulatestructinit( triinloc     );
                MEMORY triangulatestructinit( trioutloc    );

                triinloc->numberofpointattributes = 0;
                triinloc->numberofpoints = numpoints;
                MEMORY ps_malloc( &(triinloc->pointlist), numpoints*2 );

                PS_SIT index=0;
                for(PS_SIT r=0; r<data_->ndp; r++)
                    if(data_->mmimages[src][r]==img+1)
                    {
                        triinloc->pointlist[index*2  ] = vars_->newloc[r*2  ];
                        triinloc->pointlist[index*2+1] = vars_->newloc[r*2+1];
                        index++;
			if (vars_->newloc[r*2  ]<-data_->imgx || vars_->newloc[r*2  ]>2*data_->imgx || 
			    vars_->newloc[r*2+1]<-data_->imgy || vars_->newloc[r*2+1]>2*data_->imgy)
			{
			    if (numpoints>=4)
			{
			    --index;
			    --numpoints;
			    --triinloc->numberofpoints;
			}
			    else
			    {
				triinloc->pointlist[(index-1)*2  ] = OPERA randomgaussian(cdata_->ps_gsl_ran_r);
				triinloc->pointlist[(index-1)*2+1] = OPERA randomgaussian(cdata_->ps_gsl_ran_r);
			    }
			}
                    }

                // also add border pixels
                for(PS_SIT vert=0; vert<(PS_SIT)data_->mmborder[src][img][0]; ++vert)
                {
                    triinloc->pointlist[index*2  ] = data_->mmborder_defl[src][img][vert*2+0+1];
                    triinloc->pointlist[index*2+1] = data_->mmborder_defl[src][img][vert*2+1+1];
                    index++;
			if (data_->mmborder_defl[src][img][vert*2+0+1]<-data_->imgx || 
			    data_->mmborder_defl[src][img][vert*2+0+1]>2*data_->imgx || 
			    data_->mmborder_defl[src][img][vert*2+1+1]<-data_->imgy || 
			    data_->mmborder_defl[src][img][vert*2+1+1]>2*data_->imgy)
			{
			    if (numpoints>=4)
			{
			    --index;
			    --numpoints;
			    --triinloc->numberofpoints;
			}
			    else
			    {
				triinloc->pointlist[(index-1)*2  ] = OPERA randomgaussian(cdata_->ps_gsl_ran_r);
				triinloc->pointlist[(index-1)*2+1] = OPERA randomgaussian(cdata_->ps_gsl_ran_r);
			    }
			}
                }

                GEOM getconvexhull( triinloc,trioutloc,&(numverts[img]),
                                    &(convexhulls[img]),NULL );

                MEMORY triangulatestructdestruct( triinloc  );
                MEMORY triangulatestructdestruct( trioutloc );
                MEMORY ps_free(                   triinloc  );
                MEMORY ps_free(                   trioutloc );

                // compute areas
                areas[img] = GEOM areapoly(convexhulls[img],numverts[img]);
            }
            if( data_->debug ||
                data_->printvec==1 || data_->printvec==3 )
            {
                for(PS_SIT img=0; img<numimgs; img++)
                {
                    double **temp_hull;
                    MEMORY ps_malloc( &temp_hull, numverts[img]+1, 2 );

                    for( PS_SIT j=0; j<numverts[img]; ++j )
                    {
                        temp_hull[j][0] = convexhulls[img][j*2  ];
                        temp_hull[j][1] = convexhulls[img][j*2+1];
                    }
                    if( numverts[img] )
                    {
                        temp_hull[numverts[img]][0] = convexhulls[img][0];
                        temp_hull[numverts[img]][1] = convexhulls[img][1];
                    }

                    PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name,
                                             true, "mmimages."+OPERA tostring(data_->penaltynames[src])+"."+
                                             OPERA tostring(img)+".dat", temp_hull,numverts[img]+1,2, 0, data_->precision, NULL );

                    MEMORY ps_free( temp_hull, numverts[img]+1 );
                }
            }
        }

        if(/*!vars_->terminatelensing && */data_->penaltymatrix[src][3][0])
        {
            double steps = CONSTANT mm_stepsize_area; // steps in arcseconds
            steps *= data_->arc2pix;

            // compare pairs of images
            for(PS_SIT img1=0; img1<numimgs; img1++)
            {
                for(PS_SIT img2=img1+1; img2<numimgs; img2++)
                {
                    // I fill in the img1 and img2 polygons with uniformly distributed
                    // points so that the area calculation works better. (a bad hack)
                    PS_SIT /*numfillin1 = 0, numfillin2 = 0,*/ numfillin = 0;
                    double minimgx1, maximgx1, minimgy1, maximgy1;
                    double minimgx2, maximgx2, minimgy2, maximgy2;
                    double minx, miny, maxx, maxy;
                    OPERA assign_p_infinity( &minimgx1 );
                    OPERA assign_n_infinity( &maximgx1 );
                    OPERA assign_p_infinity( &minimgy1 );
                    OPERA assign_n_infinity( &maximgy1 );
                    OPERA assign_p_infinity( &minimgx2 );
                    OPERA assign_n_infinity( &maximgx2 );
                    OPERA assign_p_infinity( &minimgy2 );
                    OPERA assign_n_infinity( &maximgy2 );

                    for( PS_SIT j=0; j<numverts[img1]; ++j )
                    {
                        if( convexhulls[img1][j*2  ] < minimgx1 )
                            minimgx1 = convexhulls[img1][j*2  ];
                        if( convexhulls[img1][j*2  ] > maximgx1 )
                            maximgx1 = convexhulls[img1][j*2  ];
                        if( convexhulls[img1][j*2+1] < minimgy1 )
                            minimgy1 = convexhulls[img1][j*2+1];
                        if( convexhulls[img1][j*2+1] > maximgy1 )
                            maximgy1 = convexhulls[img1][j*2+1];
                    }
                    for( PS_SIT j=0; j<numverts[img2]; ++j )
                    {
                        if( convexhulls[img2][j*2  ] < minimgx2 )
                            minimgx2 = convexhulls[img2][j*2  ];
                        if( convexhulls[img2][j*2  ] > maximgx2 )
                            maximgx2 = convexhulls[img2][j*2  ];
                        if( convexhulls[img2][j*2+1] < minimgy2 )
                            minimgy2 = convexhulls[img2][j*2+1];
                        if( convexhulls[img2][j*2+1] > maximgy2 )
                            maximgy2 = convexhulls[img2][j*2+1];
                    }

                    /*
                      minimgx1 = std::ceil( minimgx1 );
                      minimgy1 = std::ceil( minimgy1 );
                      maximgx1 = std::floor( maximgx1 );
                      maximgy1 = std::floor( maximgy1 );
                      minimgx2 = std::ceil( minimgx2 );
                      minimgy2 = std::ceil( minimgy2 );
                      maximgx2 = std::floor( maximgx2 );
                      maximgy2 = std::floor( maximgy2 );
                    */

                    minx = std::min (minimgx1,minimgx2) - steps;
                    miny = std::min (minimgy1,minimgy2) - steps;
                    maxx = std::max (maximgx1,maximgx2) + steps;
                    maxy = std::max (maximgy1,maximgy2) + steps;

                    /*
                      for( PS_SIT x=minimgx1; x<=maximgx1; ++x )
                      for( PS_SIT y=minimgy1; y<=maximgy1; ++y )
                      if( GEOM isinpoly( x, y, convexhulls[img1], numverts[img1] ) &&
                      !GEOM isinpoly( x, y, convexhulls[img2], numverts[img2] ) )
                      ++numfillin1;
                      for( PS_SIT x=minimgx2; x<=maximgx2; ++x )
                      for( PS_SIT y=minimgy2; y<=maximgy2; ++y )
                      if( GEOM isinpoly( x, y, convexhulls[img2], numverts[img2] ) &&
                      !GEOM isinpoly( x, y, convexhulls[img1], numverts[img1] ) )
                      ++numfillin2;
                    */

                    for( double x=minx; x<=maxx; x+=steps )
                        for( double y=miny; y<=maxy; y+=steps )
                            if( GEOM isinpoly( x, y, convexhulls[img1], numverts[img1] ) ||
                                GEOM isinpoly( x, y, convexhulls[img2], numverts[img2] ) )
                                ++numfillin;

                    double areaoverlap=0;

                    // construct triangulation of union of hulls
                    struct triangulateio *triinloc, *trioutloc;
                    MEMORY ps_malloc(            &triinloc , 1 );
                    MEMORY ps_malloc(            &trioutloc, 1 );
                    MEMORY triangulatestructinit( triinloc     );
                    MEMORY triangulatestructinit( trioutloc    );

                    triinloc->numberofpointattributes = 0;
                    triinloc->numberofpoints = numverts[img1] + numverts[img2] + numfillin;
                    //numfillin1 + numfillin2;
                    MEMORY ps_malloc( &triinloc->pointlist, triinloc->numberofpoints*2 );

                    // add points
                    for(PS_SIT j=0; j<numverts[img1]; j++)
                    {
                        triinloc->pointlist[j*2  ] = convexhulls[img1][j*2  ];
                        triinloc->pointlist[j*2+1] = convexhulls[img1][j*2+1];
                    }
                    for(PS_SIT j=0; j<numverts[img2]; j++)
                    {
                        triinloc->pointlist[(numverts[img1]+j)*2  ] = convexhulls[img2][j*2  ];
                        triinloc->pointlist[(numverts[img1]+j)*2+1] = convexhulls[img2][j*2+1];
                    }
                    PS_SIT currpointlistind = numverts[img1] + numverts[img2];
                    for( double x=minx; x<=maxx; x+=steps )
                        for( double y=miny; y<=maxy; y+=steps )
                            if( GEOM isinpoly( x, y, convexhulls[img1], numverts[img1] ) ||
                                GEOM isinpoly( x, y, convexhulls[img2], numverts[img2] ) )
                            {
                                triinloc->pointlist[currpointlistind*2  ] = x;
                                triinloc->pointlist[currpointlistind*2+1] = y;
                                ++currpointlistind;
                            }
                    /*
                      for( PS_SIT x=minimgx1; x<=maximgx1; ++x )
                      for( PS_SIT y=minimgy1; y<=maximgy1; ++y )
                      if( GEOM isinpoly( x, y, convexhulls[img1], numverts[img1] ) &&
                      !GEOM isinpoly( x, y, convexhulls[img2], numverts[img2] ) )
                      {
                      triinloc->pointlist[currpointlistind*2  ] = x;
                      triinloc->pointlist[currpointlistind*2+1] = y;
                      ++currpointlistind;
                      }
                      for( PS_SIT x=minimgx2; x<=maximgx2; ++x )
                      for( PS_SIT y=minimgy2; y<=maximgy2; ++y )
                      if( GEOM isinpoly( x, y, convexhulls[img2], numverts[img2] ) &&
                      !GEOM isinpoly( x, y, convexhulls[img1], numverts[img1] ) )
                      {
                      triinloc->pointlist[currpointlistind*2  ] = x;
                      triinloc->pointlist[currpointlistind*2+1] = y;
                      ++currpointlistind;
                      }
                    */

                    // triangulate points
                    ps_tri_triangulate((char*)CONSTANT triswitchesnominangle,
                                       triinloc, trioutloc, (struct triangulateio*)NULL);

                    // test triangles to see if contained in both
                    double pos[2][3];
                    for(PS_SIT tri = 0; tri < trioutloc->numberoftriangles; tri++)
                    {
                        for(PS_SIT f=0; f<3; f++)
                        {
                            pos[0][f] = triinloc->pointlist[trioutloc->trianglelist[tri*3+f]*2  ];
                            pos[1][f] = triinloc->pointlist[trioutloc->trianglelist[tri*3+f]*2+1];
                        }
                        double xctr = (pos[0][0]+pos[0][1]+pos[0][2]) / 3.0;
                        double yctr = (pos[1][0]+pos[1][1]+pos[1][2]) / 3.0;

                        if(GEOM isinpoly(xctr,yctr,convexhulls[img1],numverts[img1]) &&
                           GEOM isinpoly(xctr,yctr,convexhulls[img2],numverts[img2]))
                        {
                            double coordlist[6] = { pos[0][0], pos[1][0],
                                                    pos[0][1], pos[1][1],
                                                    pos[0][2], pos[1][2] };
                            areaoverlap += GEOM areatriangle(coordlist);
                        }
                    }

                    MEMORY triangulatestructdestruct( triinloc  );
                    MEMORY triangulatestructdestruct( trioutloc );
                    MEMORY ps_free(                   triinloc  );
                    MEMORY ps_free(                   trioutloc );

                    double frac1 = areaoverlap / areas[img1];
                    double frac2 = areaoverlap / areas[img2];

                    /*
                      if( frac1 > 1 )
                      frac1 = 1;
                      if( frac2 > 1 )
                      frac2 = 1;
                    */
                    if (OPERA equalszero (frac1-1))
                        frac1=1;
                    if (OPERA equalszero (frac2-1))
                        frac2=1;
                    if (frac1>=1 || frac2>=1)
                        frac1 = frac2 = 1.0;

                    if(!areaoverlap)
                    {
                        OPERA assign_p_infinity( &frac1 );

                        for(PS_SIT i=0; i<numverts[img1]; ++i)
                        {
                            for(PS_SIT j=0; j<numverts[img2]; ++j)
                            {
                                double dist2 = OPERA distance2(convexhulls[img1][i*2  ],
                                                               convexhulls[img1][i*2+1],
                                                               convexhulls[img2][j*2  ],
                                                               convexhulls[img2][j*2+1]);
                                if(dist2 < frac1)
                                    frac1 = dist2;
                            }
                        }

                        frac1 *= -data_->pix2arc;
                        frac2 = frac1;
                    }


                    // penalize
                    penaltypenalizer(data_,cdata_,vars_, src,3,frac1,"mismatched images");
                    penaltypenalizer(data_,cdata_,vars_, src,3,frac2,"mismatched images");

                    //if(vars_->terminatelensing)
                    //    break;
                }

                //if(vars_->terminatelensing)
                //    break;
            }
        }

        if(/*!vars_->terminatelensing && */data_->penaltymatrix[src][2][0])
        {
            for(PS_SIT img=0; img<numimgs; img++)
            {
                double maxdist = GEOM getmaxdist(convexhulls[img],numverts[img]);
                // q = a/b = a/(area/(pi*a)) = pi*a^2/area = pi*maxdist^2/(4*area)
                double axisratio = CONSTANT pi*maxdist*maxdist/(4.0*areas[img]);
                if(axisratio < 1)
                    axisratio = 1;

                penaltypenalizer(data_,cdata_,vars_, src,2,axisratio,"axis ratio");

                //if(vars_->terminatelensing)
                //    break;
            }
        }

        if(/*!vars_->terminatelensing && */data_->penaltymatrix[src][4][0])
            // delete this later
            //if(src==3)
        {
            //std::cout << "src = " << src << " -- numimages = " << numimgs << std::endl;
            struct triangulateio *triinloc, *trioutloc;
            MEMORY ps_malloc(            &triinloc , 1 );
            MEMORY ps_malloc(            &trioutloc, 1 );
            MEMORY triangulatestructinit( triinloc     );
            MEMORY triangulatestructinit( trioutloc    );

            if(makesepconvexhulls)
            {
                PS_SIT numpts = 0;
                for(PS_SIT img=0; img<numimgs; img++)
                    numpts += numverts[img];

                triinloc->numberofpointattributes = 0;
                triinloc->numberofpoints = numpts;
                MEMORY ps_malloc( &(triinloc->pointlist), numpts*2 );

                for(PS_SIT img=0; img<numimgs; img++)
                {
                    PS_SIT startpos = 0;
                    for(PS_SIT img2=0; img2<img; img2++)
                        startpos += numverts[img2];

                    for(PS_SIT j=0; j<numverts[img]; j++)
                    {
                        triinloc->pointlist[(startpos+j)*2  ] = convexhulls[img][j*2  ];
                        triinloc->pointlist[(startpos+j)*2+1] = convexhulls[img][j*2+1];
                    }
                }
            }
            else
            {
                PS_SIT numpts=0;
                for(PS_SIT r=0; r<data_->ndp; r++)
                    if(data_->mmimages[src][r])
                        numpts++;

                // also add border pixels
                for (PS_SIT b=0; b<numimgs; ++b)
                {
                    numpts += (PS_SIT)data_->mmborder[src][b][0];
                }

                triinloc->numberofpointattributes = 0;
                triinloc->numberofpoints = numpts;
                MEMORY ps_malloc( &(triinloc->pointlist), numpts*2 );

                PS_SIT index=0;
                for(PS_SIT r=0; r<data_->ndp; r++)
                    if(data_->mmimages[src][r])
                    {
                        triinloc->pointlist[index*2  ] = vars_->newloc[r*2  ];
                        triinloc->pointlist[index*2+1] = vars_->newloc[r*2+1];
                        index++;
                    }

                // also add border pixels
                for (PS_SIT b=0; b<numimgs; ++b)
                {
                    for(PS_SIT vert=0; vert<(PS_SIT)data_->mmborder[src][b][0]; ++vert)
                    {
                        triinloc->pointlist[index*2  ] = data_->mmborder_defl[src][b][vert*2+0+1];
                        triinloc->pointlist[index*2+1] = data_->mmborder_defl[src][b][vert*2+1+1];
                        index++;
                    }
                }
            }

            PS_SIT numv;
            double *chull;
            GEOM getconvexhull( triinloc,trioutloc,&numv,&chull,NULL );

            MEMORY triangulatestructdestruct( triinloc  );
            MEMORY triangulatestructdestruct( trioutloc );
            MEMORY ps_free(                   triinloc  );
            MEMORY ps_free(                   trioutloc );

            double area = GEOM areapoly(chull,numv)*(data_->wcsinfo[4]*data_->wcsinfo[4] +
                                                     data_->wcsinfo[5]*data_->wcsinfo[5])*3600*3600;

            MEMORY ps_free( chull );

            penaltypenalizer(data_,cdata_,vars_, src,4,area,"convex hull size");

            //if(vars_->terminatelensing)
            //    break;
        }

        if(makesepconvexhulls)
        {
            MEMORY ps_free( convexhulls, numimgs );

            //if(vars_->terminatelensing)
            //    break;
        }
    }

    //if( !data_->doreconstruction )
    //    vars_->terminatelensing = 1;
    if( vars_->terminatelensing )
        vars_->endedearly=1;

    if(data_->verbose != 1 && !vars_->terminatelensing)
    {
        for (PS_SIT src=0; src<data_->extlengths[11]; ++src)
        {
            if(data_->penaltymatrix[src][3][0]==2)
            {
                PRINTER print2screen(data_->print2screenname,
                                     OPERA tostring (data_->penaltynames[src]) + ": " +
                                     "model passed mismatched images test",
                                     cdata_->print2screenmutex);
            }
            else if( data_->penaltymatrix[src][3][0] )
            {
                PRINTER print2screen(data_->print2screenname,
                                     OPERA tostring (data_->penaltynames[src]) + ": " +
                                     "mismatched images penalty = " +
                                     OPERA tostring(vars_->penalties[3]),
                                     cdata_->print2screenmutex);
            }
            if(data_->penaltymatrix[src][2][0]==2)
            {
                PRINTER print2screen(data_->print2screenname,
                                     OPERA tostring (data_->penaltynames[src]) + ": " +
                                     "model passed axis ratio test",
                                     cdata_->print2screenmutex);
            }
            else if( data_->penaltymatrix[src][2][0] )
            {
                PRINTER print2screen(data_->print2screenname,
                                     OPERA tostring (data_->penaltynames[src]) + ": " +
                                     "axis ratio penalty = " +
                                     OPERA tostring(vars_->penalties[2]),
                                     cdata_->print2screenmutex);
            }
            if(data_->penaltymatrix[src][4][0]==2)
            {
                PRINTER print2screen(data_->print2screenname,
                                     OPERA tostring (data_->penaltynames[src]) + ": " +
                                     "model passed convex hull size test",
                                     cdata_->print2screenmutex);
            }
            else if( data_->penaltymatrix[src][4][0] )
            {
                PRINTER print2screen(data_->print2screenname,
                                     OPERA tostring (data_->penaltynames[src]) + ": " +
                                     "convex hull size penalty = " +
                                     OPERA tostring(vars_->penalties[4]),
                                     cdata_->print2screenmutex);
            }
        }
    }
}

void pixsrc_common::sourcereconstructions(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if (1!=data_->verbose)
        PRINTER print2screen(data_->print2screenname,
                             "setting up matrices for de-lensing",
                             cdata_->print2screenmutex);

    // handle uv-plane data
    bool is_uv = data_->is_uvdata;
    PS_SIT sh_ndp = data_->uv_ndp;
    MATRIX *uv_lo_ptr=NULL, *uv_lo=NULL;
    MATRIX *uv_tmplo_ptr=NULL, *uv_tmplo=NULL;
    VECTOR *uv_data_ptr=NULL, *uv_data=NULL;
    VECTOR *uv_model_ptr=NULL, *uv_model=NULL;
    vars_->uv_residualptr=NULL;
    vars_->uv_residual=NULL;

    MEMORY ps_malloc (&vars_->dataptr     , 1);
    MEMORY ps_malloc (&vars_->mpsptr      , 1);
    MEMORY ps_malloc (&vars_->lensedmpsptr, 1);
    MEMORY ps_malloc (&vars_->residualptr , 1);
    vars_->data      = new (vars_->dataptr)      VECTOR (cdata_, data_, vars_->lonr);
    vars_->mps       = new (vars_->mpsptr)       VECTOR (cdata_, data_, vars_->lonc);
    vars_->lensedmps = new (vars_->lensedmpsptr) VECTOR (cdata_, data_, vars_->lonr);
    vars_->residual  = new (vars_->residualptr)  VECTOR (cdata_, data_, vars_->lonr);

    // if using image-plane data
    PS_SIT counter = 0;
    for(PS_SIT x = 0; x < data_->ndp; x++)
        if(vars_->r4r[x]!=-1)
            vars_->data->set(counter++,data_->data[x]);

    // if using uv-plane data
    if (is_uv)
    {
        if (1==data_->transmode)
        {}
        else
        {
            // allocate
            PS_SIT initsize1 = !data_->numgpu2use ? vars_->lensingoperator->get_nnz() : 0;
            MEMORY ps_malloc (&uv_lo_ptr, 1);
            MEMORY ps_malloc (&uv_data_ptr, 1);
            MEMORY ps_malloc (&uv_model_ptr, 1);
            MEMORY ps_malloc (&vars_->uv_residualptr, 1);
            uv_lo        = new (uv_lo_ptr)        MATRIX (cdata_, data_, sh_ndp*2,
                                                          vars_->lonc, initsize1, 1, data_);
            uv_data      = new (uv_data_ptr)      VECTOR (cdata_, data_, sh_ndp*2);
            uv_model     = new (uv_model_ptr)     VECTOR (cdata_, data_, sh_ndp*2);
            vars_->uv_residual = new (vars_->uv_residualptr) VECTOR (cdata_, data_, sh_ndp*2);

            // copy uv data
            for (PS_SIT x=0; x<sh_ndp*2; ++x)
                uv_data->set (x, data_->uv_data[x]);

            // transform lensing operator into uv space
            ((MATRIX*)data_->uv_transform_mat)->mult (vars_->lensingoperatornobo, uv_lo, 0, 0, cdata_->numthreads);
        }
    }

    // set pointer and var
    if (is_uv && 1==data_->transmode)
    {}
    else
    {
        vars_->lonr4stat     = is_uv ? sh_ndp*2           : vars_->lonr;
        vars_->data4stat     = is_uv ? uv_data            : vars_->data;
        vars_->lo4stat       = is_uv ? uv_lo              : vars_->lensingoperator;
        vars_->residual4stat = is_uv ? vars_->uv_residual : vars_->residual;
    }

    if( !data_->usersetsrc || data_->reg )
    {
        MEMORY ps_malloc( &(vars_->b1ptr ), 1 );
        MEMORY ps_malloc( &(vars_->lodptr), 1 );

        PS_SIT initsize1 = ( !data_->numgpu2use ) ? vars_->lensingoperator->get_nnz() : 0;

        vars_->b1  = new (vars_->b1ptr ) MATRIX( cdata_, data_, vars_->lonc,
                                                 vars_->lonc, initsize1, 1,  data_ );

        vars_->lod = new (vars_->lodptr) VECTOR( cdata_, data_, vars_->lonc);

        if (!is_uv)
        {
            vars_->lensingoperator->mult( vars_->lensingoperator,vars_->b1,
                                          1, 0, cdata_->numthreads       );
            vars_->b1->mult( 1.0 / vars_->datavariance );
        }
        else
        {
            if (1==data_->transmode)
            {
                MATRIX uvtmpmat (cdata_, data_, vars_->lonc, vars_->lonr, -1, -1, data_);
                vars_->lensingoperatornobo->mult ((MATRIX*)data_->uv_transform_mat, &uvtmpmat,
                                                  1, 0, cdata_->numthreads);
                uvtmpmat.mult (vars_->lensingoperatornobo, vars_->b1, 0, 0, cdata_->numthreads);
            }
            else
            {
                // uv_tmplo = lensingoperator^T . noisecovariance^-1
                MEMORY ps_malloc (&uv_tmplo_ptr, 1);
                uv_tmplo = new (uv_tmplo_ptr) MATRIX (cdata_, data_, vars_->lonc, sh_ndp*2,
                                                      uv_lo->get_nnz(), 1, data_);
                for (PS_SIT m=0; m<sh_ndp*2; ++m)
                    for (PS_SIT n=0; n<vars_->lonc; ++n)
                        uv_tmplo->set (n, m, uv_lo->get(m,n)/data_->uv_newvariance[m/2]);
                uv_tmplo->mult (uv_lo, vars_->b1, 0, 0, cdata_->numthreads);
            }
        }

        if( data_->debug )
        {
            PRINTER print( vars_->tracker, cdata_->basename, data_->name, true,
                           "b1.MATRIX", vars_->lonc, vars_->lonc, vars_->b1, 0);
        }

        if (!is_uv)
        {
            vars_->lensingoperator->mult( vars_->data, vars_->lod,
                                          1, cdata_->numthreads          );
        }
        else
        {
            if (1==data_->transmode)
            {
                vars_->lensingoperatornobo->mult ((VECTOR*)data_->uv_transform_vec, vars_->lod,
                                                  1, cdata_->numthreads);
            }
            else
            {
                uv_tmplo->mult (uv_data, vars_->lod, 0, cdata_->numthreads);
            }
        }

    }

    // if user defined source or analytic source optimization needed
    ANALYTICSRC getsource( data_, cdata_, vars_ );

    // wait for regularization to return
    PS_SIT dim = 1;
    PS_SIT waitlist[1] = {0};
    OPERA pthreadswait(vars_,dim,waitlist);

    if( !data_->usersetsrc || data_->reg )
    {
        if( data_->traceparams[2]>0 )
        {
            vars_->lambda1 = data_->traceparams[0] *
                std::pow(data_->traceparams[1],vars_->tracker);
        }
        else if( data_->lambdaguess<0 )
        {
            vars_->lambda1 = -data_->lambdaguess;
        }
        else if( data_->lambdaguess>0 )
        {
            findlambda( data_, cdata_, vars_ );
        }
        else
        {
            vars_->lambda1 = 0;
        }
    }
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->lambdaguess && ( !data_->usersetsrc || data_->reg ) )
    {
        MEMORY ps_malloc( &(vars_->a1ptr), 1 );

        PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->b1->get_nnz() : 0;

        vars_->a1 = new (vars_->a1ptr)
            MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 1, data_ );
        vars_->b1->plus( vars_->c1,vars_->a1, 0, 0, 1,
                         vars_->lambda1, cdata_->numthreads );

    }

    pthread_create(&(vars_->pthreads[2]), cdata_->attrdetached, STATISTIC computedeta, vars_->dav);

    if( !data_->usersetsrc || data_->reg )
    {
        if( data_->debug && data_->lambdaguess )
        {
            PRINTER print( vars_->tracker, cdata_->basename, data_->name, true,
                           "a1.MATRIX", vars_->lonc, vars_->lonc, vars_->a1, 0);
        }

        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "de-lensing source",
                                 cdata_->print2screenmutex);
    }

    if (!data_->hacksrc)
    {
    if( data_->lambdaguess && ( !data_->usersetsrc || data_->reg ) )
    {
        /*********************************************************************************/
        /*********************************START DELENSING*********************************/
        /*********************************************************************************/

        // START most straightforward method

        /*
          void *clock_ptr;
          double elapsed_time;
          CUDA ps_start_clock( &clock_ptr );
        */
        if(trieqon&&vars_->a1->are_we_using_cuda)
        {

#ifdef __USE_PIXSRC_CUDA__
            CUDA tri_eq_solve( vars_->a1, vars_->mps, vars_->lod ,1);
#endif

        }
        else
        {
            vars_->a1->linequationsolve(vars_->mps,vars_->lod, 1);
        }
        /*
          CUDA ps_stop_clock( clock_ptr, &elapsed_time );
          string et = "pixsrc took " + OPERA tostring(elapsed_time) + " ms to lineqsolve\n";
          //PRINTER print2screen("", et, cdata_->print2screenmutex);
          */

        // END most straightforward method

        // START splitting up most straightforward method into many steps
        // could possibly use this for multithreading
        // but this step seems somewhat fast .. test it in free time!
        // if you do choose to use this, make it general so the no-regularization
        // case can use it too
        /*
          MATRIX *alo = new MATRIX(vars_->lonc, vars_->lonr, 0);
          VECTOR *temp = new VECTOR(vars_->lonc);
          VECTOR *tempx = new VECTOR(vars_->lonc);
          for(int u=0; u<vars_->lonr; u++)
          {
          for(int h=0; h<vars_->lonc; h++)
          temp->set(h,vars_->lensingoperator->get(u,h));
          vars_->a1->linequationsolve(tempx,*temp);
          for(int g=0; g<vars_->lonc; g++)
          alo->set(g,u,tempx->get(g));
          }
          alo->mult(*(vars_->data),vars_->mps,false);
          delete alo;
          delete temp;
          delete tempx;
        */
        // END splitting up most straightforward method into many steps

        if (!data_->uv_newvariance)
            vars_->mps->mult( 1.0 / vars_->datavariance );

        /*********************************************************************************/
        /**********************************END DELENSING**********************************/
        /*********************************************************************************/
    }
    else if( !data_->usersetsrc || data_->reg )
    {
        if(trieqon&&vars_->b1->are_we_using_cuda)
        {

#ifdef __USE_PIXSRC_CUDA__
            CUDA tri_eq_solve( vars_->b1, vars_->mps, vars_->lod, 1 );
#endif

        }
        else
        {
            vars_->b1->linequationsolve(vars_->mps,vars_->lod, 1);
        }

        if (!data_->uv_newvariance)
            vars_->mps->mult( 1.0 / vars_->datavariance );

        // not using regularization is a bad idea,
        // but I try and help a little here if it messed up
        char printit = 0;
        for( PS_SIT g=0; g<vars_->lonc; ++g )
            if( !OPERA is_finite(vars_->mps->get(g)) )
            {
                printit = 1;
                vars_->mps->set( g, 0 );
            }
        if( printit )
        {
            PRINTER printwarning(data_->print2screenname,
                                 "one or more source plane pixels were not finite; "
                                 "they were set to zero flux manually.",
                                 cdata_->print2screenmutex);
        }
    }
    }
    else
    {
	if (1!=data_->verbose)
	    PRINTER print2screen(data_->print2screenname,
				 "applying user-supplied hacksrc parameters",
				 cdata_->print2screenmutex);

	// set source vector to use-supplied parameters
	for (PS_SIT g=0; g<vars_->lonc; ++g)
	    vars_->mps->set (g,data_->hacksrc[g]);
    }

    //vars_->mps->set(0,1.02317558531355492212*4.63621); // uv-gauss-norm
    /*
      vars_->mps->set(0,1);
      for (PS_SIT k=1; k<vars_->lonc; ++k)
      vars_->mps->set(k,0);
    */
    //vars_->mps->set(0,1); //gug
    /*
      for (PS_SIT c=0; c<vars_->lonc; ++c)
      if (1||std::abs(vars_->mps->get(c))<1e100)
      vars_->mps->set(c,0);
      vars_->mps->set(2*4+1,1);
      vars_->mps->set(2*4+2,1);
      vars_->mps->set(3*4+3,1);
    */
    //vars_->mps->set(28,1);
    //vars_->mps->set(1,0.35);
    //vars_->mps->set(11,0.2);



    if(data_->debug && data_->reg )
    {
        VECTOR *dummy1;
        MEMORY ps_malloc( &dummy1, 1 );

        VECTOR *reg = new (dummy1) VECTOR( cdata_, data_, vars_->lonc );
        if(vars_->h1)
        {
            vars_->h1->mult( vars_->mps,reg,false, cdata_->numthreads);
            double **vec;
            MEMORY ps_malloc( &vec, vars_->lonc, 3 );
            if(data_->gridtype)
            {
                for(PS_SIT x=0; x<vars_->lonc; x++)
                {
                    vec[x][0] = vars_->triout->pointlist[x*2  ];
                    vec[x][1] = vars_->triout->pointlist[x*2+1];
                }
            }
            else
            {
                for(PS_SIT x=0; x<vars_->lonc; x++)
                {
                    vec[x][0] = vars_->srcloc[vars_->c4cback[x]*2  ];
                    vec[x][1] = vars_->srcloc[vars_->c4cback[x]*2+1];
                }
            }

            for(PS_SIT x=0; x<vars_->lonc; x++)
                vec[x][2] = std::fabs(reg->get(x));
            PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, true,
                                     "reg.VECTOR",vec, vars_->lonc, 3, 0, data_->precision, NULL);
            MEMORY ps_free( vec, vars_->lonc );
        }
        else if(vars_->h1x && vars_->h1y)
        {
            VECTOR *dummy2, *dummy3;
            MEMORY ps_malloc( &dummy2, 1 );
            MEMORY ps_malloc( &dummy3, 1 );

            VECTOR *regx = new (dummy2) VECTOR(cdata_, data_, vars_->lonc);
            VECTOR *regy = new (dummy3) VECTOR(cdata_, data_, vars_->lonc);
            vars_->h1x->mult( vars_->mps,regx,false, cdata_->numthreads);
            vars_->h1y->mult( vars_->mps,regy,false, cdata_->numthreads);
            double **vec;
            MEMORY ps_malloc( &vec, vars_->lonc, 3 );
            for(PS_SIT x=0; x<vars_->lonc; x++)
            {
                vec[x][0] = vars_->triout->pointlist[x*2  ];
                vec[x][1] = vars_->triout->pointlist[x*2+1];
            }
            for(PS_SIT x=0; x<vars_->lonc; x++)
                vec[x][2] = std::sqrt(regx->get(x)*regx->get(x)+regy->get(x)*regy->get(x));
            PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                     "reg.VECTOR",vec, vars_->lonc, 3, 0, data_->precision, NULL);
            MEMORY ps_free( vec, vars_->lonc );

            regx->~VECTOR();
            regy->~VECTOR();
            MEMORY ps_free( dummy2 );
            MEMORY ps_free( dummy3 );
        }

        reg->~VECTOR();
        MEMORY ps_free( dummy1 );
    }

    // NOTE: magnification penalty computed in getmagnification thread

    /*
    // get unblurred, uncorrected lens model
    PS_SIT getnobo = data_->interperr || data_->printvec || data_->magparams || data_->penaltyquery[1];
    if (getnobo)
    {
    MEMORY ps_malloc (&vars_->lensedmpsnoboptr, 1);
    vars_->lensedmpsnobo = new (vars_->lensedmpsnoboptr) VECTOR (cdata_, data_, vars_->lonr);
    COMMON lensgalaxy (data_, cdata_, vars_, vars_->lensingoperatornobo, vars_->mps, vars_->lensedmpsnobo);
    }
    */

    if (data_->interperr)
    {
        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "fitting TPS to source for interpolation errors",
                                 cdata_->print2screenmutex);

        TPS *tpsfitptr, *tpsfit;
        MEMORY ps_malloc (&tpsfitptr, 1);
        tpsfit = new (tpsfitptr) TPS (cdata_, data_, vars_->triout->pointlist, vars_->triout->pointlist+1,
                                      vars_->mps->get_vec_ptr(), 2, 2, 1,
                                      vars_->triout->numberofpoints, ps_rbf_tps, 0, 0, 1, data_->transmode);
        tpsfit->get_tps_weights ();
        tpsfit->cleanup ();

        vars_->tpsfitptr = (void*)tpsfitptr;
        vars_->tpsfit    = (void*)tpsfit;

        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "done fitting TPS to source for interpolation errors",
                                 cdata_->print2screenmutex);
    }

    COMMON lensgalaxy( data_, cdata_, vars_, vars_->lensingoperator,
                       vars_->mps, vars_->lensedmps);

    vars_->data->minus(vars_->lensedmps, vars_->residual, 1, 1);

    // get and print uv model and residuals
    if (is_uv)
    {
        if (1==data_->transmode)
        {
            COMMON printuvresiduals (data_, cdata_, vars_);
        }
        else
        {
            uv_lo->mult (vars_->mps, uv_model, 0, cdata_->numthreads);
            uv_data->minus (uv_model, vars_->uv_residual, 1, 1);
            COMMON printuvplane (data_, cdata_, vars_, uv_data, uv_model, vars_->uv_residual);
            uv_lo->~MATRIX();
            uv_data->~VECTOR();
            uv_model->~VECTOR();
            MEMORY ps_free (uv_lo_ptr);
            MEMORY ps_free (uv_data_ptr);
            MEMORY ps_free (uv_model_ptr);
            if (uv_tmplo_ptr)
            {
                uv_tmplo->~MATRIX();
                MEMORY ps_free (uv_tmplo_ptr);
            }
        }
    }

    pthread_create(&(vars_->pthreads[4]), cdata_->attrdetached,
                   STATISTIC computees      , vars_->dav        );
    pthread_create(&(vars_->pthreads[5]), cdata_->attrdetached,
                   COMMON computesizepenalty, vars_->dav        );
    pthread_create(&(vars_->pthreads[6]), cdata_->attrdetached,
                   COMMON computepenalty    , vars_->dav        );
    pthread_create(&(vars_->pthreads[3]), cdata_->attrdetached,
                   STATISTIC computeed, vars_->dav              );
}

void pixsrc_common::lensgalaxy (inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                MATRIX *lo, VECTOR *gal, VECTOR *img)
{
    // if correcting flux for interpolation errors
    if (data_->interperr && 1==data_->irscheme && vars_->tpsfit)
        STATISTIC computeinterperrors (data_, cdata_, vars_);

    lo->mult( gal, img, 0, cdata_->numthreads );

    // if not modifying flux
    if (!(
            (data_->interperr && 1==data_->irscheme && vars_->tpsfit) ||
            (vars_->srcexactvec)
            ))
        return;


    VECTOR addme  (cdata_, data_, vars_->lonr);
    VECTOR addme2 (cdata_, data_, vars_->lonr);

    // if correcting flux for interpolation errors
    if (data_->interperr && 1==data_->irscheme && vars_->tpsfit)
        addme.plus (vars_->interperr, NULL, 1, 1);
    // if correcting flux using analytic source
    if (vars_->srcexactvec)
        addme.plus (vars_->srcexactvec, NULL, 1, 1);

    if (lo!=vars_->lensingoperatornobo)
    {
        vars_->blurringoperator->mult (&addme, &addme2, 0, cdata_->numthreads);
        img->plus (&addme2, NULL, 1, 1);
    }
    else
    {
        img->plus (&addme, NULL, 1, 1);
    }


}

struct amoebalambdastruct
{
    dataandvars *dav;
    MATRIX *a2, *a2d;
    VECTOR *mps2, *mps2d;
    VECTOR *reg2, *reg2d;
    VECTOR *lmps2, *lmps2d;
    VECTOR *res2, *res2d;
    PS_SIT count;
};
double getamoebalambdaevi (const gsl_vector *gvec, void *args)
{
    double lambda2 = gsl_vector_get (gvec, 0);
    lambda2 = std::pow (10.0, lambda2);

    if (lambda2<0)
        return 1e100 * ( 1 - lambda2 );

    amoebalambdastruct *als = (amoebalambdastruct*)args;
    ++als->count;
    dataandvars *dav = als->dav;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;
    ///*
    MATRIX *a2 = als->a2;
    VECTOR *mps2 = als->mps2;
    VECTOR *reg2 = als->reg2;
    VECTOR *lmps2 = als->lmps2;
    VECTOR *res2 = als->res2;
    MATRIX *a2d = als->a2d;
    VECTOR *mps2d = als->mps2d;
    VECTOR *reg2d = als->reg2d;
    VECTOR *lmps2d = als->lmps2d;
    VECTOR *res2d = als->res2d;

    if(!vars_->b1->are_we_using_cuda)
    {
        a2->~MATRIX();
        mps2->~VECTOR();
        reg2->~VECTOR();
        if (!(data_->is_uvdata && 1==data_->transmode))
        {
            lmps2->~VECTOR();
            res2->~VECTOR();
        }
        a2 = new (a2d) MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, 0, 0, data_ );
        mps2 = new (mps2d) VECTOR( cdata_, data_, vars_->lonc );
        reg2 = new (reg2d) VECTOR( cdata_, data_, vars_->lonc );
        if (!(data_->is_uvdata && 1==data_->transmode))
        {
            lmps2 = new (lmps2d) VECTOR( cdata_, data_, vars_->lonr4stat );
            res2 = new (res2d) VECTOR( cdata_, data_, vars_->lonr4stat );
        }
    }

    bool zeroit = ( als->count ) ? 0 : 1;
    //*/
    double detCdTerm=0,detATerm=0,detCTerm=0,EdTerm=0,EsTerm=0,lambdaTerm=0,PiTerm=0;
    /*
      MATRIX a2_( cdata_, data_, vars_->lonc, vars_->lonc, 0, 0, data_ );
      VECTOR mps2_( cdata_, data_, vars_->lonc );
      VECTOR reg2_( cdata_, data_, vars_->lonc );
      VECTOR lmps2_( cdata_, data_, vars_->lonr );
      VECTOR res2_( cdata_, data_, vars_->lonr );
      MATRIX *a2 = &a2_;
      VECTOR *mps2 = &mps2_;
      VECTOR *reg2 = &reg2_;
      VECTOR *lmps2 = &lmps2_;
      VECTOR *res2 = &res2_;
    */
    //void *ct;
    //float time;
    //CUDA ps_start_clock( &ct );

    vars_->b1->plus( vars_->c1, a2, 0, 0, 1, lambda2, cdata_->numthreads);
    //CUDA ps_stop_clock( ct, &time );
    //std::cout << "amoeba time-1 " << time << std::endl;
    //CUDA ps_start_clock( &ct );
    a2->linequationsolve(mps2,vars_->lod,zeroit);
    //CUDA ps_stop_clock( ct, &time );
    //std::cout << "amoeba time0 " << time << std::endl;
    //CUDA ps_start_clock( &ct );
    if (!data_->uv_newvariance)
        mps2->mult( 1.0 / vars_->datavariance );
    if (!(data_->is_uvdata && 1==data_->transmode))
    {
        vars_->lo4stat->mult( mps2, lmps2, 0, 1 );
        vars_->data4stat->minus( lmps2, res2, 1, 1 );
    }

    vars_->c1->mult( mps2, reg2, 0, 1 );

    //CUDA ps_stop_clock( ct, &time );
    //std::cout << "amoeba time1 " << time << std::endl;
    //CUDA ps_start_clock( &ct );
    detATerm = -a2->logdet()/2.0;
    //CUDA ps_stop_clock( ct, &time );
    //std::cout << "amoeba time2 " << time << std::endl;
    //CUDA ps_start_clock( &ct );
    EsTerm = -lambda2*reg2->innerproduct(mps2)/2.0;
    if (!(data_->is_uvdata && 1==data_->transmode))
    {
        if (!data_->uv_newvariance)
        {
            EdTerm = -0.5*res2->innerproduct(res2)/vars_->datavariance;
        }
        else
        {
            for (PS_SIT s=0; s<vars_->lonr4stat; ++s)
                EdTerm += res2->get(s)*res2->get(s) / data_->uv_newvariance[s/2];
            EdTerm *= -0.5;
        }
    }
    else
    {
        VECTOR uvtmpvec (cdata_, data_, vars_->lonc);
        vars_->b1->mult (mps2, &uvtmpvec, 0, cdata_->numthreads);
        EdTerm = uvtmpvec.innerproduct (mps2) +
            data_->uv_transform_scalar - 2.0 *mps2->innerproduct (vars_->lod);
        EdTerm *= -0.5;
    }
    lambdaTerm = vars_->lonc/2.0*std::log(lambda2);

    //CUDA ps_stop_clock( ct, &time );
    //std::cout << "amoeba time3 " << time << std::endl;
    //std::cout << als->count << " " << lambda2 << " " << detATerm << " " << detCTerm << " " << detCdTerm << " " << EsTerm << " " << EdTerm << " " << lambdaTerm << " " << PiTerm << std::endl;

    double chi2 = -2.0*(detATerm + detCTerm + detCdTerm +
                        EsTerm + EdTerm + lambdaTerm + PiTerm);

    if (OPERA is_finite (chi2))
        return chi2;
    else
        return 1.0e100;
}
void amoeba_lambda(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    amoebalambdastruct als;

    // allocate matrices and vectors once only, to speed up optimization
    MATRIX *a2dummy;
    MEMORY ps_malloc( &a2dummy, 1 );
    VECTOR *mps2d, *lmps2d, *res2d, *reg2d;
    MEMORY ps_malloc( &mps2d, 1 );
    MEMORY ps_malloc( &lmps2d, 1 );
    MEMORY ps_malloc( &res2d, 1 );
    MEMORY ps_malloc( &reg2d, 1 );
    MATRIX *a2 = new (a2dummy) MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, 0, 0, data_ );
    VECTOR *mps2 = new (mps2d) VECTOR( cdata_, data_, vars_->lonc );
    VECTOR *reg2 = new (reg2d) VECTOR( cdata_, data_, vars_->lonc );
    VECTOR *lmps2=NULL, *res2=NULL;

    if (!(data_->is_uvdata && 1==data_->transmode))
    {
        lmps2 = new (lmps2d) VECTOR( cdata_, data_, vars_->lonr4stat );
        res2 = new (res2d) VECTOR( cdata_, data_, vars_->lonr4stat );
    }

    als.dav = vars_->dav;
    als.a2 = a2;
    als.mps2 = mps2;
    als.reg2 = reg2;
    als.lmps2 = lmps2;
    als.res2 = res2;
    als.a2d = a2dummy;
    als.mps2d = mps2d;
    als.reg2d = reg2d;
    als.lmps2d = lmps2d;
    als.res2d = res2d;
    als.count = 0;

    PS_SIT lnvary = 1;
    gsl_vector *x  = gsl_vector_alloc (lnvary);
    gsl_vector *ss = gsl_vector_alloc (lnvary);
    gsl_vector_set_all (x,  0);
    gsl_vector_set_all (ss, 0);
    gsl_vector_set     (x,  0, std::log10(data_->lambdaguess));
    gsl_vector_set     (ss, 0, 2);

    OPERA multidimmin (data_, cdata_, lnvary, data_->regftol, "log regularization strength",
                       x, ss, &als, getamoebalambdaevi, NULL, 50, 1);

    vars_->lambda1 = gsl_vector_get (x, 0);
    vars_->lambda1 = std::pow (10.0, vars_->lambda1);
    gsl_vector_free (x);
    gsl_vector_free (ss);

    a2->~MATRIX();
    mps2->~VECTOR();
    if (!(data_->is_uvdata && 1==data_->transmode))
    {
        lmps2->~VECTOR();
        res2->~VECTOR();
    }
    reg2->~VECTOR();
    MEMORY ps_free( a2dummy );
    MEMORY ps_free( mps2d );
    MEMORY ps_free( lmps2d );
    MEMORY ps_free( res2d );
    MEMORY ps_free( reg2d );
}

void pixsrc_common::findlambda(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if(data_->verbose!=1)
        PRINTER print2screen(data_->print2screenname,
                             "finding regularization strength",
                             cdata_->print2screenmutex);

    if (1==data_->optfinder)
    {
        amoeba_lambda( data_, cdata_, vars_ );
        return;
    }


    double scalar;
    double prevX0[2];
    double prevX1[2];
    double twosmallest[2][2];
    double term[5];
    std::fill(prevX0,prevX0+2,0);
    std::fill(prevX1,prevX1+2,0);
    twosmallest[0][0] = 0;
    twosmallest[1][0] = 0;
    OPERA assign_p_infinity( &twosmallest[0][1] );
    OPERA assign_p_infinity( &twosmallest[1][1] );

    PS_SIT counter=0;
    PS_SIT numLoopsTotal=0;

    MATRIX *dummy1;
    VECTOR *dummy2, *dummy3;
    MEMORY ps_malloc( &dummy1, 1 );
    MEMORY ps_malloc( &dummy2, 1 );
    MEMORY ps_malloc( &dummy3, 1 );

    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->b1->get_nnz() : 0;

    MATRIX *a        = new (dummy1)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 1, data_ );

    VECTOR *smp      = new (dummy2) VECTOR(cdata_, data_, vars_->lonc);
    VECTOR *csmp     = new (dummy3) VECTOR(cdata_, data_, vars_->lonc);

    bool need_to_destroy_a = 1;

    for(;;)
    {
        if( numLoopsTotal==50  ||
            !(OPERA is_finite(prevX0[1]) && OPERA is_finite(prevX0[0]) ) ||
            (counter!=0&&prevX0[0]==prevX1[0]&&prevX0[1]==prevX1[1])  )
        {
            counter=-1;
            OPERA assign_p_infinity( &vars_->lambda1 );

            PRINTER printwarning(data_->print2screenname,
                                 "Giving up on trying to marginalize over lambda. "
                                 "Penalizing model severely.",
                                 cdata_->print2screenmutex);
            break;
        }

        if(counter==1)
            vars_->lambda1=data_->lambdaguess;

        if(counter==0)
        {
            vars_->lambda1=0;
            scalar=-vars_->lonc;
        }
        else
        {
            std::fill(term,term+5,0);

            //if(counter >  1)
            if( !a->are_we_using_cuda && counter >  1 )
            {
                a = new (dummy1)
                    MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 1, data_ );
                need_to_destroy_a = 1;

                smp->set_scalar( 1.0 );
                csmp->zeromeout();
            }

            vars_->b1->plus(vars_->c1, a, 0, 0, 1, vars_->lambda1, cdata_->numthreads );

            /*
              void *clock_ptr;
              double elapsed_time;
              CUDA ps_start_clock( &clock_ptr );
            */
            /*
              void *clock_ptr;
              double elapsed_time;
              CUDA ps_start_clock( &clock_ptr );
            */
            bool zeroit = (counter==1) ? 1 : 0;

            if(trieqon&&a->are_we_using_cuda)
            {
#ifdef __USE_PIXSRC_CUDA__

                CUDA tri_eq_solve( a, smp, vars_->lod, zeroit );
#endif
            }
            else
            {
                a->linequationsolve( smp, vars_->lod, zeroit );
            }
            /*
              CUDA ps_stop_clock( clock_ptr, &elapsed_time );
              string et = "pixsrc took " + OPERA tostring(elapsed_time) + " ms to lineqsolve\n";
              //PRINTER print2screen("", et, cdata_->print2screenmutex);
              */

            /*
              for( PS_SIT j=0; j<smp->size; ++j )
              smp->vec[j] *= 1.0 / vars_->datavariance;
            */
            if (!data_->uv_newvariance)
                smp->mult( 1.0 / vars_->datavariance );
            vars_->c1->mult( smp, csmp, 0, cdata_->numthreads );

            term[0] = vars_->lambda1*csmp->innerproduct( smp );

            /*
              CUDA ps_start_clock( &clock_ptr );
            */

            term[3] = vars_->lambda1 *
                a->trace_invA_B( vars_->c1, cdata_->numthreads, &vars_->fatalerror );

            /*
              CUDA ps_stop_clock( clock_ptr, &elapsed_time );
              et = "pixsrc took " + OPERA tostring(elapsed_time) + " ms to trace_invAB\n";
              //PRINTER print2screen("", et, cdata_->print2screenmutex);
              */


            term[4] = -vars_->lonc;

            /*
              for(PS_SIT p=0; p<3; ++p)
              std::cout << term[p] << " ";
              for(PS_SIT p=3; p<4; ++p)
              std::cout << term[p] << " ";
              for(PS_SIT p=4; p<5; ++p)
              std::cout << term[p] << " ";std::cout << std::endl;
            */

            scalar = 0;
            for( PS_SIT p=0; p<5; ++p )
                scalar += term[p];

            //if( counter >= 1 )
            if( !a->are_we_using_cuda && counter >= 1 )
            {
                a->~MATRIX();
                need_to_destroy_a = 0;
            }
        }

        if(counter==1 && scalar==prevX0[1])
        {
            vars_->lambda1 *= 10;
            PRINTER printwarning(data_->print2screenname,
                                 "bad initial guess for lambda. Most likely too small. "
                                 "Increasing lambda by x10.",
                                 cdata_->print2screenmutex);
            numLoopsTotal++;
            continue;
        }

        prevX1[0] = prevX0[0];
        prevX1[1] = prevX0[1];
        prevX0[0] = vars_->lambda1;
        prevX0[1] = scalar;

        PS_SIT max=0;
        double zer=std::fabs(twosmallest[0][1]);
        double one=std::fabs(twosmallest[1][1]);
        double tmp=std::fabs(scalar);
        if(one>zer) max=1;
        if(tmp<std::min(zer,one))
        {
            twosmallest[max][0]=vars_->lambda1;
            twosmallest[max][1]=scalar;
        }

        if(data_->verbose==3)
        {
            PRINTER print2screen(data_->print2screenname,
                                 "trying regularization " + OPERA tostring(prevX0[0]) +
                                 " => " + OPERA tostring(prevX0[1]),
                                 cdata_->print2screenmutex);
        }

        if(counter>0)
        {
            if(std::fabs((prevX0[0]-prevX1[0])/((prevX0[0]+prevX1[0])/2)) <=
               data_->regftol
               || (counter>1 && std::fabs(prevX0[1])<1 && std::fabs((prevX0[1]-prevX1[1]) / ((prevX0[1]+prevX1[1])/2)) <= data_->regftol)
               || (std::fabs(prevX0[1])<1 && std::fabs(prevX0[1]/vars_->lonc)< data_->regftol))
            {
                /*
                  std::copy( smp->vec, smp->vec + vars_->lonc, vars_->mps->vec );
                  for( PS_SIT j=0; j<vars_->lonc; ++j )
                  vars_->mps->vec[j] *= vars_->datavariance;
                */
                break;
            }

            double add2lambda;

            // following method uses Taylor series expansion
            double firstDeriv = (prevX0[1]-prevX1[1])/(prevX0[0]-prevX1[0]);
            add2lambda = -prevX0[1]/firstDeriv;

            // following method uses logarithmic Taylor series expansion
            if( add2lambda < 0 )
            {
                double logder = ( prevX0[1]-prevX1[1] ) / ( std::exp(prevX0[0])-std::exp(prevX1[0]) );
                add2lambda = -prevX0[1]/(std::exp(prevX0[0])*logder);
            }

            if(prevX0[1]>0)
                add2lambda = -0.9*prevX0[0];
            if(prevX1[1]>0)
                add2lambda = -0.9*prevX0[0];

            vars_->lambda1 += add2lambda;

            // following method uses exponential Taylor series expansion
            // not working right now
            //double expder = ( prevX0[1]-prevX1[1] ) / ( std::log(prevX0[0])-std::log(prevX1[0]) );
            //vars_->lambda1 += -prevX0[1]/(expder/prevX0[0]);
        }
        numLoopsTotal++;
        counter++;

    }

    if (need_to_destroy_a)
        a    -> ~MATRIX ();
    smp      -> ~VECTOR ();
    csmp     -> ~VECTOR ();

    MEMORY ps_free (dummy1);
    MEMORY ps_free (dummy2);
    MEMORY ps_free (dummy3);

    if(counter>1)
        vars_->lambda1=OPERA xintercept( twosmallest[0][0], twosmallest[0][1],
                                         twosmallest[1][0], twosmallest[1][1] );

    if(data_->verbose>=2)
        PRINTER print2screen(data_->print2screenname,
                             "found regularization strength = " +
                             OPERA tostring(vars_->lambda1),
                             cdata_->print2screenmutex);

    if( !OPERA is_finite(vars_->lambda1) || vars_->lambda1<=0 )
        vars_->fatalerror+=1;
}

struct lowmem_convolve_struct
{
    void *dav;
    PS_SIT start_number;
    PS_SIT r0;
    double *blurringoperator_row;
    double *row_entries;
};
void *lowmem_convolve (void *args)
{
    lowmem_convolve_struct *lmcs = (lowmem_convolve_struct*)args;
    dataandvars *dav = (dataandvars*)lmcs->dav;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    for (PS_SIT r1=lmcs->start_number; r1<vars_->lonr; r1+=cdata_->numthreads)
    {
        lmcs->row_entries[r1] = 0;
        // loop over elemens of psf and lensing operator to perform covolution
        for (PS_SIT r2=0; r2<vars_->lonr; ++r2)
            lmcs->row_entries[r1] += lmcs->blurringoperator_row[r2] * vars_->lensingoperatornobo->get (lmcs->r0,r1);
    }
    return NULL;
}
void pixsrc_common::blurit(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if(data_->debug)
        PRINTER print( vars_->tracker, cdata_->basename, data_->name, true,
                       "lensingoperatornobo.MATRIX", vars_->lonr,
                       vars_->lonc, vars_->lensingoperatornobo, 0);

    if (!data_->nopsf)
    {
        if(data_->extlengths[9]<=1 && data_->extlengths[10]<=1)
            data_->nopsf = 1;

        // if we've already created blurring operator
        if(!data_->nopsf && vars_->blurringoperator && vars_->lensingoperator &&
           vars_->blurringoperatorptr && vars_->lensingoperatorptr && !data_->lowmem)
        {
            if(data_->verbose==3)
                PRINTER print2screen(data_->print2screenname, "convolving LO with PSF",
                                     cdata_->print2screenmutex);

            vars_->blurringoperator->mult( vars_->lensingoperatornobo,
                                           vars_->lensingoperator, 0, 0, cdata_->numthreads );

            if(data_->verbose==3)
                PRINTER print2screen(data_->print2screenname, "done colvolving LO with PSF",
                                     cdata_->print2screenmutex);
        }
        else if(!data_->nopsf && !data_->lowmem)
        {
            PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonr*5 : 0;
            MEMORY ps_malloc( &vars_->blurringoperatorptr, 1 );
            vars_->blurringoperator = new (vars_->blurringoperatorptr)
                MATRIX( cdata_, data_, vars_->lonr, vars_->lonr, initsize, 0, data_ );

            // now construct blurring operator from blurrer
            PS_SIT halfSpan = (data_->extlengths[9]-1)/2;
            for(PS_SIT r0=0; r0<vars_->lonr; r0++)
            {
                PS_SIT r0Back = vars_->r4rback[r0];
                PS_SIT x0Back = data_->oldloc[r0Back*2  ];
                PS_SIT y0Back = data_->oldloc[r0Back*2+1];

                for(PS_SIT x=-halfSpan; x<=halfSpan; x++)
                {
                    PS_SIT newX = x0Back+x;
                    for( PS_SIT y=halfSpan; y>=-halfSpan; --y )
                    {
                        if(!data_->psf[x+halfSpan][y+halfSpan])
                            continue;
                        PS_SIT newY = y0Back-y;
                        if(newX>=0 && newY>=0 && newX<=data_->imgx-1 && newY<=data_->imgy-1)
                        {
                            PS_SIT newR = newX*data_->imgy + newY;
                            if(vars_->r4r[newR]!=-1)
                                if(data_->psf[x+halfSpan][y+halfSpan])
                                {
                                    vars_->blurringoperator->set( r0, vars_->r4r[newR],
                                                                  data_->psf[x+halfSpan][y+halfSpan] );
                                }
                        }
                    }
                }
            }

            if(data_->verbose==3)
                PRINTER print2screen(data_->print2screenname, "convolving LO with PSF",
                                     cdata_->print2screenmutex);

            initsize = ( !data_->numgpu2use ) ? vars_->lensingoperatornobo->get_nnz() : 0;
            MEMORY ps_malloc( &(vars_->lensingoperatorptr), 1 );
            vars_->lensingoperator = new (vars_->lensingoperatorptr)
                MATRIX( cdata_, data_, vars_->lonr, vars_->lonc, initsize, 0, data_ );

            vars_->blurringoperator->mult( vars_->lensingoperatornobo,
                                           vars_->lensingoperator, 0, 0, cdata_->numthreads );

            if(data_->verbose==3)
                PRINTER print2screen(data_->print2screenname, "done colvolving LO with PSF",
                                     cdata_->print2screenmutex);
        }
        else if(!data_->nopsf && data_->lowmem)
        {
            if(data_->verbose==3)
                PRINTER print2screen(data_->print2screenname, "convolving LO with PSF (lowmem on)",
                                     cdata_->print2screenmutex);

            PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lensingoperatornobo->get_nnz() : 0;
            MEMORY ps_malloc( &(vars_->lensingoperatorptr), 1 );
            vars_->lensingoperator = new (vars_->lensingoperatorptr)
                MATRIX( cdata_, data_, vars_->lonr, vars_->lonc, initsize, 0, data_ );

            // convolving slow (but low memory) way
            PS_SIT halfSpan = (data_->extlengths[9]-1)/2;
            double *blurringoperator_row, *row_entries;
            MEMORY ps_malloc (&blurringoperator_row, vars_->lonr);
            MEMORY ps_malloc (&row_entries, vars_->lonr);
            // loop over rows of psf operator (or lensing operator)
            for(PS_SIT r0=0; r0<vars_->lonr; r0++)
            {
                std::fill (blurringoperator_row, blurringoperator_row+vars_->lonr, 0);
                PS_SIT r0Back = vars_->r4rback[r0];
                PS_SIT x0Back = data_->oldloc[r0Back*2  ];
                PS_SIT y0Back = data_->oldloc[r0Back*2+1];

                for(PS_SIT x=-halfSpan; x<=halfSpan; x++)
                {
                    PS_SIT newX = x0Back+x;
                    for( PS_SIT y=halfSpan; y>=-halfSpan; --y )
                    {
                        if(!data_->psf[x+halfSpan][y+halfSpan])
                            continue;
                        PS_SIT newY = y0Back-y;
                        if(newX>=0 && newY>=0 && newX<=data_->imgx-1 && newY<=data_->imgy-1)
                        {
                            PS_SIT newR = newX*data_->imgy + newY;
                            if(vars_->r4r[newR]!=-1)
                                if(data_->psf[x+halfSpan][y+halfSpan])
                                {
                                    blurringoperator_row[vars_->r4r[newR]] = data_->psf[x+halfSpan][y+halfSpan];
                                }
                        }
                    }
                }
                // setup and launch threads
                pthread_t              threadshere[cdata_->numthreads];
                lowmem_convolve_struct lmcsstructs[cdata_->numthreads];
                for (PS_SIT proc=0; proc<cdata_->numthreads; ++proc)
                {
                    lmcsstructs[proc].dav = vars_->dav;
                    lmcsstructs[proc].start_number = proc;
                    lmcsstructs[proc].blurringoperator_row = blurringoperator_row;
                    lmcsstructs[proc].row_entries = row_entries;
                    lmcsstructs[proc].r0 = r0;
                    pthread_create(&threadshere[proc] ,cdata_->attrjoinable ,
                                   lowmem_convolve, &lmcsstructs[proc]);
                }
                // wait for threads to finish
                for(PS_SIT proc=0; proc<cdata_->numthreads; ++proc)
                    pthread_join (threadshere[proc],NULL);

                // set entries
                for (PS_SIT r1=0; r1<vars_->lonr; ++r1)
                    vars_->lensingoperator->set (r0,r1,row_entries[r1]);
            }
            MEMORY ps_free (blurringoperator_row);
            MEMORY ps_free (row_entries);

            if(data_->verbose==3)
                PRINTER print2screen(data_->print2screenname, "done colvolving LO with PSF (lowmem on)",
                                     cdata_->print2screenmutex);
        }
    }
    else
    {
        vars_->lensingoperator = vars_->lensingoperatornobo;
    }

    if (!data_->nopsf && data_->debug)
        PRINTER print( vars_->tracker, cdata_->basename, data_->name, 1,
                       "blurringoperator.MATRIX", vars_->lonr,
                       vars_->lonr, vars_->blurringoperator, 0);


    if( data_->printvec && !data_->nopsf )
    {
        PS_SIT len = data_->extlengths[9]*data_->extlengths[10];
        double **vec;
        MEMORY ps_malloc( &vec, len, 3 );

        for(PS_SIT x=0; x<len; x++)
        {
            vec[x][0] = x/data_->extlengths[10];
            vec[x][1] = x%data_->extlengths[10];
            vec[x][2] = data_->psf[(PS_SIT)vec[x][0]][data_->extlengths[10]-1-(PS_SIT)vec[x][1]];
        }

        if( data_->printvec==1 || data_->printvec==3 )
        {
            PRINTER print <double> ( vars_->tracker, cdata_->basename,
                                     data_->name, true,"psf.VECTOR",vec, len, 3, 0, data_->precision, NULL);
        }

        if( data_->printvec==2 || data_->printvec==3 )
        {
            PS_SIT    *dummy;
            double *temparray;
            MEMORY ps_malloc( &temparray, len );
            MEMORY ps_malloc(   &dummy  , len );

            for(PS_SIT x=0; x<len; x++)
            {
                temparray[x]=vec[x][2];
                dummy[x]=x;
            }

            PRINTER printfitsimgplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,"psf.fits",
                                       temparray, data_->wcsinfo,data_->extlengths[9],
                                       data_->extlengths[10],1,dummy,
                                       cdata_->print2screenmutex, 0, 0 );

            MEMORY ps_free( temparray );
            MEMORY ps_free(   dummy   );
        }

        MEMORY ps_free( vec, len );
    }

    if(data_->debug)
        PRINTER print( vars_->tracker, cdata_->basename, data_->name, 1,
                       "lensingoperator.MATRIX", vars_->lonr, vars_->lonc,
                       vars_->lensingoperator, 0 );
}

void pixsrc_common::printimageplanebody(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->printvec )
    {
        // setup for printing out unblurred image
        VECTOR *dummy;
        MEMORY ps_malloc( &dummy, 1 );
        VECTOR *lensedmpsnobo = new (dummy) VECTOR(cdata_, data_, vars_->lonr);
        COMMON lensgalaxy( data_, cdata_, vars_, vars_->lensingoperatornobo,
                           vars_->mps, lensedmpsnobo);

        if( data_->printvec==1 || data_->printvec==3 )
        {
            vars_->data->update_cpu();
            vars_->lensedmps->update_cpu();
            vars_->residual->update_cpu();
            lensedmpsnobo->update_cpu();

            double **vec;
            MEMORY ps_malloc( &vec, vars_->lonr, 3 );

            for(PS_SIT x=0; x<vars_->lonr; x++)
            {
                vec[x][0] = data_->oldloc[vars_->r4rback[x]*2  ];
                vec[x][1] = data_->oldloc[vars_->r4rback[x]*2+1];
            }
            for(PS_SIT x=0; x<vars_->lonr; x++)
                vec[x][2] = vars_->data->get(x);
            PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                     "data.VECTOR", vec, vars_->lonr, 3, 0, data_->precision, NULL);
            for(PS_SIT x=0; x<vars_->lonr; x++)
                vec[x][2] = vars_->lensedmps->get(x);
            PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                     "lensedmps.VECTOR", vec, vars_->lonr, 3, 0, data_->precision, NULL);
            for(PS_SIT x=0; x<vars_->lonr; x++)
                vec[x][2] = lensedmpsnobo->get(x);
            PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                     "lensedmpsnobo.VECTOR", vec, vars_->lonr, 3, 0, data_->precision, NULL);
            for(PS_SIT x=0; x<vars_->lonr; x++)
                vec[x][2] = vars_->residual->get(x);
            PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                     "residuals.VECTOR", vec, vars_->lonr, 3, 0, data_->precision, NULL);

            if( data_->interperr )
            {
                vars_->interperr->update_cpu();

                PS_SIT dim = 1;
                PS_SIT waitlist[1] = {3};
                OPERA pthreadswait(vars_,dim,waitlist);

                for(PS_SIT x=0; x<vars_->lonr; x++)
                    vec[x][2] = vars_->interperr->get(x);
                PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                         "interperr.VECTOR", vec, vars_->lonr, 3, 0, data_->precision, NULL);
            }

            MEMORY  ps_free( vec, vars_->lonr );
        }
        if( data_->printvec==2 || data_->printvec==3 )
        {
            const double *temparray;

            vars_->data     ->dissolve_scalar();
            vars_->lensedmps->dissolve_scalar();
            lensedmpsnobo   ->dissolve_scalar();
            vars_->residual ->dissolve_scalar();

            temparray = vars_->data->get_vec_ptr();
            PRINTER printfitsimgplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,"data.fits",
                                       temparray,data_->wcsinfo,data_->imgx,data_->imgy,1,
                                       vars_->r4r, cdata_->print2screenmutex, 1, 0 );
            temparray = vars_->lensedmps->get_vec_ptr();
            PRINTER printfitsimgplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,"lensedmps.fits",
                                       temparray,data_->wcsinfo,data_->imgx,data_->imgy,1,
                                       vars_->r4r, cdata_->print2screenmutex, 1, 0 );
            temparray = lensedmpsnobo->get_vec_ptr();
            PRINTER printfitsimgplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,"lensedmpsnobo.fits",
                                       temparray,data_->wcsinfo,data_->imgx,data_->imgy,1,
                                       vars_->r4r, cdata_->print2screenmutex, 1, 0 );
            temparray = vars_->residual->get_vec_ptr();
            PRINTER printfitsimgplane( vars_->tracker, cdata_->basename, data_->name,
                                       data_->print2screenname, 1,"residuals.fits",
                                       temparray,data_->wcsinfo,data_->imgx,data_->imgy,1,
                                       vars_->r4r, cdata_->print2screenmutex, 1, 0 );

            if( data_->interperr )
            {
                PS_SIT dim = 1;
                PS_SIT waitlist[1] = {3};
                OPERA pthreadswait(vars_,dim,waitlist);

                vars_->interperr->dissolve_scalar();
                temparray = vars_->interperr->get_vec_ptr();

                PRINTER printfitsimgplane( vars_->tracker, cdata_->basename, data_->name,
                                           data_->print2screenname, 1, "interperr.fits",
                                           temparray,data_->wcsinfo, data_->imgx,data_->imgy,
                                           1, vars_->r4r, cdata_->print2screenmutex, 1, 0 );
            }
        }

        lensedmpsnobo->~VECTOR();
        MEMORY ps_free( dummy );
    }
}

void pixsrc_common::printuvresiduals (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if (!data_->uv_model_filename)
        return;

    if (data_->printvec==2 || data_->printvec==3)
    {
        // read file
        char *fname;
        const char *listcc[2];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = data_->uv_model_filename;
        OPERA concatenate (listcc, 2, &fname);
        // get number of points
        double numuvpts_dbl;
        PS_SIT numuvpts;
        PRINTER readbinary (fname, &numuvpts_dbl, 1, 0);
        numuvpts = (PS_SIT)numuvpts_dbl;
        double *pos;
        MEMORY ps_malloc (&pos, numuvpts*2+1);
        PRINTER readbinary (fname, pos, numuvpts*2+1, 0);
        MEMORY ps_free (fname);
        // copy uv positions
        double **vec;
        MEMORY ps_malloc (&vec, numuvpts, 4);
        for(PS_SIT x=0; x<numuvpts; ++x)
        {
            vec[x][0] = pos[1+numuvpts*0+x];
            vec[x][1] = pos[1+numuvpts*1+x];
        }
        MEMORY ps_free (pos);

        // get spatial decomposition of image plane using RBFs
        PS_SIT num_ndp=0, numextraRBF=0;
        TPS *imgtps = NULL;
        INIT uvdata_get_rbf_mat (data_, cdata_, &imgtps, &num_ndp, 0);

        // decompose image plane
        VECTOR decompose (cdata_, data_, num_ndp);
        imgtps->tpsmat->mult (vars_->lensedmps, &decompose, 0, cdata_->numthreads);
        MEMORY ps_free (imgtps->tpsmat->mat_dense);
        imgtps->tpsmat->mat_dense = NULL;

        double thispos[2];
        MATRIX weights   (cdata_, data_, num_ndp+numextraRBF, 2, -1, -1, data_);
        VECTOR ft        (cdata_, data_, 2);
        // get Fourier transform and residuals
        PS_SIT numblocks = std::ceil(numuvpts/1000.0);
        PS_SIT uvind = 0;
        for (PS_SIT uv=0; uv<numuvpts; ++uv)
        {
            if (1!=data_->verbose)
                if (uv/1000>=uvind)
                {
                    uvind +=1;
                    PRINTER print2screen(data_->print2screenname,
                                         "processing uv block " + OPERA tostring(uvind) +
                                         " of " + OPERA tostring(numblocks),
                                         cdata_->print2screenmutex);
                }

            thispos[0] = vec[uv][0];
            thispos[1] = vec[uv][1];
            imgtps->get_fourier_weights (thispos, &weights, 1);
            weights.mult (&decompose, &ft, 1, cdata_->numthreads);
            vec[uv][2] = ft.get (0);
            vec[uv][3] = ft.get (1);
        }

        PRINTER print <double> (vars_->tracker, cdata_->basename, data_->name, 1,
                                "uvmodel.VECTOR", vec, numuvpts, 4, 0, data_->precision, NULL);

        // cleanup
        delete imgtps;
        MEMORY  ps_free (vec, numuvpts);
    }
}
void pixsrc_common::printuvplane (inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                  VECTOR *data, VECTOR *model, VECTOR *residuals)
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->printvec==2 || data_->printvec==3 )
    {
        PS_SIT sh_ndp = data_->uv_ndp;
        data->update_cpu();
        model->update_cpu();
        residuals->update_cpu();

        double **vec;
        MEMORY ps_malloc (&vec, sh_ndp, 4);

        for(PS_SIT x=0; x<sh_ndp; x++)
        {
            vec[x][0] = data_->uv_oldloc[x*2];
            vec[x][1] = data_->uv_oldloc[x*2+1];
        }
        for(PS_SIT x=0; x<sh_ndp; x++)
        {
            vec[x][2] = data->get (x*2);
            vec[x][3] = data->get (x*2+1);
        }
        PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                 "uvdata.VECTOR", vec, sh_ndp, 4, 0, data_->precision, NULL);
        for(PS_SIT x=0; x<sh_ndp; x++)
        {
            vec[x][2] = model->get (x*2);
            vec[x][3] = model->get (x*2+1);
        }
        PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                 "uvmodel.VECTOR", vec, sh_ndp, 4, 0, data_->precision, NULL);
        for(PS_SIT x=0; x<sh_ndp; x++)
        {
            vec[x][2] = residuals->get (x*2);
            vec[x][3] = residuals->get (x*2+1);
        }
        PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                 "uvresiduals.VECTOR", vec, sh_ndp, 4, 0, data_->precision, NULL);
        for(PS_SIT x=0; x<sh_ndp; x++)
        {
            vec[x][2] = residuals->get (x*2)   / std::sqrt(data_->uv_newvariance[x]);
            vec[x][3] = residuals->get (x*2+1) / std::sqrt(data_->uv_newvariance[x]);
        }
        PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                 "uvresiduals2noise.VECTOR", vec, sh_ndp, 4, 0, data_->precision, NULL);

        MEMORY  ps_free (vec, sh_ndp);
    }
}

void* pixsrc_common::printimageplane(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[8] = 1;

    COMMON printimageplanebody( data_, cdata_, vars_ );

    vars_->pthreadstracker[8] = 2;

    return NULL;
}

void pixsrc_common::findsisterimages(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if(vars_->fatalerror || !data_->findsisterimages)
    {
        return;
    }

    PS_SIT dim = 1;
    PS_SIT waitlist[1] = {7};
    OPERA pthreadswait(vars_,dim,waitlist);

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "finding sister pixels",
                             cdata_->print2screenmutex);

    for(PS_SIT r = 0; r < data_->ndp; r++)
    {
        if(vars_->r4r[r]==-1 && data_->imagemasks[r])
        {
            bool continueit = 1;
            // forcing this statement to make things faster
            if( 1 || data_->subsampling==1 )
            {
                if( vars_->pointer[r]!=-1 )
                {
                    continueit = 0;
                }
            }
            else
            {
                for(PS_SIT x=0; x<data_->subsampling; x++)
                {
                    for(PS_SIT y=0; y<data_->subsampling; y++)
                    {
                        if(vars_->sspointer[r][x*data_->subsampling+y]!=-1)
                        {
                            continueit = 0;
                            break;
                        }
                        if(!continueit)
                            break;
                    }
                    if(!continueit)
                        break;
                }
            }
            if(!continueit)
            {
                vars_->r4r[r] = 0;
            }
        }
    }

    // now we check for those pixels that will get picked up because of the psf (at FWHM only)
    PS_SIT maxx = (data_->extlengths[9] -1)/2;
    PS_SIT maxy = (data_->extlengths[10]-1)/2;
    PS_SIT xr, yr, xnew, ynew, rnew;
    // need the following to only search around pixels selected as of right now.
    PS_SIT *check_this;
    MEMORY ps_malloc (&check_this, data_->ndp);
    std::copy (vars_->r4r, vars_->r4r+data_->ndp, check_this);

    for(PS_SIT r = 0; r < data_->ndp; r++)
    {
        if (check_this[r]!=-1)
        {
            xr = r/data_->imgy;
            yr = r%data_->imgy;

            for(PS_SIT x=0; x<data_->extlengths[9]; x++)
            {
                xnew = xr-maxx+x;
                if( xnew<0 || xnew>data_->imgx-1 )
                    continue;

                for(PS_SIT y=0; y<data_->extlengths[10]; y++)
                {
                    if (!data_->psf[x][y]) // || data_->psf[x][y]<data_->fwhm_weight)
                        continue;

                    ynew = yr+maxy-y;
                    if (ynew<0 || ynew>data_->imgy-1)
                        continue;

                    rnew = xnew*data_->imgy+ynew;
                    if (!data_->imagemasks[rnew] /* || vars_->r4r[rnew]!=-1*/ )
                        continue;

                    vars_->r4r[rnew] = 0;
                }
            }
        }
    }
    MEMORY ps_free (check_this);

    PS_SIT unused=0;
    vars_->lonr = 0;
    for(PS_SIT r = 0; r < data_->ndp; r++)
    {
        if(vars_->r4r[r]!=-1)
        {
            vars_->r4r[r] = r - unused;
            vars_->r4rback[vars_->lonr++] = r;
        }
        else
        {
            unused++;
        }
    }

    // have to re-subsample now that we've picked up some pixels
    if( data_->subsampling > 1 )
        resubsample( data_, cdata_, vars_ );

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "done finding sister pixels",
                             cdata_->print2screenmutex);

    // re-subsample for interpolation errors
    COMMON resubsample_ie( vars_->dav );
}
void pixsrc_common::creater4r(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        return;
    }

    vars_->lonr = 0;
    PS_SIT unused=0;
    for(PS_SIT r = 0; r < data_->ndp; r++)
    {
        if(vars_->r4r[r]!=-1)
        {
            vars_->r4r[r] = r - unused;
            vars_->r4rback[vars_->lonr++] = r;
        }
        else
        {
            ++unused;
        }
    }
}
void pixsrc_common::resubsample(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    // find where subsamples pixels lie in source plane
    // get mins and maxs first to speed up process
    double minxlocal, minylocal, maxxlocal, maxylocal;
    OPERA assign_p_infinity( &minxlocal );
    OPERA assign_p_infinity( &minylocal );
    OPERA assign_n_infinity( &maxxlocal );
    OPERA assign_n_infinity( &maxylocal );
    double pos[2][3];

    for(PS_SIT g=0; g<vars_->lonc; g++)
    {
        if(vars_->triout->pointlist[g*2]<minxlocal)
            minxlocal = vars_->triout->pointlist[g*2];
        if(vars_->triout->pointlist[g*2+1]<minylocal)
            minylocal = vars_->triout->pointlist[g*2+1];
        if(vars_->triout->pointlist[g*2]>maxxlocal)
            maxxlocal = vars_->triout->pointlist[g*2];
        if(vars_->triout->pointlist[g*2+1]>maxylocal)
            maxylocal = vars_->triout->pointlist[g*2+1];
    }

    minxlocal -= CONSTANT smallnumber;
    minylocal -= CONSTANT smallnumber;
    maxxlocal += CONSTANT smallnumber;
    maxylocal += CONSTANT smallnumber;

    PS_SIT numvert = (PS_SIT)vars_->convexhull[0];
    double *convexhull = vars_->convexhull+1;

    {
        double pos2[2];
        double ig, defx, defy, magxx, magxy, magyx, magyy;
        double ra0,dec0;
        pthread_mutex_lock( cdata_->potdefmagmutex );
        pthread_mutex_lock( cdata_->wcsmutex       );
        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "oversampling sister pixels",
                                 cdata_->print2screenmutex);
        double xx,yy;
        double spacing = 1.0/data_->subsampling;
        for(PS_SIT r=0; r<data_->ndp; r++)
        {
            if( data_->imagemasks[r]!=0 && vars_->r4r[r]!=-1 && !vars_->sspointer[r] )
            {
                MEMORY ps_malloc (&vars_->newloc_sslo[r], data_->subsampling*data_->subsampling*2);
                MEMORY ps_malloc (&vars_->sspointer[r],   data_->subsampling*data_->subsampling);

                for(PS_SIT x=0; x<data_->subsampling; x++)
                {
                    std::fill(&vars_->sspointer[r][x*data_->subsampling],
                              &vars_->sspointer[r][(x+1)*data_->subsampling],-1);

                    xx = data_->oldloc[r*2] - 0.5 + spacing*(0.5+x);
                    for(PS_SIT y=0; y<data_->subsampling; y++)
                    {
                        yy = data_->oldloc[r*2+1] - 0.5 + spacing*(0.5+y);

                        if(!strcmp(cdata_->coordsysorig,"PIXEL"))
                            HEADER getlfroms(xx, yy, data_->px,data_->py,data_->r1,data_->r2,pos2);
                        else
                        {
                            HEADER getimgwcscoord(data_->wcs, data_->imgy,xx, yy,&ra0,&dec0);
                            HEADER getwcslfromwcss(ra0,dec0,data_->pra,data_->pdec,data_->r1,data_->r2,pos2);
                        }

#ifdef PS_HAVE_GRAVLENS
                        EXTERNAL ps_potdefmag(pos2[0], pos2[1], &ig, &defx, &defy, &magxx, &magyy, &magxy, &magyx);
#endif
#ifdef PS_HAVE_TRIAXIAL
                        EXTERNAL ps_potdefmag(pos2[0], pos2[1], &ig, &defx, &defy, &magxx, &magyy, &magxy, &magyx, cdata_->tlmparms, cdata_->tlmtime, cdata_->tlmenvirogals);
#endif

                        if(!strcmp(cdata_->coordsysorig,"PIXEL"))
                        {
                            vars_->newloc_sslo[r][(x*data_->subsampling+y)*2]   = xx - defx;
                            vars_->newloc_sslo[r][(x*data_->subsampling+y)*2+1] = yy + defy;
                        }
                        else
                        {
                            HEADER getwcssfromwcsl(-defx,-defy,ra0,dec0,0,0,pos2);
                            HEADER getimgpixcoord(data_->wcs, data_->imgy,cdata_,pos2[0],pos2[1],
                                                  &vars_->newloc_sslo[r][(x*data_->subsampling+y)*2],
                                                  &vars_->newloc_sslo[r][(x*data_->subsampling+y)*2+1]);
                        }
                    }
                }
            }
        }
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );

        // find out which pixels don't have a chance
        PS_SIT *hasachance;
        PS_SIT *hasachancecopy;
        MEMORY ps_malloc( &(hasachance    ), data_->ndp );
        MEMORY ps_malloc( &(hasachancecopy), data_->ndp );
        std::fill(hasachance,hasachance+data_->ndp,1);
        for(PS_SIT r=0; r<data_->ndp; r++)
        {
            if( data_->imagemasks[r]!=0 && vars_->r4r[r]!=-1 && !vars_->sspointer[r] )
            {
                if(vars_->newloc_sslo[r][0]>=maxxlocal ||
                   vars_->newloc_sslo[r][1]>=maxylocal ||
                   vars_->newloc_sslo[r][data_->subsampling*data_->subsampling*2-2]<=minxlocal ||
                   vars_->newloc_sslo[r][data_->subsampling*data_->subsampling*2-1]<=minylocal  )
                {
                    hasachance[r]=0;
                    continue;
                }
            }
            else
                hasachance[r]=0;
        }


        // find out which subpixels have a chance
        char ***hasachancess;
        MEMORY ps_malloc( &(hasachancess), data_->ndp );
        for(PS_SIT r=0; r<data_->ndp; r++)
        {
            if( hasachance[r] )
            {
                MEMORY ps_malloc( &(hasachancess[r]), data_->subsampling, data_->subsampling );
                for(PS_SIT d=0; d<data_->subsampling; d++)
                    std::fill(hasachancess[r][d],hasachancess[r][d]+data_->subsampling,0);

                for(PS_SIT j=0; j<data_->subsampling; j++)
                {
                    for(PS_SIT k=0; k<data_->subsampling; k++)
                    {
                        if (GEOM isinpoly (vars_->newloc_sslo[r][(j*data_->subsampling+k)*2],
                                           vars_->newloc_sslo[r][(j*data_->subsampling+k)*2+1],
                                           convexhull,numvert))
                        {
                            hasachancess[r][j][k]=1;
                        }
                    }
                }
            }
        }


        // actually find in which triangle subsampled pixels lie in source plane
        std::copy( hasachance, hasachance+data_->ndp, hasachancecopy );
        for( PS_SIT r=0; r<data_->ndp; ++r )
            hasachance[r] *= ( data_->subsampling*data_->subsampling );

        for( PS_SIT tri = 0; tri < vars_->triout->numberoftriangles; tri++ )
        {
            for(PS_SIT f=0; f<3; f++)
            {
                pos[0][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2];
                pos[1][f] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+f]*2+1];
            }
            for(PS_SIT r=0; r<data_->ndp; r++)
            {
                if( hasachance[r] )
                {
                    for(PS_SIT x=0; x<data_->subsampling; x++)
                    {
                        for(PS_SIT y=0; y<data_->subsampling; y++)
                        {
                            if(vars_->sspointer[r][x*data_->subsampling+y]==-1 &&
                               hasachancess[r][x][y] &&
                               GEOM isintri(pos,vars_->newloc_sslo[r][(x*data_->subsampling+y)*2],
                                            vars_->newloc_sslo[r][(x*data_->subsampling+y)*2+1]))
                            {
                                vars_->sspointer[r][x*data_->subsampling+y] = tri;
                                --hasachance[r];
                            }
                        }
                    }
                }
            }
        }

        for(PS_SIT r=0; r<data_->ndp; r++)
            if(hasachancecopy[r])
                MEMORY ps_free( hasachancess[r], data_->subsampling );
        MEMORY ps_free( hasachancess );
        MEMORY ps_free( hasachance                         );
        MEMORY ps_free( hasachancecopy                     );
    }
}

void pixsrc_common::printfinalmessages(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (vars_->fatalerror)
    {
        PS_SIT dim = vars_->numpthreads;
        PS_SIT waitlist[vars_->numpthreads];
        for(PS_SIT j=0; j<vars_->numpthreads; ++j)
            waitlist[j] = j;
        OPERA pthreadswaitifstarted(vars_,dim,waitlist);

        if( data_->verbose!=1 )
            PRINTER print2screen(data_->print2screenname,
                                 "fatal error: code " + OPERA tostring(vars_->fatalerror),
                                 cdata_->print2screenmutex);

        if(data_->printdetails)
        {
            PS_SIT dim2 = 3 + 6 + 1 + data_->extlengths[5] + 9 + 3;
            double writeout[dim2];
            std::fill(writeout,writeout+dim2,0);

            PS_SIT startindex=0;
            writeout[startindex+0] = vars_->fatalerror;
            writeout[startindex+1] = vars_->endedearly;
            writeout[startindex+2] = vars_->modelrejected;

            startindex += 3;
            startindex += 6;
            writeout[startindex+0] = -2.0*data_->evidence;
            PRINTER writeoutstream <double> (writeout,dim2,data_->details->stream,data_->details->lock, data_->precision, NULL);
            startindex += 1 + data_->extlengths[5] + 9 + 3;
        }

        return;
    }

    if(!vars_->endedearly)
    {
        PS_SIT dim = 1;
        PS_SIT waitlist[1] = {11};
        OPERA pthreadswait(vars_,dim,waitlist);
    }

    if(data_->printdetails)
    {
        PS_SIT dim = 3 + 6 + 1 + data_->extlengths[5] + 9 + 3;
        double writeout[dim];
        std::fill(writeout,writeout+dim,0);

        PS_SIT startindex=0;
        writeout[startindex+0] = vars_->fatalerror;
        writeout[startindex+1] = vars_->endedearly;
        writeout[startindex+2] = vars_->modelrejected;

        startindex += 3;

        if(!vars_->endedearly)
        {
            if (!data_->use_shapelets)
            {
                OPERA assign_n_infinity( &writeout[startindex+0] );
                PS_SIT maxpos=-1;
                for(PS_SIT j=0; j<vars_->lonc; j++)
                    if(vars_->mps->get(j)>writeout[startindex+0])
                    {
                        maxpos=j;
                        writeout[startindex+0] = vars_->mps->get(j);
                    }
                std::copy( &(vars_->triout->pointlist[maxpos*2]),
                           &(vars_->triout->pointlist[maxpos*2])+2, &(writeout[startindex+1]) );
            }
            //writeout[startindex+1] = vars_->srclocall[vars_->c4cback[maxpos]][0];
            //writeout[startindex+2] = vars_->srclocall[vars_->c4cback[maxpos]][1];

            pthread_mutex_lock( cdata_->wcsmutex );
            HEADER getimgwcscoord(data_->wcs, data_->imgy,writeout[startindex+1],writeout[startindex+2],
                                  &(writeout[startindex+3]),&(writeout[startindex+4]));
            HEADER getwcslfromwcss(writeout[startindex+3], writeout[startindex+4],
                                   data_->pra,data_->pdec,data_->r1,data_->r2,
                                   &(writeout[startindex+1])                            );
            pthread_mutex_unlock( cdata_->wcsmutex );

            writeout[startindex+5] = vars_->magger;
        }

        startindex += 6;

        writeout[startindex+0] = -2.0*data_->evidence;
        std::copy(vars_->penalties,vars_->penalties+data_->extlengths[5],
                  &(writeout[startindex+1]));

        startindex += 1 + data_->extlengths[5];

        if( vars_->terminatelensing )
        {}
        else if( data_->noevinochi )
        {}
        else if( data_->onlychi )
        {
            writeout[startindex+8] =  vars_->ed/(vars_->truelonr-vars_->dofreduction);
        }
        else
        {
            writeout[startindex+0] = -vars_->truelonr/2.0*std::log(vars_->datavariance);
            writeout[startindex+1] = -vars_->logdeta/2.0;
            writeout[startindex+2] =  vars_->logdetc/2.0;
            writeout[startindex+3] = -vars_->ed/2.0;
            writeout[startindex+4] = -vars_->lambda1*vars_->es/2.0;
            writeout[startindex+5] =  vars_->lonc/2.0*std::log(vars_->lambda1);
            writeout[startindex+6] = -vars_->truelonr/2.0*std::log(2*CONSTANT pi);
        }

        startindex += 9;

        writeout[startindex+0] = vars_->lonr;
        writeout[startindex+1] = vars_->lonc;
        writeout[startindex+2] = vars_->lambda1;

        startindex += 3;

        PRINTER writeoutstream <double> (writeout,dim,data_->details->stream,data_->details->lock, data_->precision, NULL);
    }

    if(!vars_->endedearly)
    {
        PS_SIT dim = 2;
        PS_SIT waitlist[2] = {8,9};
        OPERA pthreadswait(vars_,dim,waitlist);
    }
}
