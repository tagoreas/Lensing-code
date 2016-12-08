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



#include "pixsrc_shapelets.hpp"
#include "pixsrc_shapelets_operations.hpp"
#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_common_adaptive.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_triangulation.hpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_printer.hpp"
#include <algorithm>
#include <cmath>


pixsrc_shapelets::~pixsrc_shapelets(){}
pixsrc_shapelets::pixsrc_shapelets(PS_SIT imagenumber, inputdata *datavec, commoninputdata *cdata__, PS_SIT tracker_)
{
    MEMORY ps_malloc( &vars_, 1 );
    vars_->tracker = tracker_;
    vars_->imagenumber = imagenumber;
    vars_->pixsrc_class = (void*)this;
    data_  = &datavec[imagenumber];
    cdata_ = cdata__;

    vars_->num_shapelets1    = data_->num_shapelets1;
    vars_->num_shapelets2    = data_->num_shapelets2;
    vars_->numberofshapelets = data_->numberofshapelets;

    COMMON startlensing(data_, cdata_, vars_);

    if( !data_->doreconstruction )
        vars_->endedearly = 1;

    if( vars_->terminatelensing || vars_->endedearly )
    {
        double pensum[2]={0,0};
        vars_->lambda_lens = data_->lambda_lens<0 ? -data_->lambda_lens : data_->lambda_lens;

        for(PS_SIT j=0; j<data_->extlengths[5]; j++)
            pensum[0] += vars_->penalties[j];

        if (cdata_->npl && vars_->lambda_lens)
            pensum[1] += vars_->lambda_lens * NONPARAMLENS gsl_lens_penalty (data_, cdata_);

        //std::cout << "penalties: " << pensum[0] << " " << pensum[1] << std::endl;

        data_->evidence = -0.5*(pensum[0]+pensum[1]);
    }
    else
    {
        COMMON creater4r(data_, cdata_, vars_);
        createc4c();

        pthread_create(&(vars_->pthreads[0] ), cdata_->attrdetached,
                       COMMONADAPTIVE createc         , vars_->dav);
        pthread_create(&(vars_->pthreads[7] ), cdata_->attrdetached,
                       TRIANGULATION pstriangulate    , vars_->dav);

        COMMON         findsisterimages     (data_, cdata_, vars_);
        COMMONADAPTIVE createlo             (data_, cdata_, vars_);
        COMMON         sourcereconstructions(data_, cdata_, vars_);

        pthread_create(&(vars_->pthreads[8] ), cdata_->attrdetached,
                       COMMON printimageplane         , vars_->dav);
        pthread_create(&(vars_->pthreads[9] ), cdata_->attrdetached,
                       COMMONADAPTIVE printsource     , vars_->dav);
        pthread_create(&(vars_->pthreads[10]), cdata_->attrdetached,
                       COMMONADAPTIVE getmagnification, vars_->dav);
        pthread_create(&(vars_->pthreads[11]), cdata_->attrdetached,
                       STATISTIC evidencecalculator   , vars_->dav);
    }

    COMMON printfinalmessages(data_, cdata_, vars_);
    MEMORY freememory(data_, cdata_, vars_);
}
void pixsrc_shapelets::flagsrcpixels ()
{
    MEMORY ps_malloc( &(vars_->fallsinsrc), data_->ndp );
    for( PS_SIT r=0; r<data_->ndp; ++r )
    {
        //if( vars_->r4r[r] != -1 )
        {
            vars_->fallsinsrc[r] = GEOM fallsinsrcinputcircle(
                vars_->newloc[r*2], vars_->newloc[r*2+1],
                data_->srcinputcircle, data_->extlengths[2] );
        }
    }
}

void pixsrc_shapelets::createc4c ()
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "getting convex hull of ray-traced image pixels",
                             cdata_->print2screenmutex);

    // find out which image pixels are candidate source pixels
    flagsrcpixels();

    // create convex hull of qualifying pixels
    MEMORY ps_malloc( &(vars_->triin->pointlist), data_->ndp*2 );
    vars_->triin->numberofpoints = 0;
    for( PS_SIT r=0; r<data_->ndp; ++r )
    {
        if (vars_->fallsinsrc[r] && -1!=vars_->r4r[r])
        {
            vars_->triin->pointlist[vars_->triin->numberofpoints*2]   = vars_->newloc[r*2];
            vars_->triin->pointlist[vars_->triin->numberofpoints*2+1] = vars_->newloc[r*2+1];
            ++vars_->triin->numberofpoints;
        }
    }

    // don't need whole grid. So, make convex hull the grid.
    PS_SIT numverts;
    double *convexhull;
    GEOM getconvexhull( vars_->triin, vars_->triout, &numverts,
                        &convexhull, data_->triswitches );

    MEMORY triangulatestructdestruct ( vars_->triin    );
    MEMORY triangulatestructdestruct ( vars_->triout   );
    MEMORY ps_free                   ( vars_->triin    );
    MEMORY ps_free                   ( vars_->triout   );
    MEMORY ps_malloc                 (&vars_->triin, 1 );
    MEMORY ps_malloc                 (&vars_->triout, 1);
    MEMORY triangulatestructinit     ( vars_->triin    );
    MEMORY triangulatestructinit     ( vars_->triout   );

    // perturb positions of points (to keep triangulation stable)
    double centroidx, centroidy, x, y, r, cos, sin, ran;
    GEOM getcentroid (numverts, convexhull, &centroidx, &centroidy);
    for( PS_SIT j=0; j<numverts; ++j )
    {
        x = convexhull[j*2]   - centroidx;
        y = convexhull[j*2+1] - centroidy;
        r = std::sqrt (x*x + y*y);
        cos = x/r;
        sin = y/r;
        ran = std::abs(OPERA randomgaussian(cdata_->ps_gsl_ran_r));
        ran = ran<1 ? ran+1 : ran;

        convexhull[j*2]   += ran*cos;
        convexhull[j*2+1] += ran*sin;
    }

    vars_->triin->pointlist = convexhull;
    vars_->triin->numberofpoints = numverts;
    vars_->need2retriangulate = 1;

    vars_->lonc = data_->numberofshapelets;

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "done getting convex hull of ray-traced image pixels",
                             cdata_->print2screenmutex);
}
