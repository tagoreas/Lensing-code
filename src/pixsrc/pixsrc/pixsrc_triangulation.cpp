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



#include "pixsrc_triangulation.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_geometry.hpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_common.hpp"

#include <cstring>

void* pixsrc_triangulation::pstriangulate(void *args)
{
    dataandvars *dav = (dataandvars*)args;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;

    vars_->pthreadstracker[7] = 1;

    pstriangulatebody( data_, cdata_, vars_ );

    vars_->pthreadstracker[7] = 2;

    return NULL;
}

void pixsrc_triangulation::pstriangulatebody( inputdata *data_, commoninputdata *cdata_, lensvar *vars_ )
{
    if(vars_->fatalerror || data_->gridtype==0)
    {
        return;
    }

    if( vars_->need2retriangulate )
    {
        //PRINTER printtristruct(vars_->triin, cdata_->print2screenmutex);
        //PRINTER printtristruct(vars_->triout, cdata_->print2screenmutex);

        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "triangulating grid",
                                 cdata_->print2screenmutex);

        ps_tri_triangulate( data_->triswitches, vars_->triin, vars_->triout,
                            (struct triangulateio*)NULL                      );

        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname,
                                 "done triangulating grid",
                                 cdata_->print2screenmutex);

    }

    // if no Steiner points might have been added into the triangulation
    if( !data_->mintriangleangle )
    {
        vars_->triout->pointlist = vars_->triin->pointlist;
        vars_->triin ->pointlist = NULL;
    }
    else if (!data_->use_shapelets)
    {
        vars_->lonc = vars_->triout->numberofpoints;
    }

    if(data_->verbose!=1)
        PRINTER print2screen(data_->print2screenname,
                             "ray-tracing data onto grid",
                             cdata_->print2screenmutex);

    char *triused;
    MEMORY ps_malloc( &triused, vars_->triout->numberoftriangles );
    std::fill( triused, triused+vars_->triout->numberoftriangles, 0 );

    // find where subsampled pixels lie in source plane
    // get mins and maxs first to speed up process
    double minxlocal, minylocal, maxxlocal, maxylocal;
    OPERA assign_p_infinity( & minxlocal );
    OPERA assign_p_infinity( & minylocal );
    OPERA assign_n_infinity( & maxxlocal );
    OPERA assign_n_infinity( & maxylocal );

    for(PS_SIT g=0; g<vars_->triout->numberofpoints; g++)
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

    // NULLs signal no triangulation to be performed
    GEOM getconvexhull( NULL,vars_->triout,NULL,&(vars_->convexhull),NULL );
    PS_SIT numvert = (PS_SIT)vars_->convexhull[0];
    double *convexhull = vars_->convexhull+1;

    //if( data_->subsampling==1 )
    {
        // find out which pixels don't have a chance
        char *hasachance;
        MEMORY ps_malloc( &(hasachance), data_->ndp );
        std::fill(hasachance,hasachance+data_->ndp,1);

        for(PS_SIT r=0; r<data_->ndp; ++r)
        {
            if( data_->imagemasks[r]!=0 && (vars_->r4r[r]!=-1 || data_->findsisterimages) )
            {
                // if point is outside bounding rectangle
                if(vars_->newloc[r*2  ]>=maxxlocal ||
                   vars_->newloc[r*2+1]>=maxylocal ||
                   vars_->newloc[r*2  ]<=minxlocal ||
                   vars_->newloc[r*2+1]<=minylocal  )
                {
                    hasachance[r]=0;
                    continue;
                }
                // if point is outside convex hull or triangulation
                if(!GEOM isinpolypadded(vars_->newloc[r*2],vars_->newloc[r*2+1],convexhull,numvert))
                {
                    hasachance[r]=0;
                }
            }
            else
            {
                hasachance[r]=0;
            }
        }

        PS_SIT tri_seed = 0;
        MEMORY ps_malloc( &(vars_->pointer), data_->ndp );
        std::fill( vars_->pointer,vars_->pointer+data_->ndp,-1 );
        for(PS_SIT r=0; r<data_->ndp; r++)
            if( hasachance[r] )
            {
                vars_->pointer[r] = GEOM search_triangle( vars_->triout, &tri_seed,
                                                          vars_->newloc[r*2  ],
                                                          vars_->newloc[r*2+1]    );
                if( vars_->pointer[r]!=-1 )
                    triused[vars_->pointer[r]] = 1;
            }

        MEMORY ps_free( hasachance );
    }
    //else
    if( data_->subsampling > 1 )
    {
        // find subsampled pixel positions in source plane
        MEMORY ps_malloc (&vars_->newloc_sslo, data_->ndp);
        MEMORY ps_malloc (&vars_->sspointer,   data_->ndp);

        std::fill( vars_->sspointer, vars_->sspointer+data_->ndp, (PS_SIT*)NULL );
        double magxx, magxy, magyx, magyy, pot;
        if(data_->verbose==3)
            PRINTER print2screen(data_->print2screenname,
                                 "oversampling pixels",
                                 cdata_->print2screenmutex);
        double xx,yy;
        double spacing = 1.0/data_->subsampling;
        pthread_mutex_lock( cdata_->potdefmagmutex );
        pthread_mutex_lock( cdata_->wcsmutex       );
        for(PS_SIT r=0; r<data_->ndp; r++)
        {
            if( data_->imagemasks[r]!=0 && vars_->r4r[r]!=-1 )
            {
                MEMORY ps_malloc (&vars_->newloc_sslo[r],
                                  data_->subsampling*data_->subsampling*2);
                MEMORY ps_malloc (&vars_->sspointer[r],
                                  data_->subsampling*data_->subsampling);

                for(PS_SIT x=0; x<data_->subsampling; x++)
                {
                    std::fill( &vars_->sspointer[r][x*data_->subsampling],
                               &vars_->sspointer[r][(x+1)*data_->subsampling], -1 );

                    xx = data_->oldloc[r*2] - 0.5 + spacing*(0.5+x);
                    for(PS_SIT y=0; y<data_->subsampling; y++)
                    {
                        yy = data_->oldloc[r*2+1] - 0.5 + spacing*(0.5+y);
                        COMMON raytrace (data_, cdata_, xx, yy,
                                         &vars_->newloc_sslo[r][(x*data_->subsampling+y)*2],
                                         &vars_->newloc_sslo[r][(x*data_->subsampling+y)*2+1],
                                         &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
                    }
                }
            }
        }
        pthread_mutex_unlock( cdata_->potdefmagmutex );
        pthread_mutex_unlock( cdata_->wcsmutex       );

        // find out which pixels don't have a chance
        PS_SIT *hasachance;
        MEMORY ps_malloc (&hasachance, data_->ndp);
        std::fill (hasachance, hasachance+data_->ndp, 1);
        for (PS_SIT r=0; r<data_->ndp; r++)
        {
            if (data_->imagemasks[r]!=0 && vars_->r4r[r]!=-1)
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
            {
                hasachance[r]=0;
            }
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
                        if(GEOM isinpolypadded(vars_->newloc_sslo[r][(j*data_->subsampling+k)*2],
                                               vars_->newloc_sslo[r][(j*data_->subsampling+k)*2+1],convexhull,numvert))
                        {
                            hasachancess[r][j][k]=1;
                        }
                    }
                }
            }
        }


        // actually find in which triangle subsampled pixels lie in source plane
        PS_SIT tri_seed = 0;
        for(PS_SIT r=0; r<data_->ndp; r++)
            if( hasachance[r] )
                for(PS_SIT x=0; x<data_->subsampling; x++)
                    for(PS_SIT y=0; y<data_->subsampling; y++)
                        if( hasachancess[r][x][y] )
                        {
                            vars_->sspointer[r][x*data_->subsampling+y] =
                                GEOM search_triangle (vars_->triout, &tri_seed,
                                                      vars_->newloc_sslo[r][(x*data_->subsampling+y)*2],
                                                      vars_->newloc_sslo[r][(x*data_->subsampling+y)*2+1]);
                            if( vars_->sspointer[r][x*data_->subsampling+y]!=-1 )
                                triused[vars_->sspointer[r][x*data_->subsampling+y]] = 1;
                        }


        for(PS_SIT r=0; r<data_->ndp; r++)
            if(hasachance[r])
                MEMORY ps_free( hasachancess[r], data_->subsampling );
        MEMORY ps_free( hasachancess );
        MEMORY ps_free( hasachance                         );
    }

    if (!data_->use_shapelets)
    {
        char *pixelused;
        char *triremove;
        char *edgekeep;
        TRIANGULATION retriangulateholes( data_, cdata_, vars_,
                                          &triused, &pixelused, &triremove, &edgekeep );
        TRIANGULATION trimtrilist       ( data_, cdata_, vars_, triremove             );
        TRIANGULATION trimedgelist      ( data_, cdata_, vars_, edgekeep              );
        TRIANGULATION trimpointlist     ( data_, cdata_, vars_, pixelused             );
        MEMORY ps_free( pixelused );
        MEMORY ps_free( triremove );
        MEMORY ps_free( edgekeep  );
    }
    MEMORY ps_free( triused   );

    if(data_->verbose!=1)
        PRINTER print2screen(data_->print2screenname,
                             "done ray-tracing data onto grid",
                             cdata_->print2screenmutex);

    if(data_->debug)
    {
        vector<string> vec(vars_->triout->numberoftriangles*3+1);
        vec[0] = "image";
        for(PS_SIT t=0; t<vars_->triout->numberoftriangles; ++t)
        {
            PS_SIT vback=2;
            for(PS_SIT v=0; v<3; ++v)
            {
                vec[t*3+v+1] = "line(" +
                    OPERA tostring(vars_->triout->pointlist[vars_->triout->trianglelist[t*3+v    ]*2  ]) + ","
                    + OPERA tostring(vars_->triout->pointlist[vars_->triout->trianglelist[t*3+v    ]*2+1]) + ","
                    + OPERA tostring(vars_->triout->pointlist[vars_->triout->trianglelist[t*3+vback]*2  ]) + ","
                    + OPERA tostring(vars_->triout->pointlist[vars_->triout->trianglelist[t*3+vback]*2+1]) + ")";

                vback = v;
            }
        }
        PRINTER print <string> ( vars_->tracker, cdata_->basename, data_->name, true,"triangulation.reg",vec, 0, data_->precision, NULL);
    }

    return;
}

void pixsrc_triangulation::trimtrilist(inputdata *data_, commoninputdata *cdata_, lensvar *vars_, char *triremove )
{
    // removing unused triangles

    if( !triremove )
        return;

    PS_SIT notorig = vars_->triout->numberoftriangles;
    // this lil piece of code needs to be separate from the for loop that follows
    // b/c not every unused triangle is necessarily moved, and
    // triused changes as the code progresses
    for( PS_SIT t=0; t<notorig; ++t )
    {
        if( triremove[t] )
        {
            --vars_->triout->numberoftriangles;
        }
    }

    for( PS_SIT t=0; t<notorig; ++t )
    {
        // found empty spot
        if( triremove[t] )
        {
            for( PS_SIT tt=notorig-1; tt>t; --tt )
            {
                // found spot to move into empty spot
                if( !triremove[tt] )
                {
                    // move point
                    triremove[t]  = 0;
                    triremove[tt] = 1;
                    for( PS_SIT v=0; v<3; ++v )
                        vars_->triout->trianglelist[t*3+v] = vars_->triout->trianglelist[tt*3+v];

                    // correct entries in pointers
                    for( PS_SIT r=0; r<data_->ndp; ++r )
                        if( vars_->pointer[r]==tt )
                            vars_->pointer[r] = t;
                    if( data_->subsampling>1 )
                    {
                        for( PS_SIT r=0; r<data_->ndp; ++r )
                            if( vars_->sspointer[r] )
                                for( PS_SIT x=0; x<data_->subsampling; ++x )
                                    for( PS_SIT y=0; y<data_->subsampling; ++y )
                                        if( vars_->sspointer[r][x*data_->subsampling+y]==tt )
                                            vars_->sspointer[r][x*data_->subsampling+y] = t;
                    }

                    break;
                }
            }
        }
    }
}

void pixsrc_triangulation::trimpointlist( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, char *pixelused )
{
    // here i just mark which source pixels are actually used in the lensing operator
    // unused pixels are removed
    // this now voids and edgelists .. but that could be corrected

    if( !pixelused )
        return;

    PS_SIT noporig = vars_->triout->numberofpoints;
    // this lil piece of code needs to be separate from the for loop that follows
    // b/c not every unused pixel is necessarily moved, and
    // pixelused changes as the code progresses
    for( PS_SIT p=0; p<noporig; ++p )
    {
        if( !pixelused[p] )
        {
            --vars_->triout->numberofpoints;
            --vars_->lonc;
        }
    }

    for( PS_SIT c=0; c<noporig; ++c )
    {
        // found empty spot
        if( !pixelused[c] )
        {
            for( PS_SIT cc=noporig-1; cc>c; --cc )
            {
                // found spot to move into empty spot
                if( pixelused[cc] )
                {
                    // move point
                    pixelused[c]  = 1;
                    pixelused[cc] = 0;
                    vars_->triout->pointlist[c*2  ] = vars_->triout->pointlist[cc*2  ];
                    vars_->triout->pointlist[c*2+1] = vars_->triout->pointlist[cc*2+1];

                    // correct entries in trianglelist
                    for( PS_SIT t=0; t<vars_->triout->numberoftriangles*3; ++t )
                        if( vars_->triout->trianglelist[t]==cc )
                            vars_->triout->trianglelist[t] = c;

                    // correct entries in edgelist
                    for( PS_SIT e=0; e<vars_->triout->numberofedges*2; ++e )
                        if( vars_->triout->edgelist[e]==cc )
                            vars_->triout->edgelist[e] = c;

                    break;
                }
            }
        }
    }
}

void pixsrc_triangulation::trimedgelist( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, char *edgekeep )
{
    // this can be used to trim the edgelist to only include those edges
    // on the convex hull. Other edges are wrong anyway.
}

void pixsrc_triangulation::retriangulateholes( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, char **triused, char **pixelused, char **triremove, char **edgekeep )
{
    // for the pixels that are unused, i call these 'holes' in the triangulation
    // and i retriangulate these holes without the unused pixel

    if( !data_->rmpix )
    {
        *pixelused  = 0;
        *triremove  = 0;
        *edgekeep = 0;
        return;
    }

    // find unused pixels
    // initializing all pixels as unused
    MEMORY ps_malloc( &(*pixelused) , vars_->triout->numberofpoints );
    std::fill( *pixelused, *pixelused+vars_->triout->numberofpoints, 0 );
    for( PS_SIT tri=0; tri<vars_->triout->numberoftriangles; ++tri )
    {
        if( (*triused)[tri] )
        {
            for( PS_SIT v=0; v<3; ++v )
                // setting pixel status to 'used'
                (*pixelused)[vars_->triout->trianglelist[tri*3+v]] = 1;
        }
    }
    // check if there are any to remove and exit if need be
    PS_SIT num2remove = 0;
    for( PS_SIT p=0; p<vars_->triout->numberofpoints; ++p )
        if( !(*pixelused)[p] )
        {
            ++num2remove;
        }
    if( data_->verbose!=1 )
    {
        PRINTER print2screen(data_->print2screenname,
                             "removing "+OPERA tostring(num2remove) +
                             " unused pixels from triangulation",
                             cdata_->print2screenmutex);
    }
    if( !num2remove )
    {
        MEMORY ps_free( *pixelused );
        *pixelused  = 0;
        *triremove  = 0;
        *edgekeep = 0;
        return;
    }

    // prints out triangulations of filled holes
    std::ofstream *holesfile    = 0;
    std::ofstream *holesfileptr = 0;
    if( data_->debug )
    {
        MEMORY ps_malloc( &holesfileptr, 1 );
        holesfile = new (holesfileptr) std::ofstream;
        PRINTER openoutstream(cdata_->basename, data_->name, 1,
                              OPERA tostring(vars_->imagenumber)+"_triangulationholes.reg",
                              holesfile, 0);

        PRINTER writeoutstream <string> ("image", holesfile, NULL, data_->precision, NULL);
    }

    // this keeps track of triangles to remove at the end
    MEMORY ps_malloc( &(*triremove), vars_->triout->numberoftriangles );
    std::fill( *triremove, *triremove+vars_->triout->numberoftriangles, 0 );

    // this keeps track of edges to keep at the end
    MEMORY ps_malloc( &(*edgekeep), vars_->triout->numberofedges );
    std::copy( vars_->triout->edgemarkerlist, vars_->triout->edgemarkerlist +
               vars_->triout->numberofedges, *edgekeep );

    // need b/c pixelused changes as we go
    char *pixelusedorig;
    MEMORY ps_malloc( &pixelusedorig, vars_->triout->numberofpoints );
    std::copy( *pixelused, *pixelused+vars_->triout->numberofpoints, pixelusedorig );

    // declare
    PS_SIT *trisadded, *edgesadded, *tempint, *newtri2oldtri, *removelist;
    PS_SIT tricapacity, edgecapacity, numtriadded, numedgesadded, removecapacity, numremoveadded;
    PS_SIT v[2]={0,0}, capacity[2], segpos[2];
    double pos[2][2];
    char *tempchar;
    char addsegment;
    char vgood[2];
    struct triangulateio *holein,*holeout;

    tricapacity = edgecapacity = removecapacity = 10;
    MEMORY ps_malloc( &trisadded,  tricapacity*3  );
    MEMORY ps_malloc( &edgesadded, edgecapacity*2 );
    numtriadded = numedgesadded = numremoveadded = segpos[0] = segpos[1] = 0;
    removelist = 0;
    for( PS_SIT c0=0; c0<vars_->triout->numberofpoints; ++c0 )
    {
        // identify pixel to be removed
        if( !(*pixelused)[c0] )
        {
            // initialize pixellist with just this one
            MEMORY ps_malloc( &removelist, removecapacity );
            removelist[0] = c0;
            numremoveadded = 1;

            // initialize
            capacity[0] = capacity[1] = 10;
            MEMORY ps_malloc(            &holein , 1 );
            MEMORY ps_malloc(            &holeout, 1 );
            MEMORY triangulatestructinit( holein     );
            MEMORY triangulatestructinit( holeout    );
            holein->numberofpoints   = 0;
            holein->numberofsegments = 0;
            holein->numberofholes    = 0;
            holein->numberofregions  = 0;
            MEMORY ps_malloc( &(holein->pointlist  ), capacity[0]*2 );
            MEMORY ps_malloc( &(holein->segmentlist), capacity[1]*2 );
            MEMORY ps_malloc( &newtri2oldtri        , capacity[0]   );

            for( PS_SIT pix=0; pix<numremoveadded; ++pix )
            {
                if( (*pixelused)[removelist[pix]] )
                {
                    continue;
                }
                (*pixelused)[removelist[pix]] = 1;

                // scan pointlist and find surrounding pixels
                for( PS_SIT t=0; t<vars_->triout->numberoftriangles*3; ++t )
                {
                    if( vars_->triout->trianglelist[t]==removelist[pix] )
                    {
                        // if pointlist will overflow, reallocate and copy pointlist
                        if( capacity[0]-2<holein->numberofpoints )
                        {
                            capacity[0] *= 2;
                            // fix pointlisth
                            OPERA resize( &(holein->pointlist),
                                          holein->numberofpoints*2, capacity[0]*2 );
                            // fix bookmarking array
                            OPERA resize( &newtri2oldtri,
                                          holein->numberofpoints, capacity[0] );
                        }
                        // if segmentlist will overflow
                        if( capacity[1]-1<holein->numberofsegments )
                        {
                            capacity[1] *= 2;
                            OPERA resize( &(holein->segmentlist),
                                          holein->numberofsegments*2, capacity[1]*2 );
                        }

                        // mark as triangle to be removed
                        (*triremove)[t/3] = 1;

                        // find out which vertex we found and get the other 2
                        switch( t%3 )
                        {
                        case 0:
                            v[0] = vars_->triout->trianglelist[t+1];
                            v[1] = vars_->triout->trianglelist[t+2];
                            break;
                        case 1:
                            v[0] = vars_->triout->trianglelist[t-1];
                            v[1] = vars_->triout->trianglelist[t+1];
                            break;
                        case 2:
                            v[0] = vars_->triout->trianglelist[t-2];
                            v[1] = vars_->triout->trianglelist[t-1];
                            break;
                        default: PRINTER printerror( data_->print2screenname,
                                                     "the impossible happened!",
                                                     cdata_->print2screenmutex );
                            break;
                        }

                        // check if point has already been added to pointlist
                        // or should not be added (if also to be removed)
                        vgood[0] = vgood[1] = addsegment = 1;
                        for( PS_SIT j=0; j<2; ++j )
                        {
                            if( !pixelusedorig[v[j]] )
                            {
                                vgood[j] = 0;
                                addsegment = 0;

                                // need to add newly found pixel to remove list
                                // first check for overflow
                                if( removecapacity-1<numremoveadded )
                                {
                                    removecapacity *= 2;
                                    OPERA resize( &removelist,
                                                  numremoveadded, removecapacity );
                                }
                                removelist[numremoveadded++] = v[j];
                            }
                        }

                        pos[0][0] = vars_->triout->pointlist[v[0]*2  ];
                        pos[0][1] = vars_->triout->pointlist[v[0]*2+1];
                        pos[1][0] = vars_->triout->pointlist[v[1]*2  ];
                        pos[1][1] = vars_->triout->pointlist[v[1]*2+1];
                        for( PS_SIT p=0; p<holein->numberofpoints; ++p )
                        {
                            for( PS_SIT j=0; j<2; ++j)
                            {
                                if( holein->pointlist[p*2  ]==pos[j][0] &&
                                    holein->pointlist[p*2+1]==pos[j][1] )
                                {
                                    segpos[j] = p;
                                    vgood[j]  = 0;
                                }
                            }
                        }

                        // add it to triangulation
                        for( PS_SIT j=0; j<2; ++j )
                        {
                            if( vgood[j] )
                            {
                                segpos[j] = holein->numberofpoints;
                                newtri2oldtri[holein->numberofpoints] = v[j];
                                holein->pointlist[holein->numberofpoints*2  ] = pos[j][0];
                                holein->pointlist[holein->numberofpoints*2+1] = pos[j][1];
                                ++holein->numberofpoints;
                            }
                        }

                        // add segment to boundary of triangulation
                        // for cases of a concave boundary
                        if( addsegment )
                        {
                            holein->segmentlist[holein->numberofsegments*2  ] = segpos[0];
                            holein->segmentlist[holein->numberofsegments*2+1] = segpos[1];
                            ++holein->numberofsegments;
                        }
                    }
                }
            }

            // enforce the segments to be the boundary
            MEMORY ps_malloc( &(holein->segmentmarkerlist), holein->numberofsegments );
            std::fill( holein->segmentmarkerlist,
                       holein->segmentmarkerlist+holein->numberofsegments, 1 );

            // triangulate with proper switches.
            // checking for numberofpoints is equivalent to checking for the
            // the case where only one triangle on the convex hull is removed
            if( holein->numberofpoints!=2 )
            {
                ps_tri_triangulate((char*)CONSTANT triswitchesholefilling,
                                   holein, holeout, (struct triangulateio*)NULL);
            }
            else
            {
                holeout->numberoftriangles = 0;
            }

            // add new triangles to temp array of triangles
            for( PS_SIT t=0; t<holeout->numberoftriangles; ++t )
            {
                // if trisadded will overflow
                if( tricapacity-1<numtriadded )
                {
                    tricapacity *= 2;
                    OPERA resize( &trisadded,
                                  numtriadded*3, tricapacity*3 );
                }

                // add triangles to temp array
                for( PS_SIT v=0; v<3; ++v )
                {
                    trisadded[numtriadded*3+v] = newtri2oldtri[holeout->trianglelist[t*3+v]];
                }
                ++numtriadded;
            }

            // check if point(s) on convex hull were removed.
            // if a point on the hull is removed, then holein->segmentlist
            // can never form a closed polygon and triangulate() will return
            // zero triangles
            if( !holeout->numberoftriangles )
            {
                // flag old edgse to be removed
                for( PS_SIT n=0; n<numremoveadded; ++n )
                {
                    for( PS_SIT e=0; e<vars_->triout->numberofedges; ++e )
                    {
                        if( vars_->triout->edgemarkerlist[e] &&
                            ( removelist[n] == vars_->triout->edgelist[e*2  ] ||
                              removelist[n] == vars_->triout->edgelist[e*2+1] ) )
                        {
                            (*edgekeep)[e] = 0;
                        }
                    }
                }

                // add new edges
                // check for overflow first
                if( edgecapacity < numedgesadded+holein->numberofsegments )
                {
                    edgecapacity = numedgesadded+holein->numberofsegments;
                    OPERA resize( &edgesadded,
                                  numedgesadded*2, edgecapacity*2 );
                }
                for( PS_SIT e=0; e<holein->numberofsegments; ++e )
                {
                    edgesadded[numedgesadded*2  ] = newtri2oldtri[holein->segmentlist[e*2  ]];
                    edgesadded[numedgesadded*2+1] = newtri2oldtri[holein->segmentlist[e*2+1]];
                    ++numedgesadded;
                }
            }

            if( data_->debug )
            {
                PS_SIT vback;
                for( PS_SIT t=0; t<holeout->numberoftriangles; ++t)
                {
                    vback = 2;
                    for(PS_SIT v=0; v<3; ++v)
                    {
                        PRINTER writeoutstream <string> (
                            "line(" +
                            OPERA tostring(holein->pointlist[holeout->trianglelist[
                                                   t*3+v    ]*2  ]) + "," +
                            OPERA tostring(holein->pointlist[holeout->trianglelist[
                                                   t*3+v    ]*2+1]) + "," +
                            OPERA tostring(holein->pointlist[holeout->trianglelist[
                                                   t*3+vback]*2  ]) + "," +
                            OPERA tostring(holein->pointlist[holeout->trianglelist[
                                                   t*3+vback]*2+1]) + ")",
                            holesfile, NULL, data_->precision, NULL );
                        vback = v;
                    }
                }
            }
            MEMORY ps_free( newtri2oldtri );
            MEMORY triangulatestructdestruct( holein  );
            MEMORY triangulatestructdestruct( holeout );
            MEMORY ps_free(                   holein  );
            MEMORY ps_free(                   holeout );
        }
    }

    tempchar = *pixelused;
    *pixelused = pixelusedorig;
    MEMORY ps_free( tempchar );

    // fix triangulation edges
    // if no edges on the convex hull were touched, forget about it
    if( !numedgesadded )
    {
        MEMORY ps_free( *edgekeep );
        *edgekeep = 0;
    }
    else
    {
        OPERA resize( &vars_->triout->edgelist,
                      vars_->triout->numberofedges * 2,
                      (vars_->triout->numberofedges + numedgesadded) * 2 );
        OPERA resize( &vars_->triout->edgemarkerlist,
                      vars_->triout->numberofedges,
                      vars_->triout->numberofedges + numedgesadded       );

        std::copy( edgesadded, edgesadded + numedgesadded * 2,
                   vars_->triout->edgelist + vars_->triout->numberofedges * 2 );
        std::copy( *edgekeep, *edgekeep + vars_->triout->numberofedges,
                   vars_->triout->edgemarkerlist                              );

        std::fill( vars_->triout->edgemarkerlist + vars_->triout->numberofedges,
                   vars_->triout->edgemarkerlist + vars_->triout->numberofedges +
                   numedgesadded, 1                                               );

        vars_->triout->numberofedges += numedgesadded;
    }
    MEMORY ps_free( edgesadded );

    // add trisadded to 'global' trianglelist
    // first find total number of triangles
    PS_SIT totaltri = numtriadded;
    for( PS_SIT t=0; t<vars_->triout->numberoftriangles; ++t )
        if( !(*triremove)[t] )
            ++totaltri;

    // now add triangles from original triangulation
    // we have to preserve the location of used trianges
    // since they are using in vars_->pointer and other vars

    // copy triangles
    PS_SIT sizetri = std::max( totaltri, vars_->triout->numberoftriangles );
    tempint = vars_->triout->trianglelist;
    MEMORY ps_malloc( &(vars_->triout->trianglelist), sizetri*3 );
    std::copy( tempint, tempint+vars_->triout->numberoftriangles*3, vars_->triout->trianglelist );
    MEMORY ps_free( tempint );

    // copy triused list
    tempchar = *triused;
    MEMORY ps_malloc( &(*triused), sizetri );
    std::copy( tempchar, tempchar+vars_->triout->numberoftriangles, *triused );
    MEMORY ps_free( tempchar );
    // mark all newly allocated triangle spots as unused
    if( totaltri>vars_->triout->numberoftriangles )
        std::fill( *triused+vars_->triout->numberoftriangles, *triused+totaltri, (char)0 );

    // copy triremove list
    tempchar = *triremove;
    MEMORY ps_malloc( &(*triremove), sizetri );
    std::copy( tempchar, tempchar+vars_->triout->numberoftriangles, *triremove );
    MEMORY ps_free( tempchar );
    // mark all newly allocated triangle spots as to remove
    if( totaltri>vars_->triout->numberoftriangles )
        std::fill( *triremove+vars_->triout->numberoftriangles, *triremove+totaltri, (char)1 );

    // now add triangles from hole-filling
    // numtriadded could be zero if an unused pixel was on the convex hull
    if( numtriadded )
    {
        PS_SIT index = 0;
        for( PS_SIT t=0; t<sizetri; ++t )
        {
            if( (*triremove)[t] )
            {
                (*triremove)[t] = 0;
                for( PS_SIT v=0; v<3; ++v )
                {
                    vars_->triout->trianglelist[t*3+v] = trisadded[index*3+v];
                }
                ++index;
                // if we're done copy all used hole triangles
                if( index==numtriadded )
                    break;
            }
        }
    }

    MEMORY ps_free( trisadded );
    vars_->triout->numberoftriangles = sizetri;

    if( data_->debug )
    {
        PRINTER closeoutstream( holesfile );
    }

    MEMORY ps_free( holesfileptr );
}
