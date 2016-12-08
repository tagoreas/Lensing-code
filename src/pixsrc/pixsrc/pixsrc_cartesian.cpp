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



#include "pixsrc_cartesian.hpp"
#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_triangulation.hpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_regularization.hpp"
#include "pixsrc_printer_templates.cpp"

/***************************************************************/
/***************************************************************/
/*****************multithreading code follows*******************/

void* parallelcreateccartesian(void *args)
{
    pixsrc_cartesian *plc = (pixsrc_cartesian*)args;
    plc->vars_->pthreadstracker[0] = 1;
    plc->createc();
    plc->vars_->pthreadstracker[0] = 2;
    return NULL;
}
void* parallelprintsourcecartesian(void *args)
{
    pixsrc_cartesian *plc = (pixsrc_cartesian*)args;
    plc->vars_->pthreadstracker[9] = 1;
    plc->printsource();
    plc->vars_->pthreadstracker[9] = 2;
    return NULL;
}
void* parallelgetmagcartesian(void *args)
{
    pixsrc_cartesian *plc = (pixsrc_cartesian*)args;
    plc->vars_->pthreadstracker[10] = 1;
    plc->getmagnification();
    plc->vars_->pthreadstracker[10] = 2;
    return NULL;
}

/*****************multithreading code ending********************/
/***************************************************************/
/***************************************************************/


pixsrc_cartesian::~pixsrc_cartesian(){}
pixsrc_cartesian::pixsrc_cartesian( PS_SIT imagenumber,
                                    inputdata *datavec, commoninputdata *cdata__, PS_SIT tracker_ )
{
    MEMORY ps_malloc( &vars_, 1 );
    vars_->tracker     =  tracker_;
    vars_->imagenumber =  imagenumber;
    vars_->pixsrc_class = (void*)this;
    data_              = &datavec[imagenumber];
    cdata_             =  cdata__;

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
        creategrid0();
        COMMON creater4r            (data_, cdata_, vars_);
        createc4c  ();

        pthread_create(&(vars_->pthreads[0] ), cdata_->attrdetached, parallelcreateccartesian    , this      );
        pthread_create(&(vars_->pthreads[7] ), cdata_->attrdetached, TRIANGULATION pstriangulate        , vars_->dav);

        COMMON findsisterimages     (data_, cdata_, vars_);
        createlo   ();
        COMMON sourcereconstructions(data_, cdata_, vars_);

        pthread_create(&(vars_->pthreads[8] ), cdata_->attrdetached, COMMON printimageplane      , vars_->dav);
        pthread_create(&(vars_->pthreads[9] ), cdata_->attrdetached, parallelprintsourcecartesian, this      );
        pthread_create(&(vars_->pthreads[10]), cdata_->attrdetached, parallelgetmagcartesian     , this      );
        pthread_create(&(vars_->pthreads[11]), cdata_->attrdetached, STATISTIC evidencecalculator   , vars_->dav);
    }

    COMMON printfinalmessages       (data_, cdata_, vars_);
    MEMORY freememory               (data_, cdata_, vars_);
}
void pixsrc_cartesian::creategrid0()
{
    if (vars_->fatalerror)
    {
        return;
    }

    if(!data_->cartdetails)
        makecg(1,-1,-1);
    else
    {
        PS_SIT srcxtemp,srcytemp;
        srcxtemp=(PS_SIT)data_->cartdetails[4];
        srcytemp=(PS_SIT)data_->cartdetails[5];

        double redd = std::max(
            (data_->invwcsinfo[0]*data_->cartdetails[2]+data_->invwcsinfo[1] *
             data_->cartdetails[3])/3600/srcxtemp ,
            (data_->invwcsinfo[2]*data_->cartdetails[2]+data_->invwcsinfo[3] *
             data_->cartdetails[3])/3600/srcytemp );

        pthread_mutex_lock( cdata_->wcsmutex );
        HEADER getimgpixcoord( data_->wcs, data_->imgy, cdata_, data_->cartdetails[0],
                               data_->cartdetails[1], &vars_->ctrx, &vars_->ctry );
        pthread_mutex_unlock( cdata_->wcsmutex );

        makecg(redd,srcxtemp,srcytemp);
    }
}
void pixsrc_cartesian::makecg( double reduction, PS_SIT srcx, PS_SIT srcy )
{
    vars_->reduction=reduction;
    if(srcx==-1)
        vars_->srcx=(PS_SIT)(floor(vars_->rangex/reduction)+1+2);
    else
        vars_->srcx=srcx;
    if(srcy==-1)
        vars_->srcy=(PS_SIT)(floor(vars_->rangey/reduction)+1+2);
    else
        vars_->srcy=srcy;

    vars_->numsrcpoints = vars_->srcx*vars_->srcy;

    MEMORY ps_free  ( vars_->c4c      );
    MEMORY ps_free  ( vars_->c4cback  );
    MEMORY ps_free  ( vars_->srcloc   );
    MEMORY ps_malloc( &(vars_->c4c    ), vars_->numsrcpoints   );
    MEMORY ps_malloc( &(vars_->c4cback), vars_->numsrcpoints   );
    MEMORY ps_malloc( &(vars_->srcloc ), vars_->numsrcpoints*2 );
    std::fill( vars_->c4c, vars_->c4c+vars_->numsrcpoints, -1 );

    //////////////
    /// Must have already set ctrX and ctrY (in image coordinates) before calling makeCG!!!
    //////////////
    double topleft[2] = {
        topleft[0] = vars_->ctrx - (vars_->srcx-1.0)/2.0*reduction,
        topleft[1] = vars_->ctry - (vars_->srcy-1.0)/2.0*reduction };

    for(PS_SIT s = 0; s < vars_->numsrcpoints; s++)
    {
        vars_->srcloc[s*2  ] = topleft[0] + reduction*floor(s/vars_->srcy);
        vars_->srcloc[s*2+1] = topleft[1] + reduction*(s%vars_->srcy);
    }
}
void pixsrc_cartesian::createc4c()
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "building grid",
                             cdata_->print2screenmutex);

    if(!data_->srcinputcircle)
    {
        double prevnewred=-1,newred=-1,prevtotalfactor=1,factor=0.1,totalfactor;
        double trivert[6];
        while(1)
        {
            for(PS_SIT r = 0; r < data_->ndp; r++)
                if(vars_->r4r[r]!=-1)
                {
                    double x = vars_->newloc[r*2  ];
                    double y = vars_->newloc[r*2+1];

                    double xt = x - vars_->srcloc[0];
                    double yt = y - vars_->srcloc[1];

                    PS_SIT surrpix[4], whichpos[4], wpsize;
                    surrpix[1] = ( (PS_SIT)floor(xt/vars_->reduction) ) *
                        vars_->srcy+(PS_SIT)floor( yt/vars_->reduction );
                    surrpix[0] = surrpix[1]+vars_->srcy;
                    surrpix[2] = surrpix[1]+1;
                    surrpix[3] = surrpix[0]+1;

                    if(OPERA equalszero(fmod(yt,vars_->reduction)))
                    {
                        wpsize = 2;
                        whichpos[0]=0;
                        whichpos[1]=1;
                    }
                    if(OPERA equalszero(fmod(xt,vars_->reduction)))
                    {
                        wpsize = 2;
                        whichpos[0]=1;
                        whichpos[1]=2;
                    }
                    else if( OPERA equalszero(
                                 fmod(xt,vars_->reduction) - vars_->reduction / 2.0) &&
                             OPERA equalszero(
                                 fmod(yt,vars_->reduction) - vars_->reduction / 2.0 ) )
                    {
                        wpsize = 4;
                        whichpos[0]=0;
                        whichpos[1]=1;
                        whichpos[2]=2;
                        whichpos[3]=3;
                    }
                    else if(GEOM isonlinesegment(
                                vars_->srcloc[surrpix[0]*2  ],
                                vars_->srcloc[surrpix[0]*2+1],
                                vars_->srcloc[surrpix[2]*2  ],
                                vars_->srcloc[surrpix[2]*2+1], x, y ) )
                    {
                        wpsize = 2;
                        whichpos[0]=0;
                        whichpos[1]=2;
                    }
                    else if( GEOM isonlinesegment( vars_->srcloc[surrpix[1]*2  ],
                                                   vars_->srcloc[surrpix[1]*2+1],
                                                   vars_->srcloc[surrpix[3]*2  ],
                                                   vars_->srcloc[surrpix[3]*2+1], x, y ) )
                    {
                        wpsize = 2;
                        whichpos[0]=1;
                        whichpos[1]=3;
                    }
                    else
                    {
                        PS_SIT pos[] = {2,3};
                        for( PS_SIT j=0; j<3; ++j )
                        {
                            trivert[j*2  ] = vars_->srcloc[surrpix[j]*2  ];
                            trivert[j*2+1] = vars_->srcloc[surrpix[j]*2+1];
                        }
                        if(GEOM isintri(trivert,x,y))
                            pos[0]=0;

                        for( PS_SIT j=0; j<3; ++j )
                        {
                            trivert[j*2  ] = vars_->srcloc[surrpix[j+1]*2  ];
                            trivert[j*2+1] = vars_->srcloc[surrpix[j+1]*2+1];
                        }
                        if(GEOM isintri(trivert,x,y))
                            pos[1]=1;

                        PS_SIT test=pos[0]+1;
                        if(test==4) test=0;
                        bool above=false;
                        if(test==pos[1]) above=true;

                        PS_SIT diff1,diff2;
                        if(above)
                        {
                            diff1=pos[0];
                            diff2=pos[1]+2;
                            if(diff2==5) diff2=1;
                            else if(diff2==4) diff2=0;
                        }
                        else
                        {
                            diff1=pos[0]+2;
                            diff2=pos[1];
                            if(diff1==5) diff1=1;
                            else if(diff1==4) diff1=0;
                        }

                        PS_SIT posf[]={pos[0],pos[0]+1,pos[0]+2};
                        if( OPERA rank( x, y,
                                        vars_->srcloc[surrpix[diff1]*2  ],
                                        vars_->srcloc[surrpix[diff1]*2+1],
                                        vars_->srcloc[surrpix[diff2]*2  ],
                                        vars_->srcloc[surrpix[diff2]*2+1] ) == 1 )
                        {
                            posf[0]=pos[1];
                            posf[1]=pos[1]+1;
                            posf[2]=pos[1]+2;
                        }
                        if(posf[1]==4) posf[1]=0;
                        if(posf[2]==4) posf[2]=0;
                        else if(posf[2]==5) posf[2]=1;

                        wpsize = 3;
                        whichpos[0]=posf[0];
                        whichpos[1]=posf[1];
                        whichpos[2]=posf[2];
                    }

                    for(PS_SIT m=0; m<wpsize; m++)
                        vars_->c4c[surrpix[whichpos[m]]]=0;
                }

            PS_SIT currnumused=0;
            for(PS_SIT m=0; m<vars_->numsrcpoints; m++)
                if(vars_->c4c[m]==0)
                    currnumused++;

            if( std::fabs(currnumused-(vars_->lonr/2.0)) <= 0.1*vars_->lonr/2.0 )
                break;

            if(currnumused-(vars_->lonr/2.0)<0)
            {
                if(prevtotalfactor>1) //prev. tried to inc. red. or dec. NSP
                    factor=factor*0.1;
                totalfactor=1-factor;
            }
            else
            {
                if(prevtotalfactor<1)
                    factor=factor*0.1;
                totalfactor=1+factor;
            }
            prevnewred=newred;
            newred=totalfactor*vars_->reduction;

            if(prevnewred==newred && newred!=-1)
                break;

            prevtotalfactor=totalfactor;
            makecg(newred,-1,-1);
        }
    }
    else
    {
        std::fill( vars_->c4c, vars_->c4c+vars_->numsrcpoints, 0 );
    }

    if(data_->debug)
    {
        PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name, 1,
                                 "srcloc.dat",vars_->srcloc, vars_->numsrcpoints, 0, data_->precision, NULL);
    }

    flagsrcpixels();

    // c4c has been flagged. 0 means used. -1 means not used. Now, we correct the 0 entries.
    vars_->lonc=0;
    PS_SIT unused=0;
    for(PS_SIT s = 0; s < vars_->numsrcpoints; s++)
        if(vars_->c4c[s]!=-1 && vars_->fallsinsrc[s])
        {
            vars_->c4c[s] = s - unused;
            vars_->c4cback[vars_->lonc++] = s;
        }
        else
        {
            vars_->c4c[s]=-1;
            unused++;
        }
}
void pixsrc_cartesian::flagsrcpixels()
{
    MEMORY ps_malloc( &(vars_->fallsinsrc), vars_->numsrcpoints );
    for(PS_SIT r=0; r<vars_->numsrcpoints; r++)
        vars_->fallsinsrc[r] =
            GEOM fallsinsrcinputcircle( vars_->srcloc[r*2], vars_->srcloc[r*2+1],
                                        data_->srcinputcircle,data_->extlengths[2] );
}
void pixsrc_cartesian::createlo()
{
    if (vars_->fatalerror)
    {
        return;
    }

    PS_SIT dim = 1;
    PS_SIT waitlist[1] = {7};
    OPERA pthreadswait(vars_,dim,waitlist);

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "creating lensing operator",
                             cdata_->print2screenmutex);

    PS_SIT initsize = vars_->lonr*3;

    MEMORY ps_malloc( &(vars_->lensingoperatornoboptr), 1 );
    vars_->lensingoperatornobo = new (vars_->lensingoperatornoboptr)
        MATRIX( cdata_, data_, vars_->lonr,vars_->lonc, initsize, 0, data_);

    // this holds the order in which column indeces should be
    // set in the lensing operator.
    // CUDA requires row-then-column ordering.
    PS_SIT setorder[4];

    for(PS_SIT r = 0; r < data_->ndp; r++)
        if(vars_->r4r[r]!=-1)
        {
            double x = vars_->newloc[r*2  ];
            double y = vars_->newloc[r*2+1];

            if(x<vars_->srcloc[0] || y<vars_->srcloc[1] ||
               x>=vars_->srcloc[(vars_->numsrcpoints-1)*2] ||
               y>=vars_->srcloc[(vars_->numsrcpoints-1)*2+1]  )
                continue;

            double xt = x - vars_->srcloc[0];
            double yt = y - vars_->srcloc[1];

            PS_SIT xpos = (PS_SIT)floor(xt/vars_->reduction);
            PS_SIT ypos = (PS_SIT)floor(yt/vars_->reduction);

            PS_SIT surrpix[4];
            surrpix[1] = xpos*vars_->srcy+ypos;
            surrpix[0] = surrpix[1]+vars_->srcy;
            surrpix[2] = surrpix[1]+1;
            surrpix[3] = surrpix[0]+1;

            double weights[4];
            double temp1[8];
            PS_SIT temp2[4];
            for(PS_SIT v=0; v<4; v++)
            {
                temp2[v] = vars_->c4c[surrpix[v]];
                temp1[v*2  ] = vars_->srcloc[surrpix[v]*2  ];
                temp1[v*2+1] = vars_->srcloc[surrpix[v]*2+1];
            }

            OPERA linearinterpolatorcartesian( x, y, vars_->reduction, temp1, temp2, weights );

            // recording column positions in setorder
            PS_SIT setpos = 0;
            for( PS_SIT g=0; g<4; ++g )
            {
                if(!OPERA equalszero(weights[g]))
                    setorder[setpos++] = vars_->c4c[surrpix[g]];
            }

            // sorting column indices
            std::sort( setorder, setorder + setpos );

            // setting lensing operator
            for( PS_SIT sp=0; sp<setpos; ++sp )
            {
                for( PS_SIT g=0; g<4; ++g )
                    if( setorder[sp] ==  vars_->c4c[surrpix[g]] )
                    {
                        vars_->lensingoperatornobo->set(
                            vars_->r4r[r], vars_->c4c[surrpix[g]], weights[g] );
                        break;
                    }
            }
        }

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "done creating lensing operator",
                             cdata_->print2screenmutex);

    COMMON blurit(data_, cdata_, vars_);

}
void pixsrc_cartesian::createc()
{
    if (vars_->fatalerror)
    {
        return;
    }

    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : 0;

    MEMORY ps_malloc( &(vars_->c1ptr), 1 );
    vars_->c1 = new (vars_->c1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 1, data_);

    if(data_->regorder==2)
        REGULAR getdercartesian2( data_, cdata_, vars_ );
    else if(data_->regorder==1)
        REGULAR getdercartesian1( data_, cdata_, vars_ );
    else if(data_->regorder==0)
        REGULAR getdercartesian0( data_, cdata_, vars_ );

    pthread_create(&(vars_->pthreads[1]), cdata_->attrdetached,
                   STATISTIC computedetc, vars_->dav            );

    vars_->h1->mult( vars_->h1, vars_->c1,
                     1, 0, cdata_->numthreads );

    if(data_->debug)
    {
        PRINTER print( vars_->tracker, cdata_->basename, data_->name, 1,
                       "h1.MATRIX", vars_->lonc, vars_->lonc, vars_->h1, 0 );
        PRINTER print( vars_->tracker, cdata_->basename, data_->name, 1,
                       "c1.MATRIX", vars_->lonc, vars_->lonc, vars_->c1, 0 );
    }
}
void pixsrc_cartesian::getmagnification()
{
    if (vars_->fatalerror)
    {
        return;
    }

    if(data_->magparams || data_->penaltyquery[1])
    {
        double srcflux = 0, imgflux = 0;

        for(PS_SIT s=0; s<vars_->lonc; s++)
            srcflux += vars_->mps->get(s);
        srcflux *= vars_->reduction*vars_->reduction;

        for(PS_SIT s=0; s<vars_->lonr; s++)
            imgflux += vars_->data->get(s);

        vars_->magger = imgflux/srcflux;
        if(data_->verbose!=1)
            PRINTER print2screen(data_->print2screenname, "magnification #" +
                                 OPERA tostring(vars_->tracker) + " = " +
                                 OPERA tostring(vars_->magger),
                                 cdata_->print2screenmutex                     );

        double writeout[2];
        writeout[0] = vars_->lambda1;
        writeout[1] = vars_->magger;
        PRINTER writeoutstream <double> (writeout,2,data_->mags->stream,data_->mags->lock, data_->precision, NULL);

        COMMON computemagpenalty(data_, cdata_, vars_);
    }
}
void pixsrc_cartesian::printsource()
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->printvec==1 || data_->printvec==3 )
    {
        double **mps;
        MEMORY ps_malloc( &mps, vars_->lonc, 3 );
        for(PS_SIT x=0; x<vars_->lonc; x++)
        {
            mps[x][0] = vars_->srcloc[vars_->c4cback[x]*2  ];
            mps[x][1] = vars_->srcloc[vars_->c4cback[x]*2+1];
            mps[x][2] = vars_->mps->get(x);
        }

        PRINTER print <double> ( vars_->tracker, cdata_->basename, data_->name,
                                 1, "mps.VECTOR", mps, vars_->lonc, 3, 0, data_->precision, 0);
        MEMORY ps_free( mps, vars_->lonc );
    }

    if( data_->printvec==2 || data_->printvec==3 )
    {
        double *im;
        MEMORY ps_malloc( &(im), vars_->srcx*vars_->srcy );
        PS_SIT index = 0;
        for(PS_SIT y=vars_->srcy-1; y>=0; y--)
            for(PS_SIT x=0; x<vars_->srcx; x++)
            {
                PS_SIT s=x*vars_->srcy+y;
                if(vars_->c4c[s]!=-1)
                    im[index++] = (vars_->mps->get(vars_->c4c[s]));
                else
                    im[index++] = 0;
            }

        if(data_->cartdetails)
            HEADER setsrcwcs(data_, vars_, (vars_->srcx+1)/2.0-1,(vars_->srcy+1)/2.0-1,
                             data_->cartdetails[0],data_->cartdetails[1]);
        else
        {
            double one,two;
            pthread_mutex_lock( cdata_->wcsmutex );
            HEADER getimgwcscoord( data_->wcs, data_->imgy,
                                   vars_->srcloc[(((vars_->srcx-1)/2)*vars_->srcy)*2] +
                                   vars_->reduction/2*((vars_->srcx+1)%2),
                                   vars_->srcloc[((vars_->srcy-1)/2)*2+1] +
                                   vars_->reduction/2*((vars_->srcx+1)%2), &one, &two );
            pthread_mutex_unlock( cdata_->wcsmutex );
            HEADER setsrcwcs(data_, vars_, (vars_->srcx+1)/2.0-1,(vars_->srcy+1)/2.0-1,one,two);
        }

        PRINTER printfitssrcplane( vars_->tracker, cdata_->basename, data_->name,
                                   data_->print2screenname, 1,"mps.fits",vars_->srcx,
                                   vars_->srcy,im,vars_->wcsinfo, cdata_->print2screenmutex, 1, 0 );

        MEMORY ps_free( im );
    }
}
