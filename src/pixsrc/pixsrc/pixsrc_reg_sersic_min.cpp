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



#include "pixsrc_regularization.hpp"
#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer.hpp"

#include <cmath>
#include <algorithm>

void sersicminscan1( inputdata*, commoninputdata*, lensvar*  );
void sersicminscan2( inputdata*, commoninputdata*, lensvar*,
                     PS_SIT, PS_SIT*                               );
void tryallindex   ( inputdata*, commoninputdata*, lensvar*,
                     PS_SIT*, double, double, PS_SIT               );
void trysomeindex  ( inputdata*, commoninputdata*, lensvar*,
                     double*, double*, PS_SIT*, PS_SIT             );

void pixsrc_regularization::getsersicregmin(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );

    sersicminscan1( data_, cdata_, vars_ );

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
}

void sersicminscan1(inputdata *data_, commoninputdata *cdata_, lensvar *vars_ )
{
    double xx, xxx, yy, yyy;
    bool p[4];
    PS_SIT closest[4];
    double closestR2[4];

    for(PS_SIT cc = 0; cc < vars_->lonc; cc++)
    {
        xx = vars_->triout->pointlist[cc*2  ];
        yy = vars_->triout->pointlist[cc*2+1];

        std::fill( p, p + 4, 0 );
        std::fill( closest, closest + 4, -1 );

        OPERA assign_p_infinity( &closestR2[0] );
        std::fill( closestR2 + 1, closestR2 + 4, closestR2[0] );

        for(PS_SIT ccc = 0; ccc < vars_->lonc; ccc++)
            if(cc!=ccc)
            {
                xxx = vars_->triout->pointlist[ccc*2  ] - xx;
                yyy = vars_->triout->pointlist[ccc*2+1] - yy;

                if(OPERA equalszero(xxx))
                    xxx=0;
                if(OPERA equalszero(yyy))
                    yyy=0;

                double radius2 = xxx*xxx + yyy*yyy;

                if(radius2 < closestR2[0] && xxx>0 && yyy<=0)
                {
                    closestR2[0] = radius2;
                    closest[0] = ccc;
                    p[0] = 1;
                }
                if(radius2 < closestR2[1] && xxx<=0 && yyy<0)
                {
                    closestR2[1] = radius2;
                    closest[1] = ccc;
                    p[1] = 1;
                }
                if(radius2 < closestR2[2] && xxx<0 && yyy>=0)
                {
                    closestR2[2] = radius2;
                    closest[2] = ccc;
                    p[2] = 1;
                }
                if(radius2 < closestR2[3] && xxx>=0 && yyy>0)
                {
                    closestR2[3] = radius2;
                    closest[3] = ccc;
                    p[3] = 1;
                }
            }

        if(!(p[0] && p[1] && p[2] && p[3]))
        {
            sersicminscan2( data_, cdata_, vars_, cc, closest );
        }
        else
            tryallindex( data_, cdata_, vars_, closest, xx, yy, cc );
    }
}
void sersicminscan2(inputdata *data_, commoninputdata *cdata_, lensvar *vars_, PS_SIT cc, PS_SIT *closest )
{
    bool p[4];
    std::fill(p,p+4,0);
    PS_SIT total=4;
    for(PS_SIT x=0; x<4; x++)
        if(closest[x]==-1)
            total--;
        else p[x]=1;
    double xxx[4];
    double yyy[4];
    double xx = vars_->triout->pointlist[cc*2  ];
    double yy = vars_->triout->pointlist[cc*2+1];

    switch(total)
    {
    case 0: PRINTER printerror(data_->print2screenname,
                               "error in regularization. isolated pixel.",
                               cdata_->print2screenmutex                   );
        break;
    case 1: for(PS_SIT x=0; x<4; x++)
            if(p[x])
            {
                xxx[x] = vars_->triout->pointlist[closest[x]*2  ] - xx;
                yyy[x] = vars_->triout->pointlist[closest[x]*2+1] - yy;

                if(OPERA equalszero(xxx[x]))
                    xxx[x]=0;
                if(OPERA equalszero(yyy[x]))
                    yyy[x]=0;

                if(xxx[x]==0 || yyy[x]==0)
                {
                    double pos;
                    if(xxx[x]==0)
                        pos=std::fabs(yyy[x]);
                    else
                        pos=std::fabs(xxx[x]);

                    xxx[0] = pos;
                    yyy[0] = 0;
                    xxx[1] = 0;
                    yyy[1] = -pos;
                    xxx[2] = -pos;
                    yyy[2] = 0;
                    xxx[3] = 0;
                    yyy[3] = pos;
                }
                else
                {
                    double posX = xxx[x];
                    double posY = yyy[x];
                    if(posX==0) posX=1;
                    if(posY==0) posY=1;

                    switch(x)
                    {
                    case 0: xxx[1]=-posX; yyy[1]= posY;
                        xxx[2]=-posX; yyy[2]=-posY;
                        xxx[3]= posX; yyy[3]=-posY;
                        break;
                    case 1: xxx[0]=-posX; yyy[0]= posY;
                        xxx[2]= posX; yyy[2]=-posY;
                        xxx[3]=-posX; yyy[3]=-posY;
                        break;
                    case 2: xxx[1]= posX; yyy[1]=-posY;
                        xxx[0]=-posX; yyy[0]=-posY;
                        xxx[3]=-posX; yyy[3]= posY;
                        break;
                    case 3: xxx[1]=-posX; yyy[1]=-posY;
                        xxx[2]=-posX; yyy[2]= posY;
                        xxx[0]= posX; yyy[0]=-posY;
                        break;
                    }
                }
            }
        break;
    case 2: for(PS_SIT x=0; x<4; x++)
            if(p[x])
                for(PS_SIT y=0; y<4; y++)
                    if(x!=y && p[y])
                    {
                        xxx[x] = vars_->triout->pointlist[closest[x]*2  ] - xx;
                        yyy[x] = vars_->triout->pointlist[closest[x]*2+1] - yy;
                        xxx[y] = vars_->triout->pointlist[closest[y]*2  ] - xx;
                        yyy[y] = vars_->triout->pointlist[closest[y]*2+1] - yy;

                        if(OPERA equalszero(xxx[x]))
                            xxx[x]=0;
                        if(OPERA equalszero(yyy[x]))
                            yyy[x]=0;
                        if(OPERA equalszero(xxx[y]))
                            xxx[y]=0;
                        if(OPERA equalszero(yyy[y]))
                            yyy[y]=0;

                        /*
                          double avg[2] = {0,0};
                          avg[0]=(fabs(xxx[x])+fabs(xxx[y]))/2.0;
                          avg[1]=(fabs(yyy[x])+fabs(yyy[y]))/2.0;
                          for(PS_SIT b=0; b<2; b++)
                          if(avg[b]==0)
                          avg[b]=1;


                          for(PS_SIT z=0; z<4; z++)
                          if(z!=x && z!=y)
                          {
                          xxx[z]=avg[0];
                          yyy[z]=avg[1];
                          }
                          if(xxx[0]<0)xxx[0]*=-1;
                          if(yyy[0]>0)yyy[0]*=-1;
                          if(xxx[1]>0)xxx[1]*=-1;
                          if(yyy[1]>0)yyy[1]*=-1;
                          if(xxx[2]>0)xxx[2]*=-1;
                          if(yyy[2]<0)yyy[2]*=-1;
                          if(xxx[3]<0)xxx[3]*=-1;
                          if(yyy[3]<0)yyy[3]*=-1;
                        */

                        PS_SIT ind1 = x+2;
                        PS_SIT ind2 = y+2;
                        if(ind1>3)
                            ind1-=4;
                        if(ind2>3)
                            ind2-=4;
                        xxx[ind1] = -xxx[x];
                        yyy[ind1] = -yyy[x];
                        xxx[ind2] = -xxx[y];
                        yyy[ind2] = -yyy[y];
                    }
        break;
    case 3: for(PS_SIT x=0; x<4; x++)
            if(!p[x])
            {
                double avg[2] = {0,0};

                for(PS_SIT y=0; y<4; y++)
                    if(p[y])
                    {
                        xxx[y] = vars_->triout->pointlist[closest[y]*2  ] - xx;
                        yyy[y] = vars_->triout->pointlist[closest[y]*2+1] - yy;

                        if(OPERA equalszero(xxx[y]))
                            xxx[y]=0;
                        if(OPERA equalszero(yyy[y]))
                            yyy[y]=0;

                        avg[0]+=fabs(xxx[y]); avg[1]+=fabs(yyy[y]);

                    }

                /*
                  xxx[x]=avg[0]/3.0;
                  yyy[x]=avg[1]/3.0;

                  if(xxx[x]==0)xxx[x]=1;
                  if(yyy[x]==0)yyy[x]=1;
                */

                /*
                  switch(x)
                  {
                  case 0: yyy[x]*=-1;
                  break;
                  case 1: xxx[x]*=-1; yyy[x]*=-1;
                  break;
                  case 2: xxx[x]*=-1;
                  break;
                  case 3: break;
                  }
                */

                PS_SIT ind=x-2;
                if(ind<0)
                    ind+=4;
                xxx[x]=-xxx[ind];
                yyy[x]=-yyy[ind];
            }
        break;
    case 4: tryallindex( data_, cdata_, vars_, closest, xx, yy, cc );
        break;
    }
    if(total!=4)
        trysomeindex( data_, cdata_, vars_, xxx, yyy, closest, cc );
}
void tryallindex(inputdata *data_, commoninputdata *cdata_, lensvar *vars_, PS_SIT *closest, double xx, double yy, PS_SIT cc )
{

    // sort column indices
    PS_SIT order[5];
    order[0] = cc;
    std::copy( closest, closest + 4, &order[1] );
    std::sort( order, order + 5 );

    // set in order
    for( PS_SIT f=0; f<5; ++f )
    {
        if( order[f] == cc )
            vars_->h1->set(cc,cc        , (-4.0 / vars_->mps->get(cc        )) );

        else if( order[f] == closest[0] )
            vars_->h1->set(cc,closest[0], ( 1.0 / vars_->mps->get(closest[0])) );

        else if( order[f] == closest[1] )
            vars_->h1->set(cc,closest[1], ( 1.0 / vars_->mps->get(closest[1])) );

        else if( order[f] == closest[2] )
            vars_->h1->set(cc,closest[2], ( 1.0 / vars_->mps->get(closest[2])) );

        else
            vars_->h1->set(cc,closest[3], ( 1.0 / vars_->mps->get(closest[3])) );
    }

    /*
    // original
    vars_->h1->set(cc,cc        , -4.0/vars_->mps->get(cc        ) );
    vars_->h1->set(cc,closest[0],  1.0/vars_->mps->get(closest[0]) );
    vars_->h1->set(cc,closest[1],  1.0/vars_->mps->get(closest[1]) );
    vars_->h1->set(cc,closest[2],  1.0/vars_->mps->get(closest[2]) );
    vars_->h1->set(cc,closest[3],  1.0/vars_->mps->get(closest[3]) );
    */
    /*
      vars_->h1->set( cc,cc        , -4.0                                            );
      vars_->h1->set( cc,closest[0], vars_->mps->get(cc)/vars_->mps->get(closest[0]) );
      vars_->h1->set( cc,closest[1], vars_->mps->get(cc)/vars_->mps->get(closest[1]) );
      vars_->h1->set( cc,closest[2], vars_->mps->get(cc)/vars_->mps->get(closest[2]) );
      vars_->h1->set( cc,closest[3], vars_->mps->get(cc)/vars_->mps->get(closest[3]) );
    */
}
void trysomeindex(inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *xxx, double *yyy, PS_SIT *ccc, PS_SIT cc )
{
    double total = 0;
    for( PS_SIT x=0; x<4; ++x )
        if( ccc[0]!=-1 )
            ++total;

    // sort column indices
    PS_SIT order[5];
    order[0] = cc;
    std::copy( ccc, ccc + 4, &order[1] );
    std::sort( order, order + 5 );

    // set in order
    for( PS_SIT f=0; f<5; ++f )
    {
        if( order[f] == -1 )
            continue;

        if( order[f] == cc )
            vars_->h1->set(cc, cc   , (-total/ vars_->mps->get(cc    ) ));

        else if( order[f] == ccc[0] )
            vars_->h1->set(cc,ccc[0], ( 1.0  / vars_->mps->get(ccc[0]) ));

        else if( order[f] == ccc[1] )
            vars_->h1->set(cc,ccc[1], ( 1.0  / vars_->mps->get(ccc[1]) ));

        else if( order[f] == ccc[2] )
            vars_->h1->set(cc,ccc[2], ( 1.0  / vars_->mps->get(ccc[2]) ));

        else
            vars_->h1->set(cc,ccc[3], ( 1.0  / vars_->mps->get(ccc[3]) ));
    }

    /*
    // original
    vars_->h1->set(cc, cc   , -total/vars_->mps->get(cc) );
    if(ccc[0]!=-1) vars_->h1->set(cc,ccc[0],  1.0/vars_->mps->get(ccc[0]) );
    if(ccc[1]!=-1) vars_->h1->set(cc,ccc[1],  1.0/vars_->mps->get(ccc[1]) );
    if(ccc[2]!=-1) vars_->h1->set(cc,ccc[2],  1.0/vars_->mps->get(ccc[2]) );
    if(ccc[3]!=-1) vars_->h1->set(cc,ccc[3],  1.0/vars_->mps->get(ccc[3]) );
    */
    /*
      vars_->h1->set(cc, cc   , -total );
      if(ccc[0]!=-1)
      vars_->h1->set( cc, ccc[0], vars_->mps->get(cc)/vars_->mps->get(ccc[0]) );
      if(ccc[1]!=-1)
      vars_->h1->set( cc, ccc[1], vars_->mps->get(cc)/vars_->mps->get(ccc[1]) );
      if(ccc[2]!=-1)
      vars_->h1->set( cc, ccc[2], vars_->mps->get(cc)/vars_->mps->get(ccc[2]) );
      if(ccc[3]!=-1)
      vars_->h1->set( cc, ccc[3], vars_->mps->get(cc)/vars_->mps->get(ccc[3]) );
    */
}
