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
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer.hpp"
#include <cmath>
#include <algorithm>

void getCurvature ( inputdata*, commoninputdata*, lensvar*  );
void getStragglers( inputdata*, commoninputdata*, lensvar*,
                    PS_SIT, PS_SIT*                               );
void trywithindex ( inputdata*, commoninputdata*, lensvar*,
                    PS_SIT*, double, double, PS_SIT               );
void trywithcoord ( inputdata*, commoninputdata*, lensvar*,
                    double*, double*, PS_SIT*, PS_SIT             );

void pixsrc_regularization::getder12(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if( data_->regorder == 2 )
    {
        PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

        MEMORY ps_malloc( &(vars_->h1ptr), 1 );
        vars_->h1 = new (vars_->h1ptr)
            MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );
    }
    else
    {
        MEMORY ps_malloc( &(vars_->h1xptr), 1 );
        MEMORY ps_malloc( &(vars_->h1yptr), 1 );

        PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

        vars_->h1x = new (vars_->h1xptr)
            MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );
        vars_->h1y = new (vars_->h1yptr)
            MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );
    }

    getCurvature( data_, cdata_, vars_ );

    if( data_->regorder == 1 )
    {
        MATRIX *dummy1, *dummy2;
        MEMORY ps_malloc( &dummy1, 1 );
        MEMORY ps_malloc( &dummy2, 1 );

        PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : 0;

        MATRIX *h1x2 = new (dummy1)
            MATRIX( cdata_, data_, vars_->lonc,vars_->lonc, initsize, 1, data_ );
        MATRIX *h1y2 = new (dummy2)
            MATRIX( cdata_, data_, vars_->lonc,vars_->lonc, initsize, 1, data_ );
        vars_->h1x->mult( vars_->h1x,h1x2,1,0,cdata_->numthreads );
        vars_->h1y->mult( vars_->h1y,h1y2,1,0,cdata_->numthreads );

        h1x2->plus( h1y2,vars_->c1,0,0,1,1,cdata_->numthreads );

        h1x2->~MATRIX();
        h1y2->~MATRIX();

        MEMORY ps_free( dummy1 );
        MEMORY ps_free( dummy2 );
    }
    else
    {
        vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
    }
}

void getCurvature(inputdata *data_, commoninputdata *cdata_, lensvar *vars_ )
{
    double xx, yy, xxx, yyy;
    bool p[4], pp[4];
    PS_SIT closest[4];
    double closestR2[4];

    for(PS_SIT cc = 0; cc < vars_->lonc; cc++)
    {
        xx = vars_->triout->pointlist[cc*2+0];
        yy = vars_->triout->pointlist[cc*2+1];

        std::fill(p,p+4,0);
        std::fill(pp,pp+4,0);
        std::fill(closest,closest+4,-1);
        OPERA assign_p_infinity( &closestR2[0] );
        std::fill( closestR2 + 1, closestR2 + 4, closestR2[0] );

        // search triangle list for surrounding points
        vector <PS_SIT> surrpoints;
        PS_SIT foundtri;
        for (PS_SIT ss=0; ss<vars_->triout->numberoftriangles*3; ++ss)
        {
            // find triangles that contain cc and add other 2 points to list
            if (vars_->triout->trianglelist[ss]==cc)
            {
                foundtri = ss/3;
                for (PS_SIT v=0; v<3; ++v)
                    if (vars_->triout->trianglelist[foundtri*3+v] != cc)
                        surrpoints.push_back (vars_->triout->trianglelist[foundtri*3+v]);
            }
        }
        // remove duplicate points
        std::sort (surrpoints.begin(), surrpoints.end());
        surrpoints.erase (std::unique (surrpoints.begin(), surrpoints.end()), surrpoints.end());
        // find nearest points of all these (in each quadrant)
        for (PS_unsignedSIT cccc = 0; cccc < surrpoints.size(); cccc++)
        {
            PS_SIT ccc = surrpoints[cccc];
            if (cc!=ccc)
            {
                xxx = vars_->triout->pointlist[ccc*2+0] - xx;
                yyy = vars_->triout->pointlist[ccc*2+1] - yy;

                if(OPERA equalszero(xxx))
                    xxx=0;
                if(OPERA equalszero(yyy))
                    yyy=0;

                double radius2 = xxx*xxx + yyy*yyy;

                if(!pp[0] && radius2 < closestR2[0] && xxx>0 && yyy<=0)
                {
                    closestR2[0] = radius2;
                    closest[0] = ccc;
                    p[0] = 1;
                }
                if(!pp[1] && radius2 < closestR2[1] && xxx<=0 && yyy<0)
                {
                    closestR2[1] = radius2;
                    closest[1] = ccc;
                    p[1] = 1;
                }
                if(!pp[2] && radius2 < closestR2[2] && xxx<0 && yyy>=0)
                {
                    closestR2[2] = radius2;
                    closest[2] = ccc;
                    p[2] = 1;
                }
                if(!pp[3] && radius2 < closestR2[3] && xxx>=0 && yyy>0)
                {
                    closestR2[3] = radius2;
                    closest[3] = ccc;
                    p[3] = 1;
                }
            }
        }
        std::copy (p,p+4,pp);


        // if didn't get all quadrants, then search for points that surround the surrounding points
        // search triangle list for surrounding points
        if (!(p[0] && p[1] && p[2] && p[3]))
        {
            bool addthistri;
            PS_SIT numsurrpoints = surrpoints.size();
            for (PS_SIT ss=0; ss<vars_->triout->numberoftriangles*3; ++ss)
            {
                // find triangles that contain surrounding points and add other 2 points to list
                addthistri = 0;
                for (PS_SIT s3=0; s3<numsurrpoints; ++s3)
                    if (vars_->triout->trianglelist[ss]==surrpoints[s3])
                    {
                        addthistri = 1;
                        break;
                    }
                if (addthistri)
                {
                    foundtri = ss/3;
                    for (PS_SIT v=0; v<3; ++v)
                        if (vars_->triout->trianglelist[foundtri*3+v] != cc)
                            surrpoints.push_back (vars_->triout->trianglelist[foundtri*3+v]);
                }
            }
            // remove duplicate points
            std::sort (surrpoints.begin(), surrpoints.end());
            surrpoints.erase (std::unique (surrpoints.begin(), surrpoints.end()), surrpoints.end());
            // find nearest points of all these (in each quadrant)
            for (PS_unsignedSIT cccc = 0; cccc < surrpoints.size(); cccc++)
            {
                PS_SIT ccc = surrpoints[cccc];
                if (cc!=ccc)
                {
                    xxx = vars_->triout->pointlist[ccc*2+0] - xx;
                    yyy = vars_->triout->pointlist[ccc*2+1] - yy;

                    if(OPERA equalszero(xxx))
                        xxx=0;
                    if(OPERA equalszero(yyy))
                        yyy=0;

                    double radius2 = xxx*xxx + yyy*yyy;

                    if(!pp[0] && radius2 < closestR2[0] && xxx>0 && yyy<=0)
                    {
                        closestR2[0] = radius2;
                        closest[0] = ccc;
                        p[0] = 1;
                    }
                    if(!pp[1] && radius2 < closestR2[1] && xxx<=0 && yyy<0)
                    {
                        closestR2[1] = radius2;
                        closest[1] = ccc;
                        p[1] = 1;
                    }
                    if(!pp[2] && radius2 < closestR2[2] && xxx<0 && yyy>=0)
                    {
                        closestR2[2] = radius2;
                        closest[2] = ccc;
                        p[2] = 1;
                    }
                    if(!pp[3] && radius2 < closestR2[3] && xxx>=0 && yyy>0)
                    {
                        closestR2[3] = radius2;
                        closest[3] = ccc;
                        p[3] = 1;
                    }
                }
            }
            std::copy (p,p+4,pp);
        }

        // if a nearest point in a quadrant was still not found, do a brute force search
        if (!(p[0] && p[1] && p[2] && p[3]))
        {
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

                    if(!pp[0] && radius2 < closestR2[0] && xxx>0 && yyy<=0)
                    {
                        closestR2[0] = radius2;
                        closest[0] = ccc;
                        p[0] = 1;
                    }
                    if(!pp[1] && radius2 < closestR2[1] && xxx<=0 && yyy<0)
                    {
                        closestR2[1] = radius2;
                        closest[1] = ccc;
                        p[1] = 1;
                    }
                    if(!pp[2] && radius2 < closestR2[2] && xxx<0 && yyy>=0)
                    {
                        closestR2[2] = radius2;
                        closest[2] = ccc;
                        p[2] = 1;
                    }
                    if(!pp[3] && radius2 < closestR2[3] && xxx>=0 && yyy>0)
                    {
                        closestR2[3] = radius2;
                        closest[3] = ccc;
                        p[3] = 1;
                    }
                }
        }

        if( !( p[0] && p[1] && p[2] && p[3] ) )
        {
            getStragglers( data_, cdata_, vars_, cc, closest );
        }
        else
            trywithindex( data_, cdata_, vars_, closest, xx, yy, cc );
    }
}
void getStragglers(inputdata *data_, commoninputdata *cdata_, lensvar *vars_, PS_SIT cc, PS_SIT *closest )
{

    bool p[4];
    std::fill(p,p+4,false);
    PS_SIT total=4;
    for(PS_SIT x=0; x<4; x++)
        if(closest[x]==-1)
            total--;
        else p[x]=true;
    double xxx[4];
    double yyy[4];
    double xx = vars_->triout->pointlist[cc*2  ];
    double yy = vars_->triout->pointlist[cc*2+1];

    switch( total )
    {
    case 0: PRINTER printerror( data_->print2screenname,
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
    case 4: trywithindex( data_, cdata_, vars_, closest, xx, yy, cc );
        break;
    }
    if(total!=4)
        trywithcoord( data_, cdata_, vars_, xxx, yyy, closest, cc );
}
void trywithindex(inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                  PS_SIT *closest, double xx, double yy, PS_SIT cc                 )
{
    double xxx[4];
    double yyy[4];
    for(PS_SIT n=0; n<4; n++)
    {
        xxx[n] = vars_->triout->pointlist[closest[n]*2  ] - xx;
        yyy[n] = vars_->triout->pointlist[closest[n]*2+1] - yy;

        if(OPERA equalszero(xxx[n]))
            xxx[n]=0;
        if(OPERA equalszero(yyy[n]))
            yyy[n]=0;
    }

    double dEC=0,d14=0,dE4=0,dE1=0,dAC=0,d12=0,dA1=0,dA2=0,
        dBC=0,d23=0,dB2=0,dB3=0,dDC=0,d34=0,dD3=0,dD4=0;
    dEC = fabs(OPERA xintercept(xxx[0],yyy[0],xxx[3],yyy[3]));
    d14 = OPERA distance(xxx[0],yyy[0],xxx[3],yyy[3]);
    dE4 = OPERA distance(xxx[3],yyy[3],dEC,0.0);
    dE1 = OPERA distance(xxx[0],yyy[0],dEC,0.0);
    dAC = fabs(OPERA intercept(xxx[0],yyy[0],xxx[1],yyy[1]));
    d12 = OPERA distance(xxx[0],yyy[0],xxx[1],yyy[1]);
    dA1 = OPERA distance(xxx[0],yyy[0],0.0,-dAC);
    dA2 = OPERA distance(xxx[1],yyy[1],0.0,-dAC);
    dBC = fabs(OPERA xintercept(xxx[1],yyy[1],xxx[2],yyy[2]));
    d23 = OPERA distance(xxx[1],yyy[1],xxx[2],yyy[2]);
    dB2 = OPERA distance(xxx[1],yyy[1],-dBC,0.0);
    dB3 = OPERA distance(xxx[2],yyy[2],-dBC,0.0);
    dDC = fabs(OPERA intercept(xxx[2],yyy[2],xxx[3],yyy[3]));
    d34 = OPERA distance(xxx[2],yyy[2],xxx[3],yyy[3]);
    dD3 = OPERA distance(xxx[2],yyy[2],0.0,dDC);
    dD4 = OPERA distance(xxx[3],yyy[3],0.0,dDC);

    /*
      double zer = -(1/dAC + 1/dBC + 1/dDC + 1/dEC);
      double one = (dE4/d14/dEC + dA2/d12/dAC);
      double two = (dA1/d12/dAC + dB3/d23/dBC);
      double thr = (dB2/d23/dBC + dD4/d34/dDC);
      double fou = (dE1/d14/dEC + dD3/d34/dDC);
    */

    double zer ,one ,two ,thr ,fou ;
    double zerx,onex,twox,thrx,foux;
    double zery,oney,twoy,thry,fouy;

    if( data_->regorder == 2 )
    {
        zer = -(1/dAC + 1/dBC + 1/dDC + 1/dEC);
        one =  (dE4/d14/dEC + dA2/d12/dAC);
        two =  (dA1/d12/dAC + dB3/d23/dBC);
        thr =  (dB2/d23/dBC + dD4/d34/dDC);
        fou =  (dE1/d14/dEC + dD3/d34/dDC);

        if( OPERA is_finite( zer+one+two+thr+fou ) )
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
                    vars_->h1->set( cc, cc,         zer );

                else if( order[f] == closest[0] )
                    vars_->h1->set( cc, closest[0], one );

                else if( order[f] == closest[1] )
                    vars_->h1->set( cc, closest[1], two );

                else if( order[f] == closest[2] )
                    vars_->h1->set( cc, closest[2], thr );

                else
                    vars_->h1->set( cc, closest[3], fou );
            }
        }
    }
    else
    {
        zerx =  1/dEC-1/dBC;
        onex = -dE4/d14/dEC;
        twox =  dB3/d23/dBC;
        thrx =  dB2/d23/dBC;
        foux = -dE1/d14/dEC;

        zery =  1/dAC-1/dDC;
        oney = -dA2/d12/dAC;
        twoy = -dA1/d12/dAC;
        thry =  dD4/d34/dDC;
        fouy =  dD3/d34/dDC;

        if( OPERA is_finite( zerx+onex+twox+thrx+foux+zery+oney+twoy+thry+fouy ) )
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
                {
                    vars_->h1x->set( cc, cc,         zerx );
                    vars_->h1y->set( cc, cc,         zery );
                }
                else if( order[f] == closest[0] )
                {
                    vars_->h1x->set( cc, closest[0], onex );
                    vars_->h1y->set( cc, closest[0], oney );
                }
                else if( order[f] == closest[1] )
                {
                    vars_->h1x->set( cc, closest[1], twox );
                    vars_->h1y->set( cc, closest[1], twoy );
                }
                else if( order[f] == closest[2] )
                {
                    vars_->h1x->set( cc, closest[2], thrx );
                    vars_->h1y->set( cc, closest[2], thry );
                }
                else
                {
                    vars_->h1x->set( cc, closest[3], foux );
                    vars_->h1y->set( cc, closest[3], fouy );
                }
            }
        }
    }
}
void trywithcoord(inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *xxx, double *yyy, PS_SIT *ccc, PS_SIT cc )
{
    double dEC=0,d14=0,dE4=0,dE1=0,dAC=0,d12=0,dA1=0,dA2=0,
        dBC=0,d23=0,dB2=0,dB3=0,dDC=0,d34=0,dD3=0,dD4=0;
    dEC = fabs(OPERA xintercept(xxx[0],yyy[0],xxx[3],yyy[3]));
    d14 = OPERA distance(xxx[0],yyy[0],xxx[3],yyy[3]);
    dE4 = OPERA distance(xxx[3],yyy[3],dEC,0.0);
    dE1 = OPERA distance(xxx[0],yyy[0],dEC,0.0);
    dAC = fabs(OPERA intercept(xxx[0],yyy[0],xxx[1],yyy[1]));
    d12 = OPERA distance(xxx[0],yyy[0],xxx[1],yyy[1]);
    dA1 = OPERA distance(xxx[0],yyy[0],0.0,-dAC);
    dA2 = OPERA distance(xxx[1],yyy[1],0.0,-dAC);
    dBC = fabs(OPERA xintercept(xxx[1],yyy[1],xxx[2],yyy[2]));
    d23 = OPERA distance(xxx[1],yyy[1],xxx[2],yyy[2]);
    dB2 = OPERA distance(xxx[1],yyy[1],-dBC,0.0);
    dB3 = OPERA distance(xxx[2],yyy[2],-dBC,0.0);
    dDC = fabs(OPERA intercept(xxx[2],yyy[2],xxx[3],yyy[3]));
    d34 = OPERA distance(xxx[2],yyy[2],xxx[3],yyy[3]);
    dD3 = OPERA distance(xxx[2],yyy[2],0.0,dDC);
    dD4 = OPERA distance(xxx[3],yyy[3],0.0,dDC);

    /*
      double zer = -(1/dAC + 1/dBC + 1/dDC + 1/dEC);
      double one = (dE4/d14/dEC + dA2/d12/dAC);
      double two = (dA1/d12/dAC + dB3/d23/dBC);
      double thr = (dB2/d23/dBC + dD4/d34/dDC);
      double fou = (dE1/d14/dEC + dD3/d34/dDC);
    */

    double zer ,one ,two ,thr ,fou ;
    double zerx,onex,twox,thrx,foux;
    double zery,oney,twoy,thry,fouy;

    if( data_->regorder == 2 )
    {
        zer = -(1/dAC + 1/dBC + 1/dDC + 1/dEC);
        one =  (dE4/d14/dEC + dA2/d12/dAC);
        two =  (dA1/d12/dAC + dB3/d23/dBC);
        thr =  (dB2/d23/dBC + dD4/d34/dDC);
        fou =  (dE1/d14/dEC + dD3/d34/dDC);

        if( OPERA is_finite( zer+one+two+thr+fou ) )
        {
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
                    vars_->h1->set( cc, cc,     zer );

                else if( order[f] == ccc[0] )
                    vars_->h1->set( cc, ccc[0], one );

                else if( order[f] == ccc[1] )
                    vars_->h1->set( cc, ccc[1], two );

                else if( order[f] == ccc[2] )
                    vars_->h1->set( cc, ccc[2], thr );

                else
                    vars_->h1->set( cc, ccc[3], fou );
            }
        }
    }
    else
    {
        zerx =  1/dEC-1/dBC;
        onex = -dE4/d14/dEC;
        twox =  dB3/d23/dBC;
        thrx =  dB2/d23/dBC;
        foux = -dE1/d14/dEC;

        zery =  1/dAC-1/dDC;
        oney = -dA2/d12/dAC;
        twoy = -dA1/d12/dAC;
        thry =  dD4/d34/dDC;
        fouy =  dD3/d34/dDC;

        if( OPERA is_finite( zerx+onex+twox+thrx+foux+zery+oney+twoy+thry+fouy ) )
        {
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
                {
                    vars_->h1x->set( cc, cc,     zerx );
                    vars_->h1y->set( cc, cc,     zery );
                }
                else if( order[f] == ccc[0] )
                {
                    vars_->h1x->set( cc, ccc[0], onex );
                    vars_->h1y->set( cc, ccc[0], oney );
                }
                else if( order[f] == ccc[1] )
                {
                    vars_->h1x->set( cc, ccc[1], twox );
                    vars_->h1y->set( cc, ccc[1], twoy );
                }
                else if( order[f] == ccc[2] )
                {
                    vars_->h1x->set( cc, ccc[2], thrx );
                    vars_->h1y->set( cc, ccc[2], thry );
                }
                else
                {
                    vars_->h1x->set( cc, ccc[3], foux );
                    vars_->h1y->set( cc, ccc[3], fouy );
                }
            }
        }
    }
}
