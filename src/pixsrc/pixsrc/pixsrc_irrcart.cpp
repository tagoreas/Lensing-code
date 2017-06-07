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



#include "pixsrc_irrcart.hpp"
#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_common_adaptive.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_triangulation.hpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_printer_templates.cpp"
#include <algorithm>
#include <cmath>

pixsrc_irrcart::~pixsrc_irrcart(){}
pixsrc_irrcart::pixsrc_irrcart(PS_SIT imagenumber, inputdata *datavec, commoninputdata *cdata__, PS_SIT tracker_)
{
    MEMORY ps_malloc( &vars_, 1 );
    vars_->tracker = tracker_;
    vars_->imagenumber=imagenumber;
    vars_->pixsrc_class = (void*)this;
    data_  = &datavec[imagenumber];
    cdata_ = cdata__;

    // sets up variables, masks, and does pre-lensing tests
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
        creategrid0();
        createc4c();

        pthread_create(&(vars_->pthreads[0] ), cdata_->attrdetached, COMMONADAPTIVE createc         , vars_->dav);
        pthread_create(&(vars_->pthreads[7] ), cdata_->attrdetached, TRIANGULATION pstriangulate    , vars_->dav);

        COMMON         findsisterimages     (data_, cdata_, vars_);
        COMMONADAPTIVE createlo             (data_, cdata_, vars_);
        COMMON         sourcereconstructions(data_, cdata_, vars_);

        pthread_create(&(vars_->pthreads[8] ), cdata_->attrdetached, COMMON printimageplane         , vars_->dav);
        pthread_create(&(vars_->pthreads[9] ), cdata_->attrdetached, COMMONADAPTIVE printsource     , vars_->dav);
        pthread_create(&(vars_->pthreads[10]), cdata_->attrdetached, COMMONADAPTIVE getmagnification, vars_->dav);
        pthread_create(&(vars_->pthreads[11]), cdata_->attrdetached, STATISTIC evidencecalculator   , vars_->dav);
    }

    COMMON printfinalmessages(data_, cdata_, vars_);
    MEMORY freememory        (data_, cdata_, vars_);
}
void pixsrc_irrcart::flagsrcpixels()
{
    MEMORY ps_malloc( &(vars_->fallsinsrc), vars_->lonc );
    for(PS_SIT r=0; r<vars_->lonc; r++)
        vars_->fallsinsrc[r] = GEOM fallsinsrcinputcircle(
            vars_->srclocall[r*2], vars_->srclocall[r*2+1],
            data_->srcinputcircle, data_->extlengths[2] );
}
void pixsrc_irrcart::creategrid0()
{
    if (vars_->fatalerror)
    {
        return;
    }

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "building grid",
                             cdata_->print2screenmutex);

    if(data_->debug)
    {
        vector< vector<double> > vec(vars_->lonr,vector<double>(3));
        for(PS_SIT cc=0; cc<vars_->lonr; cc++)
        {
            vec[cc][0] = vars_->newloc[vars_->r4rback[cc]*2  ];
            vec[cc][1] = vars_->newloc[vars_->r4rback[cc]*2+1];
            vec[cc][2] = vars_->magnification[vars_->r4rback[cc]];
        }
        PRINTER print <double> ( vars_->tracker, cdata_->basename,
                                 data_->name, true,"magnificationmap.dat",vec, 0, data_->precision, NULL);
    }

    // finding min,max,range
    OPERA assign_p_infinity( &vars_->mindensity );
    OPERA assign_n_infinity( &vars_->maxdensity );
    for(PS_SIT cc=0; cc<vars_->lonr; cc++)
    {
        if(vars_->magnification[vars_->r4rback[cc]]>vars_->maxdensity)
            vars_->maxdensity=vars_->magnification[vars_->r4rback[cc]];
        if(vars_->magnification[vars_->r4rback[cc]]<vars_->mindensity)
            vars_->mindensity=vars_->magnification[vars_->r4rback[cc]];
    }
    vars_->rangedensity = vars_->maxdensity-vars_->mindensity;


    vars_->zeroethgridsize = std::max(vars_->rangex,vars_->rangey) + CONSTANT smallnumber;

    /*
      vars_->levelmag1=1;
      while(std::pow(vars_->zeroethgridsize/std::pow(2,vars_->levelmag1+1),2)>vars_->match2mag1)
      vars_->levelmag1++;
      if(std::pow(2,vars_->levelmag1)>vars_->zeroethgridsize)
      vars_->zeroethgridsize=std::pow(2,vars_->levelmag1);
      vars_->startlevel=vars_->levelmag1;
      while(vars_->startlevel-1>0)
      if(0.5*vars_->mindensity<vars_->match2mag1*std::pow(4,vars_->startlevel-vars_->levelmag1))
      vars_->startlevel--;
      else
      break;
    */

    //vars_->levelmag1 = std::floor(2*std::log(vars_->zeroethgridsize)/std::log(4)-1-data_->levelshift);
    vars_->levelmag1 = (1+(data_->levelshift-2))/vars_->zeroethgridsize/vars_->zeroethgridsize;

    vars_->maxgridlevel=0;
    while(!OPERA equalszero(vars_->zeroethgridsize/std::pow(2.,vars_->maxgridlevel+1.)))
        vars_->maxgridlevel++;
}
void pixsrc_irrcart::createc4c()
{
    if (vars_->fatalerror)
    {
        return;
    }

    optimizegrid();
    flagsrcpixels();

    MEMORY ps_malloc( &(vars_->c4c    ), vars_->lonc );
    MEMORY ps_malloc( &(vars_->c4cback), vars_->lonc );
    PS_SIT lonctemp = vars_->lonc;
    vars_->lonc=0;
    PS_SIT unused=0;
    for( PS_SIT s=0; s<lonctemp; ++s )
    {
        if(vars_->fallsinsrc[s])
        {
            vars_->c4c[s] = s - unused;
            vars_->c4cback[vars_->lonc++] = s;
        }
        else
        {
            vars_->c4c[s]=-1;
            ++unused;
        }
    }




    // here, we triangulate the poitns and slightly move points on edge of grid
    // so that triangle.c doesn't freak out

    // examining edges
    vars_->triin->numberofpoints = vars_->lonc;
    MEMORY ps_malloc( &(vars_->triin->pointlist), vars_->lonc*2 );
    //PRINTER printtristruct(vars_->triin,  cdata_->print2screenmutex);
    //PRINTER printtristruct(vars_->triout, cdata_->print2screenmutex);
    for(PS_SIT g=0; g<vars_->lonc; g++)
    {
        vars_->triin->pointlist[g*2  ] = vars_->srclocall[vars_->c4cback[g]*2  ];
        vars_->triin->pointlist[g*2+1] = vars_->srclocall[vars_->c4cback[g]*2+1];
    }

    PS_SIT numverts;
    double *convexhull;
    GEOM getconvexhull( vars_->triin, vars_->triout, &numverts,
                        &convexhull, NULL );
    //PRINTER printtristruct(vars_->triin,  cdata_->print2screenmutex);
    //PRINTER printtristruct(vars_->triout, cdata_->print2screenmutex);
    double centroidx, centroidy;
    GEOM getcentroid( numverts, convexhull, &centroidx, &centroidy );
    MEMORY ps_free( convexhull );
    PS_SIT *gotindex;
    MEMORY ps_malloc( &(gotindex), numverts );
    PS_SIT runningindex=0;
    char skip;
    double x, y, r, cos, sin;
    double vertices[6];
    PS_SIT pvertex, ppvertex;
    char isbadtri;
    char isedge[3];

    // check for edges the triangulation missed
    // these are signaled by very skinny triangles
    for(PS_SIT tri=0; tri<vars_->triout->numberoftriangles; ++tri)
    {

        for( PS_SIT v=0; v<6; ++v )
            vertices[v] = vars_->triin->pointlist[vars_->triout->trianglelist[tri*3+v/2]*2+v%2];

        pvertex  = 2;
        ppvertex = 1;
        isbadtri = 0;
        for(PS_SIT tvertex=0; tvertex<3; ++tvertex)
        {
            // if triangle is skinny
            if(
                GEOM angle(
                    vertices[ppvertex*2],vertices[ppvertex*2+1],
                    vertices[ pvertex*2],vertices[ pvertex*2+1],
                    vertices[ tvertex*2],vertices[ tvertex*2+1]
                    )
                >= CONSTANT smalltriangleanglecutoff
                )
            {
                isbadtri = 1;
                break;
            }
            ppvertex = pvertex;
            pvertex  = tvertex;
        }

        if( isbadtri )
        {
            std::fill( isedge, isedge+3, 0 );
            for( PS_SIT v=0; v<3; ++v )
            {
                for( PS_SIT j=0; j<vars_->triout->numberofedges; ++j )
                {
                    if( vars_->triout->edgemarkerlist[j] == 1 &&
                        ( vars_->triout->edgelist[j*2  ] == vars_->triout->trianglelist[tri*3+v] ||
                          vars_->triout->edgelist[j*2+1] == vars_->triout->trianglelist[tri*3+v] ) )
                    {
                        isedge[v] = 1;
                        break;
                    }
                }
            }
            for( PS_SIT v=0; v<3; ++v )
            {
                if( !isedge[v] )
                {
                    x = vars_->triin->pointlist[vars_->triout->trianglelist[tri*3+v]*2  ] -
                        centroidx;
                    y = vars_->triin->pointlist[vars_->triout->trianglelist[tri*3+v]*2+1] -
                        centroidy;
                    r = std::sqrt( x*x + y*y );
                    cos = x/r;
                    sin = y/r;

                    vars_->triin->pointlist[vars_->triout->trianglelist[tri*3+v]*2  ] +=
                        cos;
                    vars_->triin->pointlist[vars_->triout->trianglelist[tri*3+v]*2+1] +=
                        sin;
                }
            }
        }
    }

    for( PS_SIT j=0; j<vars_->triout->numberofedges; ++j )
    {
        if( vars_->triout->edgemarkerlist[j] == 1 )
        {
            for( PS_SIT i=0; i<2; ++i )
            {
                skip = 0;
                for( PS_SIT v=0; v<runningindex; ++v )
                {
                    if( vars_->triout->edgelist[j*2+i] == gotindex[v] )
                    {
                        skip = 1;
                        break;
                    }
                }
                if( skip )
                {
                    continue;
                }

                x = vars_->triin->pointlist[vars_->triout->edgelist[j*2+i]*2  ] -
                    centroidx;
                y = vars_->triin->pointlist[vars_->triout->edgelist[j*2+i]*2+1] -
                    centroidy;
                r = std::sqrt( x*x + y*y );
                cos = x/r;
                sin = y/r;

                vars_->triin->pointlist[vars_->triout->edgelist[j*2+i]*2  ] +=
                    cos;
                vars_->triin->pointlist[vars_->triout->edgelist[j*2+i]*2+1] +=
                    sin;

                gotindex[runningindex++] = vars_->triout->edgelist[j*2+i];
            }
        }
    }
    MEMORY ps_free( gotindex );

    //PRINTER printtristruct(vars_->triin,  cdata_->print2screenmutex);
    //PRINTER printtristruct(vars_->triout, cdata_->print2screenmutex);

    // store perturbed vertices
    vars_->need2retriangulate = 1;
    double *temp = vars_->triin->pointlist;
    vars_->triin->pointlist = NULL;
    MEMORY triangulatestructdestruct(   vars_->triin      );
    MEMORY triangulatestructdestruct(   vars_->triout     );
    MEMORY ps_free                  (   vars_->triin      );
    MEMORY ps_free                  (   vars_->triout     );
    MEMORY ps_malloc                ( &(vars_->triin ), 1 );
    MEMORY ps_malloc                ( &(vars_->triout), 1 );
    MEMORY triangulatestructinit    (   vars_->triin      );
    MEMORY triangulatestructinit    (   vars_->triout     );
    vars_->triin->pointlist = temp;
    vars_->triin->numberofpoints = vars_->lonc;

    //PRINTER printtristruct(vars_->triin,  cdata_->print2screenmutex);
    //PRINTER printtristruct(vars_->triout, cdata_->print2screenmutex);

    if( data_->verbose != 1 )
        PRINTER print2screen(data_->print2screenname,
                             "done building grid",
                             cdata_->print2screenmutex);
}



//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////// code below is all adaptive Cartesian gridding ///////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



void pixsrc_irrcart::addsubgrid(vector<PS_SIT> &path)
{
    PS_SIT scanuptoinarrayofbsarrays=0;
    for(PS_unsignedSIT j=0; j<path.size(); j++)
    {
        if(!(*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][0])
        {
            (*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][0]=true;
            if(j==path.size()-1 && !(*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][1])
            {
                (*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][1]=true;
                vars_->lonc+=4;
            }

            PS_SIT indexofnullbs;

            if(vars_->grid->size()==j+1)
            {
                vars_->grid->push_back(vector< vector< vector<bool> > >());
                indexofnullbs=0;
            }
            else
                indexofnullbs = indexofnext(j,path,scanuptoinarrayofbsarrays);

            (*(vars_->grid))[j+1].insert( (*(vars_->grid))[j+1].begin()+indexofnullbs,
                                          vector< vector<bool> >
                                          (4,vector<bool>(2,false))                           );

            scanuptoinarrayofbsarrays = indexofnullbs;
        }
        else if(j==path.size()-1 && !(*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][1])
        {
            (*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][1]=true;
            vars_->lonc+=4;
        }
        else if(j!=path.size()-1)
            scanuptoinarrayofbsarrays = indexofnext(j,path,scanuptoinarrayofbsarrays);
    }
}
void pixsrc_irrcart::getpathfrompositioninarrays(PS_SIT *pos, vector<PS_SIT> &result)
{
    result.resize(pos[0]+1);
    result[pos[0]] = pos[2];

    PS_SIT targetIndex=pos[1];
    PS_SIT level = pos[0]-1;
    while(1)
    {
    outerloop:
        if(!(level>=0))
            break;

        PS_SIT nextIndex=0;
        for(PS_unsignedSIT index=0; index<(*(vars_->grid))[level].size(); index++)
            for(PS_SIT posInBitArray=0; posInBitArray<4; posInBitArray++)
                if((*(vars_->grid))[level][index][posInBitArray][0])
                {
                    if(nextIndex == targetIndex)
                    {
                        targetIndex=index;
                        result[level]=posInBitArray;

                        // rigging a continue outer loop statement
                        level--;
                        goto outerloop;
                    }
                    nextIndex++;
                }
        level--;
    }
}
void pixsrc_irrcart::getpositioninarraysfrompath(vector<PS_SIT> &path, PS_SIT *result)
{
    PS_SIT scanuptoinarrayofbsarrays=0;

    for(PS_unsignedSIT j=0; j<path.size()-1; j++)
        scanuptoinarrayofbsarrays = indexofnext(j,path,scanuptoinarrayofbsarrays);

    result[0] = path.size()-1;
    result[1] = scanuptoinarrayofbsarrays;
    result[2] = path[path.size()-1];
}
PS_SIT pixsrc_irrcart::indexofnext(PS_SIT subgridLevel, vector<PS_SIT> &path, PS_SIT stopAt)
{
    PS_SIT indexofnext=0;
    for(PS_SIT j2=0; j2<=stopAt; j2++)
    {
        PS_SIT scanuptoinbitsetarray=4;
        if(j2==stopAt)
            scanuptoinbitsetarray=path[subgridLevel];

        for(PS_SIT k=0; k<scanuptoinbitsetarray; k++)
            if((*(vars_->grid))[subgridLevel][j2][k][0])
                indexofnext++;
    }
    return indexofnext;
}
bool pixsrc_irrcart::nodeexists(vector<PS_SIT> &path)
{
    if(path[0]==-1)
        return false;
    PS_SIT scanuptoinarrayofbsarrays=0;
    for(PS_unsignedSIT j=0; j<path.size(); j++)
    {
        try
        {
            if(!(*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][0])
                return false;
            if(j==path.size()-1 && !(*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][1])
                return false;
        }
        catch(PS_SIT e)
        {
            return false;
        }
        scanuptoinarrayofbsarrays=indexofnext(j,path,scanuptoinarrayofbsarrays);
    }
    return true;
}
bool pixsrc_irrcart::nodepossiblyexistshereorbelow(vector<PS_SIT> &path)
{
    if(path[0]==-1)
        return false;
    PS_SIT scanuptoinarrayofbsarrays=0;
    for(PS_unsignedSIT j=0; j<path.size(); j++)
    {
        try
        {
            if(!(*(vars_->grid))[j][scanuptoinarrayofbsarrays][path[j]][0])
                return false;
        }
        catch(PS_SIT e)
        {
            return false;
        }
        scanuptoinarrayofbsarrays=indexofnext(j,path,scanuptoinarrayofbsarrays);
    }
    return true;
}
void pixsrc_irrcart::optimizegrid()
{
    vars_->lonc=0;
    makegrid();
/*
  for(PS_unsignedSIT level=0; level<vars_->grid->size(); level++)
  {
  for(PS_unsignedSIT m=0; m<(*(vars_->grid))[level].size(); m++)
  {
  for(PS_SIT m2=0; m2<4; m2++)
  cout << (*(vars_->grid))[level][m][m2][0];
  cout << " ";

  }
  cout << endl;
  for(PS_unsignedSIT m=0; m<(*(vars_->grid))[level].size(); m++)
  {
  for(PS_SIT m2=0; m2<4; m2++)
  cout << (*(vars_->grid))[level][m][m2][1];
  cout << " ";

  }
  cout << endl;
  cout << endl;
  }
*/
    setsrclocall();
}
void pixsrc_irrcart::setsrclocall()
{
    PS_SIT position=0;
    MEMORY ps_malloc( &(vars_->srclocall), vars_->lonc*2 );

    vars_->gridpointer->clear();
    vars_->gridpointer->resize(vars_->grid->size());
    for(PS_unsignedSIT level=0; level<vars_->grid->size(); level++)
    {
        (*(vars_->gridpointer))[level].resize((*(vars_->grid))[level].size());
        for(PS_unsignedSIT m=0; m<(*(vars_->grid))[level].size(); m++)
        {
            (*(vars_->gridpointer))[level][m].resize(4);
            for(PS_SIT m2=0; m2<4; m2++)
            {
                if((*(vars_->grid))[level][m][m2][1])
                {
                    (*(vars_->gridpointer))[level][m][m2].resize(4);
                    std::fill((*(vars_->gridpointer))[level][m][m2].begin(),
                              (*(vars_->gridpointer))[level][m][m2].end(),-1);
                }
                else
                    (*(vars_->gridpointer))[level][m][m2].resize(0);
            }
        }
    }

    for(PS_SIT level=0; level<(PS_SIT)vars_->gridpointer->size(); level++)
        for(PS_SIT m=0; m<(PS_SIT)(*(vars_->gridpointer))[level].size(); m++)
            for(PS_SIT m2=0; m2<4; m2++)
            {
                PS_SIT poss[3] = {level,m,m2};
                vector<PS_SIT> path;
                getpathfrompositioninarrays(poss, path);
                PS_SIT length = path.size();
                double posAtCtr[2];
                positionofnode(path, posAtCtr);
                double finSpacing = vars_->zeroethgridsize/std::pow(2.,length+1.);

                for(PS_unsignedSIT v=0; v<(*(vars_->gridpointer))[level][m][m2].size(); v++)
                    if((*(vars_->gridpointer))[level][m][m2][v]==-1)
                    {
                        double posOfVert[2];
                        switch(v)
                        {
                        case 0: posOfVert[0] = posAtCtr[0] + finSpacing;
                            posOfVert[1] = posAtCtr[1] - finSpacing;
                            break;
                        case 1: posOfVert[0] = posAtCtr[0] - finSpacing;
                            posOfVert[1] = posAtCtr[1] - finSpacing;
                            break;
                        case 2: posOfVert[0] = posAtCtr[0] - finSpacing;
                            posOfVert[1] = posAtCtr[1] + finSpacing;
                            break;
                        case 3: posOfVert[0] = posAtCtr[0] + finSpacing;
                            posOfVert[1] = posAtCtr[1] + finSpacing;
                            break;
                        default:PRINTER printerror(data_->print2screenname,
                                                   "accessing nonexistent quadrant",
                                                   cdata_->print2screenmutex);
                            // this is just to suppress compiler warning
                            posOfVert[0] = posOfVert[1] = 0;
                            break;
                        }
                        vars_->srclocall[position*2  ] = posOfVert[0];
                        vars_->srclocall[position*2+1] = posOfVert[1];

                        vector< vector<PS_SIT> > surrPath;
                        findsurroundingboxes(path,v, surrPath);

                        // START examine same level
                        for(PS_SIT b=0; b<4; b++)
                            if(nodeexists(surrPath[b]))
                            {
                                //printPath(surrPath[b]);
                                PS_SIT pos[3];
                                getpositioninarraysfrompath(surrPath[b], pos);
                                PS_SIT vNew=b-2;
                                if(vNew<0)
                                    vNew+=4;

                                (*(vars_->gridpointer))[pos[0]][pos[1]][pos[2]][vNew] = position;
                            }
                        // END examine same level

                        // START examine sublevels
                        for(PS_SIT b=0; b<4; b++)
                        {
                            PS_SIT vNew=b-2;
                            if(vNew<0)
                                vNew+=4;
                            vector<PS_SIT> sublevel(surrPath[b].begin(),surrPath[b].begin()+length);
                            PS_SIT extra=1;
                            while(nodepossiblyexistshereorbelow(sublevel))
                            {
                                sublevel.resize(length+extra);
                                for(PS_SIT h=0; h<extra; h++)
                                    sublevel[length+h]=vNew;
                                if(nodeexists(sublevel))
                                {
                                    PS_SIT pos[3];
                                    getpositioninarraysfrompath(sublevel, pos);
                                    (*(vars_->gridpointer))[pos[0]][pos[1]][pos[2]][vNew] = position;
                                }
                                extra++;
                            }
                        }
                        // END examine sublevels

                        // START examine prelevels
                        for(PS_SIT b=0; b<4; b++)
                        {
                            PS_SIT vNew=b-2;
                            if(vNew<0)
                                vNew+=4;
                            PS_SIT preLength=length;
                            vector<PS_SIT> preLevel(surrPath[b].begin(),surrPath[b].begin()+preLength);
                            while(preLevel[preLength-1]==vNew)
                            {
                                preLength--;
                                if(preLength==0)
                                    break;
                                preLevel.resize(preLength);

                                if(nodeexists(preLevel))
                                {
                                    PS_SIT pos[3];
                                    getpositioninarraysfrompath(preLevel, pos);
                                    (*(vars_->gridpointer))[pos[0]][pos[1]][pos[2]][vNew] = position;
                                }
                            }
                        }
                        // END examine prelevels
                        position++;
                    }
            }
}
void pixsrc_irrcart::makegrid()
{
    vars_->grid->insert( vars_->grid->begin(), vector< vector< vector<bool> > >
                         ( 1,vector< vector<bool> >
                           ( 4,vector<bool>
                             ( 2,false ))));

    for(PS_SIT cc=0; cc<vars_->lonr; cc++)
        testthisnodewithgivendensityandaddtogrid(vars_->newloc[vars_->r4rback[cc]*2  ],
                                                 vars_->newloc[vars_->r4rback[cc]*2+1],
                                                 vars_->magnification[vars_->r4rback[cc]] );

    testforoverlappinggridpoints();
}
void pixsrc_irrcart::testforoverlappinggridpoints()
{
    vars_->gridpointer->resize(vars_->grid->size());
    for(PS_unsignedSIT level=0; level<vars_->grid->size(); level++)
    {
        (*(vars_->gridpointer))[level].resize((*(vars_->grid))[level].size());
        for(PS_unsignedSIT m=0; m<(*(vars_->grid))[level].size(); m++)
        {
            (*(vars_->gridpointer))[level][m].resize(4);
            for(PS_SIT m2=0; m2<4; m2++)
            {
                if((*(vars_->grid))[level][m][m2][1])
                {
                    (*(vars_->gridpointer))[level][m][m2].resize(4,-1);
                    std::fill((*(vars_->gridpointer))[level][m][m2].begin(),
                              (*(vars_->gridpointer))[level][m][m2].end(),-1);
                }
                else
                    (*(vars_->gridpointer))[level][m][m2].resize(0);
            }
        }
    }

    for(PS_SIT level=0; level<(PS_SIT)vars_->gridpointer->size(); level++)
        for(PS_SIT m=0; m<(PS_SIT)(*(vars_->gridpointer))[level].size(); m++)
            for(PS_SIT m2=0; m2<4; m2++)
            {
                PS_SIT poss[3] = {level,m,m2};
                vector<PS_SIT> path;
                getpathfrompositioninarrays(poss, path);
                PS_SIT length = path.size();

                /*
                  if(length!=2)
                  continue;
                  if(path[0]!=2 || path[1]!=2)
                  continue;
                */

                for(PS_unsignedSIT v=0; v<(*(vars_->gridpointer))[level][m][m2].size(); v++)
                    if((*(vars_->gridpointer))[level][m][m2][v]==-1)
                    {
                        vector< vector<PS_SIT> > surrPath;
                        findsurroundingboxes(path,v, surrPath);

                        // START examine same level
                        for(PS_SIT b=0; b<4; b++)
                            if(nodeexists(surrPath[b]))
                            {
                                PS_SIT pos[3];
                                getpositioninarraysfrompath(surrPath[b], pos);
                                PS_SIT vNew=b-2;
                                if(vNew<0)
                                    vNew+=4;

                                vars_->lonc--;
                                (*(vars_->gridpointer))[pos[0]][pos[1]][pos[2]][vNew] = 0;
                            }
                        vars_->lonc++;
                        // END examine same level

                        // START examine sublevels
                        for(PS_SIT b=0; b<4; b++)
                        {
                            PS_SIT vNew=b-2;
                            if(vNew<0)
                                vNew+=4;

                            vector<PS_SIT> sublevel(surrPath[b].begin(),surrPath[b].begin()+length);
                            PS_SIT extra=1;
                            while(nodepossiblyexistshereorbelow(sublevel))
                            {
                                sublevel.resize(length+extra);
                                for(PS_SIT h=0; h<extra; h++)
                                    sublevel[length+h]=vNew;

                                if(nodeexists(sublevel))
                                {
                                    PS_SIT pos[3];
                                    getpositioninarraysfrompath(sublevel, pos);

                                    vars_->lonc--;
                                    (*(vars_->gridpointer))[pos[0]][pos[1]][pos[2]][vNew] = 0;
                                }
                                extra++;
                            }
                        }
                        // END examine sublevels

                        // START examine prelevels
                        for(PS_SIT b=0; b<4; b++)
                        {
                            PS_SIT vNew=b-2;
                            if(vNew<0)
                                vNew+=4;
                            PS_SIT preLength=length;
                            vector<PS_SIT> preLevel(surrPath[b].begin(),surrPath[b].begin()+preLength);

                            while(preLevel[preLength-1]==vNew)
                            {
                                preLength--;
                                if(preLength==0)
                                    break;
                                preLevel.resize(preLength);

                                if(nodeexists(preLevel))
                                {
                                    PS_SIT pos[3];
                                    getpositioninarraysfrompath(preLevel, pos);

                                    vars_->lonc--;
                                    (*(vars_->gridpointer))[pos[0]][pos[1]][pos[2]][vNew] = 0;
                                }
                            }
                        }
                        // END examine prelevels
                    }
            }
}
void pixsrc_irrcart::findsurroundingboxes(vector<PS_SIT> &path, PS_SIT lastIndex, vector< vector<PS_SIT> > &result)
{
    PS_SIT length = path.size();
    result.resize(4,vector<PS_SIT>(length));

    if(fabs(lastIndex-path[length-1])==2)
        for(PS_SIT b=0; b<4; b++)
        {
            std::copy(path.begin(),path.end()-1,result[b].begin());
            result[b][length-1] = b;
        }
    else
    {
        PS_SIT maxIncrement = OPERA round(std::pow(2.,length-1.));
        PS_SIT distX=0, distY=0;
        for(PS_SIT b=0; b<length; b++)
        {
	  PS_SIT increment=maxIncrement/OPERA round(std::pow(2.,b+0.));
            if(path[b]==0 || path[b]==3)
                distX+=increment;
            if(path[b]==2 || path[b]==3)
                distY+=increment;
        }

        PS_SIT distSurr[4][2];
        for(PS_SIT b=0; b<4; b++)
        {
            distSurr[b][0]=distX;
            distSurr[b][1]=distY;
        }
        switch(lastIndex)
        {
        case 0: distSurr[0][0]++;
            distSurr[0][1]--;
            distSurr[1][1]--;
            distSurr[3][0]++;
            break;
        case 1: distSurr[0][1]--;
            distSurr[1][0]--;
            distSurr[1][1]--;
            distSurr[2][0]--;
            break;
        case 2: distSurr[1][0]--;
            distSurr[2][0]--;
            distSurr[2][1]++;
            distSurr[3][1]++;
            break;
        case 3: distSurr[0][0]++;
            distSurr[2][1]++;
            distSurr[3][0]++;
            distSurr[3][1]++;
            break;
        default: ;//PRINTER printerror("trying to find surrounding boxes. ");
        }

        bool pseudoPaths[4][length][2];
        for(PS_SIT g=0; g<4; g++)
            for(PS_SIT m=0; m<length; m++)
                std::fill(pseudoPaths[g][m],pseudoPaths[g][m]+2,false);

        for(PS_SIT b=0; b<4; b++)
        {
            if(distSurr[b][0]<0 || distSurr[b][0]>=2*maxIncrement ||
               distSurr[b][1]<0 || distSurr[b][1]>=2*maxIncrement)
            {
                result[b][0]=-1;
                continue;
            }
            if(distSurr[b][0]==distX&&distSurr[b][1]==distY)
            {
                result[b]=path;
                continue;
            }

            for(PS_SIT level=0; level<length; level++)
            {
                PS_SIT increment=maxIncrement/OPERA round(std::pow(2.,level+0.));
                if(distSurr[b][0]>=increment)
                {
                    distSurr[b][0]-=increment;
                    pseudoPaths[b][level][0]=true;
                }
                if(distSurr[b][1]>=increment)
                {
                    distSurr[b][1]-=increment;
                    pseudoPaths[b][level][1]=true;
                }
            }
            for(PS_SIT level=0; level<length; level++)
            {
                if(pseudoPaths[b][level][0] && pseudoPaths[b][level][1])
                    result[b][level]=3;
                else if(pseudoPaths[b][level][0] && !pseudoPaths[b][level][1])
                    result[b][level]=0;
                else if(!pseudoPaths[b][level][0] && pseudoPaths[b][level][1])
                    result[b][level]=2;
                else
                    result[b][level]=1;
            }
        }
    }
}
void pixsrc_irrcart::testthisnodewithgivendensityandaddtogrid(double x, double y, double density)
{
    //double z2 = vars_->zeroethgridsize*vars_->zeroethgridsize;
    PS_SIT size=0;

    while(1)
    {
        //if(density>OPERA round(std::pow(4.0,size-vars_->levelmag1)/z2))
        if( density >= std::pow(4.0,size+0.) * vars_->levelmag1 )
            size++;
        else
            break;
    }


    if( size > vars_->maxgridlevel )
        size = vars_->maxgridlevel;

    if( !size )
        size = 1;

    vector<PS_SIT> path;
    path.resize(size);
    vector<PS_SIT> pathTemp;

    for(PS_SIT level=0; level<size; level++)
    {
        pathTemp.resize(level);
        std::copy(path.begin(),path.begin()+level,pathTemp.begin());
        double position[2];
        positionofnode(pathTemp, position);

        if(x>position[0])
        {
            if(y>position[1])
                path[level]=3;
            else
                path[level]=0;
        }
        else if(x<position[0])
        {

            if(y<position[1])
                path[level]=1;
            else
                path[level]=2;
        }
        else
        {
            if(y>position[1])
                path[level]=3;
            else if(y<position[1])
                path[level]=1;
            else
                path[level]=0;
        }
    }

    if( OPERA is_exactly_p_inf( data_->srcinputcircle[0][2] ) )
    {
        addsubgrid(path);
    }
    else
    {
        double position[2];
        positionofnode(path, position);
        double spacing = vars_->zeroethgridsize/std::pow(2.,path.size()+1.);
        if( GEOM fallsinsrcinputcircle(
                position[0]+spacing,position[1]-spacing,
                data_->srcinputcircle,data_->extlengths[2]) ||
            GEOM fallsinsrcinputcircle(
                position[0]-spacing,position[1]-spacing,
                data_->srcinputcircle,data_->extlengths[2]) ||
            GEOM fallsinsrcinputcircle(
                position[0]-spacing,position[1]+spacing,
                data_->srcinputcircle,data_->extlengths[2]) ||
            GEOM fallsinsrcinputcircle(
                position[0]+spacing,position[1]+spacing,
                data_->srcinputcircle,data_->extlengths[2]) )
            addsubgrid(path);
    }
}
void pixsrc_irrcart::positionofnode(vector<PS_SIT> &path, double *result)
{
    result[0] = vars_->minx+vars_->xoffset;
    result[1] = vars_->miny+vars_->yoffset;

    for(PS_unsignedSIT m=0; m<=path.size(); m++)
    {
        double spacing = vars_->zeroethgridsize/std::pow(2.,m+1.);


        if(m==path.size())
        {
            result[0]+=spacing;
            result[1]+=spacing;
        }
        else
        {
            if(path[m]==0)
                result[0] += spacing;
            else if(path[m]==2)
                result[1] += spacing;
            else if(path[m]==3)
            {
                result[0] += spacing;
                result[1] += spacing;
            }
        }
    }
}
