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
#include "pixsrc_geometry.hpp"
#include "pixsrc_memory_templates.cpp"
#include <cmath>

// see
// "Development of Irregular-Grid Finite Difference Method (IFDM)
// for Governing Equations in Strong Form"
// by
// GEORGE XU @ Institute of High Performance Computing
// for a generalization of the method here

// the derivatives here are constructed on a grid where
// right points in an increasing direction
// and down points in an increasing direction

void pixsrc_regularization::greengaussgetder2(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    PS_SIT dim = 1;
    PS_SIT waitlist[1] = {7};
    OPERA pthreadswait(vars_,dim,waitlist);

    // contains deltax,deltay,reference position(midpoint) for 3 vertices of each triangle
    // it connected the point at [tri][vertex] to the vertex preceding it.. -1 wraps around to 2
    double **dS;
    MEMORY ps_malloc( &(dS), vars_->triout->numberoftriangles, 12 );
    char isbadtri;
    PS_SIT pvertex, ppvertex, nextvertex;
    double centroidx, centroidy, vertices[6];
    for(PS_SIT tri=0; tri<vars_->triout->numberoftriangles; ++tri)
    {
        for( PS_SIT v=0; v<6; ++v )
            vertices[v] = vars_->triout->pointlist[vars_->triout->trianglelist[tri*3+v/2]*2+v%2];

        // here we attempt to find bad triangles
        isbadtri = 0;
        pvertex  = 2;
        ppvertex = 1;
        for(PS_SIT tvertex=0; tvertex<3; ++tvertex)
        {
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

        if(isbadtri)
        {
            std::fill(dS[tri],dS[tri]+12,0);
            continue;
        }

        centroidx = 0;
        centroidy = 0;
        for(PS_SIT vertex=0; vertex<3; ++vertex)
        {
            centroidx += vertices[vertex*2  ];
            centroidy += vertices[vertex*2+1];
        }
        centroidx /= 3.0;
        centroidy /= 3.0;

        nextvertex=2;
        for(PS_SIT vertex=0; vertex<3; ++vertex)
        {
            dS[tri][vertex*4+2] = (vertices[vertex*2  ] + vertices[nextvertex*2  ]) / 2.0;
            dS[tri][vertex*4+3] = (vertices[vertex*2+1] + vertices[nextvertex*2+1]) / 2.0;
            dS[tri][vertex*4  ] = centroidx - dS[tri][vertex*4+2];
            dS[tri][vertex*4+1] = centroidy - dS[tri][vertex*4+3];

            nextvertex = vertex;
        }
    }

    // figure out to which triangles each point belongs
    vector< vector<PS_SIT> > surrtrilist(vars_->triout->numberofpoints,vector<PS_SIT>());
    for(PS_SIT t=0; t<vars_->triout->numberoftriangles*3; ++t)
    {
        surrtrilist[vars_->triout->trianglelist[t]].push_back(t/3);
    }

    double domainareas[vars_->triout->numberofpoints];
    std::fill( domainareas,domainareas+vars_->triout->numberofpoints,0 );
    double **dSsums = (double**)malloc(vars_->triout->numberofpoints*sizeof(double*));

    // these replace h1x and h1y here
    vector< vector<double> > derx(vars_->triout->numberofpoints,vector<double>());
    vector< vector<double> > dery(vars_->triout->numberofpoints,vector<double>());
    double tempx,tempy;
    // compute derivative
    for(PS_SIT point=0; point<vars_->triout->numberofpoints; ++point)
    {
        PS_SIT numsurrtri = (PS_SIT)surrtrilist[point].size();
        dSsums[point] = (double*)malloc(6*numsurrtri*sizeof(double));
        std::fill( dSsums[point],dSsums[point]+6*numsurrtri,0 );

        // first index hold index of point. second holds value
        double weights[surrtrilist[point].size()*16];

        for(PS_SIT tri=0; tri<(PS_SIT)surrtrilist[point].size(); ++tri)
        {
            isbadtri=1;
            for(PS_SIT j=0; j<12;++j)
            {
                if(dS[surrtrilist[point][tri]][j])
                {
                    isbadtri = 0;
                    break;
                }
            }

            if(isbadtri)
            {
                for(PS_SIT j=1; j<16; j+=2)
                    weights[tri*16+j] = 0;
                continue;
            }

            PS_SIT weightindex = 0;

            // vertex2 and vertex3 connect to the central point (vertex1)
            // vertex1 and vertex2 are to be referenced in dS
            PS_SIT vertex1 = -1;
            for(PS_SIT vertex=0; vertex<3; ++vertex)
            {
                if(vars_->triout->trianglelist[surrtrilist[point][tri]*3+vertex] == point)
                {
                    vertex1=vertex;
                    break;
                }
            }

            PS_SIT vertex2 = (vertex1<2) ? vertex1+1 : 0;
            PS_SIT vertex3 = (vertex1>0) ? vertex1-1 : 2;

            // central point is shifted to coordinate systems where midpoints are origin
            double centralshifted[4];
            centralshifted[0] = vars_->triout->pointlist[point*2  ]-dS[surrtrilist[point][tri]][4*vertex1+2];
            centralshifted[1] = vars_->triout->pointlist[point*2+1]-dS[surrtrilist[point][tri]][4*vertex1+3];
            centralshifted[2] = vars_->triout->pointlist[point*2  ]-dS[surrtrilist[point][tri]][4*vertex2+2];
            centralshifted[3] = vars_->triout->pointlist[point*2+1]-dS[surrtrilist[point][tri]][4*vertex2+3];

            // this block gets vertex1 x
            {
                // now we find out which way the normal unit vector points
                double negatefacx = 1.0;
                double negatefacy = 1.0;
                double slope     = -OPERA slope(0.0,0.0,
                                                dS[surrtrilist[point][tri]][4*vertex1  ],
                                                dS[surrtrilist[point][tri]][4*vertex1+1]);

                // if the line formed by this dS is above central point
                if( centralshifted[1] < slope*centralshifted[0] )
                {
                    if(slope > 0.0)
                        negatefacx = -1.0;
                }
                else
                {
                    negatefacy = -1.0;
                    if(slope < 0.0)
                        negatefacx = -1.0;
                }

                tempx = negatefacx * 0.5 * std::fabs(dS[surrtrilist[point][tri]][vertex1*4+1]);
                tempy = negatefacy * 0.5 * std::fabs(dS[surrtrilist[point][tri]][vertex1*4  ]);

                dSsums[point][tri*6+vertex3*2  ] += tempx;
                dSsums[point][tri*6+vertex3*2+1] += tempy;

                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                weights[tri*16+weightindex++] = tempx;
                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex3];
                weights[tri*16+weightindex++] = tempx;

                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                weights[tri*16+weightindex++] = tempy;
                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex3];
                weights[tri*16+weightindex++] = tempy;
            }

            // this block gets vertex2 x
            {
                // now we find out which way the normal unit vector points
                double negatefacx = 1.0;
                double negatefacy = 1.0;
                double slope     = -OPERA slope(0.0,0.0,
                                                dS[surrtrilist[point][tri]][4*vertex2],
                                                dS[surrtrilist[point][tri]][4*vertex2+1]);

                // if the line formed by this dS is above central point
                if( centralshifted[3] < slope*centralshifted[2] )
                {
                    if(slope > 0.0)
                        negatefacx = -1.0;
                }
                else
                {
                    negatefacy = -1.0;
                    if(slope < 0.0)
                        negatefacx = -1.0;
                }

                tempx = negatefacx * 0.5 * std::fabs(dS[surrtrilist[point][tri]][vertex2*4+1]);
                tempy = negatefacy * 0.5 * std::fabs(dS[surrtrilist[point][tri]][vertex2*4  ]);

                dSsums[point][tri*6+vertex2*2  ] += tempx;
                dSsums[point][tri*6+vertex2*2+1] += tempy;

                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex2];
                weights[tri*16+weightindex++] = tempx;
                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                weights[tri*16+weightindex++] = tempx;

                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex2];
                weights[tri*16+weightindex++] = tempy;
                weights[tri*16+weightindex++] = vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                weights[tri*16+weightindex++] = tempy;
            }

            // this is not ideal but ..
            // the 4 points that form the quadrilateral that overlaps with the
            // domain area for this triangle have coordinates here with their
            // origins at the midpoint connected to vertex1. So here are those
            // coordinate transformations ..
            double coords[8] =
                {
                    centralshifted[0],centralshifted[1],
                    centralshifted[0]-centralshifted[2],centralshifted[1]-centralshifted[3],
                    dS[surrtrilist[point][tri]][4*vertex1],dS[surrtrilist[point][tri]][4*vertex1+1],
                    0.0,0.0
                };
            domainareas[point] += GEOM areapoly(coords,4);

            // this here prints out grid points and assosiated centroids
            // I stream-edited this to make a ds9 reg file.
            //              cout << "positions " << vars_->triout->pointlist[point*2] << " " << vars_->triout->pointlist[point*2+1] << " " << dS[surrtrilist[point][tri]][4*vertex1+2]+dS[surrtrilist[point][tri]][4*vertex1] << " " << dS[surrtrilist[point][tri]][4*vertex1+3]+dS[surrtrilist[point][tri]][4*vertex1+1] << endl;
        }

        MEMORY ps_free( dS, vars_->triout->numberoftriangles );

        char foundweight;
        for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
        {
            if( weights[w+1] )
            {
                foundweight = 0;
                if( ! ((w/4)%2) )
                {
                    for( PS_SIT m=0; m<(PS_SIT)derx[point].size(); m+=2 )
                    {
                        if( derx[point][m] == (PS_SIT)weights[w] )
                        {
                            derx[point][m+1] += weights[w+1]/domainareas[point];
                            foundweight = 1;
                            break;
                        }
                    }
                    if( ! foundweight )
                    {
                        derx[point].push_back(weights[w]);
                        derx[point].push_back(weights[w+1]/domainareas[point]);
                    }
                }
                else
                {
                    for( PS_SIT m=0; m<(PS_SIT)dery[point].size(); m+=2 )
                    {
                        if( dery[point][m] == (PS_SIT)weights[w] )
                        {
                            dery[point][m+1] += weights[w+1]/domainareas[point];
                            foundweight = 1;
                            break;
                        }
                    }
                    if( ! foundweight )
                    {
                        dery[point].push_back(weights[w]);
                        dery[point].push_back(weights[w+1]/domainareas[point]);
                    }
                }
            }
        }
    }

    // for testing purposes against der1 only
    /*
      vars_->h1x = new MATRIX(vars_->lonc,vars_->lonc,vars_->lonc);
      vars_->h1y = new MATRIX(vars_->lonc,vars_->lonc,vars_->lonc);
      for(PS_SIT p=0; p<vars_->triout->numberofpoints; ++p)
      {
      for(PS_SIT j=0; j<derx[p].size(); j+=2)
      vars_->h1x->set(p,(PS_SIT)derx[p][j],derx[p][j+1]);
      for(PS_SIT j=0; j<dery[p].size(); j+=2)
      vars_->h1y->set(p,(PS_SIT)dery[p][j],dery[p][j+1]);
      }
      MATRIX *h1x2 = new MATRIX( vars_->lonc,vars_->lonc,vars_->lonc );
      MATRIX *h1y2 = new MATRIX( vars_->lonc,vars_->lonc,vars_->lonc );
      vars_->h1x->mult(vars_->h1x,h1x2,1,0,cdata_->numthreads);
      vars_->h1y->mult(vars_->h1y,h1y2,1,0,cdata_->numthreads);
      h1x2->plus( h1y2,vars_->c1,0,0,1,1,cdata_->numthreads );
      delete h1x2;
      delete h1y2;
      return;
    */

    // first derivatives, domain areas, and dS's computed
    // starting second derivative computation

    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc,vars_->lonc, initsize, 0, data_ );

    for( PS_SIT point=0; point<vars_->triout->numberofpoints; ++point )
    {
        // records surrounding points already analyzed
        vector<PS_SIT> gotitalready;
        PS_SIT vpoint;
        // looper is here because I haven't made sure the vertexes captured
        // are the most counter-clockwise lying
        // w.r.t the central point
        for( PS_SIT looper=0; looper<2; ++looper )
            for( PS_SIT tri=0; tri<(PS_SIT)surrtrilist[point].size(); ++tri )
            {
                PS_SIT tvertex=-1;
                PS_SIT tpoint;
                for( PS_SIT vertex=0; vertex<3; ++vertex )
                {
                    vpoint = vars_->triout->trianglelist[surrtrilist[point][tri]*3+vertex];
                    if( vpoint != point )
                    {
                        char itsgood = 1;
                        for( PS_SIT j=0; j<(PS_SIT)gotitalready.size(); ++j )
                        {
                            if( gotitalready[j] == vpoint )
                            {
                                itsgood = 0;
                                break;
                            }
                        }
                        if( itsgood )
                        {
                            tvertex = vertex;
                            tpoint = vpoint;
                        }
                    }
                    if( tvertex!=-1 )
                        break;
                }

                if( tvertex == -1 )
                    continue;

                gotitalready.push_back( tpoint );

                for( PS_SIT d=0; d<(PS_SIT)derx[point].size(); d+=2 )
                {
                    vars_->h1->set( point, (PS_SIT)derx[point][d],
                                    (
                                        derx[point][d+1] *
                                        dSsums[point][tri*6+tvertex*2  ] / domainareas[point] ) );
                }
                for( PS_SIT d=0; d<(PS_SIT)dery[point].size(); d+=2 )
                {
                    vars_->h1->set( point, (PS_SIT)dery[point][d],
                                    (
                                        dery[point][d+1] *
                                        dSsums[point][tri*6+tvertex*2+1] / domainareas[point] ) );
                }

                for( PS_SIT d=0; d<(PS_SIT)derx[tpoint].size(); d+=2 )
                {
                    vars_->h1->set( tpoint, (PS_SIT)derx[tpoint][d],
                                    (
                                        derx[tpoint][d+1] *
                                        dSsums[point][tri*6+tvertex*2  ] / domainareas[point] ) );
                }
                for( PS_SIT d=0; d<(PS_SIT)dery[tpoint].size(); d+=2 )
                {
                    vars_->h1->set( tpoint, (PS_SIT)dery[tpoint][d],
                                    (
                                        dery[tpoint][d+1] *
                                        dSsums[point][tri*6+tvertex*2+1] / domainareas[point] ) );
                }

            }

        free( dSsums[point] );
    }
    free( dSsums );

    /*
      for(PS_SIT l1=0; l1<vars_->lonc; ++l1){
      for(PS_SIT l2=0; l2<vars_->lonc; ++l2){
      if(vars_->h1->get(l1,l2))
      cout << l1 << " " << l2 << " " << vars_->h1->get(l1,l2) << " ";
      }
      cout << endl;
      }
    */

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
}
