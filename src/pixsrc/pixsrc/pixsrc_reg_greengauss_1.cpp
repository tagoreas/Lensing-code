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
#include "pixsrc_geometry.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations_templates.cpp"
#include <cmath>
#include <algorithm>

// see
// "Development of Irregular-Grid Finite Difference Method (IFDM)
// for Governing Equations in Strong Form"
// by
// GEORGE XU @ Institute of High Performance Computing
// for a generalization of the method here

// the derivatives here are constructed on a grid where
// right points in an increasing direction
// and down points in an increasing direction

void pixsrc_regularization::greengaussgetder1
(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{

    {
        PS_SIT dim = 1;
        PS_SIT waitlist[1] = {7};
        OPERA pthreadswait(vars_,dim,waitlist);
    }

    MEMORY ps_malloc( &(vars_->h1xptr) , 1 );
    MEMORY ps_malloc( &(vars_->h1yptr) , 1 );

    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    vars_->h1x = new (vars_->h1xptr)
        MATRIX( cdata_, data_, vars_->lonc,vars_->lonc, initsize, 0, data_ );
    vars_->h1y = new (vars_->h1yptr)
        MATRIX( cdata_, data_, vars_->lonc,vars_->lonc, initsize, 0, data_ );

    // contains deltax,deltay,reference position(midpoint) for 3 vertices of each triangle
    // it connected the point at [tri][vertex] to the vertex preceding it.. -1 wraps around to 2
    double **dS;
    MEMORY ps_malloc( &(dS), vars_->triout->numberoftriangles, 12 );

    // lists the triangles that surround each point in triangulation
    vector< vector<PS_SIT> > surrtrilist(vars_->triout->numberofpoints,vector<PS_SIT>());

    // holds areas of gradient smoothing domain surrounding each point
    double domainareas[vars_->triout->numberofpoints];
    std::fill(domainareas,domainareas+vars_->triout->numberofpoints,0);

    // fill in dS for computations later
    {
        char isbadtri;
        PS_SIT pvertex, ppvertex, nextvertex;
        double centroidx, centroidy, vertices[6];
        for(PS_SIT tri=0; tri<vars_->triout->numberoftriangles; ++tri)
        {
            for( PS_SIT v=0; v<6; ++v )
                vertices[v] = vars_->triout->pointlist[
                    vars_->triout->trianglelist[tri*3+v/2] * 2 + v%2 ];

            // here we attempt to find bad triangles
            isbadtri = 0;
            pvertex  = 2;
            ppvertex = 1;
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
            if(isbadtri)
            {
                std::fill(dS[tri],dS[tri]+12,0);
                continue;
            }

            // find centroid of triangle
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
                // reference is midpoint between two vertices
                dS[tri][vertex*4+2] = (vertices[vertex*2  ] + vertices[nextvertex*2  ]) / 2.0;
                dS[tri][vertex*4+3] = (vertices[vertex*2+1] + vertices[nextvertex*2+1]) / 2.0;
                // centroid stored in this reference frame
                dS[tri][vertex*4  ] = centroidx - dS[tri][vertex*4+2];
                dS[tri][vertex*4+1] = centroidy - dS[tri][vertex*4+3];

                nextvertex = vertex;
            }
        }
    }

    // figure out to which triangles each point belongs
    for(PS_SIT t=0; t<vars_->triout->numberoftriangles*3; ++t)
    {
        surrtrilist[vars_->triout->trianglelist[t]].push_back(t/3);
    }

    // compute dS-contained derivatives
    {
        char isbadtri;
        PS_SIT vertex1, vertex2, vertex3, weightindex;
        double negatefacx, negatefacy, tempx, tempy, slope, centralshifted[4], coords[8];

        PS_SIT setpos_x, setpos_y, weights_capacity = 0;
        double *weights = 0, *weights_for_sort_x = 0, *weights_for_sort_y = 0;


        for(PS_SIT point=0; point<vars_->triout->numberofpoints; ++point)
        {
            // first index hold index of point. second holds value
            if( !weights )
            {
                weights_capacity = surrtrilist[point].size();

                MEMORY ps_malloc( &weights,            weights_capacity*16 );
                MEMORY ps_malloc( &weights_for_sort_x, weights_capacity*4  );
                MEMORY ps_malloc( &weights_for_sort_y, weights_capacity*4  );
            }

            // get weights and areas from all surrounding triangles
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
                    std::fill( &weights[tri*16], &weights[tri*16] + 16, 0 );
                    continue;
                }

                weightindex = 0;

                // vertex2 and vertex3 connect to the central point (vertex1)
                // vertex1 and vertex2 are to be referenced in dS
                vertex1 = -1;
                for(PS_SIT vertex=0; vertex<3; ++vertex)
                {
                    if(vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex] == point)
                    {
                        vertex1=vertex;
                        break;
                    }
                }

                vertex2 = (vertex1<2) ? vertex1+1 : 0;
                vertex3 = (vertex1>0) ? vertex1-1 : 2;

                // central point is shifted to coordinate systems where midpoints are origin
                centralshifted[0] = vars_->triout->pointlist[point*2  ] -
                    dS[surrtrilist[point][tri]][4*vertex1+2];
                centralshifted[1] = vars_->triout->pointlist[point*2+1] -
                    dS[surrtrilist[point][tri]][4*vertex1+3];
                centralshifted[2] = vars_->triout->pointlist[point*2  ] -
                    dS[surrtrilist[point][tri]][4*vertex2+2];
                centralshifted[3] = vars_->triout->pointlist[point*2+1] -
                    dS[surrtrilist[point][tri]][4*vertex2+3];

                // this block gets vertex1
                {
                    // now we find out which way the normal unit vector points
                    negatefacx = 1.0;
                    negatefacy = 1.0;
                    slope = -OPERA slope( 0.0, 0.0, dS[surrtrilist[point][tri]][4*vertex1],
                                          dS[surrtrilist[point][tri]][4*vertex1+1]         );

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
                    if( ! OPERA is_finite( slope ) )
                    {
                        negatefacx = ( centralshifted[0] < 0.0 ) ? 1 : -1;
                    }

                    tempx = negatefacx * 0.5 *
                        std::fabs(dS[surrtrilist[point][tri]][vertex1*4+1]);
                    tempy = negatefacy * 0.5 *
                        std::fabs(dS[surrtrilist[point][tri]][vertex1*4  ]);

                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                    weights[tri*16+weightindex++] = tempx;
                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex3];
                    weights[tri*16+weightindex++] = tempx;

                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                    weights[tri*16+weightindex++] = tempy;
                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex3];
                    weights[tri*16+weightindex++] = tempy;
                }

                // this block gets vertex2
                {
                    // now we find out which way the normal unit vector points
                    negatefacx = 1.0;
                    negatefacy = 1.0;
                    slope = -OPERA slope( 0.0, 0.0, dS[surrtrilist[point][tri]][4*vertex2],
                                          dS[surrtrilist[point][tri]][4*vertex2+1]         );

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
                    if( ! OPERA is_finite( slope ) )
                    {
                        negatefacx = ( centralshifted[2] < 0.0 ) ? 1 : -1;
                    }

                    tempx = negatefacx * 0.5 *
                        std::fabs(dS[surrtrilist[point][tri]][vertex2*4+1]);
                    tempy = negatefacy * 0.5 *
                        std::fabs(dS[surrtrilist[point][tri]][vertex2*4  ]);

                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex2];
                    weights[tri*16+weightindex++] = tempx;
                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                    weights[tri*16+weightindex++] = tempx;

                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex2];
                    weights[tri*16+weightindex++] = tempy;
                    weights[tri*16+weightindex++] =
                        vars_->triout->trianglelist[3*surrtrilist[point][tri]+vertex1];
                    weights[tri*16+weightindex++] = tempy;
                }

                // this is not ideal but ..
                // the 4 points that form the quadrilateral that overlaps with the
                // domain area for this triangle have coordinates here with their
                // origins at the midpoint connected to vertex1. So here I do some
                // coordinate transformations ..
                coords[0] = centralshifted[0];
                coords[1] = centralshifted[1];
                coords[2] = centralshifted[0]-centralshifted[2];
                coords[3] = centralshifted[1]-centralshifted[3];
                coords[4] = dS[surrtrilist[point][tri]][4*vertex1  ];
                coords[5] = dS[surrtrilist[point][tri]][4*vertex1+1];
                coords[6] = 0.0;
                coords[7] = 0.0;
                domainareas[point] += GEOM areapoly(coords,4);
            }

            // remove duplicates in weights
            for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
            {
                if( weights[w] == -1 || (w/4)%2 )
                    continue;

                for( PS_SIT w2=w+2; w2<(PS_SIT)surrtrilist[point].size()*16; w2+=2 )
                {
                    if( (w2/4)%2 )
                        continue;

                    if( weights[w] == weights[w2] )
                    {
                        weights[w2] = -1;
                        weights[w+1] += weights[w2+1];
                    }
                }
            }
            for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
            {
                if( weights[w] == -1 || !((w/4)%2) )
                    continue;

                for( PS_SIT w2=w+2; w2<(PS_SIT)surrtrilist[point].size()*16; w2+=2 )
                {
                    if( !((w2/4)%2) )
                        continue;

                    if( weights[w] == weights[w2] )
                    {
                        weights[w2] = -1;
                        weights[w+1] += weights[w2+1];
                    }
                }
            }

            // store all weights into an array for sorting
            setpos_x = setpos_y = 0;
            for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
            {
                if( weights[w+1] && weights[w]!=-1 )
                {
                    if( ! ((w/4)%2) )
                        weights_for_sort_x[setpos_x++] = weights[w];

                    else
                        weights_for_sort_y[setpos_y++] = weights[w];
                }
            }

            std::sort( weights_for_sort_x, weights_for_sort_x + setpos_x );
            std::sort( weights_for_sort_y, weights_for_sort_y + setpos_y );

            // actually enter weights into x-direction derivative operators
            for( PS_SIT g=0; g<setpos_x; ++g )
                for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
                {
                    if( ! ((w/4)%2) &&
                        weights_for_sort_x[g] == weights[w] )
                    {
                        vars_->h1x->set( point, (PS_SIT)weights[w],
                                         (weights[w+1]/domainareas[point]) );
                        break;
                    }
                }

            // actually enter weights into y-direction derivative operators
            for( PS_SIT g=0; g<setpos_y; ++g )
                for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
                    if( (w/4)%2  &&
                        weights_for_sort_y[g] == weights[w] )
                    {
                        vars_->h1y->set( point, (PS_SIT)weights[w],
                                         (weights[w+1]/domainareas[point]) );
                        break;
                    }
            /*
            // original (before row-ajor, then column sorting of matrix
            for( PS_SIT w=0; w<(PS_SIT)surrtrilist[point].size()*16; w+=2 )
            {
            if( weights[w+1] )
            {
            if( ! ((w/4)%2) )
            vars_->h1x->set( point, (PS_SIT)weights[w],
            weights[w+1]/domainareas[point] );
            else
            vars_->h1y->set( point, (PS_SIT)weights[w],
            weights[w+1]/domainareas[point] );
            }
            }
            */

            // if we just analyzed the last point or we need more memory for the
            // next point
            if( point == vars_->triout->numberofpoints-1 ||
                (PS_SIT)surrtrilist[point+1].size() > weights_capacity )
            {
                MEMORY ps_free( weights            );
                MEMORY ps_free( weights_for_sort_x );
                MEMORY ps_free( weights_for_sort_y );
                weights            = 0;
                weights_for_sort_x = 0;
                weights_for_sort_y = 0;
            }
        }
    }

    // complete line integral around domain area for edge triangles
    {
        char gotit1, gotit2;
        PS_SIT vertex1, vertex2, vertex3, point1, point2, point3, tri, tpoint;
        double shiftedx, shiftedy, tempx, tempy, negatefacx, negatefacy, slope;

        double *adjust_x, *adjust_y;
        // must be at least 4. i chose 10
        PS_SIT adjust_capacity_x = 10, adjust_capacity_y = 10;
        PS_SIT adjust_num_x = 0,       adjust_num_y = 0;
        MEMORY ps_malloc( &adjust_x, adjust_capacity_x*3 );
        MEMORY ps_malloc( &adjust_y, adjust_capacity_y*3 );

        for( PS_SIT e=0; e<vars_->triout->numberofedges; ++e )
        {
            // skip if not on edge
            if( vars_->triout->edgemarkerlist[e] != 1 )
                continue;

            point1 = vars_->triout->edgelist[e*2  ];
            point2 = vars_->triout->edgelist[e*2+1];

            // find vertex1, vertex2, and triangle these belong to
            vertex1 = vertex2 = tri = -1;
            for( PS_SIT t=0; t<(PS_SIT)surrtrilist[point1].size(); ++t )
            {
                gotit1 = gotit2 = 0;
                for( PS_SIT v=0; v<3; ++v )
                {
                    tpoint = vars_->triout->trianglelist[surrtrilist[point1][t]*3+v];

                    if( tpoint == point1 )
                    {
                        vertex1 = v;
                        gotit1  = 1;
                    }
                    else if( tpoint == point2 )
                    {
                        vertex2 = v;
                        gotit2  = 1;
                    }
                    if( gotit1 && gotit2 )
                    {
                        tri = surrtrilist[point1][t];
                        break;
                    }
                }
                if( tri != -1 )
                    break;
            }

            vertex3 = 3 - vertex2 - vertex1;
            point3  = vars_->triout->trianglelist[tri*3+vertex3];

            // I could use the midpoint reference coordinate stored in dS
            // here, but it looks like it'll actually be computationally
            // more work. Plus, it's just more annoying bookkeeping
            // I'll compute it internally to this block.

            // here I get the weight for point1
            // the point2 weight is easily relatable to this

            shiftedx = vars_->triout->pointlist[point2*2  ] -
                vars_->triout->pointlist[point1*2  ];
            shiftedy = vars_->triout->pointlist[point2*2+1] -
                vars_->triout->pointlist[point1*2+1];

            negatefacx = 1.0;
            negatefacy = 1.0;
            slope = -OPERA slope( 0.0, 0.0, shiftedx, shiftedy );

            // have to figure out what side of the triangle we're on
            // so I examine where the 3rd point lies w.r.t. point1 and point2
            if( ( vars_->triout->pointlist[point3*2+1] -
                  vars_->triout->pointlist[point1*2+1] ) <
                slope * ( vars_->triout->pointlist[point3*2] -
                          vars_->triout->pointlist[point1*2] ) )
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
            if( ! OPERA is_finite( slope ) )
            {
                negatefacx = ( vars_->triout->pointlist[point3*2] -
                               vars_->triout->pointlist[point1*2] < 0.0 ) ? 1 : -1;
            }

            tempx = negatefacx * std::fabs( shiftedy ) / 2.0;
            tempy = negatefacy * std::fabs( shiftedx ) / 2.0;

            // set weights in temporary adjustment arrays

            // check for overflow
            if( adjust_capacity_x < adjust_num_x + 4 )
            {
                OPERA resize( &adjust_x, adjust_capacity_x*3, adjust_capacity_x*2*3 );
                adjust_capacity_x *= 2;
            }
            if( adjust_capacity_y < adjust_num_y + 4 )
            {
                OPERA resize( &adjust_y, adjust_capacity_y*3, adjust_capacity_y*2*3 );
                adjust_capacity_y *= 2;
            }

            // set weights x
            adjust_x[adjust_num_x*3  ] = point1;
            adjust_x[adjust_num_x*3+1] = point1;
            adjust_x[adjust_num_x*3+2] = 0.75 * tempx/domainareas[point1];
            ++adjust_num_x;
            adjust_x[adjust_num_x*3  ] = point1;
            adjust_x[adjust_num_x*3+1] = point2;
            adjust_x[adjust_num_x*3+2] = 0.25 * tempx/domainareas[point1];
            ++adjust_num_x;
            adjust_x[adjust_num_x*3  ] = point2;
            adjust_x[adjust_num_x*3+1] = point2;
            adjust_x[adjust_num_x*3+2] = 0.75 * tempx/domainareas[point2];
            ++adjust_num_x;
            adjust_x[adjust_num_x*3  ] = point2;
            adjust_x[adjust_num_x*3+1] = point1;
            adjust_x[adjust_num_x*3+2] = 0.25 * tempx/domainareas[point2];
            ++adjust_num_x;

            // set weights y
            adjust_y[adjust_num_y*3  ] = point1;
            adjust_y[adjust_num_y*3+1] = point1;
            adjust_y[adjust_num_y*3+2] = 0.75 * tempy/domainareas[point1];
            ++adjust_num_y;
            adjust_y[adjust_num_y*3  ] = point1;
            adjust_y[adjust_num_y*3+1] = point2;
            adjust_y[adjust_num_y*3+2] = 0.25 * tempy/domainareas[point1];
            ++adjust_num_y;
            adjust_y[adjust_num_y*3  ] = point2;
            adjust_y[adjust_num_y*3+1] = point2;
            adjust_y[adjust_num_y*3+2] = 0.75 * tempy/domainareas[point2];
            ++adjust_num_y;
            adjust_y[adjust_num_y*3  ] = point2;
            adjust_y[adjust_num_y*3+1] = point1;
            adjust_y[adjust_num_y*3+2] = 0.25 * tempy/domainareas[point2];
            ++adjust_num_y;

            /*
            // original before row-major then column sorting
            // set weights in derivative operators
            // point1
            vars_->h1x->set( point1, point1, 0.75 * tempx/domainareas[point1] );
            vars_->h1x->set( point1, point2, 0.25 * tempx/domainareas[point1] );
            vars_->h1y->set( point1, point1, 0.75 * tempy/domainareas[point1] );
            vars_->h1y->set( point1, point2, 0.25 * tempy/domainareas[point1] );
            // point2
            vars_->h1x->set( point2, point2, 0.75 * tempx/domainareas[point2] );
            vars_->h1x->set( point2, point1, 0.25 * tempx/domainareas[point2] );
            vars_->h1y->set( point2, point2, 0.75 * tempy/domainareas[point2] );
            vars_->h1y->set( point2, point1, 0.25 * tempy/domainareas[point2] );
            */
        }

        //vars_->h1x->remove_duplicates();
        //vars_->h1y->remove_duplicates();

        vars_->h1x->atasbc( adjust_x, adjust_num_x );
        vars_->h1y->atasbc( adjust_y, adjust_num_y );

        MEMORY ps_free( adjust_x );
        MEMORY ps_free( adjust_y );
    }


    MATRIX *dummy1, *dummy2;
    MEMORY ps_malloc( &dummy1 , 1 );
    MEMORY ps_malloc( &dummy2 , 1 );

    initsize = ( !data_->numgpu2use ) ? vars_->lonc : 0;

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

    MEMORY ps_free( dS, vars_->triout->numberoftriangles );

}
