#ifndef PS_FPT
#ifdef SINGLE_PRECISION
#define PS_FPT float
#endif
#ifdef DOUBLE_PRECISION
#define PS_FPT double
#endif
#endif



#include "pixsrc_common_analytic.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_geometry.hpp"
#include "pixsrc_operations.hpp"
#include <cmath>

void pixsrc_common_analytic::setusersource( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *source, double *penalty, double **bounds, int numbounds, double *src0, VECTOR *fluxhere )
{
    double *sourceloc = 0;
    COMMONANALYTIC expandsrcwithflags( data_, cdata_, vars_, source, &sourceloc, src0 );

    // test bounds
    if( bounds && penalty )
    {
        double sum;
        double basepenalty = 1.0e100;
        int totmt = sourceloc[0]*8;

        for( int n=0; n<numbounds; ++n )
        {
            sum = 0;
            for( int j=0; j<totmt; ++j )
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

    COMMONANALYTIC srcflux( data_, cdata_, vars_, sourceloc, fluxhere, 1,
                            vars_->triout->pointlist, NULL, 1, NULL);

    MEMORY ps_free( sourceloc );
}

void pixsrc_common_analytic::expandsrc( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *src, double **srcreturn, double *src0 )
{
    int length = src0[0]*9+1;
    /*
      if( *srcreturn )
      MEMORY ps_free( *srcreturn );
    */
    if( !*srcreturn )
        MEMORY ps_malloc( &(*srcreturn), length );

    // expanding the abridged source parameters that amoeba uses for computeed
    int trackerloc = 0;
    (*srcreturn)[0] = src0[0];
    for( int j=0; j<src0[0]; ++j )
    {
        // copy user initial params
        std::copy( src0 + 1+j*17, src0 + 1+j*17+9,
                   *srcreturn + 1+j*9 );

        // copy newly varied params
        for( int g=0; g<8; ++g )
        {
            if( src0[10+j*17+g] )
            {
                (*srcreturn)[1+j*9+1+g] = src[trackerloc++];
            }
        }
    }
}

void pixsrc_common_analytic::expandsrcwithflags( inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *src, double **srcreturn, double *src0 )
{
    int length = src0[0]*17+1;
    /*
      if( *srcreturn )
      MEMORY ps_free( *srcreturn );
    */
    if( !*srcreturn )
        MEMORY ps_malloc( &(*srcreturn), length );

    // expanding the abridged source parameters that amoeba uses
    std::copy( src0, src0+length, *srcreturn );

    // putting in newly varied parameters
    int index = 0;
    for( int s=0; s<src0[0]; ++s )
        for( int p=s*17+1+9; p<s*17+1+17; ++p )
            if( src0[p] )
                (*srcreturn)[p-8] = src[index++];
}

// multithreading code .. decided not worth it
/*
  struct srcfluxstruct
  {
  inputdata *data_;
  int srcindex;
  double xpos, ypos, val;
  double a, b, c, d;
  };

  void* pixsrc_common_analytic::srcfluxthread( void* )
  {
  srcfluxstruct sfs = (srcfluxstruct*)void;

  sfs->val = COMMONANALYTIC sersicflux( sfs->a, xpos, ypos, sfs->b, sfs->c, sfs->d );
  }
*/

void pixsrc_common_analytic::srcflux( inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                      double *sourceloc_, VECTOR *source, int sourcetype,
                                      double *pos, int *posflags, int subsample, double **posss )
{
    // sourceloc is the source
    // source is where flux will be stored
    // sourcetype indicates whether sourceloc contains vary flags or not
    // pos holds positions where flux is to be calculated (can be NULL if subsampling)
    // posflags indicate whether to skip a certain position in source (NULL means use all)
    // subsample indicates whether to integrate flux over region or get flux at center
    // posss holds (source plane) positions of subsampled pixels (can be NULL if not subsampling)

    // GUG: write multithreaded code here. get all subsampled positions (or not subsamples) in one array and split it up across cpus/gpus.

    source->zeromeout ();

    double ss2 = subsample*subsample;
    double x=0,y=0,xx,yy,val,qaxis;
    double cos1,sin1,xrot,yrot;

    int numsrcs = sourceloc_[0];
    double *sourceloc;
    MEMORY ps_malloc( &sourceloc, 1 + numsrcs*9 + sourcetype*numsrcs*8  );
    std::copy(       sourceloc_,
                     sourceloc_ + 1 + numsrcs*9 + sourcetype*numsrcs*8,
                     sourceloc                                          );

    double *srcdata;

    // convert sources to pixel coordinates
    pthread_mutex_lock( cdata_->wcsmutex );
    for( int src=0; src<numsrcs; ++src )
    {
        // current source
        srcdata = &(sourceloc[1+src*9+sourcetype*src*8]);

        switch ((int)srcdata[0])
        {
            // if sersic
        case 0:
        {
            // converting center position from arcseconds to pixels
            double poss[2];
            HEADER getwcssfromwcsl( srcdata[2], srcdata[3], data_->pra, data_->pdec,
                                    data_->r1, data_->r2, poss                       );
            HEADER getimgpixcoord( data_, cdata_, poss[0], poss[1],
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
                int vecindex = srcdata[0]-1001;
                // convert shift from arcsec to pixel with a little hack

                // first find image position of first pixel (unshifted)
                double pos1[2] = {cdata_->vector_src_pos[vecindex][0],
                                  cdata_->vector_src_pos[vecindex][1]};
                double pos_orig[2], pos_shift[2], pos2[2];
                HEADER getimgpixcoord (data_, cdata_, pos1[0], pos1[1],
                                       &pos_orig[0], &pos_orig[1]);
                // now apply srcdata shift at position of first pixel and get its image position
                HEADER getwcssfromwcsl (srcdata[2], srcdata[3], pos1[0], pos1[1], 0, 0, pos2);
                HEADER getimgpixcoord (data_, cdata_, pos2[0], pos2[1],
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

    int postracker;

    // vector sync occurs here

    for( int src=0; src<numsrcs; ++src )
    {
        // this source
        srcdata = &(sourceloc[1+src*9+sourcetype*src*8]);

        switch ((int)srcdata[0])
        {
            // if sersic
        case 0:
        {
            // convert angle to one measure from +x axis
            qaxis = 1.0-srcdata[4];
            cos1 = std::cos(srcdata[5]*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2);
            sin1 = std::sin(srcdata[5]*CONSTANT deg2rad - data_->rollangle + CONSTANT piby2);

            postracker = 0;
            for(int s=0; s<source->get_size(); ++s)
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

                    val = COMMONANALYTIC sersicflux (srcdata[1], xrot, yrot,
                                                     qaxis, srcdata[6], srcdata[8]);
                    source->set( s, source->get(s) + val );
                }
                else
                {
                    val = 0;

                    // was gonna use this for multithreading
                    //if( cdata_->numthreads == 1 )
                    {
                        for( int ss=0; ss<ss2; ++ss )
                        {
                            xx = posss[postracker][ss*2  ] - srcdata[2];
                            yy = posss[postracker][ss*2+1] - srcdata[3];

                            xrot = xx*cos1 - yy*sin1;
                            yrot = xx*sin1 + yy*cos1;

                            val += COMMONANALYTIC sersicflux (srcdata[1], xrot, yrot,
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
        default:
        {
            // if vector image
            if (srcdata[0]>1000)
            {
                int vecindex = srcdata[0]-1001;
                int numpoints = cdata_->vector_src_numpix[vecindex];

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
                for (int m=0; m<4; ++m)
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
                for (int s=0; s<numpoints; ++s)
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

                postracker = 0;
                int triseed = 0;
                for(int s=0; s<source->get_size(); ++s)
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

                        val = COMMONANALYTIC vectorflux (srcdata[1], x, y, trioutloc,
                                                         cdata_->vector_src_flux[vecindex], &triseed);
                        source->set( s, source->get(s) + val );
                    }
                    else
                    {
                        val = 0;

                        for( int ss=0; ss<ss2; ++ss )
                        {
                            xx = posss[postracker][ss*2  ] - srcdata[2];
                            yy = posss[postracker][ss*2+1] - srcdata[3];

                            val += COMMONANALYTIC vectorflux (srcdata[1], x, y, trioutloc,
                                                              cdata_->vector_src_flux[vecindex], &triseed);
                        }

                        source->set( s, source->get(s) + val );
                    }

                    if( posflags )
                        ++postracker;
                }

                // restore original pixel positions to triangle struct and destroy transformed pixels
                trioutloc->pointlist = vertices_orig;
                MEMORY ps_free (vertices_transformed);
            }

            break;
        }
        }
    }

    if (subsample>1)
        source->set_scalar( 1.0 / ss2 );

    MEMORY ps_free( sourceloc );
}

double pixsrc_common_analytic::sersicflux( double I0, double xrot, double yrot,
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

double pixsrc_common_analytic::vectorflux (double I0, double x, double y,
                                           struct triangulateio *trioutloc, double *flux, int *triseed)
{
    // find where pixel x,y lies in the triangulation and return flux
    int index = GEOM search_triangle( trioutloc, triseed, x, y);

    if (index==-1)
        return 0;

    // find vertices of triangle
    double pos[2][3], vals[3], result;
    for (int f=0; f<3; ++f)
    {
        pos[0][f] = trioutloc->pointlist[trioutloc->trianglelist[index*3+f]*2  ];
        pos[1][f] = trioutloc->pointlist[trioutloc->trianglelist[index*3+f]*2+1];
        vals[f]   = flux[trioutloc->trianglelist[index*3+f]];
    }
    OPERA planarvalueinterpolation (pos, vals, x, y, &result);

    return result*I0;
}
