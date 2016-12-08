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



#include "pixsrc_wcs.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_memory_templates.cpp"
#include <cmath>                                        // for cos
#include <cstdlib>                                      // for atoi
#include <iomanip>

#include "pixsrc_constants.hpp"
#include "pixsrc_operations.hpp"

void pixsrc_wcs::getimgpixcoord( ps_WorldCoor *wcs, PS_SIT imgy, commoninputdata *cdata_,
                                 double ra, double dec, double *x, double *y )
{
    PS_SIT offscale=0;
    EXTERNAL ps_wcsc2pix(wcs,ra,dec,cdata_->coordsysfinal,x,y,&offscale);
    *x -= 1;
    *y  = imgy-*y;
}

void pixsrc_wcs::getimgwcscoord( ps_WorldCoor *wcs, PS_SIT imgy,
                                 double x, double y, double *ra, double *dec )
{
    x+=1;
    y=imgy-y;

    /*
      char wcstring[150];
      EXTERNAL ps_pix2wcst (data_->wcs, x, y, wcstring, (char)64);

      char **vec;
      PS_SIT numvec;
      OPERA split( wcstring, " ", &vec, &numvec );

      *ra =  OPERA todouble(vec[0]);
      *dec = OPERA todouble(vec[1]);

      MEMORY ps_free( vec, numvec );
    */

    EXTERNAL ps_pix2wcs (wcs, x, y, ra, dec );

    /*
      while(*ra<0)
      *ra+=360;
      while(*ra>360)
      *ra-=360;
      while(*dec<-90)
      *dec+=180;
      while(*dec>90)
      *dec-=180;
      */
}

void pixsrc_wcs::getwcslfromwcss( double ra, double dec, double pra, double pdec,
                                  double r1, double r2, double *result            )
{
    result[0] = (ra -pra ) * cos(CONSTANT deg2rad*(dec+pdec)/2.0) * 3600.0 + r1;
    result[1] = (dec-pdec)                                        * 3600.0 + r2;
}

void pixsrc_wcs::getwcssfromwcsl( double ra, double dec, double pra, double pdec,
                                  double r1, double r2, double *result            )
{
    result[1] = (dec-r2)                                              / 3600.0 + pdec;
    result[0] = (ra -r1) / cos(CONSTANT deg2rad*(result[1]+pdec)/2.0) / 3600.0 + pra;
}

void pixsrc_wcs::getlfroms( double sx, double sy, double px, double py,
                            double rx, double ry, double *result        )
{
    // s is a point in image coordinates
    // r in cartesian. p in image coordinates
    result[0] =  sx+rx-px;
    result[1] = -sy+ry+py;
}

void pixsrc_wcs::getsfroml( double Lx, double Ly, double px, double py,
                            double rx, double ry, double *result        )
{
    // L is a point in Chuck's coordinates
    result[0] =  Lx-rx+px;
    result[1] = -Ly+ry+py;
}

void pixsrc_wcs::setsrcwcs( inputdata *data_, lensvar *vars_,
                            double crpix1, double crpix2, double crval1, double crval2 )
{
    vars_->wcsinfo[0]=crpix1+1;
    vars_->wcsinfo[1]=vars_->srcy-crpix2;
    vars_->wcsinfo[2]=crval1;
    vars_->wcsinfo[3]=crval2;
    vars_->wcsinfo[4]=data_->wcsinfo[4]*vars_->reduction;
    vars_->wcsinfo[5]=data_->wcsinfo[5]*vars_->reduction;
    vars_->wcsinfo[6]=data_->wcsinfo[6]*vars_->reduction;
    vars_->wcsinfo[7]=data_->wcsinfo[7]*vars_->reduction;
    vars_->wcsinfo[8]=data_->wcsinfo[8];
}
