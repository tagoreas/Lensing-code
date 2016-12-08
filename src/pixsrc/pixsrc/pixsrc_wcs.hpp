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



#ifndef PIXSRC_WCS_HPP_
#define PIXSRC_WCS_HPP_

#define HEADER pixsrc_wcs::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_wcs
{

public:

    static void getimgpixcoord ( ps_WorldCoor*, PS_SIT, commoninputdata*, double, double, double*, double* );
    static void getimgwcscoord ( ps_WorldCoor*, PS_SIT, double, double, double*, double* );
    static void setsrcwcs      ( inputdata*, lensvar*,         double, double, double,  double  );
    static void getwcslfromwcss(       double, double, double, double, double, double,  double* );
    static void getwcssfromwcsl(       double, double, double, double, double, double,  double* );
    static void getlfroms      (       double, double, double, double, double, double,  double* );
    static void getsfroml      (       double, double, double, double, double, double,  double* );

};

#endif
