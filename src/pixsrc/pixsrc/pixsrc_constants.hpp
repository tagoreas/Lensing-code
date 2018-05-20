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



#ifndef PIXSRC_CONSTANTS_HPP_
#define PIXSRC_CONSTANTS_HPP_

#define CONSTANT pixsrc_constants::

#include <cstring>

class pixsrc_constants
{

public:

    const static size_t uintmax;
    const static long   longmax;
    const static PS_SIT    regspeed;
    const static PS_SIT    longstringsize;
    const static PS_SIT    shortwaitms;
    const static PS_SIT    longwaitms;

    const static double d_inf;
    const static double d_epsilon;
    const static float  s_inf;
    const static float  s_epsilon;

    const static double smallnumber;
    const static double pi;
    const static double twopi;
    const static double piby2;
    const static double threepiby2;
    const static double sqrtpi;
    const static double smalltriangleanglecutoff;
    const static double deg2rad;
    const static double arc2rad;
    const static double rad2arc;
    const static char   dir_in[];
    const static char   dir_out[];
    const static char   bnseparator[];
    const static char   imagelistfile[];
    const static char   parameterfile[];
    const static char   shapeintmaskfilereg[];
    const static char   imagemaskfilereg[];
    const static char   imagemaskfileascii[];
    const static char   magmaskfilereg[];
    const static char   magmaskfileascii[];
    const static char   badpixelmaskfilereg[];
    const static char   badpixelmaskfileascii[];
    const static char   chi2maskfilereg[];
    const static char   chi2maskfileascii[];
    const static char   srcmaskfilereg[];
    const static char   srcmaskfileascii[];
    const static char   imagesfilestart[];
    const static char   imagesfileendreg[];
    const static char   imagesfileendascii[];
    const static char   fitsfile[];
    const static char   psffile[];
    const static double cutoff;
    const static double fwhm2sigma;
    const static double mm_stepsize;
    const static double mm_stepsize_area;
    const static PS_SIT    precision0;

    const static PS_SIT    shapelet_min_pixelsplit;
    const static PS_SIT    shapelet_max_pixelsplit;
    const static PS_SIT    shapelet_erf_approx_order;
    const static double shapelet_norm[251];

    const static double shapelet_max[51];
    const static double laguerre_norm[];
    const static double laguerre_min[];
    const static double laguerre_max[];

    // Jonathan Richard Shewchuk's Triangle code
    const static char triswitchesnominangle[];
    const static char triswitchesholefilling[];

};

#endif
