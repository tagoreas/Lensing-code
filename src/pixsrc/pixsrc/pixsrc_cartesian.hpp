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



#ifndef PIXSRC_CARTESIAN_HPP_
#define PIXSRC_CARTESIAN_HPP_

#define CARTESIAN pixsrc_cartesian::

#include "pixsrc_common.hpp"
#include <cmath>
#include <algorithm> // for max

class pixsrc_cartesian
{

public:

    pixsrc_cartesian(PS_SIT, inputdata*, commoninputdata*, PS_SIT );
    ~pixsrc_cartesian();

    void createc();
    void printsource();
    void getmagnification();

    lensvar *vars_;

private:

    inputdata *data_;
    commoninputdata *cdata_;

    void creategrid0();
    void makecg( double, PS_SIT, PS_SIT );
    void createc4c();
    void flagsrcpixels();
    void createlo();

    pixsrc_cartesian( const pixsrc_cartesian& );
    pixsrc_cartesian operator=( const pixsrc_cartesian& );

};

#endif
