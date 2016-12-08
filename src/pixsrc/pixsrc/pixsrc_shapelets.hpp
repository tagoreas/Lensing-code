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



#ifndef PIXSRC_SHAPELETS_HPP_
#define PIXSRC_SHAPELETS_HPP_

#define SHAPELETS pixsrc_shapelets::

#include "pixsrc_common.hpp"

class pixsrc_shapelets
{

public:

    pixsrc_shapelets(PS_SIT imagenumber, inputdata datavec[], commoninputdata *cdata__, PS_SIT magtracker_);
    ~pixsrc_shapelets();

private:

    inputdata *data_;
    commoninputdata *cdata_;
    lensvar *vars_;

    void flagsrcpixels();
    void createc4c();
    pixsrc_shapelets(const pixsrc_shapelets &ps);
    pixsrc_shapelets operator=(const pixsrc_shapelets &ps);

};

#endif
