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



#ifndef PIXSRC_ADAPTIVE_HPP_
#define PIXSRC_ADAPTIVE_HPP_

#define ADAPTIVE pixsrc_adaptive::

#include "pixsrc_common.hpp"

class pixsrc_adaptive
{

public:

    pixsrc_adaptive(PS_SIT imagenumber, inputdata datavec[], commoninputdata *cdata__, PS_SIT magtracker_);
    ~pixsrc_adaptive();

private:

    inputdata *data_;
    commoninputdata *cdata_;
    lensvar *vars_;

    void flagsrcpixels();
    void createc4c();

    pixsrc_adaptive(const pixsrc_adaptive &ps);
    pixsrc_adaptive operator=(const pixsrc_adaptive &ps);

};

#endif
