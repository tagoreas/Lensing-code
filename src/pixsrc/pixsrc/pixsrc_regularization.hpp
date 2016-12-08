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



#ifndef PIXSRC_REGULARIZATION_HPP_
#define PIXSRC_REGULARIZATION_HPP_

#define REGULAR pixsrc_regularization::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_regularization
{
public:

    static void getsersicregmin  ( inputdata*, commoninputdata*, lensvar* );
    static void getsersicregmax  ( inputdata*, commoninputdata*, lensvar* );
    static void greengaussgetder1( inputdata*, commoninputdata*, lensvar* );
    static void greengaussgetder2( inputdata*, commoninputdata*, lensvar* );
    static void getder0          ( inputdata*, commoninputdata*, lensvar* );
    static void getder12         ( inputdata*, commoninputdata*, lensvar* );
    static void getdercartesian0 ( inputdata*, commoninputdata*, lensvar* );
    static void getdercartesian1 ( inputdata*, commoninputdata*, lensvar* );
    static void getdercartesian2 ( inputdata*, commoninputdata*, lensvar* );
    static void regshapelets     ( inputdata*, commoninputdata*, lensvar* );

private:

    pixsrc_regularization(const pixsrc_regularization&);
    pixsrc_regularization operator=(const pixsrc_regularization&);
};

#endif
