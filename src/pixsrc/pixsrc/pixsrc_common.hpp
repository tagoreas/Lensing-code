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



#ifndef PIXSRC_COMMON_HPP_
#define PIXSRC_COMMON_HPP_

#define COMMON pixsrc_common::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_geometry.hpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_constants.hpp"
#include "pixsrc_matrix.hpp"

class pixsrc_common
{

public:

    static void* computesizepenalty     ( void*                                        );
    static void* computepenalty         ( void*                                        );
    static void* printimageplane        ( void*                                        );

    static void  resubsample            ( inputdata*, commoninputdata*, lensvar*       );
    static void  computemagpenalty      ( inputdata*, commoninputdata*, lensvar*       );
    static void  startlensing           ( inputdata*, commoninputdata*, lensvar*       );
    static void  testformismatchedimages( inputdata*, commoninputdata*, lensvar*       );
    static void  sourcereconstructions  ( inputdata*, commoninputdata*, lensvar*       );
    static void  findlambda             ( inputdata*, commoninputdata*, lensvar*       );
    static void  blurit                 ( inputdata*, commoninputdata*, lensvar*       );
    static void  creater4r              ( inputdata*, commoninputdata*, lensvar*       );
    static void  findsisterimages       ( inputdata*, commoninputdata*, lensvar*       );
    static void  printfinalmessages     ( inputdata*, commoninputdata*, lensvar*       );
    static void  lensgalaxy             ( inputdata*, commoninputdata*, lensvar*,
                                          MATRIX*, VECTOR*, VECTOR*                    );
    static void  penaltypenalizer       ( inputdata*, commoninputdata*, lensvar*,
                                          PS_SIT, PS_SIT, double, string                     );
    static void  raytrace               ( inputdata*, commoninputdata*,
                                          double , double , double*, double*,
                                          double*, double*, double*, double*, double*, PS_SIT, PS_SIT);
    static void  raytrace_s2i           ( inputdata*, commoninputdata*,
                                          double , double , PS_SIT*, double***            );

private:

    static void  subsample_ie           ( void*                                   );
    static void  resubsample_ie         ( void*                                   );

    static void  computesizepenaltybody ( inputdata*,commoninputdata*,lensvar*    );
    static void  computepenaltybody     ( inputdata*,commoninputdata*,lensvar*    );
    static void  printimageplanebody    ( inputdata*,commoninputdata*,lensvar*    );
    static void  printuvplane           ( inputdata*,commoninputdata*,lensvar*,
                                          VECTOR*, VECTOR*, VECTOR*);
    static void  printuvresiduals       ( inputdata*,commoninputdata*,lensvar*);

    pixsrc_common(const pixsrc_common&);
    pixsrc_common operator=(const pixsrc_common&);

};

#endif
