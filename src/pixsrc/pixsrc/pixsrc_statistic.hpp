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



#ifndef PIXSRC_STATISTIC_HPP_
#define PIXSRC_STATISTIC_HPP_

#define STATISTIC pixsrc_statistic::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_statistic
{

public:

    static void* computees              ( void* );
    static void* computeed              ( void* );
    static void* computedeta            ( void* );
    static void* computedetc            ( void* );
    static void* evidencecalculator     ( void* );

    static void  actualcomputeed        ( inputdata*, commoninputdata*, lensvar*);
    static void  computeinterperrors    ( inputdata*, commoninputdata*, lensvar*);

private:

    static void   evidencecalculatorbody( inputdata*, commoninputdata*, lensvar*                    );
    static void   computedetabody       ( inputdata*, commoninputdata*, lensvar*                    );
    static void   computedetcbody       ( inputdata*, commoninputdata*, lensvar*                    );
    static void   computeesbody         ( inputdata*, commoninputdata*, lensvar*                    );
    static void   computeedbody         ( inputdata*, commoninputdata*, lensvar*                    );

    pixsrc_statistic( const pixsrc_statistic& );
    pixsrc_statistic operator=( const pixsrc_statistic& );

};

#endif
