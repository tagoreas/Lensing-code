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



#ifndef PIXSRC_ANALYTIC_SOURCE_HPP_
#define PIXSRC_ANALYTIC_SOURCE_HPP_

#define ANALYTICSRC pixsrc_analytic_source::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"
#include "pixsrc_vector.hpp"

class pixsrc_analytic_source
{

public:

    static void   getsource               (inputdata*,commoninputdata*,lensvar*);

private:

    static double lenssourceandreturnchi2 ( const gsl_vector*, void*);
    static double sersicflux              ( double, double, double, double, double, double           );
    static double vectorflux              ( double, double, double, struct triangulateio*, double*, PS_SIT* );
    static void   setusersource           ( inputdata*, commoninputdata*, lensvar*, double*, double*, double**, PS_SIT, double*, VECTOR* );
    static void   expandsrc               ( inputdata*,commoninputdata*,lensvar*, double*, double**, double*  );
    static void   expandsrcwithflags      ( inputdata*,commoninputdata*,lensvar*, double*, double**, double*  );
    static void   srcflux                 ( inputdata*,commoninputdata*,lensvar*, double*, VECTOR*,
                                            PS_SIT, double*);
    static void   imgflux                 ( inputdata*,commoninputdata*,lensvar*, double*,
                                            PS_SIT, double*, PS_SIT*, PS_SIT, double**                        );

    pixsrc_analytic_source(const pixsrc_analytic_source&);
    pixsrc_analytic_source operator=(const pixsrc_analytic_source&);

};

#endif
