#ifndef PS_FPT
#ifdef SINGLE_PRECISION
#define PS_FPT float
#endif
#ifdef DOUBLE_PRECISION
#define PS_FPT double
#endif
#endif



#ifndef PIXSRC_COMMON_ANALYTIC_HPP_
#define PIXSRC_COMMON_ANALYTIC_HPP_

#define COMMONANALYTIC pixsrc_common_analytic::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"
#include "pixsrc_vector.hpp"

class pixsrc_common_analytic
{

public:

    static void   setusersource           ( inputdata*, commoninputdata*, lensvar*, double*, double*, double**, int, double*, VECTOR* );
    static void   expandsrc               ( inputdata*,commoninputdata*,lensvar*, double*, double**, double*  );
    static void   expandsrcwithflags      ( inputdata*,commoninputdata*,lensvar*, double*, double**, double*  );
    static void   srcflux                 ( inputdata*,commoninputdata*,lensvar*, double*, VECTOR*,
                                            int, double*, int*, int, double**                        );
    static double sersicflux              ( double, double, double, double, double, double           );
    static double vectorflux              ( double, double, double, struct triangulateio*, double*, int* );

private:

    pixsrc_common_analytic(const pixsrc_common_analytic&);
    pixsrc_common_analytic operator=(const pixsrc_common_analytic&);

};

#endif
