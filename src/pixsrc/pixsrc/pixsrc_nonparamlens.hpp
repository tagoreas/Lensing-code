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



#ifndef PIXSRC_NONPARAMLENS_HPP_
#define PIXSRC_NONPARAMLENS_HPP_

#define NONPARAMLENS pixsrc_nonparamlens::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_nonparamlens
{

public:

    static double npl_penalty       (const gsl_vector*, void*);
    static double gsl_src_penalty   (inputdata*, commoninputdata*, lensvar*);
    static double gsl_lens_penalty  (inputdata*, commoninputdata*);
    static void   minimize          (inputdata*, commoninputdata*);
    static void   print_model       (inputdata*, commoninputdata*);
    static void   init_npl          (inputdata*, commoninputdata*);
    static void   get_tps_weights   (inputdata*, commoninputdata*, double);
    static void   get_alpha         (inputdata*, commoninputdata*,
                                     double, double, double*, double*);
    static double get_pot           (inputdata*, commoninputdata*,
                                     double, double);
    static double get_kappa         (inputdata*, commoninputdata*,
                                     double, double);
    static double get_grad_kappa    (inputdata*, commoninputdata*,
                                     double, double);

};

#endif
