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



#ifndef PIXSRC_COMMON_ADAPTIVE_HPP_
#define PIXSRC_COMMON_ADAPTIVE_HPP_

#define COMMONADAPTIVE pixsrc_common_adaptive::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_common_adaptive
{
public:
    static void* createc          ( void* );
    static void* printsource      ( void* );
    static void* getmagnification ( void* );
    static void  createlo         ( inputdata*, commoninputdata*, lensvar*);

private:

    static void createlo_s2i           (inputdata*, commoninputdata*, lensvar*);
    static void createlo_i2s           (inputdata*, commoninputdata*, lensvar*);
    static void getmagnificationbody   (inputdata*, commoninputdata*, lensvar*);
    static void createcbody            (inputdata*, commoninputdata*, lensvar*);
    static void getmagfromtriangulation(inputdata*, commoninputdata*, lensvar*);
    static void actualprintsource      (inputdata*, commoninputdata*, lensvar*, VECTOR*);
    static void printsourcebody        (inputdata*, commoninputdata*, lensvar*);

    pixsrc_common_adaptive(const pixsrc_common_adaptive&);
    pixsrc_common_adaptive operator=(const pixsrc_common_adaptive&);
};

#endif
