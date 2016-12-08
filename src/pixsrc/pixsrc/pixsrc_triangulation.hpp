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



#ifndef PIXSRC_TRIANGULATION_HPP_
#define PIXSRC_TRIANGULATION_HPP_

#define TRIANGULATION pixsrc_triangulation::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_triangulation
{

public:

    static void* pstriangulate          (void*);

private:

    static void  pstriangulatebody     (inputdata*,commoninputdata*,lensvar*);
    static void  trimtrilist            (inputdata*,commoninputdata*,lensvar*,char*);
    static void  trimpointlist          (inputdata*,commoninputdata*,lensvar*,char*);
    static void  trimedgelist          (inputdata*,commoninputdata*,lensvar*,char*);
    static void  retriangulateholes     (inputdata*,commoninputdata*,lensvar*,char**,char**,char**,char**);

    pixsrc_triangulation(const pixsrc_triangulation&);
    pixsrc_triangulation operator=(const pixsrc_triangulation&);

};

#endif
