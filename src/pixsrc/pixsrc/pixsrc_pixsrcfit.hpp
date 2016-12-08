#ifndef PS_FPT
#ifdef SINGLE_PRECISION
#define PS_FPT float
#endif
#ifdef DOUBLE_PRECISION
#define PS_FPT double
#endif
#endif



#ifndef PIXSRC_PIXSRCFIT_HPP_
#define PIXSRC_PIXSRCFIT_HPP_

#define PIXSRCFIT pixsrc_pixsrcfit::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_pixsrcfit
{

public:

    static void pixsrcfit( inputdata*,commoninputdata*,lensvar* );

private:

    static double setsourceandreturnchi2 ( void*, double*       );

    pixsrc_pixsrcfit(const pixsrc_pixsrcfit&);
    pixsrc_pixsrcfit operator=(const pixsrc_pixsrcfit&);

};

#endif
