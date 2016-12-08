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



#include "pixsrc_regularization.hpp"
#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_d2matrix.cpp"
#include "pixsrc_printer.hpp"
#include <cmath>

// this could be a lot more memory efficient (remove surrpointlist)

void pixsrc_regularization::getsersicregmax(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );

    D2MATRIX <PS_SIT> *dummy;
    MEMORY ps_malloc( &dummy, 1 );
    D2MATRIX <PS_SIT> *surrpoints = new  (dummy) D2MATRIX<PS_SIT> ( vars_->lonc );

    PS_SIT thisc, thisc1, thisc2;
    char addit;
    for( PS_SIT t=0; t<vars_->triout->numberoftriangles; ++t )
    {
        for( PS_SIT v=0; v<3; ++v )
        {
            thisc = vars_->triout->trianglelist[t*3+v];
            thisc1 = v+1;
            thisc2 = v+2;
            if(thisc1>2)
                thisc1 -= 3;
            if(thisc2>2)
                thisc2 -= 3;
            thisc1 = vars_->triout->trianglelist[t*3+thisc1];
            thisc2 = vars_->triout->trianglelist[t*3+thisc2];

            addit = 1;
            for( PS_SIT cc=0; cc<surrpoints->size(thisc); ++cc )
                if( surrpoints->get(thisc,cc)==thisc1 )
                {
                    addit = 0;
                    break;
                }

            if( addit )
                surrpoints->pushback( thisc, thisc1 );

            addit = 1;
            for( PS_SIT cc=0; cc<surrpoints->size(thisc); ++cc )
                if( surrpoints->get(thisc,cc)==thisc2 )
                {
                    addit = 0;
                    break;
                }

            if( addit )
                surrpoints->pushback( thisc, thisc2 );

        }
    }

    for( PS_SIT c=0; c<vars_->lonc; ++c )
        surrpoints->pushback( c, c );

    PS_SIT dummyint;

    for( PS_SIT c=0; c<vars_->lonc; ++c )
    {
        surrpoints->sort( c );

        for( PS_SIT cc=0; cc<surrpoints->size(c); ++cc )
        {
            dummyint = surrpoints->get( c, cc );

            if( c == dummyint )
                vars_->h1->set( c, c, (-(surrpoints->size(c)-1)/vars_->mps->get(c) ));
            else
                vars_->h1->set( c, dummyint, (1.0/vars_->mps->get(dummyint) ));
        }
    }

    surrpoints->~D2MATRIX();
    MEMORY ps_free( dummy );

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
}
