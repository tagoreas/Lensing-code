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
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer.hpp"
#include <algorithm>

void pixsrc_regularization::getdercartesian2( inputdata *data_,
                                              commoninputdata *cdata_, lensvar *vars_ )
{
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );

    PS_SIT setorder[5];

    for(PS_SIT cc = 0; cc < vars_->lonc; cc++)
    {
        PS_SIT c0=vars_->c4cback[cc];
        double x = vars_->srcloc[c0*2  ];
        double y = vars_->srcloc[c0*2+1];

        setorder[0] = cc;

        setorder[1] = (x>vars_->srcloc[0] && vars_->c4c[c0-vars_->srcy]!=-1) ?
            vars_->c4c[c0-vars_->srcy] : -1;

        setorder[2] = (x<vars_->srcloc[(vars_->numsrcpoints-1)*2] &&
                       vars_->c4c[c0+vars_->srcy]!=-1) ?
            vars_->c4c[c0+vars_->srcy] : -1;

        setorder[3] = (y>vars_->srcloc[1] && vars_->c4c[c0-1]!=-1) ?
            vars_->c4c[c0-1] : -1;

        setorder[4] = (y<vars_->srcloc[(vars_->numsrcpoints-1)*2+1] &&
                       vars_->c4c[c0+1]!=-1) ?
            vars_->c4c[c0+1] : -1;

        std::sort( setorder, setorder + 5 );

        for( PS_SIT g=0; g<5; ++g )
        {
            if( setorder[g] == -1 )
                continue;

            if(      setorder[g] == vars_->c4c[c0-vars_->srcy]    )
            {
                vars_->h1->set( cc, vars_->c4c[c0-vars_->srcy], 1 );
            }
            else if( setorder[g] == vars_->c4c[c0+vars_->srcy]    )
            {
                vars_->h1->set( cc, vars_->c4c[c0+vars_->srcy], 1 );
            }
            else if( setorder[g] == vars_->c4c[c0-1]    )
            {
                vars_->h1->set( cc, vars_->c4c[c0-1], 1 );
            }
            else if( setorder[g] == vars_->c4c[c0+1]    )
            {
                vars_->h1->set( cc, vars_->c4c[c0+1], 1 );
            }
            else
            {
                vars_->h1->set( cc, cc, (-4) );
            }
        }

        /*
          if(x>vars_->srcloc[0] && vars_->c4c[c0-vars_->srcy]!=-1)
          vars_->h1->set(cc,vars_->c4c[c0-vars_->srcy],1);
          if(x<vars_->srcloc[(vars_->numsrcpoints-1)*2] && vars_->c4c[c0+vars_->srcy]!=-1)
          vars_->h1->set(cc,vars_->c4c[c0+vars_->srcy],1);
          if(y>vars_->srcloc[1] && vars_->c4c[c0-1]!=-1)
          vars_->h1->set(cc,vars_->c4c[c0-1],1);
          if(y<vars_->srcloc[(vars_->numsrcpoints-1)*2+1] && vars_->c4c[c0+1]!=-1)
          vars_->h1->set(cc,vars_->c4c[c0+1],1);
          vars_->h1->set(cc,cc,-4);
        */
    }

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
}
void pixsrc_regularization::getdercartesian1( inputdata *data_,
                                              commoninputdata *cdata_, lensvar *vars_ )
{
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );

    PS_SIT setorder[3];

    for(PS_SIT cc = 0; cc < vars_->lonc; cc++)
    {
        PS_SIT c0=vars_->c4cback[cc];
        double x = vars_->srcloc[c0*2  ];
        double y = vars_->srcloc[c0*2+1];

        setorder[0] = cc;

        setorder[1] = (x<vars_->srcloc[(vars_->numsrcpoints-1)*2  ] &&
                       vars_->c4c[c0+vars_->srcy]!=-1) ? vars_->c4c[c0+vars_->srcy] : -1;

        setorder[2] = (y<vars_->srcloc[(vars_->numsrcpoints-1)*2+1] &&
                       vars_->c4c[c0+1]!=-1) ? vars_->c4c[c0+1] : -1;

        std::sort( setorder, setorder + 3 );

        for( PS_SIT g=0; g<3; ++g )
        {
            if( setorder[g] == -1 )
                continue;

            if(      setorder[g] == vars_->c4c[c0+vars_->srcy]    )
            {
                vars_->h1->set( cc, vars_->c4c[c0+vars_->srcy], 1 );
            }
            else if( setorder[g] == vars_->c4c[c0+1]    )
            {
                vars_->h1->set( cc, vars_->c4c[c0+1], 1 );
            }
            else
            {
                vars_->h1->set( cc, cc, (-2) );
            }
        }

        /*
          if(x<vars_->srcloc[(vars_->numsrcpoints-1)*2  ] && vars_->c4c[c0+vars_->srcy]!=-1)
          vars_->h1->set(cc,vars_->c4c[c0+vars_->srcy],1);
          if(y<vars_->srcloc[(vars_->numsrcpoints-1)*2+1] && vars_->c4c[c0+1]!=-1)
          vars_->h1->set(cc,vars_->c4c[c0+1],1);
          vars_->h1->set(cc,cc,-2);
        */
    }

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
}
void pixsrc_regularization::getdercartesian0( inputdata *data_,
                                              commoninputdata *cdata_, lensvar *vars_ )
{
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );

    for(PS_SIT cc = 0; cc < vars_->lonc; cc++)
        vars_->h1->set( cc, cc, 1 );

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );
}
