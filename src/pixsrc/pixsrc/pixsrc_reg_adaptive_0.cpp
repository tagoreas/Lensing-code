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
#include "pixsrc_memory_templates.cpp"
#include <cmath>

void pixsrc_regularization::getder0(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    PS_SIT initsize = ( !data_->numgpu2use ) ? vars_->lonc : vars_->lonc;

    MEMORY ps_malloc( &(vars_->h1ptr), 1 );
    vars_->h1 = new (vars_->h1ptr)
        MATRIX( cdata_, data_, vars_->lonc, vars_->lonc, initsize, 0, data_ );

    for(PS_SIT s=0; s<vars_->lonc; ++s)
        vars_->h1->set( s, s, 1 );


    /*
      for(PS_unsignedSIT level=0; level<vars_->gridpointer->size(); level++)
      {
      //double weight = 1.0/pow(2,level+1);
      for(PS_unsignedSIT m=0; m<(*(vars_->gridpointer))[level].size(); m++)
      for(PS_SIT m2=0; m2<4; m2++)
      for(PS_unsignedSIT v=0; v<(*(vars_->gridpointer))[level][m][m2].size(); v++)
      if(vars_->c4c[(*(vars_->gridpointer))[level][m][m2][v]]!=-1)
      //vars_->h1->set(vars_->c4c[(*(vars_->gridpointer))[level][m][m2][v]],vars_->c4c[(*(vars_->gridpointer))[level][m][m2][v]],weight);
      vars_->h1->set(vars_->c4c[(*(vars_->gridpointer))[level][m][m2][v]],vars_->c4c[(*(vars_->gridpointer))[level][m][m2][v]],1);
      }
    */

    vars_->h1->mult( vars_->h1,vars_->c1, 1, 0, cdata_->numthreads );

}
