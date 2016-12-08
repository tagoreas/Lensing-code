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
#include "pixsrc_printer_templates.cpp"
#include <cmath>

void set_shapelet_laplacian (lensvar *vars_, PS_SIT s1, PS_SIT s2, PS_SIT numterms, PS_SIT *term, double **c1temp)
{
    PS_SIT ind1, ind2, index[5];
    double val1, val2, value[5];
    index[0] = (s1+0)*vars_->num_shapelets2+(s2+0);
    index[1] = (s1-2)*vars_->num_shapelets2+(s2+0);
    index[2] = (s1+2)*vars_->num_shapelets2+(s2+0);
    index[3] = (s1+0)*vars_->num_shapelets2+(s2-2);
    index[4] = (s1+0)*vars_->num_shapelets2+(s2+2);
    value[0] = -(2*s1+1)-(2*s2+1);
    value[1] = std::sqrt ( s1   *(s1-1));
    value[2] = std::sqrt ((s1+1)*(s1+2));
    value[3] = std::sqrt ( s2   *(s2-1));
    value[4] = std::sqrt ((s2+1)*(s2+2));
    for (PS_SIT n1=0; n1<numterms; ++n1)
    {
        ind1 = index[term[n1]];
        val1 = value[term[n1]];
        for (PS_SIT n2=0; n2<numterms; ++n2)
        {
            ind2 = index[term[n2]];
            val2 = value[term[n2]];
            c1temp[ind1][ind2] += val1*val2;
        }
    }
}

void pixsrc_regularization::regshapelets(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    if (data_->use_shapelets==1)
    {
        PS_SIT ind1, ind2;
        double **c1temp;
        MEMORY ps_malloc (&c1temp, vars_->lonc, vars_->lonc);
        for (PS_SIT s1=0; s1<vars_->lonc; ++s1)
            std::fill (c1temp[s1], c1temp[s1]+vars_->lonc, 0.0);

        if (data_->regorder==0)
        {
            // zeroth order derivative regularization
            for(PS_SIT s=0; s<vars_->lonc; ++s)
                c1temp[s][s] += 1.0;
        }
        else if (data_->regorder==1)
        {
            // 1st derivative regularization

            // 1st derivative in x direction (squared)
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; ++s2)
            {
                ind1 = 1*vars_->num_shapelets2+s2;
                c1temp[ind1][ind1] += 1.0;

                ind1 = (vars_->num_shapelets1-2)*vars_->num_shapelets2+s2;
                c1temp[ind1][ind1] += vars_->num_shapelets1-1;

                ind1 = (vars_->num_shapelets1-1)*vars_->num_shapelets2+s2;
                c1temp[ind1][ind1] += vars_->num_shapelets1;

                for (PS_SIT s1=1; s1<vars_->num_shapelets1-1; ++s1)
                {
                    ind1 = (s1-1)*vars_->num_shapelets2+s2;
                    ind2 = (s1+1)*vars_->num_shapelets2+s2;
                    c1temp[ind1][ind1] += s1;
                    c1temp[ind2][ind2] += s1+1;
                    c1temp[ind1][ind2] += -std::sqrt(s1*(s1+1));
                    c1temp[ind2][ind1] += -std::sqrt(s1*(s1+1));
                }
            }

            // 1st derivative in y direction (squared)
            for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
            {
                ind1 = s1*vars_->num_shapelets2+1;
                c1temp[ind1][ind1] += 1.0;

                ind1 = s1*vars_->num_shapelets2+(vars_->num_shapelets2-2);
                c1temp[ind1][ind1] += vars_->num_shapelets2-1;

                ind1 = s1*vars_->num_shapelets2+(vars_->num_shapelets2-1);
                c1temp[ind1][ind1] += vars_->num_shapelets2;

                for (PS_SIT s2=1; s2<vars_->num_shapelets2-1; ++s2)
                {
                    ind1 = s1*vars_->num_shapelets2+(s2-1);
                    ind2 = s1*vars_->num_shapelets2+(s2+1);
                    c1temp[ind1][ind1] += s2;
                    c1temp[ind2][ind2] += s2+1;
                    c1temp[ind1][ind2] += -std::sqrt(s2*(s2+1));
                    c1temp[ind2][ind1] += -std::sqrt(s2*(s2+1));
                }
            }
        }
        else if (data_->regorder==2)
        {
            // Laplacian (squared)

            // first 9 terms
            // see shapelets paper for equations
            // numterms is the number of terms inside the ()^2 .. maximum is 5
            // term stores the indices of the (max 5) terms that contribute.
            // limit1/2 are limits for the sums over shapelet indices
            // outer1/2 loop over the first 9 sums
            PS_SIT numterms, term[5], limit1[2], limit2[2];
            term[0] = 0;
            for (PS_SIT outer1=0; outer1<3; ++outer1)
                for (PS_SIT outer2=0; outer2<3; ++outer2)
                {
                    numterms = 1 + (outer1==2 ? 1 : outer1+1) + (outer2==2 ? 1 : outer2+1);
                    if (outer1==0)
                    {
                        limit1[0] = 0;
                        limit1[1] = 2;
                        term[1] = 2;
                    }
                    else if (outer1==1)
                    {
                        limit1[0] = 2;
                        limit1[1] = vars_->num_shapelets1-2;
                        term[1] = 1;
                        term[2] = 2;
                    }
                    else
                    {
                        limit1[0] = vars_->num_shapelets1-2;
                        limit1[1] = vars_->num_shapelets1;
                        term[1] = 1;
                    }
                    if (outer2==0)
                    {
                        limit2[0] = 0;
                        limit2[1] = 2;
                        term[numterms-1] = 4;
                    }
                    else if (outer2==1)
                    {
                        limit2[0] = 2;
                        limit2[1] = vars_->num_shapelets2-2;
                        term[numterms-1] = 3;
                        term[numterms-2] = 4;
                    }
                    else
                    {
                        limit2[0] = vars_->num_shapelets2-2;
                        limit2[1] = vars_->num_shapelets2;
                        term[numterms-1] = 3;
                    }
                    for (PS_SIT s1=limit1[0]; s1<limit1[1]; ++s1)
                        for (PS_SIT s2=limit2[0]; s2<limit2[1]; ++s2)
                            set_shapelet_laplacian (vars_, s1, s2, numterms, term, c1temp);
                }

            // last 2 terms
            for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
                for (PS_SIT s2=vars_->num_shapelets2; s2<vars_->num_shapelets2+1; ++s2)
                {
                    ind1 = s1*vars_->num_shapelets2+s2-2;
                    c1temp[ind1][ind1] += s2*(s2-1);
                }
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; ++s2)
                for (PS_SIT s1=vars_->num_shapelets1; s1<vars_->num_shapelets1+1; ++s1)
                {
                    ind1 = (s1-2)*vars_->num_shapelets2+s2;
                    c1temp[ind1][ind1] += s1*(s1-1);
                }
        }
        else if (data_->regorder==-1)
        {
            // regularize by size of source (prefers smaller sources)

            // integral over x coordinate
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; ++s2)
            {
                ind1 = 1*vars_->num_shapelets2+s2;
                c1temp[ind1][ind1] += 1.0;

                ind1 = (vars_->num_shapelets1-2)*vars_->num_shapelets2+s2;
                c1temp[ind1][ind1] += vars_->num_shapelets1-1;

                ind1 = (vars_->num_shapelets1-1)*vars_->num_shapelets2+s2;
                c1temp[ind1][ind1] += vars_->num_shapelets1;

                for (PS_SIT s1=1; s1<vars_->num_shapelets1-1; ++s1)
                {
                    ind1 = (s1-1)*vars_->num_shapelets2+s2;
                    ind2 = (s1+1)*vars_->num_shapelets2+s2;
                    c1temp[ind1][ind1] += s1;
                    c1temp[ind2][ind2] += s1+1;
                    c1temp[ind1][ind2] += std::sqrt(s1*(s1+1));
                    c1temp[ind2][ind1] += std::sqrt(s1*(s1+1));
                }
            }

            // integral over y coordinate
            for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
            {
                ind1 = s1*vars_->num_shapelets2+1;
                c1temp[ind1][ind1] += 1.0;

                ind1 = s1*vars_->num_shapelets2+(vars_->num_shapelets2-2);
                c1temp[ind1][ind1] += vars_->num_shapelets2-1;

                ind1 = s1*vars_->num_shapelets2+(vars_->num_shapelets2-1);
                c1temp[ind1][ind1] += vars_->num_shapelets2;

                for (PS_SIT s2=1; s2<vars_->num_shapelets2-1; ++s2)
                {
                    ind1 = s1*vars_->num_shapelets2+(s2-1);
                    ind2 = s1*vars_->num_shapelets2+(s2+1);
                    c1temp[ind1][ind1] += s2;
                    c1temp[ind2][ind2] += s2+1;
                    c1temp[ind1][ind2] += std::sqrt(s2*(s2+1));
                    c1temp[ind2][ind1] += std::sqrt(s2*(s2+1));
                }
            }
        }

        for (PS_SIT s1=0; s1<vars_->lonc; ++s1)
            for (PS_SIT s2=0; s2<vars_->lonc; ++s2)
            {
                vars_->c1->set (s1, s2, c1temp[s1][s2]);
            }

        MEMORY ps_free (c1temp, vars_->lonc);

    }
    else if (data_->use_shapelets==2)
    {
        PRINTER printerror(data_->print2screenname,
                           "polar shapelets disabled",
                           cdata_->print2screenmutex);
        /*
          PS_SIT index[2] = {0,0};
          double min_val = 0.01;
          for(PS_SIT s=0; s<vars_->lonc; ++s)
          {
          vars_->c1->set (s, s, index[0]*index[0]+index[1]*index[1]+min_val);

          // go to next radial mode if ready
          if (index[0]==-index[1])
          {
          ++index[0];
          index[1] = index[0]%2;
          }
          // otherwise increment angular mode if zero
          else if (index[1]==0)
          {
          index[1]+=2;
          }
          // otherwise increment angular mode if negative
          else if (index[1]<0)
          {
          index[1] = -index[1]+2;
          }
          // otherwise get imaginary angular mode if positive
          else
          {
          index[1] *= -1;
          }
          }
        */
    }
}
