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



#ifndef PIXSRC_IRRCART_HPP_
#define PIXSRC_IRRCART_HPP_

#define IRRCART pixsrc_irrcart::

#include "pixsrc_common.hpp"

class pixsrc_irrcart
{

public:

    pixsrc_irrcart(PS_SIT imagenumber, inputdata datavec[], commoninputdata*, PS_SIT magtracker_);
    ~pixsrc_irrcart();

private:

    inputdata *data_;
    commoninputdata *cdata_;
    lensvar *vars_;

    void flagsrcpixels();
    void creategrid0();
    void createc4c();
    void quadrantsearch(double xx, double yy, double result[2][4]);
    void addsubgrid(vector<PS_SIT> &path);
    void getpathfrompositioninarrays(PS_SIT *pos, vector<PS_SIT> &result);
    void getpositioninarraysfrompath(vector<PS_SIT> &path, PS_SIT *result);
    PS_SIT  indexofnext(PS_SIT subgridLevel, vector<PS_SIT> &path, PS_SIT stopAt);
    bool nodeexists(vector<PS_SIT> &path);
    bool nodepossiblyexistshereorbelow(vector<PS_SIT> &path);
    void optimizegrid();
    void setsrclocall();
    void makegrid();
    void testforoverlappinggridpoints();
    void findsurroundingboxes(vector<PS_SIT> &path, PS_SIT lastIndex, vector< vector<PS_SIT> > &result);
    void testthisnodewithgivendensityandaddtogrid(double x, double y, double density);
    void positionofnode(vector<PS_SIT> &path, double *result);
    pixsrc_irrcart(const pixsrc_irrcart &ps);
    pixsrc_irrcart operator=(const pixsrc_irrcart &ps);

};

#endif
