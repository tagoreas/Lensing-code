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



#include "pixsrc_constants.hpp"
#include <climits> // integral type limits
#include <cfloat>  // float type limits
#include <cmath>

const size_t pixsrc_constants::uintmax = (size_t)UINT_MAX;
const long   pixsrc_constants::longmax = LONG_MAX;

const float  pixsrc_constants::s_inf = FLT_MAX/100.0f;
const float  pixsrc_constants::s_epsilon = FLT_EPSILON*100.0f;
const double pixsrc_constants::d_inf = DBL_MAX/100.0;
const double pixsrc_constants::d_epsilon = DBL_EPSILON*100.0;

const double pixsrc_constants::smallnumber = d_epsilon*100.0;
const double pixsrc_constants::pi = 3.1415926535897932384626433832795028841971693993751;
const double pixsrc_constants::twopi = 2.0*pi;
const double pixsrc_constants::piby2 = pi/2.0;
const double pixsrc_constants::threepiby2 = 3.0*piby2;
const double pixsrc_constants::sqrtpi = std::sqrt (pi);
const double pixsrc_constants::smalltriangleanglecutoff = pi-smallnumber*100;
const double pixsrc_constants::deg2rad = pi/180.0;
const double pixsrc_constants::arc2rad = pi/(180.0*3600.0);
const double pixsrc_constants::rad2arc = (180.0*3600.0)/pi;
const char   pixsrc_constants::dir_in[] = "pixsrc_in/";
const char   pixsrc_constants::dir_out[] = "pixsrc_out/";
const char   pixsrc_constants::bnseparator[] = "_";
const char   pixsrc_constants::imagelistfile[] = ".include";
const char   pixsrc_constants::parameterfile[] = ".parameters";
const char   pixsrc_constants::shapeintmaskfilereg[] = ".shapeintmask.reg";
const char   pixsrc_constants::imagemaskfilereg[] = ".datamask.reg";
const char   pixsrc_constants::imagemaskfileascii[] = ".datamask.ascii";
const char   pixsrc_constants::magmaskfilereg[] = ".magmask.reg";
const char   pixsrc_constants::magmaskfileascii[] = ".magmask.ascii";
const char   pixsrc_constants::badpixelmaskfilereg[] = ".badpixelmask.reg";
const char   pixsrc_constants::badpixelmaskfileascii[] = ".badpixelmask.ascii";
const char   pixsrc_constants::chi2maskfilereg[] = ".chi2mask.reg";
const char   pixsrc_constants::chi2maskfileascii[] = ".chi2mask.ascii";
const char   pixsrc_constants::srcmaskfilereg[] = ".srcmask.reg";
const char   pixsrc_constants::srcmaskfileascii[] = ".srcmask.ascii";
const char   pixsrc_constants::imagesfilestart[] = ".mmimages.";
const char   pixsrc_constants::imagesfileendreg[] = ".reg";
const char   pixsrc_constants::imagesfileendascii[] = ".ascii";
const char   pixsrc_constants::fitsfile[] = ".fits";
const char   pixsrc_constants::psffile[] = ".psf.fits";
// integrating gaussian up to this multiple of sigma leads to .99 integral
const double pixsrc_constants::cutoff    = 2.575829303548901;
// conversion from FWHM of gaussian dist. to sigma
const double pixsrc_constants::fwhm2sigma= 2.3548200450309493;
const double pixsrc_constants::mm_stepsize = 0.5;
const double pixsrc_constants::mm_stepsize_area = 0.5;
const PS_SIT    pixsrc_constants::precision0 = 7;
const PS_SIT    pixsrc_constants::regspeed = 1;
const PS_SIT    pixsrc_constants::longstringsize = 1000;
const PS_SIT    pixsrc_constants::shortwaitms = 1000;
const PS_SIT    pixsrc_constants::longwaitms = 100000;

// Jonathan Richard Shewchuk's Triangle code
// these are flags for the delaunay triangulation
const char   pixsrc_constants::triswitchesnominangle[] = "zQenN";
const char   pixsrc_constants::triswitchesholefilling[] = "zQNp";
