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



#ifndef PIXSRC_LENSVAR_HPP_
#define PIXSRC_LENSVAR_HPP_

#define LENSVAR pixsrc_lensvar::

#include <vector>
#include "pixsrc_matrix.hpp"
#include "pixsrc_triangle.h"
#include "pixsrc_inputdata.hpp"
#include <pthread.h>

using std::vector;

// forward declaration for dataandvars
struct lensvar;

// dataandvars contains both inputdata and lensvar so that the struct can be passed
// to multihreaded functions, which only takes one argument of type void.
struct dataandvars
{
    inputdata *data;
    commoninputdata *cdata;
    lensvar *vars;
};

// for multithreading
enum ps_thread_name
{
    ps_thread_createc, ps_thread_detc, ps_thread_deta, ps_thread_ed,
    ps_thread_es, ps_thread_sizepenalty, ps_thread_penalty, ps_thread_triangulate,
    ps_thread_printimageplane, ps_thread_printsource, ps_thread_getmag, ps_thread_evidence,
    ps_thread_numpthreads
};

// contains image-specific data
struct lensvar
{
    /**********************************************/
    ////////////// START used by all ///////////////
    /**********************************************/

    // the pixsrc instance
    void *pixsrc_class;

    const static PS_SIT numpthreads = ps_thread_numpthreads;
    pthread_mutex_t pthreadslock[numpthreads];
    char pthreadstracker[numpthreads];
    pthread_t pthreads[numpthreads];

    // used for calculating interpolation errors
    double *bestanalyticsrc;
    void *tpsfitptr;
    void *tpsfit;

    char endedearly;
    char modelrejected;
    char srcopt;
    bool shapeletopt;
    bool terminatelensing;

    dataandvars *dav;

    double dofreduction; // number of effective free parameters in the source inversion

    PS_SIT imagenumber;
    PS_SIT tracker;
    double magger; // magnification
    double minx, maxx, miny, maxy;
    double rangex, rangey;
    double lambda1, ed, es, logdeta, logdetc;
    double lambda_lens;
    double variance;
    double datavariance; // may be different from "variance" for uv-plane data
    double* penalties;
    PS_SIT lonr;
    PS_SIT truelonr;
    PS_SIT lonc;
    PS_SIT *r4r;
    PS_SIT *r4rback;
    char *fallsinsrc;
    double *newloc;
    double *srcloc;

    // uv plane stuff
    PS_SIT lonr4stat;
    pixsrc_vector *residual4stat;
    pixsrc_vector *data4stat;
    pixsrc_matrix *lo4stat;

    // vectors and matrices
    pixsrc_vector *data;
    pixsrc_vector *mps;
    pixsrc_vector *lensedmps;
    pixsrc_vector *residual;
    pixsrc_vector *uv_residual;
    pixsrc_vector *lod;
    pixsrc_vector *interperr;
    pixsrc_vector *ieanalytic;
    pixsrc_vector *psf_analytic;
    pixsrc_vector *lensedmpsnobo;
    pixsrc_vector *srcexactvec;
    pixsrc_vector *ncv; // noise covariance vector
    pixsrc_matrix *lensingoperator;
    pixsrc_matrix *lensingoperatornobo;
    pixsrc_matrix *blurringoperator;
    pixsrc_matrix *b1;
    pixsrc_matrix *h1;
    pixsrc_matrix *h1x;
    pixsrc_matrix *h1y;
    pixsrc_matrix *c1;
    pixsrc_matrix *a1;

    // pointers for memory of vectors and matrices
    // these must of same type as class instantiated on memory (and not void*)
    // so that MEMORY ps_malloc can simply be passed the number
    // of objects to be instantiated instead of object size
    pixsrc_vector *dataptr;
    pixsrc_vector *mpsptr;
    pixsrc_vector *lensedmpsptr;
    pixsrc_vector *residualptr;
    pixsrc_vector *uv_residualptr;
    pixsrc_vector *lodptr;
    pixsrc_vector *interperrptr;
    pixsrc_vector *ieanalyticptr;
    pixsrc_vector *psf_analyticptr;
    pixsrc_vector *lensedmpsnoboptr;
    pixsrc_vector *srcexactvecptr;
    pixsrc_vector *ncvptr;
    pixsrc_matrix *lensingoperatorptr;
    pixsrc_matrix *lensingoperatornoboptr;
    pixsrc_matrix *blurringoperatorptr;
    pixsrc_matrix *b1ptr;
    pixsrc_matrix *h1ptr;
    pixsrc_matrix *h1xptr;
    pixsrc_matrix *h1yptr;
    pixsrc_matrix *c1ptr;
    pixsrc_matrix *a1ptr;


    double wcsinfo[9];
    double **newloc_sslo; // new locations subsampled
    PS_SIT **sspointer; // points from subsampled data pixels to the source triangle containing them.
    PS_SIT *pointer; // for subsample=1 case

    double **newloc_ssas; // for analytic source subsampling

    // fatal error flags
    // 0  => everything ok
    // 1  => error finding regularization strength
    // Other errors below yet to be filled out (mainly using reports from UMFPACK)
    // 2  =>
    // 4  =>
    // 8  =>
    // 16 =>
    PS_SIT fatalerror;

    // nonfatal error flags
    // 0  => everything okay
    PS_SIT nonfatalerror;

    /**********************************************/
    ////////////// END used by all /////////////////
    /**********************************************/

    /**********************************************/
    ////////////// START used by multiple //////////
    /**********************************************/

    // START used by adaptive and irrcart

    // Jonathan Richard Shewchuk's Triangle code
    struct triangulateio *triin;
    struct triangulateio *triout;

    double *convexhull;
    char need2retriangulate;

    // END used by adaptive and irrcart

    // START used by cart and irrcart

    PS_SIT srcx, srcy, numsrcpoints;
    double reduction;
    double zeroethgridsize;
    PS_SIT *c4c;
    PS_SIT *c4cback;

    // END used by cart and irrcart

    /**********************************************/
    ////////////// END used by multiple ////////////
    /**********************************************/

    /**********************************************/
    ////////////// START from lensingcartesian /////
    /**********************************************/

    double ctrx, ctry;

    /**********************************************/
    ////////////// END from lensingcartesian ///////
    /**********************************************/

    /**********************************************/
    ////////////// START from lensingirrcart ///////
    /**********************************************/

    double *srclocall;
    double maxdensity ,mindensity, rangedensity;
    //double match2mag1;
    PS_SIT maxgridlevel;
    double levelmag1;
    double *magnification;
    vector< vector< vector< vector<bool> > > > *grid;
    vector< vector< vector< vector<PS_SIT > > > > *gridpointer;
    double xoffset,yoffset;

    /**********************************************/
    ////////////// END from lensingirrcart /////////
    /**********************************************/

    /**********************************************/
    ////////////// START from lensingadaptive //////
    /**********************************************/

    //vector<PS_SIT> r4c; // which data position corresponds to a given source position
    //vector<PS_SIT> c4r; // which source position corresponds to a given data potision

    /**********************************************/
    ////////////// END from lensingadaptive ////////
    /**********************************************/

    double shapelet_ctr[2];
    double shapelet_scale;
    PS_SIT num_shapelets1;
    PS_SIT num_shapelets2;
    PS_SIT numberofshapelets;

};

#endif
