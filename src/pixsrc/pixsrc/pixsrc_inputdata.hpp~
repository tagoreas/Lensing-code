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



#ifndef PIXSRC_INPUTDATA_HPP_
#define PIXSRC_INPUTDATA_HPP_

#define INPUTDATA pixsrc_inputdata::

#include "pixsrc_external_wcs.hpp"
#include <stdbool.h>
#include <fstream>
#include <pthread.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_rng.h>
#include <string>
#include <vector>

struct gsl_mdm_params;
struct inputdata;

typedef struct
{
    std::ofstream *stream;
    std::ofstream *streamptr;
    pthread_mutex_t *lock;
    pthread_mutex_t *lockptr;

} outstreamstruct;

typedef struct
{
    char *name;
    PS_SIT busid;
    PS_SIT locid;
    double gpubenchmark;
    double cpubenchmark;
    void *handle;
    void *handle_cublas;
    pthread_mutex_t *cusparse_lock;
    pthread_mutex_t *cublas_lock;
    void *descr;
    void *descr_sym;
    void *descr_tri_upper;
    bool ignoreme;
    PS_SIT num_kernels;
    //pthread_mutex_t memcpy_mutex;

} pixsrccudastruct;

typedef struct
{
    PS_SIT *printvec0;

    struct gsl_mdm_params *params;
    gsl_vector *x;

    double meandist2;
    double reg;
    PS_SIT    numvary;
    PS_SIT    dim;
    PS_SIT    nrhs;
    PS_SIT    LDA;
    PS_SIT    LDB;
    PS_SIT    info;
    int/*PS_SIT*/    *ipiv;
    double *tps_grid;
    double *tps_mat;
    double *tps_mat_fac;
    double *b_vec;
    double *x_vec;
    char    trans;
    double pert_region[8];
    PS_SIT *index;

    // reference perturbations
    PS_SIT *dim_ref;
    PS_SIT *numvary_ref;
    double **x_vec_ref;
    double **tps_grid_ref;

} npl_struct;


typedef struct
{
    PS_SIT numcall;  // counts number of times pixsrc has been run
    PS_SIT numimages; // total number of images
    PS_SIT numthreads; // total number of computer processors to use (suggestion really)
    char *basename;
    char *coordsysorig;
    char *coordsysfinal;
    pthread_mutex_t *potdefmagmutex;
    pthread_mutex_t *wcsmutex;
    pthread_mutex_t *print2screenmutex;
    pthread_attr_t *attrjoinable;
    pthread_attr_t *attrdetached;

    pixsrccudastruct *gpudevices;
    PS_SIT numgpudevices;

    // non-parametric lens model stuff
    PS_SIT npl;
    double npl_reg_parms[2]; // npl reg strength and stepsize
    double npl_stepsize;
    double npl_ftolsize;
    PS_SIT    npl_num_ref_pot;
    char **npl_filename;
    double npl_ctr[2];
    double npl_size[2];
    PS_SIT npl_num_pts[2];

    // GNU GSL random number generator
    const gsl_rng_type *ps_gsl_ran_T;
    gsl_rng *ps_gsl_ran_r;
    
    // lens model parameters (if not linking with lensmodel and using triaxial power laws)
    double **tlmparms;
    time_t **tlmtime;
    double **tlmenvirogals;

} commoninputdata;

struct inputdata
{
    //// extlengths records dimensions of arrays, where marked with 'x' (even if I think I know them)
    ////
    //// 0  ->
    //// 1  -> cartdetails[x]
    //// 2  -> srcinputcircle[x][]
    //// 3  -> srcinputcircle[][x]
    //// 4  -> traceparams[x]
    //// 5  -> penaltymatrix[][x][]
    //// 6  -> penaltymatrix[][][x]
    //// 7  -> invwcsinfo[x]
    //// 8  -> data[x]
    //// 9  -> psf[x][]
    //// 10 -> psf[][x]
    //// 11 -> mmimages[x][]
    //// 12 -> mmimages[][x]
    //// 13 -> penaltymatrix[x][][]
    //// 14 -> usersetsrc[x]
    //// 15 -> wcsinfo[x]
    //// 16 -> imagemasks[x]
    //// 17 -> shapeintmasks[x][]
    //// 18 -> shapeintmasks[][x]
    //// 19 -> chi2mask[x]

    double arc2pix; // pixels per arcsecond
    double pix2arc; // arcseconds per pixel

    PS_SIT numgpu2use;
    PS_SIT *gpu2use;

    char   *hacksrcfn; // hacksrc filename
    double *hacksrc;   // user-supplied values for source vector

    double *usersetsrc; // source set by user and not by inversion of lens equation
    PS_SIT srcrestart;
    PS_SIT srcnvary;
    PS_SIT numsrcbounds;
    double **srcbounds;
    double *srcstepsizes;
    PS_SIT guess_src_position;
    PS_SIT num_vec_src;
    PS_SIT *vector_src_numpix;
    double **vector_src_pos;
    double **vector_src_flux;

    double srcftol;
    char irscheme;
    PS_SIT srcexact;

    PS_SIT *extlengths;

    outstreamstruct *mags;
    outstreamstruct *magsptr;
    outstreamstruct *traces;
    outstreamstruct *tracesptr;
    outstreamstruct *details;
    outstreamstruct *detailsptr;

    bool doreconstruction;
    bool nopsf;
    bool psffromfile;
    PS_SIT psf_oversample;
    double fwhm_weight;
    bool debug;
    bool fatalwarn;
    PS_SIT  printvec;
    bool printdetails;
    bool onlychi;
    bool noevinochi;
    bool useall;
    bool findsisterimages;
    char rmpix;

    // 0 = bad pixel, 1 = outside mask, 2 = inside mask
    char *imagemasks;
    char *chi2mask;
    char **mmimages;
    double **shapeintmask;
    double **shapeintmaskminmax;

    // Jonathan Richard Shewchuk's Triangle code
    // this does delaunay triangulation, gives triangle vertex positions,
    // gives segment positions, and gives boundary segments
    char *triswitches;

    PS_SIT gridtype;
    char reg;
    PS_SIT interperr;

    PS_SIT optfinder;
    PS_SIT regaccuracy;
    PS_SIT regorder;
    PS_SIT verbose;
    PS_SIT subsampling;
    PS_SIT *oldloc;
    double *oldloc_wcs;
    double *oldloc_arc;

    PS_SIT *num_mmimages;
    double ***mmborder;
    double ***mmborder_defl;

    npl_struct *npl_st;
    double levelshift;
    double *fillbadpix;
    char fillbadpixnoise;
    double majoraxis;
    double minoraxis;
    double angleaxis;
    double lambdaguess;
    double lambda_lens;
    double myvariance0;
    double *wcsinfo;
    double rollangle;
    double pra, pdec;
    double px, py, r1, r2;
    double mintriangleangle; // minimum angle in any triangle in triangulation
    PS_SIT magparams;

    double ***penaltymatrix;
    PS_SIT *penaltyquery;
    char **penaltynames;
    char **penaltyfilenames;
    // order of penalties is:
    // 0 brightness-weighted source size
    // 1 magnification
    // 2 axis ratio
    // 3 mismatched images
    // 4 convex hull size
    // 5 flux outside data images

    double *traceparams;
    double *cartdetails;
    double **srcinputcircle;
    PS_SIT imgx,imgy;
    PS_SIT ndp;
    double *data;
    char *name;
    char *print2screenname;
    double evidence;
    double **psf;
    PS_SIT noisemap;
    PS_SIT lowmem;

    // shapelet stuff
    PS_SIT num_shapelets1;
    PS_SIT num_shapelets2;
    PS_SIT numberofshapelets;
    PS_SIT use_shapelets;
    double fixedshapeletparms[4];
    PS_SIT shapeletpixsplit[2];

    // rotate and translate coordinate system
    double rottrans[6];

    // construct lensing operation by ray-tracing from lens to source plane
    // or source to lens plane
    PS_SIT raydirection;
    double regftol; // regularization tolerance
    PS_SIT precision; // precision for printing doubles to files
    PS_SIT fullmag; // compute magnification using all pixels in image plane
    double *srcpixelscale; // source pixel scale (for printing FITS files)
    PS_SIT srcpixelscale0; // default value

    // interferometric data flag
    PS_SIT is_uvdata;
    char *uv_filename; // measurement sets
    char *uv_model_filename; // measurement sets
    double *uv_oldloc; // data positions
    double *uv_data;   // data
    void   *uv_transform_vec, *uv_transform_vec_ptr; // convert image plane to uv plane
    void   *uv_transform_mat, *uv_transform_mat_ptr; // convert image plane to uv plane
    double uv_transform_scalar;                      // convert image plane to uv plane
    PS_SIT    uv_ndp;
    double uv_pointing[2];
    double uv_cutoff[2];
    PS_SIT uv_mode;
    double *uv_newvariance;
    PS_SIT uv_ndporig; // number of uv data pts, before uv cutoff or combining
    PS_SIT *uv_idx; // keeps track of original index of uv data points
    double uv_numnoisesamples; // double because negative numbers fix noise
    double uv_deltapts; // number of data points per TPS splines
    PS_SIT uv_maxintndp; // maximum number of integrated data points
    double uv_taper; // taper scale
    PS_SIT transmode;
    double uv_newpixelscale;
    PS_SIT uv_matvecsca_exist;
    PS_SIT uv_padzeros;
    PS_SIT uv_rbftype;
    double uv_rbfscale;
    PS_SIT uv_del1, uv_del2; // submatrix multiliication params
    PS_SIT uv_pbeam; // apply primary beam 


    // vector source stuff (similar to analytic source stuff)
    struct triangulateio **vector_src_triout;
    double *rotation_axis;

    double *invwcsinfo;

    // ps_Worldcoor comes from "pixsrc_external_wcs.hpp"
    ps_WorldCoor *wcs;

    // *inputdata data below points to first image
    inputdata *data_;
};

struct gsl_mdm_params
{
    inputdata *data;
    commoninputdata *cdata;
};

struct gsl_mdm
{
    gsl_mdm() : T(gsl_multimin_fminimizer_nmsimplex2rand),
                s(NULL), ss(NULL), iter(0) {}

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *ss;
    gsl_multimin_function minex_func;
    size_t iter;
    PS_SIT status;
    double size;
};

// these class members will hold the index of each input parameter
class ps_parms_ind_class
{
public:
    PS_SIT ps_parm_ind_debug;
    PS_SIT ps_parm_ind_regorder;
    PS_SIT ps_parm_ind_regstrength;
    PS_SIT ps_parm_ind_coorsys;
    PS_SIT ps_parm_ind_grid;
    PS_SIT ps_parm_ind_psf;
    PS_SIT ps_parm_ind_sisterpix;
    PS_SIT ps_parm_ind_oversample;
    PS_SIT ps_parm_ind_noise;
    PS_SIT ps_parm_ind_penalty6;
    PS_SIT ps_parm_ind_statistic;
    PS_SIT ps_parm_ind_magnification;
    PS_SIT ps_parm_ind_verbosity;
    PS_SIT ps_parm_ind_images;
    PS_SIT ps_parm_ind_varystrength;
    PS_SIT ps_parm_ind_threads;
    PS_SIT ps_parm_ind_penalty1;
    PS_SIT ps_parm_ind_penalty2;
    PS_SIT ps_parm_ind_penalty3;
    PS_SIT ps_parm_ind_penalty4;
    PS_SIT ps_parm_ind_penalty5;
    PS_SIT ps_parm_ind_details;
    PS_SIT ps_parm_ind_fillbadpix;
    PS_SIT ps_parm_ind_regaccuracy;
    PS_SIT ps_parm_ind_minangle;
    PS_SIT ps_parm_ind_reg;
    PS_SIT ps_parm_ind_rmpix;
    PS_SIT ps_parm_ind_hacksrc;
    PS_SIT ps_parm_ind_src;
    PS_SIT ps_parm_ind_srcbound;
    PS_SIT ps_parm_ind_srcstepsize;
    PS_SIT ps_parm_ind_srcrestart;
    PS_SIT ps_parm_ind_interperr;
    PS_SIT ps_parm_ind_irscheme;
    PS_SIT ps_parm_ind_srcftol;
    PS_SIT ps_parm_ind_cuda;
    PS_SIT ps_parm_ind_optfinder;
    PS_SIT ps_parm_ind_nonparamlens;
    PS_SIT ps_parm_ind_reglens;
    PS_SIT ps_parm_ind_nplstepsize;
    PS_SIT ps_parm_ind_srcctrguess;
    PS_SIT ps_parm_ind_nplftol;
    PS_SIT ps_parm_ind_npltpsreg;
    PS_SIT ps_parm_ind_shapelets;
    PS_SIT ps_parm_ind_shapeparms;
    PS_SIT ps_parm_ind_shapepixsplit;
    PS_SIT ps_parm_ind_penaltyname;
    PS_SIT ps_parm_ind_dataname;
    PS_SIT ps_parm_ind_rottrans;
    PS_SIT ps_parm_ind_noisemap;
    PS_SIT ps_parm_ind_raydirection;
    PS_SIT ps_parm_ind_regftol;
    PS_SIT ps_parm_ind_outprecision;
    PS_SIT ps_parm_ind_fullmag;
    PS_SIT ps_parm_ind_pixelscalesrc;
    PS_SIT ps_parm_ind_uvdata;
    PS_SIT ps_parm_ind_uvmodelpos;
    PS_SIT ps_parm_ind_lowmem;
    PS_SIT ps_parm_ind_uvnoise;
    PS_SIT ps_parm_ind_uvtaper;
    PS_SIT ps_parm_ind_uvcutoff;
    PS_SIT ps_parm_ind_uvnoisesample;
    PS_SIT ps_parm_ind_uvnewpixelscale;
    PS_SIT ps_parm_ind_uvrbf;
    PS_SIT ps_parm_ind_uvpadzeros;
    PS_SIT ps_parm_ind_uvmatrixsize;
    PS_SIT ps_parm_ind_uvpbeam;
    PS_SIT ps_parm_ind_srcexact;
    PS_SIT ps_parm_ind_fatalwarn;
};

// these enum members will hold the category of each input parameter
enum ps_parms_cat_enum
{
    ps_parm_cat_general,
    ps_parm_cat_dataNgridding,
    ps_parm_cat_shapelets,
    ps_parm_cat_regularization,
    ps_parm_cat_analyticsrc,
    ps_parm_cat_uvdata,
    ps_parm_cat_gpu,
    ps_parm_cat_penalty,
    ps_parm_cat_potentialperturbation,
    ps_parm_cat_interperr,
    ps_parm_cat_numcategories
};

struct ps_parms_struct
{
    ps_parms_ind_class *pindex;
    std::string scategory[10];
    ps_parms_cat_enum   category; // category
    std::string          sname;    // name
    std::vector <double> fvalue;   // default value(s)
    std::string          svalue;   // default value
    std::vector <std::string> qdescr;   // quick decription
    std::vector <std::string> notes;    // detailed description
    std::vector< std::vector< std::vector<std::string> > > entries;  // valid inputs
};

#endif
