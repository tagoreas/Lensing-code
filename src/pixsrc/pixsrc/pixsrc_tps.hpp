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


#ifndef PIXSRC_TPS_HPP_
#define PIXSRC_TPS_HPP_

#define TPS pixsrc_tps

// This class creates a thin plate spline

#include "pixsrc_matrix.hpp"
#include "pixsrc_vector.hpp"

typedef enum
{
    ps_rbf_tps,
    ps_rbf_gaussian,
    ps_rbf_tps_weighted,
    ps_rbf_linear_weighted
} ps_rbf_type;

class pixsrc_tps
{
public:
    pixsrc_tps      (commoninputdata*, inputdata*, const double *x, const double *y, const double *z,
                     PS_SIT xstride, PS_SIT ystride, PS_SIT zstride, PS_SIT numpts_, ps_rbf_type, double parm, PS_SIT integrateme, PS_SIT addlinearterms, PS_SIT transmode);
    pixsrc_tps      (commoninputdata*, inputdata*, const double *x, const double *y, const double *z,
                     PS_SIT xstride, PS_SIT ystride, PS_SIT zstride, PS_SIT numtotalpts_, bool *is_control, ps_rbf_type, double parm, PS_SIT integrateme, PS_SIT addlinearterms, PS_SIT transmode);
    void init_exact (const double *x, const double *y, const double *z,
                     PS_SIT xstride, PS_SIT ystride, PS_SIT zstride, PS_SIT numpts_, PS_SIT addlinearterms);

    void set_regularization (double reg);
    void get_tps_weights ();
    double interpolate (double x, double y);
    void   interpolate (double *pos,  PS_SIT *r4r, VECTOR *img);
    void   interpolate (double **pos, PS_SIT *r4r, VECTOR *img, PS_SIT subsample);
    void   precompute_fourier (double *pos, PS_SIT ndp, double *fac2);
    void   get_fourier_weights (double *pos, MATRIX *weights, PS_SIT ndp);
    void   get_fourier_weights_col (PS_SIT col, double *pos, double *weights, PS_SIT ndp, double *fac2, PS_SIT start, PS_SIT end);
    void   cleanup ();
    PS_SIT    get_num_ctrl_pts ();
    ~pixsrc_tps();

    MATRIX *tpsmat, *tpsmatptr;
    MATRIX *tmat, *tmatptr;
    double fac1;
    double *xcoords, *ycoords;

private:

    static double ps_rbf_integrate_y (double y, void *x_);
    static double ps_rbf_integrate_x (double x, void *parms);

    commoninputdata *cdata_;
    inputdata *data_;
    PS_SIT integrate_quadrant_query (PS_SIT quadrant, double *lim);
    void tps_init (commoninputdata*, inputdata*, ps_rbf_type, double rbf_parm, PS_SIT integrateme, PS_SIT addlinearterms, PS_SIT transmode);
    //void set_precision (float*);
    void set_precision (double*);
    double kernel (double);
    double kernel_integral (double*);
    void fourier_kernel (double x, double y, double u, double v, double *wts, double *fac2);

    // wrappers for BLAS functions (to deal with single/double precision)
    typedef void (*ps_tps_gemm) (const char *transa, const char *transb,
                                 PS_SIT *l, PS_SIT *n, PS_SIT *m, double *alpha,
                                 const double *a, PS_SIT *lda, double *b, PS_SIT *ldb,
                                 double *beta, double *c, PS_SIT *ldc);
    typedef void (*ps_tps_sytrf) (char* UPLO, PS_SIT* N, double* A, PS_SIT* LDA,
                                  int*/*PS_SIT**/ IPIV, double* WORK, PS_SIT* lwork, PS_SIT* INFO);
    typedef void (*ps_tps_sytrs) (char *UPLO, PS_SIT *N, PS_SIT *NRHS, double *A,
                                  PS_SIT* LDA, int*/*PS_SIT**/ IPIV, double *B, PS_SIT *LDB, PS_SIT *INFO);
    ps_tps_gemm  my_pstpsgemm;
    ps_tps_sytrf my_pstpssytrf;
    ps_tps_sytrs my_pstpssytrs;
    ps_rbf_type  psrbftype;
    double gauss_sigma, ft_gauss_sigma, dxdyfac;

    void *gsl_workspace_tps_x, *gsl_workspace_tps_y;

    PS_SIT is_exact, integrate_rbf, uselinearterms;
    PS_SIT numpts, dim, numtotalpts;
    int/*PS_SIT*/ *ipiv;
    PS_SIT transmode;
    double meandist2, reg;
    double *tps_mat_fac, *x_vec, *b_vec;
    double *a_mat;
    double *x_rad_offset, *y_rad_offset; // WCS offsets (in radians) -- for Fourier kernel
    pixsrc_tps operator=(const pixsrc_tps &mdm);
    pixsrc_tps(const pixsrc_tps &mdm);
};

#endif
