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



//
// contains c++ wrappers for external functions called from within pixsrc
//

#ifndef PIXSRC_EXTERNAL_HPP_
#define PIXSRC_EXTERNAL_HPP_

#define EXTERNAL pixsrc_external::

#include "pixsrc_external_wcs.hpp"
#include <ctime>

class pixsrc_external
{

public:

    // from blas/lapack

    static double            ps_ddot_         ( PS_SIT *N, double*, PS_SIT *incx, double*, PS_SIT *incy);
    static void              ps_dgesvd_       ( const char* jobu, const char* jobvt,
                                                const PS_SIT* M, const PS_SIT* N, double* A,
                                                const PS_SIT* lda, double* S, double* U,
                                                const PS_SIT* ldu, double* VT, const PS_SIT* ldvt,
                                                double* work, const PS_SIT *lwork, const PS_SIT* info );
    static void              ps_dgetrf_       ( PS_SIT* M, PS_SIT *N, double* A, PS_SIT* lda,
                                                int *IPIV, PS_SIT* INFO                            );
    static void              ps_dsytrf_       ( char* UPLO, PS_SIT* N, double* A, PS_SIT* LDA,
                                                int *IPIV, double* WORK, PS_SIT* lwork, PS_SIT* INFO  );
    static void              ps_ssytrf_       ( char* UPLO, PS_SIT* N, float* A, PS_SIT* LDA,
                                                int *IPIV, float* WORK, PS_SIT* lwork, PS_SIT* INFO  );
    static void              ps_dgetri_       ( PS_SIT* N, double* A, PS_SIT* lda, int *IPIV,
                                                double* WORK, PS_SIT* lwork, PS_SIT* INFO             );

    static void              ps_dgesv_        (PS_SIT* n, PS_SIT* nrhs, double* a, PS_SIT* lda,
                                               int *IPIV, double* b, PS_SIT* ldb, PS_SIT* info);

    static void              ps_dgetrs_       (char *TRANS, PS_SIT *N, PS_SIT *NRHS, double *A,
                                               PS_SIT *LDA, int *IPIV, double *B, PS_SIT *LDB, PS_SIT *INFO);
    static void              ps_dsytrs_       (char *UPLO, PS_SIT *N, PS_SIT *NRHS, double *A,
                                               PS_SIT* LDA, int *IPIV, double *B, PS_SIT *LDB, PS_SIT *INFO);
    static void              ps_ssytrs_       (char *UPLO, PS_SIT *N, PS_SIT *NRHS, float *A,
                                               PS_SIT* LDA, int *IPIV, float *B, PS_SIT *LDB, PS_SIT *INFO);
    static void              ps_dgemm_        ( const char *transa, const char *transb,
                                                PS_SIT *l, PS_SIT *n, PS_SIT *m, double *alpha,
                                                const double *a, PS_SIT *lda, double *b, PS_SIT *ldb,
                                                double *beta, double *c, PS_SIT *ldc                 );
    static void              ps_sgemm_        ( const char *transa, const char *transb,
                                                PS_SIT *l, PS_SIT *n, PS_SIT *m, float *alpha,
                                                const float *a, PS_SIT *lda, float *b, PS_SIT *ldb,
                                                float *beta, float *c, PS_SIT *ldc                 );

    static void              ps_dgemv_        ( const char *transa, PS_SIT *m, PS_SIT *n, double *alpha,
                                                const double *a, PS_SIT *lda, double *x, PS_SIT *incx,
                                                double *beta, double *y, PS_SIT *incy                );

#ifdef PS_HAVE_GRAVLENS
    // from gravlens
    static void              ps_potdefmag    ( double, double,  double*, double*,
                                               double*, double*, double*, double*, double*   );
    static void              ps_find_img     (double u, double v, PS_SIT *numimg, double ***imgs);
#endif
#ifdef PS_HAVE_TRIAXIAL
    // these use triaxial halo models (when not linking with gravlens/lensmodel)
    static void              ps_potdefmag    ( double, double,  double*, double*,
                                               double*, double*, double*, double*, double*, double**, time_t**, double**, PS_SIT);
    static void              ps_find_img     (double u, double v, PS_SIT *numimg, double ***imgs);
#endif

    // from WCS library

    const static PS_SIT         ps_WCS_J2000;
    const static PS_SIT         ps_WCS_B1950;
    const static PS_SIT         ps_WCS_ECLIPTIC;
    const static PS_SIT         ps_WCS_GALACTIC;
    static void              ps_setsys       ( PS_SIT                                           );
    static char*             ps_GetFITShead  ( char*, PS_SIT                                    );

    // the ps_WorldCoor struct is defined in "pixsrc_external_wcs.hpp"

    static struct ps_WorldCoor* ps_GetFITSWCS   ( char*, char*, PS_SIT, double*, double*, double*,
                                                  double*, double*, PS_SIT*, PS_SIT*, PS_SIT*, double*   );
    static void              ps_wcsc2pix     ( ps_WorldCoor*, double, double, char*,
                                               double*, double*, PS_SIT*                        );
    static void              ps_pix2wcs      ( ps_WorldCoor*, double, double, double*, double*  );
    static PS_SIT               ps_nowcs        ( ps_WorldCoor*                                    );
    static PS_SIT               ps_iswcs        ( ps_WorldCoor*                                    );
    static void              ps_wcsininit   ( ps_WorldCoor*, char*                             );
    static void              ps_wcsoutinit   ( ps_WorldCoor*, char*                             );
    static void              ps_wcscsys      ( char*                                         );
    static ps_WorldCoor*     ps_wcsinit      ( char*                                         );
    static void              ps_wcsfree      (ps_WorldCoor*);
    static ps_WorldCoor*     ps_wcskinit     (PS_SIT nxpix, PS_SIT nypix, char *ctype1, char *ctype2, double crpix1, double crpix2, double crval1, double crval2, double *cd, double cdelt1, double cdelt2, double crota, PS_SIT equinox, double epoch);

};

#endif
