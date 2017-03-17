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
// This helps future coders keep track of where functions, definitions,
// variables, etc are all coming from.
//

#include "pixsrc_external.hpp"

// from blas/lapack
extern "C" void              dgesvd_       ( const char* jobu, const char* jobvt,
                                            const int* M, const int* N, double* A,
                                            const int* lda, double* S, double* U,
                                            const int* ldu, double* VT, const int* ldvt,
                                            double* work, const int *lwork, const int* info );
extern "C" double            ddot_         ( int *N, double*, int *incx, double*, int *incy);
extern "C" void              dgetrf_       ( int* M, int *N, double* A, int* lda,
                                             int* IPIV, int* INFO                            );
extern "C" void              dsytrf_       ( char* UPLO, int *N, double* A, int* LDA,
                                             int* IPIV, double* WORK, int* lwork, int* INFO  );
extern "C" void              ssytrf_       ( char* UPLO, int *N, float* A, int* LDA,
                                             int* IPIV, float* WORK, int* lwork, int* INFO  );
extern "C" void              dgetri_       ( int* N, double* A, int* lda, int* IPIV,
                                             double* WORK, int* lwork, int* INFO             );

extern "C" void              dgesv_        ( int* n, int* nrhs, double* a, int* lda,
                                             int* ipiv, double* b, int* ldb, int* info);

extern "C" void              dgetrs_(char *TRANS, int *N, int *NRHS, double *A,
                                     int *LDA, int *IPIV, double *B, int *LDB, int *INFO );
extern "C" void              dsytrs_(char *UPLO, int *N, int *NRHS, double *A, int* LDA,
                                     int *IPIV, double *B, int *LDB, int *INFO );
extern "C" void              ssytrs_(char *UPLO, int *N, int *NRHS, float *A, int* LDA,
                                     int *IPIV, float *B, int *LDB, int *INFO );
extern "C" void              dgemm_        ( const char *transa, const char *transb,
                                             int *l, int *n, int *m, double *alpha,
                                             const double *a, int *lda, double *b, int *ldb,
                                             double *beta, double *c, int *ldc                  );
extern "C" void              sgemm_        ( const char *transa, const char *transb,
                                             int *l, int *n, int *m, float *alpha,
                                             const float *a, int *lda, float *b, int *ldb,
                                             float *beta, float *c, int *ldc                  );

extern "C" void              dgemv_        ( const char *transa, int *m, int *n, double *alpha,
                                             const double *a, int *lda, double *x, int *incx,
                                             double *beta, double *y, int *incy                 );

#ifdef PS_HAVE_GRAVLENS
// from Chuck's code
extern "C" void              potdefmag    ( double x, double y,  double *ig1,
                                            double *defx, double *defy, double *ig2,
                                            double *ig3, double *ig4, double *ig5        );
extern "C" void              find_img     ( double u, double v,  int *numimg, double ***imgs);
#endif
#ifdef PS_HAVE_TRIAXIAL
// these use triaxial halo models (when not linking with gravlens/lensmodel)
extern "C" void              potdefmag    ( double x, double y,  double *ig1,
                                            double *defx, double *defy, double *ig2,
                                            double *ig3, double *ig4, double *ig5, double **parms, time_t **tlmtime, double **envirogals, PS_SIT imagenumber);
extern "C" void              find_img     ( double u, double v,  int *numimg, double ***imgs);
#endif

// from WCS library

extern "C" char             *GetFITShead  ( char*, int                                   );
extern "C" void              setsys       ( int                                          );
extern "C" struct ps_WorldCoor *GetFITSWCS   ( char*, char*, int,
                                               double*, double*, double*,
                                               double*, double*, int*,
                                               int*, int*, double*                          );
extern "C" void              wcsc2pix     ( ps_WorldCoor*, double, double, char*,
                                            double*, double*, int*                       );
extern "C" void              pix2wcs      ( ps_WorldCoor*, double, double, double*, double* );
extern "C" int               nowcs        ( ps_WorldCoor*                                   );
extern "C" int               iswcs        ( ps_WorldCoor*                                   );
extern "C" void              wcsininit   ( ps_WorldCoor*, char*                            );
extern "C" void              wcsoutinit   ( ps_WorldCoor*, char*                            );
extern "C" int               wcscsys      ( char*                                        );
extern "C" ps_WorldCoor*     wcsinit      ( char*                                        );
extern "C" void              wcsfree      ( ps_WorldCoor*                                   );
extern "C" ps_WorldCoor*     wcskinit     (int nxpix, int nypix, char *ctype1, char *ctype2, double crpix1, double crpix2, double crval1, double crval2, double *cd, double cdelt1, double cdelt2, double crota, int equinox, double epoch);


//
// pixsrc wrappers
//

double pixsrc_external::ps_ddot_   (PS_SIT *N, double *x, PS_SIT *incx, double *y, PS_SIT *incy)
{
    int Ni    = (int)(*N);
    int incxi = (int)(*incx);
    int incyi = (int)(*incy);
    return ddot_   (&Ni, x, &incxi, y, &incyi);
}
void pixsrc_external::ps_dgetrs_   (char *TRANS, PS_SIT *N, PS_SIT *NRHS, double *A,
                                    PS_SIT *LDA, int *IPIV, double *B, PS_SIT *LDB, PS_SIT *INFO )
{
    int Ni    = (int)(*N);
    int NRHSi = (int)(*NRHS);
    int LDAi  = (int)(*LDA);
    int LDBi  = (int)(*LDB);
    int INFOi = (int)(*INFO);
    dgetrs_   (TRANS, &Ni, &NRHSi, A, &LDAi, IPIV, B, &LDBi, &INFOi );
    *INFO = (PS_SIT)INFOi;
}
void pixsrc_external::ps_dsytrs_   (char *UPLO, PS_SIT *N, PS_SIT *NRHS, double *A, PS_SIT* LDA,
                                    int *IPIV, double *B, PS_SIT *LDB, PS_SIT *INFO )
{
    int Ni    = (int)(*N);
    int NRHSi = (int)(*NRHS);
    int LDAi  = (int)(*LDA);
    int LDBi  = (int)(*LDB);
    int INFOi = (int)(*INFO);
    dsytrs_   (UPLO, &Ni, &NRHSi, A, &LDAi, IPIV, B, &LDBi, &INFOi );
    *INFO = (PS_SIT)INFOi;
}
void pixsrc_external::ps_ssytrs_   (char *UPLO, PS_SIT *N, PS_SIT *NRHS, float *A, PS_SIT* LDA,
                                    int *IPIV, float *B, PS_SIT *LDB, PS_SIT *INFO )
{
    int Ni    = (int)(*N);
    int NRHSi = (int)(*NRHS);
    int LDAi  = (int)(*LDA);
    int LDBi  = (int)(*LDB);
    int INFOi = (int)(*INFO);
    ssytrs_   (UPLO, &Ni, &NRHSi, A, &LDAi, IPIV, B, &LDBi, &INFOi );
    *INFO = (PS_SIT)INFOi;
}

void pixsrc_external::ps_dgemm_    ( const char *transa, const char *transb,
                                     PS_SIT *l, PS_SIT *n, PS_SIT *m, double *alpha,
                                     const double *a, PS_SIT *lda, double *b, PS_SIT *ldb,
                                     double *beta, double *c, PS_SIT *ldc                )
{
    int li   = (int)(*l);
    int ni   = (int)(*n);
    int mi   = (int)(*m);
    int ldai = (int)(*lda);
    int ldbi = (int)(*ldb);
    int ldci = (int)(*ldc);
    dgemm_ (transa, transb, &li, &ni, &mi, alpha, a, &ldai, b, &ldbi, beta, c, &ldci);
}
void pixsrc_external::ps_sgemm_    ( const char *transa, const char *transb,
                                     PS_SIT *l, PS_SIT *n, PS_SIT *m, float *alpha,
                                     const float *a, PS_SIT *lda, float *b, PS_SIT *ldb,
                                     float *beta, float *c, PS_SIT *ldc                )
{
    int li   = (int)(*l);
    int ni   = (int)(*n);
    int mi   = (int)(*m);
    int ldai = (int)(*lda);
    int ldbi = (int)(*ldb);
    int ldci = (int)(*ldc);
    sgemm_ (transa, transb, &li, &ni, &mi, alpha, a, &ldai, b, &ldbi, beta, c, &ldci);
}

void pixsrc_external::ps_dgemv_    ( const char *transa, PS_SIT *m, PS_SIT *n, double *alpha,
                                     const double *a, PS_SIT *lda, double *x, PS_SIT *incx,
                                     double *beta, double *y, PS_SIT *incy               )
{
    int mi    = (int)(*m);
    int ni    = (int)(*n);
    int ldai  = (int)(*lda);
    int incxi = (int)(*incx);
    int incyi = (int)(*incy);
    dgemv_ (transa, &mi, &ni, alpha, a, &ldai, x, &incxi, beta, y, &incyi);
}

void pixsrc_external::ps_dgesvd_   ( const char* jobu, const char* jobvt,
				     const PS_SIT* M, const PS_SIT* N, double* A,
				     const PS_SIT* lda, double* S, double* U,
				     const PS_SIT* ldu, double* VT, const PS_SIT* ldvt,
				     double* work, const PS_SIT *lwork, const PS_SIT* info )
{
    int mi     = (int)(*M);
    int ni     = (int)(*N);
    int ldai   = (int)(*lda);
    int ldui   = (int)(*ldu);
    int ldvti  = (int)(*ldvt);
    int lworki = (int)(*lwork);
    int infoi  = (int)(*info);
    dgesvd_( jobu, jobvt, &mi, &ni, A, &ldai, S, U, &ldui, VT, &ldvti, work, &lworki, &infoi );
}


void pixsrc_external::ps_dgetrf_    ( PS_SIT* M, PS_SIT *N, double* A, PS_SIT* lda,
                                      int* IPIV, PS_SIT* INFO                            )
{
    int Mi    = (int)(*M);
    int Ni    = (int)(*N);
    int ldai  = (int)(*lda);
    int INFOi = (int)(*INFO);
    dgetrf_( &Mi, &Ni, A, &ldai, IPIV, &INFOi );
    *INFO = (PS_SIT)INFOi;
}
void pixsrc_external::ps_dsytrf_    ( char* UPLO, PS_SIT *N, double* A, PS_SIT* LDA,
                                      int* IPIV, double* work, PS_SIT *lwork, PS_SIT* INFO  )
{
    int Ni     = (int)(*N);
    int LDAi   = (int)(*LDA);
    int lworki = (int)(*lwork);
    int INFOi  = (int)(*INFO);
    dsytrf_( UPLO, &Ni, A, &LDAi, IPIV, work, &lworki, &INFOi );
    *INFO = (PS_SIT)INFOi;
}
void pixsrc_external::ps_ssytrf_    ( char* UPLO, PS_SIT *N, float* A, PS_SIT* LDA,
                                      int* IPIV, float* work, PS_SIT *lwork, PS_SIT* INFO  )
{
    int Ni     = (int)(*N);
    int LDAi   = (int)(*LDA);
    int lworki = (int)(*lwork);
    int INFOi  = (int)(*INFO);
    ssytrf_( UPLO, &Ni, A, &LDAi, IPIV, work, &lworki, &INFOi );
    *INFO = (PS_SIT)INFOi;
}
void pixsrc_external::ps_dgetri_    ( PS_SIT* N, double* A, PS_SIT* lda, int* IPIV,
                                      double* WORK, PS_SIT* lwork, PS_SIT* INFO             )
{
    int Ni     = (int)(*N);
    int ldai   = (int)(*lda);
    int lworki = (int)(*lwork);
    int INFOi  = (int)(*INFO);
    dgetri_( &Ni, A, &ldai, IPIV, WORK, &lworki, &INFOi );
    *INFO = (PS_SIT)INFOi;
}

void pixsrc_external::ps_dgesv_    ( PS_SIT* n, PS_SIT* nrhs, double* a, PS_SIT* lda,
                                     int* ipiv, double* b, PS_SIT* ldb, PS_SIT* info)
{
    int ni    = (int)(*n);
    int nrhsi = (int)(*nrhs);
    int ldai  = (int)(*lda);
    int ldbi  = (int)(*ldb);
    int infoi = (int)(*info);
    dgesv_ (&ni, &nrhsi, a, &ldai, ipiv, b, &ldbi, &infoi);
    *info = (PS_SIT)infoi;
}

#ifdef PS_HAVE_GRAVLENS
void pixsrc_external::ps_potdefmag( double x, double y,  double *ig1,
                                    double *defx, double *defy, double *ig2,
                                    double *ig3, double *ig4, double *ig5          )
{
    potdefmag( x, y, ig1, defx, defy, ig2, ig3, ig4, ig5 );
}
#endif
#ifdef PS_HAVE_TRIAXIAL
void pixsrc_external::ps_potdefmag( double x, double y,  double *ig1,
                                    double *defx, double *defy, double *ig2,
                                    double *ig3, double *ig4, double *ig5, double **parms, time_t **tlmtime, double **envirogals, PS_SIT imagenum)
{
    potdefmag( x, y, ig1, defx, defy, ig2, ig3, ig4, ig5, parms, tlmtime, envirogals, imagenum);
}
#endif

void pixsrc_external::ps_find_img (double u, double v, PS_SIT *numimg, double ***imgs)
{
    int tmp;
    find_img (u, v, &tmp, imgs);
    *numimg = (PS_SIT)tmp;
}

// the following WCS keywords are macros in wcs.h .. I've copied them here
// instead of including wcs.h for reasons listed in "pixsrc_external_wcs.hpp"
// the variables "ps_external_wcs_*" are found in "pixsrc_external_wcs.hpp"
const PS_SIT pixsrc_external::ps_WCS_J2000    = (PS_SIT)ps_external_wcs_j2000;
const PS_SIT pixsrc_external::ps_WCS_B1950    = (PS_SIT)ps_external_wcs_b1950;
const PS_SIT pixsrc_external::ps_WCS_ECLIPTIC = (PS_SIT)ps_external_wcs_ecliptic;
const PS_SIT pixsrc_external::ps_WCS_GALACTIC = (PS_SIT)ps_external_wcs_galactic;

char* pixsrc_external::ps_GetFITShead( char *a, PS_SIT b )
{
    return GetFITShead( a, (int)b );
}
void pixsrc_external::ps_setsys( PS_SIT a )
{
    setsys( (int)a );
}
struct ps_WorldCoor* pixsrc_external::ps_GetFITSWCS( char *a, char *b, PS_SIT c,
                                                     double *d, double *e, double *f,
                                                     double *g, double *h, PS_SIT *i,
                                                     PS_SIT *j, PS_SIT *k, double *l          )
{
    int ii = (int)(*i);
    int ji = (int)(*j);
    int ki = (int)(*k);
    struct ps_WorldCoor* returnme = GetFITSWCS( a, b, (int)c, d, e, f, g, h, &ii, &ji, &ki, l );
    *i = (PS_SIT)ii;
    *j = (PS_SIT)ji;
    *k = (PS_SIT)ki;
    return returnme;
}
struct ps_WorldCoor* pixsrc_external::ps_wcsinit (char *a)
{
    return wcsinit (a);
}
void pixsrc_external::ps_wcsc2pix( ps_WorldCoor *a, double b, double c, char* d,
                                   double* e, double* f, PS_SIT* g                      )
{
    int gi = (int)(*g);
    wcsc2pix( a, b, c, d, e, f, &gi );
    *g = (PS_SIT)gi;
}
void pixsrc_external::ps_pix2wcs ( ps_WorldCoor *a, double b, double c, double *d, double *e )
{
    pix2wcs( a, b, c, d, e );
}
PS_SIT pixsrc_external::ps_nowcs     ( ps_WorldCoor *a                                   )
{
    return (PS_SIT) nowcs( a );
}
PS_SIT pixsrc_external::ps_iswcs     ( ps_WorldCoor *a                                   )
{
    return (PS_SIT) iswcs( a );
}
void pixsrc_external::ps_wcsfree     ( ps_WorldCoor *a                                   )
{
    return wcsfree (a);
}
void pixsrc_external::ps_wcsoutinit( ps_WorldCoor *a, char *b                             )
{
    wcsoutinit( a, b );
}
void pixsrc_external::ps_wcsininit( ps_WorldCoor *a, char *b                             )
{
    wcsininit( a, b );
}
void pixsrc_external::ps_wcscsys   ( char *a                                         )
{
    wcscsys( a );
}
ps_WorldCoor* pixsrc_external::ps_wcskinit (PS_SIT nxpix, PS_SIT nypix, char *ctype1, char *ctype2, double crpix1, double crpix2, double crval1, double crval2, double *cd, double cdelt1, double cdelt2, double crota, PS_SIT equinox, double epoch)
{
    return wcskinit ((int)nxpix, (int)nypix, ctype1, ctype2, crpix1, crpix2, crval1, crval2, cd, cdelt1, cdelt2, crota, (int)equinox, epoch);
}
