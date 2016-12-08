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



#ifndef PIXSRC_SHAPELETS_OPERATIONS_HPP_
#define PIXSRC_SHAPELETS_OPERATIONS_HPP_

#define SHAPELETSOPERA pixsrc_shapelets_operations::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

class pixsrc_shapelets_operations
{

public:

    static void createlo_shapelets     (inputdata*, commoninputdata*, lensvar*);
    // overloaded - extra argument indicates basis functions are integrated over
    // the source plane
    static double get_hermite_sb       (inputdata*, commoninputdata*, lensvar*,
                                        double, double);
/*
  static void get_hermite_sb       (inputdata*, commoninputdata*, lensvar*,
  double, double, double, double *, double*);
*/
    static void  get_binomials      (inputdata*, commoninputdata*, lensvar*,
                                     double*, PS_SIT);
    static double get_flux_one_pixel (inputdata*, commoninputdata*, lensvar*,
                                      PS_SIT, double**);
    static double integrate_square (inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                    double *hermvals_square, PS_SIT numsh_1, PS_SIT numsh_2,
                                    double sh_ctrx, double sh_ctry, double sh_scale,
                                    double a, double b, double c, double d, VECTOR*, double*);
    static void hermite                (PS_SIT, double, double*);
/*
  static double get_laguerre_sb      (inputdata*, commoninputdata*, lensvar*,
  double, double);
*/

private:

    //static double src_recon_iter       (inputdata*, commoninputdata*, lensvar*);
    static void lo_set_hermite         (inputdata*, commoninputdata*, lensvar*);
    static void* lo_set_hermite_thread (void*);
/*
  static void laguerre               (PS_SIT, PS_SIT, double, double**);
  static void lo_set_laguerre        (inputdata*, commoninputdata*, lensvar*);
*/

    static void integrate_shapelets (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *tri, double *hermvals, double *hermvals_split_tri);
    static void integrate_triangle (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double c, double d, double m1, double m2, double b1, double b2, double *hermvals, double *hermvals_split_tri);
    static double i00 (inputdata*, commoninputdata*, double m, double c, double d, double b, double q);
    static double i00_calc (inputdata*, commoninputdata*, double m, double u1, double u2, double b, double q, double sign);
    static double j00 (double m, double x1, double x2, double b, double q);
    static double j10 (double m, double x1, double x2, double b, double q);
    static double j01 (double m, double x1, double x2, double b, double q);
    static double jn0 (PS_SIT nx, double c, double d, double m, double b, double q, double **jnn, double *phi_c, double *phi_d);
    static double jn0_uv (PS_SIT nx, double m, double x, double b, double q, double *phi);
    static double jn0_vdu (PS_SIT nx, double m, double b, double q, double **jnn);
    static double j0n (PS_SIT ny, double c, double d, double m, double b, double q, double **jnn, double *phi_mcb, double *phi_mdb);
    static double j0n_uv (PS_SIT ny, double m, double x, double b, double q, double *phi);
    static double j0n_vdu (PS_SIT ny, double m, double b, double q, double **jnn);
    static double j1n (PS_SIT ny, double m, double b, double q, double **jnn);
    static double jnxny (PS_SIT nx, PS_SIT ny, double c, double d, double m, double b, double q, double **jnn, double *phi_c, double *phi_d, double *phi_mcb, double *phi_mdb);
    static double jnxny_uv (PS_SIT nx, PS_SIT ny, double m, double x, double b, double q, double *phi, double *phi_mxb);
    static double jnxny_vdu (PS_SIT nx, PS_SIT ny, double m, double b, double q, double **jnn);
    static double in0 (PS_SIT nx, double c, double d, double m, double b, double q, double **inn, double *phi_c, double *phi_d);
    static double in0_uv (PS_SIT nx, double m, double x, double b, double q, double *phi);
    static double in0_vdu (PS_SIT nx, double m, double b, double q, double **inn);
    static double inxny (PS_SIT nx, PS_SIT ny, double **inn, double **jnn);
    static double i10 (double m, double x1, double x2, double b, double q);

};

#endif
