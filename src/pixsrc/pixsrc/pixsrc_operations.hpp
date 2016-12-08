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



#ifndef PIXSRC_OPERATIONS_HPP_
#define PIXSRC_OPERATIONS_HPP_

#define OPERA pixsrc_operations::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"
#include "pixsrc_constants.hpp"
#include <vector>
#include <string>

using std::vector;
using std::string;

class pixsrc_operations
{

public:

    // I'm unsure about using this function with multithreading without
    // using locks, but so far I'm only using it in one place, so no
    // lock needed :)
    // it'd be better just to make multiple gseeds
    static double randomgaussian (const gsl_rng*);


    // c-string handling functionsplanar
    static void   trim                       ( char*, PS_SIT*                                       );
    static PS_SIT    sizestring                 ( const char*                                       );
    static void   split                      ( const char*, const char*, char***, PS_SIT*           );
    static void   readfile                   ( const char*, char***, PS_SIT*, pthread_mutex_t *lock       );
    static bool   fileexists                 (string);
    static void   concatenate                ( const char**, PS_SIT, char**                         );
    static bool   isanumberint               ( string                                            );
    static bool   isanumberfloat             ( string                                            );

    static PS_SIT    coords2coordl              ( PS_SIT, PS_SIT                                          );
    static PS_SIT    coordl2coords              ( PS_SIT, PS_SIT                                          );
    static PS_SIT    rank                       ( double, double, double, double, double, double    );
    static void   linearinterpolatorcartesian( double, double, double, double*, PS_SIT*, double*    );
    static double getfitscardvalue           ( string                                            );
    static void   planarvalueinterpolation   ( double[2][3], double[3], double, double, double*  );
    static void   planarinterpolation3pts    ( double, double, double[2][3], double[3]           );
    static bool   equalsnpos                 ( PS_SIT                                               );
    static void   pthreadswait               ( lensvar*, PS_SIT, PS_SIT*                               );
    static void   pthreadswaitifstarted      ( lensvar*, PS_SIT, PS_SIT*                               );
    static void   assign_p_infinity          ( double*                                           );
    static void   assign_p_infinity          ( float*                                            );
    static bool   equalszero                 ( double                                            );
    static bool   equalszero                 ( float                                             );
    static double distance                   ( double, double, double, double                    );
    static double distance2                  ( double, double, double, double                    );
    static void multidimmin                  (inputdata*, commoninputdata*,
                                              PS_SIT, double, string, gsl_vector*, gsl_vector*,
                                              void*, double(*)(const gsl_vector*,void*),
                                              double*, PS_SIT, PS_SIT);
    static void shuffle                      (PS_SIT *arr, PS_SIT size);

    /*
     * These are functions that are no longer used in the code
     * but I keep them around just in case
     *
     * static void   quadfit                    ( vector <double>, vector <double>, vector<double>& );
     * static void   intersection               ( double, double, double, double, double,
     * double, double, double, vector<double>&           );
     * static double getweight                  ( PS_SIT, PS_SIT, double**, PS_SIT, bool                     );
     * static double rank                       ( double, double, double, double, double,
     *                                        double, double, double                            );
     * static void   planarinterpolation4pts    ( double, double, double[2][4], double[4]           );
     */

    // templates are defined in pixsrc_operations_templates.cpp

    template <typename T>
    static void   resize                     ( T**, PS_SIT, PS_SIT );
    template <typename T>
    static string tostring                   ( T             );

    // slope, intercept, and xintercept are templates because
    // the regularization strength finding routine might require
    // float versions of these functions
    template <typename T>
    static T      slope                      ( T, T, T, T    );
    template <typename T>
    static T      intercept                  ( T, T, T, T    );
    template <typename T>
    static T      xintercept                 ( T, T, T, T    );

    template <typename T>
    static T      convert_string             ( string        );
    template <typename T>
    static bool   is_finite                   ( T             );
    template <typename T>
    static PS_SIT    round                      ( T             );
    template <typename T>
    static void   assign_n_infinity          ( T*            );
    template <typename T>
    static bool   is_infinite                ( T             );
    template <typename T>
    static bool   is_exactly_p_inf           ( T             );
    template <typename T>
    static T      minarr                     ( T*, PS_SIT, PS_SIT  );
    template <typename T>
    static T      maxarr                     ( T*, PS_SIT, PS_SIT  );
    template <typename T>
    static void   swap                       ( T*, T*);

};

#endif
