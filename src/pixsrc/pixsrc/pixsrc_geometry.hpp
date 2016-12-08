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



#ifndef PIXSRC_GEOMETRY_HPP_
#define PIXSRC_GEOMETRY_HPP_

#define GEOM pixsrc_geometry::

class pixsrc_geometry
{
public:

    static void   getcentroid          ( PS_SIT, double*, double*, double*                  );
    static void   getconvexhull        ( struct triangulateio*, struct triangulateio*,
                                         PS_SIT*, double**, char*                           );
    static PS_SIT    search_triangle      ( struct triangulateio*, PS_SIT*, double, double     );
    static double areatriangle         ( double*                                         );
    static double areapoly             ( double**, PS_SIT                                   );
    static double areapoly             ( double*, PS_SIT                                    );
    static void   clippoly             ( double*, PS_SIT, double*, PS_SIT, double**, PS_SIT*      );
    static double getmaxdist           ( double*, PS_SIT                                    );
    static bool   isintricross         ( double*                                         );
    static bool   isintri              ( double*, double, double                         );
    static bool   isintri              ( double[2][3], double, double                    );
    static bool   isintripadded        ( double*, double, double                         );
    static bool   isintripadded        ( double[2][3], double, double                    );
    static bool   isonlinesegment      ( double, double, double, double, double, double  );
    static double angle                ( double, double, double, double, double, double  );
    static bool   allindiffquadrants   ( double, double, double, double, double, double,
                                         double, double                                  );
    static void   norm_parallelogram_ic( double, double, double, double, double, double,
                                         double, double, bool, bool, bool, bool, PS_SIT*    );
    static bool   isinpolycalc         ( double*, PS_SIT                                    );
    static bool   isinpoly             ( double, double, double*, PS_SIT                    );
    static bool   isinpolypadded       ( double, double, double*, PS_SIT                    );
    static bool   fallsinsrcinputcircle( double, double, double**, double                );
    static void   create_border        ( double*, PS_SIT, double, double**);

    /* These functions aren't used anymore but I keep them around
     *
     * static double areaquad             ( double, double, double, double,
     *                                      double, double, double, double                  );
     * static void   normalparallelogram  ( PS_SIT, PS_SIT, double, double, double, double,
     *                                      double, double, double, double, PS_SIT, bool,
     *                                      bool, bool, bool, PS_SIT*                          );
     * static void   abnormalparallelogram( PS_SIT, PS_SIT, double, double, double, double,
     *                                      double, double, double, double, PS_SIT, bool,
     *                                      bool, bool, bool, PS_SIT*                          );
     * static void   trianglecolumnsknown ( PS_SIT, PS_SIT, double, double, double, double,
     *                                      double, double, double, double, PS_SIT, bool,
     *                                      bool, bool, bool, PS_SIT*                          );
     * static void   trianglerowsknown    ( PS_SIT, PS_SIT, double, double, double, double,
     *                                      double, double, double, double, PS_SIT, bool,
     *                                      bool, bool, bool, PS_SIT*                          );
     */
};

#endif
