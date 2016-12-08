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



#ifndef PIXSRC_INIT_HPP_
#define PIXSRC_INIT_HPP_

#define INIT pixsrc_init::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_tps.hpp"

class pixsrc_init
{

public:

    static void setuplibwcs                ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void setup_one_libwcs           (commoninputdata *cdata_,
					    string fn, ps_WorldCoor **wcs, char *coordsysfinal,
					    double *px, double *py, double *pra, double *pdec, PS_SIT imgy);
    static void settingup                  ( char*, char**, char**, inputdata*, commoninputdata*,
                                             PS_SIT*, char***, char***                         );
    static void readparameters             ( char*, char**, char**, inputdata*, commoninputdata*,
                                             PS_SIT , char** , char**                          );
    static void readimages                 ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void read_one_image             (inputdata *data_, commoninputdata *cdata_, string fn,
					    PS_SIT *wcsinfolength,    double **wcsinfo_,
					    PS_SIT *invwcsinfolength, double **invwcsinfo_,
					    PS_SIT *imgdatalength,    double **imgdata_,
					    PS_SIT *imgx, PS_SIT *imgy, double *rollangle,
					    double *pix2arc, double *arc2pix, PS_SIT *ndp);
    static void setdefaultparameters       ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void readmasks                  ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void readmmmasks                ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void setuppsfs                  ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void finalizeNcounterindications( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void uvdata                     ( char*, char**, char**,
                                             inputdata*, commoninputdata*                   );
    static void printsettings              ( inputdata*, commoninputdata*                   );
    static void parmscreator               ( char***, ps_parms_struct**, PS_SIT*               );
    static void setsource                  ( char*, inputdata*, commoninputdata*            );
    static void setsourcebounds            ( char*, inputdata*, commoninputdata*            );
    static void setsourcestepsizes         ( char*, inputdata*, commoninputdata*            );
    static void clearsource                (        inputdata*, commoninputdata*            );
    static void uvdata_get_rbf_mat         (inputdata *data_, commoninputdata *cdata_,
                                            TPS **imgtps_, PS_SIT *num_ndp_, PS_SIT transmode);

private:

    static void setgpus                    ( char*, inputdata*, commoninputdata*            );
    static void printgpudevices            (                    commoninputdata*            );
    static void readmaskascii              ( PS_SIT* ,PS_SIT** , char**,
                                             inputdata*, commoninputdata*                   );
    static void readmaskpolygon            ( PS_SIT*, PS_SIT*, double***, char**,
                                             inputdata*, commoninputdata*                   );
    static void readmaskcircle             ( PS_SIT*, PS_SIT*, double***, char**,
                                             inputdata*, commoninputdata*                   );
    static PS_SIT  readarr                    ( char*, const char*, const char*, double**,
                                                commoninputdata*                               );
    static void create_mm_border           ( double***, PS_SIT, PS_SIT, PS_SIT, PS_SIT,
                                             inputdata*, commoninputdata*                   );
    static void uvdata_read                ( inputdata*, commoninputdata*                   );
    static void uvdata_img2uv              ( inputdata*, commoninputdata*                   );
    static void uvdata_adjustimgplane      ( inputdata*, commoninputdata*                   );

};


#endif
