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



#include "pixsrc_init.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_printer.hpp"
#include <cstring>
#include <cmath>

void pixsrc_init::parmscreator(char ***parmslist, ps_parms_struct **pstructlist, PS_SIT *size)
{

    // NOTE: Every time a new parameter is added, the pstruct entry for it needs to be created below.
    //       The default value must also be added in the setdefaultparameters function.
    //       Lastly, parameter must be read in in the readparameters function.

    // create structs holding parameter details
    *size = 70;
    ps_parms_struct *pstruct = new ps_parms_struct [*size];

    pstruct[0].pindex = new ps_parms_ind_class();
    for (PS_SIT i=0; i<*size; ++i)
    {
        pstruct[i].scategory[0] = "GENERAL";
        pstruct[i].scategory[1] = "DATA AND GRIDDING";
        pstruct[i].scategory[2] = "SHAPELETS";
        pstruct[i].scategory[3] = "REGULARIZATION";
        pstruct[i].scategory[4] = "ANALYTIC SOURCE";
        pstruct[i].scategory[5] = "UV PLANE MODELLING";
        pstruct[i].scategory[6] = "GPU CODING";
        pstruct[i].scategory[7] = "PENALTIES";
        pstruct[i].scategory[8] = "LENS POTENTIAL PERTURBATION";
        pstruct[i].scategory[9] = "INTERPOLATION ERRORS";

        if (i)
            pstruct[i].pindex = NULL;
    }

    PS_SIT ind = -1;
    // fill in details for parameters
    ++ind;
    pstruct->pindex->ps_parm_ind_debug = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "debug:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Turn debug mode on(1) or off(0).");
    pstruct[ind].notes.push_back  ("Debug mode will print lots of internal variables to file.");
    pstruct[ind].notes.push_back  ("It can be slow and require lots of disk space.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Turn debug mode on.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Turn debug mode off.");

    ++ind;
    pstruct->pindex->ps_parm_ind_fatalwarn = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "fatalwarn:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Turn fatal warnings during initialization on(1) or off(0).");
    pstruct[ind].notes.push_back  ("Turn fatal warnings during initialization on(1) or off(0).");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Turn fatal warnings on.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Turn fatal warnings off.");

    ++ind;
    pstruct->pindex->ps_parm_ind_regorder = ind;
    pstruct[ind].category = ps_parm_cat_regularization;
    pstruct[ind].sname    = "regorder:";
    pstruct[ind].svalue   = "2";
    pstruct[ind].fvalue.push_back (2);
    pstruct[ind].qdescr.push_back ("Select form of regularization.");
    for (PS_SIT i=0; i<4; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" -1");
    pstruct[ind].entries[0][1].push_back ("Penalize large source sizes (only valid when using shapelets).");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Penalize sources with large surface brightnesses.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[2][1].push_back ("Penalize sources with large gradients of surface brightness.");
    pstruct[ind].entries[3][0].push_back (pstruct[ind].sname+" 2");
    pstruct[ind].entries[3][1].push_back ("Penalize sources with large Laplacians of surface brightness.");

    ++ind;
    pstruct->pindex->ps_parm_ind_regstrength = ind;
    pstruct[ind].category = ps_parm_cat_regularization;
    pstruct[ind].sname    = "regstrength:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Guess or fix source regularization strength.");
    pstruct[ind].notes.push_back  ("A value of zero will turn regularization off.");
    pstruct[ind].notes.push_back  ("However, it is numerically more stable to fix the regularization strength to a small value.");
    pstruct[ind].notes.push_back  ("The order of magnitude of the small value will depend on the noise level in the data.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" -\\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Fix regularization strength to \\textless{}value\\textgreater{}, where \\textless{}value\\textgreater{} is positive.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[1][1].push_back ("Provide initial guess, where \\textless{}value\\textgreater{} is positive.");

    ++ind;
    pstruct->pindex->ps_parm_ind_coorsys = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "coorsys:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Set coordinate system and reference point.");
    pstruct[ind].notes.push_back  ("The positions of lenses in the lens model are set using offsets");
    pstruct[ind].notes.push_back  ("in right ascension and declination, measured in arcseconds.");
    pstruct[ind].notes.push_back  ("This parameter defines an origin with respect to which these offsets are applied.");
    pstruct[ind].notes.push_back  ("IMPORTANT: When specifying lens models, angles are measured west of north.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}SYS\\textgreater{} \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} \\textless{}D\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("SYS is the coordinate system you want to use and can be either");
    pstruct[ind].entries[0][1].push_back ("J2000, B1950, ecliptic (2000), galactic (2000).");
    pstruct[ind].entries[0][1].push_back ("For ecliptic or galactic do not type \"(2000)\".");
    pstruct[ind].entries[0][1].push_back ("A is the right ascension (RA) of the reference position (degrees).");
    pstruct[ind].entries[0][1].push_back ("B is the declination (Dec) of the reference position (degrees).");
    pstruct[ind].entries[0][1].push_back ("C is the RA offset of the reference position (arcseconds).");
    pstruct[ind].entries[0][1].push_back ("D is the Dec offset of the reference position (arcseconds).");

    ++ind;
    pstruct->pindex->ps_parm_ind_grid = ind;
    pstruct[ind].category = ps_parm_cat_dataNgridding;
    pstruct[ind].sname    = "grid:";
    pstruct[ind].svalue   = "1 2";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].fvalue.push_back (2);
    pstruct[ind].fvalue.push_back (10);
    pstruct[ind].qdescr.push_back ("Specify grid on which the source will be reconstructed.");
    pstruct[ind].notes.push_back  ("If shapelets are turned on, this option will be ignored.");
    pstruct[ind].notes.push_back  ("The fully adaptive grid is created by ray-tracing a subset of the data pixels to the source plane.");
    pstruct[ind].notes.push_back  ("pixel\\_skip specifies how many data pixels to ignore before ray-tracing a data pixel.");
    pstruct[ind].notes.push_back  ("An image plane pixel is ray-traced if (x*(imgY)+y)\\%skip\\_level==0, where x and y are");
    pstruct[ind].notes.push_back  ("the x- and y-coordinates of the pixel, and imgY is the number of pixels along the y-axis.");
    pstruct[ind].notes.push_back  ("For example, skip\\_level=2 selects every other pixel.");
    pstruct[ind].notes.push_back  ("For irregular Cartesian grids, the zeroth level grid is created using five pixels, spanning the source region.");
    pstruct[ind].notes.push_back  ("Smaller \"nested\" Cartesian grids are created if the magnification in those regions is sufficiently high.");
    pstruct[ind].notes.push_back  ("Because the size of the initial grid is arbritray, the magnification threshold for subgridding is also arbritrary.");
    pstruct[ind].notes.push_back  ("shift\\_level ascribes a magnification to the zeroth level grid. It must be between 1 and infinity, exclusive.");
    pstruct[ind].notes.push_back  ("In general, as shift\\_level increases, the number of source pixels decreases.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1 [pixel\\_skip]");
    pstruct[ind].entries[0][1].push_back ("Select the fully adaptive grid.");
    pstruct[ind].entries[0][1].push_back ("pixel\\_skip controls which image pixels are used to create grid (default "+OPERA tostring(pstruct[ind].fvalue[1])+").");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 2 [shift\\_level]");
    pstruct[ind].entries[1][1].push_back ("Select irregular Cartesian grid.");
    pstruct[ind].entries[1][1].push_back ("shift\\_level sets the reference magnification for the zeroth level grid (default "+OPERA tostring(pstruct[ind].fvalue[2])+").");

    ++ind;
    pstruct->pindex->ps_parm_ind_psf = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "psf:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (101);
    pstruct[ind].qdescr.push_back ("Specify the point spread function (PSF).");
    pstruct[ind].notes.push_back  ("The PSF is computed only once when pixsrc is initialized.");
    pstruct[ind].notes.push_back  ("During convolution, less than 1\\% of the flux is lost.");
    pstruct[ind].notes.push_back  ("Input PSF files must be square with odd dimensions.");
    for (PS_SIT i=0; i<3; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Turn off PSF convolution.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Read PSF from file named \"pixsrc\\_in/\\textless{}imagename\\textgreater{}.psf.fits\".");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 2 \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} [D]");
    pstruct[ind].entries[2][1].push_back ("Set PSF as an elliptical Gaussian.");
    pstruct[ind].entries[2][1].push_back ("\\textless{}A\\textgreater{} is the FWHM along the major axis (arcseconds).");
    pstruct[ind].entries[2][1].push_back ("\\textless{}B\\textgreater{} is the FWHM along the minor axis (arcseconds).");
    pstruct[ind].entries[2][1].push_back ("\\textless{}C\\textgreater{} is the position angle of the major axis (degrees, east of north).");
    pstruct[ind].entries[2][1].push_back ("[D] is the oversampling factor of the PSF in x \\& y directions (default 101).");

    ++ind;
    pstruct->pindex->ps_parm_ind_sisterpix = ind;
    pstruct[ind].category = ps_parm_cat_dataNgridding;
    pstruct[ind].sname    = "sisterpix:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Turn search for sister pixels on(1) or off(0).");
    pstruct[ind].notes.push_back ("The source region is mapped out using image pixels.");
    pstruct[ind].notes.push_back ("Then, image pixels not flagged as bad pixels are ray-traced to the source plane.");
    pstruct[ind].notes.push_back ("If the image pixel lies inside the source region or within one FWHM of the PSF, it is included in the analysis.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Turn on search for sister pixels.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Turn off search for sister pixels.");

    ++ind;
    pstruct->pindex->ps_parm_ind_oversample = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "oversample:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Specify the image pixel oversampling factor, for creating the lensing operator, in x \\& y directions.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set oversampling factor to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_noise = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "noise:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Set the 1-sigma of Gaussian noise in the data.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the 1-sigma of Gaussian noise in the data to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penalty6 = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penalty6:";
    pstruct[ind].svalue   = "0 0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Disabled.");

    ++ind;
    pstruct->pindex->ps_parm_ind_statistic = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "statistic:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Select which statistic will be computed for lens model ranking.");
    pstruct[ind].notes.push_back ("Any additional penalty functions (see penalty? commands) will be added to this statistic.");
    for (PS_SIT i=0; i<3; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Compute Bayesian evidence.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Compute traditional $\\chi^{2}$: $[(data-model)/sigma]^{2}$.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 2");
    pstruct[ind].entries[2][1].push_back ("Compute nothing.");

    ++ind;
    pstruct->pindex->ps_parm_ind_magnification = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "magnification:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Turn magnification calculation on(1) or off(0).");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Do not compute magnification.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Compute magnification.");

    ++ind;
    pstruct->pindex->ps_parm_ind_mag_uncertainty = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "maguncertainty:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Turn magnification uncertainty calculation on(1) or off(0).");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Do not compute magnification uncertainties.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Compute magnification uncertainties.");

    ++ind;
    pstruct->pindex->ps_parm_ind_fullsrccov = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "fullsrccov:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Use full source covariance matrix for magnification errors yes(1) or no(0).");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Use full source covariance matrix for calculating mag errors.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Use only diagonals of covariance matrix for calculating mag errors. Not recommended but might be necessary to avoid large matrices.");

    ++ind;
    pstruct->pindex->ps_parm_ind_magsamples = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "magsamples:";
    pstruct[ind].svalue   = "1000";
    pstruct[ind].fvalue.push_back (1000);
    pstruct[ind].qdescr.push_back ("How many realizations of the source brightness distribution to sample for calculating mag errors.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1000");
    pstruct[ind].entries[0][1].push_back ("Sample 1000 realizations of the source brightness distribution to sample for calculating mag errors.");

    ++ind;
    pstruct->pindex->ps_parm_ind_verbosity = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "verbosity:";
    pstruct[ind].svalue   = "3";
    pstruct[ind].fvalue.push_back (3);
    pstruct[ind].qdescr.push_back ("Control how often pixsrc prints to the screen.");
    pstruct[ind].notes.push_back ("Currently, \"verbosity: 2\" is identical to \"verbosity: 3\" and prints as much as possible.");
    for (PS_SIT i=0; i<3; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Print very little.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 2");
    pstruct[ind].entries[1][1].push_back ("Print only significant progress.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 3");
    pstruct[ind].entries[2][1].push_back ("Print as much as possible.");

    ++ind;
    pstruct->pindex->ps_parm_ind_images = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "images:";
    pstruct[ind].svalue   = "3";
    pstruct[ind].fvalue.push_back (3);
    pstruct[ind].qdescr.push_back ("Select which data and/or images will be written to file.");
    pstruct[ind].notes.push_back ("Files to be written include (but not limited to): data, models, residuals, source reconstructions.");
    for (PS_SIT i=0; i<4; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Write nothing to file.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Write text files only.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 2");
    pstruct[ind].entries[2][1].push_back ("Write FITS files only.");
    pstruct[ind].entries[3][0].push_back (pstruct[ind].sname+" 2");
    pstruct[ind].entries[3][1].push_back ("Write text and FITS files.");

    ++ind;
    pstruct->pindex->ps_parm_ind_varystrength = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "varystrength:";
    pstruct[ind].svalue   = "0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Trace the lens model statistic and/or magnification as a function of regularization strength.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Trace the statistic and/or magnification from \\textless{}A\\textgreater{} to \\textless{}B\\textgreater{}, taking \\textless{}C\\textgreater{} logarithmic steps.");

    ++ind;
    pstruct->pindex->ps_parm_ind_threads = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "threads:";
    pstruct[ind].svalue   = "8";
    pstruct[ind].fvalue.push_back (2);
    pstruct[ind].qdescr.push_back ("Set the number of CPU threads per input data file.");
    pstruct[ind].notes.push_back ("This parameter specifies the maximum number of CPU threads per input data file.");
    pstruct[ind].notes.push_back ("It also specifies the maximum number of input data files to simultaneoulsy analyze.");
    pstruct[ind].notes.push_back ("Thus, the maximum possible number of CPU threads is the square of this parameter.");
    pstruct[ind].notes.push_back ("Some optimized, multithreaded basic linear algebra subroutines (BLAS) libraries will use some pre-configured number of threads, and pixsrc cannot control this.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the maximum number of CPU threads per data file to A.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penalty1 = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penalty1:";
    pstruct[ind].svalue   = "0 0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Penalize the brightness-weighted source size (square arcseconds) of the reconstructed source.");
    pstruct[ind].qdescr.push_back ("The weighting function is the absolute value of the surface brightness.");
    pstruct[ind].notes.push_back ("A number of penalty functions (priors) can be placed on the source reconstruction.");
    pstruct[ind].notes.push_back ("In mode 1 (see below), the penalty is an additional $\\chi^{2}$ term based on the difference between what is expected from the prior and what is calculated.");
    pstruct[ind].notes.push_back ("This expectation is referred to as the best estimate.");
    pstruct[ind].notes.push_back ("This additional $\\chi^{2}$ is added to the primary statistic being computed ($\\chi^{2}$, evidence, or neither).");
    pstruct[ind].notes.push_back ("In mode 2, hard lower and upper limits on the prior are specified.");
    pstruct[ind].notes.push_back ("If the calculated value is below the lower limit or above the upper limit, then pixsrc is terminated and the $\\chi^{2}$ term is computed and returned, without computing the primary statistic.");
    pstruct[ind].notes.push_back ("Otherwise, the primary statistic is computed as usual.");
    pstruct[ind].notes.push_back ("There is also a hybrid mode (mode 3) in which mode 1 operates if the calculated value is within the lower/upper limits and mode 2 operates otherwise.");
    pstruct[ind].notes.push_back ("By default, a \"penalty?\" command will apply to all input data files and to all \"mmimages\" region files.");
    pstruct[ind].notes.push_back ("To change this behavior, use the \"dataname\" and \"penaltyname\" commands.");
    pstruct[ind].notes.push_back ("For all penalty functions, the following arguments apply:");
    pstruct[ind].notes.push_back ("arg \\textless{}A\\textgreater{}: mode: 0(off), 1(best estimate), 2(lower/upper bounds), 3(1 \\& 2)");
    pstruct[ind].notes.push_back ("arg \\textless{}B\\textgreater{}: best estimate");
    pstruct[ind].notes.push_back ("arg \\textless{}C\\textgreater{}: uncertainty in best estimate");
    pstruct[ind].notes.push_back ("arg \\textless{}D\\textgreater{}: lower limit");
    pstruct[ind].notes.push_back ("arg \\textless{}E\\textgreater{}: upper limit");
    pstruct[ind].notes.push_back ("arg \\textless{}F\\textgreater{}: uncertainty outside lower/upper");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} \\textless{}D\\textgreater{} \\textless{}E\\textgreater{} \\textless{}F\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("See notes on \"penalty1\" for usage.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penalty2 = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penalty2:";
    pstruct[ind].svalue   = "0 0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Penalize the magnification.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} \\textless{}D\\textgreater{} \\textless{}E\\textgreater{} \\textless{}F\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("See notes on \"penalty1\" for usage.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penalty3 = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penalty3:";
    pstruct[ind].svalue   = "0 0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Penalize the axis ratio (major/minor) of the \"mmimages\" (computed before source reconstruction).");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} \\textless{}D\\textgreater{} \\textless{}E\\textgreater{} \\textless{}F\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("See notes on \"penalty1\" for usage.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penalty4 = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penalty4:";
    pstruct[ind].svalue   = "0 0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Penalize the fractional overlap of the \"mmimages\" in the source plane (computed before source reconstruction).");
    pstruct[ind].qdescr.push_back ("If there is no overlap between two images, then a negative value is returned, whose magnitude is the minimum distance between the two (ray-traced) images.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} \\textless{}D\\textgreater{} \\textless{}E\\textgreater{} \\textless{}F\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("See notes on \"penalty1\" for usage.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penalty5 = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penalty5:";
    pstruct[ind].svalue   = "0 0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Penalize the size (square arcseconds) of the union of all ray-traced \"mmimages\" in the source plane (computed before source reconstruction).");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}A\\textgreater{} \\textless{}B\\textgreater{} \\textless{}C\\textgreater{} \\textless{}D\\textgreater{} \\textless{}E\\textgreater{} \\textless{}F\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("See notes on \"penalty1\" for usage.");

    ++ind;
    pstruct->pindex->ps_parm_ind_details = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "details:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Write detailed log file for each pixsrc run: on(1), off(0).");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Write log file.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Do not write log file.");

    ++ind;
    pstruct->pindex->ps_parm_ind_fillbadpix = ind;
    pstruct[ind].category = ps_parm_cat_dataNgridding;
    pstruct[ind].sname    = "fillbadpix:";
    pstruct[ind].svalue   = "0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Change bad pixels to good pixels and change the pixel value.");
    for (PS_SIT i=0; i<3; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0 0");
    pstruct[ind].entries[0][1].push_back ("Do not change bad pixels.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1 noise");
    pstruct[ind].entries[1][1].push_back ("Change bad pixels and replace pixel values with noise.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 1 \\textless{}value\\textgreater{}");
    pstruct[ind].entries[2][1].push_back ("Change bad pixels and replace pixel values with \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_regaccuracy = ind;
    pstruct[ind].category = ps_parm_cat_regularization;
    pstruct[ind].sname    = "regaccuracy:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Compute higher accuracy regularization: on(1), off(0).");
    pstruct[ind].notes.push_back ("Higher accuracy will require more computational effort.");
    pstruct[ind].notes.push_back ("This option is only available for some combinations of grids and regularization.");
    pstruct[ind].notes.push_back ("This option has no effect on shapelets, which are analytic and exact.");
    pstruct[ind].notes.push_back ("For derivative-based regularization schemes, higher accuracy means using the divergence theorem to compute derivatives.");
    pstruct[ind].notes.push_back ("Higher accuracy is not yet available for curvature (Laplacian-based) regularization.");
    pstruct[ind].notes.push_back ("For analytic source regularization (ASR), higher accuracy means that all pixels connected to a given pixel are used in regularizing it, as opposed to only using the nearest four pixels.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Compute higher accuracy regularization.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Do no compute higher accuracy regularization.");

    ++ind;
    pstruct->pindex->ps_parm_ind_minangle = ind;
    pstruct[ind].category = ps_parm_cat_dataNgridding;
    pstruct[ind].sname    = "minangle:";
    pstruct[ind].svalue   = "0.1"; // this is fed into triangle code
    pstruct[ind].fvalue.push_back (0.1);
    pstruct[ind].qdescr.push_back ("Set the minimum angle in any triangle in the pseudo-Delaunay triangulation of the source plane grid.");
    pstruct[ind].notes.push_back ("Small angles can lead to problems in calculating derivatives and surface brightness.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the minimum angle to \\textless{}value\\textgreater{}, given in degrees.");

    ++ind;
    pstruct->pindex->ps_parm_ind_reg = ind;
    pstruct[ind].category = ps_parm_cat_regularization;
    pstruct[ind].sname    = "reg:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Turn regularization on(1) or off(0).");
    pstruct[ind].notes.push_back ("Turning regularization off is necessary to use strictly analytic source(s) to model the data.");
    pstruct[ind].notes.push_back ("Otherwise, to turn regularization off, it is numerically more stable to fix the regularization strength to a small value.");
    pstruct[ind].notes.push_back ("The order of magnitude of the small value will depend on the noise level in the data.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Regularize the source.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Do not regularize the source.");

    ++ind;
    pstruct->pindex->ps_parm_ind_rmpix = ind;
    pstruct[ind].category = ps_parm_cat_dataNgridding;
    pstruct[ind].sname    = "rmpix:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Remove source pixels not constrained by lensing operator: on(1), off(0).");
    pstruct[ind].notes.push_back ("Not removing unconstrained pixels may lead to numerical instability and/or more computational effort.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Remove unconstrained pixels.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Do not remove unconstrained pixels.");

    ++ind;
    pstruct->pindex->ps_parm_ind_src = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "src:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Fit one or more analytic sources to the data.");
    pstruct[ind].notes.push_back ("If regularization is turned on, analytic source regularization (ASR) will be used.");
    pstruct[ind].notes.push_back ("If regularization is turned off, the analytic source(s) will be used to model the data.");
    pstruct[ind].notes.push_back ("Three source profile types are implemented: \"none\", \"sersic\", and \"vector\\_\\textless{}vecfilename\\textgreater{}\".");
    pstruct[ind].notes.push_back ("\"none\" denotes a source with zero surface brightness everywhere.");
    pstruct[ind].notes.push_back ("\"sersic\" denotes a S\\'{e}rsic light profile with spherically symmetric light profile, $I(r)$, given by $I(r) = I_{0}\\exp(-(r/R)^{1/n})$.");
    pstruct[ind].notes.push_back ("If there are N sources, the first N lines of the file should each contain the source profile type, followed by 8 source parameters.");
    pstruct[ind].notes.push_back ("The next N lines should contain 8 parameter flags which denote whether the corresponding source parameter will(1) or will not(0) be optimized.");
    pstruct[ind].notes.push_back ("For the \"sersic\" profile type, the 8 source parameters correspond to: $I_{0}$, ra, dec, e, PA, R, Z, k; where $I_{0}$ is the normalization, ra and dec denote the source center in R.A. and Dec. offsets, e and PA are the ellipticity and position angle (east of north) of the major axis, R is the scale radius, Z is an ignored parameter, and n is the S\\'{e}rsic index.");
    pstruct[ind].notes.push_back ("\"vector\\_\\textless{}vecfilename\\textgreater{}\" allows an arbritrary surface brightness distribution to be lensed.");
    pstruct[ind].notes.push_back ("The file \\textless{}vecfilename\\textgreater{} contains 3 columns: right ascension, declination, surface brightness, where the positions are given in arcseconds offsets, using the coordinate system defined using \"coorsys:\".");
    pstruct[ind].notes.push_back ("Like the \"sersic\" profile, there can be multiple sources.");
    pstruct[ind].notes.push_back ("The syntax is: vector\\_\\textless{}vecfilename\\textgreater{} $I_{0}$ ra dec m11 m12 m21 m22 Z.");
    pstruct[ind].notes.push_back ("$I_{0}$ is an overall scaling of the surface brightness. ra and dec are position offsets in arcseconds. Z is ignored. m?? denote elements of a linear transformation of the coordinates, so that the new coordinates are given by");
    pstruct[ind].notes.push_back ("\\begin{equation}");
    pstruct[ind].notes.push_back ("\\bigg({ra' \\atop dec'}\\bigg) = \\bigg({m11~m12 \\atop m21~m22}\\bigg) \\bigg({ra \\atop dec}\\bigg)");
    pstruct[ind].notes.push_back ("\\end{equation}");
    pstruct[ind].notes.push_back ("This allows for rotations, stretches, reflections, etc. of the coordinates to be optimized.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Read analytic source parameters from \"pixsrc\\_in/\\textless{}filename\\textgreater{}\".");

    ++ind;
    pstruct->pindex->ps_parm_ind_srcbound = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "srcbound:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Specify bounds for analytic source parameters.");
    pstruct[ind].notes.push_back ("A bound consists of a combination of linear terms and a lower and upper bound for the combination..");
    pstruct[ind].notes.push_back ("If there are N bounds, there are N lines. For a particular bound, if there are M terms, there are 3*M+2 entries per line: isrc iparm ai  jsrc jparm aj  lower\\_bound upper\\_bound, where isrc>=1 and 1<=iparm<=8.");
    pstruct[ind].notes.push_back ("The bound is given by: lower\\_bound $<=$ ai*pi + aj*pj + ... + aN*pN $<=$ upper\\_bound, where pj is parameter number jparm of source number jsrc.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Read analytic source bounds from \"pixsrc\\_in/\\textless{}filename\\textgreater{}\".");

    ++ind;
    pstruct->pindex->ps_parm_ind_hacksrc = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "hacksrc:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Specify values of source parameters (as stored in code/C++ array).");
    pstruct[ind].notes.push_back ("Specify values of source parameters (as stored in code/C++ array). Not applicable to analytic source.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0 [filename]");
    pstruct[ind].entries[0][1].push_back ("Do not hack the source vector in the code.");
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1 \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Do not solve for most probable source. Instead, use source parameters found in \\textless{}filename\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_srcstepsize = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "srcstepsize:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Specify initial stepsizes for analytic source parameters.");
    pstruct[ind].notes.push_back ("If there are N sources, there are N lines containing 8 entries each, and each entry specifies the stepsize for the corresponding analytic source parameter.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Read analytic source stepsizes from \"pixsrc\\_in/\\textless{}filename\\textgreater{}\".");

    ++ind;
    pstruct->pindex->ps_parm_ind_srcrestart = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "srcrestart:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Set the number of optimizations of the analytic source parameters to perform.");
    pstruct[ind].notes.push_back ("Selecting the best of several optimizations helps ensure a global minimum is found.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Performs \\textless{}value\\textgreater{} optimizations of the analytic source parameters.");

    ++ind;
    pstruct->pindex->ps_parm_ind_srcexact = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "srcexact:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Set method for lensing analytic source: exact(1), approx(0).");
    pstruct[ind].notes.push_back ("Exact method does not use lensing matrix but ray-traces to source plane directly. Approximate methos uses lensing matrix and source grid, which can be quicker.");
    pstruct[ind].notes.push_back ("If using exact method the source will not appear in output source-plane images, and analytic source regularization will be disabled.");
    pstruct[ind].notes.push_back ("If using approximate method, oversampling can still be used.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Do not use lensing matrix and source grid.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Use lensing matrix and source grid.");

    ++ind;
    pstruct->pindex->ps_parm_ind_interperr = ind;
    pstruct[ind].category = ps_parm_cat_interperr;
    pstruct[ind].sname    = "interperr:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Compute interpolation errors from an analytic function: on(>0), off(0).");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Do not compute interpolation errors.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[1][1].push_back ("Compute interpolation errors, oversampling by a factor of \\textless{}value\\textgreater{} in x \\& y dimensions.");

    ++ind;
    pstruct->pindex->ps_parm_ind_irscheme = ind;
    pstruct[ind].category = ps_parm_cat_interperr;
    pstruct[ind].sname    = "irscheme:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Specify what to do with interpolation errors.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Add interpolation errors to noise covariance matrix.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("\"Correct\" the lensed and PSF-smeared image-plane surface brightness.");

    ++ind;
    pstruct->pindex->ps_parm_ind_srcftol = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "srcftol:";
    pstruct[ind].svalue   = "1e-3";
    pstruct[ind].fvalue.push_back (1e-3);
    pstruct[ind].qdescr.push_back ("Set convergence criterion for analytic source optimization.");
    pstruct[ind].notes.push_back ("This parameter is passed to the GSL multidimensional minimizer.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the convergence criterion for analytic source parameters to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_cuda = ind;
    pstruct[ind].category = ps_parm_cat_gpu;
    pstruct[ind].sname    = "cuda:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Turn GPU computing on/off and/or benchmark GPU(s).");
    for (PS_SIT i=0; i<4; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" help");
    pstruct[ind].entries[0][1].push_back ("List the available GPUs detected by pixsrc/CUDA and their names.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" benchmark");
    pstruct[ind].entries[1][1].push_back ("Test the GPUs and their speeds.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" ignore \\textless{}gpulist\\textgreater{}");
    pstruct[ind].entries[2][1].push_back ("Ignore the GPUs in the white-space separated list, \\textless{}gpulist\\textgreater{}.");
    pstruct[ind].entries[3][0].push_back (pstruct[ind].sname+" \\textless{}gpuname\\textgreater{}");
    pstruct[ind].entries[3][1].push_back ("Use the GPU called \\textless{}gpuname\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_optfinder = ind;
    pstruct[ind].category = ps_parm_cat_regularization;
    pstruct[ind].sname    = "optfinder:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Select how the regularization strength is found.");
    pstruct[ind].notes.push_back ("Maximizing the evidence can be faster; however, the root finder is more robust.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Maximize the evidence to find optimal regularization strength.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Use a derivative-based root finder.");

    ++ind;
    pstruct->pindex->ps_parm_ind_nonparamlens = ind;
    pstruct[ind].category = ps_parm_cat_potentialperturbation;
    pstruct[ind].sname    = "nonparamlens:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Use non-parametric lens potential perturbations.");
    pstruct[ind].notes.push_back ("FILL IN THIS SECTION AFTER INPUT FORMAT IS CLEANED UP.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Turn off potential perturbations.");

    ++ind;
    pstruct->pindex->ps_parm_ind_reglens = ind;
    pstruct[ind].category = ps_parm_cat_potentialperturbation;
    pstruct[ind].sname    = "reglens:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Set the regularization strength for the penalty on lens potential perturbations.");
    pstruct[ind].notes.push_back ("The lens potential penalty function is the magnitude of the gradient of the convergence.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the lens potential perturbation regularization strength to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_nplstepsize = ind;
    pstruct[ind].category = ps_parm_cat_potentialperturbation;
    pstruct[ind].sname    = "nplstepsize:";
    pstruct[ind].svalue   = "2";
    pstruct[ind].fvalue.push_back (2);
    pstruct[ind].qdescr.push_back ("Set the stepsize for lens potential perturbations.");
    pstruct[ind].notes.push_back ("The stepsize is passed to the GSL multidimensional minimizer.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the lens potential perturbation stepsize to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_srcctrguess = ind;
    pstruct[ind].category = ps_parm_cat_analyticsrc;
    pstruct[ind].sname    = "srcctrguess:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Estimate center of source for each lens model: on(1), off(0).");
    pstruct[ind].notes.push_back ("For each lens model, the average of the brightness-weighted positions of the ray-traced image-plane pixels is chosen as an initial guess for the source center.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Estimate center of source for each lens model.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Do not estimate center of source (use user-specified center).");

    ++ind;
    pstruct->pindex->ps_parm_ind_nplftol = ind;
    pstruct[ind].category = ps_parm_cat_potentialperturbation;
    pstruct[ind].sname    = "nplftol:";
    pstruct[ind].svalue   = "1e-3";
    pstruct[ind].fvalue.push_back (1e-3);
    pstruct[ind].qdescr.push_back ("Set convergence criterion for lens potential optimization.");
    pstruct[ind].notes.push_back ("This parameter is passed to the GSL multidimensional minimizer.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the convergence criterion for lens potential perturbation to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_npltpsreg = ind;
    pstruct[ind].category = ps_parm_cat_potentialperturbation;
    pstruct[ind].sname    = "npltpsreg:";
    pstruct[ind].svalue   = "0 3";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (3);
    pstruct[ind].qdescr.push_back ("Set TPS regularization strength and stepsize.");
    pstruct[ind].notes.push_back ("Set TPS regularization strength and stepsize.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}strength\\textgreater{} \\textless{}ss\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the regularization strength to $10^{<strength>}$ and the stepsize (of the log regularization strength) to \\textless{}ss\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_shapelets = ind;
    pstruct[ind].category = ps_parm_cat_shapelets;
    pstruct[ind].sname    = "shapelets:";
    pstruct[ind].svalue   = "0 10 10";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (10);
    pstruct[ind].fvalue.push_back (10);
    pstruct[ind].qdescr.push_back ("Use shapelets (instead of using a source grid) to parameterize source.");
    pstruct[ind].notes.push_back ("Although grid parameters are ignored, a sparse grid is constructed to compute the convex hull.");
    pstruct[ind].notes.push_back ("Source regularization can still be used.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0 [nx] [ny]");
    pstruct[ind].entries[0][1].push_back ("Do not use shapelets.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1 \\textless{}nx\\textgreater{} \\textless{}ny\\textgreater");
    pstruct[ind].entries[1][1].push_back ("Use \\textless{}nx\\textgreater{} $\\times$ \\textless{}ny\\textgreater shapelets to reconstruct the source.");

    ++ind;
    pstruct->pindex->ps_parm_ind_shapeparms = ind;
    pstruct[ind].category = ps_parm_cat_shapelets;
    pstruct[ind].sname    = "shapeparms:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Fix shapelets center and/or scale.");
    for (PS_SIT i=0; i<4; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0 [ra] [dec] [sc]");
    pstruct[ind].entries[0][1].push_back ("Automatically estimate source center and scale.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1 \\textless{}ra\\textgreater{} \\textless{}dec\\textgreater{} [sc]");
    pstruct[ind].entries[1][1].push_back ("Fix shapelets center to (\\textless{}ra\\textgreater{},\\textless{}dec\\textgreater{}), given as R.A. and Dec. offsets.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 2 [ra] [dec] \\textless{}sc\\textgreater{}");
    pstruct[ind].entries[2][1].push_back ("Fix shapelets scale to \\textless{}sc\\textgreater{}, given in arcseconds.");
    pstruct[ind].entries[3][0].push_back (pstruct[ind].sname+" 3 \\textless{}ra\\textgreater{} \\textless{}dec\\textgreater{} \\textless{}sc\\textgreater{}");
    pstruct[ind].entries[3][1].push_back ("Fix shapelets center to (\\textless{}ra\\textgreater{},\\textless{}dec\\textgreater{}), given as R.A. and Dec. offsets, and shapelets scale to \\textless{}sc\\textgreater{}, given in arcseconds.");

    ++ind;
    pstruct->pindex->ps_parm_ind_shapepixsplit = ind;
    pstruct[ind].category = ps_parm_cat_shapelets;
    pstruct[ind].sname    = "shapepixsplit:";
    pstruct[ind].svalue   = "1 8";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].fvalue.push_back (8);
    pstruct[ind].qdescr.push_back ("Set minimum and maximum pixel splitting in image plane.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}min\\textgreater{} \\textless{}max\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set minimum and maximum pixel splitting in image plane to \\textless{}min\\textgreater{} and \\textless{}max\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_penaltyname = ind;
    pstruct[ind].category = ps_parm_cat_penalty;
    pstruct[ind].sname    = "penaltyname:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Set the name of the \"mmimages\" region file to which subsequent \"penalty?\" calls will apply.");
    pstruct[ind].notes.push_back ("The argument to this command can be the entire filename (excluding directory prefixes) or just the string between \"mmimages.\" and \".reg\".");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}name\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the name of the \"mmimages\" region file to which subsequent \"penalty?\" calls apply to \\textless{}name\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_dataname = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "dataname:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].qdescr.push_back ("Set the name of the input data file to which subsequent commands will apply.");
    pstruct[ind].notes.push_back ("The argument to this command can be the entire filename (excluding directory prefixes) or the filename with the \".fits\" extension omitted.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}name\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the name of the input data file to which subsequent commands apply to \\textless{}name\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_rottrans = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "rottrans:";
    pstruct[ind].svalue   = "0 0 0 0 0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Rotate and translate input data, relative to coordinate system specified by \"coorsys\".");
    pstruct[ind].notes.push_back ("This option is useful for registering multiple data files.");
    pstruct[ind].notes.push_back ("For example, the following two sets of commands would produce the same results:");
    pstruct[ind].notes.push_back ("\"coorsys: J2000 160 25 0 0\" and \"rottrans: 0 0 0  0  0\"");
    pstruct[ind].notes.push_back ("\"coorsys: J2000 160 25 1 1\" and \"rottrans: 0 0 0 -1 -1\"");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}rot\\_ra\\textgreater{} \\textless{}rot\\_dec\\textgreater{} \\textless{}theta\\textgreater{} \\textless{}ra\\textgreater{} \\textless{}dec\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Rotate the data by an angle \\textless{}theta\\textgreater{} (given in degrees) about the point (\\textless{}rot\\_ra\\textgreater{},\\textless{}rot\\_dec\\textgreater{}) and translate by \\textless{}ra\\textgreater{},\\textless{}dec\\textgreater{} in R.A. and Dec. (All positions are given in arcsecond offsets.).");

    ++ind;
    pstruct->pindex->ps_parm_ind_noisemap = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "noisemap:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Compute noise and S/N ratio maps in the source plane: on(1), off(0).");
    pstruct[ind].notes.push_back ("This can be a long and memory intensive process.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Compute noise and S/N ratio maps in source plane.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Do not compute noise and S/N ratio maps in source plane.");

    ++ind;
    pstruct->pindex->ps_parm_ind_raydirection = ind;
    pstruct[ind].category = ps_parm_cat_dataNgridding;
    pstruct[ind].sname    = "raydirection:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Constrain image pixels by lensing source pixels forward (0) or by ray-tracing image pixels backwards (1).");
    pstruct[ind].notes.push_back ("As coded, constraining image pixels by lensing source pixels forward does not work very well.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Constrain image pixels by ray-tracing image pixels backwards.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Constrain image pixels by lensing source pixels forward.");

    ++ind;
    pstruct->pindex->ps_parm_ind_regftol = ind;
    pstruct[ind].category = ps_parm_cat_regularization;
    pstruct[ind].sname    = "regftol:";
    pstruct[ind].svalue   = "1e-3";
    pstruct[ind].fvalue.push_back (1e-3);
    pstruct[ind].qdescr.push_back ("Set convergence criterion for source regularization strength.");
    pstruct[ind].notes.push_back ("This parameter is passed to the GSL multidimensional minimizer if evidence is maximized.");
    pstruct[ind].notes.push_back ("This parameter is the fractional change needed for convergence if using the derivative-based solver.");
    pstruct[ind].notes.push_back ("See \"optfinder\" for more details about finding regularization strength.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the convergence criterion for source regularization strength to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_outprecision = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "outprecision:";
    pstruct[ind].svalue   = "7";
    pstruct[ind].fvalue.push_back (7);
    pstruct[ind].qdescr.push_back ("Set the decimal precision for writing floating-point numbers to file.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the decimal precision for writing floating-point numbers to file to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_fullmag = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "fullmag:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Select whether model surface brightness across all image pixels or only those masked as good pixels will be used to compute magnification.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[0][1].push_back ("Use model surface brightness across all image pixels to compute magnification.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[1][1].push_back ("Use model surface brightness across only masked image pixels to compute magnification.");

    ++ind;
    pstruct->pindex->ps_parm_ind_pixelscalesrc = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "pixelscalesrc:";
    pstruct[ind].svalue   = "1/4 img scale";
    pstruct[ind].fvalue.push_back (4); // 1/4 of image pixel scale
    pstruct[ind].qdescr.push_back ("Set the source pixel scale used to write FITS files.");
    pstruct[ind].notes.push_back ("The default value is 1/4 the image pixel scale.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the source pixel scale used to write FITS files to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvdata = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvdata:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Specify a file containing visibility measurements for uv-plane modelling.");
    pstruct[ind].notes.push_back ("The file containing visibility measurements have 4 columns: u, v, real\\_vis, imag\\_vis.");
    pstruct[ind].notes.push_back ("The file containing integrated visibility measurements is produced by pixsrc. It has 6 columns: u, v, real\\_integral, imag\\_integral, real\\_noise, imag\\_noise.");
    for (PS_SIT i=0; i<3; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0 [filename]");
    pstruct[ind].entries[0][1].push_back ("Do not do uv-plane modelling.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1 \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[1][1].push_back ("Read visibility measurements from pixsrc\\_in/\\textless{}filename\\textgreater{}.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" 2 \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[2][1].push_back ("Read integrated visibility measurements from pixsrc\\_in/\\textless{}filename\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvmodelpos = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvmodelpos:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Specify a file containing uv positions where model visibilities should be calculated.");
    pstruct[ind].notes.push_back ("This should be a binary file, which can be loaded into a C++ double array. First number is the number of uv points.");
    pstruct[ind].notes.push_back ("The remaining numbers should be all u points, followed by all v points.");
    pstruct[ind].notes.push_back ("This is valid if uvdata is turned on with option 1. This is a slow, non-optimized calculation and is only really intended to be performed on the final lens and source model.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0 [filename]");
    pstruct[ind].entries[0][1].push_back ("Do not produce model visibilities.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1 \\textless{}filename\\textgreater{}");
    pstruct[ind].entries[1][1].push_back ("Read uv positions from pixsrc\\_in/\\textless{}filename\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvmatrixsize = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvmatrixsize:";
    pstruct[ind].svalue   = "100000 2500";
    pstruct[ind].fvalue.push_back (100000);
    pstruct[ind].fvalue.push_back (2500);
    pstruct[ind].qdescr.push_back ("Specify dimensions of submatrices in \"divide and conquer\" multiplication.");
    pstruct[ind].notes.push_back ("Two matrices will be allocated in memory. Each contains 2$\\times$\\textless{}dim1\\textgreater{}$\\times$\\textless{}dim2\\textgreater{} double elements.");
    pstruct[ind].notes.push_back ("The default setting requires approximately 8Gb of memory, unless the number of uv data points is less than \\textless{}dim1\\textgreater{} and/or the number of image pixels is less than \\textless{}dim2\\textgreater{}.");
    pstruct[ind].notes.push_back ("It may be more efficient to keep \\textless{}dim1\\textgreater{} greater than \\textless{}dim2\\textgreater{}.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}dim1\\textgreater{} \\textless{}dim2\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Process \\textless{}dim1\\textgreater{} complex uv data points and \\textless{}dim2\\textgreater{} image pixels at a time.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvpbeam = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvpbeam:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Turn primary beam application on or off.");
    pstruct[ind].notes.push_back ("If primary beam is being applied, the code will look for \"pixsrc\\_in/\\textless{}basename\\textgreater{}\\_\\textless{}imagename\\textgreater{}.pbeam.dat\".");
    pstruct[ind].notes.push_back ("If the file has only one non-blank line containing the alphanumeric string \\textless{}descr\\textgreater{}, then the primary beam found in \"pixsrc\\_in/\\textless{}imagename\\textgreater{}.pbeam.\\textless{}descr\\textgreater{}.fits\" will be applied to all uv data points.");
    pstruct[ind].notes.push_back ("If the file is empty or does not exist, the primary beam found in \"pixsrc\\_in/\\textless{}imagename\\textgreater{}.pbeam.fits\" will be applied to all uv data points.");
    pstruct[ind].notes.push_back ("The primary beam can, however, vary between baselines. To achieve this, the file should contain as many lines as there are complex uv data points.");
    pstruct[ind].notes.push_back ("Each line should, as before, contain one alphanumeric string \\textless{}descr\\textgreater{}, which corresponds to the FITS file, \"pixsrc\\_in/\\textless{}imagename\\textgreater{}.pbeam.\\textless{}descr\\textgreater{}.fits\".");
    pstruct[ind].notes.push_back ("NOTE: Baseline-specific primary beams are not yet implemented.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Do not apply primary beam.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Apply primary beam.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvnewpixelscale = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvpixelscale:";
    pstruct[ind].svalue   = "-0.5";
    pstruct[ind].fvalue.push_back (-0.5);
    pstruct[ind].qdescr.push_back ("Set the image pixel scale for uv modelling.");
    pstruct[ind].notes.push_back ("");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("If \\textless{}value\\textgreater{} is negative, the pixel scale is set to 1/$b_{max}\\times206265\\times$(-\\textless{}value\\textgreater{}), where $b_{max}$ is the maximum baseline. Otherwise, the pixel scale is set to \\textless{}value\\textgreater{}.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvrbf = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvrbf:";
    pstruct[ind].svalue   = "0 1";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Set RBF type and scale length.");
    pstruct[ind].notes.push_back ("");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}type\\textgreater{} \\textless{}scale\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set the RBF type to Gaussian (type=0). Set the scale length to scale.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvpadzeros = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvpadzeros:";
    pstruct[ind].svalue   = "1";
    pstruct[ind].fvalue.push_back (1);
    pstruct[ind].qdescr.push_back ("Assume pixels surrounding the max contain zero surface brightness and use them to constrain RBF weights.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Do not use surrounding pixels.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Use surrounding pixels.");

    ++ind;
    pstruct->pindex->ps_parm_ind_lowmem = ind;
    pstruct[ind].category = ps_parm_cat_general;
    pstruct[ind].sname    = "lowmem:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Disabled -- Use algorithms tuned for low memory systems: on(1), off(0).");
    pstruct[ind].notes.push_back ("Disabled.");
    for (PS_SIT i=0; i<2; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Use normal algorithms.");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" 1");
    pstruct[ind].entries[1][1].push_back ("Use low memory algorithms.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvtaper = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvtaper:";
    pstruct[ind].svalue   = "0";
    pstruct[ind].fvalue.push_back (0);
    pstruct[ind].qdescr.push_back ("Apply uv taper; i.e., increase noise at longer baselines.");
    for (PS_SIT i=0; i<3; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" 0");
    pstruct[ind].entries[0][1].push_back ("Do not apply taper");
    pstruct[ind].entries[1][0].push_back (pstruct[ind].sname+" -1");
    pstruct[ind].entries[1][1].push_back ("Automatically determine taper scale.");
    pstruct[ind].entries[2][0].push_back (pstruct[ind].sname+" \\textless{}value\\textgreater{}");
    pstruct[ind].entries[2][1].push_back ("Apply taper: exp[(b/\\textless{}value\\textgreater{})${}^2$], where b is baseline lensgth.");

    ++ind;
    pstruct->pindex->ps_parm_ind_uvcutoff = ind;
    pstruct[ind].category = ps_parm_cat_uvdata;
    pstruct[ind].sname    = "uvcutoff:";
    pstruct[ind].svalue   = "none";
    pstruct[ind].fvalue.push_back (-1); // -1 is a flag
    pstruct[ind].fvalue.push_back (-1); // -1 is a flag
    pstruct[ind].qdescr.push_back ("Set lower and upper cutoffs for uv-plane data.");
    for (PS_SIT i=0; i<1; ++i)
        pstruct[ind].entries.push_back (vector< vector <string> > (2, vector<string>()));
    pstruct[ind].entries[0][0].push_back (pstruct[ind].sname+" \\textless{}lower\\textgreater{} \\textless{}upper\\textgreater{}");
    pstruct[ind].entries[0][1].push_back ("Set lower and upper cutoffs for uv-plane data to \\textless{}lower\\textgreater{} and \\textless{}upper\\textgreater{}.");


    if (ind+1!=*size)
        PRINTER printerror("pixsrc", "internal error: bad num parms. set: " +
                           OPERA tostring (*size) +". actual: "+OPERA tostring (ind+1), NULL);

    if (parmslist)
    {
        // create list of parameters
        char **plist;
        MEMORY ps_malloc( &plist, *size, 100 );
        *parmslist=plist;
        // copy parm names from structs into list
        for (PS_SIT i=0; i<*size; ++i)
            strcpy (plist[i], pstruct[i].sname.c_str());
    }

    if (pstructlist)
    {
        *pstructlist = pstruct;
    }
    else
    {
        delete pstruct->pindex;
        delete [] pstruct;
    }
}

void pixsrc_init::setdefaultparameters(char *bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_)
{
    MEMORY ps_malloc( &(cdata_->coordsysorig ), 100 );
    MEMORY ps_malloc( &(cdata_->coordsysfinal), 100 );
    strcpy(cdata_->coordsysorig,"J2000");

    cdata_->numgpudevices = 0;
    cdata_->gpudevices    = 0;

    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        data_[g].numgpu2use = 0;
        data_[g].gpu2use    = 0;

        data_[g].fwhm_weight = 0;
        data_[g].useall = 1;
        data_[g].srcpixelscale = NULL;
        data_[g].shapeintmask = NULL;
        OPERA assign_n_infinity( &data_[g].evidence );
        data_[g].extlengths[17] = data_[g].extlengths[18] = 0;

        data_[g].pra = -1;
        data_[g].pdec = -1;
        data_[g].px = 1;
        data_[g].py = 1;
        data_[g].r1 = 1;
        data_[g].r2 = 1;

        data_[g].cartdetails = NULL;
        data_[g].srcinputcircle = NULL;
        data_[g].psf = NULL;
        data_[g].mmimages = NULL;
        data_[g].chi2mask = NULL;
        data_[g].mags = NULL;
        data_[g].traces = NULL;
        data_[g].details = NULL;
        data_[g].detailsptr = NULL;

        data_[g].usersetsrc        = 0;
        data_[g].srcbounds         = 0;
        data_[g].srcstepsizes      = 0;
        data_[g].srcnvary          = 0;
        data_[g].numsrcbounds      = 0;
        data_[g].num_vec_src       = 0;
        data_[g].vector_src_numpix = 0;
        data_[g].vector_src_pos    = 0;
        data_[g].vector_src_flux   = 0;
        data_[g].vector_src_triout = NULL;
        data_[g].rotation_axis     = NULL;

        // interferometric data stuff
        data_[g].uv_model_filename  = NULL;
        data_[g].uv_oldloc          = NULL;
        data_[g].uv_data            = NULL;
        data_[g].uv_idx             = NULL;
        data_[g].uv_transform_vec_ptr   = NULL;
        data_[g].uv_transform_mat_ptr   = NULL;
        data_[g].uv_transform_vec   = NULL;
        data_[g].uv_transform_mat   = NULL;
        data_[g].uv_newvariance     = NULL;
        data_[g].uv_ndp             = 0;
        data_[g].uv_deltapts        = 30000;
        data_[g].uv_maxintndp       = -10000;
        data_[g].transmode          = 1;
        data_[g].uv_matvecsca_exist = 0;

        // penalty variables
        data_[g].extlengths[13] = 1;
        data_[g].extlengths[5] = 6;
        data_[g].extlengths[6] = 6;
        data_[g].penaltymatrix = 0;
        data_[g].penaltyfilenames = 0;
        data_[g].penaltynames = 0;
        MEMORY ps_malloc (&data_[g].penaltyquery, data_[g].extlengths[5]);
        std::fill (data_[g].penaltyquery, data_[g].penaltyquery + data_[g].extlengths[5], 0);

        // set masks to null
        data_[g].imagemasks = 0;
        data_[g].chi2mask = 0;
        data_[g].num_mmimages = 0;
        data_[g].mmborder = 0;
        data_[g].mmborder_defl = 0;
        data_[g].mmimages = 0;
        data_[g].shapeintmaskminmax = 0;
        data_[g].srcinputcircle = 0;

	// hack src
	data_[g].hacksrcfn = NULL;
	data_[g].hacksrc   = NULL;
    }


    // set more defaults
    PS_SIT size;
    ps_parms_struct *pstruct;
    INIT parmscreator (NULL, &pstruct, &size);
    ps_parms_ind_class *pind = pstruct->pindex;

    cdata_->numthreads         = (PS_SIT)pstruct[pind->ps_parm_ind_threads].fvalue[0];
    cdata_->npl_stepsize       =      pstruct[pind->ps_parm_ind_nplstepsize].fvalue[0];
    cdata_->npl_ftolsize       =      pstruct[pind->ps_parm_ind_nplftol].fvalue[0];

    cdata_->npl             = 0;
    cdata_->npl_num_ref_pot = 0;
    cdata_->npl_filename    = NULL;
    cdata_->npl_ctr[0]      = cdata_->npl_ctr[1]     = 0;
    cdata_->npl_size[0]     = cdata_->npl_size[1]    = 0;
    cdata_->npl_num_pts[0]  = cdata_->npl_num_pts[1] = 0;
    cdata_->npl_reg_parms[0]   = pstruct[pind->ps_parm_ind_npltpsreg].fvalue[0];
    cdata_->npl_reg_parms[1]   = pstruct[pind->ps_parm_ind_npltpsreg].fvalue[1];
    if (pstruct[pind->ps_parm_ind_nonparamlens].fvalue[0])
        PRINTER printerror("",pstruct[pind->ps_parm_ind_nonparamlens].sname +" bad default value.",cdata_->print2screenmutex);

    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        data_[g].debug              = (PS_SIT)pstruct[pind->ps_parm_ind_debug].fvalue[0];
        data_[g].fatalwarn          = (PS_SIT)pstruct[pind->ps_parm_ind_fatalwarn].fvalue[0];
        data_[g].regorder           = (PS_SIT)pstruct[pind->ps_parm_ind_regorder].fvalue[0];
        data_[g].lambdaguess        =      pstruct[pind->ps_parm_ind_regstrength].fvalue[0];
        data_[g].findsisterimages   = (PS_SIT)pstruct[pind->ps_parm_ind_sisterpix].fvalue[0];
        data_[g].subsampling        = (PS_SIT)pstruct[pind->ps_parm_ind_oversample].fvalue[0];
        data_[g].myvariance0        =         pstruct[pind->ps_parm_ind_noise].fvalue[0]
	    *pstruct[pind->ps_parm_ind_noise].fvalue[0];
        data_[g].magparams          = (PS_SIT)pstruct[pind->ps_parm_ind_magnification].fvalue[0];
        data_[g].maguncertainty     = (PS_SIT)pstruct[pind->ps_parm_ind_mag_uncertainty].fvalue[0];
        data_[g].fullsrccov         = (PS_SIT)pstruct[pind->ps_parm_ind_fullsrccov].fvalue[0];
        data_[g].magsamples         = (PS_SIT)pstruct[pind->ps_parm_ind_magsamples].fvalue[0];
        data_[g].verbose            = (PS_SIT)pstruct[pind->ps_parm_ind_verbosity].fvalue[0];
        data_[g].printvec           = (PS_SIT)pstruct[pind->ps_parm_ind_images].fvalue[0];
        data_[g].printdetails       = (PS_SIT)pstruct[pind->ps_parm_ind_details].fvalue[0];
        data_[g].regaccuracy        = (PS_SIT)pstruct[pind->ps_parm_ind_regaccuracy].fvalue[0];
        data_[g].reg                = (PS_SIT)pstruct[pind->ps_parm_ind_reg].fvalue[0];
        data_[g].rmpix              = (PS_SIT)pstruct[pind->ps_parm_ind_rmpix].fvalue[0];
        data_[g].interperr          = (PS_SIT)pstruct[pind->ps_parm_ind_interperr].fvalue[0];
        data_[g].irscheme           = (PS_SIT)pstruct[pind->ps_parm_ind_irscheme].fvalue[0];
        data_[g].srcftol            =      pstruct[pind->ps_parm_ind_srcftol].fvalue[0];
        data_[g].srcexact           = (PS_SIT)pstruct[pind->ps_parm_ind_srcexact].fvalue[0];
        data_[g].optfinder          = (PS_SIT)pstruct[pind->ps_parm_ind_optfinder].fvalue[0];
        data_[g].lambda_lens        =      pstruct[pind->ps_parm_ind_reglens].fvalue[0];
        data_[g].noisemap           = (PS_SIT)pstruct[pind->ps_parm_ind_noisemap].fvalue[0];
        data_[g].raydirection       = (PS_SIT)pstruct[pind->ps_parm_ind_raydirection].fvalue[0];
        data_[g].regftol            =      pstruct[pind->ps_parm_ind_regftol].fvalue[0];
        data_[g].srcpixelscale0     = (PS_SIT)pstruct[pind->ps_parm_ind_pixelscalesrc].fvalue[0];
        data_[g].precision          = (PS_SIT)pstruct[pind->ps_parm_ind_outprecision].fvalue[0];
        data_[g].fullmag            = (PS_SIT)pstruct[pind->ps_parm_ind_fullmag].fvalue[0];
        data_[g].lowmem             = (PS_SIT)pstruct[pind->ps_parm_ind_lowmem].fvalue[0];
        data_[g].uv_taper           =      pstruct[pind->ps_parm_ind_uvtaper].fvalue[0];
        data_[g].uv_newpixelscale   =      pstruct[pind->ps_parm_ind_uvnewpixelscale].fvalue[0];
        data_[g].uv_padzeros        = (PS_SIT)pstruct[pind->ps_parm_ind_uvpadzeros].fvalue[0];
        data_[g].uv_rbftype         = (PS_SIT)pstruct[pind->ps_parm_ind_uvrbf].fvalue[0];
        data_[g].uv_rbfscale        =      pstruct[pind->ps_parm_ind_uvrbf].fvalue[1];
        data_[g].uv_del1            = (PS_SIT)pstruct[pind->ps_parm_ind_uvmatrixsize].fvalue[0];
        data_[g].uv_del2            = (PS_SIT)pstruct[pind->ps_parm_ind_uvmatrixsize].fvalue[1];
        data_[g].uv_pbeam           = (PS_SIT)pstruct[pind->ps_parm_ind_uvpbeam].fvalue[0];
        data_[g].srcrestart         = (PS_SIT)pstruct[pind->ps_parm_ind_srcrestart].fvalue[0];
        data_[g].guess_src_position = (PS_SIT)pstruct[pind->ps_parm_ind_srcctrguess].fvalue[0];



        data_[g].gridtype       = (PS_SIT)pstruct[pind->ps_parm_ind_grid].fvalue[0];
        if      (1==data_[g].gridtype) data_[g].levelshift = pstruct[pind->ps_parm_ind_grid].fvalue[1];
        else if (2==data_[g].gridtype) data_[g].levelshift = pstruct[pind->ps_parm_ind_grid].fvalue[2];
        else PRINTER printerror("",pstruct[pind->ps_parm_ind_grid].sname +" bad default value.",cdata_->print2screenmutex);

        switch ((PS_SIT)pstruct[pind->ps_parm_ind_psf].fvalue[0])
        {
        case 0:
            data_[g].nopsf       = 1;
            data_[g].psffromfile = 0;
            break;
        case 1:
            data_[g].nopsf       = 0;
            data_[g].psffromfile = 1;
            break;
        default:
            PRINTER printerror("",pstruct[pind->ps_parm_ind_psf].sname +" bad default value.",cdata_->print2screenmutex);
            break;
        }
        data_[g].psf_oversample   = (PS_SIT)pstruct[pind->ps_parm_ind_psf].fvalue[1];
        data_[g].majoraxis = data_[g].minoraxis = data_[g].angleaxis = -1;

        switch ((PS_SIT)pstruct[pind->ps_parm_ind_statistic].fvalue[0])
        {
        case 0:
            data_[g].onlychi = 0;
            data_[g].noevinochi = 0;
            break;
        case 1:
            data_[g].onlychi = 1;
            data_[g].noevinochi = 0;
            break;
        case 2:
            data_[g].onlychi = 0;
            data_[g].noevinochi = 1;
            break;
        default:
            PRINTER printerror("",pstruct[pind->ps_parm_ind_statistic].sname +" bad default value.",cdata_->print2screenmutex);
        }

        data_[g].mintriangleangle = pstruct[pind->ps_parm_ind_minangle].fvalue[0];
        MEMORY ps_malloc (&data_[g].triswitches, 100);
        std::fill (data_[g].triswitches, data_[g].triswitches+100, 0);
        strcpy(data_[g].triswitches,string("zQenq"+pstruct[pind->ps_parm_ind_minangle].svalue).c_str());

        data_[g].extlengths[4] = 3;
        MEMORY ps_malloc (&data_[g].traceparams, data_[g].extlengths[4]);
        std::copy (pstruct[pind->ps_parm_ind_varystrength].fvalue.begin(),
                   pstruct[pind->ps_parm_ind_varystrength].fvalue.end(), data_[g].traceparams);
        if (data_[g].traceparams[2])
            PRINTER printerror("",pstruct[pind->ps_parm_ind_varystrength].sname +" bad default value.",cdata_->print2screenmutex);

        data_[g].fillbadpix = NULL;
        data_[g].fillbadpixnoise = (PS_SIT)pstruct[pind->ps_parm_ind_fillbadpix].fvalue[0];
        if (data_[g].fillbadpixnoise)
            PRINTER printerror("",pstruct[pind->ps_parm_ind_fillbadpix].sname +" bad default value.",cdata_->print2screenmutex);

        data_[g].use_shapelets  = (PS_SIT)pstruct[pind->ps_parm_ind_shapelets].fvalue[0];
        data_[g].num_shapelets1 = (PS_SIT)pstruct[pind->ps_parm_ind_shapelets].fvalue[1];
        data_[g].num_shapelets2 = (PS_SIT)pstruct[pind->ps_parm_ind_shapelets].fvalue[2];
        data_[g].numberofshapelets = data_[g].num_shapelets1*data_[g].num_shapelets2;

        std::copy (pstruct[pind->ps_parm_ind_shapeparms].fvalue.begin(),
                   pstruct[pind->ps_parm_ind_shapeparms].fvalue.end(), data_[g].fixedshapeletparms);
        if (data_[g].fixedshapeletparms[0])
            PRINTER printerror("",pstruct[pind->ps_parm_ind_shapeparms].sname +" bad default value.",cdata_->print2screenmutex);

	data_[g].shapeletpixsplit[0] = (PS_SIT)pstruct[pind->ps_parm_ind_shapepixsplit].fvalue[0];
	data_[g].shapeletpixsplit[1] = (PS_SIT)pstruct[pind->ps_parm_ind_shapepixsplit].fvalue[1];

        data_[g].rottrans[0] = 1; // cosine term
        std::fill (data_[g].rottrans+1, data_[g].rottrans+6, 0);
        for (PS_SIT h=0; h<(PS_SIT)pstruct[pind->ps_parm_ind_rottrans].fvalue.size(); ++h)
            if (pstruct[pind->ps_parm_ind_rottrans].fvalue[h])
                PRINTER printerror("",pstruct[pind->ps_parm_ind_rottrans].sname +" bad default value.",cdata_->print2screenmutex);

        data_[g].is_uvdata      = 0;
        data_[g].uv_mode        = 0;
        data_[g].uv_pointing[0] = 0;
        data_[g].uv_pointing[1] = 0;
        data_[g].uv_filename    = NULL;
        if (pstruct[pind->ps_parm_ind_uvdata].fvalue[0])
            PRINTER printerror("",pstruct[pind->ps_parm_ind_uvdata].sname +" bad default value.",cdata_->print2screenmutex);

        std::copy (pstruct[pind->ps_parm_ind_uvcutoff].fvalue.begin(),
                   pstruct[pind->ps_parm_ind_uvcutoff].fvalue.end(), data_[g].uv_cutoff);
    }

    delete pstruct->pindex;
    delete [] pstruct;


    // read mm-image masks (need to this now in order to read pixsrc parameters)
    INIT readmmmasks (bn, name, namewithext, data_, cdata_);
}

void pixsrc_init::readparameters(char *bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_, PS_SIT dimil, char **ignorelist, char **ignorelistwithext)
{
    // read parameter list

    char **parmslist, **vec, *fname;
    PS_SIT parmslistsize, vecsize;
    const char *listcc[3];
    listcc[0] = CONSTANT dir_in;
    listcc[1] = bn;
    listcc[2] = CONSTANT parameterfile;
    OPERA concatenate( listcc, 3, &fname );

    ps_parms_struct *pstruct;
    INIT parmscreator( &parmslist, &pstruct, &parmslistsize );
    OPERA readfile( fname, &vec, &vecsize, cdata_->print2screenmutex);
    MEMORY ps_free( fname );

    PS_SIT imageindex = -1; // applies to all images
    PS_SIT penaltyindex = -1;
    double *input = NULL;
    PS_SIT dim;
    ps_parms_ind_class *pind = pstruct->pindex;
    for(PS_SIT line=0; line<vecsize; line++)
    {
        if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_regorder],
                      OPERA sizestring(parmslist[pind->ps_parm_ind_regorder]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_regorder]),
                          "i",pstruct[pind->ps_parm_ind_regorder].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].regorder = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_dataname],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_dataname]) ) )
        {
            OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_dataname]), NULL );

            // figure out if ignoring this image
            PS_SIT ignoring = 0;
            for(PS_SIT j=0; j<dimil; ++j)
            {
                if( !strcmp(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_dataname]), ignorelist[j]) ||
                    !strcmp(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_dataname]), ignorelistwithext[j]) )
                {
                    imageindex = -2;
                    ignoring = 1;
                    break;
                }
            }
            if (!ignoring)
            {
                // figure out which images the parameters apply to:
                PS_SIT foundone=0;
                for(PS_SIT g=0; g<cdata_->numimages; g++)
                {
                    if( !strcmp( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_dataname]), namewithext[g] ) ||
                        !strcmp( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_dataname]), name[g]) )
                    {
                        imageindex = g;
                        foundone = 1;
                        break;
                    }
                }
                if (!foundone)
                    PRINTER printerror("","data name not found: " +
                                       OPERA tostring (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_dataname])),
                                       cdata_->print2screenmutex);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penaltyname],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penaltyname]) ) )
        {
            if (-1==imageindex)
                PRINTER printerror("","You must set an image before you set penalty.",
                                   cdata_->print2screenmutex);
            // if image hasn't been set, read mm-image names from first image
            //PS_SIT imgind = -1==imageindex ? 0 : imageindex;

            OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penaltyname]), NULL );
            PS_SIT imgind = imageindex;
            PS_SIT foundit = 0;
            for (PS_SIT penind=0; penind<data_[imgind].extlengths[11]; ++penind)
            {
                char *longname;
                const char *listcc[6];
                listcc[0] = bn;
                listcc[1] = CONSTANT bnseparator;
                listcc[2] = name[imgind];
                listcc[3] = CONSTANT imagesfilestart;
                listcc[4] = data_[imgind].penaltynames[penind];
                listcc[5] = CONSTANT imagesfileendreg;
                OPERA concatenate (listcc, 6, &longname);

                if (!strcmp (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penaltyname]),
                             data_[imgind].penaltynames[penind]) ||
                    !strcmp (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penaltyname]),
                             longname))
                {
                    penaltyindex = penind;
                    foundit = 1;
                    break;
                }
                MEMORY ps_free (longname);
            }
            if (!foundit)
                PRINTER printerror("","Penalty file name not found: " +
                                   OPERA tostring (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penaltyname])),
                                   cdata_->print2screenmutex);
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_nonparamlens],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]) ) )
        {
            if(imageindex != -1)
                PRINTER printerror("","You cannot specify "
                                   "non-parametric lens model for just one image.\n"
                                   "Must specify for all images.",
                                   cdata_->print2screenmutex);

	    /*
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]),
                          "i s d d d i i Z",pstruct[pind->ps_parm_ind_nonparamlens].sname.c_str(),
                          &input,cdata_);
            cdata_->npl = (PS_SIT)input[0];

            // read npl grid
            if (cdata_->npl==1 || cdata_->npl==2)
            {
                cdata_->npl_ctr[0] = input[1];
                cdata_->npl_ctr[1] = input[2];
                cdata_->npl_size[0] = input[3];
                cdata_->npl_size[1] = input[4];
                cdata_->npl_num_pts[0] = (PS_SIT)input[5];
                cdata_->npl_num_pts[1] = (PS_SIT)input[6];
            }
            // read filename of previous npl run
            PS_SIT index_string = 2==cdata_->npl ? 7 : 1;
            if (cdata_->npl==2 || cdata_->npl==3)
            {
                OPERA trim (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]), NULL);

                // copy filenames
                char *filenames;
                PS_SIT strlength = OPERA sizestring (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]));
                MEMORY ps_malloc (&filenames, strlength+1);
                strcpy (filenames, vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]));

                // find first filename
                char *ptr_start = vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]);
                char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]), " " );
                for (PS_SIT n=0; n<index_string; ++n)
                    ptr = strtok( NULL, " " );

                // parse filenames
                MEMORY ps_free (cdata_->npl_filename);
                OPERA split (filenames+(ptr-ptr_start), " ",
                             &cdata_->npl_filename, &cdata_->npl_num_ref_pot);
                MEMORY ps_free (filenames);
            }
	    */


            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]),
                          "i Z",pstruct[pind->ps_parm_ind_nonparamlens].sname.c_str(),
                          &input,cdata_);
            cdata_->npl = (PS_SIT)input[0];

            // read npl grid
            if (cdata_->npl==1)
            {
		dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]),
			      "i s d d d i i Z",pstruct[pind->ps_parm_ind_nonparamlens].sname.c_str(),
			      &input,cdata_);
		if (dim<7)
                PRINTER printerror("","Must specify grid for non-parameteric lens potential.",
                                   cdata_->print2screenmutex);

                cdata_->npl_ctr[0] = input[1];
                cdata_->npl_ctr[1] = input[2];
                cdata_->npl_size[0] = input[3];
                cdata_->npl_size[1] = input[4];
                cdata_->npl_num_pts[0] = (PS_SIT)input[5];
                cdata_->npl_num_pts[1] = (PS_SIT)input[6];
            }
            // read filename of previous npl run
            PS_SIT index_string = 1==cdata_->npl ? 7 : 1;
            if ((cdata_->npl==1 && dim>7) || cdata_->npl==2) //cdata_->npl==2 || cdata_->npl==3)
            {
                OPERA trim (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]), NULL);

                // copy filenames
                char *filenames;
                PS_SIT strlength = OPERA sizestring (vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]));
                MEMORY ps_malloc (&filenames, strlength+1);
                strcpy (filenames, vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]));

                // find first filename
                char *ptr_start = vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]);
                char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nonparamlens]), " " );
                for (PS_SIT n=0; n<index_string; ++n)
                    ptr = strtok( NULL, " " );

                // parse filenames
                MEMORY ps_free (cdata_->npl_filename);
                OPERA split (filenames+(ptr-ptr_start), " ",
                             &cdata_->npl_filename, &cdata_->npl_num_ref_pot);
		if (2==cdata_->npl)
		    --cdata_->npl_num_ref_pot;
                MEMORY ps_free (filenames);
            }
	    else 
	    {
		MEMORY ps_free (cdata_->npl_filename);
		cdata_->npl_filename = NULL;
		cdata_->npl_num_ref_pot = 0;
	    }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_nplstepsize],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_nplstepsize]) ) )
        {
            if(imageindex != -1)
                PRINTER printerror("","Currently, you cannot specify non-parametric lens model step size for just one image.\n"
                                   "Must specify for all images.",
                                   cdata_->print2screenmutex);

            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nplstepsize]),"d",pstruct[pind->ps_parm_ind_nplstepsize].sname.c_str(),&input,cdata_);
            cdata_->npl_stepsize = input[0];
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_nplftol],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_nplftol]) ) )
        {
            if(imageindex != -1)
                PRINTER printerror("","You cannot specify simplex size for convergence for just one image.\n"
                                   "Must specify for all images.",
                                   cdata_->print2screenmutex);

            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_nplftol]),"d",pstruct[pind->ps_parm_ind_nplftol].sname.c_str(),&input,cdata_);
            cdata_->npl_ftolsize = input[0];
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_npltpsreg],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_npltpsreg]) ) )
        {
            if(imageindex != -1)
                PRINTER printerror("","You cannot specify regularization strength parameters for just one image.\n"
                                   "Must specify for all images.",
                                   cdata_->print2screenmutex);

            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_npltpsreg]),"d d",pstruct[pind->ps_parm_ind_npltpsreg].sname.c_str(),&input,cdata_);
            cdata_->npl_reg_parms[0] = input[0];
            cdata_->npl_reg_parms[1] = input[1];
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_cuda],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_cuda]) ) )
        {
            OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_cuda]), NULL );

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;

#ifdef __USE_PIXSRC_CUDA__
                CUDA detectgpus (cdata_);
                INIT setgpus( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_cuda]), &data_[g], cdata_);
#else
                PRINTER printerror("","This version of pixsrc does not "
                                   "support CUDA (GPU computing).\n"
                                   "The \"USE_PIXSRC_CUDA\" environmental variable "
                                   "must be set to \"YES\" to compile CUDA code.",
                                   cdata_->print2screenmutex);
#endif

            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_src],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_src]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_src]),"s",pstruct[pind->ps_parm_ind_src].sname.c_str(),&input,cdata_);
            OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_src]), NULL );
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                INIT setsource( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_src]), &data_[g], cdata_);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_srcctrguess],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_srcctrguess]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcctrguess]),"i",pstruct[pind->ps_parm_ind_srcctrguess].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].guess_src_position = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_srcbound],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_srcbound]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcbound]),"s",pstruct[pind->ps_parm_ind_srcbound].sname.c_str(),&input,cdata_);
            OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcbound]), NULL );
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                INIT setsourcebounds( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcbound]), &data_[g], cdata_);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_srcstepsize],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_srcstepsize]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcstepsize]),"s",pstruct[pind->ps_parm_ind_srcstepsize].sname.c_str(),&input,cdata_);
            OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcstepsize]), NULL );
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                INIT setsourcestepsizes( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcstepsize]), &data_[g], cdata_);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_reg],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_reg]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_reg]),"i",pstruct[pind->ps_parm_ind_reg].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].reg = (char)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_noisemap],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_noisemap]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_noisemap]),"i",pstruct[pind->ps_parm_ind_noisemap].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].noisemap = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_raydirection],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_raydirection]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_raydirection]),"i",pstruct[pind->ps_parm_ind_raydirection].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].raydirection = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_outprecision],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_outprecision]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_outprecision]),"i",pstruct[pind->ps_parm_ind_outprecision].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].precision = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_fullmag],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_fullmag]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_fullmag]),"i",pstruct[pind->ps_parm_ind_fullmag].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].fullmag = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_pixelscalesrc],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_pixelscalesrc]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_pixelscalesrc]),"d",pstruct[pind->ps_parm_ind_pixelscalesrc].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                if (!data_[g].srcpixelscale)
                    MEMORY ps_malloc (&data_[g].srcpixelscale, 1);
                data_[g].srcpixelscale[0] = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_regftol],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_regftol]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_regftol]),"d",pstruct[pind->ps_parm_ind_regftol].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].regftol = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_shapelets],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_shapelets]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_shapelets]),"i i s",pstruct[pind->ps_parm_ind_shapelets].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].use_shapelets = (PS_SIT)input[0];
                if (data_[g].use_shapelets==1)
                {
                    data_[g].num_shapelets1 = (PS_SIT)input[1];
                    data_[g].num_shapelets2 = (PS_SIT)input[2];
                    data_[g].numberofshapelets = data_[g].num_shapelets1 * data_[g].num_shapelets2;
                }
                else if (data_[g].use_shapelets==2)
                {
                    PRINTER printerror("pixsrc",
                                       "polar shapelets disabled",
                                       cdata_->print2screenmutex);
                    /*
                      data_[g].num_shapelets1 = (PS_SIT)input[1];
                      data_[g].numberofshapelets = data_[g].num_shapelets1 * (data_[g].num_shapelets1+1)/2;
                    */
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_shapeparms],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_shapeparms]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_shapeparms]),"i d d d",pstruct[pind->ps_parm_ind_shapeparms].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].fixedshapeletparms[0] = input[0];

                if (data_[g].fixedshapeletparms[0]==1 || data_[g].fixedshapeletparms[0]==3)
                {
                    data_[g].fixedshapeletparms[1] = input[1];
                    data_[g].fixedshapeletparms[2] = input[2];
                }
                if (data_[g].fixedshapeletparms[0]==2 || data_[g].fixedshapeletparms[0]==3)
                {
                    data_[g].fixedshapeletparms[3] = input[3];
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_shapepixsplit],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_shapepixsplit]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_shapepixsplit]),"i i",pstruct[pind->ps_parm_ind_shapepixsplit].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;

                data_[g].shapeletpixsplit[0] = (PS_SIT)input[0];
                data_[g].shapeletpixsplit[1] = (PS_SIT)input[1];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_srcrestart],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_srcrestart]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcrestart]),"i",pstruct[pind->ps_parm_ind_srcrestart].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].srcrestart = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_rmpix],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_rmpix]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_rmpix]),"i",pstruct[pind->ps_parm_ind_rmpix].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].rmpix = (char)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_interperr],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_interperr]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_interperr]),"i",pstruct[pind->ps_parm_ind_interperr].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].interperr = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_irscheme],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_irscheme]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_irscheme]),"i",pstruct[pind->ps_parm_ind_irscheme].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].irscheme = (char)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_srcftol],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_srcftol]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcftol]),"d",pstruct[pind->ps_parm_ind_srcftol].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].srcftol = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_srcexact],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_srcexact]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_srcexact]),"i",pstruct[pind->ps_parm_ind_srcexact].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].srcexact = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_regaccuracy],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_regaccuracy]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_regaccuracy]),"i",pstruct[pind->ps_parm_ind_regaccuracy].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].regaccuracy = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_verbosity],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_verbosity]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_verbosity]),"i",pstruct[pind->ps_parm_ind_verbosity].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].verbose = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_optfinder],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_optfinder]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_optfinder]),"i",pstruct[pind->ps_parm_ind_optfinder].sname.c_str(),&input,cdata_);
            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].optfinder = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_grid],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_grid]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_grid]),"i d d d d d i i",pstruct[pind->ps_parm_ind_grid].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].gridtype = (PS_SIT)input[0];
            }

            if(input[0]==0 && ((dim==2 && input[1]==1) || (dim==8 && input[1]==0)) )
            {
                bool autocc = (PS_SIT)input[1];

                for(PS_SIT g=0; g<cdata_->numimages; g++)
                {
                    if(imageindex!=-1 && imageindex!=g)
                        continue;
                    data_[g].extlengths[1]=0;
                    if(!autocc)
                    {
                        data_[g].extlengths[1]=6;
                        MEMORY ps_free( data_[g].cartdetails );
                        MEMORY ps_malloc( &(data_[g].cartdetails), data_[g].extlengths[1] );
                        std::copy( input+2, input+2+data_[g].extlengths[1], data_[g].cartdetails );
                    }
                }
            }
            else if(input[0]==1 && dim==1)
            {
                for(PS_SIT g=0; g<cdata_->numimages; g++)
                {
                    if(imageindex!=-1 && imageindex!=g)
                        continue;
                    data_[g].levelshift = pstruct[pind->ps_parm_ind_grid].fvalue[1];
                }
            }
            else if(input[0]==2 && dim==1)
            {
                for(PS_SIT g=0; g<cdata_->numimages; g++)
                {
                    if(imageindex!=-1 && imageindex!=g)
                        continue;
                    data_[g].levelshift = pstruct[pind->ps_parm_ind_grid].fvalue[2];
                }
            }
            else if( ( input[0]==1 || input[0]==2 ) && dim==2 )
            {
                for(PS_SIT g=0; g<cdata_->numimages; g++)
                {
                    if(imageindex!=-1 && imageindex!=g)
                        continue;
                    data_[g].levelshift=input[1];
                }
            }
            else
            {
                PRINTER printerror("","invalid grid details!", cdata_->print2screenmutex);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_psf],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_psf]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_psf]),"i d d d i",pstruct[pind->ps_parm_ind_psf].sname.c_str(),&input,cdata_);

            double majoraxis=-1, minoraxis=-1, angleaxis=-1;
            PS_SIT oversample=11;
            bool nopsf=1, psffromfile=0;

            if(input[0]==0 && dim==1)
            {}
            else if(input[0]==1 && dim==1)
            {
                nopsf=0;
                psffromfile=1;
            }
            else if(input[0]==2 && dim>=4)
            {
                nopsf=0;
                majoraxis=input[1];
                minoraxis=input[2];
                angleaxis=input[3];
            }
            else
            {
                PRINTER printerror("","invalid psf details!", cdata_->print2screenmutex);
            }

            if (input[0]==2 && dim==5)
            {
                oversample = (PS_SIT)input[4];
                if (oversample<=0)
                    PRINTER printerror("","psf: oversample factor must be a positive integer.", cdata_->print2screenmutex);
            }

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].nopsf = nopsf;
                data_[g].psffromfile = psffromfile;
                if (input[0]==2 && dim>=4)
                {
                    data_[g].majoraxis = majoraxis;
                    data_[g].minoraxis = minoraxis;
                    data_[g].angleaxis = angleaxis;
                }
                if (input[0]==2 && dim==5)
                    data_[g].psf_oversample = oversample;
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_coorsys],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_coorsys]) ) )
        {
            if(imageindex != -1)
                PRINTER printerror("","Coordinate system is a global property!\n"
                                   "You cannot define it for a particular image.\n"
                                   "To shift a particular image, use the \"rottrans\" parameter",
                                   cdata_->print2screenmutex);

            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_coorsys]),
                          "s d d d d",pstruct[pind->ps_parm_ind_coorsys].sname.c_str(),&input,cdata_);

            char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_coorsys]), " " );
            strcpy( cdata_->coordsysorig, ptr );

            double p1=-1, p2=-1, r1=-1, r2=-1;

            if(!strcmp(cdata_->coordsysorig,"PIXEL") && dim==1)
            {
                p1 = p2 = r1 = r2 = 1;
            }
            else if(dim==5)
            {
                p1 = input[1];
                p2 = input[2];
                r1 = input[3];
                r2 = input[4];
            }
            else
            {
                PRINTER printerror("","invalid coordinate system defined!",
                                   cdata_->print2screenmutex);
            }

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                data_[g].px = data_[g].pra = p1;
                data_[g].py = data_[g].pdec = p2;
                data_[g].r1 = r1;
                data_[g].r2 = r2;
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_regstrength],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_regstrength]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_regstrength]),"d",pstruct[pind->ps_parm_ind_regstrength].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].lambdaguess = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_reglens],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_reglens]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_reglens]),"d",pstruct[pind->ps_parm_ind_reglens].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].lambda_lens = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_oversample],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_oversample]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_oversample]),"i",pstruct[pind->ps_parm_ind_oversample].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].subsampling = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_debug],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_debug]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_debug]),"i",pstruct[pind->ps_parm_ind_debug].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].debug = (input[0] != 0);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_fatalwarn],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_fatalwarn]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_fatalwarn]),"i",pstruct[pind->ps_parm_ind_fatalwarn].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].fatalwarn = (input[0] != 0);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_details],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_details]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_details]),"i",pstruct[pind->ps_parm_ind_details].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].printdetails = (input[0] != 0);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_images],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_images]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_images]),"i",pstruct[pind->ps_parm_ind_images].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].printvec = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_sisterpix],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_sisterpix]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_sisterpix]),"i",pstruct[pind->ps_parm_ind_sisterpix].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].findsisterimages = ( input[0] != 0);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_fillbadpix],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_fillbadpix]) ) )
        {
            dim =  readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_fillbadpix]),
                           "i s",pstruct[pind->ps_parm_ind_fillbadpix].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                if(input[0])
                {
                    MEMORY ps_free( data_[g].fillbadpix );
                    MEMORY ps_malloc( &(data_[g].fillbadpix), 1 );

                    // copy g=0 image to other images
                    if (imageindex==-1 && g)
                    {
                        data_[g].fillbadpixnoise = data_[0].fillbadpixnoise;
                        data_[g].fillbadpix[0] = data_[0].fillbadpix[0];
                    }
                    else
                    {
                        char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_fillbadpix]), " " );
                        ptr = strtok( NULL, " " );
                        data_[g].fillbadpixnoise = 0;
                        if( !strcmp(ptr,"noise") )
                            data_[g].fillbadpixnoise = 1;
                        else
                            data_[g].fillbadpix[0] = input[1];
                    }
                }
                else if(data_[g].fillbadpix)
                {
                    MEMORY ps_free( data_[g].fillbadpix );
                    data_[g].fillbadpix = NULL;
                    data_[g].fillbadpixnoise = 0;
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvdata],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvdata]) ) )
        {
            dim =  readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvdata]),
                           "i s d d",pstruct[pind->ps_parm_ind_uvdata].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;

                if (dim<1)
                    PRINTER printerror("","invalid uv data input", cdata_->print2screenmutex);

                data_[g].is_uvdata = (PS_SIT)input[0];
                if(input[0])
                {
                    if (4!=dim)
                        PRINTER printerror("","invalid uv data input", cdata_->print2screenmutex);

                    MEMORY ps_free (data_[g].uv_filename);
                    MEMORY ps_malloc( &data_[g].uv_filename, OPERA sizestring(vec[line])+1 );

                    char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvdata]), " " );
                    ptr = strtok( NULL, " " );
                    strcpy (data_[g].uv_filename, ptr);

                    data_[g].uv_mode        = (PS_SIT)input[0];
                    data_[g].transmode      = 1==input[0] ? 1 : 0;
                    data_[g].uv_pointing[0] = input[2];
                    data_[g].uv_pointing[1] = input[3];
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_hacksrc],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_hacksrc]) ) )
        {
            dim =  readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_hacksrc]),
                           "i s",pstruct[pind->ps_parm_ind_hacksrc].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;

                if (dim<1)
                    PRINTER printerror("","invalid hacksrc input", cdata_->print2screenmutex);

		MEMORY ps_free (data_[g].hacksrcfn);
		data_[g].hacksrcfn = NULL;
		
                if(input[0])
                {
                    if (2!=dim)
                        PRINTER printerror("","invalid hacksrc input", cdata_->print2screenmutex);

                    MEMORY ps_malloc (&data_[g].hacksrcfn, OPERA sizestring(vec[line])+1);

                    char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvdata]), " " );
                    ptr = strtok( NULL, " " );
                    strcpy (data_[g].hacksrcfn, ptr);
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvmodelpos],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvmodelpos]) ) )
        {
            dim =  readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvmodelpos]),
                           "i s",pstruct[pind->ps_parm_ind_uvmodelpos].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;

                if (dim<1)
                    PRINTER printerror("","invalid uvmodelpos input", cdata_->print2screenmutex);

                if(input[0])
                {
                    if (2!=dim)
                        PRINTER printerror("","invalid uvmodelpos input", cdata_->print2screenmutex);

                    MEMORY ps_free (data_[g].uv_model_filename);
                    MEMORY ps_malloc (&data_[g].uv_model_filename, OPERA sizestring(vec[line])+1 );

                    char *ptr = strtok( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvmodelpos]), " " );
                    ptr = strtok( NULL, " " );
                    strcpy (data_[g].uv_model_filename, ptr);
                }
                else
                {
                    MEMORY ps_free (data_[g].uv_model_filename);
                    data_[g].uv_model_filename = NULL;
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_lowmem],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_lowmem]) ) )
        {
            dim =  readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_lowmem]),
                           "i",pstruct[pind->ps_parm_ind_lowmem].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;

                data_[g].lowmem = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_statistic],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_statistic]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_statistic]),
                          "i",pstruct[pind->ps_parm_ind_statistic].sname.c_str(),&input,cdata_);
            bool noevinochi=0, onlychi=0;

            if(input[0]==0 && dim==1)
            {}
            else if(input[0]==1 && dim==1)
            {
                onlychi=1;
            }
            else if(input[0]==2 && dim==1)
            {
                noevinochi=1;
            }
            else
            {
                PRINTER printerror("","invalid statistic answer!", cdata_->print2screenmutex);
            }

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].noevinochi = noevinochi;
                data_[g].onlychi = onlychi;
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_magnification],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_magnification]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_magnification]),"i",pstruct[pind->ps_parm_ind_magnification].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].magparams = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_mag_uncertainty],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_mag_uncertainty]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_mag_uncertainty]),"i",pstruct[pind->ps_parm_ind_mag_uncertainty].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].maguncertainty = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_fullsrccov],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_fullsrccov]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_fullsrccov]),"i",pstruct[pind->ps_parm_ind_fullsrccov].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].fullsrccov = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_magsamples],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_magsamples]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_magsamples]),"i",pstruct[pind->ps_parm_ind_magsamples].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].magsamples = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_varystrength],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_varystrength]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_varystrength]),"d d i",pstruct[pind->ps_parm_ind_varystrength].sname.c_str(),&input,cdata_);

            // change max-lambda to base of logarithm in lambda computation
            input[1] = std::pow(input[1]/input[0],1.0/input[2]);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                std::copy( input, input+data_[g].extlengths[4], data_[g].traceparams );
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penalty1],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penalty1]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penalty1]),"i d d d d d",pstruct[pind->ps_parm_ind_penalty1].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                PS_SIT penind = 0;
                std::copy( input, input+data_[g].extlengths[6], data_[g].penaltymatrix[penind][0] );
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_rottrans],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_rottrans]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_rottrans]),"d d d d d",pstruct[pind->ps_parm_ind_rottrans].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                std::copy (input+1, input+5, data_[g].rottrans+2);
                data_[g].rottrans[0] = std::cos (-input[0]*CONSTANT deg2rad);
                data_[g].rottrans[1] = std::sin (-input[0]*CONSTANT deg2rad);
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penalty2],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penalty2]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penalty2]),"i d d d d d",pstruct[pind->ps_parm_ind_penalty2].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                PS_SIT penind = 0;
                std::copy( input, input+data_[g].extlengths[6], data_[g].penaltymatrix[penind][1] );
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penalty4],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penalty4]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penalty4]),"i d d d d d",pstruct[pind->ps_parm_ind_penalty4].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                for (PS_SIT penind=0; penind<data_[g].extlengths[13]; ++penind)
                {
                    if (penaltyindex!=-1 && penaltyindex!=penind)
                        continue;
                    std::copy( input, input+data_[g].extlengths[6], data_[g].penaltymatrix[penind][3] );
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penalty5],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penalty5]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penalty5]),"i d d d d d",pstruct[pind->ps_parm_ind_penalty5].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                for (PS_SIT penind=0; penind<data_[g].extlengths[13]; ++penind)
                {
                    if (penaltyindex!=-1 && penaltyindex!=penind)
                        continue;
                    std::copy( input, input+data_[g].extlengths[6], data_[g].penaltymatrix[penind][4] );
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penalty3],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penalty3]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penalty3]),"i d d d d d",pstruct[pind->ps_parm_ind_penalty3].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                for (PS_SIT penind=0; penind<data_[g].extlengths[13]; ++penind)
                {
                    if (penaltyindex!=-1 && penaltyindex!=penind)
                        continue;
                    std::copy( input, input+data_[g].extlengths[6], data_[g].penaltymatrix[penind][2] );
                }
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_penalty6],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_penalty6]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_penalty6]),"i d d d d d",pstruct[pind->ps_parm_ind_penalty6].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                PS_SIT penind = 0;
                std::copy( input, input+data_[g].extlengths[6], data_[g].penaltymatrix[penind][5] );
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_noise],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_noise]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_noise]),"d",pstruct[pind->ps_parm_ind_noise].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].myvariance0 = input[0]*input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvtaper],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvtaper]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvtaper]),"d",pstruct[pind->ps_parm_ind_uvtaper].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_taper = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvrbf],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvrbf]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvrbf]),"i d",pstruct[pind->ps_parm_ind_uvrbf].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_rbftype  = (PS_SIT)input[0];
                data_[g].uv_rbfscale = input[1];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvmatrixsize],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvmatrixsize]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvmatrixsize]),"i i",pstruct[pind->ps_parm_ind_uvmatrixsize].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_del1 = (PS_SIT)input[0];
                data_[g].uv_del2 = (PS_SIT)input[1];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvpbeam],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvpbeam]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvpbeam]),"i",pstruct[pind->ps_parm_ind_uvpbeam].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_pbeam = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvnewpixelscale],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvnewpixelscale]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvnewpixelscale]),"d",pstruct[pind->ps_parm_ind_uvnewpixelscale].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_newpixelscale = input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvpadzeros],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvpadzeros]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvpadzeros]),"d",pstruct[pind->ps_parm_ind_uvpadzeros].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_padzeros = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_uvcutoff],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_uvcutoff]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_uvcutoff]),"d d",pstruct[pind->ps_parm_ind_uvcutoff].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].uv_cutoff[0] = input[0];
                data_[g].uv_cutoff[1] = input[1];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_threads],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_threads]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_threads]),"i",pstruct[pind->ps_parm_ind_threads].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                cdata_->numthreads = (PS_SIT)input[0];
            }
        }
        else if( !strncmp( vec[line], parmslist[pind->ps_parm_ind_minangle],
                           OPERA sizestring(parmslist[pind->ps_parm_ind_minangle]) ) )
        {
            dim = readarr(vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_minangle]),"d",pstruct[pind->ps_parm_ind_minangle].sname.c_str(),&input,cdata_);

            for(PS_SIT g=0; g<cdata_->numimages; g++)
            {
                if(imageindex!=-1 && imageindex!=g)
                    continue;
                data_[g].mintriangleangle = input[0];

                if( !input[0] )
                {
                    *(data_[g].triswitches+4) = 'N';
                    *(data_[g].triswitches+5) = 0;
                }
                else
                {
                    *(data_[g].triswitches+4) = 'q';
                    OPERA trim( vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_minangle]), NULL );
                    strcpy( data_[g].triswitches+5, vec[line]+OPERA sizestring(parmslist[pind->ps_parm_ind_minangle]) );
                }
            }
        }
        else
        {
            PRINTER printerror("","unrecognized option in parameters file: " +
                               string(vec[line]), cdata_->print2screenmutex);
        }
    }

    MEMORY ps_free( input                    );
    MEMORY ps_free( parmslist, parmslistsize );
    MEMORY ps_free( vec, vecsize             );
    delete pstruct->pindex;
    delete [] pstruct;
}
