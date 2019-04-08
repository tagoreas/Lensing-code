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
#include "pixsrc_geometry.hpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_constants.hpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_external.hpp"
#include "pixsrc_shapelets_operations.hpp"
#include <dirent.h> // for file/directory access
#include <cstring>
#include <algorithm>
#include <cmath>
#include <fitsio.h>
#include <time.h>

void pixsrc_init::finalizeNcounterindications(char *bn, char **name, char **namewithext,
                                              inputdata *data_, commoninputdata *cdata_ )
{
    // for checking for possible setting configuration user did not really want
    bool killpixsrc = 0;

    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        data_[g].data_ = data_;

        // set penalty flags to true if true for any source
        std::fill (data_[g].penaltyquery, data_[g].penaltyquery + data_[g].extlengths[5], 0);
        // loop over penalties
        for (PS_SIT j=0; j<data_[g].extlengths[5]; ++j)
            // loop over sources
            for (PS_SIT s=0; s<data_[g].extlengths[13]; ++s)
                data_[g].penaltyquery[j] = data_[g].penaltyquery[j] || (PS_SIT)data_[g].penaltymatrix[s][j][0];

        // figure out if there's a need to de-lens
        data_[g].doreconstruction =
            !data_[g].noevinochi ||
            data_[g].penaltyquery[0] ||
            data_[g].penaltyquery[1] ||
            data_[g].penaltyquery[5] ||
            data_[g].magparams ||
            data_[g].printvec ? 1 : 0;


        if( data_[g].magparams && data_[g].traceparams[2] &&
            data_[g].magparams != ( data_[g].traceparams[2]+1 ) )
        {
	       data_[g].magparams = (PS_SIT)data_[g].traceparams[2]+1;
            /*
              PRINTER printerror("","cannot trace magnification and regularization "
              "strength simultaneously if sampling is different!",
              cdata_->print2screenmutex                            );
            */
        }
        if(data_[g].magparams)
        {
            MEMORY ps_malloc( &(data_[g].magsptr        ), 1 );
            data_[g].mags = new (data_[g].magsptr) outstreamstruct;

            MEMORY ps_malloc( &(data_[g].mags->streamptr), 1 );
            MEMORY ps_malloc( &(data_[g].mags->lockptr  ), 1 );
            data_[g].mags->stream = new (data_[g].mags->streamptr) std::ofstream;
            data_[g].mags->lock   = new (data_[g].mags->lockptr  ) pthread_mutex_t;

            pthread_mutex_init(data_[g].mags->lock,NULL);

            PRINTER openoutstream( cdata_->basename,
                                   data_[g].name, 1,"magnifications.dat",
                                   data_[g].mags->stream, 0               );
            PRINTER writeoutstream( "# regularization_const. magnification",
                                    data_[g].mags->stream, data_[g].mags->lock, data_[g].precision, NULL);
        }
        if(data_[g].maguncertainty)
        {
	  
            MEMORY ps_malloc( &(data_[g].mag_errsptr        ), 1 );
            data_[g].mag_errs = new (data_[g].mag_errsptr) outstreamstruct;

            MEMORY ps_malloc( &(data_[g].mag_errs->streamptr), 1 );
            MEMORY ps_malloc( &(data_[g].mag_errs->lockptr  ), 1 );
            data_[g].mag_errs->stream = new (data_[g].mag_errs->streamptr) std::ofstream;
            data_[g].mag_errs->lock   = new (data_[g].mag_errs->lockptr  ) pthread_mutex_t;

            pthread_mutex_init(data_[g].mag_errs->lock,NULL);

            PRINTER openoutstream( cdata_->basename,
                                   data_[g].name, 1,"magnification_uncertainties.dat",
                                   data_[g].mag_errs->stream, 0               );
            PRINTER writeoutstream( "# 0% 2.5% 5% 16% 30% 40% 50% 60% 70% 84% 95% 97.5% 100% quantiles",
                                    data_[g].mag_errs->stream, data_[g].mag_errs->lock, data_[g].precision, NULL);
            PRINTER writeoutstream( "",
                                    data_[g].mag_errs->stream, data_[g].mag_errs->lock, data_[g].precision, NULL);
        }
        if(data_[g].traceparams[2])
        {
            MEMORY ps_malloc( &(data_[g].tracesptr        ), 1 );
            data_[g].traces = new (data_[g].tracesptr) outstreamstruct;

            MEMORY ps_malloc( &(data_[g].traces->streamptr), 1 );
            MEMORY ps_malloc( &(data_[g].traces->lockptr  ), 1 );
            data_[g].traces->stream = new (data_[g].traces->streamptr) std::ofstream;
            data_[g].traces->lock   = new (data_[g].traces->lockptr  ) pthread_mutex_t;

            pthread_mutex_init(data_[g].traces->lock,NULL);

            PRINTER openoutstream(cdata_->basename, data_[g].name, 1,
                                  "stat.vs.strength.dat",data_[g].traces->stream, 0);
            PRINTER writeoutstream(
                "# regularization_const statistic "
                "detCd_term detA_term detC_term Ed_term Es_term lambda_term pi_term "
                "brightness-weighted_size_penalty mag_penalty axis_ratio_penalty "
                "mismatched_images_penalty conxev_hull_size_penalty flux_outside_data_penalty",
                data_[g].traces->stream, data_[g].traces->lock, data_[g].precision, NULL );
        }
        if(data_[g].printdetails)
        {
            MEMORY ps_malloc( &(data_[g].detailsptr        ), 1 );
            data_[g].details = new (data_[g].detailsptr) outstreamstruct;

            MEMORY ps_malloc( &(data_[g].details->streamptr), 1 );
            MEMORY ps_malloc( &(data_[g].details->lockptr  ), 1 );
            data_[g].details->stream = new (data_[g].details->streamptr) std::ofstream;
            data_[g].details->lock   = new (data_[g].details->lockptr  ) pthread_mutex_t;

            pthread_mutex_init(data_[g].details->lock,NULL);

            PRINTER openoutstream(cdata_->basename, data_[g].name, 1,
                                  "details.dat",data_[g].details->stream, 0);

            PRINTER writeoutstream(
                "# peak src_peak_ra dec(arcs) src_peak_ra dec(absolute) "
                "magnifcation statistic penalty_1 2 3 4 5 6 "
                "Cd A C Ed Es lambda pi chi2",
                data_[g].details->stream, data_[g].details->lock, data_[g].precision, NULL );
        }
        if(data_[g].penaltyquery[5])
        {
            data_[g].findsisterimages = true;
        }

        if( data_[g].debug && !data_[g].printvec )
        {
            data_[g].printvec = 1;
        }

        pthread_mutex_lock (cdata_->wcsmutex);
        MEMORY ps_malloc (&data_[g].oldloc,     data_[g].ndp*2);
        MEMORY ps_malloc (&data_[g].oldloc_wcs, data_[g].ndp*2);
        MEMORY ps_malloc (&data_[g].oldloc_arc, data_[g].ndp*2);
        double pos0[2], pos1[2];
        for(PS_SIT r = 0; r < data_[g].ndp; ++r)
        {
            data_[g].oldloc[r*2]   = r / data_[g].imgy;
            data_[g].oldloc[r*2+1] = r % data_[g].imgy;
            HEADER getimgwcscoord (data_[g].wcs, data_[g].imgy,
                                   data_[g].oldloc[r*2], data_[g].oldloc[r*2+1],
                                   &pos0[0], &pos0[1]);
            data_[g].oldloc_wcs[r*2  ] = pos0[0];
            data_[g].oldloc_wcs[r*2+1] = pos0[1];
            HEADER getwcslfromwcss (pos0[0], pos0[1], data_[g].pra, data_[g].pdec,
                                    data_[g].r1, data_[g].r2, pos1);
            data_[g].oldloc_arc[r*2  ] = pos1[0];
            data_[g].oldloc_arc[r*2+1] = pos1[1];
        }
        pthread_mutex_unlock (cdata_->wcsmutex);

        if( !data_[g].reg || !data_[g].lambdaguess )
        {
            data_[g].reg = 0;
            data_[g].lambdaguess = 0;
        }


        // for analytic sources using vector source, I triangluate the vertices here.
        // technically, the triangulation can change if the vertices are stretches/skewed,
        // but for small stretches/skews, it shouldn't be a big deal
        if (data_[g].num_vec_src && !data_[g].vector_src_triout)
            MEMORY ps_malloc (&(data_[g].vector_src_triout), data_[g].num_vec_src);
        if (data_[g].num_vec_src && !data_[g].rotation_axis)
            MEMORY ps_malloc (&(data_[g].rotation_axis), data_[g].num_vec_src*2);
        if (data_[g].usersetsrc)
        {
            PS_SIT vecindex, numpoints;
            for (PS_SIT s=0; s<(PS_SIT)data_[g].usersetsrc[0]; ++s)
            {
                // if this is a vector source, convert coordinates from WCS to pixel coordinates
                // and triangulate them
                if ((PS_SIT)data_[g].usersetsrc[1+s*17]>1000)
                {
                    vecindex = (PS_SIT)data_[g].usersetsrc[1+s*17]-1001;
                    numpoints = data_[g].vector_src_numpix[vecindex];
                    // convert input pixel positions to image coordinates
                    double *vertices_orig;
                    MEMORY ps_malloc( &(vertices_orig), numpoints*2 );
                    double pos3[2];
                    for (PS_SIT s=0; s<numpoints; ++s)
                    {
                        HEADER getwcssfromwcsl (data_[g].vector_src_pos[vecindex][s*2],
                                                data_[g].vector_src_pos[vecindex][s*2+1],
                                                data_[g].pra, data_[g].pdec,
                                                data_[g].r1, data_[g].r2, pos3);
                        data_[g].vector_src_pos[vecindex][s*2]   = pos3[0];
                        data_[g].vector_src_pos[vecindex][s*2+1] = pos3[1];
                        HEADER getimgpixcoord (data_[g].wcs, data_[g].imgy, cdata_, pos3[0], pos3[1],
                                               &vertices_orig[s*2], &vertices_orig[s*2+1]);
                    }

                    // create triangulation structures
                    struct triangulateio *triinloc, *trioutloc;
                    MEMORY ps_malloc(            &triinloc , 1 );
                    MEMORY ps_malloc(            &trioutloc, 1 );
                    MEMORY triangulatestructinit( triinloc     );
                    MEMORY triangulatestructinit( trioutloc    );
                    triinloc->numberofpointattributes = 0;
                    triinloc->numberofpoints = numpoints;
                    triinloc->pointlist = vertices_orig;

                    // triangulate points
                    ps_tri_triangulate ((char*)CONSTANT triswitchesnominangle,
                                        triinloc, trioutloc, (struct triangulateio*)NULL);
                    trioutloc->numberofpoints = numpoints;
                    trioutloc->pointlist = triinloc->pointlist;
                    triinloc->pointlist = 0;

                    // destroy input triangle struct
                    MEMORY triangulatestructdestruct( triinloc  );
                    MEMORY ps_free(                   triinloc  );
                    // store triangulation permanently
                    data_[g].vector_src_triout[vecindex] = trioutloc;

                    // now, find center of source to determine the rotation
                    // axis for linear transformations of the pixel coordinates
                    double wt, norm=0, sctr[2] = {0,0};
                    for (PS_SIT s=0; s<numpoints; ++s)
                    {
                        wt = data_[g].vector_src_flux[vecindex][s];
                        wt *= wt;
                        wt *= wt;
                        sctr[0] += vertices_orig[s*2]   *wt;
                        sctr[1] += vertices_orig[s*2+1] *wt;
                        norm += wt;
                    }
                    data_[g].rotation_axis[vecindex*2]   = sctr[0] /norm;
                    data_[g].rotation_axis[vecindex*2+1] = sctr[1] /norm;

                    // now, shift pixel coordinates of triangulation so that
                    // the centroid of the convex hull is at the origin.
                    for (PS_SIT s=0; s<numpoints; ++s)
                    {
                        data_[g].vector_src_triout[vecindex]->pointlist[s*2  ] -=
                            data_[g].rotation_axis[vecindex*2];
                        data_[g].vector_src_triout[vecindex]->pointlist[s*2+1] -=
                            data_[g].rotation_axis[vecindex*2+1];
                    }
                }
            }
        }
	
	// load user-specified source parameters
	if (data_[g].hacksrcfn)
	{
	    // read file
	    char **vec=NULL, **vecsplit=NULL;
	    PS_SIT vecsize=0,  vecsplitsize=0;
	    string fn = string(CONSTANT dir_in) + string(data_[g].hacksrcfn);
	    OPERA readfile (fn.c_str(), &vec, &vecsize, cdata_->print2screenmutex);
	    
	    // parse data
	    MEMORY ps_malloc (&data_[g].hacksrc, vecsize);
	    for (PS_SIT j=0; j<vecsize; ++j)
	    {
		MEMORY ps_free (vecsplit, vecsplitsize);
		vecsplit = NULL;
		OPERA split (vec[j], " ", &vecsplit, &vecsplitsize);
		
		if (1!=vecsplitsize)
		    PRINTER printerror("", "error in \"" +
				       fn + "\" file", cdata_->print2screenmutex);

		data_[g].hacksrc[j] = OPERA convert_string <double> (vecsplit[0]);
	    }	    
	    MEMORY ps_free (vec,      vecsize);
	    MEMORY ps_free (vecsplit, vecsplitsize);
	}

	
#include "pixsrc_init_safeguard.cpp"
	
	
    }
    
    if (killpixsrc)
	PRINTER printerror("pixsrc", 
			   "Found settings that are uncommon (but possibly okay) or will lead to problems.\n"
			   "Exiting now.. To allow these settings, set the \"fatalwarn:\" parameter to \"0\".", 
			   cdata_->print2screenmutex);
    
    
#ifdef PIXSRC_MEMORYDEBUG
    MEMORY nummallocs = 0;
    MEMORY numfrees   = 0;
    pthread_mutex_init(& MEMORY  freemutex ,NULL);
    pthread_mutex_init(& MEMORY mallocmutex,NULL);
#endif
    
}

void pixsrc_init::setup_one_libwcs (commoninputdata *cdata_,
				    string fn, ps_WorldCoor **wcs, char *coordsysfinal,
				    double *px, double *py, double *pra, double *pdec, PS_SIT imgy)
{
    char b1950[] = "B1950";
    char j2000[] = "J2000";
    strcpy(coordsysfinal,cdata_->coordsysorig);
    pthread_mutex_lock( cdata_->wcsmutex );
    PS_SIT pixelsys = 0;
    if(!strcmp(cdata_->coordsysorig,"PIXEL"))
	pixelsys=1;
    
    char fnc[fn.size()+1];
    strcpy (fnc, fn.c_str());
    char *header = EXTERNAL ps_GetFITShead (fnc, 0);
    double cra, cdec, dra, ddec, secpix, eqout=0.0;
    PS_SIT wp, hp, sysout=0;
    
    if(!strcmp(cdata_->coordsysorig,"J2000"))
    {
	EXTERNAL ps_setsys(EXTERNAL ps_WCS_J2000);
	sysout = EXTERNAL ps_WCS_J2000;
	eqout = 2000.0;
    }
    else if(!strcmp(cdata_->coordsysorig,"B1950"))
    {
	EXTERNAL ps_setsys(EXTERNAL ps_WCS_B1950);
	sysout = EXTERNAL ps_WCS_B1950;
	eqout = 1950.0;
    }
    else if(!strcmp(cdata_->coordsysorig,"ecliptic"))
    {
	EXTERNAL ps_setsys(EXTERNAL ps_WCS_ECLIPTIC);
	sysout = EXTERNAL ps_WCS_ECLIPTIC;
	eqout = 2000.0;
    }
    else if(!strcmp(cdata_->coordsysorig,"galactic"))
    {
	EXTERNAL ps_setsys(EXTERNAL ps_WCS_GALACTIC);
	sysout = EXTERNAL ps_WCS_GALACTIC;
	eqout = 2000.0;
    }
    
    *wcs = EXTERNAL ps_GetFITSWCS ( fnc, header, 0, &cra, &cdec,&dra, &ddec,
				    &secpix,&wp, &hp, &sysout, &eqout );
    if( EXTERNAL ps_nowcs(*wcs) )
	PRINTER printerror("","no WCS in file " + fn, cdata_->print2screenmutex);
    
    MEMORY ps_free( header );
    
    if((*wcs)->sysout == EXTERNAL ps_WCS_B1950)
	EXTERNAL ps_wcsoutinit (*wcs, b1950);
    if((*wcs)->sysout == EXTERNAL ps_WCS_J2000)
	EXTERNAL ps_wcsoutinit (*wcs, j2000);
    
    if(pixelsys)
    {
	strcpy(coordsysfinal,(*wcs)->radecin);
    }
    else
    {
	EXTERNAL ps_wcscsys(cdata_->coordsysorig);
	EXTERNAL ps_wcsoutinit (*wcs, cdata_->coordsysorig);
    }
    
    (*wcs)->degout = 1;
    (*wcs)->ndec = 10;
    
    if( pixelsys )
    {
	*px -= 1;
	*py = imgy-*py;
	HEADER getimgwcscoord (*wcs, imgy, *px, *py, pra, pdec);
    }
    else
    {
	HEADER getimgpixcoord (*wcs, imgy, cdata_, *pra, *pdec, px, py);
    }
    pthread_mutex_unlock( cdata_->wcsmutex );
}

void pixsrc_init::setuplibwcs(char *bn, char **name, char **namewithext,
                              inputdata *data_, commoninputdata *cdata_)
{
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
	string fn = string(CONSTANT dir_in) + string(namewithext[g]);
	INIT setup_one_libwcs (cdata_, fn, &data_[g].wcs, cdata_->coordsysfinal,
			       &data_[g].px, &data_[g].py, &data_[g].pra, &data_[g].pdec, data_[g].imgy);
    }
}

void pixsrc_init::settingup(char * bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_, PS_SIT *dimil, char ***ignorelist, char ***ignorelistwithext)
{
    // remove existing files or create directory
    DIR *dir = NULL;
    struct dirent *ent;
    dir = opendir ( CONSTANT dir_out );
    if (dir != NULL)
    {
        while ((ent = readdir (dir)) != NULL)
        {
            string namefile(ent->d_name);
            string bnstring(bn);
            string comparetobn = bnstring+"_";
            if(!namefile.substr(0,bnstring.size()+1).compare(comparetobn))
            {
                remove( ( string(CONSTANT dir_out) + namefile ).c_str() );
            }
        }
        closedir (dir);
    }
    else
    {
        for(PS_SIT g=0; g<cdata_->numimages; g++)
            if( 1 || data_[g].debug || data_[g].printvec || data_[g].printdetails )
            {
                // make this more portable in free time!
                string cmd = "mkdir -p " + string(CONSTANT dir_out);
                if (system(cmd.c_str()))
                    PRINTER printerror("pixsrc","failed to run: " + cmd,
                                       cdata_->print2screenmutex);
                break;
            }
    }

    // first lets do some image-independent setup

    // initialize GSL random number generator
    gsl_rng_env_setup();
    cdata_->ps_gsl_ran_T = gsl_rng_default;
    cdata_->ps_gsl_ran_r = gsl_rng_alloc (cdata_->ps_gsl_ran_T);

    PS_SIT extlengthextention = 20;
    cdata_->numcall=0;
    MEMORY ps_malloc( &(cdata_->basename), OPERA sizestring(bn)+1 );
    strcpy( cdata_->basename, bn );

    MEMORY ps_malloc( &(cdata_->print2screenmutex), 1    );
    MEMORY ps_malloc( &(cdata_->potdefmagmutex   ), 1    );
    MEMORY ps_malloc( &(cdata_->wcsmutex         ), 1    );
    MEMORY ps_malloc( &(cdata_->fitsiomutex      ), 1    );
    pthread_mutex_init( cdata_->print2screenmutex , NULL );
    pthread_mutex_init( cdata_->potdefmagmutex    , NULL );
    pthread_mutex_init( cdata_->wcsmutex          , NULL );
    pthread_mutex_init( cdata_->fitsiomutex       , NULL );

    MEMORY ps_malloc( &(cdata_->attrjoinable), 1 );
    MEMORY ps_malloc( &(cdata_->attrdetached), 1 );
    pthread_attr_init(  cdata_->attrjoinable     );
    pthread_attr_init(  cdata_->attrdetached     );
    pthread_attr_setdetachstate(cdata_->attrjoinable, PTHREAD_CREATE_JOINABLE);
    pthread_attr_setdetachstate(cdata_->attrdetached, PTHREAD_CREATE_DETACHED);

    // get names of images
    char *fname;
    const char *listcc[3];
    listcc[0] = CONSTANT dir_in;
    listcc[1] = bn;
    listcc[2] = CONSTANT imagelistfile;
    OPERA concatenate( listcc, 3, &fname );

    char **vec, ***vecsplit;
    PS_SIT vecsize,  *vecsplitsize;

    OPERA readfile( fname, &vec, &vecsize, cdata_->print2screenmutex);
    MEMORY ps_free( fname );

    // these will contain parsed lines of vec
    MEMORY ps_malloc( &vecsplit    , vecsize );
    MEMORY ps_malloc( &vecsplitsize, vecsize );

    // find number of images and those to ignore and do error checking
    PS_SIT numignore = 0;
    PS_SIT numgood   = 0;
    for( PS_SIT j=0; j<vecsize; ++j )
    {
        OPERA split( vec[j], " ", &(vecsplit[j]), &(vecsplitsize[j]) );

        if( !strcmp(vecsplit[j][0], "-ignore") )
        {
            if( vecsplitsize[j] != 2 )
            {
                PRINTER printerror("", "\"-ignore\" command used improperly in \"" +
                                   string(CONSTANT imagelistfile) +
                                   "\" file: \"-ignore\" takes only one argument.\n"
                                   "Use multiple -ignore commands to ignore multiple images.",
                                   cdata_->print2screenmutex);
            }
            ++numignore;
            // signals ignore image
            vecsplitsize[j] = 0;
        }
        else
        {
            if( vecsplitsize[j] != 1 )
            {
                PRINTER printerror("", "error in \"" +
                                   string(CONSTANT imagelistfile) + "\" file", cdata_->print2screenmutex);
            }
            ++numgood;
            // signals include image
            vecsplitsize[j] = 1;
        }
    }

    PS_SIT numimages = ( !numgood ) ? 1 : numgood;
    cdata_->numimages = numimages;

    *dimil = numignore;
    MEMORY ps_malloc( &(*ignorelist)       , *dimil );
    MEMORY ps_malloc( &(*ignorelistwithext), *dimil );

    PS_SIT badindex  = 0;
    PS_SIT goodindex = 0;
    PS_SIT sizevec;
    PS_SIT sizeconst = OPERA sizestring(CONSTANT fitsfile);
    for( PS_SIT jj=0; jj<vecsize; ++jj )
    {
        // if image good
        if( vecsplitsize[jj] == 1 )
        {
            sizevec = OPERA sizestring( vecsplit[jj][0] );
            MEMORY ps_malloc( &(namewithext[goodindex]), sizevec+1 );
            strcpy( namewithext[goodindex], vecsplit[jj][0] );

            if( sizevec >= sizeconst )
            {
                // copy into test the end of the file to ignore
                char test[sizevec+1];
                strcpy( test, vecsplit[jj][0]+sizevec-sizeconst );

                std::transform(test,test+sizeconst+1,test,::tolower);

                // if last part of filename is CONSTANT fitsfile, then chop it off
                if( !strcmp(test, CONSTANT fitsfile) )
                    vecsplit[jj][0][sizevec-sizeconst] = 0;
            }

            sizevec = OPERA sizestring( vecsplit[jj][0] );
            MEMORY ps_malloc( &(name[goodindex]), sizevec+1 );
            strcpy( name[goodindex], vecsplit[jj][0] );

            ++goodindex;
        }
        else
        {
            sizevec = OPERA sizestring(vecsplit[jj][1]);
            MEMORY ps_malloc( &((*ignorelistwithext)[badindex]), sizevec+1 );
            strcpy( (*ignorelistwithext)[badindex], vecsplit[jj][1] );

            if( sizevec >= sizeconst )
            {
                // copy into test the end of the file to ignore
                char test[sizevec+1];
                strcpy( test, vecsplit[jj][1]+sizevec-sizeconst );

                std::transform(test,test+sizeconst+1,test,::tolower);

                // if last part of filename is CONSTANT fitsfile, then chop it off
                if( !strcmp(test, CONSTANT fitsfile) )
                    vecsplit[jj][1][sizevec-sizeconst] = 0;
            }

            sizevec = OPERA sizestring( vecsplit[jj][1] );
            MEMORY ps_malloc( &((*ignorelist)[badindex]), sizevec+1 );
            strcpy( (*ignorelist)[badindex], vecsplit[jj][1] );

            ++badindex;
        }
    }

    // free data that was read in and split
    for( PS_SIT j=0; j<vecsize; ++j )
        MEMORY ps_free( vecsplit[j], vecsplitsize[j] );
    MEMORY ps_free( vecsplit          );
    MEMORY ps_free( vec, vecsize      );
    MEMORY ps_free( vecsplitsize      );

    // if no images listed to be included use basename as imagename
    if( !numgood )
    {
        MEMORY ps_malloc( &(name[0]), OPERA sizestring(bn)+1 );
        strcpy( name[0], bn );
        const char *listcc[2];
        listcc[0] = bn;
        listcc[1] = CONSTANT fitsfile;
        OPERA concatenate( listcc, 2, &(namewithext[0]) );
    }

    PS_SIT minnamelength = OPERA sizestring("pixsrc");
    PS_SIT sizename;
    for( PS_SIT g=0; g<numimages; ++g )
    {
        sizename = OPERA sizestring(name[g]);
        MEMORY ps_malloc( &(data_[g].name), sizename+1 );
        strcpy( data_[g].name, name[g] );

        MEMORY ps_malloc( &(data_[g].extlengths), extlengthextention );
        std::fill( data_[g].extlengths, data_[g].extlengths+extlengthextention, 0 );

        if( sizename > minnamelength )
            minnamelength = sizename;
    }
    for( PS_SIT g=0; g<numimages; ++g )
    {
        PS_SIT size = OPERA sizestring(name[g]);
        MEMORY ps_malloc( &(data_[g].print2screenname), minnamelength+1 );
        for( PS_SIT let=0; let<size; ++let )
            data_[g].print2screenname[let] = name[g][let];
        for( PS_SIT let=size; let<minnamelength; ++let )
            data_[g].print2screenname[let] = ' ';
        data_[g].print2screenname[minnamelength] = 0;
    }

    for( PS_SIT g=0; g<numimages; ++g )
        for( PS_SIT g2=g+1; g2<numimages; ++g2 )
        {
            string test1(name[g]);
            std::transform(test1.begin(),test1.end(),test1.begin(),::tolower);
            string test2(name[g2]);
            std::transform(test2.begin(),test2.end(),test2.begin(),::tolower);
            if( !test1.compare(test2) )
                PRINTER printerror("","multiple files found with similar or same name: \"" + string(name[g]) + "\"" , cdata_->print2screenmutex);
        }
}

void pixsrc_init::read_one_image(inputdata *data_, commoninputdata *cdata_, string fn,
				 PS_SIT *wcsinfolength,    double **wcsinfo_,
				 PS_SIT *invwcsinfolength, double **invwcsinfo_,
				 PS_SIT *imgdatalength,    double **imgdata_,
				 PS_SIT *imgx, PS_SIT *imgy, double *rollangle,
				 double *pix2arc, double *arc2pix, PS_SIT *ndp)
{
    *wcsinfolength = 10;
    MEMORY ps_malloc (&(*wcsinfo_), *wcsinfolength);
    *invwcsinfolength = 4;
    MEMORY ps_malloc (&(*invwcsinfo_), *invwcsinfolength);
    double *wcsinfo    = *wcsinfo_;
    double *invwcsinfo = *invwcsinfo_;
    
    fitsfile *fptr;
    
    std::fill(wcsinfo, wcsinfo+*wcsinfolength, 0);
    char card[FLEN_CARD];
    int/*PS_SIT*/ status = 0, nkeys;
    
    fits_open_file(&fptr, fn.c_str(), READONLY, &status);
    fits_get_hdrspace(fptr, &nkeys, NULL, &status);
    for(PS_SIT ii=1; ii<=nkeys; ii++)
    {
	fits_read_record(fptr, ii, card, &status);
	string temp(card);
	
	if(temp.find("CRPIX1") == 0)
	    wcsinfo[0] = OPERA getfitscardvalue(temp);
	else if(temp.find("CRPIX2") == 0)
	    wcsinfo[1] = OPERA getfitscardvalue(temp);
	else if(temp.find("CRVAL1") == 0)
	    wcsinfo[2] = OPERA getfitscardvalue(temp);
	else if(temp.find("CRVAL2") == 0)
	    wcsinfo[3] = OPERA getfitscardvalue(temp);
	else if(temp.find("CD1_1") == 0)
	    wcsinfo[4] = OPERA getfitscardvalue(temp);
	else if(temp.find("CD1_2") == 0)
	    wcsinfo[5] = OPERA getfitscardvalue(temp);
	else if(temp.find("CD2_1") == 0)
	    wcsinfo[6] = OPERA getfitscardvalue(temp);
	else if(temp.find("CD2_2") == 0)
	    wcsinfo[7] = OPERA getfitscardvalue(temp);
	else if(temp.find("CTYPE1") == 0 && !OPERA equalsnpos(temp.find("RA")))
	    wcsinfo[8] = 1;
	else if(temp.find("NAXIS1") == 0)
	    *imgx = (PS_SIT)OPERA getfitscardvalue(temp);
	else if(temp.find("NAXIS2") == 0)
	    *imgy = (PS_SIT)OPERA getfitscardvalue(temp);
    }
    
    // if rotation matrix keywords not found
    if(wcsinfo[4] == 0 && wcsinfo[5] == 0 &&
       wcsinfo[6] == 0 && wcsinfo[7] == 0)
	PRINTER printerror(data_->print2screenname,
			   "CD?_? keywords not found in FITS header",
			   cdata_->print2screenmutex);
    
    // if coordinate systems is right-handed (positive determinant)
    if (wcsinfo[4]*wcsinfo[7] -
	wcsinfo[5]*wcsinfo[6] > 0)
	PRINTER printerror(data_->print2screenname,
			   "right-handed coordinate systems may cause problems",
			   cdata_->print2screenmutex);
    
    // angle from north to y-axis, east positive
    *rollangle = -std::atan2 (-wcsinfo[5], wcsinfo[7]);
    
    for(PS_SIT ii=1; ii<=nkeys; ii++)
    {
	fits_read_record(fptr, ii, card, &status);
	string temp(card);
	if(temp.find("NAXIS3") == 0)
	{
	    double naxis3 = OPERA getfitscardvalue(temp);
	    if(naxis3!=0)
		PRINTER printerror(data_->print2screenname,
				   "multidimensional fits files currently not supported, even if NAXIS3=1."
				   "\nNAXIS3 must be removed. Try re-saving the file. Sorry!", cdata_->print2screenmutex);
	}
    }
    
    double factor=wcsinfo[4]*wcsinfo[7]-wcsinfo[5]*wcsinfo[6];
    invwcsinfo[0]=wcsinfo[7]/factor;
    invwcsinfo[1]=-wcsinfo[5]/factor;
    invwcsinfo[2]=-wcsinfo[6]/factor;
    invwcsinfo[3]=wcsinfo[4]/factor;
    
    *arc2pix = (
	std::sqrt( invwcsinfo[0] * invwcsinfo[0] +
		   invwcsinfo[1] * invwcsinfo[1] ) +
	std::sqrt( invwcsinfo[2] * invwcsinfo[2] +
		   invwcsinfo[3] * invwcsinfo[3] ) ) / 2.0 / 3600.0;
    *pix2arc = 1.0 / *arc2pix;
    
    *ndp = *imgx**imgy;
    *imgdatalength = *ndp;
    MEMORY ps_malloc (&(*imgdata_), *imgdatalength);
    double *imgdata = *imgdata_;
    
    int/*PS_SIT*/ naxis, bitpix;
    long naxes[2] = {1,1}, fpixel[2] = {1,1};
    fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);
    double *pixels;
    MEMORY ps_malloc( &(pixels), naxes[0] );
    for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
    {
	fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status);
	
	for(PS_SIT ii=0; ii<naxes[0]; ii++)
	    imgdata[*imgy-fpixel[1] + ii**imgy] = pixels[ii];
    }
    
    MEMORY ps_free( pixels );
    fits_close_file(fptr, &status);    
}

void pixsrc_init::readimages(char *bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_)
{
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
	string filename (CONSTANT dir_in + string(namewithext[g]));
	INIT read_one_image (&data_[g], cdata_, filename,
			     &data_[g].extlengths[15], &data_[g].wcsinfo,
			     &data_[g].extlengths[7],  &data_[g].invwcsinfo,
			     &data_[g].extlengths[8],  &data_[g].data,
			     &data_[g].imgx, &data_[g].imgy, &data_[g].rollangle,
			     &data_[g].pix2arc, &data_[g].arc2pix, &data_[g].ndp);
    }
}

void pixsrc_init::readmmmasks(char *bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_)
{
    // read image masks for quick rejection of lens models
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        // get names of files
        PS_SIT sizeendfile = OPERA sizestring( CONSTANT imagesfileendreg );
        PS_SIT numfiles = 0;
        char **mmimagemaskfn;

        char *namefilecompare;
        const char *listcc[4];
        listcc[0] = bn;
        listcc[1] = CONSTANT bnseparator;
        listcc[2] = name[g];
        listcc[3] = CONSTANT imagesfilestart;
        OPERA concatenate( listcc, 4, &namefilecompare );

        DIR *dir;
        struct dirent *ent;
        dir = opendir ( CONSTANT dir_in );
        if( dir != NULL )
        {
            while ((ent = readdir (dir)) != NULL)
            {
                if( strstr( ent->d_name, namefilecompare ) == ent->d_name &&
                    !strcmp( ent->d_name+OPERA sizestring(ent->d_name)-sizeendfile,
                             CONSTANT imagesfileendreg ) )
                {
                    ++numfiles;
                }
            }
            closedir( dir );
        }

        data_[g].extlengths[13] = 0==numfiles ? 1 : numfiles;
        data_[g].extlengths[5] = 6;
        data_[g].extlengths[6] = 6;
        MEMORY ps_malloc (&data_[g].penaltynames, numfiles);
        MEMORY ps_malloc (&data_[g].penaltymatrix,
                          data_[g].extlengths[13], data_[g].extlengths[5], data_[g].extlengths[6]);
        for (PS_SIT s=0; s<data_[g].extlengths[13]; ++s)
        {
            for(PS_SIT j=0; j<data_[g].extlengths[5]; j++)
            {
                std::fill( data_[g].penaltymatrix[s][j],
                           data_[g].penaltymatrix[s][j] + data_[g].extlengths[6]-1, 0 );

                OPERA assign_p_infinity (&data_[g].penaltymatrix[s][j][data_[g].extlengths[6]-1]);
                data_[g].penaltymatrix[s][j][data_[g].extlengths[6]-1] =
                    std::sqrt (data_[g].penaltymatrix[s][j][data_[g].extlengths[6]-1]/1e3);
            }
        }

        MEMORY ps_malloc (&mmimagemaskfn, numfiles);
        PS_SIT fileindex = 0;

        dir = opendir ( CONSTANT dir_in );
        if( dir != NULL )
        {
            while ((ent = readdir (dir)) != NULL)
            {
                if( strstr( ent->d_name, namefilecompare ) == ent->d_name &&
                    !strcmp( ent->d_name+OPERA sizestring(ent->d_name)-sizeendfile,
                             CONSTANT imagesfileendreg ) )
                {
                    // copy name of mmimage file into penaltynames
                    // copy filename
                    MEMORY ps_malloc (&data_[g].penaltynames[fileindex],
                                      OPERA sizestring (ent->d_name)+1);
                    PS_SIT sind = 0;
                    char *ptr = ent->d_name+OPERA sizestring (namefilecompare);
                    while (ptr!=ent->d_name+OPERA sizestring(ent->d_name)-sizeendfile)
                    {
                        data_[g].penaltynames[fileindex][sind++] = *ptr;
                        ++ptr;
                    }
                    data_[g].penaltynames[fileindex][sind] = 0;
                    listcc[0] = CONSTANT dir_in;
                    listcc[1] = ent->d_name;
                    OPERA concatenate( listcc, 2, &(mmimagemaskfn[fileindex]) );
                    ++fileindex;
                }
            }
            closedir( dir );
        }

        data_[g].penaltyfilenames = mmimagemaskfn;
        data_[g].extlengths[11] = numfiles;
        data_[g].extlengths[12] = data_[g].ndp;

        MEMORY ps_free( namefilecompare );
    }
}
void pixsrc_init::readmasks(char *bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_)
{
    // namewithext is unused and is NULL sometimes

    for (PS_SIT g=0; g<cdata_->numimages; ++g)
    {
        data_[g].extlengths[16] = data_[g].ndp;
        if (data_[g].imagemasks)
            MEMORY ps_free (data_[g].imagemasks);
        MEMORY ps_malloc (&data_[g].imagemasks, data_[g].extlengths[16]);
        std::fill (data_[g].imagemasks,data_[g].imagemasks+data_[g].extlengths[16],1);

        if (data_[g].magmasks)
            MEMORY ps_free (data_[g].magmasks);
        MEMORY ps_malloc (&data_[g].magmasks, data_[g].extlengths[16]);
        std::fill (data_[g].magmasks,data_[g].magmasks+data_[g].extlengths[16],-1);
    }

    // read chi2 mask
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        PS_SIT dim1,dim2;
        double **polys = NULL;

        char *file;
        const char *listcc[5];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = bn;
        listcc[2] = CONSTANT bnseparator;
        listcc[3] = name[g];
        listcc[4] = CONSTANT chi2maskfilereg;
        OPERA concatenate( listcc, 5, &file );

        readmaskpolygon( &dim1,&dim2,&polys, &file, &data_[g], cdata_);
        MEMORY ps_free( file );

        if(dim1)
        {
            data_[g].extlengths[19] = data_[g].extlengths[8];
            if (data_[g].chi2mask)
                MEMORY ps_free (data_[g].chi2mask);
            MEMORY ps_malloc( &(data_[g].chi2mask), data_[g].extlengths[19] );
            std::fill(data_[g].chi2mask,data_[g].chi2mask+data_[g].extlengths[19],0);
        }

        // find masked pixels
        for(PS_SIT r=0; r<data_[g].extlengths[19]; r++)
        {
            PS_SIT x=r/data_[g].imgy;
            PS_SIT y=r%data_[g].imgy;
            for(PS_SIT s=0; s< dim1; s++)
            {
                if(GEOM isinpoly(x,y,polys[s],dim2/2))
                {
                    data_[g].chi2mask[r]=1;
                    break;
                }
            }
        }

        MEMORY ps_free( polys, dim1 );
    }

    // read bad pixel mask
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {

        // read reg file
        PS_SIT dim1,dim2;
        double **polys = NULL;

        char *file;
        const char *listcc[5];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = bn;
        listcc[2] = CONSTANT bnseparator;
        listcc[3] = name[g];
        listcc[4] = CONSTANT badpixelmaskfilereg;
        OPERA concatenate( listcc, 5, &file );

        readmaskpolygon( &dim1,&dim2,&polys, &file, &data_[g], cdata_);
        MEMORY ps_free( file );

        if(dim1)
            data_[g].useall = 0;

        // find bad pixels
        for(PS_SIT r=0; r<data_[g].extlengths[16]; r++)
        {
            // if value for this pixel is NaN, it's a bad pixel
            if (data_[g].data[r]!=data_[g].data[r])
            {
                data_[g].imagemasks[r]=0;
                data_[g].useall = 0;
                continue;
            }
            for(PS_SIT s=0; s< dim1; s++)
            {
                PS_SIT x=r/data_[g].imgy;
                PS_SIT y=r%data_[g].imgy;
                if(GEOM isinpoly(x,y,polys[s],dim2/2))
                {
                    data_[g].imagemasks[r]=0;
                    break;
                }
            }
        }

        MEMORY ps_free( polys, dim1 );

        // read ascii file
        file = NULL;
        listcc[4] = CONSTANT badpixelmaskfileascii;
        OPERA concatenate( listcc, 5, &file );

        PS_SIT *list=NULL;
        dim1=0;
        readmaskascii(&dim1,&list,&file,&(data_[g]),cdata_);
        MEMORY ps_free( file );

        if(dim1)
            data_[g].useall=false;
        for(PS_SIT j=0; j<dim1; j++)
        {
            PS_SIT rpos = list[2*j]*data_[g].imgy+list[2*j+1];
            data_[g].imagemasks[rpos] = 0;
        }

        MEMORY ps_free( list );

        // check for fill bad pixel option
        if(data_[g].fillbadpix)
        {
            for(PS_SIT r=0; r<data_[g].ndp; r++)
            {
                double sigma = std::sqrt( data_[g].myvariance0 );
                if(!data_[g].imagemasks[r])
                {
                    if( !data_[g].fillbadpixnoise )
                        data_[g].data[r] = *(data_[g].fillbadpix);
                    else
                        data_[g].data[r] = sigma*OPERA randomgaussian(cdata_->ps_gsl_ran_r);
                }
            }

            data_[g].useall=true;
            std::fill(data_[g].imagemasks,data_[g].imagemasks+data_[g].ndp,1);
        }
    }

    // read data mask
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        char datamaskispresent = 0;

        PS_SIT dim1,dim2;
        double **polys = NULL;

        char *file;
        const char *listcc[5];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = bn;
        listcc[2] = CONSTANT bnseparator;
        listcc[3] = name[g];
        listcc[4] = CONSTANT imagemaskfilereg;
        OPERA concatenate( listcc, 5, &file );

        readmaskpolygon( &dim1,&dim2,&polys, &file, &data_[g], cdata_);
        MEMORY ps_free( file );

        if(dim1)
        {
            data_[g].useall=false;
            datamaskispresent = 1;
        }

        // find masked pixels
        for(PS_SIT r=0; r<data_[g].extlengths[16]; r++)
        {
            if(data_[g].imagemasks[r]==1)
            {
                PS_SIT x=r/data_[g].imgy;
                PS_SIT y=r%data_[g].imgy;
                for(PS_SIT s=0; s< dim1; s++)
                {
                    if(GEOM isinpoly(x,y,polys[s],dim2/2))
                    {
                        data_[g].imagemasks[r]=2;
                        break;
                    }
                }
            }
        }

        MEMORY ps_free( polys, dim1 );

        if(!datamaskispresent)
        {
            for(PS_SIT r=0; r<data_[g].extlengths[16]; r++)
                if(data_[g].imagemasks[r]==1)
                    data_[g].imagemasks[r] = 2;
        }
    }

    // read magnification bootstrap mask
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        PS_SIT dim1,dim2;
        double **polys = NULL;

        char *file;
        const char *listcc[5];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = bn;
        listcc[2] = CONSTANT bnseparator;
        listcc[3] = name[g];
        listcc[4] = CONSTANT magmaskfilereg;
        OPERA concatenate( listcc, 5, &file );

        readmaskpolygon( &dim1,&dim2,&polys, &file, &data_[g], cdata_);
        MEMORY ps_free( file );

        // find masked pixels
        for(PS_SIT r=0; r<data_[g].extlengths[16]; r++)
        {
            if(data_[g].magmasks[r]== -1)
            {
                PS_SIT x=r/data_[g].imgy;
                PS_SIT y=r%data_[g].imgy;
                for(PS_SIT s=0; s< dim1; s++)
                {
                    if(GEOM isinpoly(x,y,polys[s],dim2/2))
                    {
                        data_[g].magmasks[r]=s;
                        break;
                    }
                }
            }
        }

        MEMORY ps_free( polys, dim1 );
    }

    // read mismatched image masks
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        // correct the number of data points
        data_[g].extlengths[12] = data_[g].ndp;

        if (data_[g].num_mmimages)
            MEMORY ps_free (data_[g].num_mmimages);
        if (data_[g].mmborder)
            MEMORY ps_free (data_[g].mmborder);
        if (data_[g].mmborder_defl)
            MEMORY ps_free (data_[g].mmborder_defl);
        if (data_[g].mmimages)
            MEMORY ps_free (data_[g].mmimages, data_[g].extlengths[11]);
        MEMORY ps_malloc (&(data_[g].num_mmimages),  data_[g].extlengths[11]);
        MEMORY ps_malloc (&(data_[g].mmborder),      data_[g].extlengths[11]);
        MEMORY ps_malloc (&(data_[g].mmborder_defl), data_[g].extlengths[11]);
        MEMORY ps_malloc (&(data_[g].mmimages),      data_[g].extlengths[11],
                          data_[g].extlengths[12]);
        for(PS_SIT src=0; src<data_[g].extlengths[11]; src++)
            std::fill (data_[g].mmimages[src],
                       data_[g].mmimages[src]+data_[g].extlengths[12], 0);

        for(PS_SIT src=0; src<data_[g].extlengths[11]; src++)
        {
            PS_SIT dim1,dim2;
            double **polys = NULL;

            readmaskpolygon( &dim1,&dim2,&polys, &(data_[g].penaltyfilenames[src]), &data_[g], cdata_);
            data_[g].num_mmimages[src]=dim1;

            INIT create_mm_border (&polys, dim1, dim2/2, src, g, data_, cdata_);

            // find masked pixels
            for(PS_SIT r=0; r<data_[g].extlengths[12]; r++)
            {
                PS_SIT x=r/data_[g].imgy;
                PS_SIT y=r%data_[g].imgy;
                for(PS_SIT s=0; s< dim1; s++)
                {
                    if(GEOM isinpoly(x,y,polys[s],dim2/2))
                    {
                        data_[g].mmimages[src][r]=s+1;
                        break;
                    }
                }
            }

            MEMORY ps_free( polys, dim1 );
        }
    }

    // read shapelet integration mask
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        PS_SIT dim1,dim2;
        double **polys = NULL;

        char *file;
        const char *listcc[5];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = bn;
        listcc[2] = CONSTANT bnseparator;
        listcc[3] = name[g];
        listcc[4] = CONSTANT shapeintmaskfilereg;
        OPERA concatenate( listcc, 5, &file );

        readmaskpolygon( &dim1,&dim2,&polys, &file, &data_[g], cdata_);
        MEMORY ps_free( file );

        if(dim1)
        {
            data_[g].shapeintmask = polys;
            data_[g].extlengths[17] = dim1;
            data_[g].extlengths[18] = dim2;

            // find mins and maxs for each mask
            if (data_[g].shapeintmaskminmax)
                MEMORY ps_free (data_[g].shapeintmaskminmax, data_[g].extlengths[17]);
            MEMORY ps_malloc (&data_[g].shapeintmaskminmax, data_[g].extlengths[17], 4);
            for (PS_SIT m=0; m<data_[g].extlengths[17]; ++m)
            {
                OPERA assign_p_infinity (&data_[g].shapeintmaskminmax[m][0]);
                OPERA assign_p_infinity (&data_[g].shapeintmaskminmax[m][1]);
                OPERA assign_n_infinity (&data_[g].shapeintmaskminmax[m][2]);
                OPERA assign_n_infinity (&data_[g].shapeintmaskminmax[m][3]);
                for (PS_SIT s=0; s<data_[g].extlengths[18]/2; ++s)
                {
                    if (data_[g].shapeintmask[m][s*2]  <data_[g].shapeintmaskminmax[m][0])
                        data_[g].shapeintmaskminmax[m][0] = data_[g].shapeintmask[m][s*2];
                    if (data_[g].shapeintmask[m][s*2]  >data_[g].shapeintmaskminmax[m][2])
                        data_[g].shapeintmaskminmax[m][2] = data_[g].shapeintmask[m][s*2];
                    if (data_[g].shapeintmask[m][s*2+1]<data_[g].shapeintmaskminmax[m][1])
                        data_[g].shapeintmaskminmax[m][1] = data_[g].shapeintmask[m][s*2+1];
                    if (data_[g].shapeintmask[m][s*2+1]>data_[g].shapeintmaskminmax[m][3])
                        data_[g].shapeintmaskminmax[m][3] = data_[g].shapeintmask[m][s*2+1];
                }
            }
        }
        else
        {
            MEMORY ps_free( polys, dim1 );
        }
    }

    // read source mask
    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        char *file;
        const char *listcc[5];
        listcc[0] = CONSTANT dir_in;
        listcc[1] = bn;
        listcc[2] = CONSTANT bnseparator;
        listcc[3] = name[g];
        listcc[4] = CONSTANT srcmaskfilereg;
        OPERA concatenate( listcc, 5, &file );

        if (data_[g].srcinputcircle)
            MEMORY ps_free (data_[g].srcinputcircle, data_[g].extlengths[2]);
        readmaskcircle( &data_[g].extlengths[2],  &data_[g].extlengths[3],
                        &data_[g].srcinputcircle, &file, &data_[g], cdata_ );
        MEMORY ps_free( file );
    }
}

void pixsrc_init::create_mm_border (double ***polys_, PS_SIT dim1, PS_SIT dim2, PS_SIT src, PS_SIT g, inputdata *data_, commoninputdata *cdata_)
{
    double steps = CONSTANT mm_stepsize; // steps in pixel units

    double **polys = *polys_;
    MEMORY ps_malloc (&(data_[g].mmborder[src]), dim1);
    MEMORY ps_malloc (&(data_[g].mmborder_defl[src]), dim1);

    // iterate over images to get number of total points
    for (PS_SIT i=0; i<dim1; ++i)
    {
        GEOM create_border (polys[i], dim2, steps, &data_[g].mmborder[src][i]);

        MEMORY ps_malloc (&(data_[g].mmborder_defl[src][i]), (PS_SIT)(data_[g].mmborder[src][i][0])*2+1);
        data_[g].mmborder_defl[src][i][0] = data_[g].mmborder[src][i][0];
    }
}

void pixsrc_init::setuppsfs(char *bn, char **name, char **namewithext, inputdata *data_, commoninputdata *cdata_)
{
    // set up psf's

    double psf_flux_cutoff = 0.99;//0.95;//0.99;
    double psf_flux_half_cutoff = 0.50;

    for(PS_SIT g=0; g<cdata_->numimages; g++)
    {
        if(data_[g].nopsf)
        {
            data_[g].extlengths[9] = 0;
            data_[g].extlengths[10] = 0;
        }
        else if(data_[g].psffromfile)
        {
            char *filename;
            const char *listcc[3];
            listcc[0] = CONSTANT dir_in;
            listcc[1] = name[g];
            listcc[2] = CONSTANT psffile;
            OPERA concatenate( listcc, 3, &filename );

            fitsfile *fptr;
            int/*PS_SIT*/ status = 0;
            fits_open_file(&fptr, filename, READONLY, &status);
            MEMORY ps_free( filename );
            int/*PS_SIT*/ naxis, bitpix;
            long naxes[2] = {1,1}, fpixel[2] = {1,1};
            fits_get_img_param(fptr, 2, &bitpix, &naxis, naxes, &status);

            if(naxes[0]%2==0 || naxes[1]%2==0 || naxes[0]!=naxes[1])
                PRINTER printerror( "","Code does not currently support"
                                    "even-dimension or non-square PSFs",
                                    cdata_->print2screenmutex            );

            data_[g].extlengths[9] = naxes[0];
            data_[g].extlengths[10] = naxes[1];
            MEMORY ps_malloc( &(data_[g].psf), naxes[0], naxes[1] );

            double *pixels;
            MEMORY ps_malloc( &(pixels), naxes[0] );
            for (fpixel[1] = naxes[1]; fpixel[1] >= 1; fpixel[1]--)
            {
                fits_read_pix(fptr, TDOUBLE, fpixel, naxes[0], NULL, pixels, NULL, &status);
                for(PS_SIT ii=0; ii<naxes[0]; ii++)
                    data_[g].psf[ii][fpixel[1]-1] = pixels[ii];
            }
            MEMORY ps_free( pixels );
            fits_close_file(fptr, &status);

            PS_SIT numtotal=naxes[0]*naxes[1];
            double vec[numtotal];
            PS_SIT counter=0;
            for(PS_SIT z1=0; z1<naxes[0]; z1++)
                for(PS_SIT z2=0; z2<naxes[1]; z2++)
                    vec[counter++]=data_[g].psf[z2][z1];
            std::sort(vec, vec+numtotal);
            double sum=0,total_sum=0;
	    PS_SIT zblah;
            PS_SIT fwhm_flag = 0;
            for(zblah=numtotal-1; zblah>=0; zblah--)
	        total_sum += vec[zblah];
            for(zblah=numtotal-1; zblah>=0; zblah--)
            {
                if(sum>=psf_flux_half_cutoff*total_sum && !fwhm_flag )
                {
                    data_[g].fwhm_weight = vec[zblah+1];
                    fwhm_flag = 1;
                }
                if(sum>=psf_flux_cutoff*total_sum)
                    break;
                sum+=vec[zblah];
            }
            double minweight = vec[zblah+1];

            for(PS_SIT z1=0; z1<naxes[0]; z1++)
                for(PS_SIT z2=0; z2<naxes[1]; z2++)
                    if(data_[g].psf[z2][z1]<minweight)
                        data_[g].psf[z2][z1]=0;
        }
        else
        {
            // convert arcseconds to pixel coordinates
            data_[g].majoraxis *= data_[g].arc2pix;
            data_[g].minoraxis *= data_[g].arc2pix;
            // get position angle w.r.t. to x-y axes
            data_[g].angleaxis *= CONSTANT deg2rad;
            data_[g].angleaxis -= data_[g].rollangle;

            double sum = 0;
            double integrallimit = CONSTANT cutoff;
            data_[g].extlengths[9] = data_[g].extlengths[10] = 0;

            double sigmax = data_[g].majoraxis/CONSTANT fwhm2sigma;//*pow(1.1,counter);
            double sigmay = data_[g].minoraxis/CONSTANT fwhm2sigma;//*pow(1.1,counter);

            double cos1 = std::cos(data_[g].angleaxis);
            double sin1 = std::sin(data_[g].angleaxis);
            double fac1 = 1/(2.0*CONSTANT pi*sigmax*sigmay);
            double facx = -1/(2.0*sigmax*sigmax);
            double facy = -1/(2.0*sigmay*sigmay);

            while(sum<psf_flux_cutoff)
            {
                MEMORY ps_free( data_[g].psf, data_[g].extlengths[9] );

                PS_SIT minmax = (PS_SIT)std::ceil(integrallimit*sigmax);

                double foci = integrallimit*std::sqrt(sigmax*sigmax-sigmay*sigmay);
                double f1x = foci*std::fabs(sin1);
                double f1y = foci*cos1;
                if(data_[g].angleaxis>=0)
                    f1x *= -1;
                double f2x = -f1x;
                double f2y = -f1y;

                data_[g].extlengths[9] = data_[g].extlengths[10] = 2*minmax+1;
                MEMORY ps_malloc( &(data_[g].psf), data_[g].extlengths[9], data_[g].extlengths[10] );

                PS_SIT xty = data_[g].extlengths[9]*data_[g].extlengths[10];
                double *psf1d;
                MEMORY ps_malloc( &psf1d, xty );

                // for oversampling of PSF
                // spacing = spacing between oversampled points
                // start_position = start position for oversampling, relative to the pixel center
                PS_SIT oversample = data_[g].psf_oversample;
                double spacing = 1.0/oversample;
                double start_position = -(oversample-1)/2.0*spacing;

                bool goodpixel;
                for(PS_SIT xx=0; xx<data_[g].extlengths[9]; xx++)
                {
                    for(PS_SIT yy=0; yy<data_[g].extlengths[10]; yy++)
                    {
                        double distf1, distf2, subx, suby;
                        double x = xx-minmax;
                        double y = yy-minmax;

                        goodpixel = 1;

                        double sumflux, xrot, yrot;

                        if( goodpixel )
                        {
                            sumflux = 0;
                            for( PS_SIT r=0; r<oversample; ++r )
                            {
                                for( PS_SIT d=0; d<oversample; ++d )
                                {
                                    // find positions to oversample inside this pixel
                                    // check if position lies inside integration region
                                    // rotate the position according to orientation of major axis
                                    subx = x + start_position + r*spacing;
                                    suby = y + start_position + d*spacing;
                                    distf1 = OPERA distance(subx,suby,f1x,f1y);
                                    distf2 = OPERA distance(subx,suby,f2x,f2y);
                                    if( distf1+distf2 <= 2.0*sigmax*integrallimit )
                                    {
                                        xrot = subx*sin1-suby*cos1;
                                        yrot = subx*cos1+suby*sin1;
                                        sumflux += fac1*exp(facx*xrot*xrot)*exp(facy*yrot*yrot);
                                    }
                                }
                            }
                            data_[g].psf[xx][yy] = sumflux / (oversample*oversample);
                        }
                        else
                        {
                            data_[g].psf[xx][yy] = 0;
                        }
                    }
                }

                sum=0;
                PS_SIT psf1d_ind = 0;
                for(PS_SIT xx=0; xx<data_[g].extlengths[9]; xx++)
                    for(PS_SIT y=0; y<data_[g].extlengths[10]; y++)
                    {
                        sum += data_[g].psf[xx][y];
                        psf1d[psf1d_ind++] = data_[g].psf[xx][y];
                    }
                integrallimit *= 1.1;

                std::sort (psf1d, psf1d+xty);
                double sum_fwhm_weight=0;
                PS_SIT zblah;
                for(zblah=xty-1; zblah>=0; zblah--)
                {
                    if (sum_fwhm_weight>=psf_flux_half_cutoff)
                    {
                        data_[g].fwhm_weight = psf1d[zblah+1];
                        break;
                    }
                    sum_fwhm_weight+=psf1d[zblah];
                }
                MEMORY ps_free (psf1d);
            }
        }
    }
}

void pixsrc_init::printsettings(inputdata *data_, commoninputdata *cdata_)
{
    PS_SIT leave = 1;
    for(PS_SIT g=0; g<cdata_->numimages; g++)
        if(data_[g].verbose>1)
        {
            leave = 0;
            break;
        }
    if (leave) return;

    std::cout << std::endl;
    PRINTER print2screen("pixsrc", "pixsrc data and settings loaded", cdata_->print2screenmutex);
    std::cout << std::endl;
}
