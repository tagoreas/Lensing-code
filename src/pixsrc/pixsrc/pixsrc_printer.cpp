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



#include "pixsrc_printer.hpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_memory_templates.cpp"
#include <fitsio.h>
#include <cstdlib>                                      // for exit command
#include <cstring>
#include <iomanip>

using std::cout;
using std::endl;
using std::cerr;

void pixsrc_printer::openoutstream(string basename, string imagename, bool num,
                                   string filename, std::ofstream *stream, bool nonstandardname)
{
    if (!nonstandardname)
    {
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    stream->open(filename.c_str());
}
void pixsrc_printer::closeoutstream(std::ofstream *stream)
{
    stream->close();
}

void pixsrc_printer::print(PS_SIT tracker, string basename, string imagename, bool num,
                           string filename, PS_SIT nrow, PS_SIT ncol, MATRIX *source, bool nonstandardname)
{
    source->print_me( tracker, basename, imagename, num, filename, nonstandardname );
}

void pixsrc_printer::printfitssrcplane(PS_SIT tracker, string basename, string imagename,
                                       string print2screenname, bool num, string filenametemp,
                                       PS_SIT imgx, PS_SIT imgy, double *im, double *wcsinfo,
                                       pthread_mutex_t *lock, bool writehdr, bool nonstandardname)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filenametemp = imagename + "_" + filenametemp;
        filenametemp = basename + "_" + filenametemp;
        filenametemp = CONSTANT dir_out + filenametemp;
    }


    fitsfile *fptr;
    int status = 0; //PS_SIT status=0;
    long  fpixel, nelements;

    /* initialize FITS image parameters */
    char filename[filenametemp.size()+1];
    strcpy(filename,filenametemp.c_str());

    PS_SIT bitpix   =  DOUBLE_IMG;

    long naxis    =   2;  /* 2-dimensional image */
    /* image is im.size() pixels wide by "other" rows */
    long naxes[2] = { imgx, imgy };
    /* allocate memory for the whole image */
    double **array;
    MEMORY ps_malloc( &array, naxes[1] );

    /* initialize pointers to the start of each row of the image */
    for(PS_SIT ii=0; ii<naxes[1]; ii++)
        array[ii] = im + naxes[0]*ii;

    remove(filename);               /* Delete old file if it already exists */

    if(fits_create_file(&fptr, filename, &status)) /* create new FITS file */
    {
        printwarning( print2screenname, "error writing " + filenametemp, lock);
        MEMORY ps_free (array);
        return;
    }

    if(fits_create_img(fptr,  bitpix, naxis, naxes, &status))
    {
        printwarning( print2screenname, "error writing " + filenametemp, lock);
        MEMORY ps_free (array);
        return;
    }

    fpixel = 1;                               /* first pixel to write      */
    nelements = naxes[0] * naxes[1];          /* number of pixels to write */

    /* write the array of unsigned integers to the FITS file */

    if(fits_write_img(fptr, TDOUBLE, fpixel, nelements, array[0], &status))
    {
        PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
        MEMORY ps_free (array);
        return;
    }

    MEMORY ps_free (array);

    if(writehdr)
    {
        char ra[9], dec[9], ctype1[7], ctype2[7], empty[1], crpix1[7], crpix2[7],
            crval1[7], crval2[7], cd11[6], cd22[6], cd12[6], cd21[6];
        empty[0] = '\0';
        strcpy( ra    , string("RA---TAN").c_str() );
        strcpy( dec   , string("DEC--TAN").c_str() );
        strcpy( ctype1, string("CTYPE1"  ).c_str() );
        strcpy( ctype2, string("CTYPE2"  ).c_str() );
        strcpy( crval1, string("CRVAL1"  ).c_str() );
        strcpy( crval2, string("CRVAL2"  ).c_str() );
        strcpy( crpix1, string("CRPIX1"  ).c_str() );
        strcpy( crpix2, string("CRPIX2"  ).c_str() );
        strcpy( cd11  , string("CD1_1"   ).c_str() );
        strcpy( cd12  , string("CD1_2"   ).c_str() );
        strcpy( cd21  , string("CD2_1"   ).c_str() );
        strcpy( cd22  , string("CD2_2"   ).c_str() );

        if(wcsinfo[8]==1)
        {
            if (fits_update_key(fptr, TSTRING, ctype1, ra, empty, &status))
            {
                PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
                return;
            }
            if (fits_update_key(fptr, TSTRING, ctype2, dec, empty, &status))
            {
                PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
                return;
            }
        }
        else if(wcsinfo[8]==0)
        {
            if (fits_update_key(fptr, TSTRING, ctype2, ra, empty, &status))
            {
                PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
                return;
            }
            if (fits_update_key(fptr, TSTRING, ctype1, dec, empty, &status))
            {
                PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
                return;
            }
        }

        // perhaps we can still include double-precision WCS even if
        // data is single precision .. check this later

        if (fits_update_key(fptr, TDOUBLE, crpix1, &wcsinfo[0],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, crpix2, &wcsinfo[1],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, crval1, &wcsinfo[2],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, crval2, &wcsinfo[3],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, cd11  , &wcsinfo[4],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, cd12  , &wcsinfo[5],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, cd21  , &wcsinfo[6],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
        if (fits_update_key(fptr, TDOUBLE, cd22  , &wcsinfo[7],  empty, &status))
        {
            PRINTER printwarning( print2screenname, "error writing " + filenametemp, lock);
            return;
        }
    }

    if ( fits_close_file(fptr, &status) )                /* close the file */
        printwarning( print2screenname, "error writing " + filenametemp, lock);
}
void pixsrc_printer::printfitsimgplane(PS_SIT tracker, string basename, string imagename,
                                       string print2screenname, bool num, string filename,
                                       const double *vec, double *wcsinfo, PS_SIT imgx, PS_SIT imgy,
                                       double scalefactor, PS_SIT *r4r, pthread_mutex_t *lock,
                                       bool writehdr, bool nonstandardname)
{
    double *im;
    MEMORY ps_malloc( &im, imgx*imgy );
    PS_SIT index=0;
    for(PS_SIT y=imgy-1; y>=0; y--)
        for(PS_SIT x=0; x<imgx; x++)
        {
            PS_SIT s=x*imgy+y;
            if(r4r[s]!=-1)
                im[index++]=vec[r4r[s]]*scalefactor;
            else
                im[index++]=0;
        }

    printfitssrcplane( tracker, basename, imagename, print2screenname, num,
                       filename,imgx,imgy,im,wcsinfo, lock, writehdr, nonstandardname );

    MEMORY ps_free( im );
}
void pixsrc_printer::print2screen(string imagename, string message, pthread_mutex_t *lock)
{
    struct tm current;
    time_t now;

    if (lock)
        pthread_mutex_lock(lock);

    time(&now);
    current = *localtime(&now);
    string hour,min,sec="";
    if(current.tm_hour<10)
        hour="0";
    if(current.tm_min<10)
        min="0";
    if(current.tm_sec<10)
        sec="0";
    cout << imagename << " @ " << hour << current.tm_hour << ":" << min << current.tm_min <<
        ":" << sec << current.tm_sec << " :: " << message << endl;

    if (lock)
        pthread_mutex_unlock(lock);
}

void pixsrc_printer::printwarning(string imagename, string message, pthread_mutex_t *lock)
{
    struct tm current;
    time_t now;

    if (lock)
        pthread_mutex_lock(lock);

    time(&now);
    current = *localtime(&now);
    string hour,min,sec="";
    if(current.tm_hour<10)
        hour="0";
    if(current.tm_min<10)
        min="0";
    if(current.tm_sec<10)
        sec="0";
    cerr << imagename << " @ " << hour << current.tm_hour << ":" << min << current.tm_min <<
        ":" << sec << current.tm_sec << " :: WARNING: " << message << endl;

    if (lock)
        pthread_mutex_unlock(lock);
}
void pixsrc_printer::printerror(string imagename, string message, pthread_mutex_t *lock)
{
    struct tm current;
    time_t now;

    if (lock)
        pthread_mutex_lock(lock);

    time(&now);
    current = *localtime(&now);
    string hour,min,sec="";
    if(current.tm_hour<10)
        hour="0";
    if(current.tm_min<10)
        min="0";
    if(current.tm_sec<10)
        sec="0";
    cerr << imagename << " @ " << hour << current.tm_hour << ":" << min << current.tm_min <<
        ":" << sec << current.tm_sec << " :: ERROR: " << message << endl;

    exit(1);

    if (lock)
        pthread_mutex_unlock(lock);
}

void pixsrc_printer::printtristruct(struct triangulateio *tri, pthread_mutex_t *lock)
{
    // this fucntion is mainly for debugging

    if (lock)
        pthread_mutex_lock(lock);

    cout << endl;
    cout << "triangulateio struct  -> " << tri                        << endl;
    cout << "num of points  --------> " << tri->numberofpoints        << endl;
    cout << "num of tris  ----------> " << tri->numberoftriangles     << endl;
    cout << "num of segments  ------> " << tri->numberofsegments      << endl;
    cout << "pointlist  ------------> " << tri->pointlist             << endl;
    cout << "pointattributelist  ---> " << tri->pointattributelist    << endl;
    cout << "pointmarkerlist  ------> " << tri->pointmarkerlist       << endl;
    cout << "trianglelist  ---------> " << tri->trianglelist          << endl;
    cout << "triangleattributelist -> " << tri->triangleattributelist << endl;
    cout << "trianglearealist  -----> " << tri->trianglearealist      << endl;
    cout << "neighborlist  ---------> " << tri->neighborlist          << endl;
    cout << "segmentlist  ----------> " << tri->segmentlist           << endl;
    cout << "segmentmarkerlist  ----> " << tri->segmentmarkerlist     << endl;
    cout << "holelist  -------------> " << tri->holelist              << endl;
    cout << "regionlist  -----------> " << tri->regionlist            << endl;
    cout << "edgelist  -------------> " << tri->edgelist              << endl;
    cout << "edgemarkerlist  -------> " << tri->edgemarkerlist        << endl;
    cout << "normlist  -------------> " << tri->normlist              << endl;
    cout << endl;

    if (lock)
        pthread_mutex_unlock(lock);
}

PS_SIT pixsrc_printer::remove_file(PS_SIT tracker, string basename, string imagename, bool num,
                                   string filename, bool nonstandardname)
{
    if (!nonstandardname)
    {
        imagename = imagename + "_" + OPERA tostring(tracker);
        if(num)
            filename = imagename + "_" + filename;
        filename = CONSTANT dir_out + basename + "_" + filename;
    }

    return std::remove (filename.c_str());
}
