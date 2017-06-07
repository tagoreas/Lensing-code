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
#include "pixsrc_geometry.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_constants.hpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_tps.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_shapelets_operations.hpp"
#include "pixsrc_tps.hpp"
#include "pixsrc_d2matrix.cpp"
#include <cstring>
#include <algorithm>
#include <cmath>
#include <gsl/gsl_integration.h>

// check r.a. increase leftwards and y increasing downwards stuff .. TPS has a -1 fudge factor that seems to work fr now .. but it coul dbreak if pointing center shifted

void ps_test_tps(commoninputdata *cdata_, inputdata *data_);

void pixsrc_init::uvdata(char *bn, char **name, char **namewithext,
                         inputdata *data_, commoninputdata *cdata_ )
{
    bool atleast1 = 0;
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
        if (data_[g].is_uvdata)
            if (1!=data_[g].verbose)
                atleast1 = 1;

    if (atleast1)
        PRINTER print2screen("pixsrc",
                             "\n\n·´¯`·.¸¸.·´¯`·.¸><(((º>\n"
                             "                        <º)))><¸.·´¯`·.¸¸.·´¯`·\n"
                             "\nNote: uv-plane pre-processing can be a cpu/memory "
                             "intensive step and take a long time."
                             "\nHowever, they are one-time calculations, and "
                             "the results can be saved and loaded for next time."
                             "\nUse the help command or see the documentation "
                             "for help on saving/loading uv-plane parameters.\n"
                             "\n·´¯`·.¸¸.·´¯`·.¸><(((º>\n"
                             "                        <º)))><¸.·´¯`·.¸¸.·´¯`·\n\n",
                             cdata_->print2screenmutex);


    // for testing only
    //ps_test_tps (cdata_, data_);

    // check for saved binary files
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
        if (data_[g].is_uvdata && 1==data_[g].transmode)
        {
            // check if image to uv plane matrices already exist
            string fn1 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                + string(data_[g].name) + ".uvmat1.bin";
            string fn2 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                + string(data_[g].name) + ".uvmat2.bin";
            string fn3 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                + string(data_[g].name) + ".uvmat3.bin";
            data_[g].uv_matvecsca_exist = OPERA fileexists (fn1)
                && OPERA fileexists(fn2) && OPERA fileexists(fn3);
        }

    // get data
    INIT uvdata_read (data_, cdata_);

    // adjust image plane data to match uv data
    INIT uvdata_adjustimgplane (data_, cdata_);

    // re-read image plane masks
    char **names;
    MEMORY ps_malloc (&names, cdata_->numimages);
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
        names[g] = data_[g].name;
    INIT readmasks (cdata_->basename, names, 0, data_, cdata_);
    MEMORY ps_free (names);

    // compute shapelet decomposition of image plane and transformation to uv plane
    INIT uvdata_img2uv (data_, cdata_);
}

void printstruct (ps_WorldCoor* s)
{
    std::cout << s->xref << " ";
    std::cout << s->yref << " ";
    std::cout << s->xrefpix << " ";
    std::cout << s->yrefpix << " ";
    std::cout << s->xinc << " ";
    std::cout << s->yinc << " ";
    std::cout << s->rot << " ";
    std::cout << s->cd[0] << " ";
    std::cout << s->cd[1] << " ";
    std::cout << s->cd[2] << " ";
    std::cout << s->cd[3] << " ";
    std::cout << s->dc[0] << " ";
    std::cout << s->dc[1] << " ";
    std::cout << s->dc[2] << " ";
    std::cout << s->dc[3] << " ";
    std::cout << s->equinox << " ";
    std::cout << s->epoch << " ";
    std::cout << s->nxpix << " ";
    std::cout << s->nypix << " ";
    std::cout << s->plate_ra << " ";
    std::cout << s->plate_dec << " ";
    std::cout << s->plate_scale << " ";
    std::cout << s->x_pixel_offset << " ";
    std::cout << s->y_pixel_offset << " ";
    std::cout << s->x_pixel_size << " ";
    std::cout << s->y_pixel_size << " ";
    std::cout << s->eqin << " ";
    std::cout << s->eqout << " ";
    std::cout << s->sysin << " ";
    std::cout << s->sysout << " ";
    std::cout << s->syswcs << " ";
    std::cout << std::endl;
}

void pixsrc_init::uvdata_adjustimgplane (inputdata *data_, commoninputdata *cdata_)
{
    pthread_mutex_lock (cdata_->wcsmutex);
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
    {
        if (data_[g].is_uvdata)
        {
            if (1!=data_[g].verbose)
                PRINTER print2screen(data_[g].print2screenname,
                                     "adjusting image plane pixel scale",
                                     cdata_->print2screenmutex);

            double newscale;
            if (data_[g].uv_matvecsca_exist)
            {
                double sumsave[3];
                string fn3 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                    + string(data_[g].name) + ".uvmat3.bin";
                PRINTER readbinary (fn3, sumsave, 3, cdata_->print2screenmutex);
                newscale = sumsave[2];
            }
            else
            {
                if (data_[g].uv_newpixelscale<0)
                {
                    // get max baseline
                    double bl, maxbaseline;
                    OPERA assign_n_infinity (&maxbaseline);
                    for (PS_SIT r=0; r<data_[g].uv_ndp; ++r)
                    {
                        bl = data_[g].uv_oldloc[r*2+0] *data_[g].uv_oldloc[r*2+0] +
                            data_[g].uv_oldloc[r*2+1] *data_[g].uv_oldloc[r*2+1];
                        if (bl>maxbaseline)
                            maxbaseline = bl;
                    }
                    maxbaseline = std::sqrt (maxbaseline);
                    newscale = 1.0 /maxbaseline *CONSTANT rad2arc *-data_[g].uv_newpixelscale;
                }
                else
                {
                    newscale = data_[g].uv_newpixelscale;
                }
            }
            data_[g].uv_newpixelscale = newscale;

            // get new (adjusted) data parameters
            PS_SIT newimgx, newimgy, newndp;
            double *newdata, newcrpix[2], newcrval[2], newcdmat[4];

            if (1!=data_[g].verbose)
                PRINTER print2screen(data_[g].print2screenname,
                                     "old, new pixel scales: " + OPERA tostring(data_[g].pix2arc) +
                                     ", " + OPERA tostring(newscale),
                                     cdata_->print2screenmutex);

            //newscale = data_[g].pix2arc;
            //newscale = 0.1458512432174512800/4;

            // make new grid slightly larger than old one (new grid is centered on old grid)
            newimgx = (PS_SIT)std::ceil (data_[g].imgx *data_[g].pix2arc /newscale) +1;
            newimgy = (PS_SIT)std::ceil (data_[g].imgy *data_[g].pix2arc /newscale) +1;

            //newimgx = data_[g].imgx;
            //newimgy = data_[g].imgy;

            newndp = newimgx *newimgy;
            // set reference pixel at center of grid and get reference WCS value
            newcrpix[0] = 0.5*(newimgx-1);
            newcrpix[1] = 0.5*(newimgy-1);
            HEADER getimgwcscoord (data_[g].wcs, data_[g].imgy, 0.5*(data_[g].imgx-1), 0.5*(data_[g].imgy-1),
                                   &newcrval[0], &newcrval[1]);

            // create new WCS struct (using a fake FITS header)
            std::copy (data_[g].wcsinfo+4, data_[g].wcsinfo+8, newcdmat);
            std::transform (newcdmat, newcdmat+4, newcdmat,
                            std::bind1st(std::multiplies<double>(), newscale/data_[g].pix2arc));
            // fill header with spaces and null terminate
            char *fakeheader;
            // each header entry must be exactly 80 characters wide
            PS_SIT numhlines = 16, linewidth=80, headind=0;
            MEMORY ps_malloc (&fakeheader, numhlines*linewidth+1);
            std::fill (fakeheader, fakeheader+numhlines*linewidth, ' ');
            fakeheader[numhlines*linewidth] = 0;
            // fill in values
            strcpy (fakeheader+linewidth*headind++, "SIMPLE  = T");
            strcpy (fakeheader+linewidth*headind++, "NAXIS   = 2");
            strcpy (fakeheader+linewidth*headind++, "EQUINOX = 2000");
            strcpy (fakeheader+linewidth*headind++, "CTYPE1  = 'RA---TAN'");
            strcpy (fakeheader+linewidth*headind++, "CTYPE2  = 'DEC--TAN'");
            strcpy (fakeheader+linewidth*headind++, "NAXIS1  = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%10i",  (int)newimgx);
            strcpy (fakeheader+linewidth*headind++, "NAXIS2  = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%10i",  (int)newimgy);
            strcpy (fakeheader+linewidth*headind++, "CD1_1   = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcdmat[0]);
            strcpy (fakeheader+linewidth*headind++, "CD2_2   = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcdmat[3]);
            strcpy (fakeheader+linewidth*headind++, "CD1_2   = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcdmat[1]);
            strcpy (fakeheader+linewidth*headind++, "CD2_1   = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcdmat[2]);
            strcpy (fakeheader+linewidth*headind++, "CRPIX1  = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcrpix[0]+1);
            strcpy (fakeheader+linewidth*headind++, "CRPIX2  = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newimgy-newcrpix[1]);
            strcpy (fakeheader+linewidth*headind++, "CRVAL1  = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcrval[0]);
            strcpy (fakeheader+linewidth*headind++, "CRVAL2  = ");
            snprintf (fakeheader+linewidth*(headind-1)+10, 40, "%.18E", newcrval[1]);
            strcpy (fakeheader+linewidth*headind++, "END");
            // remove null characters that strcpy copied and create WCS struct
            for (PS_SIT hs=0; hs<linewidth*numhlines; ++hs)
                if (!fakeheader[hs])
                    fakeheader[hs] = ' ';
            ps_WorldCoor *newwcs = EXTERNAL ps_wcsinit (fakeheader);
            //std::cout << fakeheader << std::endl;
            MEMORY ps_free (fakeheader);

            // interpolate original data onto new data
            MEMORY ps_malloc (&newdata, newndp);
	    std::fill (newdata, newdata+newndp, 0);
            PS_SIT surrpos[4];
            double norm, flux, wtx, wty, oldx, newx, oldy, newy;
            double dxdx = newscale /data_[g].pix2arc;
            for (PS_SIT x=0; x<newimgx; ++x)
            {
                newx = x - 0.5*(newimgx-1);
                oldx = newx *dxdx + 0.5*(data_[g].imgx-1);
                for (PS_SIT y=0; y<newimgy; ++y)
                {
                    newy = y - 0.5*(newimgy-1);
                    oldy = newy *dxdx + 0.5*(data_[g].imgy-1);

                    // corners of old grid pixels that surround new pixel
                    surrpos[0] = (PS_SIT)std::floor (oldx);
                    surrpos[1] = (PS_SIT)std::floor (oldy);
                    surrpos[2] = (PS_SIT)std::ceil  (oldx);
                    surrpos[3] = (PS_SIT)std::ceil  (oldy);
                    flux = norm = 0;
                    for (PS_SIT xx=0; xx<2; ++xx)
                    {
                        if (surrpos[xx*2+0]<0 || surrpos[xx*2+0]>=data_[g].imgx)
                            continue;
                        for (PS_SIT yy=0; yy<2; ++yy)
                        {
                            if (surrpos[yy*2+1]<0 || surrpos[yy*2+1]>=data_[g].imgy)
                                continue;

                            wtx   = 1 - std::abs (surrpos[xx*2+0]-oldx);
                            wty   = 1 - std::abs (surrpos[yy*2+1]-oldy);
                            norm += wtx*wty;
                            flux += wtx*wty
                                *data_[g].data[surrpos[xx*2+0]*data_[g].imgy+surrpos[yy*2+1]];
                        }
                    }
		    if (norm)
			newdata[x*newimgy+y] = flux /norm;
                }
            }
            //printstruct(newwcs);
            //printstruct(data_[g].wcs);

            // overwrite original image plane data
            EXTERNAL ps_wcsfree (data_[g].wcs);
            MEMORY ps_free (data_[g].data);
            data_[g].wcs = newwcs;
            data_[g].data = newdata;
            data_[g].imgx = newimgx;
            data_[g].imgy = newimgy;
            data_[g].ndp = newndp;
            //data_[g].px = newcrpix[0];
            //data_[g].py = newcrpix[1];
            //data_[g].pra = newcrval[0];
            //data_[g].pdec = newcrval[1];
            data_[g].wcs->degout = 1;
            data_[g].wcs->ndec = 10;
            data_[g].wcsinfo[0] = newcrpix[0]+1;
            data_[g].wcsinfo[1] = newimgy-newcrpix[1];
            data_[g].wcsinfo[2] = newcrval[0];
            data_[g].wcsinfo[3] = newcrval[1];
            std::copy (newcdmat, newcdmat+4, data_[g].wcsinfo+4);
            data_[g].wcsinfo[8] = 1;
            data_[g].rollangle = -std::atan2 (-data_[g].wcsinfo[5], data_[g].wcsinfo[7]);
            double factor=data_[g].wcsinfo[4]*data_[g].wcsinfo[7]-data_[g].wcsinfo[5]*data_[g].wcsinfo[6];
            data_[g].invwcsinfo[0]=data_[g].wcsinfo[7]/factor;
            data_[g].invwcsinfo[1]=-data_[g].wcsinfo[5]/factor;
            data_[g].invwcsinfo[2]=-data_[g].wcsinfo[6]/factor;
            data_[g].invwcsinfo[3]=data_[g].wcsinfo[4]/factor;
            data_[g].arc2pix = (
                std::sqrt( data_[g].invwcsinfo[0] * data_[g].invwcsinfo[0] +
                           data_[g].invwcsinfo[1] * data_[g].invwcsinfo[1] ) +
                std::sqrt( data_[g].invwcsinfo[2] * data_[g].invwcsinfo[2] +
                           data_[g].invwcsinfo[3] * data_[g].invwcsinfo[3] ) ) / 2.0 / 3600.0;
            data_[g].pix2arc = 1.0 / data_[g].arc2pix;
            data_[g].ndp = data_[g].imgx*data_[g].imgy;
            data_[g].extlengths[8] = data_[g].ndp;
            HEADER getimgpixcoord (data_[g].wcs, data_[g].imgy, cdata_, data_[g].pra, data_[g].pdec,
                                   &data_[g].px, &data_[g].py);
        }
    }
    pthread_mutex_unlock (cdata_->wcsmutex);
}

void pixsrc_init::uvdata_get_rbf_mat (inputdata *data_, commoninputdata *cdata_, TPS **imgtps, PS_SIT *num_ndp_, PS_SIT transmode)
{
    // count the number of image pixels being used to create the lensing operator
    PS_SIT num_ndp = 0;
    if (data_->useall)
        num_ndp = data_->ndp;
    else
        for (PS_SIT r=0; r<data_->ndp; ++r)
            if (2==data_->imagemasks[r])
                ++num_ndp;
    *num_ndp_ = num_ndp;

    if (!imgtps)
        return;

    // flag pixels as either masked or surrounding masked pixels
    PS_SIT searchscale = (PS_SIT)(std::ceil (4*data_->uv_rbfscale));
    char *imgflag;
    PS_SIT xpos, ypos, rnew, numtotalpts4tps=num_ndp;
    MEMORY ps_malloc (&imgflag, data_->ndp);
    std::fill (imgflag, imgflag+data_->ndp, 0);
    // flag masked pixels
    for (PS_SIT r=0; r<data_->ndp; ++r)
        if (2==data_->imagemasks[r])
            imgflag[r] = 1;
    // search pixels that surround masked pixels and flag them
    // if they are not masked themselves
    if (data_->uv_padzeros)
    {
        for (PS_SIT r=0; r<data_->ndp; ++r)
        {
            if (1==imgflag[r])
            {
                xpos = r/data_->imgy;
                ypos = r%data_->imgy;
                for (PS_SIT xnew=xpos-searchscale; xnew<=xpos+searchscale; ++xnew)
                {
                    if (xnew<0 || xnew>data_->imgx-1)
                        continue;
                    for (PS_SIT ynew=ypos-searchscale; ynew<=ypos+searchscale; ++ynew)
                    {
                        if (ynew<0 || ynew>data_->imgy-1)
                            continue;

                        rnew = xnew*data_->imgy+ynew;
                        if (0==imgflag[rnew])
                        {
                            imgflag[rnew] = 2;
                            ++numtotalpts4tps;
                        }
                    }
                }
            }
        }
    }

    // copy image positions
    double *imgpos;
    bool *ctrpts = NULL;
    PS_SIT x, y, tracker = 0;
    MEMORY ps_malloc (&imgpos, (PS_SIT)numtotalpts4tps*2);
    MEMORY ps_malloc (&ctrpts, (PS_SIT)numtotalpts4tps);
    std::fill (ctrpts, ctrpts+numtotalpts4tps, 0);
    for (PS_SIT r=0; r<data_->ndp; ++r)
    {
        x = r/data_->imgy;
        y = r%data_->imgy;
        if (1==imgflag[r])
            ctrpts[tracker] = 1;
        if (imgflag[r])
        {
            imgpos[tracker*2+0] = x;
            imgpos[tracker*2+1] = y;
            ++tracker;
        }
    }
    MEMORY ps_free (imgflag);

    if (1!=data_->verbose)
        PRINTER print2screen (data_->print2screenname,
                              "getting RBF decomposition of image plane",
                              cdata_->print2screenmutex);

    // create TPS (matrix only) for approximating image plane
    ps_rbf_type myrbf;
    switch (data_->uv_rbftype)
    {
    case 0:
        myrbf = ps_rbf_gaussian;
        break;
    default:
        PRINTER printerror(data_->print2screenname, "RBF type not recognized",
                           cdata_->print2screenmutex);
        break;
    }

    *imgtps = new TPS (cdata_, data_,
                       imgpos, imgpos+1, NULL,
                       2, 2, -1, numtotalpts4tps, ctrpts, myrbf, data_->uv_rbfscale, 1, -1, transmode);
    MEMORY ps_free (imgpos);
    MEMORY ps_free (ctrpts);
}

struct ps_uv_ft_getcol_struct
{
    inputdata *data_;
    TPS *imgtps;
    PS_SIT imgstart, imgend, uvstart, uvend, uvlen;
    PS_SIT threadnumber, numthreads;
    double *col, *fac2;
};
void *ps_uv_multithread_ft_getcol (void *args)
{
    ps_uv_ft_getcol_struct *parms = (ps_uv_ft_getcol_struct*) args;
    inputdata *data_ = parms->data_;
    TPS *imgtps      = parms->imgtps;
    PS_SIT imgstart  = parms->imgstart;
    PS_SIT imgend    = parms->imgend;
    PS_SIT uvstart   = parms->uvstart;
    PS_SIT uvend     = parms->uvend;
    PS_SIT uvlen     = parms->uvlen;
    PS_SIT threadnumber = parms->threadnumber;
    PS_SIT numthreads   = parms->numthreads;
    double *col      = parms->col;
    double *fac2     = parms->fac2;

    for (PS_SIT m=imgstart+threadnumber; m<imgend; m+=numthreads)
    {
        imgtps->get_fourier_weights_col (m, data_->uv_oldloc,
                                         col+(m-imgstart)*uvlen,
                                         data_->uv_ndp, fac2,
                                         uvstart/2, uvend/2);
    }
    return NULL;
}

void pixsrc_init::uvdata_img2uv (inputdata *data_, commoninputdata *cdata_)
{
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
    {
        if (data_[g].is_uvdata)
        {
            if (data_[g].verbose!=1)
                PRINTER print2screen(data_[g].print2screenname,
                                     "getting image plane to uv plane transformation (one-time calculation)",
                                     cdata_->print2screenmutex);

            // check if image to uv plane matrices already exist
            string fn1 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                + string(data_[g].name) + ".uvmat1.bin";
            string fn2 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                + string(data_[g].name) + ".uvmat2.bin";
            string fn3 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
                + string(data_[g].name) + ".uvmat3.bin";
            bool doexist = data_[g].uv_matvecsca_exist;

            // may need these
            TPS *imgtps = NULL;
            MATRIX *qtcdinvq = NULL;

            // get spatial decomposition of image plane using RBFs
            PS_SIT num_ndp = 0;
            PS_SIT numextraRBF = 0;
            if (!doexist)
                INIT uvdata_get_rbf_mat (&data_[g], cdata_, &imgtps, &num_ndp, data_[g].transmode);
            else
                INIT uvdata_get_rbf_mat (&data_[g], cdata_, NULL, &num_ndp, data_[g].transmode);

            if (!doexist)
                if (1!=data_[g].verbose)
                    PRINTER print2screen (data_[g].print2screenname,
                                          "getting image- to uv-plane matrices",
                                          cdata_->print2screenmutex);

            if (1==data_[g].transmode)
            {
                // uv matrices
                MATRIX *uv_transform_mat_ptr, *uv_transform_mat;
                VECTOR *uv_transform_vec_ptr, *uv_transform_vec;
                MEMORY ps_malloc (&uv_transform_vec_ptr, 1);
                MEMORY ps_malloc (&uv_transform_mat_ptr, 1);
                double sum, sumsave[3];

                if (!doexist)
                {
                    // get Q^t.Cd^-1.Q
                    // get Q^t.Cd^-1.d
                    qtcdinvq = new MATRIX (cdata_, &data_[g], num_ndp+numextraRBF,
                                           num_ndp+numextraRBF, -1, -1, &data_[g]);
                    VECTOR qtcdinvd (cdata_, &data_[g], num_ndp+numextraRBF);
                    PS_SIT lastuvblockind=-1;

                    // encapsulated to free memory quickly
                    {
                        // parms
                        // find minimum number of submat's needed (for memory allocation)
                        PS_SIT del1=data_[g].uv_del1, del2=data_[g].uv_del2;
                        del1 *=2;

                        del1    = std::min (del1, data_[g].uv_ndp*2);
                        del2    = std::min (del2, num_ndp+numextraRBF);
                        
                        // get submatrices for divide and conquer multiplication
                        MATRIX submat1 (cdata_, &data_[g], del1, del2, -1, -1, &data_[g]);
                        MATRIX submat2 (cdata_, &data_[g], del1, del2, -1, -1, &data_[g]);
                        double *col1 = submat1.mat_dense;
                        double *col2 = submat2.mat_dense;
                        PS_SIT start1, start2, start3, end1, end2, end3, len1, len2, len3;

                        // pre-compute some exponentials
                        double *fac2;
                        MEMORY ps_malloc (&fac2, (PS_SIT)data_[g].uv_ndp);
                        imgtps->precompute_fourier (data_[g].uv_oldloc, data_[g].uv_ndp, fac2);

                        if (1!=data_[g].verbose)
                            PRINTER print2screen (data_[g].print2screenname,
                                                  "Note: time needed is linear in number of uv blocks "
                                                  "and quadratic in number of image blocks",
                                                  cdata_->print2screenmutex);

                        // loop over blocks of uv data points
                        PS_SIT blockind=1;
                        PS_SIT numuvblocks  = (PS_SIT)std::ceil (data_[g].uv_ndp*2/(double)(del1));
                        PS_SIT numimgblocks = (PS_SIT)std::ceil ((num_ndp+numextraRBF)/(double)(del2));
                        string uvstr;
                        for (PS_SIT smblock=0; smblock<data_[g].uv_ndp*2; smblock+=del1)
                        {
                            start1 = smblock;
                            end1   = smblock + del1;
                            end1   = end1>data_[g].uv_ndp*2 ? data_[g].uv_ndp*2 : end1;
                            len1   = end1 - start1;

                            PS_SIT imgblockind = 1;
                            uvstr = "working on uv block " + OPERA tostring (blockind)
                                + " of " + OPERA tostring (numuvblocks) + ", image block ";

                            // loop over basis functions (left matrix)
                            for (PS_SIT cm=0; cm<num_ndp+numextraRBF; cm+=del2)
                            {
                                start2 = cm;
                                end2   = cm + del2;
                                end2   = end2>num_ndp+numextraRBF ? num_ndp+numextraRBF : end2;
                                len2   = end2 - start2;

                                if (1!=data_[g].verbose)
                                    PRINTER print2screen (data_[g].print2screenname,
                                                          uvstr + OPERA tostring (imgblockind) +
                                                          " of " + OPERA tostring (numimgblocks),
                                                          cdata_->print2screenmutex);

                                // get weights of basis functions at uv positions
                                // multithreaded
                                pthread_t *img1_threads;
                                ps_uv_ft_getcol_struct *img1_structs;
                                MEMORY ps_malloc (&img1_threads, cdata_->numthreads);
                                MEMORY ps_malloc (&img1_structs, cdata_->numthreads);
                                for (PS_SIT p=0; p<cdata_->numthreads; ++p)
                                {
                                    img1_structs[p].data_        = &data_[g];
                                    img1_structs[p].imgtps       = imgtps;
                                    img1_structs[p].imgstart     = start2;
                                    img1_structs[p].imgend       = end2;
                                    img1_structs[p].uvstart      = start1;
                                    img1_structs[p].uvend        = end1;
                                    img1_structs[p].uvlen        = len1;
                                    img1_structs[p].threadnumber = p;
                                    img1_structs[p].numthreads   = cdata_->numthreads;
                                    img1_structs[p].col          = col1;
                                    img1_structs[p].fac2         = fac2;
                                    pthread_create (&img1_threads[p], cdata_->attrjoinable,
                                                    ps_uv_multithread_ft_getcol, &img1_structs[p]);
                                }
                                for (PS_SIT p=0; p<cdata_->numthreads; ++p)
                                    pthread_join (img1_threads[p], NULL);
                                MEMORY ps_free (img1_threads);
                                MEMORY ps_free (img1_structs);

                                // multiply inverse covariance matrix
                                for (PS_SIT d=start1; d<end1; ++d)
                                    for (PS_SIT m=0; m<len2; ++m)
                                        col1[m*len1+d-start1] /= data_[g].uv_newvariance[d/2];

                                // loop over basis functions (right matrix)
                                for (PS_SIT cn=0; cn<num_ndp+numextraRBF; cn+=del2)
                                {
                                    start3 = cn;
                                    end3   = cn + del2;
                                    end3   = end3>num_ndp+numextraRBF ? num_ndp+numextraRBF : end3;
                                    len3   = end3 - start3;

                                    // only get approx half of symmetric matrix
                                    if (start3>=end2)
                                        continue;

                                    // get weights of basis functions at uv positions
                                    // multithreaded
                                    pthread_t *img2_threads;
                                    ps_uv_ft_getcol_struct *img2_structs;
                                    MEMORY ps_malloc (&img2_threads, cdata_->numthreads);
                                    MEMORY ps_malloc (&img2_structs, cdata_->numthreads);
                                    for (PS_SIT p=0; p<cdata_->numthreads; ++p)
                                    {
                                        img2_structs[p].data_        = &data_[g];
                                        img2_structs[p].imgtps       = imgtps;
                                        img2_structs[p].imgstart     = start3;
                                        img2_structs[p].imgend       = end3;
                                        img2_structs[p].uvstart      = start1;
                                        img2_structs[p].uvend        = end1;
                                        img2_structs[p].uvlen        = len1;
                                        img2_structs[p].threadnumber = p;
                                        img2_structs[p].numthreads   = cdata_->numthreads;
                                        img2_structs[p].col          = col2;
                                        img2_structs[p].fac2         = fac2;
                                        pthread_create (&img2_threads[p], cdata_->attrjoinable,
                                                        ps_uv_multithread_ft_getcol, &img2_structs[p]);
                                    }
                                    for (PS_SIT p=0; p<cdata_->numthreads; ++p)
                                        pthread_join (img2_threads[p], NULL);
                                    MEMORY ps_free (img2_threads);
                                    MEMORY ps_free (img2_structs);

                                    // do (partial) submatrix multiplications
                                    submat1.submultadd (&submat2, qtcdinvq, 1, 0, len1, len2, len1, len3,
                                                        len2, len3, start2, start3, cdata_->numthreads);
                                }
                                ++imgblockind;
                            }

                            // save current snapshot and remove previous snapshot
                            if (1!=data_[g].verbose)
                                PRINTER print2screen(data_[g].print2screenname,
                                                     "saving current state of uv-plane calculations",
                                                     cdata_->print2screenmutex);

                            PS_SIT status=0, uvstat=1;
                            if (!status)
                                status += PRINTER print <PS_SIT> (-1, cdata_->basename, data_->name,
                                                                  1, "snapshot."+OPERA tostring(blockind)+".stat",
                                                                  &uvstat, 1, 0, data_->precision, 0);
                            if (!status)
                                status += PRINTER printbinary <double> (-1, cdata_->basename, data_[g].name,
                                                                        1, "snapshot."+OPERA tostring(blockind)+".bin",
                                                                        qtcdinvq->mat_dense, qtcdinvq->nnz, 0);
                            uvstat = 0;
                            if (!status)
                                status += PRINTER print <PS_SIT> (-1, cdata_->basename, data_->name,
                                                                  1, "snapshot."+OPERA tostring(blockind)+".stat",
                                                                  &uvstat, 1, 0, data_->precision, 0);
                            if (!status && 1!=blockind)
                            {
                                PRINTER remove_file (-1, cdata_->basename, data_->name, 1,
                                                     "snapshot."+OPERA tostring(blockind-1)+".stat", 0);
                                PRINTER remove_file (-1, cdata_->basename, data_->name, 1,
                                                     "snapshot."+OPERA tostring(blockind-1)+".bin", 0);
                            }
                            if (status)
                                PRINTER printwarning(data_[g].print2screenname,
                                                     "error in writing current state of calculations to file",
                                                     cdata_->print2screenmutex);

                            lastuvblockind = blockind;
                            ++blockind;
                        }

                        if (1!=data_[g].verbose)
                            PRINTER print2screen (data_[g].print2screenname,
                                                  "few more calculations: 1 of 2",
                                                  cdata_->print2screenmutex);
                        // gug: merge this with above loop
                        double *col3;
                        MEMORY ps_malloc (&col3, (PS_SIT)data_[g].uv_ndp*2);
                        for (PS_SIT cm=0; cm<num_ndp+numextraRBF; ++cm)
                        {
                            // get Fourier weights for this basis function and rescale by noise
                            imgtps->get_fourier_weights_col (cm, data_[g].uv_oldloc, col3,
                                                             data_[g].uv_ndp, fac2, 0, data_[g].uv_ndp);
                            for (PS_SIT d=0; d<data_[g].uv_ndp*2; ++d)
                                col3[d] /= data_[g].uv_newvariance[d/2];

                            // get qtcdinvd
                            PS_SIT num1 = data_[g].uv_ndp*2, num2 = 1;
                            sum = EXTERNAL ps_ddot_ (&num1, col3, &num2, data_[g].uv_data, &num2);
                            qtcdinvd.set (cm, sum);
                        }
                        MEMORY ps_free (col3);

                        // cleanup
                        MEMORY ps_free (fac2);
                    }

                    // set other half of symmetric matrix
                    for (PS_SIT cm=0; cm<num_ndp+numextraRBF; ++cm)
                        for (PS_SIT cn=cm+1; cn<num_ndp+numextraRBF; ++cn)
                            qtcdinvq->set (cm, cn, qtcdinvq->get (cn, cm));

                    if (1!=data_[g].verbose)
                        PRINTER print2screen (data_[g].print2screenname,
                                              "few more calculations: 2 of 2",
                                              cdata_->print2screenmutex);

		    double *pbeam;
		    MEMORY ps_malloc (&pbeam, num_ndp);
		    if (data_[g].uv_pbeam)
		    {
			if (1!=data_[g].verbose)
			    PRINTER print2screen (data_[g].print2screenname,
						  "reading primary beam data",
						  cdata_->print2screenmutex);
			
			// read in primary beam data file
			string fn1 = string(CONSTANT dir_in) + string(cdata_->basename) + "_"
			    + string(data_[g].name) + ".pbeam.dat";
			PS_SIT nlines;
			char **pbdata;
			OPERA readfile (fn1.c_str(), &pbdata, &nlines, cdata_->print2screenmutex);

			string fn2;
			if (0==nlines)
			    fn2 = string(CONSTANT dir_in) + string(data_[g].name) + ".pbeam.fits";
			else if (1==nlines)
			    fn2 = string(CONSTANT dir_in) + string(data_[g].name)
				+ ".pbeam." + string(pbdata[0]) + ".fits";
			else if (data_[g].uv_ndporig==nlines)
			    PRINTER printerror (data_[g].name, "Baseline-specific primary beams are not yet implemented.", 
					       cdata_->print2screenmutex);
			else
			    PRINTER printerror (data_[g].name, fn1+" contains an invalid number of lines", 
						cdata_->print2screenmutex);
			
			if (!OPERA fileexists (fn2))
			    PRINTER printerror (data_[g].name, fn2+" does not exist", cdata_->print2screenmutex);

			// read data
			PS_SIT wcsinfolength, invwcsinfolength, imgdatalength, imgx, imgy, ndp;
			double *wcsinfo=NULL, *invwcsinfo=NULL, *imgdata=NULL, rollangle, pix2arc, arc2pix;
			INIT read_one_image (&data_[g], cdata_, fn2, 
					     &wcsinfolength, &wcsinfo,
					     &invwcsinfolength, &invwcsinfo,
					     &imgdatalength, &imgdata,
					     &imgx, &imgy, &rollangle,
					     &pix2arc, &arc2pix, &ndp);
			// setup WCS stuff
			ps_WorldCoor *wcs=NULL;
			char coorsysfinal[100];
			double px, py, pra=data_[g].pra, pdec=data_[g].pdec;
			INIT setup_one_libwcs (cdata_, fn2, &wcs, coorsysfinal,
					       &px, &py, &pra, &pdec, imgy);
			
			// loop over image pixels (lensed sky plane data)
			double img_pos_xy[2], img_pos_radec[2], pb_pos_xy[2];
			PS_SIT surrpix[2];
			pthread_mutex_lock (cdata_->wcsmutex);
			for (PS_SIT i=0; i<num_ndp; ++i)
			{
			    // get primary beam pixel coordinate for this sky plane pixel coorindate
			    // kind of a map between data WCS and primary beam WCS
			    img_pos_xy[0] = imgtps->xcoords[i];
			    img_pos_xy[1] = imgtps->ycoords[i];
			    HEADER getimgwcscoord (data_[g].wcs, data_[g].imgy, 
						   img_pos_xy[0], img_pos_xy[1], 
						   &img_pos_radec[0], &img_pos_radec[1]);
			    HEADER getimgpixcoord (wcs, imgy, cdata_,
						   img_pos_radec[0], img_pos_radec[1],
						   &pb_pos_xy[0], &pb_pos_xy[1]);

			    // find surrounding pixels and interpolate
			    surrpix[0] = (PS_SIT)std::floor (pb_pos_xy[0]);
			    surrpix[1] = (PS_SIT)std::floor (pb_pos_xy[1]);
			    pbeam[i] = 0;
			    for (PS_SIT x0=surrpix[0]; x0<=surrpix[0]+1; ++x0)
				if (x0>=0 && x0<=imgx-1)
				    for (PS_SIT y0=surrpix[1]; y0<=surrpix[1]+1; ++y0)
					if (y0>=0 && y0<=imgy-1)
					    pbeam[i] += std::abs(surrpix[0]-x0) *std::abs(surrpix[1]-y0)
						*imgdata[x0*imgy+y0];
			}
			pthread_mutex_unlock (cdata_->wcsmutex);

			// cleanup
			MEMORY ps_free (wcsinfo);
			MEMORY ps_free (invwcsinfo);
			MEMORY ps_free (imgdata);
			EXTERNAL ps_wcsfree (wcs);
			MEMORY ps_free (pbdata, nlines);
		    }

                    // get final scalar
                    sumsave[0] = sumsave[1] = 0;
                    sumsave[2] = data_[g].uv_newpixelscale;
                    for (PS_SIT mn=0; mn<data_[g].uv_ndp*2; ++mn)
                    {
                        sumsave[0] += data_[g].uv_data[mn]*data_[g].uv_data[mn] /data_[g].uv_newvariance[mn/2];
                        sumsave[1] += std::log (data_[g].uv_newvariance[mn/2]);
                    }

                    // get final vector
                    {
                        uv_transform_vec = new (uv_transform_vec_ptr)
                            VECTOR (cdata_, &data_[g], num_ndp+numextraRBF);
                        if (!data_[g].uv_padzeros)
                            imgtps->tpsmat->linequationsolve (uv_transform_vec, &qtcdinvd, 0);
                        else
                            imgtps->tpsmat->mult (&qtcdinvd, uv_transform_vec, 1, cdata_->numthreads);
                        uv_transform_vec->size -= numextraRBF;
			
			// apply primary beam
		    if (data_[g].uv_pbeam)
			for (PS_SIT n=0; n<uv_transform_vec->size; ++n)
			    uv_transform_vec->vec[n] *= pbeam[n];
                    }


                    // get final matrix
                    {
                        if (!data_[g].uv_padzeros)
                        {
                            // solve two linear equations, in place
                            // solution is -- r^-1 q^t cd ^-1 q r^-1
                            imgtps->tpsmat->inplace_leq (qtcdinvq);
                            // transpose to get -- q^t cd^-1 q r^-1
                            qtcdinvq->inplace_transpose ();
                            // now we get our answer -- r^-1 q^t cd^-1 q r^-1
                            imgtps->tpsmat->inplace_leq (qtcdinvq);
                        }
                        else
                        {
                            MATRIX *tmpmat  = new MATRIX (cdata_, &data_[g],
                                                          num_ndp+numextraRBF, num_ndp+numextraRBF,
                                                          -1, -1, &data_[g]);
                            imgtps->tpsmat->mult (qtcdinvq, tmpmat, 1, 0, cdata_->numthreads);
                            tmpmat->mult (imgtps->tpsmat, qtcdinvq, 0, 0, cdata_->numthreads);
                        }

                        // cleanup
                        MEMORY ps_free (imgtps->tpsmat);
                        imgtps->tpsmat = NULL;

                        uv_transform_mat = new (uv_transform_mat_ptr)
                            MATRIX (cdata_, &data_[g], 1, 1, -1, -1, &data_[g]);
                        MEMORY ps_free (uv_transform_mat->mat_dense);
                        uv_transform_mat->mat_dense = qtcdinvq->mat_dense;
                        uv_transform_mat->ncol = qtcdinvq->ncol;
                        uv_transform_mat->nrow = qtcdinvq->nrow;
                        uv_transform_mat->nnz  = qtcdinvq->nnz;
                        qtcdinvq->mat_dense = NULL;

                        // throw away zero-spacing uv coordinate
                        if (numextraRBF)
                            uv_transform_mat->remove_last_row_col ();
			
			// apply primary beam -- P^t (r^-1 q^t cd^-1 q r^-1) P
		    if (data_[g].uv_pbeam)
			for (PS_SIT m=0; m<num_ndp; ++m)
			    for (PS_SIT n=0; n<num_ndp; ++n)
				uv_transform_mat->mat_dense[n*num_ndp+m] 
				    *= pbeam[m]*pbeam[n];
                    }
		    
		    MEMORY ps_free (pbeam);

                    // print to file
                    if (1!=data_[g].verbose)
                        PRINTER print2screen (data_[g].print2screenname,
                                              "writing matrices to file",
                                              cdata_->print2screenmutex);
                    PS_SIT status = 0;
                    status += PRINTER printbinary <double> (-1, cdata_->basename, data_[g].name,
                                                            1, "uvmat1.bin",
                                                            uv_transform_mat->mat_dense,
                                                            uv_transform_mat->nnz, 0);
                    status += PRINTER printbinary <double> (-1, cdata_->basename, data_[g].name,
                                                            1, "uvmat2.bin",
                                                            uv_transform_vec->get_vec_ptr(),
                                                            uv_transform_vec->get_size(), 0);
                    status += PRINTER printbinary <double> (-1, cdata_->basename, data_[g].name,
                                                            1, "uvmat3.bin", sumsave, 3, 0);
                    if (!status)
                    {
                        PRINTER remove_file (-1, cdata_->basename, data_->name, 1,
                                             "snapshot."+OPERA tostring(lastuvblockind)+".stat", 0);
                        PRINTER remove_file (-1, cdata_->basename, data_->name, 1,
                                             "snapshot."+OPERA tostring(lastuvblockind)+".bin", 0);
                    }
                    else
                    {
                        PRINTER printwarning(data_[g].print2screenname,
                                             "error in writing final matrices to file",
                                             cdata_->print2screenmutex);
                    }
                }
                else
                {
                    // allocate memory for uv matrices
                    uv_transform_mat = new (uv_transform_mat_ptr)
                        MATRIX (cdata_, &data_[g], num_ndp, num_ndp, -1, -1, &data_[g]);
                    uv_transform_vec = new (uv_transform_vec_ptr)
                        VECTOR (cdata_, &data_[g], num_ndp);

                    // read in uv matrices
                    PRINTER readbinary (fn1, uv_transform_mat->mat_dense,
                                        uv_transform_mat->nnz, cdata_->print2screenmutex);
                    PRINTER readbinary (fn2, uv_transform_vec->vec,
                                        uv_transform_vec->get_size(), cdata_->print2screenmutex);
                    PRINTER readbinary (fn3, sumsave, 2, cdata_->print2screenmutex);
                }

                // store final stuff
                data_[g].uv_transform_scalar  = sumsave[0];
                data_[g].uv_transform_vec_ptr = (void*) uv_transform_vec_ptr;
                data_[g].uv_transform_mat_ptr = (void*) uv_transform_mat_ptr;
                data_[g].uv_transform_vec     = (void*) uv_transform_vec;
                data_[g].uv_transform_mat     = (void*) uv_transform_mat;

                // free some memory
                MEMORY ps_free (data_[g].uv_data);
                MEMORY ps_free (data_[g].uv_oldloc);
                MEMORY ps_free (data_[g].uv_newvariance);
                data_[g].uv_data = NULL;
                data_[g].uv_oldloc = NULL;
                data_[g].uv_newvariance = NULL;
                MEMORY ps_malloc (&data_[g].uv_newvariance, 1);
                data_[g].uv_newvariance[0] = sumsave[1];
            }
            else
            {
                // get weights of (Fourier transformed) TPS evaluated at uv plane positions
                MATRIX weights (cdata_, &data_[g], num_ndp+numextraRBF, data_[g].uv_ndp*2,
                                -1, -1, &data_[g]);
                imgtps->get_fourier_weights (data_[g].uv_oldloc, &weights, data_[g].uv_ndp);

                // store final matrix for repeated use
                MATRIX *uv_transform_ptr, *uv_transform;
                MEMORY ps_malloc (&uv_transform_ptr, 1);
                uv_transform = new (uv_transform_ptr)
                    MATRIX (cdata_, &data_[g], data_[g].uv_ndp*2, num_ndp, -1, -1, &data_[g]);
                weights.mult (imgtps->tpsmat, uv_transform, 1, 0, cdata_->numthreads);
                data_[g].uv_transform_mat_ptr = (void*) uv_transform_ptr;
                data_[g].uv_transform_mat     = (void*) uv_transform;
            }

            delete imgtps;
            delete qtcdinvq;
        }
    }
}

// comparator, for sorting uv data by u, then by v
struct uv_compare
{
    double *oldloc;

    uv_compare(double *oldloc_):
        oldloc(oldloc_) {}

    bool operator() (PS_SIT i, PS_SIT j)
        {
            if (oldloc[i*2]==oldloc[j*2])
                return oldloc[i*2+1]<oldloc[j*2+1];
            return oldloc[i*2]<oldloc[j*2];
        }
};

void pixsrc_init::uvdata_read (inputdata *data_, commoninputdata *cdata_)
{
    for (PS_SIT g=0; g<cdata_->numimages; ++g)
    {
        if (data_[g].is_uvdata)
        {
            if (1==data_[g].transmode)
            {
                // check if image to uv plane matrices already exist
                bool doexist = data_[g].uv_matvecsca_exist;

                if (doexist)
                {
                    data_[g].uv_data = NULL;
                    data_[g].uv_oldloc = NULL;
                    data_[g].uv_newvariance = NULL;
                    continue;
                }
            }

            PS_SIT gcount = 0;

            // read file
            char *fname;
            const char *listcc[2];
            listcc[0] = CONSTANT dir_in;
            listcc[1] = data_[g].uv_filename;
            OPERA concatenate( listcc, 2, &fname );

            // check if this is a binary file
            PS_SIT isbinary = 0;
            PS_SIT sizefn = OPERA sizestring (fname);
            if (sizefn>=4 && !strcmp (fname+sizefn-4,".bin"))
                isbinary = 1;

            if (isbinary)
            {
                // get number of uv points
                double numuvpts_dbl;
                PS_SIT numuvpts;
                PRINTER readbinary (fname, &numuvpts_dbl, 1, 0);
                data_[g].uv_ndp = data_[g].uv_ndporig = numuvpts = (PS_SIT)numuvpts_dbl;
		double *uvbinfile;

                MEMORY ps_malloc (&uvbinfile, 1+numuvpts*6);
                PRINTER readbinary (fname, uvbinfile, 1+numuvpts*6, 0);
                MEMORY ps_malloc (&data_[g].uv_oldloc,      numuvpts*2);
                MEMORY ps_malloc (&data_[g].uv_data,        numuvpts*2);
                MEMORY ps_malloc (&data_[g].uv_newvariance, numuvpts);
                MEMORY ps_malloc (&data_[g].uv_idx,         numuvpts);
                MEMORY ps_free (fname);

                // store data
                for (PS_SIT ind=0; ind<numuvpts; ++ind)
                {
		    data_[g].uv_idx[ind]         = ind;
                    data_[g].uv_oldloc[ind*2]    = uvbinfile[1+numuvpts*0+ind];
                    data_[g].uv_oldloc[ind*2+1]  = uvbinfile[1+numuvpts*1+ind];
                    data_[g].uv_data[ind*2]      = uvbinfile[1+numuvpts*2+ind];
                    data_[g].uv_data[ind*2+1]    = uvbinfile[1+numuvpts*3+ind];
                    data_[g].uv_newvariance[ind] = uvbinfile[1+numuvpts*4+ind]
                        *uvbinfile[1+numuvpts*4+ind];
                }
                MEMORY ps_free (uvbinfile);

                // apply uv cutoff
                if (-1!=data_[g].uv_cutoff[0] || -1!=data_[g].uv_cutoff[1])
                {
                    PS_SIT sortind = 0;
                    double cut[2] = {data_[g].uv_cutoff[0]*data_[g].uv_cutoff[0],
                                     data_[g].uv_cutoff[1]*data_[g].uv_cutoff[1]};

                    for (PS_SIT i=0; i<data_[g].uv_ndp; ++i)
                    {
                        double uv2 = data_[g].uv_oldloc[i*2+0] *data_[g].uv_oldloc[i*2+0]
                            + data_[g].uv_oldloc[i*2+1] *data_[g].uv_oldloc[i*2+1];
                        if (uv2>=cut[0] && uv2<=cut[1])
                        {
                            data_[g].uv_idx[sortind]         = data_[g].uv_idx[i];
                            data_[g].uv_oldloc[sortind*2+0]  = data_[g].uv_oldloc[i*2+0];
                            data_[g].uv_oldloc[sortind*2+1]  = data_[g].uv_oldloc[i*2+1];
                            data_[g].uv_data[sortind*2+0]    = data_[g].uv_data[i*2+0];
                            data_[g].uv_data[sortind*2+1]    = data_[g].uv_data[i*2+1];
                            data_[g].uv_newvariance[sortind] = data_[g].uv_newvariance[i];
                            ++sortind;
                        }
                    }
                    if (data_[g].verbose!=1)
                        PRINTER print2screen(data_[g].print2screenname,
                                             "removed "+ OPERA tostring(data_[g].uv_ndp-sortind) +
                                             " of " + OPERA tostring(data_[g].uv_ndp) +
                                             " uv data points due to uv cutoff",
                                             cdata_->print2screenmutex);
                    data_[g].uv_ndp = sortind;
                    MEMORY ps_realloc (&data_[g].uv_idx,         data_[g].uv_ndp);
                    MEMORY ps_realloc (&data_[g].uv_oldloc,      data_[g].uv_ndp*2);
                    MEMORY ps_realloc (&data_[g].uv_data,        data_[g].uv_ndp*2);
                    MEMORY ps_realloc (&data_[g].uv_newvariance, data_[g].uv_ndp);
                }
            }
            else
            {
                char **vec, ***vecsplit;
                PS_SIT vecsize,  *vecsplitsize;
                OPERA readfile( fname, &vec, &vecsize, cdata_->print2screenmutex);

                // these will contain parsed lines of vec
                MEMORY ps_malloc( &vecsplit    , vecsize );
                MEMORY ps_malloc( &vecsplitsize, vecsize );
		data_[g].uv_ndporig = vecsize;

                for( PS_SIT j=0; j<vecsize; ++j )
                {
                    OPERA split( vec[j], " ", &(vecsplit[j]), &(vecsplitsize[j]) );

                    if (5!=vecsplitsize[j])
                    {
                        PRINTER printerror("", "error in \"" +
                                           string(fname) + "\" file", cdata_->print2screenmutex);
                    }

                    double x = OPERA convert_string <double> (vecsplit[j][0]);
                    double y = OPERA convert_string <double> (vecsplit[j][1]);

                    // apply uv cutoff
                    if (
                        (-1==data_[g].uv_cutoff[1] ||
                         x*x + y*y <= data_[g].uv_cutoff[1]*data_[g].uv_cutoff[1])
                        &&
                        (-1==data_[g].uv_cutoff[0] ||
                         x*x + y*y >= data_[g].uv_cutoff[0]*data_[g].uv_cutoff[0])
                        )
                        ++gcount;
                }

                if (-1!=data_[g].uv_cutoff[0] || -1!=data_[g].uv_cutoff[1])
                    if (data_[g].verbose!=1)
                        PRINTER print2screen(data_[g].print2screenname,
                                             "removed "+OPERA tostring(vecsize-gcount) +
                                             " of " + OPERA tostring(vecsize) +
                                             " uv data points due to uv cutoff",
                                             cdata_->print2screenmutex);
                // allocate memory for data
                data_[g].uv_ndp = gcount;
                MEMORY ps_malloc (&data_[g].uv_oldloc,      gcount*2);
                MEMORY ps_malloc (&data_[g].uv_data,        gcount*2);
                MEMORY ps_malloc (&data_[g].uv_newvariance, gcount);
                MEMORY ps_malloc (&data_[g].uv_idx,         gcount);

                gcount = 0;
                for( PS_SIT j=0; j<vecsize; ++j )
                {
                    double x = OPERA convert_string <double> (vecsplit[j][0]);
                    double y = OPERA convert_string <double> (vecsplit[j][1]);
                    double r = OPERA convert_string <double> (vecsplit[j][2]);
                    double i = OPERA convert_string <double> (vecsplit[j][3]);
                    double nx = OPERA convert_string <double> (vecsplit[j][4]);

                    // apply uv cutoff
                    if (!(
                            (-1==data_[g].uv_cutoff[1] ||
                             x*x + y*y <= data_[g].uv_cutoff[1]*data_[g].uv_cutoff[1])
                            &&
                            (-1==data_[g].uv_cutoff[0] ||
                             x*x + y*y >= data_[g].uv_cutoff[0]*data_[g].uv_cutoff[0])
                            ))
                        continue;

		    data_[g].uv_idx[gcount]        = j;
		    data_[g].uv_oldloc[gcount*2]   = x;
                    data_[g].uv_oldloc[gcount*2+1] = y;
                    data_[g].uv_data[gcount*2]     = r;
                    data_[g].uv_data[gcount*2+1]   = i;
                    data_[g].uv_newvariance[gcount]   = nx*nx;
                    ++gcount;
                }

                // free data that was read in and split
                MEMORY ps_free( fname );
                for( PS_SIT j=0; j<vecsize; ++j )
                    MEMORY ps_free( vecsplit[j], vecsplitsize[j] );
                MEMORY ps_free (vecsplit    );
                MEMORY ps_free (vec, vecsize);
                MEMORY ps_free (vecsplitsize);
            }

            if (0!=data_[g].uv_taper)
            {
                double uv2, expfac;
                if (-1==data_[g].uv_taper)
                    expfac = CONSTANT twopi *CONSTANT pi
                        *data_[g].pix2arc *CONSTANT arc2rad
                        *data_[g].pix2arc *CONSTANT arc2rad;
                else
                    expfac = 1.0 /(data_[g].uv_taper *data_[g].uv_taper);

                if (data_[g].verbose!=1)
                    PRINTER print2screen(data_[g].print2screenname,
                                         "uv taper scale (lambda): " +
                                         OPERA tostring (std::sqrt(1.0/expfac)),
                                         cdata_->print2screenmutex);

                for (PS_SIT i=0; i<data_[g].uv_ndp; ++i)
                {
                    uv2 = OPERA distance2 (0, 0,
                                           data_[g].uv_oldloc[i*2+0],
                                           data_[g].uv_oldloc[i*2+1]);
                    data_[g].uv_newvariance[i] *= std::exp (uv2 *expfac);
                }
            }
	    
	    if (1==data_[g].transmode)
	    {
		PRINTER print2screen(data_[g].print2screenname,
				     "combining duplicate uv data points",
				     cdata_->print2screenmutex);
		// sort data by u and then v
		PS_SIT *tmpindices, *tmpvec_i;
		double *tmpvec;
		MEMORY ps_malloc (&tmpindices, data_[g].uv_ndp);
		MEMORY ps_malloc (&tmpvec,     data_[g].uv_ndp);
		MEMORY ps_malloc (&tmpvec_i,   data_[g].uv_ndp);
		// label indices of unsorted vector
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i)
		    tmpindices[i] = i;
		// sort indices by uv distance
		std::sort (tmpindices, tmpindices+data_[g].uv_ndp, uv_compare(data_[g].uv_oldloc));
		// sort data arrays
		// uv coordinate
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) tmpvec[i] = data_[g].uv_oldloc[i*2];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) data_[g].uv_oldloc[i*2] = tmpvec[tmpindices[i]];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) tmpvec[i] = data_[g].uv_oldloc[i*2+1];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) data_[g].uv_oldloc[i*2+1] = tmpvec[tmpindices[i]];
		// data
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) tmpvec[i] = data_[g].uv_data[i*2];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) data_[g].uv_data[i*2] = tmpvec[tmpindices[i]];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) tmpvec[i] = data_[g].uv_data[i*2+1];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) data_[g].uv_data[i*2+1] = tmpvec[tmpindices[i]];
		// noise
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) tmpvec[i] = data_[g].uv_newvariance[i];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) data_[g].uv_newvariance[i] = tmpvec[tmpindices[i]];
		// indexing
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) tmpvec_i[i] = data_[g].uv_idx[i];
		for (PS_SIT i=0; i<data_[g].uv_ndp; ++i) data_[g].uv_idx[i] = tmpvec_i[tmpindices[i]];
		
		// cleanup
		MEMORY ps_free (tmpindices);
		MEMORY ps_free (tmpvec_i);
		MEMORY ps_free (tmpvec);
		
		double *sortedoldloc, *sorteddata, *sortedvariance;
		PS_SIT *sortedidx;
		MEMORY ps_malloc (&sortedoldloc,   data_[g].uv_ndp*2);
		MEMORY ps_malloc (&sorteddata,     data_[g].uv_ndp*2);
		MEMORY ps_malloc (&sortedvariance, data_[g].uv_ndp);
		MEMORY ps_malloc (&sortedidx,      data_[g].uv_ndp);
		
		PS_SIT sindex = 0, numremoved=0, numduplicates=0;
		
		// get first index, if not a duplicate
		if (!(data_[g].uv_oldloc[0*2+0]==data_[g].uv_oldloc[1*2+0] &&
                          data_[g].uv_oldloc[0*2+1]==data_[g].uv_oldloc[1*2+1]))
		{
		    sortedoldloc[sindex*2+0] = data_[g].uv_oldloc[0*2+0];
		    sortedoldloc[sindex*2+1] = data_[g].uv_oldloc[0*2+1];
		    sorteddata[sindex*2+0]   = data_[g].uv_data[0*2+0];
		    sorteddata[sindex*2+1]   = data_[g].uv_data[0*2+1];
		    sortedvariance[sindex]   = data_[g].uv_newvariance[0];
		    sortedidx[sindex]        = data_[g].uv_idx[0];
		    ++sindex;
		}
		// remove duplicate uv points (inverse variance weighting)
		for (PS_SIT i=1; i<data_[g].uv_ndp; ++i)
		{
		    if (data_[g].uv_oldloc[i*2+0]==data_[g].uv_oldloc[(i-1)*2+0] &&
			data_[g].uv_oldloc[i*2+1]==data_[g].uv_oldloc[(i-1)*2+1])
		    {
			++numduplicates;
		    }
		    // if we found duplicates (but this uv point is not one of them)
		    else if (numduplicates)
		    {
			PS_SIT index1 = i-1-numduplicates;
			double brightness[2]={0,0}, invvarsum=0;
			for (PS_SIT j=0; j<numduplicates+1; ++j)
			{
			    brightness[0] += data_[g].uv_data[(index1+j)*2+0]
				/data_[g].uv_newvariance[index1+j];
			    brightness[1] += data_[g].uv_data[(index1+j)*2+1]
				/data_[g].uv_newvariance[index1+j];
			    invvarsum += 1.0 /data_[g].uv_newvariance[index1+j];
			}
			
			// set first data point to combined result
			sortedoldloc[sindex*2+0] = data_[g].uv_oldloc[index1*2+0];
			sortedoldloc[sindex*2+1] = data_[g].uv_oldloc[index1*2+1];
			sorteddata[sindex*2+0]   = brightness[0] /invvarsum;
			sorteddata[sindex*2+1]   = brightness[1] /invvarsum;
			sortedvariance[sindex]   = 1.0 /invvarsum;
			sortedidx[sindex]        = data_[g].uv_idx[index1];
			++sindex;
			
			// adjust some things
			numremoved += numduplicates;
			numduplicates = 0;
		    }
		    else
		    {
			sortedoldloc[sindex*2+0] = data_[g].uv_oldloc[i*2+0];
			sortedoldloc[sindex*2+1] = data_[g].uv_oldloc[i*2+1];
			sorteddata[sindex*2+0]   = data_[g].uv_data[i*2+0];
			sorteddata[sindex*2+1]   = data_[g].uv_data[i*2+1];
			sortedvariance[sindex]   = data_[g].uv_newvariance[i];
			sortedidx[sindex]        = data_[g].uv_idx[i];
			++sindex;
		    }
		}
		// check final block of uv data that was scanned, if duplicates exist
		if (numduplicates)
		{
		    PS_SIT index1 = data_[g].uv_ndp-1-numduplicates;
		    double brightness[2]={0,0}, invvarsum=0;
		    for (PS_SIT j=0; j<numduplicates+1; ++j)
		    {
			brightness[0] += data_[g].uv_data[(index1+j)*2+0]
			    /data_[g].uv_newvariance[index1+j];
			brightness[1] += data_[g].uv_data[(index1+j)*2+1]
			    /data_[g].uv_newvariance[index1+j];
			invvarsum += 1.0 /data_[g].uv_newvariance[index1+j];
		    }
		    
		    // set first data point to combined result
		    sortedoldloc[sindex*2+0] = data_[g].uv_oldloc[index1*2+0];
		    sortedoldloc[sindex*2+1] = data_[g].uv_oldloc[index1*2+1];
		    sorteddata[sindex*2+0]   = brightness[0] /invvarsum;
		    sorteddata[sindex*2+1]   = brightness[1] /invvarsum;
		    sortedvariance[sindex]   = 1.0 /invvarsum;
		    sortedidx[sindex]        = data_[g].uv_idx[index1];
		    ++sindex;
		    // adjust some things
		    numremoved += numduplicates;
		    numduplicates = 0;
		}
		
		data_[g].uv_ndp -= numremoved;
		PRINTER print2screen(data_[g].print2screenname,
				     "removed " + OPERA tostring(numremoved) +
				     " of " + OPERA tostring(data_[g].uv_ndp+numremoved) +
				     " duplicate uv data points using inverse variance weighting",
				     cdata_->print2screenmutex);
		if (sindex!=data_[g].uv_ndp)
		    PRINTER printerror(data_[g].print2screenname,
				       "internal error: uv data: remove duplicates: sindex!=data_->uv_ndp",
				       cdata_->print2screenmutex);
		
		MEMORY ps_free (data_[g].uv_oldloc);
		MEMORY ps_free (data_[g].uv_data);
		MEMORY ps_free (data_[g].uv_newvariance);
		MEMORY ps_free (data_[g].uv_idx);
		MEMORY ps_realloc (&sortedoldloc,   sindex*2);
		MEMORY ps_realloc (&sorteddata,     sindex*2);
		MEMORY ps_realloc (&sortedvariance, sindex);
		MEMORY ps_realloc (&sortedidx,      sindex);
		data_[g].uv_oldloc      = sortedoldloc;
		data_[g].uv_data        = sorteddata;
		data_[g].uv_newvariance = sortedvariance;
		data_[g].uv_idx         = sortedidx;
	    }
        }
    }
}




// for testing
// used in paper
/*
  PS_SIT ps_test_size=51;
  double ps_test_scale = 10;
  double ps_test_func (double x, double y)
  {
  PS_SIT size = ps_test_size;
  double scale = ps_test_scale;
  return std::exp(-std::pow(std::sqrt((x-size/2.)*(x-size/2.)
  +(y-size/2.)*(y-size/2.))/(size/1.5),1/.25))
  *(
  std::cos((x-size/2.)/scale*CONSTANT twopi) *std::sin((y-size/2.)/scale*CONSTANT twopi)
  *std::cos((x-size/2.)/scale/1.5*(y-size/2.)/scale/1.5*CONSTANT twopi)
  );
  return std::exp(-std::pow(std::sqrt((x-size/2.)*(x-size/2.)
  +(y-size/2.)*(y-size/2.))/(size/4.5),1/.25))
  *(
  std::cos((x-size/2.)/scale*CONSTANT twopi) *std::sin((y-size/2.)/scale*CONSTANT twopi)
  *std::cos((x-size/2.)/scale/1.5*(y-size/2.)/scale/1.5*CONSTANT twopi)
  );
  }
  double ps_test_rbf_integrate_y (double y, void *parms)
  {
  double x = *((double*)parms);
  return ps_test_func (x, y);
  }
  double ps_test_rbf_integrate_x (double x, void *parms)
  {
  double integral, interr;
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (10000);
  gsl_function F1;
  F1.function = &ps_test_rbf_integrate_y;
  F1.params = &x;
  PS_SIT y = *((PS_SIT*)parms);
  gsl_integration_qags (&F1, y-0.5, y+0.5, 0, 1.0e-6, 10000, ws, &integral, &interr);
  gsl_integration_workspace_free (ws);
  return integral;
  }
  void ps_test_tps (commoninputdata *cdata_, inputdata *data_)
  {
  // parms
  PS_SIT size = ps_test_size;
  PS_SIT size2 = size*size;
  double samplewidth=0.1;

  // utility variables
  PS_SIT r;
  double *imgpos, *img, err, maxerr=0, trueval, interpval;

  // create fake data -- 2d cosine
  MEMORY ps_malloc (&img, size2);
  MEMORY ps_malloc (&imgpos, size2*2);
  gsl_integration_workspace *ws = gsl_integration_workspace_alloc (10000);
  for (PS_SIT x=0; x<size; ++x)
  {
  for (PS_SIT y=0; y<size; ++y)
  {
  r = x*size + y;
  imgpos[r*2+0] = x;
  imgpos[r*2+1] = y;
  //img[r] = ps_test_func(x,y);
  //continue;
  //if (x==0&&y==0) {img[r]=0;continue;}
  double interr;
  gsl_function F1;
  F1.function = &ps_test_rbf_integrate_x;
  F1.params = &y;
  gsl_integration_qags (&F1, x-0.5, x+0.5, 0, 1.0e-6, 10000, ws, &img[r], &interr);
  }
  }

  TPS *imgtps;

  if (1)
  {
  PS_SIT nummore = 1;
  PS_SIT newnum = size2 + (4*4*size + 4*16)*nummore;
  double *newimg, *newimgpos;
  bool *ctrl;
  MEMORY ps_malloc (&ctrl, newnum);
  MEMORY ps_malloc (&newimg, newnum);
  MEMORY ps_malloc (&newimgpos, newnum*2);
  std::fill (ctrl, ctrl+size2, 1);
  std::fill (ctrl+size2, ctrl+newnum, 0);
  std::copy (img, img+size2, newimg);
  std::copy (imgpos, imgpos+size2*2, newimgpos);
  PS_SIT ind = size2;
  for (PS_SIT pp=0; pp<nummore; ++pp)
  {
  for (PS_SIT x=-4; x<0; ++x)
  {
  for (PS_SIT y=0; y<size; ++y)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT x=size; x<size+4; ++x)
  {
  for (PS_SIT y=0; y<size; ++y)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT y=-4; y<0; ++y)
  {
  for (PS_SIT x=0; x<size; ++x)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT y=size; y<size+4; ++y)
  {
  for (PS_SIT x=0; x<size; ++x)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT y=-4; y<0; ++y)
  {
  for (PS_SIT x=-4; x<0; ++x)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT y=-4; y<0; ++y)
  {
  for (PS_SIT x=size; x<size+4; ++x)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT y=size; y<size+4; ++y)
  {
  for (PS_SIT x=-4; x<0; ++x)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  for (PS_SIT y=size; y<size+4; ++y)
  {
  for (PS_SIT x=size; x<size+4; ++x)
  {
  newimgpos[ind*2+0] = x;
  newimgpos[ind*2+1] = y;
  newimg[ind] = 0;
  ++ind;
  }
  }
  }

  //for (PS_SIT b1=0; b1<newnum; ++b1)
  //    std::cout << "b1 " << newimgpos[b1*2] << " " << newimgpos[b1*2+1] << " " << newimg[b1] << std::endl;
  //exit(01);

  // do RBF fitting
  imgtps = new TPS(cdata_, &data_[0],
  newimgpos, newimgpos+1, newimg,

  2, 2, 1, newnum, ctrl, ps_rbf_gaussian, 1, 1, -1, 0);
  }
  else
  {
  // do RBF fitting
  imgtps = new TPS (cdata_, &data_[0],
  imgpos, imgpos+1, img,
  2, 2, 1, size2, ps_rbf_gaussian, 1, 1, -1, 0);
  }

  imgtps->get_tps_weights ();

  // get max error
  double pad = 5;
  for (double x=0-pad; x<size+pad; x+=samplewidth)
  {
  for (double y=0-pad; y<size+pad; y+=samplewidth)
  {
  if (x<0||x>size-1||y<0||y>size-1)
  trueval = 0;
  else
  trueval = ps_test_func (x,y);
  interpval = imgtps->interpolate (x,y);
  err = interpval-trueval;
  maxerr = std::abs(err)>std::abs(maxerr) ? err : maxerr;
  std::cout << std::setprecision(15) << "err: " << x << " " << y << " "
  << trueval << " " << interpval << " " << err << " " << (err-trueval)/trueval << std::endl;
  }
  }

  // talk to me
  std::cout << "Maximum error in TPS interpolation: " << maxerr << std::endl;
  exit(0);
  }
*/
