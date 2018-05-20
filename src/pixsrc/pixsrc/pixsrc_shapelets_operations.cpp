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



#include "pixsrc_shapelets_operations.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_common.hpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_printer.hpp"
#include "pixsrc_external.hpp"
#include <algorithm>
#include <cmath>
#include <gsl/gsl_sf_erf.h>
#include <unistd.h>

void pixsrc_shapelets_operations::createlo_shapelets (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    // if the user is fixing the shapelets center
    if (data_->fixedshapeletparms[0]==1 || data_->fixedshapeletparms[0]==3)
    {
        double pos[2];
        pthread_mutex_lock (cdata_->wcsmutex);
        HEADER getwcssfromwcsl (data_->fixedshapeletparms[1], data_->fixedshapeletparms[2],
                                data_->pra, data_->pdec,
                                data_->r1, data_->r2, pos);
        HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, pos[0], pos[1],
                               &vars_->shapelet_ctr[0], &vars_->shapelet_ctr[1]);
        pthread_mutex_unlock (cdata_->wcsmutex);
    }
    else
    {
        // estimate source center
        double sigfac = 3;
        double wt, weight, norm;
        do
        {
            norm = 0;
            vars_->shapelet_ctr[0] = vars_->shapelet_ctr[1] = 0;
            for( PS_SIT r=0; r<data_->ndp; ++r )
                if (-1!=vars_->r4r[r] && data_->data[r]>sigfac*std::sqrt(vars_->variance))
                {
                    wt = data_->data[r]*data_->data[r];
                    weight = wt*wt;
                    vars_->shapelet_ctr[0] += vars_->newloc[r*2  ] * weight;
                    vars_->shapelet_ctr[1] += vars_->newloc[r*2+1] * weight;
                    norm += weight;
                }
            vars_->shapelet_ctr[0] /= norm;
            vars_->shapelet_ctr[1] /= norm;
            sigfac *= 0.75;
        }
        while (norm==0);
    }

    // if the user is fixing the shapelet scale
    if (data_->fixedshapeletparms[0]==2 || data_->fixedshapeletparms[0]==3)
    {
        vars_->shapelet_scale = data_->fixedshapeletparms[3]*data_->arc2pix;
    }
    else
    {
        // estimate source scale
        double sigfac = 3;
        double weight, norm;
        do
        {
            norm = vars_->shapelet_scale = 0;
            for( PS_SIT r=0; r<data_->ndp; ++r )
                if (-1!=vars_->r4r[r] && data_->data[r]>sigfac*std::sqrt(vars_->variance))
                {
                    weight = std::sqrt(data_->data[r]*data_->data[r]);
                    weight *= weight;
                    vars_->shapelet_scale += OPERA distance
                        (vars_->shapelet_ctr[0], vars_->shapelet_ctr[1],
                         vars_->newloc[r*2], vars_->newloc[r*2+1])* weight;
                    norm += weight;
                }
            vars_->shapelet_scale /= norm;
            sigfac *= 0.75;
        }
        while (norm==0);
    }

    if (data_->use_shapelets==1)
        vars_->numberofshapelets = vars_->num_shapelets1*vars_->num_shapelets2;
    else if (data_->use_shapelets==2)
        vars_->numberofshapelets = vars_->num_shapelets1*(vars_->num_shapelets1+1)/2;
    vars_->lonc = vars_->numberofshapelets;

    vars_->shapeletopt = 1;

    if(data_->verbose==3)
    {
        double pos1[2], pos2[2];
        pthread_mutex_lock (cdata_->wcsmutex);
        HEADER getimgwcscoord (data_->wcs, data_->imgy, vars_->shapelet_ctr[0], vars_->shapelet_ctr[1],
                               &pos1[0], &pos1[1]);
        HEADER getwcslfromwcss (pos1[0], pos1[1], data_->pra, data_->pdec,
                                data_->r1, data_->r2, pos2);
        pthread_mutex_unlock (cdata_->wcsmutex);
        PRINTER print2screen(data_->print2screenname,
                             "shapelets center, scale: (" +
                             OPERA tostring(pos2[0]) + ", " +
                             OPERA tostring(pos2[1]) + "), " +
                             OPERA tostring (vars_->shapelet_scale*data_->pix2arc),
                             cdata_->print2screenmutex);
    }

    // get and set weights in lensing operator
    if (data_->use_shapelets==1)
    {
        SHAPELETSOPERA lo_set_hermite  (data_, cdata_, vars_);
    }
    else if (data_->use_shapelets==2)
    {
        PRINTER printerror(data_->print2screenname,
                           "polar shapelets disabled",
                           cdata_->print2screenmutex);
        //SHAPELETSOPERA lo_set_laguerre (data_, cdata_, vars_);
    }
}

struct herm_integrate
{
    inputdata *data_;
    commoninputdata *cdata_;
    lensvar *vars_;
    PS_SIT startr;
    double *hermvals;
    double *hermvals_split_tri;
};

void* pixsrc_shapelets_operations::lo_set_hermite_thread (void *args)
{
    // calculate hermite polynomials

    herm_integrate *pt_args    = (herm_integrate*)args;
    inputdata *data_           = pt_args->data_;
    commoninputdata *cdata_    = pt_args->cdata_;
    lensvar *vars_             = pt_args->vars_;
    PS_SIT startr                 = pt_args->startr;
    double *hermvals           = pt_args->hermvals;
    double *hermvals_split_tri = pt_args->hermvals_split_tri;

    PS_SIT rback;
    PS_SIT num_mag_steps;
    double mag_step_size, startpos[2];
    double tri[6], tricopy[6];
    double magxx, magxy, magyx, magyy, pot;
    double **newpos;

    for (PS_SIT r=startr; r<vars_->lonr; r+=cdata_->numthreads)
    {
        rback = vars_->r4rback[r];
        startpos[0] = data_->oldloc[rback*2]  -0.5;
        startpos[1] = data_->oldloc[rback*2+1]-0.5;

        // compute pixel splitting
        double logmag = std::max
            (0.0, std::log (vars_->magnification[rback])/std::log(4));
        num_mag_steps =
            OPERA round (std::pow (2, std::floor (logmag))) +
            data_->shapeletpixsplit[0]-1;
        num_mag_steps = std::min (data_->shapeletpixsplit[1],num_mag_steps);
        num_mag_steps = std::max (data_->shapeletpixsplit[0],num_mag_steps);

        mag_step_size = 1.0 / num_mag_steps;

        MEMORY ps_malloc (&newpos, 2, (num_mag_steps+1)*2);
        std::fill (hermvals, hermvals+vars_->numberofshapelets, 0);

        // get deflections for first TWO x positions
        pthread_mutex_lock (cdata_->potdefmagmutex);
        pthread_mutex_lock (cdata_->wcsmutex);
        // check if shapelet integration mask exists and apply it
        if (data_->shapeintmask)
        {
            for (PS_SIT d=0; d<2; ++d)
                std::fill (newpos[d], newpos[d]+(num_mag_steps+1)*2, -1e200);
            // loop over different masks and raytrace
            for(PS_SIT m=0; m<data_->extlengths[17]; ++m)
            {
                PS_SIT xistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][0]-startpos[0])/mag_step_size);
                PS_SIT xiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][2]-startpos[0])/mag_step_size);
                PS_SIT yistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][1]-startpos[1])/mag_step_size);
                PS_SIT yiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][3]-startpos[1])/mag_step_size);
                xistart = std::max (xistart,(PS_SIT)0);
                xiend   = std::min (xiend,(PS_SIT)2);
                yistart = std::max (yistart,(PS_SIT)0);
                yiend   = std::min (yiend,(PS_SIT)(num_mag_steps+1));

                for (PS_SIT xs=xistart; xs<xiend; ++xs)
                {
                    for (PS_SIT ys=yistart; ys<yiend; ++ys)
                    {
                        double thisx = startpos[0]+xs*mag_step_size;
                        double thisy = startpos[1]+ys*mag_step_size;
                        if (GEOM isinpoly (thisx,thisy,data_->shapeintmask[m],data_->extlengths[18]/2))
                        {
                            COMMON raytrace (data_, cdata_,
                                             thisx, thisy,
                                             &newpos[xs][ys*2], &newpos[xs][ys*2+1],
                                             &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
                        }
                    }
                }
            }
        }
        else
        {
            // get positions of all ray-traced triangles (for first TWO x coordinates)
            for (PS_SIT xs=0; xs<2; ++xs)
                for (PS_SIT ys=0; ys<num_mag_steps+1; ++ys)
                    COMMON raytrace (data_, cdata_,
                                     startpos[0]+xs*mag_step_size, startpos[1]+ys*mag_step_size,
                                     &newpos[xs][ys*2], &newpos[xs][ys*2+1],
                                     &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
        }
        pthread_mutex_unlock (cdata_->potdefmagmutex);
        pthread_mutex_unlock (cdata_->wcsmutex);

        PS_SIT minind, maxind;
        PS_SIT xind1=1, xind2=0;
        for (PS_SIT xs=0; xs<num_mag_steps; ++xs)
        {
            // as we loop over xs, the index of the 'newer' x coordinate
            // flips between 0 and 1.
            xind1 = 1-xind1;
            xind2 = 1-xind2;
            // only ray-trace another x-coordinate if we've looped through at least once
            if (xs)
            {
                pthread_mutex_lock (cdata_->potdefmagmutex);
                pthread_mutex_lock (cdata_->wcsmutex);
                // check if shapelet integration mask exists and apply it
                if (data_->shapeintmask)
                {
                    std::fill (newpos[xind2], newpos[xind2]+(num_mag_steps+1)*2, -1e200);
                    // loop over different masks and raytrace
                    for(PS_SIT m=0; m<data_->extlengths[17]; ++m)
                    {
                        PS_SIT xistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][0]-startpos[0])/mag_step_size);
                        PS_SIT xiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][2]-startpos[0])/mag_step_size);
                        PS_SIT yistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][1]-startpos[1])/mag_step_size);
                        PS_SIT yiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][3]-startpos[1])/mag_step_size);
                        xistart = std::max (xistart,(PS_SIT)(xs+1));
                        xiend   = std::min (xiend,(PS_SIT)(xs+2));
                        yistart = std::max (yistart,(PS_SIT)0);
                        yiend   = std::min (yiend,(PS_SIT)(num_mag_steps+1));

                        for (PS_SIT xss=xistart; xss<xiend; ++xss)
                        {
                            for (PS_SIT yss=yistart; yss<yiend; ++yss)
                            {
                                double thisx = startpos[0]+xss*mag_step_size;
                                double thisy = startpos[1]+yss*mag_step_size;
                                if (GEOM isinpoly (thisx,thisy,data_->shapeintmask[m],data_->extlengths[18]/2))
                                {
                                    COMMON raytrace (data_, cdata_,
                                                     thisx, thisy,
                                                     &newpos[xind2][yss*2], &newpos[xind2][yss*2+1],
                                                     &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
                                }
                            }
                        }
                    }
                }
                else
                {
                    // get positions of all ray-traced triangles (for this x coordinate)
                    for (PS_SIT xss=xs+1; xss<xs+2; ++xss)
                        for (PS_SIT yss=0; yss<num_mag_steps+1; ++yss)
                            COMMON raytrace (data_, cdata_,
                                             startpos[0]+xss*mag_step_size, startpos[1]+yss*mag_step_size,
                                             &newpos[xind2][yss*2], &newpos[xind2][yss*2+1],
                                             &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
                }
                pthread_mutex_unlock (cdata_->potdefmagmutex);
                pthread_mutex_unlock (cdata_->wcsmutex);
            }

            for (PS_SIT ys=0; ys<num_mag_steps; ++ys)
            {
                // get positions of lower triangle vertices
                tri[0] = (newpos[xind1][ys*2]       -vars_->shapelet_ctr[0]);
                tri[2] = (newpos[xind2][(ys+1)*2]   -vars_->shapelet_ctr[0]);
                tri[4] = (newpos[xind2][ys*2]       -vars_->shapelet_ctr[0]);
                tri[1] = (newpos[xind1][ys*2+1]     -vars_->shapelet_ctr[1]);
                tri[3] = (newpos[xind2][(ys+1)*2+1] -vars_->shapelet_ctr[1]);
                tri[5] = (newpos[xind2][ys*2+1]     -vars_->shapelet_ctr[1]);
                std::copy (tri, tri+6, tricopy);
                // sort x-coordinates
                minind = tri[0]     <tri[2] ? 0      : 2;
                minind = tri[minind]<tri[4] ? minind : 4;
                maxind = tri[0]     >tri[2] ? 0      : 2;
                maxind = tri[maxind]>tri[4] ? maxind : 4;
                tri[0] = tricopy[minind];
                tri[1] = tricopy[minind+1];
                tri[2] = tricopy[6-maxind-minind];
                tri[3] = tricopy[6-maxind-minind+1];
                tri[4] = tricopy[maxind];
                tri[5] = tricopy[maxind+1];
                SHAPELETSOPERA integrate_shapelets (data_, cdata_, vars_, tri, hermvals, hermvals_split_tri);

                // get positions of upper triangle vertices
                tri[0] = (newpos[xind1][ys*2]       -vars_->shapelet_ctr[0]);
                tri[2] = (newpos[xind2][(ys+1)*2]   -vars_->shapelet_ctr[0]);
                tri[4] = (newpos[xind1][(ys+1)*2]   -vars_->shapelet_ctr[0]);
                tri[1] = (newpos[xind1][ys*2+1]     -vars_->shapelet_ctr[1]);
                tri[3] = (newpos[xind2][(ys+1)*2+1] -vars_->shapelet_ctr[1]);
                tri[5] = (newpos[xind1][(ys+1)*2+1] -vars_->shapelet_ctr[1]);
                std::copy (tri, tri+6, tricopy);
                // sort x-coordinates
                minind = tri[0]     <tri[2] ? 0      : 2;
                minind = tri[minind]<tri[4] ? minind : 4;
                maxind = tri[0]     >tri[2] ? 0      : 2;
                maxind = tri[maxind]>tri[4] ? maxind : 4;
                tri[0] = tricopy[minind];
                tri[1] = tricopy[minind+1];
                tri[2] = tricopy[6-maxind-minind];
                tri[3] = tricopy[6-maxind-minind+1];
                tri[4] = tricopy[maxind];
                tri[5] = tricopy[maxind+1];
                SHAPELETSOPERA integrate_shapelets (data_, cdata_, vars_, tri, hermvals, hermvals_split_tri);
            }
        }

        for (PS_SIT s=0; s<vars_->numberofshapelets; ++s)
        {
            vars_->lensingoperatornobo->set (r, s, hermvals[s]/(num_mag_steps*num_mag_steps*2.0));
        }
        MEMORY ps_free (newpos, 2);
    }

    return NULL;
}

void pixsrc_shapelets_operations::lo_set_hermite (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    // structures for multi-threading
    pthread_t *herm_threads;
    herm_integrate *herm_structs;
    double **hermvals, **hermvals_split_tri;
    MEMORY ps_malloc (&herm_structs, cdata_->numthreads);
    MEMORY ps_malloc (&herm_threads, cdata_->numthreads);
    MEMORY ps_malloc (&hermvals, cdata_->numthreads, vars_->numberofshapelets);
    MEMORY ps_malloc (&hermvals_split_tri, cdata_->numthreads, vars_->numberofshapelets);

    // fill in structures and launch threads
    for (PS_SIT p=0; p<cdata_->numthreads; ++p)
    {
        herm_structs[p].data_              = data_;
        herm_structs[p].cdata_             = cdata_;
        herm_structs[p].vars_              = vars_;
        herm_structs[p].startr             = p;
        herm_structs[p].hermvals           = hermvals[p];
        herm_structs[p].hermvals_split_tri = hermvals_split_tri[p];

        pthread_create (&herm_threads[p] ,cdata_->attrjoinable,
                        SHAPELETSOPERA lo_set_hermite_thread, &herm_structs[p]);
    }

    // wait for threads to finish
    for (PS_SIT p=0; p<cdata_->numthreads; ++p)
        pthread_join (herm_threads[p], NULL);

    // cleanup
    MEMORY ps_free (herm_structs);
    MEMORY ps_free (herm_threads);
    MEMORY ps_free (hermvals, cdata_->numthreads);
    MEMORY ps_free (hermvals_split_tri, cdata_->numthreads);
}

void pixsrc_shapelets_operations::integrate_shapelets (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *tri, double *hermvals, double *hermvals_split_tri)
{

    if (data_->shapeintmask &&
        (tri[0]==-1e200||tri[1]==-1e200||tri[2]==-1e200||tri[3]==-1e200||tri[4]==-1e200||tri[5]==-1e200))
        return;

    // for each triangle, two integrations are performed
    // 1. from leftmost vertex to middle vertex
    // 2. from middle vertex to rightmost vertex
    double m1 = -OPERA slope (tri[0],tri[1],tri[4],tri[5]);
    double m2 = -OPERA slope (tri[0],tri[1],tri[2],tri[3]);
    double m3 = -OPERA slope (tri[2],tri[3],tri[4],tri[5]);
    if (std::abs(m1)<1e-5) m1 = m1<0 ? -1e-5 : 1e-5;
    if (std::abs(m2)<1e-5) m2 = m2<0 ? -1e-5 : 1e-5;
    if (std::abs(m3)<1e-5) m3 = m3<0 ? -1e-5 : 1e-5;

    double b1 = -OPERA intercept (tri[0],tri[1],tri[4],tri[5]);
    double b2 = -OPERA intercept (tri[0],tri[1],tri[2],tri[3]);
    double b3 = -OPERA intercept (tri[2],tri[3],tri[4],tri[5]);
    double intpoint[2] = {tri[2],m1*tri[2]+b1};
    double intfac = intpoint[1]>tri[3] ? 1 : -1;
    double triarea = 0;
    std::fill (hermvals_split_tri, hermvals_split_tri+vars_->numberofshapelets, 0);

    if (std::abs(m1)<1e5&&std::abs(m2)<1e5)
        SHAPELETSOPERA integrate_triangle (data_, cdata_, vars_, tri[0], intpoint[0], m2, m1, b2, b1, hermvals_split_tri, &triarea);
    if (std::abs(m1)<1e5&&std::abs(m3)<1e5)
        SHAPELETSOPERA integrate_triangle (data_, cdata_, vars_, intpoint[0], tri[4], m3, m1, b3, b1, hermvals_split_tri, &triarea);

    if (triarea)
        for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
            for (PS_SIT s2=0; s2<vars_->num_shapelets2; ++s2)
                hermvals[s1*vars_->num_shapelets2+s2] += intfac/triarea
                    *hermvals_split_tri[s1*vars_->num_shapelets2+s2];
}

void pixsrc_shapelets_operations::integrate_triangle (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double c, double d, double m1, double m2, double b1, double b2, double *hermvals_split_tri, double *triarea)
{
    PS_SIT max_num_sh = std::max (vars_->num_shapelets1, vars_->num_shapelets2);
    double scale_fac = std::pow (vars_->shapelet_scale, -0.5);
    double q = vars_->shapelet_scale;
    *triarea += std::abs (0.5*(m2-m1)*(d*d-c*c) + (b2-b1)*(d-c));

    double *phi_c, *phi_d, *phi_mcb1, *phi_mdb1, *phi_mcb2, *phi_mdb2;
    double **jnn1, **jnn2, **inn1, **inn2;
    MEMORY ps_malloc (&phi_c, max_num_sh);
    MEMORY ps_malloc (&phi_d, max_num_sh);
    MEMORY ps_malloc (&phi_mcb1, max_num_sh);
    MEMORY ps_malloc (&phi_mdb1, max_num_sh);
    MEMORY ps_malloc (&phi_mcb2, max_num_sh);
    MEMORY ps_malloc (&phi_mdb2, max_num_sh);
    MEMORY ps_malloc (&jnn1, vars_->num_shapelets1, vars_->num_shapelets2);
    MEMORY ps_malloc (&jnn2, vars_->num_shapelets1, vars_->num_shapelets2);
    MEMORY ps_malloc (&inn1, vars_->num_shapelets1, vars_->num_shapelets2);
    MEMORY ps_malloc (&inn2, vars_->num_shapelets1, vars_->num_shapelets2);

    // phi terms
    SHAPELETSOPERA hermite (max_num_sh-1, c/q, phi_c);
    SHAPELETSOPERA hermite (max_num_sh-1, d/q, phi_d);
    SHAPELETSOPERA hermite (max_num_sh-1, (m1*c+b1)/q, phi_mcb1);
    SHAPELETSOPERA hermite (max_num_sh-1, (m1*d+b1)/q, phi_mdb1);
    SHAPELETSOPERA hermite (max_num_sh-1, (m2*c+b2)/q, phi_mcb2);
    SHAPELETSOPERA hermite (max_num_sh-1, (m2*d+b2)/q, phi_mdb2);
    for (PS_SIT s1=0; s1<max_num_sh; ++s1)
    {
        phi_c[s1] *= CONSTANT shapelet_norm[s1] *scale_fac *std::exp (-c*c/(2*q*q));
        phi_d[s1] *= CONSTANT shapelet_norm[s1] *scale_fac *std::exp (-d*d/(2*q*q));
        phi_mcb1[s1] *= CONSTANT shapelet_norm[s1] *
            scale_fac *std::exp (-(m1*c+b1)*(m1*c+b1)/(2*q*q));
        phi_mdb1[s1] *= CONSTANT shapelet_norm[s1] *
            scale_fac *std::exp (-(m1*d+b1)*(m1*d+b1)/(2*q*q));
        phi_mcb2[s1] *= CONSTANT shapelet_norm[s1] *
            scale_fac *std::exp (-(m2*c+b2)*(m2*c+b2)/(2*q*q));
        phi_mdb2[s1] *= CONSTANT shapelet_norm[s1] *
            scale_fac *std::exp (-(m2*d+b2)*(m2*d+b2)/(2*q*q));
    }

    // j terms
    if (max_num_sh>1)
    {
        jnn1[0][0] = SHAPELETSOPERA j00 (m1,c,d,b1,q);
        jnn2[0][0] = SHAPELETSOPERA j00 (m2,c,d,b2,q);
    }
    if (vars_->num_shapelets2>1)
    {
        jnn1[0][1] = SHAPELETSOPERA j01 (m1,c,d,b1,q);
        jnn2[0][1] = SHAPELETSOPERA j01 (m2,c,d,b2,q);
    }
    if (vars_->num_shapelets1>1)
    {
        jnn1[1][0] = SHAPELETSOPERA j10 (m1,c,d,b1,q);
        jnn2[1][0] = SHAPELETSOPERA j10 (m2,c,d,b2,q);
    }
    for (PS_SIT s1=2; s1<vars_->num_shapelets1; ++s1)
    {
        jnn1[s1][0] = SHAPELETSOPERA jn0 (s1,c,d,m1,b1,q,jnn1,phi_c,phi_d);
        jnn2[s1][0] = SHAPELETSOPERA jn0 (s1,c,d,m2,b2,q,jnn2,phi_c,phi_d);
    }
    for (PS_SIT s1=2; s1<vars_->num_shapelets1; ++s1)
    {
        jnn1[0][s1] = SHAPELETSOPERA j0n (s1,c,d,m1,b1,q,jnn1,phi_mcb1,phi_mdb1);
        jnn2[0][s1] = SHAPELETSOPERA j0n (s1,c,d,m2,b2,q,jnn2,phi_mcb2,phi_mdb2);
    }
    if (vars_->num_shapelets1>1)
    {
        for (PS_SIT s1=1; s1<vars_->num_shapelets1-1; ++s1)
        {
            jnn1[1][s1] = SHAPELETSOPERA j1n (s1,m1,b1,q,jnn1);
            jnn2[1][s1] = SHAPELETSOPERA j1n (s1,m2,b2,q,jnn2);
        }
    }
    for (PS_SIT s1=2; s1<vars_->num_shapelets1; ++s1)
    {
        for (PS_SIT s2=1; s2<vars_->num_shapelets2; ++s2)
        {
            jnn1[s1][s2] = SHAPELETSOPERA jnxny (s1,s2,c,d,m1,b1,q,jnn1,
                                                 phi_c, phi_d, phi_mcb1,phi_mdb1);
            jnn2[s1][s2] = SHAPELETSOPERA jnxny (s1,s2,c,d,m2,b2,q,jnn2,
                                                 phi_c, phi_d, phi_mcb2,phi_mdb2);
        }
    }

    // i terms
    inn1[0][0] = SHAPELETSOPERA i00 (data_,cdata_, m1,c,d,b1,q);
    inn2[0][0] = SHAPELETSOPERA i00 (data_,cdata_, m2,c,d,b2,q);
    if (vars_->num_shapelets1>1)
    {
        inn1[1][0] = SHAPELETSOPERA i10 (m1,c,d,b1,q);
        inn2[1][0] = SHAPELETSOPERA i10 (m2,c,d,b2,q);
    }
    if (vars_->num_shapelets2>1)
    {
        for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
        {
            inn1[s1][1] = jnn1[s1][0];
            inn2[s1][1] = jnn2[s1][0];
        }
    }
    for (PS_SIT s1=0; s1<std::min((PS_SIT)2,vars_->num_shapelets1); ++s1)
    {
        for (PS_SIT s2=2; s2<vars_->num_shapelets2; ++s2)
        {
            inn1[s1][s2] = SHAPELETSOPERA inxny (s1,s2,inn1,jnn1);
            inn2[s1][s2] = SHAPELETSOPERA inxny (s1,s2,inn2,jnn2);
        }
    }
    for (PS_SIT s1=2; s1<vars_->num_shapelets1; ++s1)
    {
        inn1[s1][0] = SHAPELETSOPERA in0 (s1,c,d,m1,b1,q,inn1,phi_c,phi_d);
        inn2[s1][0] = SHAPELETSOPERA in0 (s1,c,d,m2,b2,q,inn2,phi_c,phi_d);
    }
    for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
    {
        for (PS_SIT s2=2; s2<vars_->num_shapelets2; ++s2)
        {
            inn1[s1][s2] = SHAPELETSOPERA inxny (s1,s2,inn1,jnn1);
            inn2[s1][s2] = SHAPELETSOPERA inxny (s1,s2,inn2,jnn2);
        }
    }

    for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
        for (PS_SIT s2=0; s2<vars_->num_shapelets2; ++s2)
        {
            hermvals_split_tri[s1*vars_->num_shapelets2+s2] += inn2[s1][s2]-inn1[s1][s2];
        }

    MEMORY ps_free (phi_c);
    MEMORY ps_free (phi_d);
    MEMORY ps_free (phi_mcb1);
    MEMORY ps_free (phi_mdb1);
    MEMORY ps_free (phi_mcb2);
    MEMORY ps_free (phi_mdb2);
    MEMORY ps_free (jnn1, vars_->num_shapelets1);
    MEMORY ps_free (jnn2, vars_->num_shapelets1);
    MEMORY ps_free (inn1, vars_->num_shapelets1);
    MEMORY ps_free (inn2, vars_->num_shapelets1);
}

double pixsrc_shapelets_operations::i00 (inputdata *data_, commoninputdata *cdata_, double m, double c, double d, double b, double q)
{
    double lower = (m*c+b)/(std::sqrt(2.0)*q);
    double upper = (m*d+b)/(std::sqrt(2.0)*q);

    if (lower<=0 && upper<=0)
        return -i00_calc (data_, cdata_, m,-upper, -lower,b,q,-1.0);
    else if (lower>=0 && upper>=0)
        return  i00_calc (data_, cdata_, m, lower, upper,b,q, 1.0);
    else
        return -i00_calc (data_, cdata_, m,0,-lower,b,q,-1.0) +
            i00_calc (data_, cdata_, m,0,upper,b,q,1.0);
}
double pixsrc_shapelets_operations::i00_calc (inputdata *data_, commoninputdata *cdata_, double m, double u1, double u2, double b, double q, double sign)
{
    if (CONSTANT shapelet_erf_approx_order==2)
    {
        // The following is working code, for a second order polynomial
        // fit to the erf function
        // a1,a2 are constants. d1-6 are fitted to the error function.
        // prefac multiplies the integral
        double a2 = 1.0/(m*m);
        double a1 = sign * std::sqrt(2.0)*b*a2/q;

        double d1 = -1.4618753656656642;
        double d2 = -0.86400152258553;
        double d3 = 0.33424917797640685;
        double d4 = 0.27570810605268997;
        double expbarg = -b*b*a2/(2.0*q*q);
        double prefac = q/m/8.0;

        // some commonly occuring terms
        double SPI = std::sqrt (CONSTANT pi);

        double A12 = a1*a1;
        double SA2 = std::sqrt (a2);
        double A1PD1 = a1 + d1;
        double A1PD12 = A1PD1*A1PD1;
        double A2MD2 = a2 - d2;
        double A2MD2_12 = std::sqrt (A2MD2);
        double A2MD2_M52 = std::pow (A2MD2, -5.0/2.0);

        double expfac1_1 = 4.0*SPI*gsl_sf_erf(2*a2*u1-a1)/SA2;
        double expfac2_1 = A2MD2_M52*2.0*A2MD2_12*(A1PD1*d4+2.0*A2MD2*(d3+d4*u1));
        double expfac3_1 = -A2MD2_M52*(2.0*A2MD2*(2*A2MD2+A1PD1*d3)+(A1PD12+2.0*A2MD2)*d4)*SPI*gsl_sf_erf((2*A2MD2*u1-A1PD1)/(2.0*A2MD2_12));
        double expfac1_2 = 4.0*SPI*gsl_sf_erf(2*a2*u2-a1)/SA2;
        double expfac2_2 = A2MD2_M52*2.0*A2MD2_12*(A1PD1*d4+2.0*A2MD2*(d3+d4*u2));
        double expfac3_2 = -A2MD2_M52*(2.0*A2MD2*(2*A2MD2+A1PD1*d3)+(A1PD12+2.0*A2MD2)*d4)*SPI*gsl_sf_erf((2*A2MD2*u2-A1PD1)/(2.0*A2MD2_12));
        expfac1_1 *= -1;
        expfac2_1 *= -1;
        expfac3_1 *= -1;

        // arguments of the exponentials
        double exparg1 = A12/(4.0*a2) + expbarg;
        double exparg2_1 = u1*(A1PD1-A2MD2*u1) + expbarg;
        double exparg2_2 = u2*(A1PD1-A2MD2*u2) + expbarg;
        double exparg3 = A1PD12/(4.0*A2MD2) + expbarg;

        // add exponentials, while trying to avoid overflow (for small slopes m)
        double maxarg = std::max (exparg1,std::max(exparg3,std::max(exparg2_1,exparg2_2)));
        double logarg =
            (expfac1_1+expfac1_2) * std::exp (exparg1-maxarg) +
            (expfac3_1+expfac3_2) * std::exp (exparg3-maxarg) +
            expfac2_1 * std::exp (exparg2_1-maxarg) +
            expfac2_2 * std::exp (exparg2_2-maxarg);

        double sgnmult = logarg>0 ? 1 : (logarg<0 ? -1 : 0);

        return sgnmult*prefac *
            std::exp (std::log (sgnmult*logarg) + maxarg);
    }
    else if (CONSTANT shapelet_erf_approx_order!=2)
    {
        // The following is working code, for a fourth order polynomial
        // fit to the erf function
        // a1,a2 are constants. d1-6 are fitted to the error function.
        // prefac multiplies the integral
        double a2 = 1.0/(m*m);
        double a1 = sign * std::sqrt(2.0)*b*a2/q;

        double d1 = -1.6093661701173445;
        double d2 = -0.9144569742603699;
        double d3 = 0.48105948656562003;
        double d4 = 0.3927744807707903;
        double d5 = 0.05250001720254402;
        double d6 = 0.036174940195734834;
        double expbarg = -b*b*a2/(2.0*q*q);
        double prefac = q/m;

        // some commonly occuring terms
        double SPI = std::sqrt (CONSTANT pi);
        double D1D2 = d1*d2;
        double D1D4 = d1*d4;
        double D1D5 = d1*d5;
        double D1D3 = d1*d3;
        double D1D6 = d1*d6;
        double D2D3 = d2*d3;
        double D2D4 = d2*d4;
        double D2D5 = d2*d5;
        double D2D6 = d2*d6;
        double D12 = d1*d1;
        double D22 = d2*d2;
        double D23 = D22*d2;
        double D13 = D12*d1;

        double A12 = a1*a1;
        double A22 = a2*a2;
        double A23 = A22*a2;
        double A13 = A12*a1;
        double A1PD1 = a1 + d1;
        double A1PD12 = A1PD1*A1PD1;
        double A1PD13 = A1PD12*A1PD1;
        double A2MD2 = a2 - d2;
        double A2MD22 = A2MD2*A2MD2;
        double A2MD24 = A2MD22*A2MD22;
        double A2MD292 = std::pow (A2MD2,9.0/2.0);
        double LONG2 = A12*A12*d6 + D12*D12*d6;

        // u terms for lower (1) and upper (2) bounds
        double U2_1 = u1*u1;
        double D5PD6U_1 = d5 + d6*u1;
        double D4Pd5pd6uU_1 = d4 + D5PD6U_1*u1;
        double LONG1_1 = a1*D4Pd5pd6uU_1 + d1*D4Pd5pd6uU_1 + 2*d5 + 3*d6*u1;
        double LONG3_1 = (a1*D5PD6U_1 + d1*D5PD6U_1 + 5*d6)*A1PD1;
        double LONG4_1 = d3 + D4Pd5pd6uU_1*u1;
        double U2_2 = u2*u2;
        double D5PD6U_2 = d5 + d6*u2;
        double D4Pd5pd6uU_2 = d4 + D5PD6U_2*u2;
        double LONG1_2 = a1*D4Pd5pd6uU_2 + d1*D4Pd5pd6uU_2 + 2*d5 + 3*d6*u2;
        double LONG3_2 = (a1*D5PD6U_2 + d1*D5PD6U_2 + 5*d6)*A1PD1;
        double LONG4_2 = d3 + D4Pd5pd6uU_2*u2;

        // factors multiplying the exponentials
        double expfac1_1 = (2*(A1PD13*d6 + 4*D22*LONG1_1 - 2*d2*LONG3_1 + 8*A23*LONG4_1 - 8*D23*LONG4_1 + 2*a2*(-4*d2*LONG1_1 + LONG3_1 + 12*D22*LONG4_1) + 4*A22*(D1D4 + a1*D4Pd5pd6uU_1 + 2*d5 - 6*d2*LONG4_1 + D1D5*u1 + 3*d6*u1 + D1D6*U2_1)))/A2MD24/32;
        double expfac2_1 = -((12*D1D5*D22 - 8*D1D3*D23 - 2*D13*D2D5 - 12*D12*D2D6 + 4*D12*D22*d4 - 8*D23*d4 + 8*A23*(D1D3 + d4) + 2*A13*(2*D1D6 - D2D5 + a2*d5) + 12*D22*d6 + 4*A22*(3*D1D5 - 6*D2D4 - 6*D1D2*d3 + D12*d4 + 3*d6) + 2*a2*(12*D1D2*(D2D3 - d5) + D13*d5 + 12*d2*(D2D4 - d6) + D12*(-4*D2D4 + 6*d6)) + 2*A12*(2*A22*d4 + 2*D22*d4 + 3*D12*d6 - 3*d2*(D1D5 + 2*d6) + a2*(3*D1D5 - 4*D2D4 + 6*d6)) + 2*a1*(4*A23*d3 - 4*D23*d3 + D22*(4*D1D4 + 6*d5) + A22*(4*D1D4 - 12*D2D3 + 6*d5) + 2*D13*d6 - 3*D1D2*(D1D5 + 4*d6) + a2*(12*D22*d3 - 4*d2*(2*D1D4 + 3*d5) + 3*d1*(D1D5 + 4*d6))) + LONG2)*SPI*gsl_sf_erf((-a1 - d1 + 2*A2MD2*u1)/(2*std::sqrt(A2MD2))))/A2MD292/32;
        double expfac3_1 = 16*SPI*(gsl_sf_erf((-a1 + 2*a2*u1)/(2*std::sqrt(a2))))/std::sqrt(a2)/32;
        double expfac4_1 = -16*SPI*gsl_sf_erf((-a1 - d1 + 2*a2*u1 - 2*d2*u1)/(2*std::sqrt(A2MD2)))/std::sqrt(A2MD2)/32;
        double expfac1_2 = (2*(A1PD13*d6 + 4*D22*LONG1_2 - 2*d2*LONG3_2 + 8*A23*LONG4_2 - 8*D23*LONG4_2 + 2*a2*(-4*d2*LONG1_2 + LONG3_2 + 12*D22*LONG4_2) + 4*A22*(D1D4 + a1*D4Pd5pd6uU_2 + 2*d5 - 6*d2*LONG4_2 + D1D5*u2 + 3*d6*u2 + D1D6*U2_2)))/A2MD24/32;
        double expfac2_2 = -((12*D1D5*D22 - 8*D1D3*D23 - 2*D13*D2D5 - 12*D12*D2D6 + 4*D12*D22*d4 - 8*D23*d4 + 8*A23*(D1D3 + d4) + 2*A13*(2*D1D6 - D2D5 + a2*d5) + 12*D22*d6 + 4*A22*(3*D1D5 - 6*D2D4 - 6*D1D2*d3 + D12*d4 + 3*d6) + 2*a2*(12*D1D2*(D2D3 - d5) + D13*d5 + 12*d2*(D2D4 - d6) + D12*(-4*D2D4 + 6*d6)) + 2*A12*(2*A22*d4 + 2*D22*d4 + 3*D12*d6 - 3*d2*(D1D5 + 2*d6) + a2*(3*D1D5 - 4*D2D4 + 6*d6)) + 2*a1*(4*A23*d3 - 4*D23*d3 + D22*(4*D1D4 + 6*d5) + A22*(4*D1D4 - 12*D2D3 + 6*d5) + 2*D13*d6 - 3*D1D2*(D1D5 + 4*d6) + a2*(12*D22*d3 - 4*d2*(2*D1D4 + 3*d5) + 3*d1*(D1D5 + 4*d6))) + LONG2)*SPI*gsl_sf_erf((-a1 - d1 + 2*A2MD2*u2)/(2*std::sqrt(A2MD2))))/A2MD292/32;
        double expfac3_2 = 16*SPI*(gsl_sf_erf((-a1 + 2*a2*u2)/(2*std::sqrt(a2))))/std::sqrt(a2)/32;
        double expfac4_2 = -16*SPI*gsl_sf_erf((-a1 - d1 + 2*a2*u2 - 2*d2*u2)/(2*std::sqrt(A2MD2)))/std::sqrt(A2MD2)/32;
        expfac1_1 *= -1;
        expfac2_1 *= -1;
        expfac3_1 *= -1;
        expfac4_1 *= -1;

        // arguments of the exponentials
        double exparg1_1 = u1*(A1PD1 + (-a2 + d2)*u1) + expbarg;
        double exparg1_2 = u2*(A1PD1 + (-a2 + d2)*u2) + expbarg;
        double exparg2 = A1PD12/(4*A2MD2) + expbarg;
        double exparg3 = A12/(4*a2) + expbarg;
        double exparg4 = A1PD12/(4*A2MD2) + expbarg;

        // add exponentials, while trying to avoid overflow (for small slopes m)
        double maxarg = std::max (exparg1_1,std::max(exparg1_2,std::max(exparg2,std::max(exparg3,exparg4))));
        double logarg =
            expfac1_1*std::exp (exparg1_1-maxarg) +
            expfac1_2*std::exp (exparg1_2-maxarg) +
            (expfac2_1+expfac2_2)*std::exp (exparg2-maxarg) +
            (expfac3_1+expfac3_2)*std::exp (exparg3-maxarg) +
            (expfac4_1+expfac4_2)*std::exp (exparg4-maxarg);
        double sgnmult = logarg>0 ? 1 : (logarg<0 ? -1 : 0);

        return sgnmult*prefac *
            std::exp (std::log (sgnmult*logarg) + maxarg);
    }
    else return 0;
}

double pixsrc_shapelets_operations::j00 (double m, double x1, double x2, double b, double q)
{
    //
    // original function:
    //return -std::exp(-b*b/(2.0*(1.0+m*m)*q*q)) * q/std::sqrt(1.0+m*m) *
    //gsl_sf_erf((x+m*(b+m*x))/(std::sqrt(2.0)*std::sqrt(1+m*m)*q));
    //
    // the following prevents overflow
    //

    double expfac1 = -q/std::sqrt(1.0+m*m) *
        gsl_sf_erf((x1+m*(b+m*x1))/(std::sqrt(2.0)*std::sqrt(1+m*m)*q));
    double expfac2 = -q/std::sqrt(1.0+m*m) *
        gsl_sf_erf((x2+m*(b+m*x2))/(std::sqrt(2.0)*std::sqrt(1+m*m)*q));
    double exparg = -b*b/(2.0*(1.0+m*m)*q*q);
    expfac1 *= -1;

    double logarg = expfac1 + expfac2;
    double sgnmult = logarg>0 ? 1 : (logarg<0 ? -1 : 0);

    return sgnmult * std::exp (std::log (sgnmult*logarg) + exparg);
}

double pixsrc_shapelets_operations::j10 (double m, double x1, double x2, double b, double q)
{
    /*
    // original function:
    return (2.0*std::sqrt(1.0+m*m)*q+b*
    std::exp(std::pow(x+m*(b+m*x),2.0)/(2.0*(1.0+m*m)*q*q))*m*
    std::sqrt(2.0*CONSTANT pi)*
    gsl_sf_erf((x+m*(b+m*x))/(std::sqrt(2.0)*std::sqrt(1.0+m*m)*q)))/
    (std::exp((b*b+2.0*b*m*x+(1.0+m*m)*x*x)/(2.0*q*q))*
    std::pow(1.0+m*m,3.0/2.0)*std::sqrt(CONSTANT pi));
    //
    // the following prevents overflow
    */

    double expfac1_1 = 2.0*std::sqrt(1.0+m*m)*q /
        (std::pow(1.0+m*m,3.0/2.0)*std::sqrt(CONSTANT pi));
    double expfac2_1 = b*m*std::sqrt(2.0*CONSTANT pi)*
        gsl_sf_erf((x1+m*(b+m*x1))/(std::sqrt(2.0)*std::sqrt(1.0+m*m)*q)) /
        (std::pow(1.0+m*m,3.0/2.0)*std::sqrt(CONSTANT pi));
    double expfac1_2 = expfac1_1;
    double expfac2_2 = b*m*std::sqrt(2.0*CONSTANT pi)*
        gsl_sf_erf((x2+m*(b+m*x2))/(std::sqrt(2.0)*std::sqrt(1.0+m*m)*q)) /
        (std::pow(1.0+m*m,3.0/2.0)*std::sqrt(CONSTANT pi));
    double exparg1_1 = -(b*b+2.0*b*m*x1+(1.0+m*m)*x1*x1)/(2.0*q*q);
    double exparg2_1 = exparg1_1 + std::pow(x1+m*(b+m*x1),2.0)/(2.0*(1.0+m*m)*q*q);
    double exparg1_2 = -(b*b+2.0*b*m*x2+(1.0+m*m)*x2*x2)/(2.0*q*q);
    double exparg2_2 = exparg1_2 + std::pow(x2+m*(b+m*x2),2.0)/(2.0*(1.0+m*m)*q*q);
    expfac1_1 *= -1;
    expfac2_1 *= -1;

    double maxarg = std::max (std::max (std::max(exparg1_1,exparg2_1),exparg1_2),exparg2_2);
    double logarg =
        expfac1_1*std::exp (exparg1_1-maxarg) +
        expfac1_2*std::exp (exparg1_2-maxarg) +
        expfac2_1*std::exp (exparg2_1-maxarg) +
        expfac2_2*std::exp (exparg2_2-maxarg);
    double sgnmult = logarg>0 ? 1 : (logarg<0 ? -1 : 0);

    return sgnmult * std::exp (std::log (sgnmult*logarg) + maxarg);
}

double pixsrc_shapelets_operations::j01 (double m, double x1, double x2, double b, double q)
{
    /*
    // original funtion:
    return 1.0/(m*m+1.0)*m*q*std::sqrt(2.0/CONSTANT pi)*
    std::exp(-((m*x+b)*(m*x+b)+x*x)/2.0/(q*q)) -
    std::pow(m*m+1.0,-3.0/2.0)*b*std::exp(-b*b/(2*(m*m+1.0)*q*q))*
    gsl_sf_erf(((m*m+1.0)*x+m*b)/q/std::sqrt(2.0*(m*m+1.0)));
    //
    // the following prevents overflow
    */

    double expfac1_1 = 1.0/(m*m+1.0)*m*q*std::sqrt(2.0/CONSTANT pi);
    double expfac2_1 = -std::pow(m*m+1.0,-3.0/2.0)*b*
        gsl_sf_erf(((m*m+1.0)*x1+m*b)/q/std::sqrt(2.0*(m*m+1.0)));
    double expfac1_2 = expfac1_1;
    double expfac2_2 = -std::pow(m*m+1.0,-3.0/2.0)*b*
        gsl_sf_erf(((m*m+1.0)*x2+m*b)/q/std::sqrt(2.0*(m*m+1.0)));
    double exparg1_1 = -((m*x1+b)*(m*x1+b)+x1*x1)/2.0/(q*q);
    double exparg2_1 = -b*b/(2*(m*m+1.0)*q*q);
    double exparg1_2 = -((m*x2+b)*(m*x2+b)+x2*x2)/2.0/(q*q);
    double exparg2_2 = exparg2_1;
    expfac1_1 *= -1;
    expfac2_1 *= -1;

    double maxarg = std::max (std::max (std::max(exparg1_1,exparg2_1),exparg1_2),exparg2_2);
    double logarg =
        expfac1_1*std::exp (exparg1_1-maxarg) +
        expfac1_2*std::exp (exparg1_2-maxarg) +
        expfac2_1*std::exp (exparg2_1-maxarg) +
        expfac2_2*std::exp (exparg2_2-maxarg);
    double sgnmult = logarg>0 ? 1 : (logarg<0 ? -1 : 0);

    return sgnmult * std::exp (std::log (sgnmult*logarg) + maxarg);
}

double pixsrc_shapelets_operations::jn0 (PS_SIT nx, double c, double d, double m, double b, double q, double **jnn, double *phi_c, double *phi_d)
{
    return jn0_uv(nx,m,d,b,q,phi_d) - jn0_uv(nx,m,c,b,q,phi_c) +
        jn0_vdu(nx,m,b,q,jnn);
}
double pixsrc_shapelets_operations::jn0_uv (PS_SIT nx, double m, double x, double b, double q, double *phi)
{
    return 1.0/(std::sqrt(nx)*(m*m+1.0)) *
        (2.0*std::pow(q,3.0/2.0)*
         std::pow(CONSTANT pi,-1.0/4.0)*phi[nx-1]*
         std::exp(-(m*x+b)*(m*x+b)/(2.0*q*q)));
}
double pixsrc_shapelets_operations::jn0_vdu (PS_SIT nx, double m, double b, double q, double **jnn)
{
    return 1.0/(std::sqrt(nx)*(m*m+1.0)) *
        (-std::sqrt(2.0)*b*m/q*jnn[nx-1][0] -
         (m*m-1.0)*std::sqrt(nx-1.0)*jnn[nx-2][0]);
}

double pixsrc_shapelets_operations::j0n (PS_SIT ny, double c, double d, double m, double b, double q, double **jnn, double *phi_mcb, double *phi_mdb)
{
    return j0n_uv(ny,m,d,b,q,phi_mdb) - j0n_uv(ny,m,c,b,q,phi_mcb) +
        j0n_vdu(ny,m,b,q,jnn);
}
double pixsrc_shapelets_operations::j0n_uv (PS_SIT ny, double m, double x, double b, double q, double *phi)
{
    return 1.0/(std::sqrt(ny+1.0)*(m*m+1.0)) *
        (2*m*std::pow(q,3.0/2.0)*std::pow(CONSTANT pi,-1.0/4.0)/
         std::sqrt(ny)*phi[ny-1]*std::exp(-x*x/(2.0*q*q)));
}
double pixsrc_shapelets_operations::j0n_vdu (PS_SIT ny, double m, double b, double q, double **jnn)
{
    return 1.0/(std::sqrt(ny+1.0)*(m*m+1.0)) *
        (std::sqrt(2.0)*b/q*jnn[0][ny-1] +
         (m*m-1.0)*(ny-1.0)/std::sqrt(ny)*jnn[0][ny-2]);
}

double pixsrc_shapelets_operations::j1n (PS_SIT ny, double m, double b, double q, double **jnn)
{
    return (std::sqrt(ny+2.0)*jnn[0][ny+1] +
            ny/std::sqrt(ny+1.0)*jnn[0][ny-1] -
            std::sqrt(2.0)*b/q*jnn[0][ny])/m;
}

double pixsrc_shapelets_operations::jnxny (PS_SIT nx, PS_SIT ny, double c, double d, double m, double b, double q, double **jnn, double *phi_c, double *phi_d, double *phi_mcb, double *phi_mdb)
{
    return jnxny_uv(nx,ny,m,d,b,q,phi_d,phi_mdb) - jnxny_uv(nx,ny,m,c,b,q,phi_c,phi_mcb) +
        jnxny_vdu(nx,ny,m,b,q,jnn);
}
double pixsrc_shapelets_operations::jnxny_uv (PS_SIT nx, PS_SIT ny, double m, double x, double b, double q, double *phi, double *phi_mxb)
{
    return 1.0/(std::sqrt(nx)*(m*m+1.0))*
        (2*q*q/std::sqrt(ny+1.0)*phi[nx-1]*phi_mxb[ny]);
}
double pixsrc_shapelets_operations::jnxny_vdu (PS_SIT nx, PS_SIT ny, double m, double b, double q, double **jnn)
{
    return 1.0/(std::sqrt(nx)*(m*m+1.0))*
        (-b*m/q*std::sqrt(2)*jnn[nx-1][ny] -
         (m*m-1.0)*std::sqrt(nx-1.0)*jnn[nx-2][ny] +
         2.0*m*ny/std::sqrt(ny+1.0)*jnn[nx-1][ny-1]);
}

double pixsrc_shapelets_operations::in0 (PS_SIT nx, double c, double d, double m, double b, double q, double **inn, double *phi_c, double *phi_d)
{
    return in0_uv(nx,m,d,b,q,phi_d) - in0_uv(nx,m,c,b,q,phi_c) +
        in0_vdu(nx,m,b,q,inn);
}
double pixsrc_shapelets_operations::in0_uv (PS_SIT nx, double m, double x, double b, double q, double *phi)
{
    return 1.0/std::sqrt(nx)*
        (-std::pow(CONSTANT pi,1.0/4.0)*std::pow(q,3.0/2.0)*
         gsl_sf_erf((m*x+b)/(q*std::sqrt(2.0)))*phi[nx-1]);
}
double pixsrc_shapelets_operations::in0_vdu (PS_SIT nx, double m, double b, double q, double **inn)
{
    return 1.0/std::sqrt(nx)*
        (-m*inn[nx-1][1] + std::sqrt(nx-1.0)*inn[nx-2][0]);
}

double pixsrc_shapelets_operations::inxny (PS_SIT nx, PS_SIT ny, double **inn, double **jnn)
{
    return std::sqrt((ny-1.0)/ny)*inn[nx][ny-2] + jnn[nx][ny-1];
}

double pixsrc_shapelets_operations::i10 (double m, double x1, double x2, double b, double q)
{
    /*
    // original function:
    return (-q*gsl_sf_erf((b + m*x)/(std::sqrt(2.0)*q))*std::exp(-x*x/(2.0*q*q)) +
    q*(m*gsl_sf_erf((x + m*(b + m*x))/(std::sqrt(2.0)*std::sqrt(1.0 + m*m)*q)))*
    std::exp(-b*b/(2.0*(1.0 + m*m)*q*q)) / std::sqrt(1.0 + m*m));
    //
    // the following prevents overflow
    */

    double expfac1_1 = -q*gsl_sf_erf((b + m*x1)/(std::sqrt(2.0)*q));
    double expfac2_1 = gsl_sf_erf((x1 + m*(b + m*x1))/
                                  (std::sqrt(2.0)*std::sqrt(1.0 + m*m)*q))*
        q*m/std::sqrt(1.0 + m*m);
    double expfac1_2 = -q*gsl_sf_erf((b + m*x2)/(std::sqrt(2.0)*q));
    double expfac2_2 = gsl_sf_erf((x2 + m*(b + m*x2))/
                                  (std::sqrt(2.0)*std::sqrt(1.0 + m*m)*q))*
        q*m/std::sqrt(1.0 + m*m);
    double exparg1_1 = -x1*x1/(2.0*q*q);
    double exparg2_1 = -b*b/(2.0*(1.0 + m*m)*q*q);
    double exparg1_2 = -x2*x2/(2.0*q*q);
    double exparg2_2 = exparg2_1;
    expfac1_1 *= -1;
    expfac2_1 *= -1;

    double maxarg = std::max (std::max (std::max(exparg1_1,exparg2_1),exparg1_2),exparg2_2);
    double logarg =
        expfac1_1*std::exp (exparg1_1-maxarg) +
        expfac1_2*std::exp (exparg1_2-maxarg) +
        expfac2_1*std::exp (exparg2_1-maxarg) +
        expfac2_2*std::exp (exparg2_2-maxarg);
    double sgnmult = logarg>0 ? 1 : (logarg<0 ? -1 : 0);

    return sgnmult * std::exp (std::log (sgnmult*logarg) + maxarg);
}

/*
  void pixsrc_shapelets_operations::lo_set_laguerre (inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
  {
  // calculate laguerre polynomials
  PS_SIT rback, norm_index, sind, s2;
  double pos[2], rad, rad2, theta, val;
  double **laguerre;
  MEMORY ps_malloc (&laguerre,
  (vars_->num_shapelets1-1)/2+1, vars_->num_shapelets1);

  for (PS_SIT r=0; r<vars_->lonr; ++r)
  {
  rback = vars_->r4rback[r];
  if (vars_->pointer[rback]==-1)
  continue;

  pos[0] = (vars_->newloc[rback*2]  -vars_->shapelet_ctr[0]) / vars_->shapelet_scale;
  pos[1] = (vars_->newloc[rback*2+1]-vars_->shapelet_ctr[1]) / vars_->shapelet_scale;

  rad2 = OPERA distance2 (0,0,pos[0],pos[1]);
  rad  = std::sqrt (rad2);
  theta = std::atan2 (pos[1], pos[0]);

  for (s2=0; s2<vars_->num_shapelets1; ++s2)
  SHAPELETSOPERA laguerre (vars_->num_shapelets1-1, s2, rad2, laguerre);

  norm_index = sind = 0;
  for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
  {
  for (PS_SIT s2_0=0; s2_0<=s1; ++s2_0)
  {
  s2 = -s1+2*s2_0;
  if (s2<0)
  continue;
  if (rad>CONSTANT laguerre_max[sind] || rad<CONSTANT laguerre_min[sind])
  {
  ++norm_index;
  sind = s2==0 ? sind+1 : sind+2;
  continue;
  }

  val = 2 * CONSTANT laguerre_norm[norm_index] *
  std::pow (vars_->shapelet_scale,-(s2+1)) *
  std::pow (rad,s2) * laguerre[(s1-s2)/2][s2] *
  std::exp (-rad2/(2*vars_->shapelet_scale*vars_->shapelet_scale));
  ++norm_index;

  // real part
  vars_->lensingoperatornobo->set (r, sind, val*std::cos(s2*theta));
  ++sind;
  if (s2==0)
  continue;
  // imaginary part
  vars_->lensingoperatornobo->set (r, sind, val*std::sin(s2*theta));
  ++sind;
  }
  }
  }

  MEMORY ps_free (laguerre, (vars_->num_shapelets1-1)/2+1);
  }
*/

void pixsrc_shapelets_operations::hermite (PS_SIT n, double x, double *res)
{
    res[0] = 1;

    if (!n)
        return;

    res[1] = 2*x;

    for (PS_SIT i=2; i<=n; ++i)
    {
        res[i] = 2*x*res[i-1] - 2*(i-1)*res[i-2];
    }
}

double pixsrc_shapelets_operations::get_hermite_sb (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double x, double y)
{
    double q = vars_->shapelet_scale;
    double sb = 0;
    double pos[2] = {(x - vars_->shapelet_ctr[0]),
                     (y - vars_->shapelet_ctr[1])};
    double expf = std::exp (-(pos[0]*pos[0]+pos[1]*pos[1])/(2.0*q*q));
    double scale_fac = 1.0/q;

    double *herm1, *herm2;
    MEMORY ps_malloc (&herm1, vars_->num_shapelets1);
    MEMORY ps_malloc (&herm2, vars_->num_shapelets2);

    SHAPELETSOPERA hermite (vars_->num_shapelets1-1, pos[0]/q, herm1);
    SHAPELETSOPERA hermite (vars_->num_shapelets2-1, pos[1]/q, herm2);

    for (PS_SIT i1=0; i1<vars_->num_shapelets1; ++i1)
    {
        for (PS_SIT i2=0; i2<vars_->num_shapelets2; ++i2)
        {
            sb +=
                vars_->mps->get(i1*vars_->num_shapelets2+i2) *
                expf * scale_fac *
                CONSTANT shapelet_norm[i1] * herm1[i1] *
                CONSTANT shapelet_norm[i2] * herm2[i2];
        }
    }

    MEMORY ps_free (herm1);
    MEMORY ps_free (herm2);

    return sb;
}

double pixsrc_shapelets_operations::integrate_square (inputdata *data_, commoninputdata *cdata_, lensvar *vars_,
                                                      double *hermvals_square, PS_SIT numsh_1, PS_SIT numsh_2,
                                                      double sh_ctrx, double sh_ctry, double sh_scale,
                                                      double a, double b, double c, double d,
                                                      VECTOR *mps, double *noise)
{
    // This function computes all the integrals for a particular square.
    // It would be more efficient to separately compute the x and y integrals and
    // then multiply them afterwards for each cell/pixel.
    // But, I don't expect this function to be called repeatedly for multiple pixsrc runs.
    // Note: This function does not divide by the area of the square.
    // a and b are lower and upper x limits
    // c and d are lower and upper y limits
    // mps are shapelet weights (for computing source plane flux) and can be NULLS
    // noise is not NULL if you want to compute noise in this pixel (source plane)
    // hermvals_square is an array that will have shapelet weights upon exit (cannot be null)

    // setup
    double scale_fac = std::pow (sh_scale, -0.5);
    double q = sh_scale;
    PS_SIT numshape[2] = {numsh_1,numsh_2};
    double limits[4] = {a-sh_ctrx,b-sh_ctrx,c-sh_ctry,d-sh_ctry};
    double *phi_l, *phi_u, **integrals, *lim;
    MEMORY ps_malloc (&phi_l, numshape[0]);
    MEMORY ps_malloc (&phi_u, numshape[1]);
    MEMORY ps_malloc (&integrals, 2, std::max(numshape[0],numshape[1]));
    double flux = 0;

    // loop over x and y axes to get 1d integrals of shapelets
    for (PS_SIT ax=0; ax<2; ++ax)
    {
        // set integral limits
        lim = !ax ? limits : limits+2;

        // get shapelets evaluated at endpoints
        SHAPELETSOPERA hermite (numshape[ax]-1, lim[0]/q, phi_l);
        SHAPELETSOPERA hermite (numshape[ax]-1, lim[1]/q, phi_u);
        for (PS_SIT s=0; s<numshape[ax]; ++s)
        {
            phi_l[s] *= CONSTANT shapelet_norm[s] *scale_fac *std::exp (-lim[0]*lim[0]/(2*q*q));
            phi_u[s] *= CONSTANT shapelet_norm[s] *scale_fac *std::exp (-lim[1]*lim[1]/(2*q*q));
        }

        // get first two shapelet integrals
        integrals[ax][0] = std::sqrt(q*CONSTANT sqrtpi*0.5) *
            (gsl_sf_erf(lim[1]/(q*std::sqrt(2.0))) - gsl_sf_erf(lim[0]/(q*std::sqrt(2.0))));
        if (numshape[ax]>1)
            integrals[ax][1] = -q*std::sqrt(2.0) * (phi_u[0]-phi_l[0]);

        // get higher order integrals
        for (PS_SIT s=2; s<numshape[ax]; ++s)
            integrals[ax][s] = std::sqrt((s-1.0)/s)*integrals[ax][s-2]
                - q*std::sqrt(2.0/s)*(phi_u[s-1]-phi_l[s-1]);
    }

    // compute products of 1d shapelets to get 2d integrals
    for (PS_SIT s1=0; s1<numshape[0]; ++s1)
        for (PS_SIT s2=0; s2<numshape[1]; ++s2)
        {
            hermvals_square[s1*numshape[1]+s2] = integrals[0][s1]*integrals[1][s2];
            if (mps)
                flux +=  hermvals_square[s1*numshape[1]+s2] * mps->get (s1*numshape[1]+s2);
        }

    if (data_ && cdata_ && vars_ && data_->noisemap && noise)
    {
        VECTOR *dummy1, *dummy2, *product, *hermvec;
        MEMORY ps_malloc (&dummy1, 1);
        MEMORY ps_malloc (&dummy2, 1);
        product = new (dummy1) VECTOR (cdata_, data_, vars_->numberofshapelets);
        hermvec = new (dummy2) VECTOR (cdata_, data_, vars_->numberofshapelets);

        // copy hermvals
        for (PS_SIT s=0; s<vars_->numberofshapelets; ++s)
            hermvec->set (s,integrals[0][s/numshape[1]]*integrals[1][s%numshape[1]]);

        // multiply covariance matrix by shapelet weights
        vars_->a1->inv_dense_mult (hermvec, product, cdata_->numthreads);
        // get variance and standard deviation
        *noise = product->innerproduct (hermvec);
        *noise = std::sqrt (*noise);

        product->~VECTOR();
        hermvec->~VECTOR();
        MEMORY ps_free (dummy1);
        MEMORY ps_free (dummy2);
    }

    MEMORY ps_free (phi_l);
    MEMORY ps_free (phi_u);
    MEMORY ps_free (integrals, 2);

    return flux;
}

void pixsrc_shapelets_operations::get_binomials (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double *binom, PS_SIT max_num)
{
    // calculate square root of binomial coefficients for even shapelets

    binom[0] = 1;
    for (PS_SIT s=2; s<max_num; s+=2)
    {
        binom[s/2] = 1;
        for (PS_SIT i=1; i<=s/2; ++i)
            binom[s/2] *= (s+1-i)/(double)i;
        binom[s/2] = std::sqrt (binom[s/2]);
    }
}

double pixsrc_shapelets_operations::get_flux_one_pixel (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, PS_SIT rback, double **hermvecs, VECTOR *thismps)
{
    // get flux in one data pixel (r)
    // this function is used when computing magnification
    // It calculates flux in pixels that were masked out

    if (!thismps)
        thismps = vars_->mps;

    double sb = 0;
    double *hermvals           = hermvecs[0];
    double *hermvals_split_tri = hermvecs[1];

    PS_SIT num_mag_steps;
    double mag_step_size, startpos[2];
    double tri[6], tricopy[6];
    double magxx, magxy, magyx, magyy, pot;
    double **newpos;

    startpos[0] = data_->oldloc[rback*2]  -0.5;
    startpos[1] = data_->oldloc[rback*2+1]-0.5;

    // compute pixel splitting
    double logmag = std::max
        (0.0, std::log (vars_->magnification[rback])/std::log(4));
    num_mag_steps =
        OPERA round (std::pow (2, std::floor (logmag))) +
        data_->shapeletpixsplit[0]-1;
    num_mag_steps = std::min (data_->shapeletpixsplit[1],num_mag_steps);
    num_mag_steps = std::max (data_->shapeletpixsplit[0],num_mag_steps);

    mag_step_size = 1.0 / num_mag_steps;

    MEMORY ps_malloc (&newpos, 2, (num_mag_steps+1)*2);
    std::fill (hermvals, hermvals+vars_->numberofshapelets, 0);

    // get deflections for first TWO x positions
    pthread_mutex_lock (cdata_->potdefmagmutex);
    pthread_mutex_lock (cdata_->wcsmutex);
    // check if shapelet integration mask exists and apply it
    if (data_->shapeintmask)
    {
        for (PS_SIT d=0; d<2; ++d)
            std::fill (newpos[d], newpos[d]+(num_mag_steps+1)*2, -1e200);
        // loop over different masks and raytrace
        for(PS_SIT m=0; m<data_->extlengths[17]; ++m)
        {
            PS_SIT xistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][0]-startpos[0])/mag_step_size);
            PS_SIT xiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][2]-startpos[0])/mag_step_size);
            PS_SIT yistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][1]-startpos[1])/mag_step_size);
            PS_SIT yiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][3]-startpos[1])/mag_step_size);
            xistart = std::max (xistart,(PS_SIT)0);
            xiend   = std::min (xiend,(PS_SIT)2);
            yistart = std::max (yistart,(PS_SIT)0);
            yiend   = std::min (yiend,(PS_SIT)(num_mag_steps+1));

            for (PS_SIT xs=xistart; xs<xiend; ++xs)
            {
                for (PS_SIT ys=yistart; ys<yiend; ++ys)
                {
                    double thisx = startpos[0]+xs*mag_step_size;
                    double thisy = startpos[1]+ys*mag_step_size;
                    if (GEOM isinpoly (thisx,thisy,data_->shapeintmask[m],data_->extlengths[18]/2))
                    {
                        COMMON raytrace (data_, cdata_,
                                         thisx, thisy,
                                         &newpos[xs][ys*2], &newpos[xs][ys*2+1],
                                         &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
                    }
                }
            }
        }
    }
    else
    {
        // get positions of all ray-traced triangles (for first TWO x coordinates)
        for (PS_SIT xs=0; xs<2; ++xs)
            for (PS_SIT ys=0; ys<num_mag_steps+1; ++ys)
                COMMON raytrace (data_, cdata_,
                                 startpos[0]+xs*mag_step_size, startpos[1]+ys*mag_step_size,
                                 &newpos[xs][ys*2], &newpos[xs][ys*2+1],
                                 &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
    }
    pthread_mutex_unlock (cdata_->potdefmagmutex);
    pthread_mutex_unlock (cdata_->wcsmutex);

    PS_SIT minind, maxind;
    PS_SIT xind1=1, xind2=0;
    for (PS_SIT xs=0; xs<num_mag_steps; ++xs)
    {
        // as we loop over xs, the index of the 'newer' x coordinate
        // flips between 0 and 1.
        xind1 = 1-xind1;
        xind2 = 1-xind2;
        // only ray-trace another x-coordinate if we've looped through at least once
        if (xs)
        {
            pthread_mutex_lock (cdata_->potdefmagmutex);
            pthread_mutex_lock (cdata_->wcsmutex);
            // check if shapelet integration mask exists and apply it
            if (data_->shapeintmask)
            {
                std::fill (newpos[xind2], newpos[xind2]+(num_mag_steps+1)*2, -1e200);
                // loop over different masks and raytrace
                for(PS_SIT m=0; m<data_->extlengths[17]; ++m)
                {
                    PS_SIT xistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][0]-startpos[0])/mag_step_size);
                    PS_SIT xiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][2]-startpos[0])/mag_step_size);
                    PS_SIT yistart = (PS_SIT)std::floor ((data_->shapeintmaskminmax[m][1]-startpos[1])/mag_step_size);
                    PS_SIT yiend   = (PS_SIT)std::ceil  ((data_->shapeintmaskminmax[m][3]-startpos[1])/mag_step_size);
                    xistart = std::max (xistart,(PS_SIT)(xs+1));
                    xiend   = std::min (xiend,(PS_SIT)(xs+2));
                    yistart = std::max (yistart,(PS_SIT)0);
                    yiend   = std::min (yiend,(PS_SIT)(num_mag_steps+1));

                    for (PS_SIT xss=xistart; xss<xiend; ++xss)
                    {
                        for (PS_SIT yss=yistart; yss<yiend; ++yss)
                        {
                            double thisx = startpos[0]+xss*mag_step_size;
                            double thisy = startpos[1]+yss*mag_step_size;
                            if (GEOM isinpoly (thisx,thisy,data_->shapeintmask[m],data_->extlengths[18]/2))
                            {
                                COMMON raytrace (data_, cdata_,
                                                 thisx, thisy,
                                                 &newpos[xind2][yss*2], &newpos[xind2][yss*2+1],
                                                 &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
                            }
                        }
                    }
                }
            }
            else
            {
                // get positions of all ray-traced triangles (for this x coordinate)
                for (PS_SIT xss=xs+1; xss<xs+2; ++xss)
                    for (PS_SIT yss=0; yss<num_mag_steps+1; ++yss)
                        COMMON raytrace (data_, cdata_,
                                         startpos[0]+xss*mag_step_size, startpos[1]+yss*mag_step_size,
                                         &newpos[xind2][yss*2], &newpos[xind2][yss*2+1],
                                         &magxx, &magxy, &magyx, &magyy, &pot, -1, vars_->imagenumber);
            }
            pthread_mutex_unlock (cdata_->potdefmagmutex);
            pthread_mutex_unlock (cdata_->wcsmutex);
        }

        for (PS_SIT ys=0; ys<num_mag_steps; ++ys)
        {
            // get positions of lower triangle vertices
            tri[0] = (newpos[xind1][ys*2]       -vars_->shapelet_ctr[0]);
            tri[2] = (newpos[xind2][(ys+1)*2]   -vars_->shapelet_ctr[0]);
            tri[4] = (newpos[xind2][ys*2]       -vars_->shapelet_ctr[0]);
            tri[1] = (newpos[xind1][ys*2+1]     -vars_->shapelet_ctr[1]);
            tri[3] = (newpos[xind2][(ys+1)*2+1] -vars_->shapelet_ctr[1]);
            tri[5] = (newpos[xind2][ys*2+1]     -vars_->shapelet_ctr[1]);
            std::copy (tri, tri+6, tricopy);
            // sort x-coordinates
            minind = tri[0]     <tri[2] ? 0      : 2;
            minind = tri[minind]<tri[4] ? minind : 4;
            maxind = tri[0]     >tri[2] ? 0      : 2;
            maxind = tri[maxind]>tri[4] ? maxind : 4;
            tri[0] = tricopy[minind];
            tri[1] = tricopy[minind+1];
            tri[2] = tricopy[6-maxind-minind];
            tri[3] = tricopy[6-maxind-minind+1];
            tri[4] = tricopy[maxind];
            tri[5] = tricopy[maxind+1];
            SHAPELETSOPERA integrate_shapelets (data_, cdata_, vars_, tri, hermvals, hermvals_split_tri);

            // get positions of upper triangle vertices
            tri[0] = (newpos[xind1][ys*2]       -vars_->shapelet_ctr[0]);
            tri[2] = (newpos[xind2][(ys+1)*2]   -vars_->shapelet_ctr[0]);
            tri[4] = (newpos[xind1][(ys+1)*2]   -vars_->shapelet_ctr[0]);
            tri[1] = (newpos[xind1][ys*2+1]     -vars_->shapelet_ctr[1]);
            tri[3] = (newpos[xind2][(ys+1)*2+1] -vars_->shapelet_ctr[1]);
            tri[5] = (newpos[xind1][(ys+1)*2+1] -vars_->shapelet_ctr[1]);
            std::copy (tri, tri+6, tricopy);
            // sort x-coordinates
            minind = tri[0]     <tri[2] ? 0      : 2;
            minind = tri[minind]<tri[4] ? minind : 4;
            maxind = tri[0]     >tri[2] ? 0      : 2;
            maxind = tri[maxind]>tri[4] ? maxind : 4;
            tri[0] = tricopy[minind];
            tri[1] = tricopy[minind+1];
            tri[2] = tricopy[6-maxind-minind];
            tri[3] = tricopy[6-maxind-minind+1];
            tri[4] = tricopy[maxind];
            tri[5] = tricopy[maxind+1];
            SHAPELETSOPERA integrate_shapelets (data_, cdata_, vars_, tri, hermvals, hermvals_split_tri);
        }
    }

    for (PS_SIT s=0; s<vars_->numberofshapelets; ++s)
    {
        sb += hermvals[s]/(num_mag_steps*num_mag_steps*2.0) * thismps->get(s);
    }
    MEMORY ps_free (newpos, 2);

    return sb;
}

/*
  double pixsrc_shapelets_operations::get_laguerre_sb (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double x, double y)
  {
  PS_SIT s2, norm_index=0, sind=0;
  double val, sb = 0;
  double pos[2] = {(x - vars_->shapelet_ctr[0]) / vars_->shapelet_scale,
  (y - vars_->shapelet_ctr[1]) / vars_->shapelet_scale};

  double rad2 = OPERA distance2 (0,0,pos[0],pos[1]);
  double rad  = std::sqrt (rad2);
  double theta = std::atan2 (pos[1], pos[0]);

  double **laguerre;
  MEMORY ps_malloc (&laguerre,
  (vars_->num_shapelets1-1)/2+1, vars_->num_shapelets1);
  for (s2=0; s2<vars_->num_shapelets1; ++s2)
  SHAPELETSOPERA laguerre (vars_->num_shapelets1-1, s2, rad2, laguerre);

  for (PS_SIT s1=0; s1<vars_->num_shapelets1; ++s1)
  {
  for (PS_SIT s2_0=0; s2_0<=s1; ++s2_0)
  {
  s2 = -s1+2*s2_0;
  if (s2<0)
  continue;
  if (rad>CONSTANT laguerre_max[sind] || rad<CONSTANT laguerre_min[sind])
  {
  ++norm_index;
  sind = s2==0 ? sind+1 : sind+2;
  continue;
  }

  val = 2 * CONSTANT laguerre_norm[norm_index] *
  std::pow (vars_->shapelet_scale,-(s2+1)) *
  std::pow (rad,s2) * laguerre[(s1-s2)/2][s2] *
  std::exp (-rad2/(2*vars_->shapelet_scale*vars_->shapelet_scale));
  ++norm_index;

  sb += val * std::cos(s2*theta) * vars_->mps->get(sind++);
  if (s2!=0)
  sb += val * std::sin(s2*theta) * vars_->mps->get(sind++);
  }
  }

  MEMORY ps_free (laguerre, (vars_->num_shapelets1-1)/2+1);

  return sb;
  }
*/

/*
  void pixsrc_shapelets_operations::laguerre (PS_SIT n0, PS_SIT alpha, double x, double **res)
  {
  PS_SIT n = (n0-alpha)/2;

  res[0][alpha] = 1;
  if (n==0)
  return;

  res[1][alpha] = 1+alpha-x;
  for (PS_SIT i=2; i<=n; ++i)
  res[i][alpha] = ( (2*i-1+alpha-x)*res[i-1][alpha] -
  (i-1+alpha)*res[i-2][alpha] ) / i;
  }
*/

/*
  void pixsrc_shapelets_operations::get_hermite_sb (inputdata *data_, commoninputdata *cdata_, lensvar *vars_, double x, double y, double div, double *sb, double *noise)
  {
  PS_SIT num_mag_steps = 1;
  double mag_step_size = 1.0/div;
  double startpos[2] = {x-1.0/(div*2.0), y-1.0/(div*2.0)};
  double tri[6];
  double newpos[2][4];
  double *hermvals, *hermvals_split_tri;
  *sb = 0;

  MEMORY ps_malloc (&hermvals, vars_->numberofshapelets);
  MEMORY ps_malloc (&hermvals_split_tri, vars_->numberofshapelets);
  std::fill (hermvals, hermvals+vars_->numberofshapelets, 0);

  // get positions of all ray-traced triangles
  for (PS_SIT xs=0; xs<num_mag_steps+1; ++xs)
  for (PS_SIT ys=0; ys<num_mag_steps+1; ++ys)
  {
  newpos[xs][ys*2]   = startpos[0]+xs*mag_step_size;
  newpos[xs][ys*2+1] = startpos[1]+ys*mag_step_size;
  }

  PS_SIT xs=0, ys=0;
  // get positions of lower triangle vertices
  tri[0] = (newpos[xs][ys*2]-vars_->shapelet_ctr[0]);
  tri[2] = (newpos[xs+1][(ys+1)*2]-vars_->shapelet_ctr[0]);
  tri[4] = (newpos[xs+1][ys*2]-vars_->shapelet_ctr[0]);
  tri[1] = (newpos[xs][ys*2+1]-vars_->shapelet_ctr[1]);
  tri[3] = (newpos[xs+1][(ys+1)*2+1]-vars_->shapelet_ctr[1]);
  tri[5] = (newpos[xs+1][ys*2+1]-vars_->shapelet_ctr[1]);
  SHAPELETSOPERA integrate_shapelets (data_, cdata_, vars_, tri, hermvals, hermvals_split_tri);

  // get positions of upper triangle vertices
  tri[0] = (newpos[xs][ys*2]-vars_->shapelet_ctr[0]);
  tri[4] = (newpos[xs+1][(ys+1)*2]-vars_->shapelet_ctr[0]);
  tri[2] = (newpos[xs][(ys+1)*2]-vars_->shapelet_ctr[0]);
  tri[1] = (newpos[xs][ys*2+1]-vars_->shapelet_ctr[1]);
  tri[5] = (newpos[xs+1][(ys+1)*2+1]-vars_->shapelet_ctr[1]);
  tri[3] = (newpos[xs][(ys+1)*2+1]-vars_->shapelet_ctr[1]);
  SHAPELETSOPERA integrate_shapelets (data_, cdata_, vars_, tri, hermvals, hermvals_split_tri);

  for (PS_SIT s=0; s<vars_->numberofshapelets; ++s)
  *sb += hermvals[s] * vars_->mps->get(s);

  // if computing noise in source plane
  if (data_->noisemap && noise)
  {
  VECTOR *dummy1, *dummy2, *product, *hermvec;
  MEMORY ps_malloc (&dummy1, 1);
  MEMORY ps_malloc (&dummy2, 1);
  product = new (dummy1) VECTOR (cdata_, data_, vars_->numberofshapelets);
  hermvec = new (dummy2) VECTOR (cdata_, data_, vars_->numberofshapelets);

  // copy hermvals
  for (PS_SIT s=0; s<vars_->numberofshapelets; ++s)
  hermvec->set (s,hermvals[s]);

  // multiply covariance matrix by shapelet weights
  vars_->a1->inv_dense_mult (hermvec, product, cdata_->numthreads);
  // get variance and standard deviation
  *noise = product->innerproduct (hermvec);
  *noise = std::sqrt (*noise);

  product->~VECTOR();
  hermvec->~VECTOR();
  MEMORY ps_free (dummy1);
  MEMORY ps_free (dummy2);
  }

  MEMORY ps_free (hermvals);
  MEMORY ps_free (hermvals_split_tri);
  }
*/
