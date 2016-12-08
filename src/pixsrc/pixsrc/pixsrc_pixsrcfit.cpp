#ifndef PS_FPT
#ifdef SINGLE_PRECISION
#define PS_FPT float
#endif
#ifdef DOUBLE_PRECISION
#define PS_FPT double
#endif
#endif



#include "pixsrc_pixsrcfit.hpp"
#include "pixsrc_common_analytic.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_statistic.hpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_printer.hpp"
#include <cmath>

#include "pixsrc_common.hpp"

typedef struct
{
    dataandvars *dav;
    int *paramindex;
    int pixsrcnvary;
    int numpixsrcbounds;
    double **pixsrcbounds;
    double *pixsrc0;
} amoebastruct;

void pixsrc_pixsrcfit::pixsrcfit( inputdata* data_, commoninputdata* cdata_, lensvar* vars_ )
{
    if( !data_->pixsrcfit )
        return;

    /*
    PRINTER printerror (data_->print2screenname,
                        "estimating interpolation errors using pixsrcfit is disabled. "
                        "Update GSL minimization in pixsrc_pixsrcfit.cpp to use this function. Alternatively, you can estimate interpolation errors using ASR (analytic source regularization).",
                        cdata_->print2screenmutex);
    */

    if (data_->verbose!=1)
	PRINTER print2screen(data_->print2screenname,
			     "fitting TPS to source for computing interpolation errors",
			     cdata_->print2screenmutex);

    TPS *tpsfitptr, *tpsfit;
    MEMORY ps_malloc (&tpsfitptr, 1);
    tpsfit = new (tpsfitptr) TPS (vars_->triout->pointlist, vars_->triout->pointlist+1,
				  vars_->mps->get_vec_ptr(), 2, 2, 1,
				  vars_->triout->numberofpoints);
    tpsfit->get_tps_weights ();
    tpsfit->cleanup ();
    
    vars_->tpsfitptr = tpsfitptr;
    vars_->tpsfit    = tpsfit;

    if (data_->verbose!=1)
	PRINTER print2screen(data_->print2screenname,
			     "done fitting TPS to source for computing interpolation errors",
			     cdata_->print2screenmutex);

    /*
      if( !vars_->psf_analyticptr )
      {
      MEMORY ps_malloc( &vars_->psf_analyticptr, 1 );
      vars_->psf_analytic = new (vars_->psf_analyticptr) VECTOR( cdata_, data_, vars_->lonc );
      }


      int superrestart = 3;
      int restart = superrestart;
      int pixsrcnvary = 7;
      int numpixsrcbounds = 0;
      double *pixsrcstepsizes = NULL;
      double **pixsrcbounds = (double**)NULL;
      double *pixsrc0;
      double *bestsourceoverall;
      MEMORY ps_malloc( &bestsourceoverall, pixsrcnvary );
      double bestchioverall = 0;

      if( 1 )
      {
      double pos[2], pos2[2];
      double maxflux;
      OPERA assign_n_infinity( &maxflux );
      int maxpos=-1;
      for( int j=0; j<vars_->lonc; ++j )
      if( vars_->mps->get(j) > maxflux )
      {
      maxpos  = j;
      maxflux = vars_->mps->get(j);
      }
      std::copy( &(vars_->triout->pointlist[maxpos*2]),
      &(vars_->triout->pointlist[maxpos*2])+2, pos );

      pthread_mutex_lock( cdata_->wcsmutex );
      HEADER getimgwcscoord( data_, pos[0], pos[1],
      &(pos2[0]), &(pos2[1]) );
      HEADER getwcslfromwcss( pos2[0], pos2[1],
      data_->pra, data_->pdec, data_->r1, data_->r2, pos );
      pthread_mutex_unlock( cdata_->wcsmutex );

      MEMORY ps_malloc( &pixsrc0, 18 );
      pixsrc0[0]  = 1;
      pixsrc0[1]  = 0;
      pixsrc0[2]  = maxflux;
      pixsrc0[3]  = pos[0];
      pixsrc0[4]  = pos[1];
      pixsrc0[5]  = 0.3;
      pixsrc0[6]  = 60;
      pixsrc0[7]  = 0.15;
      pixsrc0[8]  = 0;
      pixsrc0[9]  = 0.5;
      pixsrc0[10] = 1;
      pixsrc0[11] = 1;
      pixsrc0[12] = 1;
      pixsrc0[13] = 1;
      pixsrc0[14] = 1;
      pixsrc0[15] = 1;
      pixsrc0[16] = 0;
      pixsrc0[17] = 1;
      }

      while( restart )
      {
      double *bestsource;
      double **p=0, **pdummy=0, *y=0;

      if( pixsrcnvary )
      {
      double **dither;
      int paramindex[pixsrcnvary];
      MEMORY ps_malloc( &dither, pixsrcnvary,   pixsrcnvary );
      MEMORY ps_malloc( &p     , pixsrcnvary+1, pixsrcnvary );
      MEMORY ps_malloc( &pdummy, pixsrcnvary+2                   );
      MEMORY ps_malloc( &y     , pixsrcnvary+1                   );

      // numerical recipes works on a fortran-like indexing system
      // so pdummy replaces p so that [1..nvary][1..nvary] indexing works
      for( int i=1; i<=pixsrcnvary+1; ++i )
      pdummy[i] = p[i-1]-1;

      int index = 0;
      int tracker;
      double currval;

      for( int src=0; src<pixsrc0[0]; ++src )
      {
      // switching profile type
      switch( (int)pixsrc0[1+src*17] )
      {
      // if profile type = "sersic"
      case 0:
      {
      tracker = 1;
      for( int pp=src*17+1+9; pp<src*17+1+17; ++pp )
      {
      if( pixsrc0[pp] )
      {
      paramindex[index] = tracker;
      p[0][index] = pixsrc0[pp-8];

      if( tracker==1 )
      {
      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      do
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 1.0 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      currval = p[0][index]+dither[nn][index];
      }
      while( currval<0.0 );
      }
      }
      else if( tracker==2 || tracker==3 )
      {
      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 1.0 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      }
      }
      else if( tracker==4 )
      {
      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      do
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 0.25 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      currval = p[0][index]+dither[nn][index];
      }
      while( currval<0.0 || currval>0.99 );
      }
      }
      else if( tracker==5 )
      {
      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 45 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      while( dither[nn][index] >  90 )
      dither[nn][index] -= 180;
      while( dither[nn][index] < -90 )
      dither[nn][index] += 180;
      }
      }
      else if( tracker==6 )
      {
      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      do
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 1.0 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      currval = p[0][index]+dither[nn][index];
      }
      while( currval < CONSTANT smallnumber );
      }
      }
      else if( tracker==8 )
      {
      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      do
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 3.0 * OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      currval = p[0][index]+dither[nn][index];
      }
      while( currval < CONSTANT smallnumber );
      }
      }
      ++index;
      }
      ++tracker;
      }
      break;
      }
      // if profile type = "none"
      case -1:
      {
      tracker = 1;
      for( int pp=src*17+1+9; pp<src*17+1+17; ++pp )
      {
      if( pixsrc0[pp] )
      {
      paramindex[index] = tracker;
      p[0][index] = pixsrc0[pp-8];

      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      if( pixsrcstepsizes )
      dither[nn][index] = pixsrcstepsizes[index] *
      OPERA randomgaussian (cdata_->ps_gsl_ran_r);
      else
      dither[nn][index] = 0;
      }
      ++index;
      }
      ++tracker;
      }
      break;
      }
      default:
      {

      }
      }
      }

      amoebastruct as;
      as.dav             = vars_->dav;
      as.paramindex      = paramindex;
      as.pixsrcnvary     = pixsrcnvary;
      as.numpixsrcbounds = numpixsrcbounds;
      as.pixsrcbounds    = pixsrcbounds;
      as.pixsrc0         = pixsrc0;

      if(data_->verbose!=1)
      {
      PRINTER print2screen(data_->print2screenname,
      "setting up simplex for pixelated source fitting",
      cdata_->print2screenmutex);
      }

      y[0] = PIXSRCFIT setsourceandreturnchi2( &as, pdummy[1] );

      for( int i=1; i<pixsrcnvary+1; ++i )
      {
      std::copy( p[0], p[0]+pixsrcnvary, p[i] );

      for( int nn=0; nn<pixsrcnvary; ++nn )
      {
      p[i][nn] += dither[i-1][nn];
      }

      y[i] = PIXSRCFIT setsourceandreturnchi2( &as, pdummy[i+1] );
      }

      MEMORY ps_free( dither, pixsrcnvary );

      if(data_->verbose!=1)
      {
      PRINTER print2screen(data_->print2screenname,
      "running amoeba for pixelated source fitting",
      cdata_->print2screenmutex);
      }

      int nfunc = 0;
      EXTERNAL ps_pixsrc_amoeba( &as, pdummy, y-1, pixsrcnvary, data_->srcftol,
      PIXSRCFIT setsourceandreturnchi2, &nfunc, 50000 );

      if(data_->verbose!=1)
      {
      PRINTER print2screen(data_->print2screenname,
      "amoeba finished after "+OPERA tostring(nfunc)+
      " evaluations",
      cdata_->print2screenmutex);
      }

      double bestchi = y[0];
      bestsource = p[0];
      for( int i=1; i<=pixsrcnvary; ++i )
      {
      if( y[i]<bestchi )
      {
      bestsource = p[i];
      bestchi    = y[i];
      }
      }

      // if this is the first run through or we did better
      if( restart == superrestart  ||
      (restart != superrestart && bestchi < bestchioverall) )
      {
      std::copy( bestsource, bestsource+pixsrcnvary, bestsourceoverall);
      bestchioverall = bestchi;
      }

      //vars_->residual ->zeromeout();
      //vars_->lensedmps->zeromeout();

      if( data_->verbose!=1 )
      {
      string bs = "";
      for(int j=0; j<pixsrcnvary; ++j)
      {
      bs += OPERA tostring(bestsource[j]);
      //if( j != pixsrcnvary-1 )
      bs += " ";
      }
      bs += ": " + OPERA tostring(bestchi);

      PRINTER print2screen(data_->print2screenname,
      "best source: " + bs,
      cdata_->print2screenmutex);
      }

      MEMORY ps_free( p     , pixsrcnvary+1 );
      MEMORY ps_free( pdummy          );
      MEMORY ps_free( y               );
      }

      // if on last run through while loop
      if( restart == 1 )
      {
      //COMMONANALYTIC setusersource( data_, cdata_, vars_, bestsourceoverall,
      //NULL, pixsrcbounds, numpixsrcbounds, pixsrc0                           );

      // if using analytic source to determine interpolation errors
      if( data_->interperr && data_->pixsrcfit )
      {
      COMMONANALYTIC expandsrc( data_, cdata_, vars_, bestsourceoverall,
      &vars_->bestanalyticsrc, pixsrc0         );
      }

      if( superrestart>1 && data_->verbose!=1 && pixsrcnvary )
      {
      string bs = "";
      for(int j=0; j<pixsrcnvary; ++j)
      {
      bs += OPERA tostring(bestsourceoverall[j]);
      bs += " ";
      }
      bs += ": " + OPERA tostring(bestchioverall);

      PRINTER print2screen(data_->print2screenname,
      "best source overall: " + bs,
      cdata_->print2screenmutex);
      }
      }

      --restart;
      }

      if( 1 )
      {
      MEMORY ps_free( pixsrc0 );
      }

      MEMORY ps_free( bestsourceoverall );

      double sigma = std::sqrt( vars_->variance );
      double perc_diff = 0;
      double sum, diff1, diff2, diff3;
      double fluxlimit = 2.0 * sigma;
      int numdof = 0;
      for( int s=0; s<vars_->lonc; ++s )
      {
      if( std::abs(vars_->psf_analytic->get(s)) < fluxlimit )
      continue;

      sum =  std::max( std::abs( vars_->mps->get(s) ),
      std::abs( vars_->psf_analytic->get(s) ) );

      diff1 = std::abs( vars_->mps->get(s) - vars_->psf_analytic->get(s) );
      diff2 = std::min( diff1, std::abs(diff1 -     sigma) );
      diff3 = std::min( diff2, std::abs(diff1 - 2.0*sigma) );

      perc_diff += diff3 / sum;
      ++numdof;
      }
      perc_diff /= numdof;

      if( perc_diff >= 0.15 )
      PRINTER printwarning( data_->print2screenname,
      "pixelated and analytic source differ by " +
      OPERA tostring( perc_diff*100 ) + "%. Reducing interpolation error.",
      cdata_->print2screenmutex );

      vars_->ie_reduction = ( perc_diff < 0.15 ) ? 1 : 1.0 / (1.0+2.0*perc_diff);
    */
}

double pixsrc_pixsrcfit::setsourceandreturnchi2(void *args, double *source_)
{
    amoebastruct *as = (amoebastruct*)args;
    dataandvars *dav = as->dav;
    inputdata *data_ = dav->data;
    commoninputdata *cdata_ = dav->cdata;
    lensvar *vars_ = dav->vars;
    int *paramindex = as->paramindex;
    int numpixsrcbounds = as->numpixsrcbounds;
    int pixsrcnvary = as->pixsrcnvary;
    double **pixsrcbounds = as->pixsrcbounds;
    double *pixsrc0 = as->pixsrc0;

    // degrees of freedom
    //int dof = vars_->lonc - 7 - 1;

    // this is to account for the fortran-like indexing of amoeba
    double *source = source_ + 1;

    double penalty = 0;
    double emax = 0.999;
    for( int i=0; i<pixsrcnvary; ++i )
    {
        if( paramindex[i]==1 && source[i]<0 )
            penalty += 1.0e8 * ( 1.0 + source[i]*source[i] );
        if( paramindex[i]==4 && source[i]<0 )
            penalty += 1.0e8 * ( 1.0 + source[i]*source[i] );
        if( paramindex[i]==4 && source[i]>emax )
            penalty += 1.0e8 * ( 1.0 + (source[i]-emax)*(source[i]-emax) );
        if( paramindex[i]==6 && source[i]<0 )
            penalty += 1.0e8 * ( 1.0 + source[i]*source[i] );
        if( paramindex[i]==8 && source[i]<0 )
            penalty += 1.0e8 * ( 1.0 + source[i]*source[i] );
    }
    if(penalty)
        return penalty;

    //vars_->residual->zeromeout();
    //vars_->lensedmps->zeromeout();

    COMMONANALYTIC setusersource( data_, cdata_, vars_, source, &penalty,
                                  pixsrcbounds, numpixsrcbounds, pixsrc0, vars_->psf_analytic );

    if( penalty )
        return penalty;

    VECTOR residual( cdata_, data_, vars_->lonc );

    vars_->psf_analytic->minus( vars_->mps, &residual, 1, 1 );

    return residual.innerproduct( &residual );// / ( dof * vars_->variance );
}
