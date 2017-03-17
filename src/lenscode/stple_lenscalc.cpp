#include <unistd.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_monte.h>
#include <gsl/gsl_monte_plain.h>
#include <gsl/gsl_monte_miser.h>
#include <gsl/gsl_monte_vegas.h>
#include <cmath>
#include <iostream>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_integration.h>
#include <iomanip>
//#include <mathlink.h>
#include <cstring>
#include <cfloat>
//#include <gmp.h>
#include <stdexcept>
#include <sstream>
#include <fstream>
#include <iostream>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/file.h>
#include <fcntl.h>
#include <string>
#include <vector>


// contains factorials, cosine integrals, forward declarations, etc
#include "stple_lenscalc.hpp"

// some constants
const double pi   = 3.14159265358979323846264338327950; // known and loved
const double arc2rad = pi/(180.0*3600.0);
const double pinf =  DBL_MAX;
const double ninf = -DBL_MAX;
const int min_inf_sum_iter = 5;   // minimum number of terms in infinite sum
const int max_inf_sum_iter = 149; // maximum number of terms in infinite sum
const int fatalexit        = 0;   // exit (or don't) if fatal error occurs
const int cutoffminutes    = 3;   // timeout (in minutes)

// for keeping track of error printing
static int errprint[en_num_msg];

// pixsrc constants
static void *ps_data=0, *ps_cdata=0;
static size_t ps_numimages=0;
static double *ps_tlmparms_pointer;
static time_t *ps_tlmtime_pointer;
static double *ps_tlmenvirogals_pointer;
static pthread_mutex_t mymutex;
static int numthreads;
static int thread_in_use[100];

int getlock (char const *lockName)
{
    int numloop=0, fd=-1;
    while (-1==fd)
    {
	if (numloop) usleep (100000);
	fd = open( lockName, O_RDONLY);
	if (fd>=0 && flock (fd, LOCK_EX)<0)
	{
	    close (fd);
	    fd = -1;
	}
	++numloop;
    }
    
    return fd;
}

void releaselock (int fd)
{
    if (fd>=0) close (fd);
}

extern "C"
int find_thread_index ()
{
    int fd = getlock ("ids-lock.dat");

    std::vector <std::string> originaldata;
    {
	std::string line;
	std::ifstream infile;
	infile.open ("ids.dat");
        while (!infile.eof())
        {
	    getline (infile,line);
	    if (!line.empty())
		originaldata.push_back (line);
        }
	infile.close();
    }

    int numloop = 0;
    int index = -1;
    while (-1==index)
    {
	if (numloop) usleep (100000);

	int lineno = 0;
	std::string line;
	std::ifstream infile;
	infile.open ("ids.dat");
        while (!infile.eof())
        {
	    getline (infile,line);
	    if (line[0]=='0')
	    {
		index = lineno;
		originaldata[index] = std::string("1");
		break;
	    }
	    ++lineno;
        }
	infile.close();
	++numloop;
    }
    
    std::ofstream file;
    file.open("ids.dat");
    for (int com=0; com<(int)originaldata.size(); ++com)
	file << originaldata[com] << std::endl;
    file.close();    
    
    releaselock (fd);
    return index;
}

extern "C"
void release_thread_index (int ind)
{
    int fd = getlock ("ids-lock.dat");
    
    std::vector <std::string> originaldata;
    {
	std::string line;
	std::ifstream infile;
	infile.open ("ids.dat");
        while (!infile.eof())
        {
	    getline (infile,line);
	    if (!line.empty())
		originaldata.push_back (line);
        }
	infile.close();
    }
    originaldata[ind] = std::string("0");

    std::ofstream file;
    file.open("ids.dat");
    for (int com=0; com<(int)originaldata.size(); ++com)
	file << originaldata[com] << std::endl;
    file.close();    
    
    releaselock (fd);
    return;
}
  
// struct for GSL minimization
struct gsl_mdm
{
    gsl_mdm() : T(gsl_multimin_fminimizer_nmsimplex2rand),
                s(NULL), ss(NULL), iter(0) {}

    const gsl_multimin_fminimizer_type *T;
    gsl_multimin_fminimizer *s;
    gsl_vector *ss;
    gsl_multimin_function minex_func;
    size_t iter;
    int status;
    double size;
};

// struct for chi2 minimization
struct minchist
{
    void *p;
    void *eg;
    time_t t;
};

extern "C" 
void potdefmag    ( double x, double y,  double *ig1,
		    double *defx, double *defy, double *ig2,
		    double *ig3, double *ig4, double *ig5, double **parms__, time_t **tlmtime__, double **envirogals_, int imagenumber)
{
//    std::cout << "myparms " << parms__ << std::endl;
//    usleep(10000000000000);
    double *parms = *parms__;
    double *envirogals = *envirogals_;
    time_t starttime = **tlmtime__;
  // check for timeout
//    std::cout << "pdm " << x << " " << y << " " << envirogals_ << std::endl;
//    std::cout << "pdm " << x << " " << y << " " << envirogals << std::endl;
//    std::cout << "pdm " << x << " " << y << " " << envirogals[0] << std::endl;
  if (time(NULL)-starttime>cutoffminutes*60)
  {
//      std::cout << time(NULL) << " " << starttime << std::endl;
      *defx = 0*rand()/(double)RAND_MAX-0.5;
      *defy = 0*rand()/(double)RAND_MAX-0.5;
      return;
  }  

  // setup parameters for deflection calculations
  int pind = 0;
  pind +=9; pind +=3; pind +=1; pind +=1; pind +=1; pind +=1; pind +=1; pind +=1; pind +=1; pind +=1; pind +=1;
  double newparms[pind+1+1+1+6];
  std::copy (parms,parms+pind,newparms);
  newparms[pind] = 1;
  pind +=1;
  newparms[pind] = -x / parms[16]; // convert to Nbody units (negative sign accounts for left-increasing R.A.)
  pind +=1;
  newparms[pind] =  y / parms[16]; // convert to Nbody units
  pind +=1;
  std::fill (newparms+pind,newparms+pind+6,0);
  pind +=6;
  
  // get deflections and rescale
  stple_lenscalc (newparms,envirogals,starttime,defx,defy, imagenumber);
  *defx *= -parms[16];
  *defy *= parms[16];
}
extern "C"
void find_img     ( double u, double v,  int *numimg, double ***imgs)
{
}

// optimize image-plane chi2 (callable from python)
extern "C"
void stple_imageplane_chi2 (void* params, void *envirogals_)
{
    // start timing!
    time_t starttime;
    starttime = time(NULL);
    
    // handle what the code prints and how GSL terminates
    gsl_set_error_handler_off();
//    std::cout.setstate(std::ios_base::failbit);
//    std::cerr.setstate(std::ios_base::failbit);

    double *parms = (double*)params;
    double *envirogals = (double*)envirogals_;

    // parse incoming array of parameters and results array
    int pind = 0;
    double* rotmat = parms+pind;
    pind +=9;
    double* aratio = parms+pind;
    pind +=3;
    double* s = parms+pind;
    pind +=1;
    double* sigmacr = parms+pind;
    pind +=1;
    double* MRe = parms+pind;
    pind +=1;
    double* rscale = parms+pind;
    pind +=1;
    double* rscale_arc = parms+pind;
    pind +=1;
    double* rcore = parms+pind;
    pind +=1;
    double* cutrad = parms+pind;
    pind +=1;
    double* user_PA = parms+pind;
    pind +=1;
    double* t0 = parms+pind;
    pind +=1;
    double* numimg_d = parms+pind;
    pind +=1;
    double* coor = parms+pind;
    pind +=2*(int)(*numimg_d);
    double* tdelobs = parms+pind;
    pind +=1*(int)(*numimg_d);
    double* coor_sigma = parms+pind;
    pind +=2*(int)(*numimg_d);
    double* tdelobs_sigma = parms+pind;
    pind +=1*(int)(*numimg_d);
    double* resultarr = parms+pind;
    pind +=4;

    int numimg = (int)(*numimg_d);
    
    int numeg = (int)envirogals[0];
    int lens2startind = 2+3*numeg;
    int havelens2 = (int)envirogals[lens2startind];
    double *l2parms = envirogals + lens2startind+1;
    double origl2p = l2parms[-1];
    l2parms[-1] = 0;
    

    // setup array to pass to lensing calculations 
    // structure of parms4calc array:
    //    rotmat, aratio, s, sigmacr, MRe, rscale, rcore, cutrad, numimg, coor, resultsarr(3+3*numimg elements),
    //    coor, tdel, position error, tdel error, source bounds(4 elements)
    // setup parms4calc array and some pointers (for convenience)
    int parms4calclength = pind + 3+5*numimg;
    double  parms4calc[parms4calclength];
    //                         start_of_pointer      array_length      coor2,tdel,coorerr,tdelerr,src_bounds      results        coor
    double *parms4calccoor   = parms4calc       +  parms4calclength - (2*numimg+1*numimg+2*numimg+1*numimg+4) - (3+3*numimg) - (numimg*2);
    double *parms4calcresult = parms4calccoor   + (numimg*2);
    double *parms4calccoor2  = parms4calcresult + (3+3*numimg);
    double *parms4calcbounds = parms4calccoor2  + (2*numimg+1*numimg+2*numimg+1*numimg);
    // copy values into new array
    std::copy (parms, parms + (parms4calcresult-parms4calc), parms4calc);
    std::copy (parms4calccoor, parms4calccoor+numimg*2, parms4calccoor2);
    std::copy (parms + (parms4calcresult-parms4calc), 
	       parms + (parms4calcresult-parms4calc) + (parms4calcbounds-(parms4calccoor2+numimg*2)),
	       parms4calccoor2 + numimg*2);
    
    resultarr[3] = 1.e100;

    // do one run to get initial guesses
    stple_lenscalc (parms4calc, envirogals, starttime, 0,0, -1);

    // store ellipticity and mass normalization
    std::copy (parms4calcresult, parms4calcresult+3, resultarr);
    
    // get ray-traced image positions to get source position bounding box and best guess
    double srcminmaxavg[6]={pinf,ninf,0,pinf,ninf,0}, raytracex, raytracey;
    for (int j=0; j<numimg; ++j)
    {
	raytracex = coor[j*2+0] - parms4calcresult[3+j*2+0];
	raytracey = coor[j*2+1] - parms4calcresult[3+j*2+1];

	srcminmaxavg[2] += raytracex;
	srcminmaxavg[5] += raytracey;
	if (raytracex<srcminmaxavg[0]) srcminmaxavg[0] = raytracex;
	if (raytracex>srcminmaxavg[1]) srcminmaxavg[1] = raytracex;
	if (raytracey<srcminmaxavg[3]) srcminmaxavg[3] = raytracey;
	if (raytracey>srcminmaxavg[4]) srcminmaxavg[4] = raytracey;
    }
    srcminmaxavg[2] /= numimg;
    srcminmaxavg[5] /= numimg;
    // store source bounds
    parms4calcbounds[0] = srcminmaxavg[0];
    parms4calcbounds[1] = srcminmaxavg[1];
    parms4calcbounds[2] = srcminmaxavg[3];
    parms4calcbounds[3] = srcminmaxavg[4];
    
    // optimize source position (and image positions) and get image-plane chi2
    // setup GSL minimization struct
    struct gsl_mdm mdm;
    int ndim = 2 + numimg*2;
    double tolerance = 1e-3*100;
    int maxiter = 1000;
    struct minchist mcs;
    mcs.p = parms4calc;
    mcs.eg = envirogals;
    mcs.t = starttime;
    mdm.minex_func.n      = ndim;
    mdm.minex_func.f      = get_image_plane_chi2;
    mdm.minex_func.params = &mcs;
    mdm.s                 = gsl_multimin_fminimizer_alloc (mdm.T, mdm.minex_func.n);
    // setup initial guesses and stepsizes
    gsl_vector *x  = gsl_vector_alloc (ndim);
    gsl_vector *ss = gsl_vector_alloc (ndim);
    gsl_vector_set (x,  0, srcminmaxavg[2]);
    gsl_vector_set (x,  1, srcminmaxavg[5]);
    gsl_vector_set (ss, 0, 0.1*(srcminmaxavg[1]-srcminmaxavg[0]));
    gsl_vector_set (ss, 1, 0.1*(srcminmaxavg[4]-srcminmaxavg[3]));
    for (int j=0; j<numimg; ++j)
    {
	gsl_vector_set (x,  2+j*2+0, coor[j*2+0]);
	gsl_vector_set (x,  2+j*2+1, coor[j*2+1]);
	gsl_vector_set (ss, 2+j*2+0, coor_sigma[j*2+0]*0.2);
	gsl_vector_set (ss, 2+j*2+1, coor_sigma[j*2+1]*0.2);
    }
    
    // perform iterative minimization
    gsl_multimin_fminimizer_set (mdm.s, &mdm.minex_func, x, ss);
    do
    {
        ++mdm.iter;
        mdm.status = gsl_multimin_fminimizer_iterate (mdm.s);
	
	if (mdm.status) break;

        mdm.size   = gsl_multimin_fminimizer_size (mdm.s);
        mdm.status = gsl_multimin_test_size (mdm.size, tolerance);
    }
    while (mdm.status==GSL_CONTINUE && mdm.iter<(unsigned int)maxiter);
    resultarr[3] = gsl_multimin_fminimizer_minimum (mdm.s);
    gsl_vector *best = gsl_multimin_fminimizer_x (mdm.s);
    
    double srcposfound[2] = {gsl_vector_get(mdm.s->x,0),gsl_vector_get(mdm.s->x,1)};

    gsl_multimin_fminimizer_free (mdm.s);
    gsl_vector_free (x);
    gsl_vector_free (ss);
    
    l2parms[-1] = origl2p;
    
    // get chi^2 from extended source reconstruction
    if (ps_data && ps_cdata && ps_numimages)
    {
	l2parms[1] = srcposfound[0];
	l2parms[2] = srcposfound[1];

	// copy lens model parameters
	ps_tlmparms_pointer = parms;
	ps_tlmenvirogals_pointer = envirogals;
	ps_tlmtime_pointer  = &starttime;
	// do source reconstruction
	double lE;
	ps_launch(ps_data, ps_cdata, &lE);
	resultarr[3] += -2.0*(lE-0);
    }
}

struct st4dtf
{
  void *p;
  time_t t;
};

// lensing calculations function
void stple_lenscalc (void *params, double *envirogals, time_t starttime, double *defx, double *defy, int imagenumber)
{
    double *parms = (double*)params;

    // parse incoming array of parameters and results array
    int pind = 0;
    double* rotmat = parms+pind;
    pind +=9;
    double* aratio = parms+pind;
    pind +=3;
    double* s = parms+pind;
    pind +=1;
    double* sigmacr = parms+pind;
    pind +=1;
    double* MRe = parms+pind;
    pind +=1;
    double* rscale = parms+pind;
    pind +=1;
    double* rscale_arc = parms+pind;
    pind +=1;
    double* rcore = parms+pind;
    pind +=1;
    double* cutrad = parms+pind;
    pind +=1;
    double* user_PA = parms+pind;
    pind +=1;
    double* t0 = parms+pind;
    pind +=1;
    double* numimg_d = parms+pind;
    pind +=1;
    double* coor = parms+pind;
    pind +=2*(int)(*numimg_d);
    double* resultarr = parms+pind;
    pind +=3+3*(int)(*numimg_d);

    double rc=*rcore, rt=*cutrad;
    int numimg = (int)(*numimg_d);

    double r[numimg],phi[numimg];

    // convert to polar coordinates
    for (int j=0; j<numimg; ++j)
    {
	r[j]   = std::sqrt  (coor[2*j]*coor[2*j]+coor[2*j+1]*coor[2*j+1]);
	phi[j] = std::atan2 (coor[2*j+1],coor[2*j]);
    }
    
    // get axes ratios and zero out deflections and time delays
    double q,p;
    q   = aratio[1];
    p   = aratio[2];
    double *deflection = resultarr+3;
    double *tdel = resultarr+3+numimg*2;
    std::fill (tdel,tdel+numimg,0);
    std::fill (deflection,deflection+numimg*2,0);

    // this will store parameters for calculating mass normalizations
    double massnormparms[4*2+1];
    massnormparms[8] = (*s-3)*0.5;

    // loop over density terms (inner core radiud, outer truncation radius)
    for (int rads=0; rads<2; ++rads)
    {
    double tdel4ref=0;
    // loop over number of images (or positions where lensing calculations are done)
    for (int j=0; j<numimg; ++j)
    {
	// scale radii (squared) along different axes
	double rm2,q2,p2, rcrt0,rcrt = 0==rads ? rc : rt;
	q2  = q*q*rcrt*rcrt;
	p2  = p*p*rcrt*rcrt;
	rm2 =     rcrt*rcrt;
	rcrt0 = rcrt;
	rcrt = 1;

	// Euler angle related stuff + mu
	double B1f,B2f,Af,C1f,C2f,C3f,mu,cx,cy,cxy,*rm=rotmat;
	B1f = rm[0]*rm[2]/rm2 + rm[3]*rm[5]/q2 + rm[6]*rm[8]/p2;
	B2f = rm[1]*rm[2]/rm2 + rm[4]*rm[5]/q2 + rm[7]*rm[8]/p2;
	Af  = rm[2]*rm[2]/rm2 + rm[5]*rm[5]/q2 + rm[8]*rm[8]/p2;
	C1f = rm[0]*rm[0]/rm2 + rm[3]*rm[3]/q2 + rm[6]*rm[6]/p2;
	C2f = rm[1]*rm[1]/rm2 + rm[4]*rm[4]/q2 + rm[7]*rm[7]/p2;
	C3f = rm[0]*rm[1]/rm2 + rm[3]*rm[4]/q2 + rm[6]*rm[7]/p2;
	mu  = (*s-3)*0.5;
	// my viewing angle (x-axis) is different from Chae et al. (1998) viewing angle (z-axis),
	// so some formulae are different
	// (see "norm" variable below as well)
	/*
	  cx  = C1f-B1f*B1f/Af;
	  cy  = C2f-B2f*B2f/Af;
	  cxy = 2*(C3f-B1f*B2f/Af);
	*/
	cx  = C2f-C3f*C3f/C1f;
	cy  = Af -B1f*B1f/C1f;
	cxy = 2*(B2f-C3f*B1f/C1f);

	// main variables used in calculations (related to core radius, ellipticity, PA)
	double P,S,Q;
	Q = 0.5 * sgn_pos_neg(cxy) * std::sqrt(cxy*cxy+(cx-cy)*(cx-cy));
	P = 0.5*(cx+cy);
	//S = cx==cy && 0==cxy ? 0 : std::atan2(cx-cy,cxy);
	if (Q/P<0)
	    S = -2*(*user_PA +0.5*pi -0.5*pi -0.25*pi);
	else
	    S = -2*(*user_PA -0.5*pi -0.25*pi);
	
	
	// calculate ellipticity and PA of this model (really, e*cos(2phi), e*sin(2phi))
	double tmp0,ell,pa,e1,e2;
	tmp0 = Q/P;
	pa   = 0.75*pi-0.5*S;
	if (tmp0<0) pa += 0.5*pi;
	tmp0 = std::abs(tmp0);
	ell = 1-std::sqrt((1-tmp0)/(1+tmp0));
	e1 = ell*std::cos(2*pa);
	e2 = ell*std::sin(2*pa);
	if (0==rads)
	{
	    resultarr[0] = e1;
	    resultarr[1] = e2;
	}

	// normalization for this density term
	double norm = (0==rads ? 1 : -1) *std::pow(2.,2*mu+1) *gsl_sf_beta(mu+1,mu+1)
	    *std::pow(rcrt0,-*s) /std::sqrt(/*Af*/C1f);
	// characteristic radius used in determining fastest converging formula to use for series expansion
	double reg = 1/std::sqrt(std::max(1e-13,P-std::abs(Q)));
	
	double hh,eps1,eps2,I0,del_int_I0=0,def_tdel[3] = {0,0,0};
	double tmp1,tmp2,tmp3,tmp4,tmp10,tmp11,tmp12,tmp13,tmp14;
	
	// some terms to calculate eps1,eps2 (series expansion variables) + hh
	tmp1  = rcrt*rcrt+P*r[j]*r[j];  //  rc+Prr
	tmp2  = tmp1*tmp1;              // (rc+Prr)^2
	tmp3  = Q*r[j]*r[j];            //  Qrr
	tmp4  = tmp3*tmp3;              // (Qrr)^2
	tmp10 = tmp3/tmp1;              //  Qrr/(rc+Prr)
	tmp11 = tmp10*tmp10;            // (Qrr/(rc+Prr))^2
	tmp12 = std::sqrt(1-tmp11);     // sqrt( 1 - (Qrr/(rc+Prr))^2  ) 
	tmp13 = tmp2-tmp4;              // (rc+Prr)^2 - (Qrr)^2
	tmp14 = std::sqrt(tmp13);       // sqrt( (rc+Prr)^2 - (Qrr)^2  )
	hh    = r[j] / tmp14;
	eps1  = (1-tmp12) / (1+tmp12);
	eps2  = 0.5* ( tmp1/tmp13 + 1/tmp14 );

	// extra zeros at end of parms2 array are for storing temporary variables in nested sums
	// get I0
	{
	    I0 = 0;
	    double parms2[] = {hh,eps1,eps2,P,Q,S,mu,reg,r[j],phi[j],0,0,0,0,0};
	    infinite_sum_func (5, max_inf_sum_iter, 1, 1e-5, I0_func, parms2, &I0, "I0_2m", starttime);
	    if (r[j]<=reg)
		I0 *= hh*std::pow(eps2,mu);
	    else
		I0 = (1/(r[j]*std::sqrt(P*P-Q*Q))-hh*std::pow(eps2,mu)*I0)/mu;
	}

	// get delta int_I0
	if ( ! (defx && defy) )
	{
	  double parms2[] = {P,Q,S,mu,reg,rcrt,r[0],r[j],0,0,0,0,0};
	    integrate_I0(parms2, starttime);
	    del_int_I0 = parms2[8];
	}

	// get I1,I2 / deflections,time delays
	{
	    double parms2[] = {hh,eps1,eps2,P,Q,S,mu,reg,r[j],phi[j],0,0,0,0,0};
	    struct st4dtf dtf;
	    dtf.p = parms2;
	    dtf.t = starttime;
	    infinite_sum_func (5, max_inf_sum_iter, 3, 1e-5, def_tdel_func, &dtf, def_tdel, "I1_2m/I2_2m", starttime);
	    def_tdel[0]  = norm * (I0+def_tdel[0]);
	    def_tdel[1] *= norm;
	    def_tdel[2]  = norm * (del_int_I0-0.5*r[j]*def_tdel[2]);
	}

	// add time delays from this density term
	// time of arrival of 1st image is used as reference for time delays
	if (0==j) 
	    tdel4ref = def_tdel[2];
	else
	    tdel[j] += def_tdel[2]-tdel4ref;

	// convert deflections back to Cartesian coordinates
	double cosphi,sinphi,def_x,def_y;
	cosphi = std::cos(phi[j]);
	sinphi = std::sin(phi[j]);
	def_x  = def_tdel[0]*cosphi - def_tdel[1]*sinphi;
	def_y  = def_tdel[0]*sinphi + def_tdel[1]*cosphi;

	// add deflections from this density term
	deflection[j*2+0] += def_x;
	deflection[j*2+1] += def_y;

	// store mass normalization for this density term
	massnormparms[rads*4+0] = Q;
	massnormparms[rads*4+1] = P;
	massnormparms[rads*4+2] = S;
	massnormparms[rads*4+3] = norm;
    }
    }

    // calculate mass normalization using brute force
    // mass_inf, mass_Re are in N-body nuits
    // Minf, MRe are in physical units
    double Re, Minf, mass_inf=0, mass_Re=0;
    Re = std::sqrt(*MRe/(pi**sigmacr)) / *rscale;
    {
	gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	double result, error;
	gsl_function F;
	F.function = massnorm_func;
	F.params = massnormparms;

	// 0 to cutoff radius
	gsl_integration_qags (&F, 0, rt, 0, 1e-5, 10000, w, &result, &error);
	mass_inf += 2*pi*result;
	// cutoff radius to infinity
	gsl_integration_qagiu (&F, rt, 0, 1e-5, 10000, w, &result, &error);
	mass_inf += 2*pi*result;
	// 0 to Einstein radius
	gsl_integration_qags (&F, 0, Re, 0, 1e-5, 10000, w, &result, &error);
	mass_Re += 2*pi*result;
	
	gsl_integration_workspace_free (w);
    }
    Minf = *MRe * mass_inf/mass_Re;
    resultarr[2] = mass_Re/mass_inf;
    
    // apply mass normalization and sigma_crit scaling (dimensionless units)
    double massnorm = mass_inf * *sigmacr**rscale**rscale/Minf;
    for (int j=0; j<numimg; ++j)
    {
	deflection[j*2+0] /= massnorm;
	deflection[j*2+1] /= massnorm;
	tdel[j]           /= massnorm;
    }
    
    // add deflections and delays from galaxies nearby (SIS)
    int numeg = (int)envirogals[0];
    int lens2startind = 2+3*numeg;
    int havelens2 = (int)envirogals[lens2startind];
    if (numeg>0)
    {
	double maxrE = envirogals[1];
	double newx0,newy0,newr0,newx,newy,newr,newc,news,egdefx=0,egdefy=0;
	for (int g=0; g<numeg; ++g)
	{
	    newx0  = coor[0]-envirogals[2+g*3+0]/ *rscale_arc;
	    newy0  = coor[1]-envirogals[2+g*3+1]/ *rscale_arc;
	    newr0  = std::sqrt(newx0*newx0+newy0*newy0);
	    for (int j=0; j<numimg; ++j)
	    {
		newx  = coor[j*2+0]-envirogals[2+g*3+0]/ *rscale_arc;
		newy  = coor[j*2+1]-envirogals[2+g*3+1]/ *rscale_arc;
		newr  = std::sqrt(newx*newx+newy*newy);
		newc  = newx/newr;
		news  = newy/newr;
		deflection[j*2+0] += newc *envirogals[2+g*3+2]*maxrE/ *rscale_arc;
		deflection[j*2+1] += news *envirogals[2+g*3+2]*maxrE/ *rscale_arc;
		tdel[j]           -= newr0*envirogals[2+g*3+2]*maxrE/ *rscale_arc;
		tdel[j]           += newr *envirogals[2+g*3+2]*maxrE/ *rscale_arc;
	    }
	}
    }

    if (havelens2 && 1==imagenumber)
    {
	double newx0,newy0,newr0,newx,newy,newr,newc,news;
	double *l2parms = envirogals + lens2startind+1;
	newx0  = coor[0]-l2parms[1];
	newy0  = coor[1]-l2parms[2];
	newr0  = std::sqrt(newx0*newx0+newy0*newy0);
	for (int j=0; j<numimg; ++j)
	{
	    newx  = coor[j*2+0]-0*l2parms[1];
	    newy  = coor[j*2+1]-0*l2parms[2];
	    newr  = std::sqrt(newx*newx+newy*newy);
	    newc  = newx/newr;
	    news  = newy/newr;
	    deflection[j*2+0] += newc *l2parms[0]/ *rscale_arc;
	    deflection[j*2+1] += news *l2parms[0]/ *rscale_arc;
	    deflection[j*2+0] *= l2parms[3];
	    deflection[j*2+1] *= l2parms[3];
	}
    }

    // save deflections if calculating for extended source reconstruction
    if (defx && defy)
    {
	*defx = deflection[0];
	*defy = deflection[1];
    }
}






//////////////////////////////////////////////////////////////////////
// Functions below are lower-level routines or lensing calculations //
//////////////////////////////////////////////////////////////////////




double get_image_plane_chi2 (const gsl_vector *v, void* params__)
{
  struct minchist *mcs = (struct minchist*)params__;
  void *params = mcs->p;
  double *envirogals = (double*)(mcs->eg);
  time_t starttime;
  starttime = mcs->t;
  
    // This function is passed to GSL minimizer
    
    double *parms = (double*)params;

    // parse incoming array of parameters and results array
    int pind = 0;
    double* rotmat = parms+pind;
    pind +=9;
    double* aratio = parms+pind;
    pind +=3;
    double* s = parms+pind;
    pind +=1;
    double* sigmacr = parms+pind;
    pind +=1;
    double* MRe = parms+pind;
    pind +=1;
    double* rscale = parms+pind;
    pind +=1;
    double* rscale_arc = parms+pind;
    pind +=1;
    double* rcore = parms+pind;
    pind +=1;
    double* cutrad = parms+pind;
    pind +=1;
    double* user_PA = parms+pind;
    pind +=1;
    double* t0 = parms+pind;
    pind +=1;
    double* numimg_d = parms+pind;
    pind +=1;
    double* coor = parms+pind;
    pind +=2*(int)(*numimg_d);
    double* resultarr = parms+pind;
    pind +=3+3*(int)(*numimg_d);
    double* coor0 = parms+pind;
    pind +=2*(int)(*numimg_d);
    double* tdelobs = parms+pind;
    pind +=1*(int)(*numimg_d);
    double* coor_sigma = parms+pind;
    pind +=2*(int)(*numimg_d);
    double* tdelobs_sigma = parms+pind;
    pind +=1*(int)(*numimg_d);
    double* srcbounds = parms+pind;
    pind +=4;

    int numimg = (int)(*numimg_d);
    
    // copy source position and image positions
    double srcpos[2] = {gsl_vector_get(v,0),gsl_vector_get(v,1)};
    for (int j=0; j<numimg; ++j)
    {
	coor[j*2+0] = gsl_vector_get(v,2+j*2+0);
	coor[j*2+1] = gsl_vector_get(v,2+j*2+1);
    }
    
    double delsrcposx,delsrcposy, tmp1, tmp2, tmp3, chi2 = 0;
    double chival1=0,chival2=0,chival3=0,chival4=0;
    for (int j=0; j<numimg; ++j)
    {
	// get position chi2
	tmp1  = coor[j*2+0]-coor0[j*2+0];
	if (std::abs(tmp1)>3*coor_sigma[j*2+0])
	    chival2 += (1+tmp1*tmp1)*1e10;
	else
	    chival2 += tmp1*tmp1 / (coor_sigma[j*2+0]*coor_sigma[j*2+0]);
	tmp1  = coor[j*2+1]-coor0[j*2+1];
	if (std::abs(tmp1)>3*coor_sigma[j*2+0])
	    chival3 += (1+tmp1*tmp1)*1e10;
	else
	    chival3 += tmp1*tmp1 / (coor_sigma[j*2+1]*coor_sigma[j*2+1]);
    }
    
    if (chival2+chival3>1e10)
	return chival2+chival3;

    // get deflections and potential
    stple_lenscalc (params, envirogals, starttime, 0,0, -1);
    
    // get model predictions and chi2
    for (int j=0; j<numimg; ++j)
    {
	// penalize severely if lens equation violated
	delsrcposx = coor[j*2+0] - resultarr[3+j*2+0] - srcpos[0];
	delsrcposy = coor[j*2+1] - resultarr[3+j*2+1] - srcpos[1];
	if (0&&delsrcposx*delsrcposx+delsrcposy*delsrcposy > 5e-3*5e-3)
	    chival1 += (1+delsrcposx*delsrcposx+delsrcposy*delsrcposy)*1e6;
	else chival1 += (delsrcposx*delsrcposx+delsrcposy*delsrcposy) / (coor_sigma[0]*coor_sigma[0]/10./10.);
	
	// get time delay chi2
	tmp1  = (coor[0*2+0]-srcpos[0])*(coor[0*2+0]-srcpos[0])
	    +   (coor[0*2+1]-srcpos[1])*(coor[0*2+1]-srcpos[1]);
	tmp2  = (coor[j*2+0]-srcpos[0])*(coor[j*2+0]-srcpos[0])
	    +   (coor[j*2+1]-srcpos[1])*(coor[j*2+1]-srcpos[1]);
	tmp3  = 0.5*(tmp2-tmp1) - resultarr[3+numimg*2+j];
	chival4 += (tmp3-tdelobs[j])*(tmp3-tdelobs[j]) / (tdelobs_sigma[j]*tdelobs_sigma[j]);
    }
    
    chi2 = chival1+chival2+chival3+chival4;
    
    return chi2;
}


extern "C"
void setup (char *basename, int numthreads_)
{
    // we haven't printed any error messages yet
    std::fill (errprint,errprint+en_num_msg,0);

    // setup pixsrc
    if (basename)
    {
	numthreads = numthreads_;
//	std::cout << "mutex init " << &mymutex << " " << (int)pthread_mutex_init(&mymutex,NULL) << std::endl;
	std::fill (thread_in_use,thread_in_use+100,0);
	size_t memorysizedata,memorysizecdata;
	ps_getpixsrcinfo(basename, &memorysizedata, &memorysizecdata, &ps_numimages);
	ps_data  = (void*)malloc(memorysizedata );
	ps_cdata = (void*)malloc(memorysizecdata);
	ps_loaddata(basename, ps_numimages, ps_data, ps_cdata, &ps_tlmparms_pointer, &ps_tlmtime_pointer, &ps_tlmenvirogals_pointer);
//	std::cout << &ps_tlmenvirogals_pointer << std::endl;
//	std::cout << ps_tlmenvirogals_pointer << std::endl;
    }
}







void I0_func (int m, double *sum, void *parms_)
{
    double *parms,hh,eps1,eps2,mu,P,Q,S,reg,r,phi,thissum=0;
    
    parms = (double*)parms_;
    hh   = parms[0];
    eps1 = parms[1];
    eps2 = parms[2];
    P    = parms[3];
    Q    = parms[4];
    S    = parms[5];
    mu   = parms[6];
    reg  = parms[7];
    r    = parms[8];
    phi  = parms[9];

    int k = m-1;

    if (r<=reg)
    {
	for (int el=0; el<=k; ++el)
	    thissum += (0==el%2 ? 1 : -1) * 
		safe_mul_div (zeta1_func(k,el,0,mu),std::pow(eps2,el),
			      safe_2F1 (-el-mu,-el-mu,1,eps1, 0), 
			      1,1,1);
    }
    else
    {
	thissum += std::pow(eps2,k) * safe_2F1(-k-mu,-k-mu,1,eps1, 1);
    }
    
    *sum += thissum;
}

void def_tdel_func (int m, double *sum, void *parms__)
{
    struct st4dtf *dtf = (struct st4dtf*)parms__;
    void *parms_ = dtf->p;
    time_t starttime = dtf->t;
  
    double *parms,hh,eps1,eps2,mu,P,Q,S,reg,r,phi;
    
    parms = (double*)parms_;
    hh   = parms[0];
    eps1 = parms[1];
    eps2 = parms[2];
    P    = parms[3];
    Q    = parms[4];
    S    = parms[5];
    mu   = parms[6];
    reg  = parms[7];
    r    = parms[8];
    phi  = parms[9];

    double I1[2],I2,prefac,ang,mcos,msin;
    prefac = sgnpow(-Q,m) * hh * std::pow(eps1,0.5*m);
    parms[10] = m;

    {
	I1[0] = I1[1] = 0;
	infinite_sum_func (5, max_inf_sum_iter, 2, 1e-5, I1_func, parms_, I1, "I1_2m_k", starttime);	
	double rb = 1/std::sqrt(P);
	double tmp0 = r/rb;
	double tmp1 = tmp0*tmp0;
        double t    = tmp1/(1+tmp1);

	if (tmp0>=1./2.&&tmp0<=2)
	    I1[0] *= 2.*r/pi*std::pow(-t*std::abs(Q/P),m)*std::pow(1-t,mu+1);
	else if (r<=reg)
	    I1[0] *= prefac * std::pow(eps2,mu);
	else
	    I1[0] = -prefac 
		*safe_mul_div (gsl_sf_gamma(mu-m),gsl_sf_gamma(mu+1),gsl_sf_gamma(m+1),
			       1,-1,-1)
		*(std::pow(eps2,mu)*I1[0] - std::pow(eps2,m)*I1[1]);
    }

    {
	I2 = 0;
	infinite_sum_func (5, max_inf_sum_iter, 1, 1e-5, I2_func, parms_, &I2, "I2_2m_k", starttime);
	double rb = 1/std::sqrt(P);
	double tmp0 = r/rb;

	if (tmp0>=1./2.&&tmp0<=2)
	{
	    double tmp1 = tmp0*tmp0;
	    double t    = tmp1/(1+tmp1);
	    I2 *= 2.*r/pi*std::pow(-t*std::abs(Q/P),m)*std::pow(1-t,mu+1)*std::pow(mu+m,-t);
	}
	else if (1==m && r<1./std::sqrt(P+std::abs(Q)))
	{
	    double epsQP = std::abs(Q/P);
	    if (epsQP<1e-3)
		I2  = r*(-0.5*epsQP-2*I2/pi);
	    else
		I2  = r*(-P/Q*(1-std::sqrt(1-epsQP*epsQP))-2*I2/pi);
	}
	else if (eps1<0.1 || eps2>0.9)
	{
	    // numerical integration (phi integral done analytically) done instead of using analytic result
	    // Note: find fast converging analytic result.

	    // r integral
	    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
	    double result, error;
	    gsl_function F;
	    F.function = I2_for_int_func_r;
	    F.params = parms;
	    gsl_integration_qagiu (&F, r, 0, 1e-5, 10000, w, &result, &error);
	    gsl_integration_workspace_free (w);
	    I2 = safe_mul_div (result,2*sgnpow(Q,m),pochhammer_func(-mu,m),std::pow(r,1-2*m), 
			       1,1,-1,-1);
	}
	else
	{
	    I2 *= safe_mul_div (prefac,std::pow(eps2,mu),zeta2_func(m,mu), 
				1,1,1);	
	}
    }	
    ang     = m*(2*phi-(0.5*pi-S));
    mcos    = std::cos(ang);
    msin    = std::sin(ang);
    sum[0] += mcos * (I1[0]-I2);
    sum[1] += msin * (I1[0]+I2);
    sum[2] += mcos * (I1[0]+I2) /m;
}

void I1_func (int m, double *sum, void *parms_)
{
    double *parms,hh,eps1,eps2,mu,P,Q,S,reg,r,phi,twom,thissum[2]={0,0};
    
    parms = (double*)parms_;
    hh   = parms[0];
    eps1 = parms[1];
    eps2 = parms[2];
    P    = parms[3];
    Q    = parms[4];
    S    = parms[5];
    mu   = parms[6];
    reg  = parms[7];
    r    = parms[8];
    phi  = parms[9];
    twom = parms[10];

    int k=m-1;
    m = (int)twom;
    double rb = 1/std::sqrt(P);
    double tmp0 = r/rb;

    if (tmp0>=1./2.&&tmp0<=2)
    {
	double tmp1 = tmp0*tmp0;
	double t    = tmp1/(1+tmp1);
	double icos = Icos(k+m,m);
	if (0==k%2)
	    thissum[0] += safe_mul_div (gsl_sf_gamma(k+m+mu+1),gsl_sf_gamma(mu+1),gsl_sf_gamma(k+m+1),
					(k+2*m+1.0),std::pow(t*std::abs(Q/P),k),icos,
					safe_2F1(k+m+mu+1,1,k+2*m+2,t, 3), 
					1,-1,-1,-1,1,1,1);
    }
    else if (r<=reg)
    {
	for (int el=0; el<=k; ++el)
	    thissum[0] += (0==el%2 ? 1 : -1) * 
		safe_mul_div (zeta1_func(k,el,m,mu),std::pow(eps2,el),safe_2F1(m-el-mu,-el-mu,m+1,eps1, 4), 1,1,1);
    }
    else
    {
	thissum[0] += safe_mul_div (gsl_sf_gamma(k+m+mu+1), gsl_sf_gamma(k-m+mu+1), 
				    std::pow(eps2,k), safe_2F1(m-k-mu,-k-mu,m+1,eps1, 5),
				    1,-1,1,1);
	thissum[1] += safe_mul_div (gsl_sf_gamma(k+2*m+1), gsl_sf_gamma(k+1),
				    std::pow(eps2,k), safe_2F1(-k,-k-m,m+1,eps1, 6),
				    1,-1,1,1);
    }
    
    sum[0] += thissum[0];
    sum[1] += thissum[1];
}

void I2_func (int m, double *sum, void *parms_)
{
    double *parms,hh,eps1,eps2,mu,P,Q,S,reg,r,phi,twom,thissum=0;
    
    parms = (double*)parms_;
    hh   = parms[0];
    eps1 = parms[1];
    eps2 = parms[2];
    P    = parms[3];
    Q    = parms[4];
    S    = parms[5];
    mu   = parms[6];
    reg  = parms[7];
    r    = parms[8];
    phi  = parms[9];
    twom = parms[10];

    int k=m;
    m = (int)twom;
    double rb = 1/std::sqrt(P);
    double tmp0 = r/rb;

    if (tmp0>=1./2.&&tmp0<=2)
    {
	double tmp1 = tmp0*tmp0;
	double t    = tmp1/(1+tmp1);
	double icos = Icos(k+m,m);
	if (0==k%2)
	    thissum += safe_mul_div (gsl_sf_gamma(k+m+mu+1),gsl_sf_gamma(mu+1),gsl_sf_gamma(k+m+1),
					std::pow(std::abs(Q/P),k),icos,
					safe_2F1(-k,m+mu,m+mu+1,1-t, 33), 
					1,-1,-1,1,1,1);
    }
    else if (1==m && r<1./std::sqrt(P+std::abs(Q)))
    {
	double tmpg;
	for (int el=0; el<=k; ++el)
	    if (1==el%2)
	    {
		tmpg = gsl_sf_gamma(0.5*el+1);
		thissum += safe_mul_div (std::pow(2*Q/P,el),tmpg,tmpg, 
					 gsl_sf_gamma(k-el+1),gsl_sf_gamma(el+1),gsl_sf_gamma(el+2),
					 1,1,1,-1,-1,-1);
	    }
	thissum *= std::pow(-1.,k)/((double)k)
	    *safe_mul_div (gsl_sf_gamma(mu+k+1),gsl_sf_gamma(mu+1),std::pow(P*r*r,k), 
			   1,-1,1);
    }
    else if (eps1<0.1 || eps2>0.9)
    {
	// numerical integration done instead of analytic result
	// Note: find fast converging analytic result.
    }
    else
    {
	--k;
	thissum += std::pow(eps2,k) * safe_2F1(m-k-mu,-k-mu,m+1,eps1, 7);
    }

    *sum += thissum;
}


double zeta1_func (double k, double el, double m, double mu)
{
    // Chae et al. (1998)
    return safe_mul_div (gsl_sf_gamma(k+1),       gsl_sf_gamma(k+mu+1),
		      gsl_sf_gamma(m+el+mu+1), gsl_sf_gamma(el+1),
		      gsl_sf_gamma(k-el+1),    gsl_sf_gamma(mu+1),
		      gsl_sf_gamma(el+mu+1),   gsl_sf_gamma(m+k+2),
		      1,1,1,-1,-1,-1,-1,-1);
}
double zeta2_func (double m, double mu)
{
    // Chae et al. (1998)
    return safe_mul_div (gsl_sf_gamma(m+mu),gsl_sf_gamma(m+1),gsl_sf_gamma(mu+1),
		      1,-1,-1);    
}
double associatedlegendre_func (double a1, int m, double z, double eps)
{
    // do small eps, 1st order expansion for z just bigger than unity
    if (z-1<1e-8)
    {
	double eps2;
	eps  = std::abs(eps);
	eps2 = eps*eps;
	
	return std::pow(eps,m) *
	    (
		safe_mul_div (std::pow(2.,-m),pochhammer_func(-a1, m),pochhammer_func(-a1-m,m),factorial_func(m),1,1,1,-1)
		+
		safe_mul_div (std::pow(2.,-2-m),(a1+a1*a1+m+m*m),eps2,pochhammer_func(-a1,m),pochhammer_func(-a1-m, m),
			   gsl_sf_gamma(2+m),1,1,1,1,1,-1)
		);
    }
    else
    {
	// Chae et al. (1998)
	double tmp1,tmp2,tmp3,tmp4,tmp5,tmp6;
	tmp1 = pochhammer_func (-a1,m);
	tmp2 = pochhammer_func (-a1-m,m);
	tmp3 = factorial_func (m);
	tmp4 = z-1;
	tmp5 = z+1;
	tmp6 = tmp4/tmp5;
	return safe_mul_div (tmp1,tmp2,tmp3,std::pow(tmp6,0.5*m),std::pow(0.5*tmp5,a1),
			  safe_2F1(m-a1,-a1,m+1,tmp6, 8),
			  1,1,-1,1,1,1);
    }
}



double pochhammer_func (double x, int n)
{
    return gsl_sf_poch (x,n);
}

double binomial_func (int n, int k)
{
    double ans = 1;
    for (int j=1; j<=k; ++j)
	ans *= (n+1-j)/(double)j;
    return ans;
}

double factorial_func (int n)
{
    if (n<0)
    {
	if(!errprint[en_neg_factorial])
	    std::cerr << "Warning: Factorial ((" << n << ")!) of negative number requested." << std::endl;  
	++errprint[en_neg_factorial];
        if (errprint[en_neg_factorial]>100) errprint[en_neg_factorial]=0;
	if (fatalexit) exit(1);
	return 0;
    }

    if (0==n)   return 1;
    if (n<=150) return factorial_list[n-1];

    // shouldn't get higher than 150 -- overflow could occur (64bit)
    if(!errprint[en_large_factorial])
	std::cerr << "Warning: Factorial (" << n << "!) could lead to overflow." << std::endl;  
    ++errprint[en_large_factorial];
    if (errprint[en_large_factorial]>100) errprint[en_large_factorial]=0;
    if (fatalexit) exit(1);
    
    double answer = factorial_list[150-1];
    for (int i=n;i>150;--i)
	answer *= i;
    return answer;
}

int is_negint_or_0 (double c)
{
    if (std::abs(c)<1e-3) return 1;
    if (c>0) return 0;
    if (std::abs(fmod(c,1))<1e-3) return 1;
    return 0;
}
int is_exactly_negint_or_0 (double a)
{
    if (0==a) return 1;
    if (a>0) return 0;
    if (0==std::abs(fmod(a,1))) return 1;
    return 0;
}

double safe_2F1 (double a, double b, double c, double z, int id)
{
    // This function safeguards against bad arguments (for GSL) that I think are likely to occur,
    // checks for special values, or it does a transformation for, hopefully, faster convergence.
    // Unexpected argument cause error. 
    // There is the option to link to Mathematica as a last resort, 
    // but it only works with multithreading off (if pyton's emcee is used).
    // Reference: http://functions.wolfram.com/HypergeometricFunctions/Hypergeometric2F1/17/ShowAll.html
	
    double fval=0;

    if (0==z)
	return 1;
    if (1==z && c-a-b>0)
    {
	return safe_mul_div (gsl_sf_gamma(c),gsl_sf_gamma(c-a-b),
			     gsl_sf_gamma(c-a),gsl_sf_gamma(c-b),
			     1,1,-1,-1);
    }
    if (is_exactly_negint_or_0(a) || is_exactly_negint_or_0(b))
    {
	int myint;
	double myother;
	if (is_exactly_negint_or_0(a))
	{
	    myint   = -(int)a;
	    myother = b;
	}
	else
	{
	    myint   = -(int)b;
	    myother = a;
	}

	for (int n=0; n<=myint; ++n)
	    fval += safe_mul_div (sgnpow(-1,n),binomial_func(myint,n),pochhammer_func(myother,n),
			       pochhammer_func(c,n),std::pow(z,n),  1,1,1,-1,1);
	
	return fval;
    }


    double h1=0, h2=0;
    int st1=0, st2=0;
    gsl_sf_result gslr;
    // try multiple different transformations if needed so that 0<z<=0.5
    const int numhyperfuncs = 4;
    int hyper_order[numhyperfuncs] = {-1,-1,-1,-1};
    if (z<-1)                hyper_order[0]=0;
    else if (z>=-1 && z<0)   hyper_order[0]=1;
    else if (z>=0 && z<=0.5) hyper_order[0]=2;
    else if (z>0.5 && z<=1)  hyper_order[0]=3;
 
    for (int j=0; j<numhyperfuncs; ++j)
    {
	switch (hyper_order[j])
	{
	case 0: // 1/(1-z)
	    fval = safe_mul_div (std::pow(1-z,-a),gsl_sf_gamma(c),gsl_sf_gamma(b-a),
				 gsl_sf_gamma(b),gsl_sf_gamma(c-a),safe_2F1 (a,c-b,a-b+1,1/(1-z),101),
				 1,1,1,-1,-1,1) + 
		safe_mul_div (std::pow(1-z,-b),gsl_sf_gamma(c),gsl_sf_gamma(a-b),
			      gsl_sf_gamma(a),gsl_sf_gamma(c-b),safe_2F1 (b,c-a,b-a+1,1/(1-z),102),
			      1,1,1,-1,-1,1);
	    st1 = st2 = 0;
	    break;
	case 1: // z/(z-1)
	    fval = std::pow(1-z,-a) *safe_2F1 (a,c-b,c,z/(z-1),103);
	    st1 = st2 = 0;
	    break;
	case 2: // z
	    st1  = gsl_sf_hyperg_2F1_e(a,b,c,z,&gslr);
	    h1   = gslr.val;
	    fval = h1;
	    break;
	case 3: // 1-z
	    fval = safe_mul_div (gsl_sf_gamma(c),gsl_sf_gamma(c-a-b),
				 gsl_sf_gamma(c-a),gsl_sf_gamma(c-b),safe_2F1 (a,b,a+b-c+1,1-z,104), 
				 1,1,-1,-1,1) +
		safe_mul_div (gsl_sf_gamma(c),gsl_sf_gamma(a+b-c),gsl_sf_gamma(a),
			      gsl_sf_gamma(b),std::pow(1-z,c-a-b),safe_2F1 (c-a,c-b,c-a-b+1,1-z,105), 
			      1,1,-1,-1,1,1);
	    st1 = st2 = 0;
	    break;
	default:
	    break;
	}

	if (!st1 && !st2)
	    return fval;
    }
    
    if (is_exactly_negint_or_0(c))
    {
	// do some other transformation
	
    }
    
    
    
    // if GSL failed, then check if a,b,c are too large and use recursion relations
    // force a,b,c to be positive, so that GSL can try series expansions
    // also force c=a+b at least
    
    double abs_abc[3]    = {std::abs(a),       std::abs(b),       std::abs(c)};
    int    sgn_abc[3]    = {sgn_pos_neg(a),    sgn_pos_neg(b),    sgn_pos_neg(c)};
    int    too_big[3]    = {abs_abc[0]>2||a<0, abs_abc[1]>2||b<0, abs_abc[2]>2||c<0};
    double min_abc[3]    = {fmod(a,1),         fmod(b,1),         fmod(c,1)};
    // don't allow zero values
    if (0==min_abc[0]) min_abc[0] +=sgn_abc[0];
    if (0==min_abc[1]) min_abc[1] +=sgn_abc[1];
    if (0==min_abc[2]) min_abc[2] +=sgn_abc[2];
    // force all parameters to be positive
    if (a<0) while (min_abc[0]<1.00001) ++min_abc[0];
    if (b<0) while (min_abc[1]<1.00001) ++min_abc[1];
    if (c<0) while (min_abc[2]<1.00001) ++min_abc[2];
    // force c to be large
    while (min_abc[2]<std::min(abs_abc[0],abs_abc[1])) ++min_abc[2];
    int    del_abc[3]    = {myround(std::abs(a-min_abc[0])),
			    myround(std::abs(b-min_abc[1])),
			    myround(std::abs(c-min_abc[2]))};
    int    rec_dir[3]    = {a-min_abc[0]>0,
			    b-min_abc[1]>0,
			    c-min_abc[2]>0};
    double start1_abc[3] = {too_big[0] ? min_abc[0] : a,
			    too_big[1] ? min_abc[1] : b,
			    too_big[2] ? min_abc[2] : c};
    double start2_abc[3] = {too_big[0] ? min_abc[0]+sgnpow(-1,1-rec_dir[0]) : a,
			    too_big[1] ? min_abc[1]+sgnpow(-1,1-rec_dir[1]) : b,
			    too_big[2] ? min_abc[2]+sgnpow(-1,1-rec_dir[2]) : c};
        
    double f_a1_b1_c1, f_a1_b1_c2, f_a1_b2_c1, f_a1_b2_c2;
    double f_a2_b1_c1, f_a2_b1_c2, f_a2_b2_c1, f_a2_b2_c2;
    double f_a1_b1_c, f_a1_b2_c;
    double f_a2_b1_c, f_a2_b2_c;
    double f_a1_b_c;
    double f_a2_b_c;
    double f_a_b_c;
    double f1,f2;

    // try to reduce (one by one) a/b then c, and then try all three at once if needed
    st1 = st2 = 1;

    // put largest parameter in 'a' position
    if (abs_abc[1]>abs_abc[0])
    {
	swap (abs_abc,abs_abc+1);
	swap (sgn_abc,sgn_abc+1);
	swap (too_big,too_big+1);
	swap (min_abc,min_abc+1);
	swap (del_abc,del_abc+1);
	swap (rec_dir,rec_dir+1);
	swap (start1_abc,start1_abc+1);
	swap (start2_abc,start2_abc+1);
	swap (&a,&b);
    }

    // try a
    if ((st1||st2) && too_big[0])
    {
	st1 = gsl_sf_hyperg_2F1_e (start1_abc[0],b,c,z,&gslr);
        f1  = gslr.val;
        st2 = gsl_sf_hyperg_2F1_e (start2_abc[0],b,c,z,&gslr);
        f2  = gslr.val;
        gauss_recursion (start1_abc[0], b, c, z, 0, del_abc, rec_dir, f1, f2, &f_a_b_c);
	fval = f_a_b_c;
    }
    // try b
    if ((st1||st2) && too_big[1])
    {
	st1 = gsl_sf_hyperg_2F1_e (a,start1_abc[1],c,z,&gslr);
        f1  = gslr.val;
        st2 = gsl_sf_hyperg_2F1_e (a,start2_abc[1],c,z,&gslr);
        f2  = gslr.val;
        gauss_recursion (a, start1_abc[1], c, z, 1, del_abc, rec_dir, f1, f2, &f_a_b_c);
	fval = f_a_b_c;
    }
    // try c
    if ((st1||st2) && too_big[2])
    {
	st1 = gsl_sf_hyperg_2F1_e (a, b, start1_abc[2],z,&gslr);
        f1  = gslr.val;
        st2 = gsl_sf_hyperg_2F1_e (a, b, start2_abc[2],z,&gslr);
        f2  = gslr.val;
        gauss_recursion (a, b, start1_abc[2], z, 2, del_abc, rec_dir, f1, f2, &f_a_b_c);
	fval = f_a_b_c;
    }



    // try all three at once
    // could lead to large loss in precision or terribly inacurrate result
    if (!too_big[1] && too_big[0])
    {
	swap (abs_abc,abs_abc+1);
	swap (sgn_abc,sgn_abc+1);
	swap (too_big,too_big+1);
	swap (min_abc,min_abc+1);
	swap (del_abc,del_abc+1);
	swap (rec_dir,rec_dir+1);
	swap (start1_abc,start1_abc+1);
	swap (start2_abc,start2_abc+1);
	swap (&a,&b);
    }
    if ((st1||st2) && too_big[2])
    {
	st1 = gsl_sf_hyperg_2F1_e (start1_abc[0],start1_abc[1],start1_abc[2],z,&gslr);
	f_a1_b1_c1 = gslr.val;
	st2 = gsl_sf_hyperg_2F1_e (start1_abc[0],start1_abc[1],start2_abc[2],z,&gslr);
	f_a1_b1_c2 = gslr.val;

	gauss_recursion (start1_abc[0],start1_abc[1],start1_abc[2], z, 2, del_abc, rec_dir, f_a1_b1_c1, f_a1_b1_c2, &f_a1_b1_c);

	if (too_big[1])
	{
	    st1 = gsl_sf_hyperg_2F1_e (start1_abc[0],start2_abc[1],start1_abc[2],z,&gslr);
	    f_a1_b2_c1 = gslr.val;
	    st2 = gsl_sf_hyperg_2F1_e (start1_abc[0],start2_abc[1],start2_abc[2],z,&gslr);
	    f_a1_b2_c2 = gslr.val;

	    gauss_recursion (start1_abc[0],start2_abc[1],start1_abc[2], z, 2, del_abc, rec_dir, f_a1_b2_c1, f_a1_b2_c2, &f_a1_b2_c);    
	    gauss_recursion (start1_abc[0],start1_abc[1], c, z, 1, del_abc, rec_dir, f_a1_b1_c, f_a1_b2_c, &f_a1_b_c);

	    if (too_big[0])
	    {
		st1 = gsl_sf_hyperg_2F1_e (start2_abc[0],start1_abc[1],start1_abc[2],z,&gslr);
		f_a2_b1_c1 = gslr.val;
		st2 = gsl_sf_hyperg_2F1_e (start2_abc[0],start1_abc[1],start2_abc[2],z,&gslr);
		f_a2_b1_c2 = gslr.val;

		st1 = gsl_sf_hyperg_2F1_e (start2_abc[0],start2_abc[1],start1_abc[2],z,&gslr);
		f_a2_b2_c1 = gslr.val;
		st2 = gsl_sf_hyperg_2F1_e (start2_abc[0],start2_abc[1],start2_abc[2],z,&gslr);
		f_a2_b2_c2 = gslr.val;
		
		gauss_recursion (start2_abc[0],start1_abc[1],start1_abc[2], z, 2, del_abc, rec_dir, f_a2_b1_c1, f_a2_b1_c2, &f_a2_b1_c);
		gauss_recursion (start2_abc[0],start2_abc[1],start1_abc[2], z, 2, del_abc, rec_dir, f_a2_b2_c1, f_a2_b2_c2, &f_a2_b2_c);
		gauss_recursion (start2_abc[0],start1_abc[1], c, z, 1, del_abc, rec_dir, f_a2_b1_c, f_a2_b2_c, &f_a2_b_c);
		gauss_recursion (start1_abc[0],b, c, z, 0, del_abc, rec_dir, f_a1_b_c, f_a2_b_c, &f_a_b_c);

		fval = f_a_b_c;
	    }
	    else
	    {
		fval = f_a1_b_c;
	    }
	}
	else
	{
	    fval = f_a1_b1_c;
	}
    }
    else if (st1||st2)
    {
	if (too_big[1])
	{
	    st1 = gsl_sf_hyperg_2F1_e (start1_abc[0],start1_abc[1],c,z,&gslr);
	    f_a1_b1_c = gslr.val;
	    st2 = gsl_sf_hyperg_2F1_e (start1_abc[0],start2_abc[1],c,z,&gslr);
	    f_a1_b2_c = gslr.val;

	    gauss_recursion (start1_abc[0],start1_abc[1], c, z, 1, del_abc, rec_dir, f_a1_b1_c, f_a1_b2_c, &f_a1_b_c);

	    if (too_big[0])
	    {
		st1 = gsl_sf_hyperg_2F1_e (start2_abc[0],start1_abc[1],c,z,&gslr);
		f_a2_b1_c = gslr.val;
		st2 = gsl_sf_hyperg_2F1_e (start2_abc[0],start2_abc[1],c,z,&gslr);
		f_a2_b2_c = gslr.val;
		
		gauss_recursion (start2_abc[0],start1_abc[1], c, z, 1, del_abc, rec_dir, f_a2_b1_c, f_a2_b2_c, &f_a2_b_c);
		gauss_recursion (start1_abc[0],b, c, z, 0, del_abc, rec_dir, f_a1_b_c, f_a2_b_c, &f_a_b_c);
		
		fval = f_a_b_c;
	    }
	    else
	    {
		fval = f_a1_b_c;
	    }
	}
	else
	{
	    // no recusrion relation to apply
	}
    }

    if (!st1&&!st2&&fval==fval)
    {
	if (too_big[0]||too_big[1]||too_big[2])
	    if (!errprint[en_hyper_2F1_gauss])
	    {
		std::cerr << "Warning: GSL 2F1 evaluation failed. Used recursion relations -- possible numerical errors: " 
			  << id << ": 2F1(" << a << "," << b << "," << c << "," << z << ")=" << fval << " ???" << std::endl;
		++errprint[en_hyper_2F1_gauss];
		if (errprint[en_hyper_2F1_gauss]>100) errprint[en_hyper_2F1_gauss]=0;
	    }
	
	return fval;
    }
    
    // error handling
    if (!errprint[en_hyper_2F1])
	std::cerr << "Warning: Hypergeometric function evaluation failed: " 
		  << id << " 2F1(" << a << "," << b << "," << c << "," << z << ")" << std::endl;
    ++errprint[en_hyper_2F1];
    if (errprint[en_hyper_2F1]>100) errprint[en_hyper_2F1]=0;
    if (fatalexit) exit(1);
    
    return 0;
}

void gauss_recursion (double a, double b, double c, double z, 
		      int this_index, int *del_abc, int* rec_dir,
		      double f10, double f20, double *fnew0)
{
    int ti=this_index;
    



/*
    int prec = 256;
    mpf_t f1, f2, fnew, tmp1, tmp2, tmp3, tmp4;
    mpf_init2 (f1, prec);
    mpf_init2 (f2, prec);
    mpf_init2 (fnew, prec);
    mpf_init2 (tmp1, prec);
    mpf_init2 (tmp2, prec);
    mpf_init2 (tmp3, prec);
    mpf_init2 (tmp4, prec);
    mpf_set_d (f1,f10);
    mpf_set_d (f2,f20);
    mpf_set_d (fnew,*fnew0);
    
    for (int j=2; j<=del_abc[ti]; ++j)
    {
	if (2==ti) // c parameter
	{
	    if (0==rec_dir[ti])
	    {
		mpf_set_d (tmp1, ((-(c-j)+(1-a-b+2*(c-j))*z)/((c-j)*(z-1))));
		mpf_set   (tmp2, f2);
		mpf_mul   (tmp3, tmp1, tmp2);
		mpf_set_d (tmp1, (((-1+a-(c-j))*(1-b+(c-j))*z)/((c-j)*(1+(c-j))*(z-1))));
		mpf_set   (tmp2, f1);
		mpf_mul   (tmp4, tmp1, tmp2);
		mpf_add   (fnew, tmp3, tmp4);
	    }
	    else
	    {
		mpf_set_d (tmp1, ((((c+j)-1)*(2-(c+j)-(3+a+b-2*(c+j))*z))/((1+a-(c+j))*(1+b-(c+j))*z)));
		mpf_set   (tmp2, f2);
		mpf_mul   (tmp3, tmp1, tmp2);
		mpf_set_d (tmp1, ((((c+j)-1)*((c+j)-2)*(1-z))/((1+a-(c+j))*(1+b-(c+j))*z)));
		mpf_set   (tmp2, f1);
		mpf_mul   (tmp4, tmp1, tmp2);
		mpf_add   (fnew, tmp3, tmp4);
	    }
	    mpf_set (f1,f2);
	    mpf_set (f2,fnew);
	}
	else if (0==ti) // a parameter
	{
	    if (0==rec_dir[ti])
	    {
		mpf_set_d (tmp1, ((2-c+2*(a-j)+(b-(a-j)-1)*z)/((a-j)-c+1)));
		mpf_set   (tmp2, f2);
		mpf_mul   (tmp3, tmp1, tmp2);
		mpf_set_d (tmp1, ((((a-j)+1)*(z-1))/((a-j)-c+1)));
		mpf_set   (tmp2, f1);
		mpf_mul   (tmp4, tmp1, tmp2);
		mpf_add   (fnew, tmp3, tmp4);
	    }
	    else
	    {
		mpf_set_d (tmp1, ((2-2*(a+j)+c+(-1-b+(a+j))*z)/(((a+j)-1)*(z-1))));
		mpf_set   (tmp2, f2);
		mpf_mul   (tmp3, tmp1, tmp2);
		mpf_set_d (tmp1, ((-1+(a+j)-c)/(((a+j)-1)*(z-1))));
		mpf_set   (tmp2, f1);
		mpf_mul   (tmp4, tmp1, tmp2);
		mpf_add   (fnew, tmp3, tmp4);
	    }
	    mpf_set (f1,f2);
	    mpf_set (f2,fnew);
	}
	else if (1==ti) // b parameter
	{
	    if (0==rec_dir[ti])
	    {
		mpf_set_d (tmp1, ((2-c+2*(b-j)+(a-(b-j)-1)*z)/((b-j)-c+1)));
		mpf_set   (tmp2, f2);
		mpf_mul   (tmp3, tmp1, tmp2);
		mpf_set_d (tmp1, ((((b-j)+1)*(z-1))/((b-j)-c+1)));
		mpf_set   (tmp2, f1);
		mpf_mul   (tmp4, tmp1, tmp2);
		mpf_add   (fnew, tmp3, tmp4);
	    }
	    else
	    {
		mpf_set_d (tmp1, ((2-2*(b+j)+c+(-1-a+(b+j))*z)/(((b+j)-1)*(z-1))));
                mpf_set   (tmp2, f2);
                mpf_mul   (tmp3, tmp1, tmp2);
                mpf_set_d (tmp1, ((-1+(b+j)-c))/(((b+j)-1)*(z-1)));
                mpf_set   (tmp2, f1);
                mpf_mul   (tmp4, tmp1, tmp2);
		mpf_add   (fnew, tmp3, tmp4);
	    }
	    mpf_set (f1,f2);
	    mpf_set (f2,fnew);
	}
    }

    *fnew0 = mpf_get_d (fnew);
*/

//    /*
    double f1=f10,f2=f20,fnew=0;
    for (int j=2; j<=del_abc[ti]; ++j)
    {
	if (2==ti) // c parameter
	{
	    if (0==rec_dir[ti])
		fnew = ((-(c-j)+(1-a-b+2*(c-j))*z)/((c-j)*(z-1)))*f2
		    + (((-1+a-(c-j))*(1-b+(c-j))*z)/((c-j)*(1+(c-j))*(z-1)))*f1;
	    else
		fnew = ((((c+j)-1)*(2-(c+j)-(3+a+b-2*(c+j))*z))/((1+a-(c+j))*(1+b-(c+j))*z))*f2
		    + ((((c+j)-1)*((c+j)-2)*(1-z))/((1+a-(c+j))*(1+b-(c+j))*z))*f1;
	    f1 = f2;
	    f2 = fnew;
	}
	else if (0==ti) // a parameter
	{
	    if (0==rec_dir[ti])
		fnew = ((2-c+2*(a-j)+(b-(a-j)-1)*z)/((a-j)-c+1))*f2
		    + ((((a-j)+1)*(z-1))/((a-j)-c+1))*f1;
	    else
		fnew = ((2-2*(a+j)+c+(-1-b+(a+j))*z)/(((a+j)-1)*(z-1)))*f2
		    + ((-1+(a+j)-c)/(((a+j)-1)*(z-1)))*f1;
	    f1= f2;
            f2= fnew;
	}
	else if (1==ti) // b parameter
	{
	    if (0==rec_dir[ti])
		fnew = ((2-c+2*(b-j)+(a-(b-j)-1)*z)/((b-j)-c+1))*f2
		    + ((((b-j)+1)*(z-1))/((b-j)-c+1))*f1;
	    else
		fnew = ((2-2*(b+j)+c+(-1-a+(b+j))*z)/(((b+j)-1)*(z-1)))*f2
                    + ((-1+(b+j)-c)/(((b+j)-1)*(z-1)))*f1;
	    f1= f2;
            f2= fnew;
	}
    }
    *fnew0 = fnew;
//    */
}

double Icos (int k, int m)
{
    return icos_list[k*151+m];
}

int sgn(double a)
{
    return (double(0)<a) - (a<double(0));
}
int sgn_pos_neg(double a)
{
    int mysgn = sgn(a);
    return 0==mysgn ? 1 : mysgn;
}
int sgnpow(double a, int b)
{
    return -1==sgn(a) ? (0==b%2 ? 1 : -1) : 1;
}

void infinite_sum_func (int minloops, int maxloops, int numsums, double ftol,
			void (*addends_func)(int,double*,void*), void *parms,
			double *sums, string name, time_t starttime)
{
    minloops = std::max(minloops,min_inf_sum_iter);
    maxloops = std::min(maxloops,max_inf_sum_iter);
    double prevsums3[numsums];
    double prevsums2[numsums];
    double prevsums[numsums];
    std::fill (prevsums3,prevsums3+numsums,0);
    std::fill (prevsums2,prevsums2+numsums,0);
    std::fill (prevsums,prevsums+numsums,0);
    
    int m, failed=0;
    for (m=1; m<=maxloops; ++m)
    {
      if (time(NULL)-starttime>cutoffminutes*60)
	{
	  std::fill (sums,sums+numsums,0);
	break;
	}
      /*{
	  std::cout << "KILLING BC OF TIMEOUT" << std::endl;
	  std::cerr << "KILLING BC OF TIMEOUT" << std::endl;
	    exit(1);
	}
      */
	
	std::copy (prevsums2,prevsums2+numsums,prevsums3);
	std::copy (prevsums,prevsums+numsums,prevsums2);
	std::copy (sums,sums+numsums,prevsums);
	(*addends_func) (m,sums,parms);
	
	for (int i=0; i<numsums; ++i)
	    if (!std::isfinite(sums[i]))
	    {
		std::copy (prevsums,prevsums+numsums,sums);
		failed = 1;
		break;
	    }
	if (failed) break;
	
	if (m>=minloops)
	{
	    int allgood = 1;
	    for (int i=0; i<numsums; ++i)
	    {
		if (0!=sums[i] && std::abs((prevsums[i]-sums[i])/sums[i])>=ftol)
		{
		    allgood = 0;
		    break;
		}
		if (0!=prevsums[i] && std::abs((prevsums2[i]-prevsums[i])/prevsums[i])>=ftol)
		{
		    allgood = 0;
		    break;
		}
	    }
	    if (allgood) break;
	}
    }
    
    if (m>maxloops)
    {
	if (!errprint[en_inf_sum_diverge])
	{
	    std::cerr << "Warning: Infinite sum " << name << " did not converge after " << maxloops << " iterations: ";
	    for (int j=0; j<numsums; ++j)
	    {
		std::cerr << sums[j] << "+/-" << std::max(std::abs(sums[j]-prevsums[j]),std::abs(prevsums2[j]-prevsums[j]));
		if (numsums-1!=j) std::cerr << ", " ;
	    }
	    std::cerr << std::endl;
	}
	++errprint[en_inf_sum_diverge];
        if (errprint[en_inf_sum_diverge]>100) errprint[en_inf_sum_diverge]=0;
    }
    if (failed)
    {
	if (!errprint[en_inf_sum_nan])
	{
	    std::cerr << "Warning: Infinite sum " << name << " terminated after not finite number encountered: ";
	    for (int j=0; j<numsums; ++j)
	    {
		std::cerr << sums[j] << "+/-" << std::max(std::abs(prevsums2[j]-prevsums[j]),std::abs(prevsums2[j]-prevsums3[j]));
		if (numsums-1!=j) std::cerr << ", " ;
	    }
	    std::cerr << std::endl;
	}
	++errprint[en_inf_sum_nan];
        if (errprint[en_inf_sum_nan]>100) errprint[en_inf_sum_nan]=0;
    }
    
    if (fatalexit && failed) exit(1);
}

struct st4intI0
{
  void *p;
  time_t t;
};


double I0_for_int_func (double rad, void *parms__)
{
  struct st4intI0 *s4I0 = (struct st4intI0*)parms__;
  void *parms_ = s4I0->p;
  time_t starttime;
  starttime = s4I0->t;
  
    double *parms = (double*)parms_;
    double P,Q,S,mu,reg,rcor,r1,r2,res;
    P    = parms[0];
    Q    = parms[1];
    S    = parms[2];
    mu   = parms[3];
    reg  = parms[4];
    rcor = parms[5];
    r1   = parms[6];
    r2   = parms[7];
    res  = parms[8];

    double hh,eps1,eps2;
    double tmp5,tmp6,tmp7,tmp1,tmp2,tmp3,tmp4,tmp10,tmp11,tmp12,tmp13,tmp14;
    // some terms
    tmp1  = 1+P*rad*rad;  //  rc+Prr
    tmp2  = tmp1*tmp1;          // (rc+Prr)^2
    tmp3  = Q*rad*rad;        //  Qrr
    tmp4  = tmp3*tmp3;          // (Qrr)^2
    tmp10 = tmp3/tmp1;          //  Qrr/(rc+Prr)
    tmp11 = tmp10*tmp10;        // (Qrr/(rc+Prr))^2
    tmp12 = std::sqrt(1-tmp11); // sqrt( 1 - (Qrr/(rc+Prr))^2  ) 
    tmp13 = tmp2-tmp4;          // (rc+Prr)^2 - (Qrr)^2
    tmp14 = std::sqrt(tmp13);   // sqrt( (rc+Prr)^2 - (Qrr)^2  )
    hh   = rad / tmp14;
    eps1 = (1-tmp12) / (1+tmp12);
    eps2 = 0.5* ( tmp1/tmp13 + 1/tmp14 );

    double parms2[] = {hh,eps1,eps2,P,Q,S,mu,reg,rad,0,0};

    double I0=0;
    infinite_sum_func (5, max_inf_sum_iter, 1, 1e-5, I0_func, parms2, &I0, "I0_2m_k", starttime);
    if (rad<=reg)
	I0 *= hh*std::pow(eps2,mu);
    else
	I0 = (1/(rad*std::sqrt(P*P-Q*Q))-hh*std::pow(eps2,mu)*I0)/mu;

    return I0;
}

void integrate_I0(double *parms, time_t starttime)
{
    double P,Q,S,mu,reg,rcor,r1,r2,res;
    P    = parms[0];
    Q    = parms[1];
    S    = parms[2];
    mu   = parms[3];
    reg  = parms[4];
    rcor = parms[5];
    r1   = parms[6];
    r2   = parms[7];
    res  = parms[8];
    

    if (r1==r2)
    {
	parms[8] = 0;
	return;
    }
    
    struct st4intI0 s4I0;
    s4I0.p = parms;
    s4I0.t = starttime;
    gsl_integration_workspace * w = gsl_integration_workspace_alloc (10000);
    double result, error;
    gsl_function F;
    F.function = I0_for_int_func;
    F.params = &s4I0;
    gsl_integration_qags (&F, r1, r2, 0, 1e-5, 10000, w, &result, &error);
    gsl_integration_workspace_free (w);
    
    parms[8] = result;
}




double I2_for_int_func_phi (double phivar, void *parms_)
{
    
    // Note: the cosine term is factored in by the GSL QAWO integrtion routine

    double *parms,hh,eps1,eps2,mu,P,Q,S,reg,r,phi,m,scale;
    
    parms = (double*)parms_;
    hh   = parms[0];
    eps1 = parms[1];
    eps2 = parms[2];
    P    = parms[3];
    Q    = parms[4];
    S    = parms[5];
    mu   = parms[6];
    reg  = parms[7];
    r    = parms[8];
    phi  = parms[9];
    m    = parms[10];
    scale = parms[11];

    phi = phivar;

    double tmp0 = P+Q*std::sin(2*phivar+S);

    return /*std::cos(2*m*phivar) **/ 
	( std::pow(tmp0,-(mu+1)) *safe_2F1(mu+1,mu+m,1+mu+m,-1/(r*r*tmp0), 9) )
	/scale;
}

double I2_for_int_func_r (double rvar, void *parms_)
{
    double *parms,hh,eps1,eps2,mu,P,Q,S,reg,r,phi,m,scale;
    
    parms = (double*)parms_;
    hh   = parms[0];
    eps1 = parms[1];
    eps2 = parms[2];
    P    = parms[3];
    Q    = parms[4];
    S    = parms[5];
    mu   = parms[6];
    reg  = parms[7];
    r    = parms[8];
    phi  = parms[9];
    m    = parms[10];
    scale = parms[11];

    r = rvar; 

    double tmp1,tmp2,tmp3,tmp4,tmp10,tmp11,tmp12,tmp13,tmp14;    
    tmp1  = 1+P*r*r;            //  1+Prr
    tmp2  = tmp1*tmp1;          // (1+Prr)^2
    tmp3  = Q*r*r;              //  Qrr
    tmp4  = tmp3*tmp3;          // (Qrr)^2
    tmp10 = tmp3/tmp1;          //  Qrr/(1+Prr)
    tmp11 = tmp10*tmp10;        // (Qrr/(1+Prr))^2
    tmp12 = std::sqrt(1-tmp11); // sqrt( 1 - (Qrr/(1+Prr))^2  ) 
    tmp13 = tmp2-tmp4;          // (1+Prr)^2 - (Qrr)^2
    tmp14 = std::sqrt(tmp13);   // sqrt( (1+Prr)^2 - (Qrr)^2  )

    return safe_mul_div (std::pow(r,1-2*m),std::pow(tmp13,-(mu+1)*0.5),
			 associatedlegendre_func(-(mu+1),(int)m,tmp1/tmp14,tmp3/tmp1),
			 1,1,1);
}

double massnorm_func(double r, void *parms_)
{
    double *parms = (double*)parms_;
    double Q1,Q2,P1,P2,S1,S2,N1,N2,mu;
    Q1 = parms[0];
    P1 = parms[1];
    S1 = parms[2];
    N1 = parms[3];
    Q2 = parms[4];
    P2 = parms[5];
    S2 = parms[6];
    N2 = parms[7];
    mu = parms[8];

    double A1  = 1+P1*r*r;
    double B1  = Q1*r*r;
    double mu1 = mu+1;
    double A2  = 1+P2*r*r;
    double B2  = Q2*r*r;
    double mu2 = mu+1;

    return 
	safe_mul_div (N1,std::pow(A1*A1-B1*B1,-0.5*mu1),
		      associatedlegendre_func(-mu1,0,A1/std::sqrt(A1*A1-B1*B1),B1/A1),r, 
		      1,1,1,1) +
	safe_mul_div (N2,std::pow(A2*A2-B2*B2,-0.5*mu2),
		      associatedlegendre_func(-mu2,0,A2/std::sqrt(A2*A2-B2*B2),B2/A2),r, 
		      1,1,1,1) ;
}

double safe_mul_div (double x1, double x2, double x3, double x4,
		     double x5, double x6, double x7, double x8,
		     int i1, int i2, int i3, int i4,
		     int i5, int i6, int i7, int i8)
{
    if (0==x1 || 0==x2 || 0==x3 || 0==x4 || 0==x5 || 0==x6 || 0==x7 || 0==x8)
	return 0;
    
    /*
    if (-1==i1) x1 = 1/x1;
    if (-1==i2) x2 = 1/x2;
    if (-1==i3) x3 = 1/x3;
    if (-1==i4) x4 = 1/x4;
    if (-1==i5) x5 = 1/x5;
    if (-1==i6) x6 = 1/x6;
    if (-1==i7) x7 = 1/x7;
    if (-1==i8) x8 = 1/x8;
    return x1*x2*x3*x4*x5*x6*x7*x8;
    */
    
    double ln1,ln2,ln3,ln4,ln5,ln6,ln7,ln8;
    ln1 = i1 * std::log (std::abs (x1));
    ln2 = i2 * std::log (std::abs (x2));
    ln3 = i3 * std::log (std::abs (x3));
    ln4 = i4 * std::log (std::abs (x4));
    ln5 = i5 * std::log (std::abs (x5));
    ln6 = i6 * std::log (std::abs (x6));
    ln7 = i7 * std::log (std::abs (x7));
    ln8 = i8 * std::log (std::abs (x8));
    
    int mysign = sgn(x1)*sgn(x2)*sgn(x3)*sgn(x4)*sgn(x5)*sgn(x6)*sgn(x7)*sgn(x8);
    
    return mysign * std::exp (ln1+ln2+ln3+ln4+ln5+ln6+ln7+ln8);
}
double safe_mul_div (double x1, double x2, double x3, double x4,
		     double x5, double x6, double x7,
		     int i1, int i2, int i3, int i4,
		     int i5, int i6, int i7)
{    return safe_mul_div (x1,x2,x3,x4,x5,x6,x7,1,i1,i2,i3,i4,i5,i6,i7,1); }

double safe_mul_div (double x1, double x2, double x3, double x4,
		     double x5, double x6,
		     int i1, int i2, int i3, int i4,
		     int i5, int i6)
{    return safe_mul_div (x1,x2,x3,x4,x5,x6,1,1,i1,i2,i3,i4,i5,i6,1,1); }

double safe_mul_div (double x1, double x2, double x3, double x4, double x5,
		     int i1, int i2, int i3, int i4, int i5)
{    return safe_mul_div (x1,x2,x3,x4,x5,1,1,1,i1,i2,i3,i4,i5,1,1,1); }

double safe_mul_div (double x1, double x2, double x3, double x4,
		     int i1, int i2, int i3, int i4)
{    return safe_mul_div (x1,x2,x3,x4,1,1,1,1,i1,i2,i3,i4,1,1,1,1); }

double safe_mul_div (double x1, double x2, double x3, int i1, int i2, int i3)
{    return safe_mul_div (x1,x2,x3,1,1,1,1,1,i1,i2,i3,1,1,1,1,1); }


int myround (double a)
{
    double val = fmod(a,1);
    if (val==0)
        return (int)a;
    else if (val>0.5)
        return (int)ceil(a);
    else if (a>0)
        return (int)floor(a);
    else if (std::abs(val)>0.5)
        return (int)floor(a);
    else
        return (int)ceil(a);
}
