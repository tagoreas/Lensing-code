#include <cmath>
#include <stdio.h>
#include "mem.cpp"

//typedef long long PS_SIT;
typedef int PS_SIT;

double *grid, *anggrid;
PS_SIT numpts=0,dim=0,numangles=0;
double *x_vec_avd, *x_vec_hmpr;

double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078;
double sc1=0, sc2=0, sc3=0, sc4=0, sc5=0, tr1=0, tr2=0, tr3=0, tr4=0, tr5=0;

double ang_dist (double x4, double y4, double x5, double y5)
{
    double xx4,xx5,yy4,yy5,zz4,zz5,a,b,ab,acarg,ang;
    
    x4 *= pi/180.;
    x5 *= pi/180.;
    y4 *= pi/180.;
    y5 *= pi/180.;

    xx4 = std::sqrt(1-std::cos(x5)*std::cos(x5))*std::cos(x4);
    yy4 = std::sqrt(1-std::cos(x5)*std::cos(x5))*std::sin(x4);
    zz4 = std::cos(x5);
    xx5 = std::sqrt(1-std::cos(y5)*std::cos(y5))*std::cos(y4);
    yy5 = std::sqrt(1-std::cos(y5)*std::cos(y5))*std::sin(y4);
    zz5 = std::cos(y5);

    a   = std::sqrt(xx4*xx4+yy4*yy4+zz4*zz4);
    b   = std::sqrt(xx5*xx5+yy5*yy5+zz5*zz5);
    ab  = xx4*xx5+yy4*yy5+zz4*zz5;
    acarg = ab/(a*b);
    if (acarg>1 && std::abs(acarg-1)<1e-10)
        acarg = 1;
    ang = std::acos(acarg);    
    
    return ang;
}
double r2func (double x1,double x2,double x3,double x4,double x5,
	       double y1,double y2,double y3,double y4,double y5)
{
    double r1,r2,r3,r4,r5;
    
    r1 = (x1-y1);
    r2 = (x2-y2);
    r3 = (x3-y3);
    r4 = (x4-y4);
    r5 = (x5-y5);
    
    return r1*r1*(sc1*sc1) + r2*r2*(sc2*sc2) + r3*r3*(sc3*sc3) + r4*r4*(sc4*sc4) + r5*r5*(sc5*sc5);
}

double kernel (double r2)
{
    if(r2<1e-12) 
	return (-1 + r2)* r2;
    else
	return r2*std::log(r2);
}

extern "C"
void setup_common (void *arr_, void *anggrid0_)
{
    double *arr = (double*)arr_;
    double *anggrid0 = (double*)anggrid0_;
    PS_SIT ind = 0;
    
    numangles = PS_SIT(arr[ind++]);
    numpts = PS_SIT(arr[ind++]);
    dim = numpts +5+1;
    
    ps_malloc (&grid, numpts*5);
    std::copy (arr+ind, arr+ind+5*numpts, grid);
    ind += 5*numpts;
    
    ps_malloc (&x_vec_avd,  dim*numangles);
    ps_malloc (&x_vec_hmpr, dim*numangles);
    ind += dim;

    tr1 = arr[ind++];
    tr2 = arr[ind++];
    tr3 = arr[ind++];
    tr4 = arr[ind++];
    tr5 = arr[ind++];
    sc1 = arr[ind++];
    sc2 = arr[ind++];
    sc3 = arr[ind++];
    sc4 = arr[ind++];
    sc5 = arr[ind++];
    
    for (PS_SIT vv=0; vv<numpts; ++vv)
    {
	grid[vv+numpts*0] /= sc1;
	grid[vv+numpts*0] -= tr1;
	grid[vv+numpts*1] /= sc2;
	grid[vv+numpts*1] -= tr2;
	grid[vv+numpts*2] /= sc3;
	grid[vv+numpts*2] -= tr3;
	grid[vv+numpts*3] /= sc4;
	grid[vv+numpts*3] -= tr4;
	grid[vv+numpts*4] /= sc5;
	grid[vv+numpts*4] -= tr5;
    }
    
    ps_malloc (&anggrid, numangles*2);
    std::copy (anggrid0, anggrid0+numangles*2, anggrid);
}

extern "C"
void setup_angle (int angind, void *arr1_, void *arr2_)
{
    double *arr1 = (double*)arr1_;
    double *arr2 = (double*)arr2_;
    PS_SIT ind = 0;
    ind += 2;
    ind += 5*numpts;
    std::copy (arr1+ind, arr1+ind+dim, x_vec_avd +angind*dim);
    std::copy (arr2+ind, arr2+ind+dim, x_vec_hmpr+angind*dim);
}

extern "C" 
void interpolate (void *xs_)
{
    double *xs = (double*)xs_;
    double *gridv, *gridw, *gridx, *gridy, *gridz, *t_xvec_avd, *t_xvec_hmpr;

    gridv = grid+numpts*0;
    gridw = grid+numpts*1;
    gridx = grid+numpts*2;
    gridy = grid+numpts*3;
    gridz = grid+numpts*4;
    
    double v=xs[0];
    double w=xs[1];
    double x=xs[2];
    double y=xs[3];
    double z=xs[4];
    double xyang=xs[5];
    double zang =xs[6];
    double r2, kval;
    double zval_avd[numangles];
    double zval_hmpr[numangles];
    
    for (PS_SIT ang=0; ang<numangles; ++ang)
    {
	t_xvec_avd  = x_vec_avd  + ang*dim;
	t_xvec_hmpr = x_vec_hmpr + ang*dim;
	zval_avd[ang]  = t_xvec_avd[dim-6]  + t_xvec_avd[dim-5]*(v+tr1)*sc1  + t_xvec_avd[dim-4]*(w+tr2)*sc2  + t_xvec_avd[dim-3]*(x+tr3)*sc3  + t_xvec_avd[dim-2]*(y+tr4)*sc4  + t_xvec_avd[dim-1]*(z+tr5)*sc5;
	zval_hmpr[ang] = t_xvec_hmpr[dim-6] + t_xvec_hmpr[dim-5]*(v+tr1)*sc1 + t_xvec_hmpr[dim-4]*(w+tr2)*sc2 + t_xvec_hmpr[dim-3]*(x+tr3)*sc3 + t_xvec_hmpr[dim-2]*(y+tr4)*sc4 + t_xvec_hmpr[dim-1]*(z+tr5)*sc5;
    }
    
    for (PS_SIT vv=0; vv<numpts; ++vv)
    {
	r2 = r2func(v,w,x,y,z, gridv[vv],gridw[vv],gridx[vv],gridy[vv],gridz[vv]);
	kval = kernel (r2);
	for (PS_SIT ang=0; ang<numangles; ++ang)
	{
	    t_xvec_avd  = x_vec_avd  + ang*dim;
	    t_xvec_hmpr = x_vec_hmpr + ang*dim;
	    zval_avd[ang]  += t_xvec_avd[vv]  *kval;
	    zval_hmpr[ang] += t_xvec_hmpr[vv] *kval;
	}
    }
    
    double anginterpval_avd=0,anginterpval_hmpr=0, norm=0, weight=0, dist=0;
    for (PS_SIT ang=0; ang<numangles; ++ang)
    {
	dist = ang_dist(anggrid[ang*2+0],anggrid[ang*2+1],xyang,zang);
	if (dist<1e-4)
	{
	    xs[7] = zval_avd[ang];
	    xs[8] = zval_hmpr[ang];
	    return;
	}
	weight = std::pow (dist,-4);
	norm  += weight;
	anginterpval_avd  += weight*zval_avd[ang];
	anginterpval_hmpr += weight*zval_hmpr[ang];
    }
    
    xs[7] = anginterpval_avd  /norm;
    xs[8] = anginterpval_hmpr /norm;
}


extern "C" 
void interpolatemany (void *xs_)
{
    double *xs = (double*)xs_;
    PS_SIT numpts = (PS_SIT)(xs[0]);
    double *xs1 = xs+1 + (0+numpts)*0;
    double *xs2 = xs+1 + (0+numpts)*1;
    double *xs3 = xs+1 + (0+numpts)*2;
    double *xs4 = xs+1 + (0+numpts)*3;
    double *xs5 = xs+1 + (0+numpts)*4;
    double *xs6 = xs+1 + (0+numpts)*5;
    double *xs7 = xs+1 + (0+numpts)*6;
    double *answers = xs+1 + (0+numpts)*7;
    for (PS_SIT ind=0; ind<numpts; ++ind)
    {
	double tryme[9] = {xs1[ind],xs2[ind],xs3[ind],xs4[ind],xs5[ind],xs6[ind],xs7[ind],0,0};
	interpolate (tryme);
	answers[ind]        = tryme[7];
	answers[ind+numpts] = tryme[8];
    }
}

