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



#include "pixsrc_operations_templates.cpp"
#include "pixsrc_constants.hpp"
#include "pixsrc_geometry.hpp"
#include "pixsrc_external.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer.hpp"
#include <iostream>
#include <fstream>                                      // to read in file
#include <cmath>                                        // for sqrt and fabs
#include <cctype>                                       // for isdigit
#include <sstream>
#include <unistd.h>
#include <string.h>
#include <gsl/gsl_randist.h> // for random gaussian
#include <stdlib.h>

using std::ifstream;                            // reading in file
using std::stringstream;

double pixsrc_operations::randomgaussian(const gsl_rng *r)
{
    // I'm unsure about using this function with multithreading without
    // using locks, but so far I'm using it at different times.
    return gsl_ran_gaussian (r, 1.0);
}

void pixsrc_operations::concatenate( const char **in, PS_SIT num, char **out )
{
    PS_SIT size[num];
    PS_SIT totalsize = 0;

    for( PS_SIT x=0; x<num; ++x )
    {
        size[x] = OPERA sizestring(in[x]);
        totalsize += size[x];
    }

    MEMORY ps_malloc( &(*out), totalsize+1 );

    char *startcopy = *out;
    for( PS_SIT x=0; x<num; ++x )
    {
        strcpy( startcopy, in[x] );
        startcopy += size[x];
    }
}

bool pixsrc_operations::fileexists (string filename)
{
    FILE *f = fopen (filename.c_str(), "r");
    if (NULL!=f)
    {
        fclose (f);
        return 1;
    }
    return 0;
}

void pixsrc_operations::readfile(const char *filename, char ***vec, PS_SIT *numlines, pthread_mutex_t *printlock)
{
    // removes comments signaled with a '#' and whitespace

    // open file
    ifstream file;
    file.open(filename);
    if(file.fail())
    {
        MEMORY ps_malloc( &(*vec), 0 );
        *numlines = 0;
        return;
    }

    char *buffer, *ptr;
    PS_SIT size;
    MEMORY ps_malloc( &buffer, 100000 );

    // find out how many "good" lines in file
    *numlines = 0;
    while(!file.eof())
    {
        file.getline( buffer, 100000 );

        // can't just call sizestring because need to exclude comments
        size = 0;
        ptr = buffer;

        // remove comments from string
        while( *ptr )
        {
            if( *ptr=='#' )
            {
                *ptr = 0;
                break;
            }

            ++size;
            ++ptr;
        }

        // remove white space
        OPERA trim( buffer, &size );

        if( size )
        {
            ++*numlines;
        }
    }
    file.close();

    if( !*numlines )
    {
        MEMORY ps_free( buffer );
        PRINTER printwarning("pixsrc",
                             string(filename) + " loaded but no data found",
                             printlock);
        return;
    }

    // alloctae memory for file read in
    MEMORY ps_malloc( &(*vec), *numlines );

    // rewind file buffer to beginning of file
    // file.seekg( 0, std::ios::beg );
    file.open( filename );

    // actually read in file to memory
    PS_SIT position = 0;
    while(!file.eof())
    {
        file.getline( buffer, 100000 );

        // can't just call sizestring because need to exclude comments
        size = 0;
        ptr = buffer;

        // remove comments from string
        while( *ptr )
        {
            if( *ptr=='#' )
            {
                *ptr = 0;
                break;
            }

            ++size;
            ++ptr;
        }

        // remove white space
        OPERA trim( buffer, &size );

        // write to memory
        if( size )
        {
            MEMORY ps_malloc( &((*vec)[position]), size+1 );
            strcpy( (*vec)[position], buffer );
            ++position;
        }
    }
    file.close();

    PRINTER print2screen("pixsrc",
                         string(filename) + " loaded",
                         printlock);

    MEMORY ps_free( buffer );
}
void pixsrc_operations::trim( char *line, PS_SIT *size0 )
{
    // size does NOT include terminating null character

    PS_SIT size = OPERA sizestring( line );

    // remove whitespace at end of string
    while( size && ( line[size-1]=='\n' || line[size-1]=='\r' ||
                     line[size-1]==' '  || line[size-1]=='\t'    ) )
    {
        line[size-1] = 0;
        --size;
    }

    // remove whitespace at beginning of string
    while( size && ( line[  0   ]=='\n' || line[  0   ]=='\r' ||
                     line[  0   ]==' '  || line[  0   ]=='\t'    ) )
    {
        // shift entire string over, not copying first character
        char *ptr1 = line;
        char *ptr2 = ptr1+1;
        for( PS_SIT j=0; j<size; ++j)
        {
            *ptr1 = *ptr2;
            ptr1 = ptr2;
            ++ptr2;
        }
        --size;
    }

    if( size0 )
    {
        *size0 = size;
    }
}
PS_SIT pixsrc_operations::sizestring( const char *str )
{
    // size does NOT include terminating character

    PS_SIT size = 0;
    const char *ptr = str;
    while( *ptr )
    {
        ++size;
        ++ptr;
    }
    return size;
}
void pixsrc_operations::split( const char *line0, const char *delimiter0, char ***out, PS_SIT *num )
{
    // this function will trim whitespace from ends of string before splitting it
    // a dleimiter of " " will match spaces and tabs
    const char *delimiter = strcmp (delimiter0, " ") ? delimiter0 : " \t";

    // copy string to split
    PS_SIT sizestr = OPERA sizestring( line0 );
    char line[sizestr+1];
    strcpy( line, line0 );

    // trim whitespace
    OPERA trim( line, &sizestr );

    // find number of strings result after splitting
    char *ptr;
    ptr = strtok( line, delimiter );
    PS_SIT num0 = 0;
    while( ptr )
    {
        ++num0;
        ptr = strtok( NULL, delimiter );
    }
    if( num )
        *num = num0;

    MEMORY ps_malloc( &(*out), num0 );

    // split strings

    // re-copy because strtok destroys string
    strcpy( line, line0 );
    ptr = strtok( line, delimiter );
    PS_SIT position = 0;
    while( ptr )
    {
        MEMORY ps_malloc( &((*out)[position]), OPERA sizestring(ptr)+1 );
        strcpy( (*out)[position], ptr );
        ++position;

        ptr = strtok( NULL, delimiter );
    }
}

double pixsrc_operations::distance   (double x1, double y1, double x2, double y2)
{
    return sqrt((x1-x2)*(x1-x2) + (y1-y2)*(y1-y2));
}

double pixsrc_operations::distance2  (double x1, double y1, double x2, double y2)
{
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2);
}

bool pixsrc_operations::isanumberint(string temp)
{
    PS_SIT littlee = temp.find("e");
    PS_SIT bige = temp.find("E");
    if(!equalsnpos(littlee) && !equalsnpos(bige))
        return false;
    if((!equalsnpos(littlee) || !equalsnpos(bige)))
    {
        char delimiter[2] = "e";
        if(littlee == -1)
            delimiter[0] = 'E';

        char **entries;
        PS_SIT numentries;
        OPERA split( temp.c_str(), delimiter, &entries, &numentries );

        bool answer;

        if(numentries!=2)
            answer  = false;
        if(isanumberint(entries[0]) && isanumberint(entries[1]))
            answer  = true;
        else answer = false;

        MEMORY ps_free( entries, numentries );

        return answer;
    }

    if(temp[0]=='-' || temp[0]=='+')
        temp = temp.substr(1);
    if(temp.size()==0)
        return false;

    string::iterator iter;
    string::iterator iterEnd =temp.end();
    bool afterDecimal=false;
    for(iter=temp.begin(); iter != iterEnd; iter++)
    {
        if(!isdigit(*iter) && *iter!='.')
            return false;
        else if(afterDecimal && *iter!='0')
            return false;
        else if(*iter=='.')
        {
            if(afterDecimal)
                return false;
            else
                afterDecimal = true;
        }
    }
    return true;
}
bool pixsrc_operations::isanumberfloat(string temp)
{
    PS_SIT littlee = temp.find("e");
    PS_SIT bige = temp.find("E");
    if(!equalsnpos(littlee) && !equalsnpos(bige))
        return false;
    if((!equalsnpos(littlee) || !equalsnpos(bige)))
    {
        char delimiter[2] = "e";
        if(littlee == -1)
            delimiter[0] = 'E';

        char **entries;
        PS_SIT numentries;
        OPERA split( temp.c_str(), delimiter, &entries, &numentries );

        bool answer;

        if(numentries!=2)
            answer = false;
        if(isanumberfloat(entries[0]) && isanumberfloat(entries[1]))
            answer = true;
        else answer = false;

        MEMORY ps_free( entries, numentries );

        return answer;
    }

    if(temp[0]=='-' || temp[0]=='+')
        temp = temp.substr(1);
    if(temp.size()==0)
        return false;

    string::iterator iter;
    string::iterator iterEnd =temp.end();
    bool afterDecimal=false;
    for(iter = temp.begin(); iter != iterEnd; iter++)
    {
        if(!isdigit(*iter) && *iter!='.')
            return false;
        else if(afterDecimal && *iter=='.')
            return false;
        else if(*iter=='.')
            afterDecimal=true;
    }
    return true;
}

PS_SIT pixsrc_operations::rank   (double x ,double y ,double x1,double y1,double x2,double y2)
{
    if((x1-x)*(x1-x)+(y1-y)*(y1-y)<=(x2-x)*(x2-x)+(y2-y)*(y2-y))
        return 0;
    return 1;
}

PS_SIT pixsrc_operations::coords2coordl(PS_SIT x, PS_SIT newsize)
{
    //used in blurIt()
    return x - (newsize-1)/2;
}
PS_SIT pixsrc_operations::coordl2coords(PS_SIT x, PS_SIT newsize)
{
    //used in blurIt()
    return x+(newsize-1)/2; // this is the inversion from coordS2coordL
}

void pixsrc_operations::linearinterpolatorcartesian(double x, double y, double reduction, double *positions, PS_SIT *exists, double *result)
{
    // given a point Q whose value must be interpolated from 4 points that surround
    // Q and form a square, this returns the weights given to the 4 points, ordered
    // counterclockwise starting in the first quadrant.

    std::fill( result, result+4, 0 );
    double *weights = result;

    if(positions[0*2+1]==y)
    {
        if(exists[0]!=-1 && exists[1]!=-1)
        {
            weights[0] =  distance(positions[1*2],positions[1*2+1],x,y) /
                distance(positions[0*2],positions[0*2+1],
                         positions[1*2],positions[1*2+1]);
            weights[1] = 1.0-weights[0];
        }
    }
    else if(positions[1*2]==x)
    {
        if(exists[2]!=-1 && exists[1]!=-1)
        {
            weights[1] = distance(positions[2*2],positions[2*2+1],x,y) /
                distance(positions[1*2],positions[1*2+1],
                         positions[2*2],positions[2*2+1]);
            weights[2] = 1.0-weights[1];
        }
    }
    else if(positions[1*2]+reduction/2.0==x && positions[0*2+1]+reduction/2.0==y)
    {
        if(exists[0]!=-1 && exists[1]!=-1 && exists[2]!=-1 && exists[3]!=-1)
        {
            weights[0]=0.25;
            weights[1]=0.25;
            weights[2]=0.25;
            weights[3]=0.25;
        }
    }
    else if( GEOM isonlinesegment(positions[0*2],positions[0*2+1],
                                  positions[2*2],positions[2*2+1],x,y) )
    {
        if(exists[0]!=-1 && exists[2]!=-1)
        {
            weights[0] = distance(positions[2*2],positions[2*2+1],x,y) /
                distance(positions[0*2],positions[0*2+1],
                         positions[2*2],positions[2*2+1]);
            weights[2] = 1.0-weights[0];
        }
    }
    else if( GEOM isonlinesegment(positions[1*2],positions[1*2+1],
                                  positions[3*2],positions[3*2+1],x,y) )
    {
        if(exists[3]!=-1 && exists[1]!=-1)
        {
            weights[1] = distance(positions[3*2],positions[3*2+1],x,y) /
                distance(positions[1*2],positions[1*2+1],
                         positions[3*2],positions[3*2+1]);
            weights[3] = 1.0-weights[1];
        }
    }
    else
    {
        PS_SIT pos[2] = {2,3};
        if(GEOM isintri(positions  ,x,y))
            pos[0]=0;
        if(GEOM isintri(positions+2,x,y))
            pos[1]=1;

        PS_SIT test=pos[0]+1;
        if(test==4) test=0;
        bool above=false;
        if(test==pos[1]) above=true;

        PS_SIT diff1,diff2;
        if(above)
        {
            diff1=pos[0];
            diff2=pos[1]+2;
            if(diff2==5) diff2=1;
            else if(diff2==4) diff2=0;
        }
        else
        {
            diff1=pos[0]+2;
            diff2=pos[1];
            if(diff1==5) diff1=1;
            else if(diff1==4) diff1=0;
        }

        PS_SIT posF[3]={pos[0],pos[0]+1,pos[0]+2};
        if(rank( x, y, positions[diff1*2], positions[diff1*2+1],
                 positions[diff2*2], positions[diff2*2+1]) == 1 )
        {
            posF[0]=pos[1];
            posF[1]=pos[1]+1;
            posF[2]=pos[1]+2;
        }
        if(posF[1]==4) posF[1]=0;
        if(posF[2]==4) posF[2]=0;
        else if(posF[2]==5) posF[2]=1;

        double x1=positions[posF[0]*2  ];
        double y1=positions[posF[0]*2+1];
        double x2=positions[posF[1]*2  ];
        double y2=positions[posF[1]*2+1];
        double x3=positions[posF[2]*2  ];
        double y3=positions[posF[2]*2+1];


        double del23 = x2-x3;
        double del31 = x3-x1;
        double del12 = x1-x2;
        double denom = del23*y1 + del31*y2 + del12*y3;
        double w1=(del23*y + (x3-x)*y2 + (x-x2)*y3)/denom;
        double w3=(del12*y + (x2-x)*y1 + (x-x1)*y2)/denom;
        double w2=1.0-w1-w3;

        if(exists[posF[0]]!=-1 && exists[posF[1]]!=-1 && exists[posF[2]]!=-1)
        {
            weights[posF[0]] = w1;
            weights[posF[1]] = w2;
            weights[posF[2]] = w3;
        }
    }
}
double pixsrc_operations::getfitscardvalue(string s)
{
    char **entries1,  **entries2;
    PS_SIT numentries1, numentries2;

    OPERA split( s.c_str()  , "=", &entries1, &numentries1 );
    OPERA split( entries1[1], "/", &entries2, &numentries2 );

    double answer = OPERA convert_string <double> ( entries2[0] );

    MEMORY ps_free( entries1, numentries1 );
    MEMORY ps_free( entries2, numentries2 );

    return answer;
}

void pixsrc_operations::planarinterpolation3pts( double x, double y,
                                                 double p[2][3], double weights[3] )
{
    std::fill( weights, weights + 3, 0 );

    double del23 = p[0][1] - p[0][2];
    double del31 = p[0][2] - p[0][0];
    double del12 = p[0][0] - p[0][1];
    double denom = del23*p[1][0] + del31*p[1][1] + del12*p[1][2];

    weights[0] = (del23*y + (p[0][2]-x)*p[1][1] + (x-p[0][1])*p[1][2])/denom;
    weights[2] = (del12*y + (p[0][1]-x)*p[1][0] + (x-p[0][0])*p[1][1])/denom;
    weights[1] = 1.0 - weights[2] - weights[0];
}

void pixsrc_operations::planarvalueinterpolation( double pos[2][3], double val[3],
                                                  double x, double y, double *result )
{
    double del23 = pos[0][1]-pos[0][2];
    double del31 = pos[0][2]-pos[0][0];
    double del12 = pos[0][0]-pos[0][1];
    double denom = del23*pos[1][0] + del31*pos[1][1] + del12*pos[1][2];

    double weights[3];
    weights[0] = (del23*y + (pos[0][2]-x)*pos[1][1] + (x-pos[0][1])*pos[1][2])/denom;
    weights[2] = (del12*y + (pos[0][1]-x)*pos[1][0] + (x-pos[0][0])*pos[1][1])/denom;
    weights[1] = 1.0 - weights[2] - weights[0];

    *result = weights[0]*val[0] + weights[1]*val[1] + weights[2]*val[2];
}

bool pixsrc_operations::equalsnpos(PS_SIT checkit)
{
    if(checkit==-1 || checkit>100000)
        return 1;
    return 0;
}

void pixsrc_operations::assign_p_infinity( double* inf )
{
    *inf = CONSTANT d_inf;
}

void pixsrc_operations::assign_p_infinity( float* inf )
{
    *inf = CONSTANT s_inf;
}

bool pixsrc_operations::equalszero( double a )
{
    return ( std::abs(a) <= CONSTANT d_epsilon ) ? 1 : 0;
}

bool pixsrc_operations::equalszero( float a )
{
    return ( std::abs(a) <= CONSTANT s_epsilon ) ? 1 : 0;
}

void pixsrc_operations::pthreadswait(lensvar *vars_, PS_SIT dim, PS_SIT *waitlist)
{
    char compound;
    for(PS_SIT i=0; i<dim; ++i)
    {
        do
        {
            pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));
            compound = vars_->pthreadstracker[waitlist[i]];
            pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));

            if( compound!=2 )
            {
                // sleep 1 millisecond
                usleep( CONSTANT shortwaitms );
            }
        }
        while( compound!=2 );
    }


    /*
      char compound;
      do
      {
      compound = 1;
      for(PS_SIT i=0; i<dim; ++i)
      {
      pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));
      compound = compound && vars_->pthreadstracker[waitlist[i]];
      pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));
      }
      }
      while(!compound);

      for(PS_SIT i=0; i<dim; ++i)
      {
      pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));

      if( vars_->pthreadstracker[waitlist[i]] == 1 )
      {
      pthread_join(vars_->pthreads[waitlist[i]],NULL);
      vars_->pthreadstracker[waitlist[i]] = 2;
      }

      pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));
      }
    */
}
void pixsrc_operations::pthreadswaitifstarted(lensvar *vars_, PS_SIT dim, PS_SIT *waitlist)
{
    char compound = 1;
    for(PS_SIT i=0; i<dim; ++i)
    {
        pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));
        compound = compound && vars_->pthreadstracker[waitlist[i]];
        pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));
    }

    if( !compound )
        usleep( CONSTANT longwaitms );

    for(PS_SIT i=0; i<dim; ++i)
    {
        do
        {
            pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));
            compound = vars_->pthreadstracker[waitlist[i]];
            pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));

            if( compound==1 )
            {
                usleep( CONSTANT shortwaitms );
            }
        }
        while( compound==1 );
    }


    /*
      char compound = 1;
      for(PS_SIT i=0; i<dim; ++i)
      {
      pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));
      compound = compound && vars_->pthreadstracker[waitlist[i]];
      pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));
      }

      if( !compound )
      usleep( CONSTANT longwaitms );

      for(PS_SIT i=0; i<dim; ++i)
      {
      if( vars_->pthreadstracker[waitlist[i]] )
      {
      pthread_mutex_lock(&(vars_->pthreadslock[waitlist[i]]));

      if( vars_->pthreadstracker[waitlist[i]] == 1 )
      {
      pthread_join(vars_->pthreads[waitlist[i]],NULL);
      vars_->pthreadstracker[waitlist[i]] = 2;
      }

      pthread_mutex_unlock(&(vars_->pthreadslock[waitlist[i]]));
      }
      }
    */
}

void pixsrc_operations::multidimmin (inputdata *data_, commoninputdata *cdata_,
                                     PS_SIT nvary, double tolerance, string optname,
                                     gsl_vector *x, gsl_vector *ss,
                                     void *params, double (*penfunc)(const gsl_vector*,void*),
                                     double *stat, PS_SIT maxiter, PS_SIT printfreq)
{
    // nvary is number of parameters being optimized
    // tolerance is passed to minimizer to check for convergence
    // optname is a short description of the minimization
    // x contains initial guess
    // ss are stepsizes
    // params will be passed to the minimizing function
    // penfunc is the minimizing function
    // stat will contain the lowest function evaluation upon exit (can be null)
    // maxiter is the maximum number of iterations
    // printfeq is the frequency (in units of 1/iterations) at which the status is printed

    // struct that contains default pixsrc setup for GSL minimizer
    struct gsl_mdm mdm;

    mdm.minex_func.n      = nvary;
    mdm.minex_func.f      = penfunc;
    mdm.minex_func.params = params;
    mdm.s                 = gsl_multimin_fminimizer_alloc (mdm.T, mdm.minex_func.n);

    gsl_multimin_fminimizer_set (mdm.s, &mdm.minex_func, x, ss);
    do
    {
        if ((PS_SIT)mdm.iter>maxiter)
        {
            PRINTER printwarning (data_->print2screenname,
                                  optname + ": exceeded max " +
                                  OPERA tostring(maxiter)+" iterations",
                                  cdata_->print2screenmutex);
            break;
        }

        ++mdm.iter;
        mdm.status = gsl_multimin_fminimizer_iterate (mdm.s);

        if (mdm.status)
            break;

        mdm.size   = gsl_multimin_fminimizer_size (mdm.s);
        mdm.status = gsl_multimin_test_size (mdm.size, tolerance);

        if (data_->verbose!=1 && (printfreq && (1==printfreq||mdm.iter%printfreq==1)))
            PRINTER print2screen (data_->print2screenname,
                                  optname + ": iter "
                                  + OPERA tostring (mdm.iter) + ", fval "
                                  + OPERA tostring (mdm.s->fval) + ", size "
                                  + OPERA tostring(mdm.size),
                                  cdata_->print2screenmutex);
    }
    while (mdm.status == GSL_CONTINUE);

    if(data_->verbose!=1)
    {
        if (mdm.status == GSL_SUCCESS)
            PRINTER print2screen (data_->print2screenname,
                                  optname + " converged after "
                                  + OPERA tostring(mdm.iter)
                                  + " iter's with size " + OPERA tostring(mdm.size),
                                  cdata_->print2screenmutex);
        else
            PRINTER printwarning (data_->print2screenname,
                                  optname + ": optimization failed",
                                  cdata_->print2screenmutex);
    }

    // copy best value into x
    gsl_vector *best = gsl_multimin_fminimizer_x (mdm.s);
    for (PS_SIT nv=0; nv<(PS_SIT)mdm.minex_func.n; ++nv)
        gsl_vector_set (x, nv, gsl_vector_get (best, nv));
    if (stat)
        *stat = gsl_multimin_fminimizer_minimum (mdm.s);

    // print it
    if(data_->verbose!=1)
    {
        string beststr = "";
        for (PS_SIT nv=0; nv<(PS_SIT)mdm.minex_func.n; ++nv)
            beststr += OPERA tostring (gsl_vector_get(x, nv)) + " ";
        PRINTER print2screen (data_->print2screenname,
                              optname + ": " + beststr,
                              cdata_->print2screenmutex);
    }

    gsl_multimin_fminimizer_free (mdm.s);
}

void pixsrc_operations::shuffle (PS_SIT *array, PS_SIT n)
{
    if (n > 1) 
    {
        for (PS_SIT i=0; i<n-1; ++i) 
        {
	    PS_SIT j = i + rand() / (RAND_MAX / (n-i)+1);
	    PS_SIT t = array[j];
	    array[j] = array[i];
	    array[i] = t;
        }
    }
}

/*
  void pixsrc_operations::intersection(double ax,double ay,double bx,double by
  double cx,double cy,double dx,double dy, vector<double> &result)
  {
  result[0]=result[1]=CONSTANT inf;

  double m0 = (by-ay)/(bx-ax);
  double m1 = (dy-cy)/(dx-cx);
  double b0 = ay-m0*ax;
  double b1 = cy-m1*cx;

  if((equalszero(bx-ax) && equalszero(dx-cx)) || equalszero(m0-m1))
  {
  result[0] = CONSTANT inf;
  result[1] = CONSTANT inf;
  }
  else if(equalszero(bx-ax))
  {
  result[0] = bx;
  result[1] = m1*result[0]+b1;

  }
  else if(equalszero(dx-cx))
  {
  result[0] = cx;
  result[1] = m0*result[0]+b0;
  }
  else
  {
  result[0] = (b1-b0)/(m0-m1);
  result[1] = m0*result[0]+b0;
  }
  }
*/

/*
  void pixsrc_operations::quadfit(vector<double> ex, vector<double> why, vector<double> &result)
  {
  result[0]=result[1]=result[2]=CONSTANT inf;
  PS_SIT N=why.size();
  if(N<3)
  return;
  double sumx=0, sumy=0, sumxy=0, sumx2=0, sumx3=0, sumx4=0, sumx2y=0;

  for(PS_SIT x = 0; x < N; x++)
  {
  sumx += ex[x];
  sumy += why[x];
  sumxy += ex[x]*why[x];
  double squ = ex[x]*ex[x];
  sumx2 += squ;
  sumx3 += squ*ex[x];
  sumx4 += squ*squ;
  sumx2y += squ*why[x];
  }

  double denom = N*sumx2*sumx4 + 2*sumx*sumx2*sumx3 - sumx2*sumx2*sumx2 - sumx*sumx*sumx4 - N*sumx3*sumx3;
  if(equalszero(denom))
  return;
  result[0] = (N*sumx2*sumx2y + sumx*sumx3*sumy + sumx*sumx2*sumxy - sumx2*sumx2*sumy - sumx*sumx*sumx2y - N*sumx3*sumxy)/denom;
  result[1] = (N*sumx4*sumxy + sumx*sumx2*sumx2y + sumx2*sumx3*sumy - sumx2*sumx2*sumxy - sumx*sumx4*sumy - N*sumx3*sumx2y)/denom;
  result[2] = (sumx2*sumx4*sumy + sumx2*sumx3*sumxy + sumx*sumx3*sumx2y - sumx2*sumx2*sumx2y - sumx*sumx4*sumxy - sumx3*sumx3*sumy)/denom;
  }
*/

/*
  void pixsrc_operations::planarinterpolation4pts(double x, double y, double p[2][4], double weights[4])
  {
  std::fill(weights,weights+4,0);

  bool isInTriangle[4];
  isInTriangle[0] = GEOM isintriangle(p[0][0],p[1][0],p[0][1],p[1][1],p[0][3],p[1][3],x,y); // top right triangle
  isInTriangle[1] = GEOM isintriangle(p[0][0],p[1][0],p[0][2],p[1][2],p[0][3],p[1][3],x,y); // bottom right triangle
  isInTriangle[2] = GEOM isintriangle(p[0][2],p[1][2],p[0][1],p[1][1],p[0][3],p[1][3],x,y); // bottom left triangle
  isInTriangle[3] = GEOM isintriangle(p[0][0],p[1][0],p[0][1],p[1][1],p[0][2],p[1][2],x,y); // top left triangle

  PS_SIT surrPix[3];
  if(isInTriangle[0] || isInTriangle[1] || isInTriangle[2] || isInTriangle[3])
  GEOM norm_parallelogram_ic(p[0][0],p[1][0],p[0][1],p[1][1],p[0][2],p[1][2],p[0][3],p[1][3],isInTriangle[0],isInTriangle[1],isInTriangle[2],isInTriangle[3],&surrPix[0]);
  else
  return;

  double del23 = p[0][surrPix[1]] - p[0][surrPix[2]];
  double del31 = p[0][surrPix[2]] - p[0][surrPix[0]];
  double del12 = p[0][surrPix[0]] - p[0][surrPix[1]];
  double denom = del23*p[1][surrPix[0]] + del31*p[1][surrPix[1]] + del12*p[1][surrPix[2]];

  weights[surrPix[0]] = (del23*y + (p[0][surrPix[2]]-x)*p[1][surrPix[1]] + (x-p[0][surrPix[1]])*p[1][surrPix[2]])/denom;
  weights[surrPix[2]] = (del12*y + (p[0][surrPix[1]]-x)*p[1][surrPix[0]] + (x-p[0][surrPix[0]])*p[1][surrPix[1]])/denom;
  weights[surrPix[1]] = 1.0-weights[surrPix[2]]-weights[surrPix[0]];
  }
*/

/*
  double pixsrc_operations::getweight(PS_SIT x, PS_SIT y, double **w, PS_SIT max, bool fromfile)
  {
  // used in blurIt()
  if(fromfile)
  return w[x+max][y+max];

  if(y>=0)
  return w[x+max][y];
  return w[-x+max][-y];
  }
*/

/*
  double pixsrc_operations::rank(double x0,double y0,double x1,double y1,
  double x2,double y2,double x3,double y3)
  {
  double inter[4];
  inter[0] =  intercept(x0,y0,x1,y1);
  inter[1] = xintercept(x1,y1,x2,y2);
  inter[2] =  intercept(x2,y2,x3,y3);
  inter[3] = xintercept(x3,y3,x0,y0);

  return inter[0]*inter[0]+inter[1]*inter[1]+inter[2]*inter[2]+inter[3]*inter[3];
  //return Math.max(Math.abs(inter[0]/inter[2]),Math.abs(inter[2]/inter[0])) + Math.max(Math.abs(inter[1]/inter[3]),Math.abs(inter[3]/inter[1]));
  }
*/
