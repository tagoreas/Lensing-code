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



#ifndef PIXSRC_OPERATIONS_TEMPLATES_CPP_
#define PIXSRC_OPERATIONS_TEMPLATES_CPP_

#include "pixsrc_operations.hpp"
#include "pixsrc_memory.hpp"
#include <sstream>
#include <cmath>

template <typename T>
void pixsrc_operations::resize( T **ptr, PS_SIT oldsize, PS_SIT newsize )
{
    T *temp = *ptr;
    MEMORY ps_malloc( &(*ptr), newsize );
    std::copy( temp,
               temp + oldsize,
               *ptr );
    MEMORY ps_free( temp );
}

template <typename T>
string pixsrc_operations::tostring( T a )
{
    std::stringstream out;
    out << a;
    return out.str();
}

template <typename T>
bool pixsrc_operations::is_infinite( T num )
{
    T inf;
    OPERA assign_p_infinity( &inf );

    if( num==num && ( num >= inf || num <= -inf ) )
        return 1;
    else
        return 0;
}

template <typename T>
bool pixsrc_operations::is_exactly_p_inf( T num )
{
    T inf;
    OPERA assign_p_infinity( &inf );

    return ( num==inf ) ? 1 : 0;
}

template <typename T>
void pixsrc_operations::assign_n_infinity( T* inf )
{
    OPERA assign_p_infinity( inf );
    *inf = -*inf;
}

template <typename T>
T pixsrc_operations::slope (T x1, T y1, T x2, T y2)
{
    //returns infinity if slope is infinite.
    //slope is calculated for standard coordinate system given java's coordinate points.

    if(equalszero(x1-x2))
    {
        T inf;
        OPERA assign_p_infinity( &inf );
        return inf;
    }
    else
        return (y1-y2)/(x2-x1);
}

template <typename T>
T pixsrc_operations::intercept  (T x1, T y1, T x2, T y2)
{
    T inf;
    OPERA assign_p_infinity( &inf );

    //returns not a number if slope is infinite.
    //intercept calculated for standard coordinate system too, given java cooridnates.
    if(equalszero(x1))
        return y1;
    else if(equalszero(x2))
        return y2;
    else if(equalszero(x1-x2))
        return inf;
    else
        return -(y2+slope(x1,y1,x2,y2)*x2);
}

template <typename T>
T pixsrc_operations::xintercept (T x1, T y1, T x2, T y2)
{
    T inf;
    OPERA assign_p_infinity( &inf );

    if(equalszero(x1-x2))
        return x1;
    else if(equalszero(y1))
        return x1;
    else if(equalszero(y2))
        return x2;
    else if(equalszero(y1-y2))
        return inf;
    else
        return -intercept(x1,y1,x2,y2)/slope(x1,y1,x2,y2);
}

template <typename T>
T pixsrc_operations::convert_string( string a )
{
    T i;
    std::stringstream out(a);
    out >> i;
    return i;
}

template <typename T>
bool pixsrc_operations::is_finite( T a )
{
    T inf;
    OPERA assign_p_infinity( &inf );

    if( a==a && a < inf && a > -inf )
        return 1;
    else
        return 0;
}

template <typename T>
PS_SIT pixsrc_operations::round( T a )
{
    double val = fmod(a,1);
    if (val==0)
        return (PS_SIT)a;
    else if (val>0.5)
        return (PS_SIT)ceil(a);
    else if (a>0)
        return (PS_SIT)floor(a);
    else if (std::abs(val)>0.5)
        return (PS_SIT)floor(a);
    else
        return (PS_SIT)ceil(a);
}

template <typename T>
T pixsrc_operations::minarr (T *a, PS_SIT len, PS_SIT skip)
{
    T val;
    OPERA assign_p_infinity (&val);
    for (PS_SIT i=0; i<len; i+=skip)
        if (a[i]<val)
            val = a[i];
    return val;
}

template <typename T>
T pixsrc_operations::maxarr (T *a, PS_SIT len, PS_SIT skip)
{
    T val;
    OPERA assign_n_infinity (&val);
    for (PS_SIT i=0; i<len; i+=skip)
        if (a[i]>val)
            val = a[i];
    return val;
}

template <typename T>
void pixsrc_operations::swap (T *a, T *b)
{
    T tmp = *a;
    *a = *b;
    *b = tmp;
}

#endif
