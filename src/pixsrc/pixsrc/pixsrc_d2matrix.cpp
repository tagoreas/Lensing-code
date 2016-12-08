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




// see pixsrc_d2matrix.hpp for info

#include "pixsrc_d2matrix.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"

#include <algorithm>

template <class T>
pixsrc_d2matrix<T>::pixsrc_d2matrix()
{
    capacity1  = 0;
    extension1 = 0;

    mdm      = (T**)NULL;
    capacity = 0;
    sizearr  = 0;
}

template <class T>
pixsrc_d2matrix<T>::pixsrc_d2matrix( PS_SIT dim1 )
{
    capacity1  = dim1;
    extension1 = dim1;

    MEMORY ps_malloc( &mdm     , capacity1 );
    MEMORY ps_malloc( &capacity, capacity1 );
    MEMORY ps_malloc( &sizearr , capacity1 );
    std::fill( mdm     , mdm      + capacity1, (T*)NULL );
    std::fill( capacity, capacity + capacity1,     0    );
    std::fill( sizearr , sizearr  + capacity1,     0    );
}

template <class T>
pixsrc_d2matrix<T>::~pixsrc_d2matrix()
{
    /*
      for( PS_SIT j=0; j<extension1; ++j )
      if( mdm[j] )
      MEMORY ps_free( mdm[j] );
    */
    MEMORY ps_free( mdm, capacity1 );
    MEMORY ps_free( sizearr  );
    MEMORY ps_free( capacity );
}

template <class T>
void pixsrc_d2matrix<T>::sort( PS_SIT index )
{
    std::sort( mdm[index], mdm[index] + sizearr[index] );
}

template <class T>
void pixsrc_d2matrix<T>::pushback( PS_SIT path, T val)
{
    // check for overflow of first dimension
    if( path >= capacity1 )
    {

        if( !capacity1 )
        {
            capacity1 = (path+1)*2;
            MEMORY ps_malloc( &mdm     , capacity1 );
            MEMORY ps_malloc( &capacity, capacity1 );
            MEMORY ps_malloc( &sizearr , capacity1 );
            std::fill( mdm     , mdm      + capacity1, (T*)NULL );
            std::fill( capacity, capacity + capacity1,     0    );
            std::fill( sizearr , sizearr  + capacity1,     0    );
        }
        else
        {
            PS_SIT oldsizearr = capacity1;
            capacity1 = (path+1)*2;
            OPERA resize( &mdm     , oldsizearr, capacity1 );
            OPERA resize( &capacity, oldsizearr, capacity1 );
            OPERA resize( &sizearr , oldsizearr, capacity1 );
            std::fill( mdm      + oldsizearr, mdm      + capacity1, (T*)NULL );
            std::fill( capacity + oldsizearr, capacity + capacity1,     0    );
            std::fill( sizearr  + oldsizearr, sizearr  + capacity1,     0    );
        }

        extension1 = path+1;
    }
    else if( path >= extension1 )
    {
        extension1 = path+1;
    }

    // check if anything's been set in this index yet
    if( !mdm[path] )
    {
        // defualt size for vector is 3
        capacity[path] = 3;
        MEMORY ps_malloc( &(mdm[path]), capacity[path] );
        std::fill( mdm[path], mdm[path]+capacity[path], 0 );
    }

    // check for overflow of second dimension
    if( sizearr[path] >= capacity[path] )
    {
        PS_SIT oldsizearr = capacity[path];
        capacity[path] = (sizearr[path]+1)*2;
        OPERA resize( &(mdm[path]), oldsizearr, capacity[path] );
        std::fill( mdm[path] + oldsizearr, mdm[path] + capacity[path], 0 );
    }

    mdm[path][ sizearr[path] ] = val;
    ++sizearr[path];
}

template <class T>
void pixsrc_d2matrix<T>::set( PS_SIT *path, T val)
{
    // check for overflow of first dimension
    if( path[0] >= capacity1 )
    {
        if( !capacity1 )
        {
            capacity1 = (path[0]+1)*2;
            MEMORY ps_malloc( &mdm     , capacity1 );
            MEMORY ps_malloc( &capacity, capacity1 );
            MEMORY ps_malloc( &sizearr , capacity1 );
            std::fill( mdm     , mdm      + capacity1, (T*)NULL );
            std::fill( capacity, capacity + capacity1,     0    );
            std::fill( sizearr , sizearr  + capacity1,     0    );
        }
        else
        {
            PS_SIT oldsizearr = capacity1;
            capacity1 = (path[0]+1)*2;
            OPERA resize( &mdm     , oldsizearr, capacity1 );
            OPERA resize( &capacity, oldsizearr, capacity1 );
            OPERA resize( &sizearr , oldsizearr, capacity1 );
            std::fill( mdm      + oldsizearr, mdm      + capacity1, (T*)NULL );
            std::fill( capacity + oldsizearr, capacity + capacity1,     0    );
            std::fill( sizearr  + oldsizearr, sizearr  + capacity1,     0    );
        }

        extension1 = path[0]+1;
    }
    else if( path[0] >= extension1 )
    {
        extension1 = path[0]+1;
    }

    // check if anything's been set in this index yet
    if( !mdm[path[0]] )
    {
        capacity[path[0]] = (path[1]+1)*2;
        MEMORY ps_malloc( &(mdm[path[0]]), capacity[path[0]] );
        std::fill( mdm[path[0]], mdm[path[0]]+capacity[path[0]], 0 );
    }

    // check for overflow of second dimension
    if( path[1] >= capacity[path[0]] )
    {
        PS_SIT oldsizearr = capacity[path[0]];
        capacity[path[0]] = (path[1]+1)*2;
        sizearr[path[0]] = path[1]+1;
        OPERA resize( &(mdm[path[0]]), oldsizearr, capacity[path[0]] );
        std::fill( mdm[path[0]]+oldsizearr, mdm[path[0]]+capacity[path[0]], 0 );
    }
    else if( path[1] >= sizearr[path[0]] )
    {
        sizearr[path[0]] = path[1]+1;
    }

    mdm[ path[0] ][ path[1] ] = val;
}

template <class T>
void pixsrc_d2matrix<T>::set( PS_SIT x, PS_SIT y, T val )
{
    PS_SIT path[2] = { x, y };
    set( path, val );
}

template <class T>
T pixsrc_d2matrix<T>::get( PS_SIT x, PS_SIT y )
{
    return mdm[x][y];
}

template <class T>
T pixsrc_d2matrix<T>::get( PS_SIT *path )
{
    return mdm[ path[0] ][ path[1] ];
}
template <class T>
T* pixsrc_d2matrix<T>::get_pointer( PS_SIT x )
{
    return mdm[x];
}

template <class T>
PS_SIT pixsrc_d2matrix<T>::size( PS_SIT path )
{
    // path==-1 requests 1st dimensions size
    if( path==-1 )
        return extension1;
    return sizearr[path];
}
