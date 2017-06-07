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



#ifndef PIXSRC_MEMORY_TEMPLATES_CPP_
#define PIXSRC_MEMORY_TEMPLATES_CPP_

#include "pixsrc_memory.hpp"
#include "pixsrc_constants.hpp"
#include <stdlib.h>
#include <iostream>
#include <unistd.h>

template <typename T>
void pixsrc_memory::ps_malloc( T **ptr, size_t dim1 )
{
    size_t total = (size_t)dim1;

    *ptr = (T*)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)
    {
        return;
    }

    do
    {
        *ptr = (T*)malloc( dim1 * sizeof(T) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(double)(dim1 *sizeof(T)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep ((size_t)1e7);
        }
        else if(0)
        {
            try
            {
                MEMORY check_ptr (*ptr, dim1);
            }
            catch (...)
            {
                std::cerr << "pixsrc :: WARNING: malloc returned non-NULL pointer, but memory is invalid\n"
                          << "operating system is probably over-committing memory.\n"
                          << "trying to malloc " << 1.0e-9*(double)(dim1 *sizeof(T)) << " Gbytes\n"
                          << "waiting for memory to be available .." << std::endl;
                MEMORY ps_free (*ptr);
                *ptr = NULL;
            }
        }


#ifdef PIXSRC_MEMORYDEBUG
        pthread_mutex_lock  (&mallocmutex);
        if( *ptr )
            ++nummallocs;
        pthread_mutex_unlock(&mallocmutex);
#endif

    }
    while( !(*ptr) );
}

template <typename T>
void pixsrc_memory::ps_malloc( T ***ptr, size_t dim1 , size_t dim2 )
{
    size_t total = (size_t)(dim1*dim2);

    *ptr = (T**)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)//if( !total || total > CONSTANT uintmax/sizeof(T*) )
    {
        return;
    }

    do
    {
        *ptr = (T**)malloc( dim1 * sizeof(T*) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T*)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep ((size_t)1e7);
        }

#ifdef PIXSRC_MEMORYDEBUG
        pthread_mutex_lock  (&mallocmutex);
        if( *ptr )
            ++nummallocs;
        pthread_mutex_unlock(&mallocmutex);
#endif

    }
    while( !(*ptr) );

    for(size_t j=0; j<dim1; ++j)
    {
        ps_malloc<T>( &((*ptr)[j]), dim2 );
    }
}

template <typename T>
void pixsrc_memory::ps_malloc( T ****ptr, size_t dim1 , size_t dim2 , size_t dim3 )
{
    size_t total = (size_t)(dim1*dim2*dim3);

    *ptr = (T***)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)//if( !total || total > CONSTANT uintmax/sizeof(T**) )
    {
        return;
    }

    do
    {
        *ptr = (T***)malloc( dim1 * sizeof(T**) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T**)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep ((size_t)1e7);
        }

#ifdef PIXSRC_MEMORYDEBUG
        pthread_mutex_lock  (&mallocmutex);
        if( *ptr )
            ++nummallocs;
        pthread_mutex_unlock(&mallocmutex);
#endif

    }
    while( !(*ptr) );

    for (size_t j=0; j<dim1; ++j)
    {
        ps_malloc<T>( &((*ptr)[j]), dim2, dim3 );
    }
}

template <typename T>
void pixsrc_memory::ps_malloc( T *****ptr, size_t dim1 , size_t dim2 , size_t dim3, size_t dim4 )
{
    size_t total = (size_t)(dim1*dim2*dim3*dim4);

    *ptr = (T****)NULL;

    if(!total/* || total > CONSTANT uintmax/sizeof(T) */)//if( !total || total > CONSTANT uintmax/sizeof(T***) )
    {
        return;
    }

    do
    {
        *ptr = (T****)malloc( dim1 * sizeof(T***) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T***)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep ((size_t)1e7);
        }

#ifdef PIXSRC_MEMORYDEBUG
        pthread_mutex_lock  (&mallocmutex);
        if( *ptr )
            ++nummallocs;
        pthread_mutex_unlock(&mallocmutex);
#endif

    }
    while( !(*ptr) );

    for (size_t j=0; j<dim1; ++j)
    {
        ps_malloc<T>( &((*ptr)[j]), dim2, dim3, dim4 );
    }
}

template <typename T>
void pixsrc_memory::ps_free( T *ptr )
{

#ifdef PIXSRC_MEMORYDEBUG
    pthread_mutex_lock  (&freemutex);
    if( ptr )
        ++numfrees;
    pthread_mutex_unlock(&freemutex);
#endif

    free( ptr );
    ptr = (T*)NULL;
}

template <typename T>
void pixsrc_memory::ps_free( T **ptr, size_t dim1 )
{

#ifdef PIXSRC_MEMORYDEBUG
    pthread_mutex_lock  (&freemutex);
    if( ptr )
        ++numfrees;
    pthread_mutex_unlock(&freemutex);
#endif

    for (size_t j=0; j<dim1; ++j)
        ps_free<T>( ptr[j] );
    free( ptr );
    ptr = (T**)NULL;
}

template <typename T>
void pixsrc_memory::ps_free( T ***ptr, size_t dim1, size_t dim2 )
{

#ifdef PIXSRC_MEMORYDEBUG
    pthread_mutex_lock  (&freemutex);
    if( ptr )
        ++numfrees;
    pthread_mutex_unlock(&freemutex);
#endif

    for (size_t j=0; j<dim1; ++j)
        ps_free<T>( ptr[j], dim2 );
    free( ptr );
    ptr = (T***)NULL;
}

template <typename T>
void pixsrc_memory::ps_free( T ****ptr, size_t dim1, size_t dim2, size_t dim3 )
{

#ifdef PIXSRC_MEMORYDEBUG
    pthread_mutex_lock  (&freemutex);
    if( ptr )
        ++numfrees;
    pthread_mutex_unlock(&freemutex);
#endif

    for (size_t j=0; j<dim1; ++j)
        ps_free<T>( ptr[j], dim2, dim3 );
    free( ptr );
    ptr = (T****)NULL;
}

template <typename T>
void pixsrc_memory::ps_realloc( T **ptr, size_t dim1 )
{
    size_t total = (size_t)dim1;

    if( !total/* || total > CONSTANT uintmax/sizeof(T) */)
    {
        return;
    }

    do
    {
        *ptr = (T*)realloc( (void*)(*ptr), dim1 * sizeof(T) );

        if( !(*ptr) )
        {
            std::cerr << "pixsrc :: WARNING: trying to malloc "
                      << 1.0e-9*(dim1 *sizeof(T)) << " Gbytes\n"
                "waiting for memory to be available .." << std::endl;
            usleep ((size_t)1e7);
        }
    }
    while( !(*ptr) );
}

#endif
