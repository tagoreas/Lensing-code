//typedef long long PS_SIT;
//typedef unsigned long long PS_unsignedSIT;
typedef int PS_SIT;
typedef unsigned int PS_unsignedSIT;

typedef double PS_FPT;

#ifndef PIXSRC_MEMORY_TEMPLATES_CPP_
#define PIXSRC_MEMORY_TEMPLATES_CPP_

#include <stdlib.h>
#include <iostream>
#include <unistd.h>

template <typename T>
void ps_free( T *ptr )
{

    free( ptr );
    ptr = (T*)NULL;
}

template <typename T>
void ps_free( T **ptr, size_t dim1 )
{

    for (size_t j=0; j<dim1; ++j)
        ps_free<T>( ptr[j] );
    free( ptr );
    ptr = (T**)NULL;
}

template <typename T>
void ps_free( T ***ptr, size_t dim1, size_t dim2 )
{

    for (size_t j=0; j<dim1; ++j)
        ps_free<T>( ptr[j], dim2 );
    free( ptr );
    ptr = (T***)NULL;
}

template <typename T>
void ps_free( T ****ptr, size_t dim1, size_t dim2, size_t dim3 )
{

    for (size_t j=0; j<dim1; ++j)
        ps_free<T>( ptr[j], dim2, dim3 );
    free( ptr );
    ptr = (T****)NULL;
}

template <typename T>
void ps_realloc( T **ptr, size_t dim1 )
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

template <typename T>
void ps_malloc( T **ptr, size_t dim1 )
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
                
            }
            catch (...)
            {
                std::cerr << "pixsrc :: WARNING: malloc returned non-NULL pointer, but memory is invalid\n"
                          << "operating system is probably over-committing memory.\n"
                          << "trying to malloc " << 1.0e-9*(double)(dim1 *sizeof(T)) << " Gbytes\n"
                          << "waiting for memory to be available .." << std::endl;
                ps_free (*ptr);
                *ptr = NULL;
            }
        }


    }
    while( !(*ptr) );
}

template <typename T>
void ps_malloc( T ***ptr, size_t dim1 , size_t dim2 )
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

    }
    while( !(*ptr) );

    for(size_t j=0; j<dim1; ++j)
    {
        ps_malloc<T>( &((*ptr)[j]), dim2 );
    }
}

template <typename T>
void ps_malloc( T ****ptr, size_t dim1 , size_t dim2 , size_t dim3 )
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


    }
    while( !(*ptr) );

    for (size_t j=0; j<dim1; ++j)
    {
        ps_malloc<T>( &((*ptr)[j]), dim2, dim3 );
    }
}

template <typename T>
void ps_malloc( T *****ptr, size_t dim1 , size_t dim2 , size_t dim3, size_t dim4 )
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

    }
    while( !(*ptr) );

    for (size_t j=0; j<dim1; ++j)
    {
        ps_malloc<T>( &((*ptr)[j]), dim2, dim3, dim4 );
    }
}


#endif
