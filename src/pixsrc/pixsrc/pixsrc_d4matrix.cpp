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




// see pixsrc_d4matrix.hpp for info

#include "pixsrc_d4matrix.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"

template <class T>
pixsrc_d4matrix<T>::pixsrc_d4matrix()
{
    capacity1  = 0;
    extension1 = 0;
    capacity2  = (T*   )NULL;
    extension2 = (T*   )NULL;
    capacity3  = (T**  )NULL;
    extension3 = (T**  )NULL;
    capacity4  = (T*** )NULL;
    extension4 = (T*** )NULL;
    mdm        = (T****)NULL;
}

template <class T>
pixsrc_d4matrix<T>::pixsrc_d4matrix( PS_SIT dim1 )
{
    capacity1  = dim1;
    extension1 = dim1;
    capacity2  = (T*   )NULL;
    extension2 = (T*   )NULL;
    capacity3  = (T**  )NULL;
    extension3 = (T**  )NULL;
    capacity4  = (T*** )NULL;
    extension4 = (T*** )NULL;

    MEMORY ps_malloc( &mdm       , capacity1 );
    MEMORY ps_malloc( &capacity2 , capacity1 );
    MEMORY ps_malloc( &extension2, capacity1 );
    std::fill( mdm       , mdm        + capacity1, (T***)NULL );
    std::fill( capacity2 , capacity2  + capacity1,      0     );
    std::fill( extension2, extension2 + capacity1,      0     );
}

template <class T>
pixsrc_d4matrix<T>::~pixsrc_d4matrix()
{
    for( PS_SIT i=0; i<capacity1; ++i )
    {
        if( capacity2[i] )
        {
            MEMORY ps_free(  capacity4[i], capacity2[i] );
            MEMORY ps_free( extension4[i], capacity2[i] );

            for( PS_SIT j=0; j<capacity2[i]; ++j )
                if( capacity3[i][j] )
                    MEMORY ps_free( mdm[i][j], capacity3[i][j] );
        }
    }

    MEMORY ps_free(        mdm, capacity1 );
    MEMORY ps_free(  capacity4            );
    MEMORY ps_free( extension4            );
    MEMORY ps_free(  capacity3, capacity1 );
    MEMORY ps_free( extension3, capacity1 );
    MEMORY ps_free(  capacity2            );
    MEMORY ps_free( extension2            );
}

template <class T>
void pixsrc_d4matrix<T>::set( PS_SIT *path, T val)
{
    // START FIRST DIMENSION
    {
        // check for overflow of first dimension
        if( path[0] >= capacity1 )
        {
            if( !capacity1 )
            {
                capacity1 = (path[0]+1)*2;
                MEMORY ps_malloc( &mdm       , capacity1 );
                MEMORY ps_malloc( &capacity2 , capacity1 );
                MEMORY ps_malloc( &extension2, capacity1 );
                MEMORY ps_malloc( &capacity3 , capacity1 );
                MEMORY ps_malloc( &extension3, capacity1 );
                MEMORY ps_malloc( &capacity4 , capacity1 );
                MEMORY ps_malloc( &extension4, capacity1 );
                std::fill( mdm       , mdm        + capacity1, (T***)NULL );
                std::fill( capacity2 , capacity2  + capacity1,      0     );
                std::fill( extension2, extension2 + capacity1,      0     );
                std::fill( capacity3 , capacity3  + capacity1, (T*  )NULL );
                std::fill( extension3, extension3 + capacity1, (T*  )NULL );
                std::fill( capacity4 , capacity4  + capacity1, (T** )NULL );
                std::fill( extension4, extension4 + capacity1, (T** )NULL );
            }
            else
            {
                PS_SIT oldsizearr = capacity1;
                capacity1 = (path[0]+1)*2;
                OPERA resize( &mdm       , oldsizearr, capacity1 );
                OPERA resize( &capacity2 , oldsizearr, capacity1 );
                OPERA resize( &extension2, oldsizearr, capacity1 );
                OPERA resize( &capacity3 , oldsizearr, capacity1 );
                OPERA resize( &extension3, oldsizearr, capacity1 );
                OPERA resize( &capacity4 , oldsizearr, capacity1 );
                OPERA resize( &extension4, oldsizearr, capacity1 );
                std::fill( mdm        + oldsizearr, mdm        + capacity1, (T***)NULL );
                std::fill( capacity2  + oldsizearr, capacity2  + capacity1,      0     );
                std::fill( extension2 + oldsizearr, extension2 + capacity1,      0     );
                std::fill( capacity3  + oldsizearr, capacity3  + capacity1, (T*  )NULL );
                std::fill( extension3 + oldsizearr, extension3 + capacity1, (T*  )NULL );
                std::fill( capacity4  + oldsizearr, capacity4  + capacity1, (T** )NULL );
                std::fill( extension4 + oldsizearr, extension4 + capacity1, (T** )NULL );
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
            capacity2[path[0]] = (path[1]+1)*2;
            MEMORY ps_malloc( &(mdm[path[0]])       , capacity2[path[0]] );
            MEMORY ps_malloc( &(capacity3[path[0]]) , capacity2[path[0]] );
            MEMORY ps_malloc( &(capacity4[path[0]]) , capacity2[path[0]] );
            MEMORY ps_malloc( &(extension3[path[0]]), capacity2[path[0]] );
            MEMORY ps_malloc( &(extension4[path[0]]), capacity2[path[0]] );
            std::fill( mdm[path[0]]       , mdm[path[0]]        + capacity2[path[0]], (T**)NULL );
            std::fill( capacity3[path[0]] , capacity3[path[0]]  + capacity2[path[0]], (T  )NULL );
            std::fill( capacity4[path[0]] , capacity4[path[0]]  + capacity2[path[0]], (T* )NULL );
            std::fill( extension3[path[0]], extension3[path[0]] + capacity2[path[0]], (T  )NULL );
            std::fill( extension4[path[0]], extension4[path[0]] + capacity2[path[0]], (T* )NULL );
        }
    }
    // END FIRST DIMENSION


    // START SECOND DIMENSION
    {
        // check for overflow of second dimension
        if( path[1] >= capacity2[path[0]] )
        {
            if( !capacity2[path[0]] )
            {
                capacity2[path[0]] = (path[1]+1)*2;
                MEMORY ps_malloc( &(mdm[path[0]])       , capacity2[path[0]] );
                MEMORY ps_malloc( &(capacity3[path[0]]) , capacity2[path[0]] );
                MEMORY ps_malloc( &(extension3[path[0]]), capacity2[path[0]] );
                MEMORY ps_malloc( &(capacity4[path[0]]) , capacity2[path[0]] );
                MEMORY ps_malloc( &(extension4[path[0]]), capacity2[path[0]] );
                std::fill( mdm[path[0]]       , mdm[path[0]]        + capacity2[path[0]], (T**)NULL );
                std::fill( capacity3[path[0]] , capacity3[path[0]]  + capacity2[path[0]], (T  )NULL );
                std::fill( extension3[path[0]], extension3[path[0]] + capacity2[path[0]], (T  )NULL );
                std::fill( capacity4[path[0]] , capacity4[path[0]]  + capacity2[path[0]], (T* )NULL );
                std::fill( extension4[path[0]], extension4[path[0]] + capacity2[path[0]], (T* )NULL );
            }
            else
            {
                PS_SIT oldsizearr = capacity2[path[0]];
                capacity2[path[0]] = (path[1]+1)*2;
                OPERA resize( &(mdm[path[0]])       , oldsizearr, capacity2[path[0]] );
                OPERA resize( &(capacity3[path[0]]) , oldsizearr, capacity2[path[0]] );
                OPERA resize( &(extension3[path[0]]), oldsizearr, capacity2[path[0]] );
                OPERA resize( &(capacity4[path[0]]) , oldsizearr, capacity2[path[0]] );
                OPERA resize( &(extension4[path[0]]), oldsizearr, capacity2[path[0]] );
                std::fill( mdm[path[0]]        + oldsizearr,
                           mdm[path[0]]        + capacity2[path[0]], (T**)NULL );
                std::fill( capacity3[path[0]]  + oldsizearr,
                           capacity3[path[0]]  + capacity2[path[0]], (T  )NULL );
                std::fill( extension3[path[0]] + oldsizearr,
                           extension3[path[0]] + capacity2[path[0]], (T  )NULL );
                std::fill( capacity4[path[0]]  + oldsizearr,
                           capacity4[path[0]]  + capacity2[path[0]], (T* )NULL );
                std::fill( extension4[path[0]] + oldsizearr,
                           extension4[path[0]] + capacity2[path[0]], (T* )NULL );
            }

            extension2[path[0]] = path[1]+1;
        }
        else if( path[1] >= extension2[path[0]] )
        {
            extension2[path[0]] = path[1]+1;
        }

        // check if anything's been set in this index yet
        if( !mdm[path[0]][path[1]] )
        {
            capacity3[path[0]][path[1]] = (path[2]+1)*2;
            MEMORY ps_malloc( &(mdm[path[0]][path[1]])       , capacity3[path[0]][path[1]] );
            MEMORY ps_malloc( &(capacity4[path[0]][path[1]]) , capacity3[path[0]][path[1]] );
            MEMORY ps_malloc( &(extension4[path[0]][path[1]]), capacity3[path[0]][path[1]] );
            std::fill( mdm[path[0]][path[1]],
                       mdm[path[0]][path[1]]        + capacity3[path[0]][path[1]], (T*)NULL );
            std::fill( capacity4[path[0]][path[1]],
                       capacity4[path[0]][path[1]]  + capacity3[path[0]][path[1]], (T )NULL );
            std::fill( extension4[path[0]][path[1]],
                       extension4[path[0]][path[1]] + capacity3[path[0]][path[1]], (T )NULL );
        }
    }
    // END SECOND DIMENSION


    // START THIRD DIMENSION
    {
        // check for overflow of second dimension
        if( path[2] >= capacity3[path[0]][path[1]] )
        {
            if( !capacity3[path[0]][path[1]] )
            {
                capacity3[path[0]][path[1]] = (path[2]+1)*2;
                MEMORY ps_malloc( &(mdm[path[0]][path[1]])       , capacity3[path[0]][path[1]] );
                MEMORY ps_malloc( &(capacity4[path[0]][path[1]]) , capacity3[path[0]][path[1]] );
                MEMORY ps_malloc( &(extension4[path[0]][path[1]]), capacity3[path[0]][path[1]] );
                std::fill( mdm[path[0]][path[1]]       ,
                           mdm[path[0]][path[1]]        + capacity3[path[0]][path[1]], (T*)NULL );
                std::fill( capacity4[path[0]][path[1]] ,
                           capacity4[path[0]][path[1]]  + capacity3[path[0]][path[1]],      0    );
                std::fill( extension4[path[0]][path[1]],
                           extension4[path[0]][path[1]] + capacity3[path[0]][path[1]],      0    );
            }
            else
            {
                PS_SIT oldsizearr = capacity3[path[0]][path[1]];
                capacity3[path[0]][path[1]] = (path[2]+1)*2;
                OPERA resize( &(mdm[path[0]][path[1]])       , oldsizearr, capacity3[path[0]][path[1]] );
                OPERA resize( &(capacity4[path[0]][path[1]]) , oldsizearr, capacity3[path[0]][path[1]] );
                OPERA resize( &(extension4[path[0]][path[1]]), oldsizearr, capacity3[path[0]][path[1]] );
                std::fill( mdm[path[0]][path[1]]        + oldsizearr,
                           mdm[path[0]][path[1]]        + capacity3[path[0]][path[1]], (T*)NULL );
                std::fill( capacity4[path[0]][path[1]]  + oldsizearr,
                           capacity4[path[0]][path[1]]  + capacity3[path[0]][path[1]],     0     );
                std::fill( extension4[path[0]][path[1]] + oldsizearr,
                           extension4[path[0]][path[1]] + capacity3[path[0]][path[1]],     0     );
            }

            extension3[path[0]][path[1]] = path[2]+1;
        }
        else if( path[2] >= extension3[path[0]][path[1]] )
        {
            extension3[path[0]][path[1]] = path[2]+1;
        }

        // check if anything's been set in this index yet
        if( !mdm[path[0]][path[1]][path[2]] )
        {
            capacity4[path[0]][path[1]][path[2]] = (path[3]+1)*2;
            MEMORY ps_malloc( &(mdm[path[0]][path[1]][path[2]]), capacity4[path[0]][path[1]][path[2]] );
            std::fill( mdm[path[0]][path[1]][path[2]],
                       mdm[path[0]][path[1]][path[2]] + capacity4[path[0]][path[1]][path[2]], 0 );
        }
    }
    // END THIRD DIMENSION


    // check for overflow of fourth dimension
    if( path[3] >= capacity4[path[0]][path[1]][path[2]] )
    {
        PS_SIT oldsizearr = capacity4[path[0]][path[1]][path[2]];
        capacity4[path[0]][path[1]][path[2]] = (path[3]+1)*2;
        extension4[path[0]][path[1]][path[2]] = path[3]+1;
        OPERA resize( &(mdm[path[0]][path[1]][path[2]]), oldsizearr, capacity4[path[0]][path[1]][path[2]] );
        std::fill( mdm[path[0]][path[1]][path[2]] + oldsizearr,
                   mdm[path[0]][path[1]][path[2]] + capacity4[path[0]][path[1]][path[2]], 0 );
    }
    else if( path[3] >= extension4[path[0]][path[1]][path[2]] )
    {
        extension4[path[0]][path[1]][path[2]] = path[3]+1;
    }

    mdm[path[0]][path[1]][path[2]][path[3]] = val;
}

template <class T>
void pixsrc_d4matrix<T>::set( PS_SIT w, PS_SIT x, PS_SIT y, PS_SIT z, T val )
{
    PS_SIT path[4] = { w, x, y, z };
    set( path, val );
}

template <class T>
T pixsrc_d4matrix<T>::get( PS_SIT w, PS_SIT x, PS_SIT y, PS_SIT z )
{
    return mdm[w][x][y][z];
}

template <class T>
T pixsrc_d4matrix<T>::get( PS_SIT *path )
{
    return mdm[path[0]][path[1]][path[2]][path[3]];
}

template <class T>
PS_SIT pixsrc_d4matrix<T>::size( PS_SIT *path, PS_SIT dim )
{
    // dim==1 requests 1st dimensions size
    switch( dim )
    {
    case 0:
        return extension1;
        break;
    case 1:
        return extension2[ path[0] ];
        break;
    case 2:
        return extension3[ path[0] ][ path[1] ];
        break;
    case 3:
        return extension4[ path[0] ][ path[1] ][ path[2] ];
        break;
    default:
        return -1;
        break;
    }
}

template <class T>
PS_SIT pixsrc_d4matrix<T>::size()
{
    return extension1;
}

template <class T>
PS_SIT pixsrc_d4matrix<T>::size( PS_SIT x )
{
    return extension2[x];
}

template <class T>
PS_SIT pixsrc_d4matrix<T>::size( PS_SIT x, PS_SIT y )
{
    return extension3[x][y];
}

template <class T>
PS_SIT pixsrc_d4matrix<T>::size( PS_SIT x, PS_SIT y, PS_SIT z )
{
    return extension4[x][y][z];
}
