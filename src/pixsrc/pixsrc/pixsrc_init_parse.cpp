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



//
// this file has methods to parse files
//

#include "pixsrc_init.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_wcs.hpp"
#include "pixsrc_printer.hpp"
#include <cstring>

void pixsrc_init::readmaskpolygon(PS_SIT *ext1, PS_SIT *ext2, double ***polys, char **file, inputdata *data_, commoninputdata *cdata_)
{
    *ext1 = *ext2 = 0;

    char **vec = 0;
    PS_SIT vecsize;
    OPERA readfile( *file, &vec, &vecsize, cdata_->print2screenmutex);

    if( vecsize )
    {
        *ext1 = 0;
        PS_SIT mode=-1;
        PS_SIT maxext2 = 0;

        char **list;
        PS_SIT numlist;

        PS_SIT currentext2;
        for( PS_SIT v=0; v<vecsize; ++v )
        {
            if( !strncmp(vec[v], "polygon", 7) )
            {
                (*ext1)++;

                OPERA split( vec[v], ",", &list, &numlist );

                currentext2 = numlist/2;
                if(currentext2>maxext2)
                    maxext2 = currentext2;

                MEMORY ps_free( list, numlist );
            }
            else if( !strcmp(vec[v], "image" ) )
                mode=0;
            else if( !strcmp(vec[v],  "fk5"  ) )
                mode=1;
        }

        if(mode==-1)
            PRINTER printerror("","invalid mask/reg file",
                               cdata_->print2screenmutex);

        *ext2 = maxext2*2;

        MEMORY ps_malloc( &(*polys), *ext1, *ext2 );

        PS_SIT index1 = 0;
        char *ptr1, *ptr2;
        pthread_mutex_lock( cdata_->wcsmutex );
        for( PS_SIT v=0; v<vecsize; ++v )
        {
            if( !strncmp(vec[v], "polygon", 7) )
            {
                // cut off "polygon(" from beginning and ")" from end

                // cutt of beginning
                ptr1 = vec[v];
                ptr2 = vec[v]+8;
                while( *ptr2 )
                {
                    *ptr1 = *ptr2;
                    ++ptr1;
                    ++ptr2;
                }
                // cut off end
                *ptr1 = 0;

                OPERA split( vec[v], ",", &list, &numlist );

                if(mode==0)
                {
                    for( PS_SIT q=0; q<numlist/2; ++q )
                    {
                        (*polys)[index1][q*2  ] = OPERA convert_string <double> (list[2*q])-1;
                        (*polys)[index1][q*2+1] = data_->imgy-OPERA convert_string <double> (list[2*q+1]);
                    }
                }
                else if(mode==1)
                {
                    for( PS_SIT q=0; q<numlist/2; ++q )
                    {
                        double ra =  OPERA convert_string <double> (list[2*q  ]);
                        double dec = OPERA convert_string <double> (list[2*q+1]);
                        HEADER getimgpixcoord( data_->wcs, data_->imgy, cdata_, ra,dec,
                                               &((*polys)[index1][q*2]),
                                               &((*polys)[index1][q*2+1]) );
                    }
                }

                MEMORY ps_free( list, numlist );

                // fill remaining points with last point read in
                for( PS_SIT s=numlist/2; s<*ext2/2; ++s )
                {
                    (*polys)[index1][s*2  ] = (*polys)[index1][(s-1)*2  ];
                    (*polys)[index1][s*2+1] = (*polys)[index1][(s-1)*2+1];
                }

                index1++;
            }
        }
        pthread_mutex_unlock( cdata_->wcsmutex );
    }

    MEMORY ps_free( vec, vecsize );
}

void pixsrc_init::readmaskcircle(PS_SIT *ext1, PS_SIT *ext2, double ***circles, char **file, inputdata *data_, commoninputdata *cdata_)
{
    *ext1 = *ext2 = 0;

    char **vec = 0;
    PS_SIT vecsize;
    OPERA readfile( *file, &vec, &vecsize, cdata_->print2screenmutex);

    if( vecsize )
    {
        *ext1 = 0;
        *ext2 = 3;
        PS_SIT mode=-1;

        for( PS_SIT v=0; v<vecsize; ++v )
        {
            if( !strncmp(vec[v], "circle", 6) )
                (*ext1)++;
            else if( !strcmp(vec[v], "image") )
                mode=0;
            else if( !strcmp(vec[v], "fk5") )
                mode=1;
        }

        if( mode==-1 || !(*ext1) )
            PRINTER printerror("","invalid mask/reg file", cdata_->print2screenmutex);

        MEMORY ps_malloc( &(*circles), *ext1, *ext2 );

        PS_SIT index1 = 0;
        char *ptr1, *ptr2;
        char **list;
        PS_SIT numlist;
        pthread_mutex_lock( cdata_->wcsmutex );
        for( PS_SIT v=0; v<vecsize; ++v )
        {
            if( !strncmp(vec[v], "circle", 6) )
            {
                // cut off "circle(" from beginning and ")" from end

                // cut off beginning
                ptr1 = vec[v];
                ptr2 = vec[v]+7;
                while( *ptr2 )
                {
                    *ptr1 = *ptr2;
                    ++ptr1;
                    ++ptr2;
                }
                // cut off end
                *(ptr1-1) = 0;

                OPERA split( vec[v], ",", &list, &numlist );

                if(mode==0)
                {
                    (*circles)[index1][0] = OPERA convert_string <double> (list[0])-1;
                    (*circles)[index1][1] = data_->imgy-OPERA convert_string <double> (list[1]);
                    (*circles)[index1][2] = OPERA convert_string <double> (list[2]);
                    (*circles)[index1][2]*= (*circles)[index1][2];
                }
                else if(mode==1)
                {
                    double ra =  OPERA convert_string <double> (list[0]);
                    double dec = OPERA convert_string <double> (list[1]);
                    HEADER getimgpixcoord (data_->wcs, data_->imgy, cdata_, ra,dec,
					  &(*circles)[index1][0],&(*circles)[index1][1]);
                    double degree;
                    PS_SIT size = OPERA sizestring( list[2] );

                    if(list[2][size-1]=='\"')
                    {
                        list[2][size-1] = 0;
                        degree = OPERA convert_string <double> (list[2])/3600.0;
                    }
                    else
                    {
                        degree = OPERA convert_string <double> (list[2]);
                    }
                    (*circles)[index1][2] = degree*degree*
                        ( data_->invwcsinfo[0]*data_->invwcsinfo[0] +
                          data_->invwcsinfo[1]*data_->invwcsinfo[1] );

                }
                index1++;

                MEMORY ps_free( list, numlist );
            }
        }
        pthread_mutex_unlock( cdata_->wcsmutex );
    }
    if( !*ext1 )
    {
        *ext1=1;
        *ext2=3;
        MEMORY ps_malloc( circles, *ext1, *ext2 );
        OPERA assign_p_infinity( &(*circles)[0][2] );
    }

    MEMORY ps_free( vec, vecsize );
}

void pixsrc_init::readmaskascii(PS_SIT *dim, PS_SIT **list, char **file, inputdata *data_, commoninputdata *cdata_)
{
    char **vec = 0;
    OPERA readfile( *file, &vec, &(*dim), cdata_->print2screenmutex);

    if( *dim )
    {
        MEMORY ps_malloc( &(*list), *dim*2 );

        char **xy;
        PS_SIT xysize;

        for(PS_SIT j=0; j<*dim; j++)
        {
            OPERA split( vec[j], " ", &xy, &xysize );
            (*list)[2*j]   = OPERA convert_string <PS_SIT> (xy[0])-1;
            (*list)[2*j+1] = data_->imgy-OPERA convert_string <PS_SIT> (xy[1]);
            MEMORY ps_free( xy, xysize );
        }
    }

    MEMORY ps_free( vec, *dim );
}

PS_SIT pixsrc_init::readarr(char *in, const char *pattern0, const char *err,
                            double **result, commoninputdata *cdata_)
{
    //
    // this will return an array matching the pattern supplied.
    // it is okay if there are less arguments that in the pattern.
    // it is NOT okay if there are too many arguments.
    // if the pattern contains a string, then nothing is parsed for that argument
    // However, if the pattern specifies a string, but the value
    // is a number, it will be parsed anyway

    OPERA trim( in, NULL );

    char **entries,  **pattern;
    PS_SIT numentries, numpattern;

    OPERA split( in,       " ", &entries, &numentries);
    OPERA split( pattern0, " ", &pattern, &numpattern);

    // presence of 'Z' in string means 'in' can have any length
    if(numpattern<numentries && !strchr(pattern0,'Z'))
        PRINTER printerror("","invalid "+string(err)+": too many arguments supplied!",
                           cdata_->print2screenmutex);
    else if(!numentries)
        PRINTER printerror("","invalid "+string(err)+": no arguments supplied!",
                           cdata_->print2screenmutex);

    MEMORY ps_free  (   *result              );
    MEMORY ps_malloc( &(*result), numentries );

    for(PS_SIT j=0; j<numentries && j<numpattern; j++)
    {
        if( ( pattern[j][0]=='i' && OPERA isanumberint   (entries[j]) ) ||
            ( pattern[j][0]=='d' && OPERA isanumberfloat(entries[j]) ) )
        {
            (*result)[j] = OPERA convert_string <double> (entries[j]);
        }
        else
        {
            if( pattern[j][0]=='s' )
            {
                if (OPERA isanumberfloat (entries[j]))
                    (*result)[j] = OPERA convert_string <double> (entries[j]);
            }
            else if( pattern[j][0]=='i' )
                PRINTER printerror("","invalid "+string(err)+": must be an integer!",
                                   cdata_->print2screenmutex);
            else if( pattern[j][0]=='d' )
                PRINTER printerror("","invalid "+string(err)+": must be an number!",
                                   cdata_->print2screenmutex);
        }
    }

    MEMORY ps_free( entries, numentries );
    MEMORY ps_free( pattern, numpattern );

    return numentries;
}
