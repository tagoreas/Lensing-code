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



#include "pixsrc_init.hpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_printer.hpp"
#include <cstring>

void pixsrc_init::setsource( char *argv, inputdata *data_, commoninputdata *cdata_)
{
    char *fname;
    const char *listcc[2];
    listcc[0] = CONSTANT dir_in;
    listcc[1] = argv;
    OPERA concatenate( listcc, 2, &fname );

    char **file;
    PS_SIT filesize;
    OPERA readfile( fname, &file, &filesize, cdata_->print2screenmutex);

    if( !filesize )
    {
        PRINTER printerror( "pixsrc",
                            "no data found in file: " + string(fname),
                            cdata_->print2screenmutex );

        MEMORY ps_free( fname );

        return;
    }

    if( filesize%2 == 1 )
    {
        PRINTER printerror( "pixsrc",
                            "format of setsource file incorrect: " + string(fname),
                            cdata_->print2screenmutex );
    }

    PS_SIT totalnsrc = filesize/2;

    char **l1,  **l2;
    PS_SIT numl1, numl2;

    char **srcdata;
    // I don't know why I made srcdata have 19 elements, but whatever the reason,
    // first two elements are unused.
    MEMORY ps_malloc( &srcdata, 19 );

    // check for number of vector images
    PS_SIT num_vector_src = 0;
    for( PS_SIT src=0; src<totalnsrc; ++src )
        if (!strncmp(file[src], "vector_FILENAME", 7))
            ++num_vector_src;
    data_->num_vec_src = num_vector_src;
    MEMORY ps_malloc (&(data_->vector_src_pos),    num_vector_src);
    MEMORY ps_malloc (&(data_->vector_src_flux),   num_vector_src);
    MEMORY ps_malloc (&(data_->vector_src_numpix), num_vector_src);
    num_vector_src = 0;

    // record data from input file
    for( PS_SIT src=0; src<totalnsrc; ++src )
    {
        // split up lines with source parameters and source vary flags
        OPERA split( file[src          ], " ", &l1, &numl1 );
        OPERA split( file[src+totalnsrc], " ", &l2, &numl2 );

        if( numl1!=9 || numl2!=8 )
        {
            PRINTER printerror( "pixsrc",
                                "format of set source file incorrect: " + string(fname),
                                cdata_->print2screenmutex );
        }

        // combine the parameters and flags into one array
        for( PS_SIT j=0; j<9; ++j )
        {
            MEMORY ps_malloc( &(srcdata[2+j]), OPERA sizestring(l1[j])+1 );
            strcpy( srcdata[2+j], l1[j] );
        }
        for( PS_SIT j=0; j<8; ++j )
        {
            MEMORY ps_malloc( &(srcdata[2+9+j]), OPERA sizestring(l2[j])+1 );
            strcpy( srcdata[2+9+j], l2[j] );
        }

        MEMORY ps_free( l1, numl1 );
        MEMORY ps_free( l2, numl2 );

        // check if previous sources have been set
        double numsrc = 1;
        if( data_->usersetsrc )
        {
            numsrc += data_->usersetsrc[0];
        }

        // length of array to hold all sources
        // 17 for 1 profile type, 8 params, and 8 vary_or_not flags
        // +1 for first element that holds number of sources
        PS_SIT length = 17*numsrc+1;

        double *temp = data_->usersetsrc;
        data_->usersetsrc = NULL;
        MEMORY ps_malloc( &(data_->usersetsrc), length );
        data_->usersetsrc[0] = numsrc;

        // copy all previously set sources into array
        double *startpos = data_->usersetsrc + 1;
        if( numsrc>1 )
        {
            std::copy( temp+1, temp+1+17*((PS_SIT)numsrc-1),  startpos );
            startpos += 17*((PS_SIT)numsrc-1);
        }
        MEMORY ps_free( temp );

        // flag source profile type
        if( !strcmp(srcdata[2],"sersic") )
        {
            startpos[0] = 0;
        }
        else if( !strcmp(srcdata[2],"none") )
        {
            startpos[0] = -1;
        }
        else if (!strncmp(srcdata[2], "vector_FILENAME", 7))
        {
            startpos[0] = 1001+num_vector_src;

            // read in vector image file
            char **vecvec, *vecfname;
            PS_SIT vecvecsize;
            const char *veclistcc[2];
            veclistcc[0] = CONSTANT dir_in;
            veclistcc[1] = srcdata[2]+7;
            OPERA concatenate (veclistcc, 2, &vecfname);
            OPERA readfile (vecfname, &vecvec, &vecvecsize, cdata_->print2screenmutex);

            // allocate memory for source
            data_->vector_src_numpix[num_vector_src] = vecvecsize;
            MEMORY ps_malloc (&(data_->vector_src_pos[num_vector_src]),  vecvecsize*2);
            MEMORY ps_malloc (&(data_->vector_src_flux[num_vector_src]), vecvecsize  );

            // parse file line by line
            char *tok;
            for (PS_SIT s=0; s<vecvecsize; ++s)
            {
                tok = strtok (vecvec[s], " \t");
                if (!tok)
                {
                    PRINTER printerror("pixsrc","invalid vector source: " + string(vecfname),
                                       cdata_->print2screenmutex);
                }
                data_->vector_src_pos[num_vector_src][s*2  ] = OPERA convert_string <double> (tok);
                tok = strtok (NULL, " \t");
                if (!tok)
                {
                    PRINTER printerror("pixsrc","invalid vector source: " + string(vecfname),
                                       cdata_->print2screenmutex);
                }
                data_->vector_src_pos[num_vector_src][s*2+1] = OPERA convert_string <double> (tok);

                tok = strtok (NULL, " \t");
                if (!tok)
                {
                    PRINTER printerror("pixsrc","invalid vector source: " + string(vecfname),
                                       cdata_->print2screenmutex);
                }
                data_->vector_src_flux[num_vector_src][s] = OPERA convert_string <double> (tok);
            }

            MEMORY ps_free( vecvec, vecvecsize );
            MEMORY ps_free( vecfname );

            ++num_vector_src;
        }
        else
            PRINTER printerror("pixsrc","invalid source arguments: invalid source profile",
                               cdata_->print2screenmutex);
        ++startpos;

        // copy parameters for this source
        for( PS_SIT j=0; j<16; ++j )
        {
            startpos[j] = OPERA convert_string <double> ( srcdata[3+j] );
            if( j>7 && startpos[j] )
                ++data_->srcnvary;
        }

        for( PS_SIT j=2; j<19; ++j )
            MEMORY ps_free( srcdata[j] );

    }

    MEMORY ps_free( file, filesize );
    MEMORY ps_free( srcdata );
    MEMORY ps_free( fname );
}

void pixsrc_init::setsourcebounds( char *argv, inputdata *data_, commoninputdata *cdata_)
{
    char *fname;
    const char *listcc[2];
    listcc[0] = CONSTANT dir_in;
    listcc[1] = argv;
    OPERA concatenate( listcc, 2, &fname );

    char **file;
    PS_SIT filesize;
    OPERA readfile( fname, &file, &filesize, cdata_->print2screenmutex);

    if( !filesize )
    {
        PRINTER printerror( "pixsrc",
                            "no data found in file " + string(fname),
                            cdata_->print2screenmutex );

        MEMORY ps_free( fname );

        return;
    }

    data_->numsrcbounds = filesize;

    MEMORY ps_malloc( &(data_->srcbounds), data_->numsrcbounds, data_->usersetsrc[0]*8+2 );

    // pj = index of source parameter (0..number_of_sources*8)
    // iii = source number (>=1)
    // jjj = source parameter (1..8)
    PS_SIT pj, iii, jjj;
    char **l;
    PS_SIT llength;

    for( PS_SIT line=0; line<filesize; ++line )
    {
        OPERA split( file[line], " ", &l, &llength );

        if( (llength-2) % 3  )
        {
            PRINTER printerror( "pixsrc",
                                "invalid source bounds arguments",
                                cdata_->print2screenmutex );
        }

        std::fill( data_->srcbounds[line],
                   data_->srcbounds[line] + (PS_SIT)data_->usersetsrc[0]*8, 0 );

        // set lower and upper bounds
        data_->srcbounds[line][(PS_SIT)data_->usersetsrc[0]*8  ] =
            OPERA convert_string <double> ( l[llength-2] );
        data_->srcbounds[line][(PS_SIT)data_->usersetsrc[0]*8+1] =
            OPERA convert_string <double> ( l[llength-1] );

        // read in weights for each source parameter
        for( PS_SIT ent=0; ent<llength-2; ent+=3 )
        {
            iii = OPERA convert_string <PS_SIT> (l[ent  ]);
            jjj = OPERA convert_string <PS_SIT> (l[ent+1]);
            pj = (iii-1)*8+(jjj-1);
            data_->srcbounds[line][pj] = OPERA convert_string <double> ( l[ent+2] );
        }

        MEMORY ps_free( l, llength );
    }

    MEMORY ps_free( file, filesize );
    MEMORY ps_free( fname );
}

void pixsrc_init::setsourcestepsizes( char *argv, inputdata *data_, commoninputdata *cdata_)
{
    char *fname;
    const char *listcc[2];
    listcc[0] = CONSTANT dir_in;
    listcc[1] = argv;
    OPERA concatenate( listcc, 2, &fname );

    char **file;
    PS_SIT filesize;
    OPERA readfile( fname, &file, &filesize, cdata_->print2screenmutex);

    if( !filesize )
    {
        PRINTER printerror( "pixsrc",
                            "no data found in file " + string(fname),
                            cdata_->print2screenmutex );

        MEMORY ps_free( fname );

        return;
    }
    if( filesize != (PS_SIT)data_->usersetsrc[0] )
    {
        PRINTER printerror( "pixsrc",
                            "invalid stepsize file format",
                            cdata_->print2screenmutex );
    }

    MEMORY ps_malloc( &(data_->srcstepsizes), data_->srcnvary );

    PS_SIT index = 0;
    PS_SIT numl0;
    char **l0;
    for( PS_SIT line=0; line<filesize; ++line )
    {
        OPERA split( file[line], " ", &l0, &numl0 );
        if( numl0 != 8 )
        {
            PRINTER printerror( "pixsrc",
                                "invalid stepsize file format",
                                cdata_->print2screenmutex );
        }
        for( PS_SIT ent=0; ent<8; ++ent )
        {
            // if this parameter is to be varied
            if( data_->usersetsrc[1+line*17+9+ent] )
            {
                data_->srcstepsizes[index] = OPERA convert_string <double> ( l0[ent] );
                ++index;
            }
        }
        MEMORY ps_free( l0, numl0 );
    }

    MEMORY ps_free( file, filesize );
    MEMORY ps_free( fname );
}
