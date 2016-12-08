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



///
/// The __USE_PIXSRC_CUDA__ variable passed on to CXX
/// decides whether this source code will be compiled
///

#include "pixsrc_init.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_printer.hpp"
#include <cstring>

void pixsrc_init::setgpus( char *input, inputdata *data_, commoninputdata *cdata_ )
{
    if( !strcmp( input, "help") )
    {
        INIT printgpudevices( cdata_ );
        exit(0);
    }

    // check for ignore flag
    if( !strncmp( input, "ignore ", 7 ) )
    {
        char ** splitted;
        PS_SIT numsplit;

        OPERA split( input, " ", &splitted, &numsplit );

        for( PS_SIT j=1; j<numsplit; ++j )
        {
            for( PS_SIT dev=0; dev<cdata_->numgpudevices; ++dev )
            {
                if( !strcmp( splitted[j], cdata_->gpudevices[dev].name ) )
                {
                    cdata_->gpudevices[dev].ignoreme = 1;
                    break;
                }
                // if didn't find a match for input name
                else if( dev == cdata_->numgpudevices-1 )
                {
                    INIT printgpudevices( cdata_ );

                    PRINTER printerror("pixsrc", "invalid gpu ignore list",
                                       cdata_->print2screenmutex);
                }
            }
        }

        MEMORY ps_free( splitted, numsplit );

        return;
    }

    if( !strcmp( input, "benchmark") )
    {
        CUDA benchmark( cdata_ );
        exit(0);
    }

    PS_SIT numlist;
    char ** list;
    OPERA split( input, " ", &list, &numlist );

    data_->numgpu2use = numlist;
    MEMORY ps_malloc( &(data_->gpu2use), numlist );

    if( numlist != 1 )
    {
        PRINTER printerror( "pixsrc", "multiple gpu selection currently not supported",
                            cdata_->print2screenmutex);
    }

    char foundit;

    for( PS_SIT j=0; j<numlist; ++j )
    {
        for( PS_SIT g=0; g<cdata_->numgpudevices; ++g )
        {
            foundit = 0;
            // if name of gpu requested and name found match
            if( !strcmp( list[j], cdata_->gpudevices[g].name) )
            {
                data_->gpu2use[j] = g;
                foundit = 1;

                /*
                  string envstr = "CUDA_VISIBLE_DEVICES=" + OPERA tostring(g);
                  char envcstr[envstr.size()+1];
                  strcpy( envcstr, envstr.c_str() );
                  putenv( envcstr );
                */

                CUDA create_handle( cdata_, g );
            }
            else if( g == cdata_->numgpudevices-1 )
            {
                INIT printgpudevices( cdata_ );

                PRINTER printerror("pixsrc",
                                   "requested gpu device "+string(list[j])+" not found\n"
                                   "Note: if you've set more than 1 gpu, only the first will be visible to this application",
                                   cdata_->print2screenmutex);
            }

            if( foundit )
                break;
        }
    }

    MEMORY ps_free( list, numlist );
}

void pixsrc_init::printgpudevices( commoninputdata *cdata_ )
{
    // this function is usually followed by an exitting command

    CUDA detectgpus( cdata_ );

    PRINTER print2screen("pixsrc", OPERA tostring(cdata_->numgpudevices) + " gpu device(s) found:",
                         cdata_->print2screenmutex);

    for( PS_SIT j=0; j<cdata_->numgpudevices; ++j )
    {
        string id = OPERA tostring(cdata_->gpudevices[j].busid) + "." +
            OPERA tostring(cdata_->gpudevices[j].locid);
        PRINTER print2screen("pixsrc", string(cdata_->gpudevices[j].name),
                             cdata_->print2screenmutex);
        if( j == cdata_->numgpudevices-1 )
            id += "\n";
        PRINTER print2screen("pixsrc", "    on PCI Bus " + id,
                             cdata_->print2screenmutex);
    }
}
