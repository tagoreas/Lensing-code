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



#include "pixsrc_cartesian.hpp"
#include "pixsrc_adaptive.hpp"
#include "pixsrc_irrcart.hpp"
#include "pixsrc_shapelets.hpp"
#include "pixsrc_init.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_operations_templates.cpp"
#include "pixsrc_printer_templates.cpp"
#include "pixsrc_nonparamlens.hpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_help.hpp"
#include <unistd.h> // for sleep
#include <string.h>
#include <iostream>
#include <fstream>


//
// A lot of functions called from lensmodel are defined here.
//

double get_evidence (inputdata*, commoninputdata*);

// struct holding particular parameters for a pixsrc thread
typedef struct
{
    PS_SIT u;
    inputdata *data_;
    commoninputdata *cdata_;
    PS_SIT firstfree;
    bool *threads;
    PS_SIT tracker;
} args2launch;

// multithreading code, launching lensing routines
void* launchthat (void *args_)
{
    args2launch *args = (args2launch*)args_;

    if (args->data_[args->u].use_shapelets)
        pixsrc_shapelets ps(args->u, args->data_, args->cdata_, args->tracker);
    else if (args->data_[args->u].gridtype==0)
        pixsrc_cartesian ps(args->u, args->data_, args->cdata_, args->tracker);
    else if (args->data_[args->u].gridtype==1)
        pixsrc_adaptive  ps(args->u, args->data_, args->cdata_, args->tracker);
    else if (args->data_[args->u].gridtype==2)
        pixsrc_irrcart   ps(args->u, args->data_, args->cdata_, args->tracker);

    args->threads[args->firstfree] = false;

    return NULL;
}

extern "C" void ps_getpixsrcinfo (char *basename, size_t *memorysizedata, size_t *memorysizecdata, size_t *numimages_)
{
    PS_SIT numimages[1] = {0};
    char *file;
    const char *listcc[3];
    listcc[0] = CONSTANT dir_in;
    listcc[1] = basename;
    listcc[2] = CONSTANT imagelistfile;
    OPERA concatenate( listcc, 3, &file );

    char **list;
    PS_SIT listsize;
    OPERA readfile( file, &list, &listsize, NULL);
    MEMORY ps_free( file );
    *numimages = listsize;

    // this removes any "-ignore" entries
    char **entries;
    PS_SIT numentries;
    for(PS_SIT j=0; j<listsize; ++j)
    {
        OPERA split( list[j], " ", &entries, &numentries );

        if( !strcmp(entries[0], "-ignore") )
        {
            --*numimages;
        }

        MEMORY ps_free( entries, numentries );
    }

    MEMORY ps_free( list, listsize );

    // if no image file was present, assume one image exists
    *numimages_      = *numimages<=0 ? (size_t)1 : (size_t)(*numimages);

    *memorysizecdata =             sizeof(commoninputdata);
    *memorysizedata  = *numimages_*sizeof(      inputdata);
}

#ifdef PS_HAVE_GRAVLENS
extern "C" void ps_loaddata (char *basename, size_t numimages_, void *data__, void *cdata__)
#endif
#ifdef PS_HAVE_TRIAXIAL
    extern "C" void ps_loaddata (char *basename, size_t numimages_, void *data__, void *cdata__, double **tlmparms, time_t **tlmtime, double **tlmenvirogals)
#endif
{
    PS_SIT numimages = (PS_SIT)numimages_;
    inputdata *data_ = (inputdata*)data__;
    commoninputdata *cdata_ = (commoninputdata*)cdata__;
    char *bn = basename;
    char **name;
    char **namewithext;
    PS_SIT dimil=0;
    char **ignorelist;
    char **ignorelistwithext;

    MEMORY ps_malloc (&name       , numimages);
    MEMORY ps_malloc (&namewithext, numimages);

    INIT settingup                   (bn, name, namewithext, data_, cdata_,
                                      &dimil, &ignorelist, &ignorelistwithext);
    INIT setdefaultparameters        (bn, name, namewithext, data_, cdata_);
    INIT readparameters              (bn, name, namewithext, data_, cdata_,
                                      dimil,  ignorelist,  ignorelistwithext);
    INIT readimages                  (bn, name, namewithext, data_, cdata_);
    INIT setuplibwcs                 (bn, name, namewithext, data_, cdata_);
    INIT readmasks                   (bn, name, namewithext, data_, cdata_);
    INIT setuppsfs                   (bn, name, namewithext, data_, cdata_);
    INIT uvdata                      (bn, name, namewithext, data_, cdata_);
    INIT finalizeNcounterindications (bn, name, namewithext, data_, cdata_);
    INIT printsettings               (                       data_, cdata_);

    MEMORY ps_free (name             , numimages);
    MEMORY ps_free (namewithext      , numimages);
    MEMORY ps_free (ignorelist       , dimil    );
    MEMORY ps_free (ignorelistwithext, dimil    );
    
#ifdef PS_HAVE_TRIAXIAL
    // for triaxial halo models
    cdata_->tlmparms = tlmparms;
    cdata_->tlmtime  = tlmtime;
    cdata_->tlmenvirogals  = tlmenvirogals;
#endif
    
}

extern "C" void ps_launch (void *data__, void *cdata__, double *lE)
{
    // This hack calls a bash script to compute a statistic
    // So pixsrc is not really being used here
    if (0)
    {
        double stat = 1e15;
        if (!system("bash get-statistic.sh"))
        {
            string line;
            std::ifstream myfile ("stat.dat");
            if (myfile.is_open())
            {
                getline (myfile,line);
                stat = OPERA convert_string <double> (line);
                stat = OPERA is_finite (stat) ? stat : 1e15;
                myfile.close();
            }
        }
        *lE = -0.5*stat;
        return;
    }

    inputdata        *data_ = (      inputdata*) data__;
    commoninputdata *cdata_ = (commoninputdata*)cdata__;

    // say hello to the user (print to screen)
    string p2sn(data_[0].print2screenname);
    for(PS_unsignedSIT h=0; h<p2sn.size(); h++)
        p2sn[h]=' ';
    for(PS_SIT g=0; g<cdata_->numimages; g++)
        if(data_[g].verbose>1)
        {
            pthread_mutex_lock (cdata_->print2screenmutex);
            std::cout << std::endl;
            pthread_mutex_unlock (cdata_->print2screenmutex);

            PRINTER print2screen (p2sn,
                                  "entering pixsrc!",
                                  cdata_->print2screenmutex);
            break;
        }

    // if using non-parametric lens potential perturbations
    if( cdata_->npl )
    {
        NONPARAMLENS init_npl    (data_, cdata_);
        NONPARAMLENS minimize    (data_, cdata_);
        NONPARAMLENS print_model (data_, cdata_);
    }
    else
    {
        *lE = get_evidence( data_, cdata_ );
    }

    // say goodbye (print to screen)
    for(PS_SIT g=0; g<cdata_->numimages; g++)
        if(data_[g].verbose>1)
        {
            PRINTER print2screen (p2sn, "leaving pixsrc!\n",
                                  cdata_->print2screenmutex);
            break;
        }
}

double get_evidence (inputdata *data_, commoninputdata *cdata_)
{
    bool threads[cdata_->numthreads];
    std::fill(threads,threads+cdata_->numthreads,false);
    args2launch argsarr[cdata_->numthreads];
    pthread_t thread[cdata_->numthreads];

    for(PS_SIT u = 0; u < cdata_->numimages; u++)
    {
        PS_SIT numinnerloop = 1;
        if(data_[u].magparams > 0)
            numinnerloop = data_[u].magparams;
        if(data_[u].traceparams[2] > 0)
            numinnerloop = (PS_SIT)(data_[u].traceparams[2])+1;

        for(PS_SIT tracker=0; tracker<numinnerloop; tracker++)
        {
            bool foundthread = false;

            while(!foundthread)
            {
                PS_SIT firstfree=-1;
                for(PS_SIT x=0; x<cdata_->numthreads; x++)
                    if(!threads[x])
                    {
                        firstfree=x;
                        break;
                    }
                if(firstfree!=-1)
                {
                    threads[firstfree]           = true;
                    argsarr[firstfree].u         = u;
                    argsarr[firstfree].data_     = data_;
                    argsarr[firstfree].cdata_    = cdata_;
                    argsarr[firstfree].firstfree = firstfree;
                    argsarr[firstfree].threads   = threads;
                    argsarr[firstfree].tracker   = tracker;

                    pthread_create( &thread[firstfree], cdata_->attrdetached,
                                    launchthat,  &argsarr[firstfree]          );
                    foundthread = true;

                    //if all pixsrc threads have been spawned
                    if(u == cdata_->numimages-1 && tracker == numinnerloop-1)
                        // wait for all threads to finish
                        for(PS_SIT i=0; i<cdata_->numthreads; i++)
                            while (threads[i])
                                usleep( CONSTANT shortwaitms );
                }
                else
                    usleep( CONSTANT longwaitms );
            }
        }
    }

    for(PS_SIT u = 0; u < cdata_->numimages; u++)
    {
        if(data_[u].magparams>0)
            PRINTER closeoutstream(data_[u].mags->stream);
        if(data_[u].traceparams[2]>0)
            PRINTER closeoutstream(data_[u].traces->stream);
        //if(data_[u].printdetails)
        //PRINTER closeoutstream(data_[u].details->stream);
    }

    double logevidence=0;
    for(PS_SIT u=0; u<cdata_->numimages; u++)
    {
        if(OPERA is_finite(data_[u].evidence))
            logevidence+=data_[u].evidence;
        else
        {
            OPERA assign_n_infinity( &logevidence );
            break;
        }
    }

    cdata_->numcall++;
    return logevidence;
}
extern "C" void ps_freepixsrcmem (void *data__, void *cdata__)
{
    // this would ideally free all pixsrc inputdata and commoninputdata memory
}

extern "C" void ps_help (int argc, char **argv)
{
    if (argc==2)
    {
        HELP general ();
    }
    else if (argc==3 && !strcmp(argv[2],"parms"))
    {
        HELP parms ();
    }
    else if (argc==3 && !strcmp(argv[2],"input"))
    {
        HELP input (1);
    }
    else if (argc==3 && !strcmp(argv[2],"output"))
    {
        HELP output (1);
    }
    else if (argc==3 && !strcmp(argv[2],"template"))
    {
        HELP ps_template ();
    }
    else if (argc==3 && !strcmp(argv[2],"manual"))
    {
        HELP manual ();
    }
    else if (argc==4 && !strcmp(argv[2],"parms") && !strcmp(argv[3],"more"))
    {
        HELP parms_more ();
    }
    else
    {
        std::cout <<
            "pixsrc :: Help option not recognized. Try" << std::endl <<
            "pixsrc help"                               << std::endl <<
            std::endl;
    }
}
