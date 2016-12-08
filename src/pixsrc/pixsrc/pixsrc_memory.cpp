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




#include "pixsrc_memory_templates.cpp"
#include "pixsrc_constants.hpp"
#include "pixsrc_operations.hpp"
#include "pixsrc_tps.hpp"

#ifdef PIXSRC_MEMORYDEBUG
PS_unsignedSIT    pixsrc_memory::nummallocs;
PS_unsignedSIT    pixsrc_memory::numfrees;
pthread_mutex_t pixsrc_memory::mallocmutex;
pthread_mutex_t pixsrc_memory::freemutex;
#endif

void pixsrc_memory::free_matrices(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    // Note: residual4stat is not destroyed/freed because it only points to either
    // residual or uv_residual
    // same goes for data4stat and lo4stat

    if( vars_->data                )
        vars_->data->                ~VECTOR();
    if( vars_->mps                 )
        vars_->mps->                 ~VECTOR();
    if( vars_->lod                 )
        vars_->lod->                 ~VECTOR();
    if( vars_->lensedmps           )
        vars_->lensedmps->           ~VECTOR();
    if( vars_->residual            )
        vars_->residual->            ~VECTOR();
    if( vars_->uv_residual         )
        vars_->uv_residual->         ~VECTOR();
    if( vars_->interperr           )
        vars_->interperr->           ~VECTOR();
    if( vars_->ieanalytic          )
        vars_->ieanalytic->          ~VECTOR();
    if( vars_->psf_analytic        )
        vars_->psf_analytic->        ~VECTOR();
    if( vars_->lensedmpsnobo       )
        vars_->lensedmpsnobo->       ~VECTOR();
    if( vars_->srcexactvec         )
        vars_->srcexactvec->         ~VECTOR();
    if( vars_->ncv                 )
        vars_->ncv->                 ~VECTOR();
    if( vars_->lensingoperatornobo )
        vars_->lensingoperatornobo-> ~MATRIX();
    if( vars_->blurringoperator    )
        vars_->blurringoperator->    ~MATRIX();
    if( vars_->b1                  )
        vars_->b1->                  ~MATRIX();
    if( vars_->h1                  )
        vars_->h1->                  ~MATRIX();
    if( vars_->h1x                 )
        vars_->h1x->                 ~MATRIX();
    if( vars_->h1y                 )
        vars_->h1y->                 ~MATRIX();
    if( vars_->c1                  )
        vars_->c1->                  ~MATRIX();
    if( vars_->a1                  )
        vars_->a1->                  ~MATRIX();
    if (vars_->tpsfit)
        ((TPS*)(vars_->tpsfit))->~TPS();

    if( !data_->nopsf && vars_->lensingoperator )
        vars_->lensingoperator->~MATRIX();

    MEMORY ps_free( vars_->dataptr                );
    MEMORY ps_free( vars_->mpsptr                 );
    MEMORY ps_free( vars_->lodptr                 );
    MEMORY ps_free( vars_->lensedmpsptr           );
    MEMORY ps_free( vars_->residualptr            );
    MEMORY ps_free( vars_->uv_residualptr         );
    MEMORY ps_free( vars_->interperrptr           );
    MEMORY ps_free( vars_->ieanalyticptr          );
    MEMORY ps_free( vars_->psf_analyticptr        );
    MEMORY ps_free( vars_->lensedmpsnoboptr       );
    MEMORY ps_free( vars_->srcexactvecptr         );
    MEMORY ps_free( vars_->ncvptr                 );
    MEMORY ps_free( vars_->lensingoperatornoboptr );
    MEMORY ps_free( vars_->blurringoperatorptr    );
    MEMORY ps_free( vars_->b1ptr                  );
    MEMORY ps_free( vars_->h1ptr                  );
    MEMORY ps_free( vars_->h1xptr                 );
    MEMORY ps_free( vars_->h1yptr                 );
    MEMORY ps_free( vars_->c1ptr                  );
    MEMORY ps_free( vars_->a1ptr                  );
    MEMORY ps_free ((TPS*)vars_->tpsfitptr);

    if( !data_->nopsf )
        MEMORY ps_free( vars_->lensingoperator );
}

void pixsrc_memory::freememory(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{

    if(data_->subsampling>1 && !vars_->terminatelensing)
    {
        for(PS_SIT r=0; r<data_->ndp; r++)
        {
            if( vars_->sspointer[r] )
            {
                MEMORY ps_free (vars_-> newloc_sslo[r]);
                MEMORY ps_free (vars_->sspointer[r]);
            }
        }
    }

    if( data_->interperr > 1 && vars_->newloc_ssas )
    {
        MEMORY ps_free( vars_->newloc_ssas, data_->ndp );
    }

    MEMORY ps_free( vars_-> newloc_sslo  );
    MEMORY ps_free( vars_->sspointer     );
    MEMORY ps_free( vars_->penalties     );
    MEMORY ps_free( vars_->convexhull    );
    MEMORY ps_free( vars_->pointer       );
    MEMORY ps_free( vars_->newloc        );
    MEMORY ps_free( vars_->r4r           );
    MEMORY ps_free( vars_->r4rback       );
    MEMORY ps_free( vars_->fallsinsrc    );
    MEMORY ps_free( vars_->srcloc        );
    MEMORY ps_free( vars_->c4c           );
    MEMORY ps_free( vars_->c4cback       );
    MEMORY ps_free( vars_->srclocall     );
    MEMORY ps_free( vars_->magnification );
    MEMORY ps_free( vars_->dav           );

    MEMORY triangulatestructdestruct( vars_->triin  );
    MEMORY triangulatestructdestruct( vars_->triout );
    MEMORY ps_free(                   vars_->triin  );
    MEMORY ps_free(                   vars_->triout );

    MEMORY free_matrices( data_, cdata_, vars_ );

    if( data_->gridtype == 2 )
    {
        delete vars_->grid;
        delete vars_->gridpointer;
    }

    MEMORY ps_free( vars_->bestanalyticsrc );

    MEMORY ps_free( vars_ );

#ifdef PIXSRC_MEMORYDEBUG
    PRINTER print2screen(data_->print2screenname, "number of mallocs: " +
                         OPERA tostring(nummallocs), cdata_->print2screenmutex);
    PRINTER print2screen(data_->print2screenname, "number of frees  : " +
                         OPERA tostring( numfrees ), cdata_->print2screenmutex);
#endif

}

void pixsrc_memory::initialize(inputdata *data_, commoninputdata *cdata_, lensvar *vars_)
{
    MEMORY ps_malloc( &(vars_->penalties), data_->extlengths[5] );

    vars_->data                   = NULL;
    vars_->mps                    = NULL;
    vars_->lod                    = NULL;
    vars_->lensedmps              = NULL;
    vars_->residual               = NULL;
    vars_->uv_residual            = NULL;
    vars_->interperr              = NULL;
    vars_->ieanalytic             = NULL;
    vars_->psf_analytic           = NULL;
    vars_->lensedmpsnobo          = NULL;
    vars_->srcexactvec            = NULL;
    vars_->ncv                    = NULL;
    vars_->lensingoperator        = NULL;
    vars_->lensingoperatornobo    = NULL;
    vars_->blurringoperator       = NULL;
    vars_->b1                     = NULL;
    vars_->h1                     = NULL;
    vars_->h1x                    = NULL;
    vars_->h1y                    = NULL;
    vars_->c1                     = NULL;
    vars_->a1                     = NULL;
    vars_->tpsfit                 = NULL;
    vars_->dataptr                = NULL;
    vars_->mpsptr                 = NULL;
    vars_->lodptr                 = NULL;
    vars_->lensedmpsptr           = NULL;
    vars_->residualptr            = NULL;
    vars_->uv_residualptr         = NULL;
    vars_->interperrptr           = NULL;
    vars_->ieanalyticptr          = NULL;
    vars_->psf_analyticptr        = NULL;
    vars_->lensedmpsnoboptr       = NULL;
    vars_->srcexactvecptr         = NULL;
    vars_->ncvptr                 = NULL;
    vars_->lensingoperatorptr     = NULL;
    vars_->lensingoperatornoboptr = NULL;
    vars_->blurringoperatorptr    = NULL;
    vars_->b1ptr                  = NULL;
    vars_->h1ptr                  = NULL;
    vars_->h1xptr                 = NULL;
    vars_->h1yptr                 = NULL;
    vars_->c1ptr                  = NULL;
    vars_->a1ptr                  = NULL;
    vars_->tpsfitptr              = NULL;

    vars_->dav                    = NULL;

    vars_->newloc_sslo            = NULL;
    vars_->newloc_ssas            = NULL;
    vars_->sspointer              = NULL;
    vars_->pointer                = NULL;
    vars_->convexhull             = NULL;
    vars_->fallsinsrc             = NULL;
    vars_->srcloc                 = NULL;
    vars_->c4c                    = NULL;
    vars_->c4cback                = NULL;
    vars_->srclocall              = NULL;
    vars_->magnification          = NULL;
    vars_->bestanalyticsrc        = NULL;


    MEMORY ps_malloc(             &vars_->triin , 1 );
    MEMORY ps_malloc(             &vars_->triout, 1 );
    MEMORY triangulatestructinit(  vars_->triin     );
    MEMORY triangulatestructinit(  vars_->triout    );
}

void pixsrc_memory::triangulatestructinit(struct triangulateio *tristruct)
{
    tristruct->pointlist             = NULL;
    tristruct->pointattributelist    = NULL;
    tristruct->pointmarkerlist       = NULL;
    tristruct->trianglelist          = NULL;
    tristruct->triangleattributelist = NULL;
    tristruct->trianglearealist      = NULL;
    tristruct->neighborlist          = NULL;
    tristruct->segmentlist           = NULL;
    tristruct->segmentmarkerlist     = NULL;
    tristruct->holelist              = NULL;
    tristruct->regionlist            = NULL;
    tristruct->edgelist              = NULL;
    tristruct->edgemarkerlist        = NULL;
    tristruct->normlist              = NULL;

    tristruct->numberofpoints             = 0;
    tristruct->numberofpointattributes    = 0;
    tristruct->numberoftriangles          = 0;
    tristruct->numberofcorners            = 0;
    tristruct->numberoftriangleattributes = 0;
    tristruct->numberofsegments           = 0;
    tristruct->numberofholes              = 0;
    tristruct->numberofregions            = 0;
    tristruct->numberofedges              = 0;
}

void pixsrc_memory::triangulatestructdestruct(struct triangulateio *tristruct)
{
    if(!tristruct)
        return;

    MEMORY ps_free( tristruct->pointlist             );
    MEMORY ps_free( tristruct->pointattributelist    );
    MEMORY ps_free( tristruct->pointmarkerlist       );
    MEMORY ps_free( tristruct->trianglelist          );
    MEMORY ps_free( tristruct->triangleattributelist );
    MEMORY ps_free( tristruct->trianglearealist      );
    MEMORY ps_free( tristruct->neighborlist          );
    MEMORY ps_free( tristruct->segmentlist           );
    MEMORY ps_free( tristruct->segmentmarkerlist     );
    MEMORY ps_free( tristruct->holelist              );
    MEMORY ps_free( tristruct->regionlist            );
    MEMORY ps_free( tristruct->edgelist              );
    MEMORY ps_free( tristruct->edgemarkerlist        );
    MEMORY ps_free( tristruct->normlist              );
}

void pixsrc_memory::check_ptr (bool *ptr, size_t size) {ptr[size-1]=0;}

void pixsrc_memory::check_ptr (char *ptr, size_t size) {ptr[size-1]=0;}

void pixsrc_memory::check_ptr (signed char *ptr,          size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (signed short int *ptr,     size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (signed int *ptr,           size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (signed long int *ptr,      size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (signed long long int *ptr, size_t size) {ptr[size-1]=0;}

void pixsrc_memory::check_ptr (unsigned char *ptr,          size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (unsigned short int *ptr,     size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (unsigned int *ptr,           size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (unsigned long int *ptr,      size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (unsigned long long int *ptr, size_t size) {ptr[size-1]=0;}

void pixsrc_memory::check_ptr (float *ptr,       size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (double *ptr,      size_t size) {ptr[size-1]=0;}
void pixsrc_memory::check_ptr (long double *ptr, size_t size) {ptr[size-1]=0;}

void pixsrc_memory::check_ptr (void *ptr, size_t size) {}

