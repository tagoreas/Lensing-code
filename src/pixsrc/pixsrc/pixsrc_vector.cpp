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



#include "pixsrc_vector.hpp"
#include "pixsrc_memory_templates.cpp"
#include "pixsrc_cuda.hpp"
#include "pixsrc_inputdata.hpp"
#include "pixsrc_printer.hpp"

#include <algorithm>

// constructor
pixsrc_vector::pixsrc_vector( commoninputdata *cdata__, inputdata *data__, PS_SIT size_)
{
    cdata_ = cdata__;
    data_  = data__;

    size = size_;
    scalar = (PS_FPT)1;

    vec = (PS_FPT*)0;
    vec_to_return = 0;

    are_we_using_cuda = ( data_ && data_->numgpu2use ) ? 1 : 0;

    pthread_mutex_init (&initmutex, NULL);

    if( !are_we_using_cuda )
    {
        MEMORY ps_malloc( &vec, size_ );
        std::fill( vec, vec + size_, (PS_FPT)0 );
        up2date = pvs_cpu;
    }
    else
    {

#ifdef __USE_PIXSRC_CUDA__

        cuda_vec = (PS_FPT*)0;
        cusp_wrapper     = 0;
        cusp_wrapper_ptr = 0;

        // CUDA mutex init
        pthread_mutex_init( &pv_mutex, NULL );

        up2date = pvs_neither;

        if( are_we_using_cuda )
        {
            CUDA ps_cuda_malloc( cdata_, &cuda_vec, size, data_->gpu2use[0] );
            CUDA create_stream( this );
        }

#endif

    }
}

// destructor
pixsrc_vector::~pixsrc_vector()
{
    if( vec )
        MEMORY ps_free( vec );

    bool delete_it = 0;
#ifdef SINGLE_PRECISION
    delete_it = 1;
#endif
    if( delete_it && vec_to_return )
        MEMORY ps_free( vec_to_return );


#ifdef __USE_PIXSRC_CUDA__

    if( are_we_using_cuda )
    {
        //CUDA destroy_cusp_wrapper( this );
        CUDA ps_cuda_free( cuda_vec );
        CUDA destroy_stream( this );
    }

#endif

}

void pixsrc_vector::zeromeout()
{
    // zero out the vector
    if( !are_we_using_cuda )
    {
        std::fill( vec, vec+size, (PS_FPT)0 );
        up2date = pvs_cpu;
    }
    else
    {

#ifdef __USE_PIXSRC_CUDA__
        CUDA fill_vec( this, (PS_FPT)0 );
#endif

    }

    scalar = (PS_FPT)1;
}

const double *pixsrc_vector::get_vec_ptr()
{
    // return a pointer to the vector
    // the vector is a private variable, but this provides direct access

    this->dissolve_scalar();
    this->update_cpu();

#ifdef SINGLE_PRECISION
    if( !vec_to_return )
        MEMORY ps_malloc( &vec_to_return, size );
    double *d_ptr = vec_to_return;
    float  *s_ptr = vec;
    for( PS_SIT j=0; j<size; ++j )
    {
        *d_ptr = (double)*s_ptr;
        ++d_ptr;
        ++s_ptr;
    }
    return (const double*)vec_to_return;
#else
    return (const double*)vec;
#endif

}

void pixsrc_vector::set_status( pv_update_status stat )
{
    // update
    up2date = stat;
}

double pixsrc_vector::get_scalar()
{
    // return scalar multiplying the vector
    return (double)scalar;
}

void pixsrc_vector::mult( double val )
{
    // multiply entire vector by scalar
    scalar *= (PS_FPT)val;
}

PS_SIT pixsrc_vector::get_size()
{
    // get size of vector
    return size;
}

#ifdef __USE_PIXSRC_CUDA__
void pixsrc_vector::cuda_reset()
{
    // reset CUDA operations
    // this function does not destroy "vec" or "cuda_vec"
    // because it is assumed that the answer from some
    // CUDA operation will be placed into those two arrays

    scalar = (PS_FPT)1;
    up2date = pvs_neither;
}
#endif

/*
  void pixsrc_vector::copy(pixsrc_vector *b)
  {
  b->scalar = scalar;
  std::copy( vec, vec+size, b->vec );
  }
*/

void pixsrc_vector::init_cpu_vec()
{
    // allocate vector
    if (!vec)
    {
        if( up2date == pvs_both || up2date == pvs_cpu )
        {
            PRINTER printerror( data_->print2screenname,
                                "pixsrc_vector up2date",
                                cdata_->print2screenmutex );
        }

        pthread_mutex_lock (&initmutex);
        // test again for near-simultaneous calls to init_cpu_vec
        if (!vec)
            MEMORY ps_malloc( &vec, size );
        pthread_mutex_unlock (&initmutex);
    }
}

void pixsrc_vector::update_cpu()
{

    // update vector data on CPU

#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

    if( up2date == pvs_both || up2date == pvs_cpu )
    {
        return;
    }
    else if( up2date == pvs_gpu )
    {

#ifdef __USE_PIXSRC_CUDA__
        CUDA sync_vector( this, ccd_gpu2cpu );
#endif

    }
    else
    {
        init_cpu_vec();
    }

#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU to finish
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

}

void pixsrc_vector::update_gpu()
{

    // update vector data on GPU

#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU to finish
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

    if( up2date == pvs_both || up2date == pvs_gpu )
    {
        return;
    }
    else if( up2date == pvs_cpu )
    {

#ifdef __USE_PIXSRC_CUDA__
        CUDA sync_vector( this, ccd_cpu2gpu );
#endif

    }
    else
    {
        /*
          PRINTER printerror( data_->print2screenname,
          "pixsrc_vector: up2date error: gpu",
          cdata_->print2screenmutex );
        */
    }

#ifdef __USE_PIXSRC_CUDA__
    // wait for GPU to finish
    if( are_we_using_cuda )
        CUDA wait_for_stream( this );
#endif

}

void pixsrc_vector::set(PS_SIT y, double val)
{
    // make these calls the user has to make before doing any setting
    this->init_cpu_vec();
    up2date = pvs_cpu;

    vec[y] = (PS_FPT)val/scalar;
}

void pixsrc_vector::set_scalar( double val )
{
    scalar = (PS_FPT)val;
}

double pixsrc_vector::get(PS_SIT y)
{
    // make this a call the user has to make before doing any setting
    this->update_cpu();

    return (double)(vec[y]*scalar);
}
double pixsrc_vector::max()
{
    // find maximum value of vector

    PS_FPT max_ret=0;

    if( !are_we_using_cuda )
    {
        this->update_cpu();

        PS_FPT max = vec[0];
        if( scalar > (PS_FPT)0 )
        {
            for( PS_SIT j=1; j<size; ++j )
                if( vec[j] > max )
                    max = vec[j];
        }
        else
        {
            for( PS_SIT j=1; j<size; ++j )
                if( vec[j] < max )
                    max = vec[j];
        }

        max_ret = max*scalar;
    }
    else
    {

#ifdef __USE_PIXSRC_CUDA__
        // update is called in the ps_cuda_max function
        max_ret = CUDA ps_cuda_max( this );
#endif

    }

    return (double)max_ret;
}
void pixsrc_vector::dissolve_scalar()
{
    // multiply each element of vectcor by scalar and set scalar to 1
    if( !are_we_using_cuda )
    {
        if( scalar != (PS_FPT)1 )
        {
            this->update_cpu();

            std::transform (vec, vec+size, vec,
                            std::bind1st(std::multiplies<PS_FPT>(),scalar));
            scalar  = (PS_FPT)1;
            up2date = pvs_cpu;
        }
    }
    else
    {
#ifdef __USE_PIXSRC_CUDA__

        // update gpu is called and set in dissolve_scalar
        CUDA dissolve_scalar( this );

#endif

    }
}

double pixsrc_vector::innerproduct(pixsrc_vector *b)
{
    // compute inner product of two vectors

    PS_FPT res = 0;

    if( !are_we_using_cuda )
    {
        this->update_cpu();
        b   ->update_cpu();

        PS_FPT sum = 0;
        for(PS_SIT g=0; g<size; g++)
            sum += vec[g]*b->vec[g];
        res = scalar*b->scalar*sum;
    }
    else
    {
#ifdef __USE_PIXSRC_CUDA__

        // update is called in CUDA innerproduct
        res = CUDA innerproduct( this, b );

#endif
    }

    return (double)res;
}

void pixsrc_vector::plus(pixsrc_vector *b, pixsrc_vector *c, double scalar1_, double scalar2_)
{
    // add two vectors
    // the answer can be stored into "c" or into this vector

    PS_FPT scalar1 = (PS_FPT)scalar1_;
    PS_FPT scalar2 = (PS_FPT)scalar2_;

    if( !are_we_using_cuda )
    {
        this->update_cpu();
        b   ->update_cpu();

        // if storing result in third vector
        if( c )
        {
            c->init_cpu_vec();

            // if both vectors "vec" array are multiplied by same scalar
            if( scalar1*scalar == scalar2*b->scalar )
            {
                for(PS_SIT g=0; g<size; g++)
                    c->vec[g] = vec[g] + b->vec[g];
                c->scalar = scalar1*scalar;
            }
            else
            {
                for(PS_SIT g=0; g<size; g++)
                    c->vec[g] = scalar1*scalar*vec[g] + scalar2*b->scalar*b->vec[g];
                c->scalar = (PS_FPT)1.0;
            }

            c->up2date = pvs_cpu;
        }
        // if storing result in this vector
        else
        {
            // if both vectors "vec" array are multiplied by same scalar
            if( scalar1*scalar == scalar2*b->scalar )
            {
                for(PS_SIT g=0; g<size; g++)
                    vec[g] = vec[g] + b->vec[g];
                scalar = scalar1*scalar;
            }
            else
            {
                for(PS_SIT g=0; g<size; g++)
                    vec[g] = scalar1*scalar*vec[g] + scalar2*b->scalar*b->vec[g];
                scalar = (PS_FPT)1.0;
            }

            this->up2date = pvs_cpu;
        }
    }
    else
    {

#ifdef __USE_PIXSRC_CUDA__
        CUDA plus( this, b, c, scalar1, scalar2 );
#endif

    }
}
void pixsrc_vector::minus(pixsrc_vector *b, pixsrc_vector *c, double scalar1_, double scalar2_)
{
    // subtract two vectors
    // the answer can be stored into "c" or into this vector

    PS_FPT scalar1 = (PS_FPT)scalar1_;
    PS_FPT scalar2 = (PS_FPT)scalar2_;

    if( !are_we_using_cuda )
    {
        this->update_cpu();
        b   ->update_cpu();

        // if storing result in third vector
        if( c )
        {
            c->init_cpu_vec();

            // if both vectors "vec" array are multiplied by same scalar
            if( scalar1*scalar == scalar2*b->scalar )
            {
                for(PS_SIT g=0; g<size; g++)
                    c->vec[g] = vec[g] - b->vec[g];
                c->scalar = scalar1*scalar;
            }
            else
            {
                for(PS_SIT g=0; g<size; g++)
                    c->vec[g] = scalar1*scalar*vec[g] - scalar2*b->scalar*b->vec[g];
                c->scalar = (PS_FPT)1.0;
            }

            c->up2date = pvs_cpu;
        }
        // if storing result in this vector
        else
        {
            // if both vectors "vec" array are multiplied by same scalar
            if( scalar1*scalar == scalar2*b->scalar )
            {
                for(PS_SIT g=0; g<size; g++)
                    vec[g] = vec[g] - b->vec[g];
                scalar = scalar1*scalar;
            }
            else
            {
                for(PS_SIT g=0; g<size; g++)
                    vec[g] = scalar1*scalar*vec[g] - scalar2*b->scalar*b->vec[g];
                scalar = (PS_FPT)1.0;
            }

            this->up2date = pvs_cpu;
        }
    }
    else
    {

#ifdef __USE_PIXSRC_CUDA__
        CUDA minus( this, b, c, scalar1, scalar2 );
#endif

    }

}
