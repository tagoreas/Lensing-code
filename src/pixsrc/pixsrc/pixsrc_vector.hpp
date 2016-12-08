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



#ifndef PIXSRC_VECTOR_HPP_
#define PIXSRC_VECTOR_HPP_

#define VECTOR pixsrc_vector

#include "pixsrc_inputdata.hpp"

#include <pthread.h>

// flag for where the data is up to date
typedef enum
{
    pvs_neither = -1,
    pvs_both    = 0,
    pvs_cpu     = 1,
    pvs_gpu     = 2

} pv_update_status;


class pixsrc_vector
{

    // allow pixsrc_matrix and CUDA to access private variables directly
    friend class pixsrc_matrix;

#ifdef __USE_PIXSRC_CUDA__
    friend class pixsrc_cuda;
#endif

public:

    pixsrc_vector( commoninputdata*, inputdata*, PS_SIT);
    ~pixsrc_vector();

    bool are_we_using_cuda;

    void    init_cpu_vec   ();
    void    update_cpu     ();
    void    update_gpu     ();
    double  get_scalar     ();
    double  max            ();
    PS_SIT     get_size       ();
    void    dissolve_scalar();
    void    zeromeout      ();
    const
    double *get_vec_ptr    ();
    double  get            ( PS_SIT                                            );
    void    set_scalar     ( double                                         );
    void    mult           ( double                                         );
    void    set            ( PS_SIT, double                                    );
    double  innerproduct   ( pixsrc_vector*                                 );
    void    set_status     ( pv_update_status                               );
    void    plus           ( pixsrc_vector*, pixsrc_vector*, double, double );
    void    minus          ( pixsrc_vector*, pixsrc_vector*, double, double );

    PS_FPT *vec;
    PS_SIT size;                // # of {columns, rows, nonzeroes}
private:

    commoninputdata *cdata_;
    inputdata *data_;
    double *vec_to_return; // constant pointer to vector to return
    PS_FPT scalar;           // scalar to multiply matrix by;

    pv_update_status up2date;
    pthread_mutex_t initmutex;

#ifdef __USE_PIXSRC_CUDA__

    PS_FPT *cuda_vec;
    pthread_mutex_t pv_mutex;
    void cuda_reset();

    // CUSP object
    void *cusp_wrapper;
    void *cusp_wrapper_ptr;

    // CUDA stream
    void *stream;

#endif

    pixsrc_vector operator=( const pixsrc_vector& );
    pixsrc_vector( const pixsrc_vector& );

};

#endif
