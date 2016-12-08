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



#ifndef PIXSRC_MEMORY_HPP_
#define PIXSRC_MEMORY_HPP_

#define MEMORY pixsrc_memory::

#include "pixsrc_inputdata.hpp"
#include "pixsrc_lensvar.hpp"

#ifdef PIXSRC_MEMORYDEBUG
#include <pthread.h>
#endif


class pixsrc_memory
{

public:

    static void  triangulatestructinit    (struct triangulateio*);
    static void  triangulatestructdestruct(struct triangulateio*);
    static void  freememory               (inputdata*,commoninputdata*,lensvar*);
    static void  free_matrices            (inputdata*,commoninputdata*,lensvar*);
    static void  initialize               (inputdata*,commoninputdata*,lensvar*);

    static void check_ptr (bool*, size_t);

    static void check_ptr (char*, size_t);

    static void check_ptr (signed char*,          size_t);
    static void check_ptr (signed short int*,     size_t);
    static void check_ptr (signed int*,           size_t);
    static void check_ptr (signed long int*,      size_t);
    static void check_ptr (signed long long int*, size_t);

    static void check_ptr (unsigned char*,          size_t);
    static void check_ptr (unsigned short int*,     size_t);
    static void check_ptr (unsigned int*,           size_t);
    static void check_ptr (unsigned long int*,      size_t);
    static void check_ptr (unsigned long long int*, size_t);

    static void check_ptr (float*,       size_t);
    static void check_ptr (double*,      size_t);
    static void check_ptr (long double*, size_t);

    static void check_ptr (void*, size_t);

    template <typename T> static void ps_realloc( T**   , size_t                );
    template <typename T> static void ps_malloc( T**   , size_t                );
    template <typename T> static void ps_malloc( T***  , size_t, size_t           );
    template <typename T> static void ps_malloc( T**** , size_t, size_t, size_t      );
    template <typename T> static void ps_malloc( T*****, size_t, size_t, size_t, size_t );
    template <typename T> static void ps_free  ( T*                         );
    template <typename T> static void ps_free  ( T**   , size_t                );
    template <typename T> static void ps_free  ( T***  , size_t, size_t           );
    template <typename T> static void ps_free  ( T**** , size_t, size_t, size_t      );

#ifdef PIXSRC_MEMORYDEBUG
    static PS_unsignedSIT nummallocs;
    static PS_unsignedSIT numfrees;
    static pthread_mutex_t freemutex;
    static pthread_mutex_t mallocmutex;
#endif

};

#endif
