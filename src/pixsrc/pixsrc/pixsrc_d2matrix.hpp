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



#ifndef PIXSRC_D2MATRIX_HPP_
#define PIXSRC_D2MATRIX_HPP_

#define D2MATRIX pixsrc_d2matrix

// this is a multi-dimensional matrix.
// it mimics a nested std::vectors of dimension 2.
// no safety checking on indices is performed, so don't
// try to acces an index that hasn't been set

template <class T>
class pixsrc_d2matrix
{
public:
    pixsrc_d2matrix();
    pixsrc_d2matrix( PS_SIT dim );
    ~pixsrc_d2matrix();

    void pushback( PS_SIT, T );
    void set( PS_SIT*, T );
    void set( PS_SIT, PS_SIT, T );
    T get( PS_SIT, PS_SIT );
    T get( PS_SIT* );
    T* get_pointer (PS_SIT);
    PS_SIT size ( PS_SIT );
    void sort( PS_SIT );

private:
    // holds the actual 2d structure
    T** mdm;
    // holds value of largest element set in 1st dimension
    PS_SIT extension1;
    // holds capacity of 1st dimension
    PS_SIT capacity1;
    // holds value of largest element set in 2nd dimension
    PS_SIT *sizearr;
    // holds capacity of 2nd dimension
    PS_SIT *capacity;

    pixsrc_d2matrix operator=(const pixsrc_d2matrix &mdm);
    pixsrc_d2matrix(const pixsrc_d2matrix &mdm);
};

#endif
