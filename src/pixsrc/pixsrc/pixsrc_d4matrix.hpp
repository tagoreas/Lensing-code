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



#ifndef PIXSRC_D4MATRIX_HPP_
#define PIXSRC_D4MATRIX_HPP_

#define D4MATRIX pixsrc_d4matrix

// this is a multi-dimensional matrix.
// it mimics a nested std::vectors of dimension 4.
// no safety checking on indices is performed, so don't
// try and access an element that hasn't been set yet!!

template <class T>
class pixsrc_d4matrix
{
public:
    pixsrc_d4matrix();
    pixsrc_d4matrix( PS_SIT dim );
    ~pixsrc_d4matrix();

    void set( PS_SIT*, T );
    void set( PS_SIT, PS_SIT, PS_SIT, PS_SIT, T );
    T get( PS_SIT, PS_SIT, PS_SIT, PS_SIT );
    T get( PS_SIT* );
    PS_SIT size ( PS_SIT*, PS_SIT );
    PS_SIT size();
    PS_SIT size( PS_SIT );
    PS_SIT size( PS_SIT, PS_SIT );
    PS_SIT size( PS_SIT, PS_SIT, PS_SIT );

private:
    // holds the actual 2d structure
    T**** mdm;
    // holds value of largest element set in 1st dimension
    PS_SIT extension1;
    // holds capacity of 1st dimension
    PS_SIT capacity1;

    // holds value of largest element set in 2nd dimension
    PS_SIT *extension2;
    // holds capacity of 2nd dimension
    PS_SIT *capacity2;

    // holds value of largest element set in 3rd dimension
    PS_SIT **extension3;
    // holds capacity of 3rd dimension
    PS_SIT **capacity3;

    // holds value of largest element set in 4th dimension
    PS_SIT ***extension4;
    // holds capacity of 4th dimension
    PS_SIT ***capacity4;

    pixsrc_d4matrix operator=(const pixsrc_d4matrix &mdm);
    pixsrc_d4matrix(const pixsrc_d4matrix &mdm);
};

#endif
