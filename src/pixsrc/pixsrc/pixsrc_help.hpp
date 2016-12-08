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



#ifndef PIXSRC_HELP_HPP_
#define PIXSRC_HELP_HPP_

#define HELP pixsrc_help::

#include <vector>
#include <string>

class pixsrc_help
{

public:

    static void general ();
    static void parms ();
    static void parms_more ();
    static void ps_template ();
    static void manual ();
    static std::vector<std::string> input (PS_SIT printit);
    static std::vector<std::string> output (PS_SIT printit);

private:

    static void template_1 ();

};

#endif
