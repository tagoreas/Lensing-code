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



#ifndef PIXSRC_PRINTER_HPP_
#define PIXSRC_PRINTER_HPP_

#define PRINTER pixsrc_printer::

#include "pixsrc_matrix.hpp"
#include "pixsrc_constants.hpp"
#include <pthread.h>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>

using std::string;
using std::vector;

class pixsrc_printer
{
public:

    static void openoutstream    ( string, string, bool, string, std::ofstream *stream, bool     );
    static void closeoutstream   ( std::ofstream*                                              );
    static void print            ( PS_SIT, string, string, bool, string, PS_SIT, PS_SIT, pixsrc_matrix*, bool );
    static void printfitssrcplane( PS_SIT, string, string, string, bool, string, PS_SIT,
                                   PS_SIT, double*, double*, pthread_mutex_t*, bool, bool         );
    static void printfitsimgplane( PS_SIT, string, string, string, bool, string, const double*,
                                   double*, PS_SIT, PS_SIT, double, PS_SIT*, pthread_mutex_t*,bool, bool);
    static void print2screen     ( string, string, pthread_mutex_t*                            );
    static void printwarning     ( string, string, pthread_mutex_t*                            );
    static void printerror       ( string, string, pthread_mutex_t*                            );
    static void printtristruct   ( struct triangulateio*,  pthread_mutex_t*                    );
    static PS_SIT  remove_file      ( PS_SIT, string, string, bool, string, bool                     );

    // template functions are defined in pixsrc_printer_templates.cpp

    template <typename T>
    static PS_SIT  writeoutstream   ( T, std::ofstream*, pthread_mutex_t*, PS_SIT, vector<string>*                         );
    template <typename T>
    static PS_SIT  writeoutstream   ( T*, PS_SIT, std::ofstream*, pthread_mutex_t*, PS_SIT, vector<string>*                   );
    template <typename T>
    static PS_SIT  print            ( PS_SIT, string, string, bool, string, T, bool, PS_SIT, vector<string>*                  );
    template <typename T>
    static PS_SIT  print            ( PS_SIT, string, string, bool, string, T*,  PS_SIT, bool, PS_SIT, vector<string>*           );
    template <typename T>
    static PS_SIT  print            ( PS_SIT, string, string, bool, string, T**, PS_SIT, PS_SIT, bool, PS_SIT, vector<string>*      );
    template <typename T>
    static PS_SIT  print            ( PS_SIT, string, string, bool, string, vector<T>, bool, PS_SIT, vector<string>*          );
    template <typename T>
    static PS_SIT  print            ( PS_SIT, string, string, bool, string, vector< vector<T> >, bool, PS_SIT, vector<string>*);
    template <typename T>
    static PS_SIT  printbinary      ( PS_SIT, string, string, bool, string, const T*, size_t, bool);
    template <typename T>
    static PS_SIT  readbinary       ( string, T*, size_t, pthread_mutex_t*);

};

#endif
