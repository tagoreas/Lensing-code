
1. The libwcs library makefile has been modified. There is a conflict between global variables in lensmodel and libwcs. I've redefined one of the variables in the libwcs routines by throwing in the following compiler flag: "-Dmatinv=matinv_wcs".

2. The CUDA routines were built for CUDA 5.0. Future changes to CUDA might have to be updated in this code.

3. I made a lot changes to the triangle.c code (renamed pixsrc_triangel.c). I've allowed for 64-bit integers.

4. cfitsio must be configured with the --enable-reentrant option (for multithreaded access).
   Alternatively, disable multi-threaded file writing by putting mutex locks around the cfitsio code in pixsrc_printer.cpp

5. SuiteSpare/UMFPACK can use 'long int'. If you need 64-bit integers (most likely for uv-plane modelling)  and 'long int' isn't 64-bit, you will have to compile UMFPACK from source.
