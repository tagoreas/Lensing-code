#!/bin/csh -f
#
# Compile pixsrc
#

###############################################
### START USER MODIFIABLE SETTINGS ############
###############################################
# compile pixsrc (NO will create a dummy shared library)
setenv USE_PIXSRC          YES

# enable CUDA (GPU code)
setenv USE_PIXSRC_CUDA     NO

# create pixsrc library for use with triaxial lens code or for use with gravlens
setenv USE_PIXSRC_GRAVLENS NO
setenv USE_PIXSRC_TRIAXIAL YES

# define integer type to be long long(LL), long(L) or int(I)
# uv modelling code with large visibility sets will need 64-bit integer
#   long long(LL) requires compiling suitesparse/umfpack from source
setenv USE_PIXSRC_SIT_LL   NO
setenv USE_PIXSRC_SIT_L    YES
setenv USE_PIXSRC_SIT_I    NO

# BLAS library to use: ATLAS or openblas
setenv USE_PIXSRC_BLAS_OPENBLAS YES
setenv USE_PIXSRC_BLAS_ATLAS    NO
###############################################
#### END USER MODIFIABLE SETTINGS #############
###############################################


###############################################
########### START SETUP ENVIRONMENT ###########
###############################################

# set defaults if not yet set
if (! $?USE_PIXSRC) then
    setenv USE_PIXSRC YES
endif
if (! $?USE_PIXSRC_CUDA) then
    setenv USE_PIXSRC_CUDA YES
endif
if (! $?USE_PIXSRC_SINGLE_PRECISION) then
    setenv USE_PIXSRC_SINGLE_PRECISION NO
endif
if (! $?USE_PIXSRC_DOUBLE_PRECISION) then
    setenv USE_PIXSRC_DOUBLE_PRECISION YES
endif

# find out if the system is 32 or 64 bit
if (! $?PIXSRC_SIZE) then
    setenv PIXSRC_SIZE `uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/"`
endif
# find out if this is a Darwin or Linux system
if (! $?PIXSRC_OS) then
    setenv PIXSRC_OS `uname -s | tr '[:lower:]' '[:upper:]'`
endif

# set one last default if not set yet
if (! $?USE_PIXSRC_SIT_I) then
    if (! $?USE_PIXSRC_SIT_L) then
	if (! $?USE_PIXSRC_SIT_LL) then
	    if($PIXSRC_SIZE == 64) then
		echo setting
		setenv USE_PIXSRC_SIT_LL NO
		setenv USE_PIXSRC_SIT_L  YES
		setenv USE_PIXSRC_SIT_I  NO
	    else
		setenv USE_PIXSRC_SIT_LL NO
		setenv USE_PIXSRC_SIT_L  NO
		setenv USE_PIXSRC_SIT_I  YES
	    endif
	endif
    endif
endif

# check for CUDA compatibility
if (! $?USE_PIXSRC_CUDA) then
    if (-e "/usr/local/cuda/bin/nvcc") then
	set cudaversion = `/usr/local/cuda/bin/nvcc --version | \
			   sed -n -e 's/^.*release //p' | \
			   sed -n -e 's/\..*$//p'`
	if ($cudaversion >= 5) then
	    setenv USE_PIXSRC_CUDA YES
	else
	    echo "WARNING: CUDA found, but it is out of date (or too new, perhaps)."
	    setenv USE_PIXSRC_CUDA NO
	endif
    else
	setenv USE_PIXSRC_CUDA NO
    endif
endif

if ($USE_PIXSRC_CUDA == YES) then
    if ($USE_PIXSRC_SIT_I != YES) then
	echo 'Error: cannot define integer type as long or long long if using CUDA'
	exit
    endif
endif

###############################################
############ END SETUP ENVIRONMENT ############
###############################################


###############################################
############### START COMPILING ###############
###############################################

set BASEDIR = `pwd`
    rm -f lib/libpixsrc.so

if($USE_PIXSRC == YES) then

    # compile WCS library
    cd $BASEDIR/src/pixsrc
    set wcstarball = `ls wcstools-*.tar.gz`
    set wcsdir = `ls wcstools-*.tar.gz | sed -n -e 's/\.tar\.gz//p'`
    if (! -d wcstools) then
	rm -rf wcstools $wcsdir
	tar -xvzf $wcstarball
	mv $wcsdir wcstools
    endif
    cd $BASEDIR/src/pixsrc/wcstools/libwcs
    make CFLAGS="-g -O2 -D_FILE_OFFSET_BITS=64 -Dmatinv=matinvwcs -fPIC"

    # compile pixsrc
    cd $BASEDIR/src/pixsrc/pixsrc
    make all
    mv libpixsrc.so ../../../lib/

else

    # compile dummy pixsrc
    cd $BASEDIR/src/pixsrc/pixsrc
    g++ -fPIC -c pixsrc_dummy.cpp
    if ($PIXSRC_OS == DARWIN) then
        g++ -dynamiclib -fPIC -o libpixsrc.so pixsrc_dummy.o 
    else
        g++ -shared -fPIC -o libpixsrc.so pixsrc_dummy.o
    endif
    mv libpixsrc.so ../../../lib/

endif

###############################################
################ END COMPILING ################
###############################################
