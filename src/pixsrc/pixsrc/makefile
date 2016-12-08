
ifeq ($(USE_PIXSRC_SINGLE_PRECISION),YES)
PREC     =		-DSINGLE_PRECISION
else
PREC     =		-DDOUBLE_PRECISION
endif

ifeq ($(USE_PIXSRC_SIT_LL),YES)
PREC     +=             -DINTEGER_TYPE_LL
else
ifeq ($(USE_PIXSRC_SIT_L),YES)
PREC     +=             -DINTEGER_TYPE_L
else
PREC     +=             -DINTEGER_TYPE_I
endif
endif

ifeq ($(USE_PIXSRC_GRAVLENS),YES)
PREC     +=             -DPS_HAVE_GRAVLENS
endif
ifeq ($(USE_PIXSRC_TRIAXIAL),YES)
PREC     +=             -DPS_HAVE_TRIAXIAL
endif

CC	 =		gcc
CXX	 =		g++

CXXFLAGS = 		-O2 -g -fPIC -Wall -fmessage-length=0 $(PREC) \
			-I/usr/include/cfitsio \
			-I/usr/include/suitesparse \
			-I../../pixsrc/wcstools/libwcs \
			-I./

CFLAGS = 		-O2 -g -fPIC -Wall -fmessage-length=0 -Wno-strict-overflow $(PREC)






ifneq ($(PIXSRC_OS),DARWIN)

ifeq ($(PIXSRC_SIZE),32)
                        CFLAGS    +=	-m32
                        CXXFLAGS +=	-m32
else
                        CFLAGS    +=	-m64
                        CXXFLAGS +=	-m64
endif

else
			CFLAGS    +=    -arch i386 -arch x86_64
                        CXXFLAGS +=     -arch i386 -arch x86_64
endif






ifeq ($(USE_PIXSRC_CUDA),YES)
CXXFLAGS +=		-D__USE_PIXSRC_CUDA__

# OS Name (Linux or Darwin)
        CUDA_OSUPPER = $(shell uname -s 2>/dev/null | tr [:lower:] [:upper:])
        CUDA_OSLOWER = $(shell uname -s 2>/dev/null | tr [:upper:] [:lower:])

# Flags to detect 32-bit or 64-bit OS platform
        CUDA_OS_SIZE = $(shell uname -m | sed -e "s/i.86/32/" -e "s/x86_64/64/")
        CUDA_OS_ARCH = $(shell uname -m | sed -e "s/i386/i686/")

# These flags will override any settings
ifeq ($(i386),1)
                CUDA_OS_SIZE = 32
                CUDA_OS_ARCH = i686
endif

ifeq ($(x86_64),1)
                CUDA_OS_SIZE = 64
                CUDA_OS_ARCH = x86_64
endif

# Flags to detect either a Linux system (linux) or Mac OSX (darwin)
        CUDA_DARWIN = $(strip $(findstring DARWIN, $(CUDA_OSUPPER)))

# Location of the CUDA Toolkit binaries and libraries
        CUDA_PATH       ?= /usr/local/cuda
        CUDA_INC_PATH   ?= $(CUDA_PATH)/include
        CUDA_BIN_PATH   ?= $(CUDA_PATH)/bin
ifneq ($(CUDA_DARWIN),)
                CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib
else
ifeq ($(CUDA_OS_SIZE),32)
                        CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib
else
                        CUDA_LIB_PATH  ?= $(CUDA_PATH)/lib64
endif
endif

# Common binaries
        NVCC            ?= $(CUDA_BIN_PATH)/nvcc

# OS-specific build flags
ifneq ($(CUDA_DARWIN),)
                CUDA_LDFLAGS   := -Xlinker -rpath $(CUDA_LIB_PATH) -L$(CUDA_LIB_PATH) \
                                -framework CUDA -lcudart -lcusparse
                CUDA_CCFLAGS   := -arch $(CUDA_OS_ARCH)
else
		CUDA_LDFLAGS   := -L$(CUDA_LIB_PATH) -lcudart -lcusparse
ifeq ($(CUDA_OS_SIZE),32)
                        CUDA_CCFLAGS   := -m32
else
                        CUDA_CCFLAGS   := -m64
endif
endif

# OS-architecture specific flags
ifeq ($(CUDA_OS_SIZE),32)
                NVCCFLAGS := -m32 -arch=sm_20 $(PREC) -I../../pixsrc/cusp
else
                NVCCFLAGS := -m64 -arch=sm_20 $(PREC) -I../../pixsrc/cusp
endif

# Common includes and paths for CUDA
        CUDA_INCLUDES      := -I$(CUDA_INC_PATH)
        CXXFLAGS += $(CUDA_CCFLAGS) $(CUDA_INCLUDES)
endif















PIXSRC_LIB = 		libpixsrc.a
PIXSRC_LIB_SHARED =	libpixsrc.so

PIXSRC_OBJ = 		pixsrc_cartesian.o \
			pixsrc_adaptive.o \
			pixsrc_irrcart.o \
			pixsrc_shapelets.o \
			pixsrc_shapelets_operations.o \
			pixsrc_common.o \
			pixsrc_common_adaptive.o \
			pixsrc_reg_sersic_min.o \
			pixsrc_reg_sersic_max.o \
			pixsrc_reg_greengauss_1.o \
			pixsrc_reg_greengauss_2.o \
			pixsrc_reg_adaptive_0.o \
			pixsrc_reg_adaptive_12.o \
			pixsrc_reg_cartesian.o \
			pixsrc_reg_shapelets.o \
			pixsrc_constants.o \
			pixsrc_constants_shapelets.o \
			pixsrc_geometry.o \
			pixsrc_operations.o \
			pixsrc_wcs.o \
			pixsrc_printer.o \
			pixsrc_triangle.o \
			pixsrc_vector.o \
			pixsrc_matrix.o \
			pixsrc_communicator.o \
			pixsrc_tps.o \
			pixsrc_init_other.o \
			pixsrc_init_uvplane.o \
			pixsrc_init_parse.o \
			pixsrc_init_analyticsrc.o \
			pixsrc_init_parms.o \
			pixsrc_memory.o \
			pixsrc_external.o \
			pixsrc_analytic_source.o \
			pixsrc_triangulation.o \
			pixsrc_statistic.o \
			pixsrc_nonparamlens.o \
			pixsrc_help.o \
			pixsrc_clipper.o

ifeq ($(USE_PIXSRC_CUDA),YES)
PIXSRC_OBJ +=		pixsrc_cuda.o \
			pixsrc_cuda_helper.o \
			pixsrc_init_cuda.o \
			pixsrc_cuda_kernels.o
endif

PIXSRC_SRC_TEMPL =	pixsrc_memory_templates.cpp \
			pixsrc_operations_templates.cpp \
			pixsrc_printer_templates.cpp \
			pixsrc_d2matrix.cpp \
			pixsrc_d4matrix.cpp

PIXSRC_OBJ_SHARED_ST =	$(PIXSRC_OBJ) \
			../../pixsrc/wcstools/libwcs/libwcs.a

PIXSRC_OBJ_SHARED_SH = 	-lgsl -lumfpack -lcfitsio -lpthread
ifeq ($(USE_PIXSRC_BLAS_OPENBLAS),YES)
PIXSRC_OBJ_SHARED_SH +=	-lopenblas
endif
ifeq ($(USE_PIXSRC_BLAS_ATLAS),YES)
ifeq ($(PIXSRC_SIZE),32)
PIXSRC_OBJ_SHARED_SH +=	-L/usr/lib/atlas -ltatlas
else
PIXSRC_OBJ_SHARED_SH +=	-L/usr/lib64/atlas -ltatlas
endif
endif

all:			$(PIXSRC_LIB_SHARED)

#pixsrc_help.o:		CXXFLAGS += -fno-var-tracking

.cxx.o:
			$(CXX) $(CXXFLAGS) $(CFLAGS) -c $<

$(PIXSRC_OBJ):		$(PIXSRC_SRC_TEMPL)

#$(PIXSRC_LIB):		$(PIXSRC_OBJ)
#			ar -ruv $(PIXSRC_LIB) $(PIXSRC_OBJ)
#			ranlib  $(PIXSRC_LIB)

ifneq ($(PIXSRC_OS),DARWIN)
$(PIXSRC_LIB_SHARED):	$(PIXSRC_OBJ_SHARED_ST)
			$(CXX) -shared -fPIC -o $@ -Wl,-Bstatic $(PIXSRC_OBJ_SHARED_ST) -Wl,-Bdynamic $(PIXSRC_OBJ_SHARED_SH)
else
$(PIXSRC_LIB_SHARED):	$(PIXSRC_OBJ_SHARED_ST)
			$(CXX) -dynamiclib -fPIC -o $@ $(PIXSRC_OBJ_SHARED_ST) $(PIXSRC_OBJ_SHARED_SH) -Wl,-undefined,dynamic_lookup 
endif

ifeq ($(USE_PIXSRC_CUDA),YES)
pixsrc_cuda_kernels.o:	pixsrc_cuda_kernels.cu pixsrc_cuda_kernels.hpp
			$(NVCC) $(NVCCFLAGS) -Xcompiler -fPIC -c pixsrc_cuda_kernels.cu
endif
