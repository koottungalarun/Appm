# This is a makefile for compiling on the Cluster at OST-IET


# These are the compiler flags in MSVC. 
#       -> and its Linux counterparts
# 
# /I"C:\Program Files\eigenLinAlg"  
# /permissive- 
# /GS   -> -fstack-security-check
# /W3   -> -w3    enable remarks,warnings,and errors
# /Gy   -> -fdata-sections -ffunction-sections -W1, --gc-sections   -> remove dead code from executable
# /Zi   -> Guided auto parallelism, Optimization flag MUST be specified
# /O2   -> optimization for speed
# /Fd"x64\Release\vc141.pdb" -> Windows only: program database, for debugging purpose
# /fp:precise   -> value-safe optimizations
# /D "H5_BUILT_AS_DYNAMIC_LIB"    -> define a macro
# /Qipo -> enable inter-process optimization
# /Oi   -> enable inline expansion of intrinsic functions
# /MD   -> search for unresolved references in dynamic-link run-time library
# /std:c++17 
# /FC   -> only Visual Studio: display full path of source files passed to compiler diagnostics
# /Zc:wchar_t    -> only Microsoft: ANSI C standard conformance
# /Zc:forScope 
# /EHsc -> only Visual Studio: exception handling type
# /nologo -> suppress display of compiler version information (Visual Studio only?)
# /D "_MBCS" 
# /Fa"x64\Release\" 
# /Fo"x64\Release\" 
# /Qprof-dir "x64\Release\" 
# /Fp"x64\Release\Appm.pch" 


# Name of executable
EXE=appm

# C++ Compiler
CC=icpc



# Include path to EIGEN header files
EIGEN_INC= -I/home/HSR/rfuchs/git/eigen-3.3.7/

# Compiler flags:
# -std:c++17     C++ standard
# -w             Enable diagnostics for errors, disable warnings
# -O0, -O1, -O2  Optimization flags
# -inline-forceinline    highly recommended by EIGEN if Intel compiler is used
CFLAGS=-std=c++17 -w -O2 -inline-forceinline 

# Help on Makefile:
# https://stackoverflow.com/questions/53136024/makefile-to-compile-all-c-files-without-needing-to-specify-them/53138757

# Discover all source files
SRCS := $(wildcard *.cpp)

# Transform name of source files to name of object files
OBJS := $(patsubst %.cpp,%.o,$(SRCS))


# HDF5
HDF_INSTALL=/home/HSR/rfuchs/software/hdf5-1.12.0-linux-centos7-x86_64-shared
HDF_INC = -I$(HDF_INSTALL)/include
HDF_LIB := -L$(HDF_INSTALL)/lib
LIBSHDF := $(HDF_LIB) $(HDF_LIB)/libhdf5.so $(HDF_LIB)/libhdf5_cpp.so 

# Intel MKL
MKL_LIB=-L$(MKLROOT)/lib/intel64
MKL_INCLUDE=-I$(MKLROOT)/include
MKL_LINK=-L${MKLROOT}/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread -lm -ldl
#MKL_LINK=-L$(MKLROOT)/lib/ia32 -lmkl_intel -lmkl_sequential -lmkl_core -lpthread -lm -ldl


INCFLAGS := $(HDF_INC) $(EIGEN_INC) $(MKL_INCLUDE)
LDFLAGS := $(HDF_LIB) $(LIBSHDF) $(MKL_LIB)
LDLIBS := -lsz -lz -lhdf5 -lhdf5_cpp $(MKL_LINK) -lm 


all: $(EXE)

$(EXE): $(OBJS)
	export LD_LIBRARY_PATH=$(HDF_INSTALL)/lib:$(MKLROOT)/lib/intel64; \
	export INTEL_LICENSE_FILE=/home/HSR/rfuchs/Licences; \
	module load intelmpi/2018.1.163_64; \
	$(CC) -o $(EXE) $(LDFLAGS) $^ $(LDLIBS)

%.o: %.cpp
	export LD_LIBRARY_PATH=$(HDF_INSTALL)/lib:$(MKLROOT)/lib/intel64; \
	export INTEL_LICENSE_FILE=/home/HSR/rfuchs/Licences; \
	module load intelmpi/2018.1.163_64; \
	$(CC) $(INCFLAGS) $(CFLAGS) -c $< -o $@

clean: 
	rm *.o
	
run: 
	export LD_LIBRARY_PATH=$(HDF_INSTALL)/lib:$(MKLROOT)/lib/intel64; \
	./$(EXE)