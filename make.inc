############################################################################
#
#  Program:         SuperLU
#
#  Module:          make.inc
#
#  Purpose:         Top-level Definitions
#
#  Creation date:   May 10, 2015
#
#  Modified:	    
#		    
#
############################################################################
#
#  The name of the libraries to be created/linked to
#
SuperLUroot	= /home/xiaoye/Dropbox/Codes/SuperLU/SuperLU/build
SUPERLULIB   	= $(SuperLUroot)/SRC/libsuperlu_5.1.0.a

TMGLIB       	= libtmglib.a

# BLASLIB   	= $(SuperLUroot)/lib/libblas.a
# BLASDEF 	= -DUSE_VENDOR_BLAS

LIBS		= $(SUPERLULIB) blas

#
#  The archiver and the flag(s) to use when building archive (library)
#  If your system has no ranlib, set RANLIB = echo.
#
ARCH         = /usr/bin/ar
ARCHFLAGS    = cr
RANLIB       = /usr/bin/ranlib

CC           = /usr/bin/cc
CFLAGS 	     = -O3  -DPRNTlevel=0 -DAdd_ 
# CFLAGS       += -D
# CFLAGS       +=  
NOOPTS       = -O0
FORTRAN	     = /usr/bin/gfortran

LOADER       = $(CC)
LOADOPTS     =
