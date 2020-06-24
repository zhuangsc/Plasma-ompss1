###
#
# @file make.inc.gcc.mkl_x86_64
#
#  PLASMA is a software package provided by Univ. of Tennessee,
#  Univ. of California Berkeley and Univ. of Colorado Denver
#
# @version 2.6.0
# @author Piotr Luszczek
# @author Mathieu Faverge
# @date 2010-11-15
#
###
#
# OS:       Linux x86-64
# Compiler: GCC
# BLAS:     MKL
#
###

# PLASMA install directory (./install by default)
prefix      = /home/bsc28/bsc28687/minotauro/plasma_2.6.0/install

# Compilers, linker/loader, archiver, options
CC          = gcc
MCC			= mcc --ompss --debug
#MCC			= ~bsc15105/mcxx-workspace/install-mt/bin/mcc --ompss --debug
FC          = gfortran
LOADER      = $(FC)

ARCH        = ar
ARCHFLAGS   = cr
RANLIB      = ranlib

CFLAGS      = -O2 -DADD_
FFLAGS      = -O2
LDFLAGS     =

# Enable/Disable F90 compilation
PLASMA_F90 = 0

# BLAS library
CFLAGS += -DPLASMA_WITH_MKL -I${MKLROOT}/include -DPLASMA_WITH_SMP
LIBBLAS = -L$(MKLROOT)/lib/intel64 -lmkl_gf_lp64 -lmkl_core -lmkl_sequential -lpthread

# CBLAS library
# CBLAS is included in MKL

# LAPACK and TMG library
# LAPACK and TMG are included in MKL

# LAPACK C Wrapper library
INCCLAPACK  = -I/usr/include -I/gpfs/apps/MN3/LAPACK/SRC/lapack-3.5.0/lapacke/include
LIBCLAPACK  = -L/gpfs/apps/MN3/LAPACK/3.5.0/GCC/lib -llapacke
