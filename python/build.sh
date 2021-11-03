#!/bin/bash
# Build GCM version of RRTMG_LW 
# without Monte-Carlo Independent Column Approximation
# and using absorption coefficients from FORTRAN data file

# Base directory
export BASE=`pwd`
# Location for source files
export SRC=`pwd`/src
# Build directory
export BUILD=`pwd`/build
# Make utility
export MAKE=make
# Fortran compiler (compile-only)
export FC='gfortran -c -fPIC'
# Compilation flags (-O0 for fast test compilation)
export FFLAGS='-O0'
# f2py utility
# Note: can run f2py -c --help-fcompiler to list
# available fortran compilers
export F2PY='f2py3 -c --fcompiler=gnu95'

# Remove build directory if present
if [ -d "$BUILD" ]; then rm -r "$BUILD"; fi

# Copy required files (see README)
mkdir -p $SRC
cp ../src/rrtmg_lw_rad.nomcica.f90 $SRC/rrtmg_lw_rad.f90
cp ../src/rrtmg_lw_cldprop.f90 $SRC/
cp ../src/rrtmg_lw_init.f90 $SRC/
cp ../src/rrtmg_lw_k_g.f90 $SRC/
cp ../src/rrtmg_lw_rtrn.f90 $SRC/
cp ../src/rrtmg_lw_rtrnmr.f90 $SRC/
cp ../src/rrtmg_lw_setcoef.f90 $SRC/
cp ../src/rrtmg_lw_taumol.f90 $SRC/
cp ../modules/parkind.f90 $SRC/
cp ../modules/parrrtm.f90 $SRC/
cp ../modules/rrlw_cld.f90 $SRC/
cp ../modules/rrlw_con.f90 $SRC/
cp ../modules/rrlw_kg*.f90 $SRC/
cp ../modules/rrlw_ncpar.f90 $SRC/
cp ../modules/rrlw_ref.f90 $SRC/
cp ../modules/rrlw_tbl.f90 $SRC/
cp ../modules/rrlw_vsn.f90 $SRC/
cp ../modules/rrlw_wvn.f90 $SRC/

# Create build directory
mkdir -p $BUILD
cd $BUILD

# Generate dependency list
echo $SRC > Filepath
perl $BASE/mkSrcfiles.pl
perl $BASE/mkDepends.pl Filepath Srcfiles > Depends

# Compile
$MAKE -f $BASE/Makefile
