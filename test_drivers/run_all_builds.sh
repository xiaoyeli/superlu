#!/bin/sh
#
# This script runs various builds for SuperLU on a given machine using GCC.
#
# To run this, you must run it as:
#
#   mkdir <some_dir>
#   env TPL_BLAS_LIBRARIES=<some_dir_1>/libblas.so \
#     TPL_LAPACK_LIBRARIES=<some_dir_2>/liblapack.so \
#     <this_dir>/run_all_builds.sh
#
# If a test case passes, it is silent.  If it fails, it returns FAILED

_ABS_FILE_PATH=`readlink -f $0`
_SCRIPT_DIR=`dirname $_ABS_FILE_PATH`
#echo "_SCRIPT_DIR = $_SCRIPT_DIR"

echo "Remove all *.out files!"
rm *.out

echo ""
echo "1) Serial with default system BLAS/LAPACK"

$_SCRIPT_DIR/serial_default_blas_shared_libs.sh
$_SCRIPT_DIR/serial_default_blas_static_libs.sh

echo ""
echo "2) Serial with override of default system BLAS/LAPACK (requires TPL_BLAS_LIBRARIES and TPL_LAPACK_LIBRARIES set in env)"

$_SCRIPT_DIR/serial_override_blas_shared_libs.sh
$_SCRIPT_DIR/serial_override_blas_static_libs.sh

echo ""
echo "3) Serial with user internal CLBAS"

$_SCRIPT_DIR/serial_internal_cblas_shared_libs.sh
$_SCRIPT_DIR/serial_internal_cblas_static_libs.sh

echo ""
echo "4) MPI with override of default system BLAS/LAPACK (requires TPL_BLAS_LIBRARIES and TPL_LAPACK_LIBRARIES set in env)"

$_SCRIPT_DIR/mpi_override_blas_shared_libs.sh
$_SCRIPT_DIR/mpi_override_blas_static_libs.sh
