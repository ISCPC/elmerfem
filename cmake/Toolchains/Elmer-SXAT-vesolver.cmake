# CMake toolchain file for building with GNU compilers
#
# Author: Mikko Byckling, CSC - IT Center for Science Ltd.
# Version: 0.1

SET(CMAKE_SYSTEM_NAME Linux)
SET(CMAKE_SYSTEM_PROCESSOR x86_64)
SET(CMAKE_SYSTEM_VERSION 1)

# Specify the cross compilers (serial)
SET(CMAKE_C_COMPILER gcc)
SET(CMAKE_Fortran_COMPILER gfortran)
SET(CMAKE_CXX_COMPILER g++)

# Specify the cross compilers (parallel)
SET(MPI_C_COMPILER mpicc)
SET(MPI_CXX_COMPILER mpic++)
SET(MPI_Fortran_COMPILER mpif90)

# Compilation flags (i.e. with optimization)
SET(CMAKE_C_FLAGS "-O3 -g -m64 -ftree-vectorize -funroll-loops" CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-O3 -g -m64 -ftree-vectorize -funroll-loops" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-O3 -g -m64 -ftree-vectorize -funroll-loops -L${VESOLVER_LIBRARIES} -lvesolver" CACHE STRING "")

SET(MPI_C_LIBRARIES ${SXAT_MPI_LIBRARY_PATH}/lib64/vh/gnu/4.8.5/libmpi.so)
SET(MPI_C_INCLUDE_PATH ${SXAT_MPI_LIBRARY_PATH}/include)

SET(MPI_CXX_LIBRARIES ${SXAT_MPI_LIBRARY_PATH}/lib64/vh/gnu/4.8.5/libmpi.so)
SET(MPI_CXX_INCLUDE_PATH ${SXAT_MPI_LIBRARY_PATH}/include)

SET(MPI_Fortran_LIBRARIES "${SXAT_MPI_LIBRARY_PATH}/lib64/vh/gnu/4.8.5/libmpi.so" "${VESOLVER_LIBRARIES}/libvesolver.so")
SET(MPI_Fortran_INCLUDE_PATH "${SXAT_MPI_LIBRARY_PATH}/include" "${SXAT_MPI_LIBRARY_PATH}/lib64/vh/gnu/4.8.5/module")
SET(MPI_Fortran_LINK_FLAGS "-Wl,-Bdynamic -Wl,--enable-new-dtags -Wl,-rpath=${SXAT_MPI_LIBRARY_PATH}/lib64/vh/gnu/4.8.5 -L${SXAT_MPI_LIBRARY_PATH}/lib64/vh/gnu/4.8.5 -lmpi_mem -L${VESOLVER_LIBRARIES} -lvesolver")
