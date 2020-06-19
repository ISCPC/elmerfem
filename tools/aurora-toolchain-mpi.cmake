# takes NMPI_ROOT and NLC_ROOT from environment
SET(CMAKE_SYSTEM_NAME Linux)
#set(CMAKE_SYSTEM_NAME Aurora-VE)

#set(CMAKE_Fortran_COMPILER /opt/nec/ve/bin/nfort CACHE FILEPATH "Aurora Fortran compiler")
#set(CMAKE_CXX_COMPILER /opt/nec/ve/bin/nc++ CACHE FILEPATH "Aurora C++ compiler")
#set(CMAKE_C_COMPILER /opt/nec/ve/bin/ncc CACHE FILEPATH "Aurora C compiler")
SET(CMAKE_Fortran_COMPILER /home/P1401/bin/nfort)
SET(CMAKE_CXX_COMPILER /home/P1401/bin/nc++)
SET(CMAKE_C_COMPILER /home/P1401/bin/ncc)
set(CMAKE_LINKER /opt/nec/ve/bin/nld CACHE FILEPATH "Aurora linker")
set(CMAKE_AR /opt/nec/ve/bin/nar CACHE FILEPATH "Aurora archiver")
set(CMAKE_RANLIB /opt/nec/ve/bin/nranlib CACHE FILEPATH "Aurora ranlib")

set(OpenMP_Fortran_FLAGS "-fopenmp" CACHE STRING "Flag to enable OpenMP")
set(OpenMP_CXX_FLAGS "-fopenmp" CACHE STRING "Flag to enable OpenMP")
set(OpenMP_C_FLAGS "-fopenmp" CACHE STRING "Flag to enable OpenMP")

#set(MPI_C_COMPILER $ENV{NMPI_ROOT}/bin/mpincc CACHE FILEPATH "")
set(MPI_C_COMPILER /home/P1401/bin/mpincc CACHE FILEPATH "")
set(MPI_C_INCLUDE_PATH $ENV{NMPI_ROOT}/include CACHE FILEPATH "")
set(MPI_C_LIBRARIES $ENV{NMPI_ROOT}/lib64/ve/libmpi.a CACHE FILEPATH "")
set(MPI_C_COMPILE_FLAGS "-D_MPIPP_INCLUDE" CACHE STRING "")

#set(MPI_CXX_COMPILER $ENV{NMPI_ROOT}/bin/mpinc++ CACHE FILEPATH "")
set(MPI_CXX_COMPILER /home/P1401/bin/mpinc++ CACHE FILEPATH "")
set(MPI_CXX_INCLUDE_PATH $ENV{NMPI_ROOT}/include CACHE FILEPATH "")
set(MPI_CXX_LIBRARIES $ENV{NMPI_ROOT}/lib64/ve/libmpi++.a CACHE FILEPATH "")

#set(MPI_Fortran_COMPILER $ENV{NMPI_ROOT}/bin/mpif90 CACHE FILEPATH "")
set(MPI_Fortran_COMPILER /home/P1401/bin/mpif90 CACHE FILEPATH "")
set(MPI_Fortran_INCLUDE_PATH "$ENV{NMPI_ROOT}/lib/ve/module;$ENV{NMPI_ROOT}/include" CACHE FILEPATH "")
set(MPI_Fortran_LIBRARIES $ENV{NMPI_ROOT}/lib64/ve/libmpi.a CACHE FILEPATH "")
set(MPI_Fortran_COMPILE_FLAGS "-D_MPIPP_INCLUDE" CACHE STRING "")

#
# On the glibc enviroenment do not undefine __GNUC__!
#
#SET(CMAKE_C_FLAGS   "-U__GNUC__ -U__GNUC_MINOR__" CACHE STRING "" FORCE)
#SET(CMAKE_CXX_FLAGS "-U__GNUC__ -U__GNUC_MINOR__" CACHE STRING "" FORCE)
#SET(CMAKE_C_FLAGS   "-U__GNUC_MINOR__" CACHE STRING "" FORCE)
#SET(CMAKE_CXX_FLAGS "-U__GNUC_MINOR__" CACHE STRING "" FORCE)
SET(CMAKE_C_FLAGS "-U__GNUC_MINOR__ -O2 -g -lcblas -lblas_openmp -fopenmp"  CACHE STRING "")
SET(CMAKE_CXX_FLAGS "-U__GNUC_MINOR__ -O2 -g -lcblas -lblas_openmp -fopenmp" CACHE STRING "")
SET(CMAKE_Fortran_FLAGS "-O2 -g -llapack -lblas_openmp -fopenmp" CACHE STRING "")
SET(CMAKE_EXE_LINKER_FLAGS "-lrt" CACHE STRING "")

set(CMAKE_CROSSCOMPILING_EMULATOR "/opt/nec/ve/bin/ve_exec" CACHE FILEPATH "Command to execute VE binaries")

#set(BLAS_LIBRARIES $ENV{NLC_ROOT}/lib/libblas_sequential.a CACHE FILEPATH "BLAS library")
#set(LAPACK_LIBRARIES $ENV{NLC_ROOT}/lib/liblapack.a CACHE FILEPATH "LAPACK library")
#SET(BLAS_LIBRARIES /opt/nec/ve/nlc/2.0.0/lib/libblas_openmp.so)
SET(BLAS_LIBRARIES /opt/nec/ve/nlc/2.0.0/lib/libblas_sequential.so)
SET(LAPACK_LIBRARIES /opt/nec/ve/nlc/2.0.0/lib/liblapack.so)

set(CMAKE_C_STANDARD_COMPUTED_DEFAULT 11)
set(CMAKE_CXX_STANDARD_COMPUTED_DEFAULT 14)
set(CMAKE_Fortran_STANDARD_COMPUTED_DEFAULT 08)

# override CMake's compiler feature detection which fails for the NEC compiler, as it identifies as GNU 6.3, but does not support "-std=c++1z"
# minimal working setting here would be just setting: 
# set(CMAKE_CXX_COMPILE_FEATURES "cxx_auto_type;cxx_range_for;cxx_variadic_templates")
# a more complete list would be
set(CMAKE_CXX_COMPILE_FEATURES "cxx_template_template_parameters;cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX98_COMPILE_FEATURES "cxx_template_template_parameters")
set(CMAKE_CXX11_COMPILE_FEATURES "cxx_alias_templates;cxx_alignas;cxx_alignof;cxx_attributes;cxx_auto_type;cxx_constexpr;cxx_decltype;cxx_decltype_incomplete_return_types;cxx_default_function_template_args;cxx_defaulted_functions;cxx_defaulted_move_initializers;cxx_delegating_constructors;cxx_deleted_functions;cxx_enum_forward_declarations;cxx_explicit_conversions;cxx_extended_friend_declarations;cxx_extern_templates;cxx_final;cxx_func_identifier;cxx_generalized_initializers;cxx_inheriting_constructors;cxx_inline_namespaces;cxx_lambdas;cxx_local_type_template_args;cxx_long_long_type;cxx_noexcept;cxx_nonstatic_member_init;cxx_nullptr;cxx_override;cxx_range_for;cxx_raw_string_literals;cxx_reference_qualified_functions;cxx_right_angle_brackets;cxx_rvalue_references;cxx_sizeof_member;cxx_static_assert;cxx_strong_enums;cxx_thread_local;cxx_trailing_return_types;cxx_unicode_literals;cxx_uniform_initialization;cxx_unrestricted_unions;cxx_user_literals;cxx_variadic_macros;cxx_variadic_templates")
set(CMAKE_CXX14_COMPILE_FEATURES "${CMAKE_CXX_COMPILE_FEATURES}")
set(CMAKE_CXX17_COMPILE_FEATURES "")
set(CMAKE_CXX20_COMPILE_FEATURES "")

set(CMAKE_CXX17_STANDARD_COMPILE_OPTION "-std=c++17")
set(CMAKE_CXX17_EXTENSION_COMPILE_OPTION "-std=gnu++17")



#
# Local Fixes for building elmer
#
#SET(CMAKE_SYSTEM_NAME Linux)
#SET(CMAKE_C_COMPILER /opt/nec/ve/bin/ncc)
#SET(CMAKE_CXX_COMPILER /opt/nec/ve/bin/nc++)
#SET(CMAKE_Fortran_COMPILER /opt/nec/ve/bin/nfort)
#SET(CMAKE_Fortran_COMPILER /opt/nec/ve/bin/nfort-2.2.2)

#SET(CMAKE_C_FLAGS "-g -O2")
#SET(CMAKE_CXX_FLAGS "-g -O2")
#SET(CMAKE_Fortran_FLAGS "-g -O2")
#SET(CMAKE_C_FLAGS "-O2 -g -lcblas -lblas_openmp -fopenmp"  CACHE STRING "")
#SET(CMAKE_CXX_FLAGS "-O2 -g -lcblas -lblas_openmp -fopenmp" CACHE STRING "")
#SET(CMAKE_Fortran_FLAGS "-O2 -g -llapack -lblas_openmp -fopenmp" CACHE STRING "")
#SET(CMAKE_C_FLAGS "-O0 -g -lcblas -lblas_openmp -fopenmp"  CACHE STRING "")
#SET(CMAKE_CXX_FLAGS "-O0 -g -lcblas -lblas_openmp -fopenmp" CACHE STRING "")
#SET(CMAKE_Fortran_FLAGS "-O0 -g -llapack -lblas_openmp -fopenmp" CACHE STRING "")

#SET(CMAKE_Fortran_MODDIR_FLAG "-module")

#SET(CMAKE_C_COMPILE_OPTIONS_PIC "-fPIC")
#SET(CMAKE_C_COMPILE_OPTIONS_PIE "-fPIC")
#SET(CMAKE_CXX_COMPILE_OPTIONS_PIC "-fPIC")
#SET(CMAKE_CXX_COMPILE_OPTIONS_PIE "-fPIC")
#SET(CMAKE_Fortran_COMPILE_OPTIONS_PIC "-fPIC")
#SET(CMAKE_Fortran_COMPILE_OPTIONS_PIE "-fPIC")
#
#SET(CMAKE_SHARED_LIBRARY_C_FLAGS "-fPIC")
#SET(CMAKE_SHARED_LIBRARY_CXX_FLAGS "-fPIC")
#SET(CMAKE_SHARED_LIBRARY_Fortran_FLAGS "-fPIC")

#set(MPI_C_COMPILER /home/P1401/bin/mpincc CACHE FILEPATH "")
#set(MPI_CXX_COMPILER /home/P1401/bin/mpinc++ CACHE FILEPATH "")
#set(MPI_Fortran_COMPILER /home/P1401/bin/mpif90 CACHE FILEPATH "")

SET(BLAS_LIBRARIES /opt/nec/ve/nlc/2.0.0/lib/libblas_openmp.so)
SET(LAPACK_LIBRARIES /opt/nec/ve/nlc/2.0.0/lib/liblapack.so)
