cmake_minimum_required(VERSION 2.8)

project(Jacobi)
include_directories(/opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/include/)
aux_source_directory(. SRC_LIST)
add_executable(${PROJECT_NAME} ${SRC_LIST})
#add_executable(${PROJECT_NAME} "main.cpp" )
#ADD_LIBRARY(libmkl_intel_ilp64 STATIC /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_ilp64.a)
#ADD_LIBRARY(libmkl_gnu_thread STATIC /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_gnu_thread.a)
#ADD_LIBRARY(libmkl_core STATIC /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl//lib/intel64/libmkl_core.a)

TARGET_LINK_LIBRARIES(Jacobi -Wl,--start-group /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl//lib/intel64/libmkl_core.a -Wl,--end-group)

#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl//lib/intel64/libmkl_core.a )
#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_ilp64.a )
#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_gnu_thread.a )
#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl//lib/intel64/libmkl_core.a  )
