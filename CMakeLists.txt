cmake_minimum_required(VERSION 2.8)

project(Jacobi)
include_directories(/opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/include/)
add_executable(${PROJECT_NAME} "main.cpp" )
target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl//lib/intel64/libmkl_core.a )
#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_intel_ilp64.a )
#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl/lib/intel64/libmkl_gnu_thread.a )
#target_link_libraries(Jacobi /opt/intel/compilers_and_libraries_2017.4.196/linux/mkl//lib/intel64/libmkl_core.a  )
