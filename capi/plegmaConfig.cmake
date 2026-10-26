# CMake package for plegma's C/C++ API, from a plegma checkout built with `pixi run build-capi`:
#
#     find_package(plegma CONFIG REQUIRED PATHS /path/to/plegma/capi)
#     target_link_libraries(mytool PRIVATE plegma::plegma)
#
# plegma::plegma puts capi/ (plegma.h, plegma.hpp) on the include path and links build/libplegma.so with run paths
# to it and to the Mojo runtime libraries in plegma's pixi environment.
get_filename_component(_plegma_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(plegma_VERSION 0.2.0)
set(_plegma_library "${_plegma_root}/build/libplegma.so")
if(NOT EXISTS "${_plegma_library}")
    message(FATAL_ERROR "plegma: ${_plegma_library} not found; build it with `pixi run build-capi` in ${_plegma_root}")
endif()
if(NOT TARGET plegma::plegma)
    add_library(plegma::plegma SHARED IMPORTED)
    set_target_properties(plegma::plegma PROPERTIES
        IMPORTED_LOCATION "${_plegma_library}"
        INTERFACE_INCLUDE_DIRECTORIES "${_plegma_root}/capi"
        INTERFACE_LINK_OPTIONS "-Wl,-rpath,${_plegma_root}/build;-Wl,-rpath,${_plegma_root}/.pixi/envs/default/lib")
endif()
