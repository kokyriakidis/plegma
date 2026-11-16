# CMake package for plegma's C/C++ API, from a plegma checkout built with `pixi run build-capi`:
#
#     find_package(plegma CONFIG REQUIRED PATHS /path/to/plegma/capi)
#     target_link_libraries(mytool PRIVATE plegma::plegma)
#
# plegma::plegma puts capi/ (plegma.h, plegma.hpp) on the include path and links build/libplegma.so with a run path
# to build/, which also holds the Mojo runtime libraries it needs (copied there by build-capi).
#
# plegma::static links build/libplegma.a instead (`pixi run build-static`): no Mojo runtime at run time, usable in a
# static executable. Defined when build/libplegma.a exists.
get_filename_component(_plegma_root "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
set(plegma_VERSION 0.4.0)
set(_plegma_library "${_plegma_root}/build/libplegma.so")
set(_plegma_static_library "${_plegma_root}/build/libplegma.a")
if(NOT EXISTS "${_plegma_library}" AND NOT EXISTS "${_plegma_static_library}")
    message(FATAL_ERROR "plegma: neither ${_plegma_library} nor ${_plegma_static_library} found; build one with "
        "`pixi run build-capi` or `pixi run build-static` in ${_plegma_root}")
endif()
if(EXISTS "${_plegma_static_library}" AND NOT TARGET plegma::static)
    find_package(Threads REQUIRED)
    add_library(plegma::static STATIC IMPORTED)
    set_target_properties(plegma::static PROPERTIES
        IMPORTED_LOCATION "${_plegma_static_library}"
        INTERFACE_INCLUDE_DIRECTORIES "${_plegma_root}/capi"
        INTERFACE_LINK_LIBRARIES "Threads::Threads;m")
endif()
if(EXISTS "${_plegma_library}" AND NOT TARGET plegma::plegma)
    add_library(plegma::plegma SHARED IMPORTED)
    set_target_properties(plegma::plegma PROPERTIES
        IMPORTED_LOCATION "${_plegma_library}"
        INTERFACE_INCLUDE_DIRECTORIES "${_plegma_root}/capi"
        INTERFACE_LINK_OPTIONS "-Wl,-rpath,${_plegma_root}/build")
endif()
