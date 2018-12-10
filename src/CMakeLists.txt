cmake_minimum_required(VERSION 3.6)
project(SeqArc)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -Wno-unused-function -O2 -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS")

set(LIBBWA_SOURCE_FILES
        bwa/utils.c
        bwa/kthread.c
        bwa/kstring.c
        bwa/ksw.c
        bwa/bwt.c
        bwa/bntseq.c
        bwa/bwa.c
        bwa/bwamem.c
        bwa/bwamem_pair.c
        bwa/bwamem_extra.c
        bwa/malloc_wrap.c
        bwa/QSufSort.c
        bwa/bwt_gen.c
        bwa/rope.c
        bwa/rle.c
        bwa/is.c
        bwa/bwtindex.c
        bwa/bamlite.c
        bwa/bwape.c
        bwa/bwase.c
        bwa/bwaseqio.c
        bwa/bwashm.c
        bwa/bwt_lite.c
        bwa/bwtaln.c
        bwa/bwtgap.c
        bwa/bwtsw2_aux.c
        bwa/bwtsw2_chain.c
        bwa/bwtsw2_core.c
        bwa/bwtsw2_main.c
        bwa/bwtsw2_pair.c
        bwa/fastmap.c
        bwa/kopen.c
        bwa/maxk.c
        bwa/pemerge.c)

set(MAIN_SOURCE_FILES
        kseq.hpp
        fqzcomp.cpp
        sfh.cpp
        main.cpp)

find_package(ZLIB REQUIRED)
if (ZLIB_FOUND)
    add_definitions(-DHAVE_ZLIB)
    include_directories(${ZLIB_INCLUDE_DIRS})
endif()

add_library(bwa ${LIBBWA_SOURCE_FILES})
add_executable(SeqArc ${MAIN_SOURCE_FILES})
if (CMAKE_SYSTEM_NAME MATCHES "Linux")
    target_link_libraries(SeqArc m pthread bwa ${ZLIB_LIBRARIES} rt)
else()
    target_link_libraries(SeqArc m pthread bwa ${ZLIB_LIBRARIES})
endif()
