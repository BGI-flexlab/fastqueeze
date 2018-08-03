cmake_minimum_required(VERSION 3.6)
project(bwaTest)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_C_FLAGS "-g -Wall -Wno-unused-function -O2 -DHAVE_PTHREAD -DUSE_MALLOC_WRAPPERS")

set(LIBBWA_SOURCE_FILES
        utils.c
        kthread.c
        kstring.c
        ksw.c
        bwase.c
        bwt.c
        bntseq.c
        bwa.c
        bwamem.c
        bwamem_pair.c
        bwamem_extra.c
        malloc_wrap.c
        QSufSort.c
        bwt_gen.c
        rope.c
        rle.c
        is.c
        bwtindex.c)

set(MAIN_SOURCE_FILES
        bamlite.c
        bwaseqio.c
        bwashm.c
        bwt_lite.c
        bwtaln.c
        bwtgap.c
        bwtindex.c
        bwtsw2_aux.c
        bwtsw2_chain.c
        bwtsw2_core.c
        bwtsw2_main.c
        bwtsw2_pair.c
        fastsam.c)

add_library(bwa ${LIBBWA_SOURCE_FILES})
add_executable(fastsam ${MAIN_SOURCE_FILES})
target_link_libraries(fastsam m z pthread bwa)