cmake_minimum_required(VERSION 3.26)
project(LongBCH)

set(CMAKE_CXX_STANDARD 17)
set(GMP_INCLUDE_DIR "C:/tools/msys64/mingw64/include")
set(GMP_LIBRARY "C:/tools/msys64/mingw64/lib/libgmp.a")
set(GMPXX_LIBRARY "C:/tools/msys64/mingw64/lib/libgmpxx.a")
include_directories(${GMP_INCLUDE_DIR})
add_executable(LongBCH main.cpp
        Polinomial.cpp
        Polynomial.h
        GFq.cpp
        GFq.h
        Algorthms.h
        Algorthms.cpp
        testing.cpp
        testing.h
        Utils.h
        Utils.cpp
        Tree.h
        Tree.cpp
        BCH.h
        BCH.cpp
)
target_link_libraries(LongBCH PRIVATE ${GMPXX_LIBRARY} ${GMP_LIBRARY})