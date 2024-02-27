include(ExternalProject)

ExternalProject_Add(libtagarray
  PREFIX ${CMAKE_SOURCE_DIR}/external/libtagarray
  GIT_REPOSITORY https://github.com/foxtran/libtagarray.git
  GIT_TAG master
  BINARY_DIR ${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-build
  STAMP_DIR  ${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-stamp
  CMAKE_ARGS -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
             -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
             -DCMAKE_GENERATOR:INTERNAL=${CMAKE_GENERATOR}
             -DENABLE_FORTRAN=OFF
             -DENABLE_PYTHON=ON
             -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
             -DENABLE_TESTING=${ENABLE_TESTS}
             -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
             -DCMAKE_PROJECT_TOP_LEVEL_INCLUDES=${CMAKE_SOURCE_DIR}/external/libtagarray/src/libtagarray/cmake/conan_provider.cmake
             -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
  )
