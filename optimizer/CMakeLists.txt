cmake_minimum_required(VERSION 3.22)

project(LED
  VERSION 0.0.0
  LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)

set(BUILD_SHARED_LIBS ON)
set(PN ${PROJECT_NAME})

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)
include(CheckIPOSupported)

check_ipo_supported(RESULT IPO_AVAILABLE)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
include(CompilerFlagsSets)
include(TargetCompileCheckAddFlags)

if (NOT ${BUILD_SHARED_LIBS})
  message(FATAL_ERROR "lib${PN} requires BUILD_SHARED_LIBS ON because only shared libraries can be dynamically loaded")
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
find_package(pybind11 REQUIRED)

add_subdirectory(external)

set(SOURCES source/PyAPI.cpp)

pybind11_add_module(${PN}_py MODULE ${SOURCES})
add_dependencies(${PN}_py libtagarray)
target_include_directories(${PN}_py PRIVATE ${PROJECT_SOURCE_DIR}/include)
# libtagarray
target_link_libraries(${PN}_py PUBLIC tagarray)
target_link_directories(${PN}_py PRIVATE ${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-build)
target_include_directories(${PN}_py PRIVATE ${CMAKE_BINARY_DIR}/external/libtagarray/src/libtagarray-build/include
                                            ${CMAKE_SOURCE_DIR}/external/libtagarray/src/libtagarray/include)
# target features
target_compile_features(${PN}_py PRIVATE cxx_std_17)
target_compile_check_add_flags(${PN}_py PRIVATE CXX ${${CMAKE_CXX_COMPILER_ID}_CXX_WARNINGS})
if (IPO_AVAILABLE)
  set_target_properties(${PN}_py PROPERTIES INTERPROCEDURAL_OPTIMIZATION TRUE)
endif()
set_target_properties(${PN}_py PROPERTIES OUTPUT_NAME ${PN})

install(TARGETS ${PN}_py
        EXPORT ${PN}Targets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
        FILE_SET HEADERS DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT ${PN}Targets
    FILE ${PN}Targets.cmake
    NAMESPACE ${PN}::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PN})

write_basic_package_version_file(
    "${PN}ConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion)

install(FILES "${CMAKE_CURRENT_BINARY_DIR}/${PN}ConfigVersion.cmake" DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PN})
