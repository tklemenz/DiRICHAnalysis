cmake_minimum_required(VERSION 2.8.12 FATAL_ERROR)
project(DiRICHAnalysis)

set(CMAKE_CXX_STANDARD 11)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
set(CMAKE_MODULE_PATH $ENV{ROOTSYS}/etc/cmake)
find_package(ROOT REQUIRED COMPONENTS MathCore RIO Hist Tree Net)

set(CMAKE_BUILD_TYPE_SHARED_LIBS ON)

#find_package(fmt REQUIRED)

#set(Boost_USE_STATIC_LIBS OFF)
#set(Boost_USE_MULTITHREADED ON)
#set(Boost_USE_STATIC_RUNTIME OFF)
#find_package(Boost 1.65.0 COMPONENTS log)

add_definitions(${ROOT_CXX_FLAGS})

include_directories(${CMAKE_SOURCE_DIR} ${CMAKE_SOURCE_DIR}/src ${CMAKE_SOURCE_DIR}/include ${ROOT_INCLUDE_DIRS})# ${Boost_INCLUDE_DIRS})

set(SRCS
   src/Dummy.cxx
   )
set(SRCS
    src/Utility.cxx
   )
   
set(HEADERS
   include/Dummy.h
   )
set(HEADERS
   include/Utility.h
   )
   
ROOT_GENERATE_DICTIONARY(G__RMA ${HEADERS} LINKDEF include/LinkDef.h)
add_library(DiRICHAna SHARED ${SRCS} G__RMA)
target_link_libraries(DiRICHAna ${ROOT_LIBRARIES})
                                #fmt::fmt
                                #${Boost_LIBRARIES})

#install(TARGETS DiRICHAna ${INSTALL_OPTS})

# Create main programs using the library

set(EXE_SRCS
    macros/runDummy.C
    )

set(EXE_NAMES
    runDummy
    )

list(LENGTH EXE_SRCS count)
math(EXPR count "${count}-1")
foreach(i RANGE ${count})
  list(GET EXE_SRCS ${i} src)
  list(GET EXE_NAMES ${i} name)
  add_executable(${name} ${src})
  target_link_libraries(${name} DiRICHAna)
  set_target_properties(${name}
                      PROPERTIES
                      RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
endforeach()
