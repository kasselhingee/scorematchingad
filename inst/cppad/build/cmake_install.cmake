# Install script for directory: /home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "minsizerel")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/usr/local/include/cppad/")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/usr/local/include/cppad" TYPE DIRECTORY FILES "/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/include/cppad/" FILES_MATCHING REGEX "/[^/]*\\.hpp$" REGEX "/xrst$" EXCLUDE)
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/include/cppad/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/pkgconfig/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/cppad_lib/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/val_graph/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/example/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/introduction/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/test_more/cmake_install.cmake")
  include("/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/speed/cmake_install.cmake")

endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/kassel/Documents/professional/ANU_Compositional/scorecompdir/versioned/inst/cppad/build/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
