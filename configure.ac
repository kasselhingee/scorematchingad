# Process this file with autoconf to produce a configure script.

AC_INIT([scorematchingad], [0.0.71], [kassel.hingee@anu.edu.au])

AC_PREREQ([2.71])

# From 'writing R extensions'
: ${R_HOME=`R RHOME`}
if test -z "${R_HOME}"; then
  echo "could not determine R_HOME"
  exit 1
fi
CC=`"${R_HOME}/bin/R" CMD config CC`
CFLAGS=`"${R_HOME}/bin/R" CMD config CFLAGS`

# Writing R Extensions says C++ compiler details and language switch to C++ is needed
CXX=`"${R_HOME}/bin/R" CMD config CXX`
if test -z "$CXX"; then
  AC_MSG_ERROR([No C++ compiler is available])
fi
CXXFLAGS=`"${R_HOME}/bin/R" CMD config CXXFLAGS`
CPPFLAGS=`"${R_HOME}/bin/R" CMD config CPPFLAGS`
AC_LANG(C++)

# find cmake looking in the usual macos location too
AC_PATH_PROG([CMAKE], [cmake], [no], [$PATH$PATH_SEPARATOR/Applications/CMake.app/Contents/bin])
if test "$CMAKE" = "no"; then
  AC_MSG_ERROR([CMake is required but was not found. Please install CMake.])
fi

# Set the path to the C++ library
CPPAD_PATH="`pwd`/inst/cppad"

# Run CMake with the appropriate flags
AC_MSG_NOTICE([Running CMake...])
mkdir -p $CPPAD_PATH/build
cd $CPPAD_PATH/build
$CMAKE .. -D CMAKE_BUILD_TYPE=minsizerel -DBUILD_SHARED_LIBS:bool=OFF -DCMAKE_POSITION_INDEPENDENT_CODE=bool=ON
if test $? -ne 0; then
  AC_MSG_ERROR([CMake configuration failed])
fi
cd -

# Copy the include files to the appropriate directory
"$R_HOME/bin/Rscript" -e 'file.copy("inst/cppad/include/cppad", "inst/include/", recursive = TRUE)'
# Copy temp_file.cpp to ./src
"$R_HOME/bin/Rscript" -e 'file.copy("inst/cppad/cppad_lib/temp_file.cpp", "./src/")'
# Cleanup build directory so that CRAN doesn't see any make files
rm -r -f ./inst/cppad/build

# Create the Makevars file with the necessary compilation flags
echo "PKG_CPPFLAGS=-I../inst/include" > src/Makevars

# Output files
AC_OUTPUT
