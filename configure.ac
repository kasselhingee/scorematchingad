
AC_INIT([scorematchingad], [0.0.62], [kassel.hingee@anu.edu.au])

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
CMAKE_PATH="$PATH:/Applications/CMake.app/Contents/bin/cmake"
AC_PATH_PROG([CMAKE], [cmake], [no], [$CMAKE_PATH])
if test "$CMAKE" = "no"; then
  AC_MSG_ERROR([CMake is required but was not found. Please install CMake.])
fi

# configure cppad
cd ../inst/cppad
$CMAKE . \
	-DCMAKE_BUILD_TYPE=minsizerel \
	-DBUILD_SHARED_LIBS:bool=OFF \
	-DCMAKE_POSITION_INDEPENDENT_CODE:bool=ON
# copy the generated cppad files into a different directory
"$(R_HOME)/bin$(R_ARCH_BIN)/Rscript" -e 'file.copy("./include/cppad", "../include/", recursive = TRUE)'

# add more CPPFLAGS via the Makevars file (Could they be added earlier than this?)
echo "PKG_CPPFLAGS=-I../inst/include -Wno-ignored-attributes" > src/Makevars

AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT
