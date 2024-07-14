# CppAD: A Package for Differentiation of C++ Algorithms

*Note by Kassel Liam Hingee:*  This directory contains a modified version of CppAD 20240000.5 for use within an R package.  I have removed directories and files in the root directory that I didn't appear to need (xrst, xrst.toml, user_guide.xrst, .travis.yml, .readthedocs.yml, appveyor.yml, appendix/. batch_edit.sed, bug/ .circleci, .coin-or, .github, .gitignore) and much of bin/.  I have also replaced all use of std::cout with Rcpp::Rcout (similary for std::cerr) and added an include directive to RcppCommon for any file containing such changes.  The copyright for CppAD is held by Bradley M. Bell (Copyright (C) 2003-18). The CppAD source code is available at https://github.com/coin-or/CppAD.git .


## Documentation
[users guide](https://cppad.readthedocs.io/en/latest/user_guide.html)

## License
SPDX-License-Identifier: EPL-2.0 OR GPL-2.0-or-later

## Install

- The preferred method to test and
  [install](https://cppad.readthedocs.io/en/latest/Install.html)
  CppAD uses [cmake](https://cmake.org).

- A deprecated
  [autotools](https://cppad.readthedocs.io/en/latest/autotools.html)
  procedure can be used for this purpose, but it will eventually be removed.

- For any sub-directory *dir*,
  files of the form *dir*/`makefile.am` and *dir*/`makefile.in`
  are used to support the Autotools test and install procedure.
  In addition,
  the following files, in this directory, are also for this purpose:
  `compile`,
  `config.guess`,
  `config.sub`,
  `configure.ac`,
  `depcomp`,
  `install-sh`,
  `missing`.
