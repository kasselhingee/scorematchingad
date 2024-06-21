Issue with archived version on CRAN should be resolved. Previous attempt didn't fix it as I don't have access to an M1mac for checking. Thanks to Uwe, realised test failure on M1mac due to eigen decomposition sign identifiability. A change so that only the non-decomposed matrix is tested for equality should make tests on M1mac pass.
My sincere apologies for being slow to fix this problem.

Also error with my submission on 13 Jun seems resolved and passes checks on win_devel, and passes rhub_check() on linux, macos, windows, ubuntu-clang and more.


