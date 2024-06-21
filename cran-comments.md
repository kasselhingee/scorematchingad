Error with my submission on 13 Jun seems resolved and passes checks with on win_devel, and passes rhub_check() on all the platforms I've tried (linux, macos, windows, ubuntu-clang etc.

Issue with archived version on CRAN should also be resolved: 
Thanks to Uwe, realised test failure on M1mac due to eigen decomposition sign identifiability. A change so that only the non-decomposed matrix is tested for equality should make it pass.

My sincere apologies for being slow to fix this problem this time.

