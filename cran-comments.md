Issue reported by the CRAN package check fixed.
Error was created by a lazily-written unit test that used `expect_error()` around a `expect_equal()` to check that two vectors are different.
I've corrected with a much more direct route.

Package passes on rhub checks of linux, macos, windows, with combinations of oldrel, release and devel versions of R.


