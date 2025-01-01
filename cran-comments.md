Fixed issue with finding cmake on mac.
Also fixed a test that was sensitive to the operating system.
Passes Rhub and R CMD check on windows, mac, linux (https://github.com/kasselhingee/scorematchingad/actions/workflows/rhub.yaml and https://github.com/kasselhingee/scorematchingad/actions/workflows/R-CMD-check.yaml)

New feature of taping custom functions now works on all platforms. Which required a major change to the way tapes are accessed by R.

