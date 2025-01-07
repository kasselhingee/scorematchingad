Fixed the two issues note on CRAN:

+ finding cmake on mac.
+ test that was sensitive to the operating system.

Package passes:

+ R CMD check on windows, mac, linux (https://github.com/kasselhingee/scorematchingad/actions/runs/12644603482)
+ Rhub check with address sanitisation (https://github.com/kasselhingee/scorematchingad/actions/runs/12643534072)
+ Rhub check with valgrind finds only 'blocks are possibly lost' errors, which seems likely related to the timing of R's garbage collection (https://github.com/kasselhingee/scorematchingad/actions/runs/12626828775).

This version of the package now enables taping custom functions that works on all platforms. Which required a major change to the way tapes were accessed by R.

