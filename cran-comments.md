Dear CRAN Team,

Please find attached a revised R package for Hyvarinen score matching estimators by automatic differentiation. I have made the following changes in response to Konstanze Lauseker's review:

+ "In the description of the DESCRIPTION file ... please write year in parentheses." Thank you for picking this up. We have now followed the Author-Year style more closely with years, including years for the two references you asked about. For references without doi or similar we have switched to the [https:...]<https:...>.

+ "Please always write package names, software names and API (application programming interface) names in single quotes in title and description." I have written 'CppAD'.

+ "It seems like you have too many spaces in your description field. Probably because linebreaks count as spaces too. Please remove unecassary ones." There were a few extra spaces that I have now removed and the field contains no linebreaks. Could it be that you are referring to a rendering of the description field with fixed-width (these tend to include blank space due to the formatting of URLs)?

+ "\dontrun{} should only be used if the example really cannot be executed (e.g. because of missing additional software, missing API keys, ...) by the user. That's why wrapping examples in \dontrun{} adds the comment ("# Not run:") as a warning for the user. Please replace \dontrun with \donttest where possible. Please unwrap the examples if they are executable in < 5 sec, or replace dontrun{} with \donttest{}". The examples for customll() require additional software in the Enhances field of the DESCRIPTION - they can only be run if RcppEigen and RcppXPtrUtils is installed, and I've now added a comment as such to the help. The other use of dontrun{} has been changed to donttest{}.

---

The following are the comments I included in my first submission:

0. The RcppEigen package is listed in the Enhances section of DESCRIPTION because it is needed for calling Rcpp::cppFunction() (via RcppXPtrUtils::cppXPtr) with an RcppEigen dependency.

1. The feature involving customll() appears to only work on linux + gcc combinations. I have a provided a function for users to test if the feature works, and I hope to expand the functionality in the future.

2. Although the source tar ball is only 700KB, the installed size of the package is large due to incorporating the full CppAD source code. I have looked into avoiding this but opted otherwise because:
  + it will be very rare for CppAD to be installed on users computers, and it looks hard to check
  + the CppAD version inside the TMB R package is unsuitable for my purposes: the CppAD version is from 2015 (which misses at least one crucial feature I'm using) and contains modifications bespoke to TMB.


