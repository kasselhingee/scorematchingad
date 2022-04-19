
namespace { // begin the empty namespace

    a2type ll(const veca1 &beta,
	       const veca2 &u)
    {   size_t n  = beta.size();
        a2type y(0.);  // initialize summation
        for(size_t i = 0; i < n; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }


}

