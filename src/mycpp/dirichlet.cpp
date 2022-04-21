
namespace { // begin the empty namespace

    a1type ll(const veca1 &beta,
	       const veca1 &u)
    {   size_t n  = beta.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < n; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }


}

