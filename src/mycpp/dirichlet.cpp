
namespace { // begin the empty namespace

    a1type ll_dirichlet(const veca1 &beta,
	       const veca1 &u)
    {   size_t d  = u.size();
        a1type y(0.);  // initialize summation
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }

    a1type ll_ppi(const veca1 &beta, const veca1 & u){
        size_t d  = u.size();
        //assume the parameter vector beta is encoded as:
        //c(diag(ALs), ALs[upper.tri(ALs)], bL, betafordirichlet)
        mata1 Amat(d, d);

        a1type y(0.);  // initialize summation
        //dirichlet component
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }


}

