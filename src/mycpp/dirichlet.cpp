
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
        Amat.setZero();
        mata1 ALmat(d-1, d-1);
        ALmat.diagonal() = beta.block(0,0,d-1,1);
        //the upper triangle has (d-1)-1 rows, the rows have 1 - ()d-2) elements. Arithmetic series:
        //(d-2) * (2 + d - 3) / 2 = (d-2) * (d-1)/2
        size_t upptrisize = (d-2) * (d-1)/2;
        ALmat.triangularView<Eigen::Upper>() = beta.block(d-1, 0, upptrisize, 1);
        ALmat.triangularView<Eigen::Lower>() = beta.block(d-1, 0, upptrisize, 1);
        Amat.block(0,0,d-1,d-1) = ALmat;
        veca1 bvec(d);
        bvec.setZero();
        bvec.block(0,0,d-1,1) = beta.block(d-1 + upptrisize, 0, d-1, 1);

        a1type y(0.);  // initialize summation
        y += u.transpose()*Amat*u;
        y += b.transpose()*u;
        //dirichlet component
        for(size_t i = 0; i < d; i++)
        {   y   += beta[i] * log(u[i]);
        }
        return y;
    }


}

