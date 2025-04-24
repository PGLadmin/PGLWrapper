      FUNCTION beta(z,w)
      REAL beta,w,z
CU    USES gammln
      REAL gammln
      beta=exp(gammln(z)+gammln(w)-gammln(z+w))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j.
