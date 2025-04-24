      FUNCTION bico(n,k)
      INTEGER k,n
      REAL bico
CU    USES factln
      REAL factln
      bico=nint(exp(factln(n)-factln(k)-factln(n-k)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j.
