      FUNCTION ellf(phi,ak)
      REAL ellf,ak,phi
CU    USES rf
      REAL s,rf
      s=sin(phi)
      ellf=s*rf(cos(phi)**2,(1.-s*ak)*(1.+s*ak),1.)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j.
