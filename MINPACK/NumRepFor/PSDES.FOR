      SUBROUTINE psdes(lword,irword)
      INTEGER irword,lword,NITER
      PARAMETER (NITER=4)
      INTEGER i,ia,ib,iswap,itmph,itmpl,c1(4),c2(4)
      SAVE c1,c2
      DATA c1 /16#BAA96887,16#1E17D32C,16#03BCDC3C,16#0F33D1B2/, c2 
     */16#4B0F3B58,16#E874F0C3,16#6955C5A6, 16#55A7CA46/
      do 11 i=1,NITER
        iswap=irword
        ia=ieor(irword,c1(i))
        itmpl=iand(ia,65535)
        itmph=iand(ishft(ia,-16),65535)
        ib=itmpl**2+not(itmph**2)
        ia=ior(ishft(ib,16),iand(ishft(ib,-16),65535))
        irword=ieor(lword,ieor(c2(i),ia)+itmpl*itmph)
        lword=iswap
11    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j.
