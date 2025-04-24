      PROGRAM xsvdfit
C     driver for routine svdfit
      INTEGER NPOL,NPT
      doublePrecision SPREAD
      PARAMETER(NPT=100,SPREAD=0.02,NPOL=5)
      INTEGER i,idum,mp,np
      doublePrecision chisq,gasdev,stdErr(NPOL)
      doublePrecision x(NPT),y(NPT),sig(NPT),a(NPOL),cvm(NPOL,NPOL)
      doublePrecision u(NPT,NPOL),v(NPOL,NPOL),w(NPOL)!working space
      EXTERNAL fpoly,fleg,PolyFit
C     polynomial fit
      idum=-911
      mp=NPT
      np=NPOL
      do 11 i=1,NPT
        x(i)=0.02*i
        y(i)=1.0+x(i)*(2.0+x(i)*(3.0+x(i)*(4.0+x(i)*5.0)))
        y(i)=y(i)*(1.0+SPREAD*gasdev(idum))
        sig(i)=y(i)*SPREAD
11    continue

      call svdfit(x,y,sig,NPT,a,NPOL,u,v,w,mp,np,chisq,fpoly)
      call svdvar(v,NPOL,np,w,cvm,NPOL)
      write(*,*) 'Polynomial fit:'
      do 12 i=1,NPOL
        write(*,'(1x,f12.6,a,f10.6)') a(i),'  +-',sqrt(cvm(i,i))
12    continue
      write(*,'(1x,a,f12.6/)') 'Chi-squared',chisq
      call svdfit(x,y,sig,NPT,a,NPOL,u,v,w,mp,np,chisq,fleg)
      call svdvar(v,NPOL,np,w,cvm,NPOL)
      write(*,*) 'Legendre polynomial fit'
      do 13 i=1,NPOL
        write(*,'(1x,f12.6,a,f10.6)') a(i),'  +-',sqrt(cvm(i,i))
13    continue
      write(*,'(1x,a,f12.6/)') 'Chi-squared',chisq

      call SvdFitEz(PolyFit,x,y,mp,np,a,stdErr)
      write(*,*) 'Polynomial fit'
      do i=1,NPOL
        write(*,'(1x,f12.6,a,f10.6)') a(i),'  +-',stdErr(i)
	enddo
	stop
      END

      FUNCTION gasdev(idum)
      INTEGER idum
      doublePrecision gasdev
CU    USES ran2
      INTEGER iset
      doublePrecision fac,gset,rsq,v1,v2,ran2
      SAVE iset,gset
      DATA iset/0/
      if (iset.eq.0) then
1       v1=2.*ran2(idum)-1.
        v2=2.*ran2(idum)-1.
        rsq=v1**2+v2**2
        if(rsq.ge.1..or.rsq.eq.0.)goto 1
        fac=sqrt(-2.*log(rsq)/rsq)
        gset=v1*fac
        gasdev=v2*fac
        iset=1
      else
        gasdev=gset
        iset=0
      endif
      return
      END

      SUBROUTINE fpoly(x,p,np)
      INTEGER np
      doublePrecision x,p(np)
      INTEGER j
      p(1)=1.
      do 11 j=2,np
        p(j)=p(j-1)*x
11    continue
      return
      END

      SUBROUTINE fleg(x,pl,nl)
      INTEGER nl
      doublePrecision x,pl(nl)
      INTEGER j
      doublePrecision d,f1,f2,twox
      pl(1)=1.
      pl(2)=x
      if(nl.gt.2) then
        twox=2.*x
        f2=x
        d=1.
        do 11 j=3,nl
          f1=d
          f2=f2+twox
          d=d+1.
          pl(j)=(f2*pl(j-1)-f1*pl(j-2))/d
11      continue
      endif
      return
      END

      FUNCTION ran2(idum)
      INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
      doublePrecision ran2,AM,EPS,RNMX
      PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,IMM1=IM1-1,
     *IA1=40014,IA2=40692,IQ1=53668,IQ2=52774,IR1=12211,IR2=3791,
     *NTAB=32,NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
      INTEGER idum2,j,k,iv(NTAB),iy
      SAVE iv,iy,idum2
      DATA idum2/123456789/, iv/NTAB*0/, iy/0/
      if (idum.le.0) then
        idum=max(-idum,1)
        idum2=idum
        do 11 j=NTAB+8,1,-1

          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1
          if (idum.lt.0) idum=idum+IM1
          if (j.le.NTAB) iv(j)=idum
11      continue
        iy=iv(1)
      endif
      k=idum/IQ1
      idum=IA1*(idum-k*IQ1)-k*IR1
      if (idum.lt.0) idum=idum+IM1
      k=idum2/IQ2
      idum2=IA2*(idum2-k*IQ2)-k*IR2
      if (idum2.lt.0) idum2=idum2+IM2
      j=1+iy/NDIV
      iy=iv(j)-idum2
      iv(j)=idum
      if(iy.lt.1)iy=iy+IMM1
      ran2=min(AM*iy,RNMX)
      return
      END