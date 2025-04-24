!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      SUBROUTINE PolyFit(x,p,np)
	!polynomial regression by Numerical Recipes. np = # parameters(?) = nPolyOrder+1 if intercept != 0.
      INTEGER np
      doublePrecision x,p(np)
      INTEGER j
      p(1)=1.
      do j=2,np
		!p(j)=x**j !this is less efficient but gives the same thing
        p(j)=p(j-1)*x
	enddo
      return
      END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      SUBROUTINE POLYFIT0(x,p,np)
	!Polynomial regression with zero intercept
      INTEGER np
      doublePrecision x,p(np)
      INTEGER j
      p(1)=x
      do j=2,np
		!p(j)=x**j !this is less efficient but gives the same thing
        p(j)=p(j-1)*x
	enddo
      return
      END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
	SUBROUTINE A2Basis(eta,p,np)
	!Polynomial in numerator with (1+500*eta^4) in denominator.  For TPT A2.
	INTEGER np
	doublePrecision eta,p(np)
	INTEGER j
	eta4=eta*eta*eta*eta
	eta5=eta*eta4
	p(1)=eta/(1+500*eta4)
	do j=2,np
		p(j)=p(j-1)*eta !/(1+500*eta4)	Note: divide only to first power.
	enddo
	return
	END

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
      subroutine SvdFitEz(funcs,x,y,nData,nParms,parm,stdErr,chisq)
C     driver for routine svdfit
c     to fit a polynomial to nData points of y-data vs. x-data.
c     the polynomial order depends on the funcs specification
c		order=nParms for PolyFit0 (zero intercept assumed).
c		order=nParms-1 for PolyFit.
      INTEGER NPOL,NPT
      PARAMETER(NPT=1000,NPOL=50)
      doublePrecision chisq,cvm(NPOL,NPOL),sig(NPT)
      doublePrecision x(nData),y(nData),parm(nParms),stdErr(nParms)
      doublePrecision u(NPT,NPOL),v(NPOL,NPOL),w(NPOL)!working space
      EXTERNAL funcs
	do iPt=1,nData
		sig(iPt)=1.0 ! gives excel result for unweighted least squares.
c        sig(iPt)=y(iPt)*0.02
c        sig(iPt)=y(iPt)*0.05 ! gives %err MINimum	and err estimates based on 5% err
	enddo
	if(nParms.gt.nData-2)pause 'SvdFit warning: nParms>nData-2'
	if(nParms.gt.nData)then
		pause 'SvdFit Error: nParms>nData'
		return
	endif
      call svdfit(x,y,sig,nData,parm,nParms,u,v,w,NPT,NPOL,chisq,funcs)
      call svdvar(v,nParms,NPOL,w,cvm,NPOL)
	do iParm=1,nParms
		stdErr(iParm)=0
		if(nData.ne.nParms)then
		stdErr(iParm)=SQRT( cvm(iParm,iParm)*chisq/(nData-nParms) )
		endif
	enddo
	return
	end

      subroutine SVDFITEZ0(x,y,nData,nParms,parm,stdErr,chisq)
      INTEGER NPOL,NPT,nParms,nData,iPt,iParm
      PARAMETER(NPT=1000,NPOL=50)
      doublePrecision chisq,cvm(NPOL,NPOL),sig(NPT)
      doublePrecision x(nData),y(nData),parm(nParms),stdErr(nParms)
      doublePrecision u(NPT,NPOL),v(NPOL,NPOL),w(NPOL)!working space
      EXTERNAL POLYFIT0
	!the following are valuable for mixed language programMINg with C++. note:must comment out for debugging fortran.
	!DEC$ ATTRIBUTES  VALUE    :: nParms
	!DEC$ ATTRIBUTES  VALUE    :: nData
	!DEC$ ATTRIBUTES REFERENCE :: chisq
C     driver for routine svdfit where PolyFit0 has been assumed as the fitting function.
c     to fit a polynomial to nData points of y-data vs. x-data.
c     the polynomial order depends on the funcs specification
c		order=nParms for PolyFit0 (zero intercept assumed).
	do iPt=1,nData
		sig(iPt)=1.0
c		sig(iPt)=1.0/y(iPt) ! gives excel result for unweighted least squares.
c        sig(iPt)=y(iPt)*0.02
c        sig(iPt)=y(iPt)*0.05 ! gives %err MINimum	and err estimates based on 5% err
	enddo
      call svdfit(x,y,sig,nData,parm,nParms,u,v,w,NPT,NPOL,chisq
     &          ,POLYFIT0)
      call svdvar(v,nParms,NPOL,w,cvm,NPOL)
	do iParm=1,nParms
		stdErr(iParm)=SQRT( cvm(iParm,iParm)*chisq/(nData-nParms) )
	enddo
	return
	end

      SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,funcs)
      INTEGER ma,mp,ndata,np,NMAX,MMAX
      DOUBLE PRECISION chisq,a(ma),sig(ndata),u(mp,np),v(np,np),w(np),
     *  x(ndata),y(ndata),TOL
      EXTERNAL funcs
      PARAMETER (NMAX=1000,MMAX=50,TOL=1.d-13)
CU    USES svbksb,svdcmp
      INTEGER i,j
      DOUBLE PRECISION sum,thresh,tmp,wMAX,afunc(MMAX),b(NMAX)
      do 12 i=1,ndata
        call funcs(x(i),afunc,ma)
        tmp=1.d0/sig(i)
        do 11 j=1,ma
          u(i,j)=afunc(j)*tmp
11      continue
        b(i)=y(i)*tmp
12    continue
      call svdcmp(u,ndata,ma,mp,np,w,v)

      wMAX=0.d0
      do 13 j=1,ma
        if(w(j).gt.wMAX)wMAX=w(j)
13    continue
      thresh=TOL*wMAX
      do 14 j=1,ma
        if(w(j).lt.thresh)w(j)=0.d0
14    continue
      call svbksb(u,w,v,ndata,ma,mp,np,b,a)
      chisq=0.d0
      do 16 i=1,ndata
        call funcs(x(i),afunc,ma)
        sum=0.d0
        do 15 j=1,ma
          sum=sum+a(j)*afunc(j)
15      continue
        chisq=chisq+((y(i)-sum)/sig(i))**2
16    continue
      return
      END

      SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION b(mp),u(mp,np),v(np,np),w(np),x(np)
      PARAMETER (NMAX=500)
      INTEGER i,j,jj
      DOUBLE PRECISION s,tmp(NMAX)
      do 12 j=1,n
        s=0.d0
        if(w(j).ne.0.d0)then
          do 11 i=1,m
            s=s+u(i,j)*b(i)
11        continue
          s=s/w(j)
        endif
        tmp(j)=s
12    continue
      do 14 j=1,n
        s=0.d0
        do 13 jj=1,n
          s=s+v(j,jj)*tmp(jj)
13      continue
        x(j)=s
14    continue
      return
      END

      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
      INTEGER m,mp,n,np,NMAX
      DOUBLE PRECISION a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=500)
CU    USES pythag
      INTEGER i,its,j,jj,k,l,nm
      DOUBLE PRECISION anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
      g=0.0d0
      scale=0.0d0
      anorm=0.0d0
      do 25 i=1,n
        l=i+1
        rv1(i)=scale*g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if(i.le.m)then
          do 11 k=i,m
            scale=scale+ABS(a(k,i))
11        continue
          if(scale.ne.0.0d0)then
            do 12 k=i,m
              a(k,i)=a(k,i)/scale

              s=s+a(k,i)*a(k,i)
12          continue
            f=a(i,i)
            g=-SIGN(SQRT(s),f)
            h=f*g-s
            a(i,i)=f-g
            do 15 j=l,n
              s=0.0d0
              do 13 k=i,m
                s=s+a(k,i)*a(k,j)
13            continue
              f=s/h
              do 14 k=i,m
                a(k,j)=a(k,j)+f*a(k,i)
14            continue
15          continue
            do 16 k=i,m
              a(k,i)=scale*a(k,i)
16          continue
          endif
        endif
        w(i)=scale *g
        g=0.0d0
        s=0.0d0
        scale=0.0d0
        if((i.le.m).and.(i.ne.n))then
          do 17 k=l,n

            scale=scale+ABS(a(i,k))
17        continue
          if(scale.ne.0.0d0)then
            do 18 k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
18          continue
            f=a(i,l)
            g=-SIGN(SQRT(s),f)
            h=f*g-s
            a(i,l)=f-g
            do 19 k=l,n
              rv1(k)=a(i,k)/h
19          continue
            do 23 j=l,m
              s=0.0d0
              do 21 k=l,n
                s=s+a(j,k)*a(i,k)
21            continue
              do 22 k=l,n
                a(j,k)=a(j,k)+s*rv1(k)
22            continue
23          continue
            do 24 k=l,n

              a(i,k)=scale*a(i,k)
24          continue
          endif
        endif
        anorm=MAX(anorm,(ABS(w(i))+ABS(rv1(i))))
25    continue
      do 32 i=n,1,-1
        if(i.lt.n)then
          if(g.ne.0.0d0)then
            do 26 j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
26          continue
            do 29 j=l,n
              s=0.0d0
              do 27 k=l,n
                s=s+a(i,k)*v(k,j)
27            continue
              do 28 k=l,n
                v(k,j)=v(k,j)+s*v(k,i)
28            continue

29          continue
          endif
          do 31 j=l,n
            v(i,j)=0.0d0
            v(j,i)=0.0d0
31        continue
        endif
        v(i,i)=1.0d0
        g=rv1(i)
        l=i
32    continue
      do 39 i=MIN(m,n),1,-1
        l=i+1
        g=w(i)
        do 33 j=l,n
          a(i,j)=0.0d0
33      continue
        if(g.ne.0.0d0)then
          g=1.0d0/g
          do 36 j=l,n
            s=0.0d0
            do 34 k=l,m
              s=s+a(k,i)*a(k,j)
34          continue
            f=(s/a(i,i))*g
            do 35 k=i,m

              a(k,j)=a(k,j)+f*a(k,i)
35          continue
36        continue
          do 37 j=i,m
            a(j,i)=a(j,i)*g
37        continue
        else
          do 38 j= i,m
            a(j,i)=0.0d0
38        continue
        endif
        a(i,i)=a(i,i)+1.0d0
39    continue
      do 49 k=n,1,-1
        do 48 its=1,30
          do 41 l=k,1,-1
            nm=l-1
            if((ABS(rv1(l))+anorm).eq.anorm)  goto 2
            if((ABS(w(nm))+anorm).eq.anorm)  goto 1
41        continue
1         c=0.0d0
          s=1.0d0
          do 43 i=l,k

            f=s*rv1(i)
            rv1(i)=c*rv1(i)
            if((ABS(f)+anorm).eq.anorm) goto 2
            g=w(i)
            h=pythag(f,g)
            w(i)=h
            h=1.0d0/h
            c= (g*h)
            s=-(f*h)
            do 42 j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
42          continue
43        continue
2         z=w(k)
          if(l.eq.k)then
            if(z.lt.0.0d0)then
              w(k)=-z

              do 44 j=1,n
                v(j,k)=-v(j,k)
44            continue
            endif
            goto 3
          endif
          if(its.eq.30) pause 'no convergence in svdcmp'
          x=w(l)
          nm=k-1
          y=w(nm)
          g=rv1(nm)
          h=rv1(k)
          f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0d0*h*y)
          g=pythag(f,1.0d0)
          f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
          c=1.0d0
          s=1.0d0
          do 47 j=l,nm
            i=j+1
            g=rv1(i)
            y=w(i)
            h=s*g
            g=c*g
            z=pythag(f,h)

            rv1(j)=z
            c=f/z
            s=h/z
            f= (x*c)+(g*s)
            g=-(x*s)+(g*c)
            h=y*s
            y=y*c
            do 45 jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
45          continue
            z=pythag(f,h)
            w(j)=z
            if(z.ne.0.0d0)then
              z=1.0d0/z
              c=f*z
              s=h*z
            endif
            f= (c*g)+(s*y)

            x=-(s*g)+(c*y)
            do 46 jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
46          continue
47        continue
          rv1(l)=0.0d0
          rv1(k)=f
          w(k)=x
48      continue
3       continue
49    continue
      return
      END

      SUBROUTINE svdvar(v,ma,np,w,cvm,ncvm)
      INTEGER ma,ncvm,np,MMAX
      DOUBLE PRECISION cvm(ncvm,ncvm),v(np,np),w(np)
      PARAMETER (MMAX=20)
      INTEGER i,j,k
      DOUBLE PRECISION sum,wti(MMAX)
      do 11 i=1,ma
        wti(i)=0.d0
        if(w(i).ne.0.d0) wti(i)=1.d0/(w(i)*w(i))
11    continue
      do 14 i=1,ma
        do 13 j=1,i
          sum=0.d0
          do 12 k=1,ma
            sum=sum+v(i,k)*v(j,k)*wti(k)
12        continue
          cvm(i,j)=sum
          cvm(j,i)=sum
13      continue
14    continue
      return
      END


      FUNCTION pythag(a,b)
      doublePrecision a,b,pythag
      doublePrecision ABSa,ABSb
      ABSa=ABS(a)
      ABSb=ABS(b)
      if(ABSa.gt.ABSb)then
        pythag=ABSa*SQRT(1.+(ABSb/ABSa)**2)
      else
        if(ABSb.eq.0.)then
          pythag=0.
        else
          pythag=ABSb*SQRT(1.+(ABSa/ABSb)**2)
        endif
      endif
      return
      END
