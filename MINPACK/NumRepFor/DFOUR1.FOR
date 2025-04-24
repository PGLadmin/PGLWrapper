      SUBROUTINE four1(data,nn,isign)
      INTEGER isign,nn
      DOUBLE PRECISION data(2*nn)
      INTEGER i,istep,j,m,mmax,n
      DOUBLE PRECISION tempi,tempr
      DOUBLE PRECISION theta,wi,wpi,wpr,wr,wtemp
      n=2*nn
      j=1
c	print*,'four1: nn,n = ', nn,n
      do 11 i=1,n,2
        if(j.gt.i)then
          tempr=data(j)
          tempi=data(j+1)
          data(j)=data(i)
          data(j+1)=data(i+1)
          data(i)=tempr
          data(i+1)=tempi
        endif
        m=n/2
1       if ((m.ge.2).and.(j.gt.m)) then
          j=j-m
          m=m/2
        goto 1
        endif
        j=j+m
11    continue
c	print*,'four1: do11 done'
      mmax=2
2     if (n.gt.mmax) then
        istep=2*mmax
c	  if(istep > nn/2)print*,'istep,mmax=',istep,mmax
        theta=6.28318530717959d0/(isign*mmax)
        wpr=sin(0.5d0*theta) 
	  wpr= -2.d0*wpr*wpr
        wpi=sin(theta)
        wr=1.d0
        wi=0.d0
        do 13 m=1,mmax,2
	    !if(istep > nn/2)print*,'do13: m,n=',m,n
          do 12 i=m,n,istep
            j=i+mmax
            tempr=wr*data(j)-wi*data(j+1)
            tempi=wr*data(j+1)+wi*data(j)
            data(j)=data(i)-tempr
            data(j+1)=data(i+1)-tempi
            data(i)=data(i)+tempr
            data(i+1)=data(i+1)+tempi
12        continue
	!print*,'four1: do12 done'
          wtemp=wr
          wr=wr*wpr-wi*wpi+wr
          wi=wi*wpr+wtemp*wpi+wi
13      continue
c	  print*,'four1: do13 done. mmax=',mmax
        mmax=istep
        goto 2
      endif
c	print*,'four1: returning'
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software *j.
