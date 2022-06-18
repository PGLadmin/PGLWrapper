c
c----------------------------------------------------------
c     thermodynamics block
c
      subroutine thermo(icon,ntyp,mtyp,t00,p00,zz,zn,fg,ft,fp,fx,aux)
      implicit double precision (a-h,o-z)
      include '..\source\comflash.cmn'      
      dimension zn(*) , fg(*) , ft(*) , fp(*) , fx(pmax,*) , aux(*)
      ntemp = 0
      nder = 1
      if ( icon.gt.2 ) nder = 2
      if ( icon.eq.2 .or. icon.gt.3 .or. icon.eq.0 ) ntemp = 2
      if ( icon.gt.4 ) ntemp = 2
      call cubgen(ntyp,mtyp,t00,p00,zz,zn,fg,ft,fp,fx,aux)
      end
c
c----------------------------------------------------------
c
      subroutine vthermo(icon,t00,v00,p00,dp_dt,dp_dv,x,fg,ft,fp,fx,ierr)
c
c-------parameters of vtermo (crit. point calc)
c       icon  level of derivatives                     (input}
c       1     p and fg only
c       2     also t- and v-derivatives
c       3     also n-derivativs
c       4     all derivatives
c
c       t     temperature (k)                          (input)
c       v     total volume/r    (k*mol/mpa)            (input)
c       p     pressure    (mpa)                        (output)
c       dpdt  derivative of p/t wrt. t                 (output)
c       dpdv  derivative of p/t wrt. v                 (output)
c       x---  mixture mole numbers                     (input)
c       fg    logarithm of fugacity/mole number ratio  (output)
c       ft    t-derivative of fg (const. vol)          (output)
c       fv    vol-derivative of fg (const temp)        (output)
c       fx    comp-derivative of fg (const t. and v)   (output)
c       ierr  error indicator     (nonzero for error)  (output)
c---------------------------------------------------
      implicit double precision (a-h,o-z)
      include '..\source\comflash.cmn'
      dimension x(*) , fg(*) , ft(*) , fp(*) , fx(pmax,*)
      ntemp = 0
      nder = 1
      if ( icon.gt.2 ) nder = 2
      if ( icon.eq.2 .or. icon.gt.3 .or. icon.eq.0 ) ntemp = 2
      if ( icon.gt.4 ) ntemp = 2

      call volgen(t00,v00,p00,dp_dt,dp_dv,x,fg,ft,fp,fx,ierr)
      end
c
c
c----------------------------------------------------------
c
      subroutine cubgen(mt,ic,t,p,zz,xx,fg,ft,fp,fx,aux)
c
c     parameters:
c
c     mt:     (i):      phase type desired
c                       1 = liq, -1 = vap, 0 = min g phase
c     ic:     (o):      indicator for phase type; ic returns
c                       1:   if called with mt = 1/-1 and phase is found
c                       -1:  if called with mt = 1/-1 and phase is not found
c                       2:   if called with mt = 0 and min g phase is liquid
c                       -2:  if called with mt = 0 and min g phase is vapour
c     t:      (i):      temperature (k)
c     p:      (i):      pressure (MPa)
c     z:      (o):      calculated compressibility factor
c     xx:     (i):      mixture total composition
c     fg:     (o):      log fugacity coefficient vector
c     ft:     (o):      t-derivative of fg
c     fp:     (o):      p-derivative of fg
c     fx:     (o):      scaled composition derivative of fg
c     aux:    (o):      various residual properties
c
      implicit double precision (a-h,o-z)
c-----dummy variables
       include '..\source\comflash.cmn'
      dimension xx(*) , fg(*) , ft(*) , fp(*) , fx(pmax,*) , aux(*)
c-----local variables
      dimension x(pmax),pd(pmax),ad1(pmax),adt(pmax)
      logical pspec

      if(mt.eq.1) j = 1
      if(mt.eq.-1) j = 2
      if(mt.eq.0) then
          istab = 1
          j = 3
      end if
      n = pmax
      sumx = 0.d0
      do i = 1 , n
         sumx = sumx + xx(i)
      end do
      sumr = 1.d0/sumx
      do i = 1 , n
         x(i) = xx(i)*sumr
      end do             
      x(1:n) = xx(1:n)     
      
      t0 = t
      p0 = p*10.0d0 !Mpa ----> bar
      ider = 1
      if(isGPED) pc = pc * 10.0d0
      call eqet0
      call tempe(0)
      call regles_melange(ider,2,x)

      call phas1(j,ider,istab,p0,x,eta0,v0,g0)
      zz = z

      g0r = 0.0d0
      do i=1,nco
         g0r = g0r + x(i)*lphi(i)
      end do      

       aux(10) = dpdv / t * 1.0d-01 * rgas
       aux(8) = dpdv * 1.0d-01 * rgas   
c
c     are temperature derivatives required
c
      if ( ntemp.ne.0 ) then
         call hec_calc(x,t,v0,hec)
c------excess enthalpy/r stored in aux(3)
         aux(3) = hec / rgas
         sec = (hec - g0r*rgas*t) / t / rgas
c------excess entropy/r stored in aux(4)
         aux(4) = sec
c
c     residual cp required
c
         if ( ntemp.ge.2 ) then
            call Cpec_calc(x,t,v0,Cpec)         
c------excess heat capacity/r in aux(5)
            aux(5) = Cpec / rgas
            dvdt = -dpdt / dpdv / rgas
c---- aux(6) is pressure derivative of resid. entropy
            aux(6) = -dvdt + 1.0d0 / p
c-----aux(7) is pressure derivative of enthalpy
            aux(7) = -t*dvdt + v0/rgas
            aux(9) = dvdt
         endif
      endif

      if ( nder.ne.0 ) then
c
c     log fugacity coefficients
c
         fg(1:n) = lphi(1:n) - dlog(peq)
         if ( nder.ge.2 .or. ntemp.ne.0 ) then
c
c     composition derivatives required
c  
            if ( nder.ge.2 ) then
               do i = 1,n
                 do k = 1,n
                    fx(i,k) = dmudn(i,k) + 1.0d0 
                 enddo
                 fx(i,i) = fx(i,i) - 1.0d0/x(i) 
               enddo
            endif
c
c     temperature and pressure derivatives
c
            if ( ntemp.ne.0 ) then
               fp(1:n) = dlphidp(1:n)*10.0d0 !bar-1 -----> MPa-1
               ft(1:n) = dlphidt(1:n)
            endif
         endif
      endif

      if(isGPED) pc = pc / 10.0d0        
      if(istab.ne.2) return  
      
      end             
c
c     auxillary routines
c

      subroutine volgen(t,vv,p,dpdtv,dpdvv,x,fg,ft,fv,fx,ierr)
      implicit double precision (a-h,o-z)
      include '..\source\comflash.cmn'
      dimension x(*),fg(*),ft(*),fv(*),fx(pmax,*),x2(pmax)
      dimension dfgdv(pmax),dfgdn(pmax,pmax)
c
c-------parameters of volgen (crit. point calc)
c
c       t     temperature (k)                          (input)
c       v     total volume/r    (k*mol/mpa)            (input)
c       p     pressure    (mpa)                        (output)
c       dpdt  derivative of p/t wrt. t                 (output)
c       dpdv  derivative of p/t wrt. v                 (output)
c       x---  mixture mole numbers                     (input)
c       fg    logarithm of fugacity/mole number ratio  (output)
c       ft    t-derivative of fg (const. vol)          (output)
c       fv    vol-derivative of fg (const temp)        (output)
c       fx    comp-derivative of fg (const t. and v)   (output)
c---------------------------------------------------
c---  modified and corrected april 2021
c---
c---------------------------------------------------
c
c     sumn is total moles; volgen calculates in terms of overall moles
c     and overall volume
c
      if(isGPED) pc = pc * 10.0d0

      n = pmax
      sumn = 0.d0
      do i = 1,n
         sumn = sumn + x(i)
      enddo
      
      x2(1:n) = x(1:n) / sumn
      
      t0 = t
      call eqet0
      call tempe(0)
      ider = 1
      call regles_melange(ider,2,x2)

      ierr =0
      pn = rgas/(vv*rgas-bm)

      if (pn.lt. 0.d0) then
         ierr = 1
         if(isGPED) pc = pc / 10.0d0
         return
      endif      
c     Eta calculation       
      vv2 = vv/sumn  
      eta = bm/(vv2*rgas)

      call eqet !Calculation of z and P
      p = peq * 1.0d-01 ! bar ---> MPa 
      call fug !Calculation of ln phi

      if ( nder.ne.0 ) then
c       fg    logarithm of fugacity/mole number ratio  (output)      
         do i = 1,n
           fg(i) = lphi(i) + dlog(1.0d-01) - dlog(sumn)
         end do
         if ( nder.ge.2 .or. ntemp.ne.0 ) then
            ider = 2
            call dertp(ider,x)
c
c  volume derivative of fugacity from composition deriv. of pressure
c
            fv(1:n) = dmudv(1:n) + dpdv/peq
            dfgdv(1:n) = - fv(1:n) * bm/eta / sumn
            fv(1:n) = fv(1:n) * rgas / sumn
            dpdvv = dpdv * rgas / t * 1.0d-01 / sumn            
            if ( nder.ge.2 ) then
                do i = 1,n
                  do k = 1,n
                    som = 0.0d0
                    do j = 1,n
                       if(j == k) som = som + dmudn(i,j) * (1.0d0-x2(j))
                       if(j /= k) som = som - dmudn(i,j) * x2(j)
                    enddo      
                     dfgdn(i,k) = som / sumn
                     fx(i,k) = dfgdv(i) + dfgdn(i,k) - 1.0d0/sumn
                  enddo        
                enddo
                fx(2,1) = fx(1,2)
            endif
            if ( ntemp.gt.0 ) then
                dpdtv = dpdt / t * 1.0d-01  - p / t**2
                ft(1:n) = dmudt(1:n) + dpdt/peq
            endif
         endif
      endif

      if(isGPED) pc = pc / 10.0d0
      end


c----------------------------------------------------------
