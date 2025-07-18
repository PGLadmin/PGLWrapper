!*********************************************************
!
! SUBROUTINE GAMMA_PA
!
! Implement calls to calculate Wertheim, Physical, Combinatorial
! calculations and combine the results for Gammas and temperature
! derivative of gamma
!
!**********************************************************
SUBROUTINE GAMMA_PA(kop, Kcalc, T, P, gamma, dgamma)

use CONSTANTS
use sitenspecies
use VOLUMES
use PHYS_PARMS

implicit none

integer, INTENT(IN) :: Kcalc
real*8, INTENT(IN) :: T, P
real*8, dimension(nc), INTENT(OUT) :: gamma, dgamma

INTEGER i, j, k, l, kop(10)
! The intent of KOP is input, but the user for certain models inconsistent values
! can be specified, so the values can be overwritten below.
REAL*8 delT, rhomix, rhomixu, rhomixd
real*8, dimension(nc) :: gammares, dgammares, dgammaresd, gammacomb, &
 dgammacomb, dgammacombd, gammacombcorr, dgammacombcorr, dgammacombcorrd
real*8, dimension(nc) :: gammaw, dgammaw, dgammawd
real*8, dimension(nc) :: bvdw, VU, VD, V, rho, rhou, rhod
real*8, dimension(nc, nc) :: tau, Lambda, Vratio, VratioU, VratioD
! VrationU = volume ratio at T+delT/2
! VrationD = volume ratio at T-deltT/2

delT = 0.1D0

! Set default component molar volumes and size-related covolumes

do i = 1, nc
    V(i) = vstd(i) ! standard volume V in cm3/mol, overwrite later for T-dependent
	VU(i) = V(i) ! upper value used for T-derivative, same as V for standard volume
	VD(i) = V(i) ! down value used for T-derivative, same as V for standard volume
    bvdw(i) = r_uniquac(i)*15.17D0 !vdW covolume parameter
enddo

! these are the ln(gamma) and d(ln(gamma)/dT)
gamma = (/ (0D0, k=1,nc) /)
dgamma = (/ (0.D0, k=1,nc) /)  ! this is d(ln gamma)/dT
gammares = (/ (0D0, k=1,nc) /)
dgammares = (/ (0D0, k=1,nc) /)
gammacomb = (/ (0D0, k=1,nc) /)
dgammacomb = (/ (0D0, k=1,nc) /)
gammacombcorr = (/ (0D0, k=1,nc) /)
dgammacombcorr = (/ (0D0, k=1,nc) /)
gammaw = (/ (0D0, k=1,nc) /)
dgammaw = (/ (0D0, k=1,nc) /)
comp%x = (/ (x(k), k=1,nc) /)

! set site mole fractions
DO k = 1, nsites
    site(k)%xhost = comp(site(k)%host).x
END DO

if((kop(1).gt.1)) THEN
  write(debugfile,*) 'Sites Present: site present in call, site in file, id, name, host, host x'
  do i = 1,nsites
   write(debugfile,'(I3,2X,I4,2X,A20,2X,I3,2X,F10.5)') i, site(i)%id, site(i)%name, site(i)%host, site(i)%xhost
  enddo !i
 write(debugfile,*)'Site1, Site2, index1, index2, Keps, eps, kad'
 do l = 1, nsites
    do k = 1, nsites
        write(debugfile,'(2A20, 2I4, 3G12.5)') site(l)%name,site(k)%name,l,k,kad(l,k)*(dexp(eps(l,k)/T)-1D0), eps(l,k), kad(l,k)
    enddo !k
 enddo !l
endif


! ************** calculate the physical contribution *********************************************

if ((kop(1).gt.1).and.(mod(kcalc,2).gt.0)) write(debugfile,*) 'Component id, volume and covolume'

! KCALC=1, Calculate property only
! KCALC=2, Calculate temperature derivative of property only
! KCALC=3, Calculate property and temperature derivative.

! overwrite molar volume if kop(3) .gt. 0
if (kop(3) .gt. 0) then
    if (KCALC .gt. 1) then ! calculate values needed for T derivative
        do i = 1,nc
            CALL VLU(VU(i),i,VEQ(i),vparms(:,i),T + delT/2D0)
            CALL VLU(VD(i),i,VEQ(i),vparms(:,i),T - delT/2D0)
        enddo
    endif !kcalc > 1
	if (MOD(KCALC,2).gt.0) then	! calculate molar volume
        do i = 1,nc
            CALL VLU(V(i),i,VEQ(i),vparms(:,i),T)
        enddo
    endif !kcalc is odd
endif

if (kop(1).gt.1) then
! write components present in simulation, molar volume to be used, and user covolume.
    if (mod(kcalc,2).gt.0) then
        write(debugfile,*) ' volume and covolume'
        do i = 1, nc
            write(debugfile,'(I4, 2G12.5)') i, V(i), bvol(i)
        enddo
    endif
    if (kcalc.gt.1) then
        write(debugfile,*) ' vols at T+-deltT/2 for finite diff'
        do i = 1, nc
            write(debugfile,'(I4, 2G12.5)') i, VU(i), VD(i)
        enddo
    endif
endif

select case (kop(4))
    case (0)
		if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            call nrtl(dgammares, nc, x, aparam+bparam/(T + delT/2D0), alpha)
            call nrtl(dgammaresd, nc, x, aparam+bparam/(T - delT/2D0), alpha)
            dgammares = (dgammares-dgammaresd)/delT
		endif ! kcalc > 1
		if(MOD(kcalc,2) .gt. 0) then ! kcalc is odd
            tau = aparam+bparam/T
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMNRTL Physical tau'
                write(debugfile,'(6G12.5)') tau
            endif
            call nrtl(gammares, nc, x, tau, alpha)
        endif ! kcalc is odd

    case (1)
        do i=1,nc
            do j = 1,nc
                if (i.EQ.j) then
                    Vratio(i,j)= 1D0
					VratioU(i,j)=1D0
					VratioD(i,j)=1D0
                else
                    Vratio(i,j)= V(j)/V(i)
					! calculate upper and lower values for T-derivative
                    VratioU(i,j) = VU(j)/VU(i)
                    VratioD(i,j) = VD(j)/VD(i)
                endif
            enddo ! j
        enddo ! i
		if(kcalc .gt. 1) then ! todo implement wilson
		    ! calculate upper and lower values for T-derivative
            ! call wilson(dgammares, nc, x, VratioU*dexp(aparam+bparam/(T+delT/2D0)))
            ! call wilson(dgammaresd, nc, x, VratioD*dexp(aparam+bparam/(T-delT/2D0)))
            dgammares=(dgammares-dgammaresd)/delT
        endif ! kcalc > 1
        if( MOD(kcalc,2) .gt. 0) then ! kcalc is odd
            Lambda = Vratio*dexp(aparam+bparam/T)
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMUwilson Physical Lambda'
                write(debugfile,'(6G12.5)') Lambda
            endif
            ! call wilson(gammares, n, x, Lambda)
        endif ! kcalc is odd

    case (2) !scatchard hildebrand (eqn Elliott/Lira 12.24,12.25) todo implement SH
        if(kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            ! call scathild(dgammares, x, VU, nc, aparam*(T + delT/2D0)+bparam,R, T + delT/2D0)
            ! call scathild(dgammaresd, x, VD, nc, aparam*(T - delT/2D0)+bparam,R, T - delT/2D0)
            dgammares = (dgammares - dgammaresd)/delT
        endif ! kcalc < 1
        if (MOD(kcalc,2) .gt. 0) then !kcalc is odd
            Aij = aparam*T+bparam
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMUscathild Physical AIJ'
                write(debugfile,'(6G12.5)') AIJ
            endif
            ! call scathild(gammares, x, V, nc, Aij, R, T)
        endif ! kcalc is odd

	case (3) ! Nagata1 todo implement Nagata
        if (kcalc .gt. 1) then
		    ! calculate upper and lower values for T-derivative
            ! call nagata1(dgammares, nc, x, VU, aparam + bparam/(T + delT/2D0))
            ! call nagata1(dgammaresd, nc, x, VD, aparam + bparam/(T - delT/2D0))
            dgammares = (dgammares - dgammaresd)/delT
			WRITE(debugfile,*) 'Calculating temperature derivative of residual'
        endif
        if (MOD(kcalc,2).gt.0) then !kcalc is odd
            if(kop(1).gt.1) then
                write(debugfile,*) ' *** GMUnagata Physical aij + bij/T'
                write(debugfile,'(6G12.5)') aparam + bparam/T
            endif
            ! call nagata1(gammares, n, x, r_uniquac, aparam)
            ! call nagata1(gammares, n, x, V, aparam + bparam/T)
        endif

end select
! ********************************** end of gamma phyiscal ********************

! ********************************** Wertheim contribution ***********************

if (aspmx.eq.1) then ! only execute if association exists

	if (mod(kcalc,2) .gt. 0) then ! odd kcalc, calculate properties
		rho = 1D0/V ! molar density mol/cc
		rhomix = 1D0/(dot_product(x,V)) ! mol/cc

		CALL calc_gammaw(gammaw, kcalc, kop(1), kop(2), nc, nsites, comp, site, &
							   T, rhomix, rho, KAD, eps, bvol, &
							   PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
	endif

	if (kcalc .gt. 1) then ! need temperature derivative
		rhoU = 1D0/VU ! molar density mol/cc
		rhomixU = 1D0/(dot_product(x,VU)) ! mol/cc
		rhoD = 1D0/VD ! molar density mol/cc
		rhomixD = 1D0/(dot_product(x,VD)) ! mol/cc

		CALL calc_gammaw(dgammaw,  kcalc, kop(1),kop(2), nc, nsites, comp, site, &
							   (T + delT / 2D0), rhomixU, rhoU, KAD, eps, &
							   bvol, PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
		CALL calc_gammaw(dgammawd, kcalc, kop(1), kop(2), nc, nsites, comp, site, &
							   (T - delT / 2D0), rhomixD, rhoD, KAD, eps, &
							   bvol, PCSAFT_sigma, PCSAFT_m, PCSAFT_epsok)
		dgammaw = (dgammaw - dgammawd) / delT
	endif
endif !aspmx eq 1

! ------------ end of section if association occurs in mixture-----------

! ************************ end of Werthiem contribution ********************************

! combinatorial term always uses volume, not r
if (kop(4) .ne. 1) then  ! if not Wilson
    if(kcalc.gt.1) then
	if (kop(6) .eq. 0) then
		call flory(dgammacomb, x, VU, nc)
		call flory(dgammacombd, x, VD, nc)
		dgammacomb = (dgammacomb - dgammacombd)/delT
	elseif (kop(6) .eq. 1) then
		call flory(dgammacomb, x, VU**(2D0/3D0), nc)
		call flory(dgammacombd, x, VD**(2D0/3D0), nc)
		dgammacomb = (dgammacomb - dgammacombd)/delT
	endif
    endif !kcalc > 1
    if(mod(kcalc,2).gt.0) then
	if (kop(6) .eq. 0) then
		call flory(gammacomb, x, V, nc)
        if (kop(1).gt. 1) WRITE(debugfile,*) '*** Combinatorial Method: Flory using molar volume'
	elseif (kop(6) .eq. 1) then
		call flory(gammacomb, x, V**(2D0/3D0), nc)
        if (kop(1).gt. 1) WRITE(debugfile,*) '*** Combinatorial Method: Flory using molar volume^(2/3)'
	endif
	endif
endif

! **********Combinatorial Correction
if ((kop(4) .eq. 1) .or. (kop(6) .eq. 2)) then
    kop(5) = 1; ! force combinatorial correction off if using Wilson, or if combinatorial term is off
    ! note the gammacombcorr is set to zero at the beginning of the routine
    if(kop(1).gt.1) WRITE(debugfile,*) '*** Your selection of Comb Correction May Be Inconsistent with Other Settings'
    if(kop(1).gt.1) WRITE(debugfile,*) '*** Combinatorial Correction Method: None'
endif
if (kop(5) .ne. 1) then
    if(kcalc.gt.1) then
        select case(kop(5))
        case (0)
            call correqn1262(dgammacombcorr, x, VU, nc, bvdw)
            call correqn1262(dgammacombcorrd, x, VD, nc, bvdw)
            dgammacombcorr = (dgammacombcorr - dgammacombcorrd)/delT
        case (2)
            dgammacombcorr = 0
        case (3)
	        call correqn1262(dgammacombcorr, x, VU, nc, bvol)
            call correqn1262(dgammacombcorrd, x, VD, nc, bvol)
            dgammacombcorr = (dgammacombcorr - dgammacombcorrd)/delT
        end select
    endif
    if(mod(kcalc,2).gt.0) then
        select case(kop(5))
        case (0)
            call correqn1262(gammacombcorr, x, V, nc, bvdw)
            if (kop(1).gt. 1) WRITE(debugfile,*) '*** Combinatorial Correction Method: vdW using bvdw'
        case (2)
            call stavgugg(gammacombcorr, x, r_uniquac,q_uniquac, nc)
            if (kop(1) .gt. 1) WRITE(debugfile,*) '*** Combinatorial Correction Method: Staverman-Guggenheim usig r and q'
        case (3)
	        call correqn1262(gammacombcorr, x, V, nc, bvol)
            if (kop(1) .gt. 1) WRITE(debugfile,*) '*** Combinatorial Correction Method: vdW using COVOL'
        end select
    endif
end if

! add logs of physical and wertheim contributions
gamma = gammares + gammaw + gammacomb + gammacombcorr
dgamma = dgammares + dgammaw + dgammacomb + dgammacombcorr

if(kop(1) .gt. 0) then ! write to history and/or .csv files
	  WRITE(debugfile,*) '*** GAMMA RESULTS'
	  WRITE(debugfile,900) KCALC, nc, T, P
 900  FORMAT(' KCALC, N, T, P =', I3, 1X, I3, F10.3, 2X, E13.6)
	  if(mod(kcalc,2).gt.0) then !kcalc is odd, write gammas to history
		WRITE(debugfile,1010)
		DO 400 I = 1, nc! write gamma contributions to history
		 WRITE(debugfile,1000) i, comp(i)%name, x(i), gammares(i),gammacomb(i),gammacombcorr(i), gammaw(i), GAMMA(I), dgamma(i)
 400   CONTINUE
	! frequent opening/writing/closing file actions --> potential optimization: save in mem buffer, write at once?
	   if(kop(1) .gt. 2) then
		! write to vol.csv
		open(unit=12, file="..\..\Output\vol.csv", status='unknown', action='write', position='append')
		WRITE(12,1002) T, ',', 1D0/rhomix, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", V(i), i=1,nc)
		close(12)
		! write to gammas.csv
		open(unit=10, file="..\..\Output\gammas.csv", status='unknown', action='write', position='append')
		WRITE(10,1001) T, ',', P, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", gammares(i),",",gammacomb(i),",",gammacombcorr(i), ",",gammaw(i), ",",GAMMA(I), ",",dgamma(i), i=1,nc)
		close(10)
	   endif !kop(1).gt.2
	  endif !mod(kcalc,2)
	  if(kcalc.gt.1) then !write derivative
		WRITE(debugfile,1011)
		DO 410 I = 1, nc
		 WRITE(debugfile,1000) i, comp(i)%name, X(I), dgammares(i),dgammacomb(i),dgammacombcorr(i), dgammaw(i), GAMMA(I), dgamma(i)
 410    CONTINUE
	! frequent opening/writing/closing file actions --> potential optimization: save in mem buffer, write at once?
		if(kop(1).gt.2) then
		! write to HXS.csv
		open(unit=11, file="..\..\Output\HXS.csv", status='unknown', action='write', position='append')
		WRITE(11,1001) T, ",", P, (', ', i, ',  ', comp(i)%name, ",", X(I), ",", -dgammares(i)*R*(T**2), ",", -dgammacomb(i)*R*(T**2), ",", -dgammacombcorr(i)*R*(T**2), ",", -dgammaw(i)*R*(T**2), ",", -dgamma(i)*R*(T**2), ',', dgamma(i), i=1,nc)
		close(11)
		endif !if kcalc.gt.1
		endif !kop(1).gt.2

		WRITE(debugfile,*) 'End of function call'
	  WRITE(debugfile,"(120('*'))")
endif !kop(1).gt.0

!
!
!
!     Format statements
!
 1000   FORMAT (I3,1X, A8, 2X, 7F15.6)
 1001	FORMAT (F15.6,A1,F15.6,100(A2, I3, A3, A20, 7(A1,F15.6))) ! format for gamma, hex csv files, up to 100 components
 1002	FORMAT (F15.6,A1,F15.6,100(A2, I3, A3, A20, 2(A1,F15.6))) ! format for vol csv files, up to 100 components
 1010   FORMAT (' I     NAME          X         lnGAMMAres    lnGAMMAcomb  lnGAMMAcombcorr       lnGAMMAw       lnGAMMA      dlnGAMMA')
 1011   FORMAT (' I     NAME          X         dlnGAMMAres   dlnGAMMAcomb dlnGAMMAcombcorr     dlnGAMMAw       lnGAMMA      dlnGAMMA')
!
return
end subroutine GAMMA_PA
