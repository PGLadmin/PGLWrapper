!==============================================================================
!MODULE SECTION
!==============================================================================
module comflash
        implicit none
!c VERSION NOVEMBER 2020
        integer,parameter :: lun1 = 6
        integer,parameter :: lun2 = 55
        integer,parameter :: ncomax = 3
        integer,parameter :: phmax = 2
        integer,parameter :: nbgr_tot_max=40
        integer,parameter :: tempcmax = 500
        integer,parameter :: tempdmax = 10
        integer,parameter :: ntempmax = 9000
        integer,parameter :: parmax = 3
        integer,parameter :: nvarmax = 25

!C      List of index for the property array
        integer,parameter :: h_id = 1   !Enthalpy [J/mol]
        integer,parameter :: s_id = 2   !Entropy [J/mol/K]
        integer,parameter :: u_id = 3   !Internal energy [J/mol]
        integer,parameter :: cp_id = 4  !Isobaric heat capacity [J/mol/K]
        integer,parameter :: cv_id = 5  !Isochoric heat capacity [J/mol/K]
        integer,parameter :: qx_id = 6  !Quality (molar basis)
        integer,parameter :: dpt_id = 7 !(dP/dT)@V
        integer,parameter :: dvp_id = 8 !(dV/dP)@T
        integer,parameter :: dvt_id = 9 !(dV/dT)@P
        integer,parameter :: p_id = 10  !Pressure [bar/ Pa]
        integer,parameter :: t_id = 11  !Temperature [K]
        integer,parameter :: v_id = 12  !Bulk molar volume [m3/mol]
        integer,parameter :: vl_id = 13 !Liquid-phase molar volume [m3/mol]
        integer,parameter :: vv_id = 14 !Vapor-phase molar volume [m3/mol]
        integer,parameter :: lh_id = 17 !Heat of vaporization [J/mol]
        integer,parameter :: jt_id = 18 !Joule-Thompson coeff [K/Pa]
        integer,parameter :: b_id = 19 !Volumetric expansivity [1/K]
        integer,parameter :: ph_id = 20 !Phase region index [-]
        integer,parameter :: thc_id = 21 !Thermal conductivity [W/m/K]
        integer,parameter :: vs_id = 22 !Dynamic Viscosity [Pa.s]
        integer,parameter :: thd_id = 23 !Thermal diffusivity [cm2/s]
        integer,parameter :: kvs_id = 24 !Kinematic Viscosity [cm2/s]
        integer,parameter :: dpts_id = 25 !dP/dT (along the saturation line)

        integer, parameter :: nb_ids = 5
        integer, parameter :: nb_cstvalues = 9
        integer, parameter :: nb_tpty = 15
        integer, parameter :: nb_coeff = 7

        integer, parameter :: mw_idx = 1
        integer, parameter :: tc_idx = 2
        integer, parameter :: pc_idx = 3
        integer, parameter :: acen_idx = 4
        integer, parameter :: vc_idx = 5
        integer, parameter :: zc_idx = 6
        integer, parameter :: tpt_idx = 7
        integer, parameter :: tpp_idx = 8
        integer, parameter :: nbp_idx = 9

        ! regarder notice dippr pour plus info sur les noms
        integer, parameter :: ldn_idx = 1
        integer, parameter :: vp_idx = 2
        integer, parameter :: hvp_idx = 3
        integer, parameter :: lcp_idx = 4
        integer, parameter :: icp_idx = 5
        integer, parameter :: svr_idx = 6
        integer, parameter :: sdn_idx = 7
        integer, parameter :: svp_idx = 8
        integer, parameter :: scp_idx = 9
        integer, parameter :: lvs_idx = 10
        integer, parameter :: vvs_idx = 11
        integer, parameter :: stc_idx = 12
        integer, parameter :: ltc_idx = 13
        integer, parameter :: vtc_idx = 14
        integer, parameter :: st_idx = 15


        real(8),parameter :: r = 83.14462618d0		!83.14411d0 ! bar.cm3.mol^(-1).K^(-1) modified to PGL6ed values
        real(8),parameter :: rgas = 8.314462618d0	!8.314411d0 ! J/mol/K
        real(8) :: t_ref,v_ref

        real(8) :: t0,p0,ntt,vt,mmolt,zt,denst,dropout,zp(phmax)
        real(8) :: nt(ncomax),npt(phmax),vp(phmax),vph(phmax)
        real(8) :: densp(phmax),np(ncomax,phmax),xp(ncomax,phmax)
        real(8) :: mmolp(phmax),volat(ncomax)
        integer nco
        character*15 nom(ncomax)


        real(8) :: epseq,eps0,eps1,l0,ld(10)
        integer :: nimp,ph3,iac,isubs,cam,oa,nb,nlp,nls

        integer,parameter :: ncompmaxcp = 300
        integer,parameter :: ncompmaxcp2 = 2000

        integer :: choix_eos,choix_b,choix_ge,choix_h,choix_cP
        integer :: indice_cp(ncompmaxcp2),chemlist(ncomax)
        real(8) :: tc(ncomax),pc(ncomax),fac(ncomax)
        real(8) :: zra(ncomax),wmas(ncomax)
        real(8) :: ac(ncomax),msoave(ncomax),pcmax,parC(ncomax)
        real(8) :: parL(ncomax),parM(ncomax),parN(ncomax)
        real(8) :: qi(ncomax)
        real(8) :: tpt(ncomax),tpp(ncomax),tnbp(ncomax)
        real(8) :: href(ncomax),sref(ncomax),uref(ncomax)
        real(8) :: tabcpgp(ncompmaxcp,6),vc(ncomax),omega(ncomax)

        real(8) :: eta,peq,dpdeta,gamma,rt,v,dzdeta,qpseta,rtsb,y,z,q

        integer :: choix_flash,nombre_flash
        real(8) :: aa,bm,ee,am,amix,bmix,cm,cmix
        real(8) :: b(ncomax),asb(ncomax),bisb(ncomax),c(ncomax)
        real(8) :: aai(ncomax),bx(ncomax),cx(ncomax)
        real(8) :: e(ncomax),bc(ncomax)
        real(8) :: dEijsRT(ncomax,ncomax),d2EijsRT(ncomax,ncomax)
        real(8) :: daidx(ncomax,ncomax),en(ncomax,ncomax),k12
        real(8) :: Lij(ncomax,ncomax),kij(ncomax,ncomax),Lambda
        real(8) :: k12_const,puisb,att,at,daai_dT(ncomax)
        real(8) :: agE0(ncomax,ncomax),bgE0(ncomax,ncomax)
        real(8) :: cgE0(ncomax,ncomax)
        real(8) :: pargE(ncomax,ncomax),dpargE(ncomax,ncomax)
        real(8) :: d2pargE(ncomax,ncomax)
        real(8) :: expA(ncomax,ncomax),sigma(ncomax),mu(ncomax)
        real(8) :: alphgE(ncomax,ncomax)
        real(8) :: delta_sur_RacRT(ncomax)
        character(30) :: modele_gE
        character(50) :: tabingr_file,gE_file
        logical :: quadra_b

        real(8) :: gt,gt0
        real(8) :: etap(phmax),gp(phmax)

        real(8) :: theta1,theta2,omegab,omegaa
        real(8) :: mmol(ncomax),ss(ncomax)
        real(8) :: dpdv,dpdt
        real(8) :: lphi(ncomax),dlphideta(ncomax)
        real(8) :: dpdn(ncomax),dmudv(ncomax)
        real(8) :: dmudt(ncomax),dpdpar(parmax)
        real(8) :: dlphidp(ncomax),detadt
        real(8) :: dlphidn(ncomax,ncomax)
        real(8) :: dmudn(ncomax,ncomax),dmudpar(ncomax,parmax)
        real(8) :: dlphidt(ncomax),dlphidv(ncomax)

        real(8) :: fmaxi,l12
        integer :: iph, nitflash, m0,n0

        real(8) :: agr(nbgr_tot_max,nbgr_tot_max)
        real(8) :: bgr(nbgr_tot_max,nbgr_tot_max)
        real(8) :: srg(ncomax,nbgr_tot_max)
        real(8) :: asb0(ncomax)
        integer :: nbgr_tot,ngroupe
        integer :: knum(10),imin(10),imax(10)

        real(8) :: tmemc(tempcmax),tmemd(tempdmax)
        real(8) :: kijt(tempcmax,ncomax,ncomax)
        real(8) :: ent(tempcmax,ncomax,ncomax)
        real(8) :: dendtt(tempdmax,ncomax,ncomax)
        real(8) :: pargEt(tempcmax,ncomax,ncomax)
        real(8) :: expAt(tempcmax,ncomax,ncomax)
        integer :: itempc

        integer :: ngr,nbgr_used,liste_gr(nbgr_tot_max)
        integer :: ngr1,nomgr1(nbgr_tot_max)
        integer :: ngr2,nomgr2(nbgr_tot_max)

        real(8) :: occurgr1(nbgr_tot_max),occurgr2(nbgr_tot_max)
        real(8) :: pres_max

        integer :: ncomp1,ncomp2,codcomp1,codcomp2
        real(8) :: agr0(nbgr_tot_max,nbgr_tot_max)
        real(8) :: bgr0(nbgr_tot_max,nbgr_tot_max)
        character(15) :: comp1,comp2

        integer :: nder,ntemp
        real(8) :: ac0(ncomax),ac1(ncomax),ac2(ncomax)
        real(8) :: dlna_T(ncomax),d2lna_T(ncomax)

        integer :: choix_r1r2,n_param
        real(8) :: r1,r2,derbm(ncomax),dercm(ncomax)
        real(8) :: dalphadx(ncomax),dbmdx(ncomax),dcmdx(ncomax)
        real(8) :: der2bm(ncomax,ncomax)
        real(8) :: bij(ncomax,ncomax),aij(ncomax,ncomax)
        character(5) :: FoncAlpha
        real(8) :: dE_dT, d2E_dT2, dalpha_dT, d2alpha_dT2


        logical :: ygtx,pc_up,fl_ok,corrFL
        real*8  :: pvap(ncomax),lnphi_L(ncomax),lnphi_V(ncomax)

        logical :: kij_const = .false.
        logical :: predictif = .false.
        logical :: fitting = .false.
        logical :: simul_PRWil = .false.
        logical :: fitting_PRWil = .false.

        real(8) :: z_gauss(10),w_gauss(10)

!C      general choices (EoS, alpha function)
        integer :: choice(5)  = 0
        integer :: nb_comp

        integer, parameter  :: n_prop = 15

        real(8) :: tdep_prop(ncomax,n_prop,7)
        real(8) :: tmin_prop(ncomax,n_prop)
        real(8) :: tmax_prop(ncomax,n_prop)
        integer :: eqIDComp(ncomax,n_prop)
        character(70) :: idComp(ncomax,5)

        character(250) :: errm
contains  ! JRE20210513: "contains" of flash() in module keeps it private relative to other codes like UaWrapper.
!c*****************************************************************************
!	Added for PGL6ed
!c*****************************************************************************
!c Notes for general fugacity calls:
!c       icon  level of derivatives
!c       1     p and fg only
!c       2     also t- and v-derivatives
!c       3     also n-derivativs
!c       4     all derivatives
!c     ntyp:     (i):      phase type desired
!c                       1 = liq, -1 = vap, 0 = min g phase
!c Notes for fugacity_TP calls:
!c     mtyp:     (o):      indicator for phase type; ic returns
!c                       1:   if called with mt = 1/-1 and phase is found
!c                       -1:  if called with mt = 1/-1 and phase is not found
!c                       2:   if called with mt = 0 and min g phase is liquid
!c                       -2:  if called with mt = 0 and min g phase is vapour
!c
!c---------------------------------------------------------------------------
!c
      subroutine flash(imp,t,p00,nd)
!c
      !use comflash
      implicit none

!c
      character*80 name
      integer imp,no
      integer i,init,j
      integer j0,ider,istb,it,ph,er1,er2,l1,l2,npoints
      real*8 t,p00,nd(ncomax),eps,tempo
      real*8 eta0,v0,g0,x1u,y1u,x1s,y1s,z1old,epsmin,phi1,phi2,res
      real*8 x0(ncomax),z1
      real*8 xpold(ncomax,phmax),res2,xmin,xxx,rti,xd(ncomax)
      real*8 xinf,xsup,stepx,presqun,rien,psat1,psat2
      real*8 densmin,densmax,ecart_densite,ps_max,ps_min
      real*8 min1,min2,vliq,vgaz

      logical stable,croissant,on_inverse,traitement_special
      logical bleme_sup
      common/is_stable/stable

      real*8 tpdminimum ! modif rp tempo
      common/mapoubelle/tpdminimum ! modif rp tempo

      corrFL = .FALSE.
      nombre_flash = nombre_flash + 1
      xinf = 0.0d0
      xsup = 1.0d0

      no = 0
      if(dabs(t-t0).gt.1.0D-13) then  !entre 2 flash successifs on n'a pas
          t0 = t                      !la meme temperature.
          no = 1                      !t0 est definie ici
      endif

      if(no.ne.0) then
          call tempe(0)
          if(dabs(k12).gt.50.d0) then
             iph = 0
             name = nom(1)
             call longueur(name,l1)
             name = nom(2)
             call longueur(name,l2)
             write(lun1,131) t0,nom(1)(1:l1),nom(2)(1:l2),k12
             write(lun2,131) t0,nom(1)(1:l1),nom(2)(1:l2),k12
 131         format(1x,'A T/K= ',F7.2,' kij(',A,'-',A,') =',g17.10, &
     &              ' on ne resout pas le flash')
             return
          endif
!c
!c         la valeur de no utilisee est celle des arguments si la temperature
!c         n'a pas varie.
!c
      endif
!c
      p0 = p00
      ntt = 0.0D+0
      do i=1,nco
         ntt = ntt + nd(i)
         nt(i) = nd(i)
      end do
      xd(1:ncomax) = nd(1:ncomax) / ntt
!c
!c      print*,'VERIF DER'
!c      call verif_der ! modif rp tempo
!c      print*,'bm = ',bm
!c      stop

      if(ph3.eq.0) then
         init = 0
         call flash2(init)
      else
!c         call stab3
      endif

!c     Si une fraction mol. est si proche de 0 que son complement
!c     a 1 vaut numeriquement 1.d0, on met le complement a 1.d0 - 1.d-16
!c     pour eviter BUG dans phas1
!c     Exemple : x1 = 1.d-17 => x2 = 1 - x1 vaut numeriquement 1.d0

      if(xp(1,1).eq.1.d0 .or. xp(2,1).eq.1.d0 .or. &
     &   xp(1,2).eq.1.d0 .or. xp(2,2).eq.1.d0) then
        if(xp(1,1).eq.1.d0) then
           xp(1,1) = 1.d0 - 1.d-16
           xp(2,1) = 1.d0 - xp(1,1)
        endif
        if(xp(2,1).eq.1.d0) then
           xp(2,1) = 1.d0 - 1.d-16
           xp(1,1) = 1.d0 - xp(2,1)
        endif
        if(xp(1,2).eq.1.d0) then
           xp(1,2) = 1.d0 - 1.d-16
           xp(2,2) = 1.d0 - xp(1,2)
        endif
        if(xp(2,2).eq.1.d0) then
           xp(2,2) = 1.d0 - 1.d-16
           xp(1,2) = 1.d0 - xp(2,2)
        endif
      endif


!c     Y a-t-il une fraction mol. > 1 et tres proche de 1
!c     Si oui, on la met a 1.d0 - 1.d-16
      if(xp(1,1).gt.1.d0 .or. xp(2,1).gt.1.d0 .or. &
     &   xp(1,2).gt.1.d0 .or. xp(2,2).gt.1.d0) then
        if(xp(1,1).gt.1.d0 .AND. dabs(xp(1,1)-1.d0).lt.1.d-15) then
           xp(1,1) = 1.d0 - 1.d-16
           xp(2,1) = 1.d0 - xp(1,1)
        endif
        if(xp(2,1).gt.1.d0 .AND. dabs(xp(2,1)-1.d0).lt.1.d-15) then
           xp(2,1) = 1.d0 - 1.d-16
           xp(1,1) = 1.d0 - xp(2,1)
        endif
        if(xp(1,2).gt.1.d0 .AND. dabs(xp(1,2)-1.d0).lt.1.d-15) then
           xp(1,2) = 1.d0 - 1.d-16
           xp(2,2) = 1.d0 - xp(1,2)
        endif
        if(xp(2,2).gt.1.d0 .AND. dabs(xp(2,2)-1.d0).lt.1.d-15) then
           xp(2,2) = 1.d0 - 1.d-16
           xp(1,2) = 1.d0 - xp(2,2)
        endif
      endif

      if(choix_flash.eq.1) goto 1477

!c
!c ===========================================================
!c Analyse des solutions du flash de Peneloux/Gramajo a partir
!c    du critere du plan tangent (voir les notes de Romain)
!c   Parfois le flash est re-resolu selon methode R.P.
!c ===========================================================
!c
!c     Verification operee seulement si les solutions du flash
!c     ne sont pas trop proches des axes x1 = 0 et x1 = 1
      if(iph.eq.2) then
         if(xp(1,1).lt.0.002d0.and.xp(1,2).lt.0.002d0) goto 1477
         if(xp(1,1).gt.0.998d0.and.xp(1,2).gt.0.998d0) goto 1477
      else
         if(xp(1,1).lt.0.002d0) goto 1477
         if(xp(1,1).gt.0.998d0) goto 1477
      endif

      stable = .true.
      epsmin = 1.d-13
      xpold = xp

!c     croissant se souvient de l'ordre de rangement des solutions
!c     croissant = .TRUE.  => xp(1,1) < xp(1,2)
!c     croissant = .FALSE. => xp(1,1) > xp(1,2)

      if(xp(1,1).lt.xp(1,2)) then
         croissant = .true.
      else
         croissant = .false.
      end if

!c     on lance le TPD sur la fraction la plus eloignee des fractions x=0
!c     et x=1
      min1 = min(xp(1,1),1.d0-xp(1,1))
      min2 = min(xp(1,2),1.d0-xp(1,2))
      if(min1.lt.min2) then
         z1 = xp(1,2)
      else
         z1 = xp(1,1)
      endif
      call tpd_rp(z1,t,p00)
!c      print*,'xp = ',xp(1,1:2) ! modif rp tempo
!c      print*,'tpdmin',tpdminimum,z1,stable ! modif rp tempo
!
!c     Si la solution du TPD est monophasique stable, on verifie l'egalite
!c     des potentiels chimiques entre xp(1,1) et xp(1,2)
      if(stable.and.iph.eq.2) then

!c        Si une des fractions molaires calculees par le Flash Peneloux
!c        vaut exactement 1.d-16, c'est que le flash Peneloux
!c        renvoyait une valeur numerique de la fraction molaire
!c        complementaire egale a 1.d0 et que cette
!c        valeur a ete modifiee pour eviter des bugs dans PHAS1.
!c        Dans ce cas, la verif de l'egalite des potentiels chimiques
!c        ne sera pas concluante donc on la contourne
         if(xp(1,1).eq.1.d-16 .or. xp(2,1).eq.1.d-16 .or. &
     &      xp(1,2).eq.1.d-16 .or. xp(2,2).eq.1.d-16) then
            goto 1599
         endif

         j0    = 1
         ider  = 0
         istb  = 1

         x0(1) = xp(1,1)
         x0(2) = 1.d0 - x0(1)
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
         phi1 = dexp(lphi(1))

         x0(1) = xp(1,2)
         x0(2) = 1.d0 - x0(1)
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
         phi2 = dexp(lphi(1))

         res = dabs(xp(1,1) * phi1 - xp(1,2) * phi2)

         if(res.gt.1.d-5) then
!c           Non egalite des potentiels chimiques
!c            print*,'Mauvais calcul Peneloux a T = ',t,res

!c           On trie les solutions dans l'ordre croissant
            tempo = max(xp(1,1),xp(1,2))
            xp(1,1) = min(xp(1,1),xp(1,2))
            xp(1,2) = tempo

            eps = 0.001d0

!c           On cherche z1 le plus proche de zero tel que solution a
!c           T,P,z1 instable

!c           Pour cela on commence par cherche z1 stable
            z1 = 1.d-6
            call tpd_rp(z1,t,p00)
            if(.not.stable) then
               it = 0
               do while(.not.stable)
                  z1 = z1 / 2.d0
                  it = it + 1
                  call tpd_rp(z1,t,p00)
                  if(it.gt.1000) then
                     z1 = 1.d-12
                     stable = .true.
                     goto 1610
                  end if
               enddo
 1610          continue
            endif

!c           z1 stable trouve, on cherche z1 instable
            do while(stable)
               z1 = z1 + eps
               if(z1.ge.1.d0) goto 1598
               call tpd_rp(z1,t,p00)
            enddo

!c           On range z1 instable dans xp(1,1)
!c           A ce niveau, stable = FALSE
            xp(1,1) = z1

!c           On cherche z1 le plus proche de un tel que solution a
!c           T,P,z1 instable

!c           Pour cela on commence par cherche z1 stable

            z1 = 1.d0 - 1.d-6
            call tpd_rp(z1,t,p00)
            if(.not.stable) then
               it = 0
               do while(.not.stable)
                  z1 = z1 + (1.d0 - z1) / 2.d0
                  it = it + 1
                  call tpd_rp(z1,t,p00)
                  if(it.gt.1000) then
                     z1 = 1.d0 - 1.d-12
                     stable = .true.
                     goto 1631
                  end if
               enddo
 1631          continue
            endif

!c           z1 stable trouve, on cherche z1 instable
            do while(stable)
               z1 = z1 - eps
               if(z1.le.0.d0) goto 1598
               call tpd_rp(z1,t,p00)
            enddo

!c           On range z1 instable dans xp(1,2)
!c           A ce niveau, stable = FALSE
            xp(1,2) = z1

 1598       continue
         end if
      end if

 1599 continue

!c     Si solution stable, on sort du programme, sinon on le continue
!c      if(stable .and. iph > 0) goto 1477 ! modif du 7/10/19

      if(stable) goto 1477

      corrFL = .TRUE.
      if(xp(1,2).lt.1.d-12) then
         xp(1,2) = xp(1,1)
      end if
      if(xp(1,1).lt.1.d-12) then
         xp(1,1) = xp(1,2)
      end if

!c     x1u et y1u contiennent xp(1,1) instable et xp(1,2) instable
      eps = 1.d-3
      x1u = min(xp(1,1),xp(1,2))  ! x1 instable
      call tpd_rp(x1u,t,p00)
      if(stable) then
         do while(stable)
            x1u = x1u + eps
            if(x1u.ge.1.d0) then
               x1u = min(xp(1,1),xp(1,2))  ! x1 instable
               exit
            endif
            call tpd_rp(x1u,t,p00)
         enddo
      endif
      y1u = max(xp(1,1),xp(1,2))  ! y1 instable
      call tpd_rp(y1u,t,p00)
      if(stable) then
         do while(stable)
            y1u = y1u - eps
            if(y1u.le.0.d0) then
               y1u = max(xp(1,1),xp(1,2))  ! x1 instable
               exit
            endif
            call tpd_rp(y1u,t,p00)
         enddo
      endif

!c     recherche d'une phase stable
      stable = .false.
      z1  = x1u

      do while(.not.stable)
         z1 = z1 / 2.d0
         call tpd_rp(z1,t,p00)
      enddo

!c     x1s et y1s contiennent xp(1,1) stable et xp(1,2) stable
      x1s = z1  ! x1 stable

!c     recherche d'une phase stable
      stable = .false.
      z1     = y1u

      bleme_sup = .false.
      do while(.not.stable)
         z1 = z1 + (1.d0 - z1) / 2.d0
         if(z1.gt.(1.d0-1.d-12)) then ! Pb : le domaine est stable
!c                                       pour x1 > 0.9999999999...
            bleme_sup = .true.
            z1 = 1.d0 - 1.d-12
            goto 1830
         endif
         call tpd_rp(z1,t,p00)
      enddo
 1830 continue
      y1s = z1  ! y1 stable
!c     A ce stade, on a encadre les deux solutions que l'on cherche :
!c     x1s < xp(1,1) < x1u et y1u < xp(1,2) < y1s
!c     On met en oeuvre une dichotomie pour calculer precisement
!c     xp(1,1) et xp(1,2)

      eps = 1.d0
      it = 0
      do while(eps.gt.epsmin)
         z1old = z1
         z1 = 0.5d0 * (x1s + x1u)
         eps = dabs(z1old - z1)
         if(it.eq.0) eps = 1.d0
         call tpd_rp(z1,t,p00)
         if(stable) then
            x1s = z1
         else
            x1u = z1
         endif
         it = it + 1
      enddo

      eps = 1.d0
      it = 0
      if(bleme_sup) goto 1870
      do while(eps.gt.epsmin)
         z1old = z1
         z1 = 0.5d0 * (y1s + y1u)
         eps = dabs(z1old - z1)
         if(it.eq.0) eps = 1.d0
         call tpd_rp(z1,t,p00)
         if(stable) then
            y1s = z1
         else
            y1u = z1
         endif
         it = it + 1
      enddo

 1870 continue

      xp(1,1) = x1s
      xp(1,2) = y1s

! 1515 continue

!c     Verification de l'egalite des potentiels chimiques des
!c     solutions obtenues
!c     2 CAS SE PRESENTENT :
!c     - CAS 1 : soit xp(1,i) pas trop proche de 0 ou 1 et dans ce cas,
!c       on peut verifier l'egalite des potentiels chimiques
!c     - CAS 2 : soit ils sont vraiment proches de 0 ou 1 et dans ce cas,
!c       on perd tellement de precision dans les calculs qu'il est
!c       necessaire de choisir une solution plus robuste mais moins rapide
!c       (traitement_special = .TRUE. dans ce cas et .FALSE. sinon)

      traitement_special = .false.
      presqun = 1.d0 - 1.d-9
      if(xp(1,1).lt.1.d-9.or.xp(1,2).lt.1.d-9.or.xp(1,1).gt.presqun &
     &   .or.xp(1,2).gt.presqun) then
         traitement_special = .true.
      endif

      if(traitement_special) then
!c        on regarde s'il existe un domaine monophasique entre les
!c        compositions du diphasique
!
!c         xinf = min(xp(1,1),xp(1,2)) + 1.d0*epsmin
!c         xsup = max(xp(1,1),xp(1,2)) - 1.d0*epsmin
         xinf = min(xp(1,1),xp(1,2)) + 1.d-4
         xsup = max(xp(1,1),xp(1,2)) - 1.d-4

!c        balayage de 50 compo entre ces deux la :

         stepx = (xsup-xinf)/49.d0
         npoints = 50
         do i = 1,npoints
            z1 = xinf + dfloat(i - 1)*stepx
!         do z1 = xinf,xsup,stepx
            call tpd_rp(z1,t,p00)
            if(stable) goto 1852
         enddo

!c        si on n'a pas trouve, balayage de 1000 compo entre les 2

         stepx = (xsup-xinf)/999.d0

         do i = 1,npoints
!         do z1 = xinf,xsup,stepx
            z1 = xinf + dfloat(i - 1)*stepx
            call tpd_rp(z1,t,p00)
            if(stable) goto 1852
         enddo

!c        sinon, on a 1 diphasique classique :
         iph = 2
         goto 1477

 1852    continue

!c        a ce stade, on dispose d'une compo z1 stable entre
!c        xp(1,1) et xp(1,2)
!
!c        on recherche les phases en equilibre par dichotomie
!
!c        1er cas : on recherche la compo en eq. avec xp(1,1)
         x1s = z1
         x1u = xp(1,1)
         eps = 1.d0
         it = 0
         do while(eps.gt.epsmin)
            z1old = z1
            z1 = 0.5d0 * (x1s + x1u)
            eps = dabs(z1old - z1)
            if(it.eq.0) eps = 1.d0
            call tpd_rp(z1,t,p00)
            if(stable) then
               x1s = z1
            else
               x1u = z1
            endif
            it = it + 1
         enddo
!c         print*,'1er  eq : ',xp(1,1),z1,p00

!c        on verifie si xd(1) est compris entre z1 et xp(1,1) :
         if(z1.lt.xp(1,1)) then
            if(xd(1).lt.xp(1,1).and.xd(1).gt.z1) then
               iph = 2
               xp(1,2) = z1
               goto 302
            end if
         else
            if(xd(1).gt.xp(1,1).and.xd(1).lt.z1) then
               iph = 2
               xp(1,2) = z1
               goto 302

            end if
         endif

!c        2eme cas : on recherche la compo en eq. avec xp(1,2)
         x1s = z1
         x1u = xp(1,2)
         eps = 1.d0
         it = 0
         do while(eps.gt.epsmin)
            z1old = z1
            z1 = 0.5d0 * (x1s + x1u)
            eps = dabs(z1old - z1)
            if(it.eq.0) eps = 1.d0
            call tpd_rp(z1,t,p00)
            if(stable) then
               x1s = z1
            else
               x1u = z1
            endif
            it = it + 1
         enddo
!c         print*,'2eme eq : ',xp(1,2),z1,p00
!c        on verifie si xd(1) est compris entre z1 et xp(1,2) :
         if(z1.lt.xp(1,2)) then
            if(xd(1).lt.xp(1,2).and.xd(1).gt.z1) then
               iph = 2
               xp(1,1) = z1
               goto 302
            end if
         else
            if(xd(1).gt.xp(1,2).and.xd(1).lt.z1) then
               iph = 2
               xp(1,1) = z1
               goto 302
            end if
         endif

!c        Si on arrive ici, c'est que xd(1) n'est pas dans les
!c        2 domaines diphasiques calcules
!
!c        on verifie si le critere du plan tangent pense que
!c        la solution est monophasique stable ou pas
         z1 = xd(1)
         call tpd_rp(z1,t,p00)
         if(stable) then
            iph = 1
            xp(1,1) = xd(1)
            xp(1,2) = 0.d0
            goto 1477
         endif

!c        Autrement, ce cas demande a etre approfondi s'il est rencontre
         write(lun1,*) 'Probleme Flash. Au boulot RP2'
         write(lun2,*) 'Probleme Flash. Au boulot RP2'
!c         pause
         call prtfl
         stop
      endif

      j0    = 1
      ider  = 0
      istb  = 1

      x0(1) = xp(1,1)
      x0(2) = 1.d0 - x0(1)
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
      phi1 = dexp(lphi(1))

      x0(1) = xp(1,2)
      x0(2) = 1.d0 - x0(1)
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
      phi2 = dexp(lphi(1))

      res = dabs(xp(1,1) * phi1 - xp(1,2) * phi2)

      if(xpold(1,2).lt.1.d-12) goto 1789
      x0(1) = xpold(1,1)
      x0(2) = 1.d0 - x0(1)
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
      phi1 = dexp(lphi(1))

      x0(1) = xpold(1,2)
      x0(2) = 1.d0 - x0(1)
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
      phi2 = dexp(lphi(1))

      res2 = dabs(xpold(1,1) * phi1 - xpold(1,2) * phi2)
 1789 continue

      if(res.gt.1.d-5) then

!c        non egalite des potentiels chimiques
!c         print*,'Bleme a T = ',t,res,res2
!
!c        Dans ce cas, on recherche le y associe a xp(1,1) et le x
!c        associe a xp(1,2)
         call search_min_tpd(xp(1,1),t,p00,xmin)
         if((xp(1,1).gt.xd(1) .and. xmin.lt.xd(1)) .or. &
     &      (xp(1,1).lt.xd(1) .and. xmin.gt.xd(1))) then
             xp(1,2) = xmin
             goto 302
         end if

         call search_min_tpd(xp(1,2),t,p00,xmin)
         if((xp(1,2).gt.xd(1) .and. xmin.lt.xd(1)) .or. &
     &      (xp(1,2).lt.xd(1) .and. xmin.gt.xd(1))) then
             xp(1,1) = xmin
             goto 302
         end if

         if(dabs(xmin-xp(1,2)).lt.1.d-5) then
!c           on verifie si le critere du plan tangent pense que
!c           la solution est monophasique stable ou pas
            z1 = xd(1)
            call tpd_rp(z1,t,p00)
            if(stable) then
               iph = 1
               xp(1,1) = xd(1)
               xp(1,2) = 0.d0
               goto 1477
            endif

!c           tentative de la derniere chance
!c           on cherche manuellement les compo en equilibre avec xp(1,i)
!
!c           on commence par chercher un domaine monophasique compris
!c           entre xp(1,1) et xp(1,2)
!c           On ajoute ou retire un facteur 2*epsmin pour que les compo
!c           xinf et xsup soient dans le domaine diphasique
!            xinf = min(xp(1,1),xp(1,2)) + 1.d0*epsmin
!            xsup = max(xp(1,1),xp(1,2)) - 1.d0*epsmin
!
!c           balayage de 100 compo entre ces deux la :

            stepx = (xsup-xinf)/99.d0
            npoints = 100
!            do z1 = xinf,xsup,stepx
            do i = 1,npoints
               z1 = xinf + dfloat(i-1)*stepx
               call tpd_rp(z1,t,p00)
               if(stable) goto 1848
            enddo

!c           si on n'a pas trouve, balayage de 5000 compo entre les 2

            stepx = (xsup-xinf)/4999.d0
            npoints = 5000
!            do z1 = xinf,xsup,stepx
            do i = 1,npoints
               z1 = xinf + dfloat(i-1)*stepx
               call tpd_rp(z1,t,p00)
               if(stable) goto 1848
            enddo

!c           sinon, on abandonne :
            iph = 1
            xp(1,1) = xd(1)
            xp(1,2) = 0.d0
            goto 1477

 1848       continue

!c           a ce stade, on dispose d'une compo z1 stable entre
!c           xp(1,1) et xp(1,2)
!
!c           on recherche les phases en equilibre par dichotomie
!
!c           1er cas : on recherche la compo en eq. avec xp(1,1)
            x1s = z1
            x1u = xp(1,1)
            eps = 1.d0
            it = 0
            do while(eps.gt.epsmin)
               z1old = z1
               z1 = 0.5d0 * (x1s + x1u)
               eps = dabs(z1old - z1)
               if(it.eq.0) eps = 1.d0
               call tpd_rp(z1,t,p00)
               if(stable) then
                  x1s = z1
               else
                  x1u = z1
               endif
               it = it + 1
            enddo
!c            print*,'1er  eq : ',xp(1,1),z1,p00

!c           on verifie si xd(1) est compris entre z1 et xp(1,1) :
            if(z1.lt.xp(1,1)) then
               if(xd(1).lt.xp(1,1).and.xd(1).gt.z1) then
                  iph = 2
                  xp(1,2) = z1
                  goto 302
               end if
            else
               if(xd(1).gt.xp(1,1).and.xd(1).lt.z1) then
                  iph = 2
                  xp(1,2) = z1
                  goto 302
               end if
            endif

!c           2eme cas : on recherche la compo en eq. avec xp(1,2)
            x1s = z1
            x1u = xp(1,2)
            eps = 1.d0
            it = 0
            do while(eps.gt.epsmin)
               z1old = z1
               z1 = 0.5d0 * (x1s + x1u)
               eps = dabs(z1old - z1)
               if(it.eq.0) eps = 1.d0
               call tpd_rp(z1,t,p00)
               if(stable) then
                  x1s = z1
               else
                  x1u = z1
               endif
               it = it + 1
            enddo
!c            print*,'2eme eq : ',xp(1,2),z1,p00
!c           on verifie si xd(1) est compris entre z1 et xp(1,2) :
            if(z1.lt.xp(1,2)) then
               if(xd(1).lt.xp(1,2).and.xd(1).gt.z1) then
                  iph = 2
                  xp(1,1) = z1
                  goto 302
               end if
            else
               if(xd(1).gt.xp(1,2).and.xd(1).lt.z1) then
                  iph = 2
                  xp(1,1) = z1
                  goto 302
               end if
            endif

!c           Si on arrive ici, c'est que xd(1) n'est pas dans les
!c           2 domaines diphasiques calcules
!c           Ce cas demande a etre approfondi s'il est rencontre
            write(lun1,*) 'Probleme Flash. Au boulot RP'
            write(lun2,*) 'Probleme Flash. Au boulot RP'
!c            pause
            stop
         end if


         if(dabs(xmin-xp(1,2)).lt.1.d-3) then
            write(*,*) 'Resultat incertain, T/K = ',t,'P/bar = ',p00
            iph = 1
            xp(1,1) = xd(1)
            xp(1,2) = 0.d0
            goto 1477
         end if

!c         print*,'GROS SOUCI'
!c         print*,'xp(1,1) , xp(1,2) : ',xp(1,1),xp(1,2)
!c         print*,'p00 = ',p00
!c         print*,'tk  = ',t
!
!c         xp = xpold
!c         iph = 0
         iph = 2
         goto 1477
      end if

 302  continue

      iph = 2

!c     tri des solutions
      if(croissant.and.xp(1,1).gt.xp(1,2)) then
         tempo = xp(1,1)
         xp(1,1) = xp(1,2)
         xp(1,2) = tempo
      end if
      if(.not.croissant.and.xp(1,2).gt.xp(1,1)) then
         tempo = xp(1,1)
         xp(1,1) = xp(1,2)
         xp(1,2) = tempo
      end if

!c      print*,'P/bar = ',p00,'T/K = ',t0

 1477 continue
      if (iph.le.1) then
         xp(2,1) = 1.0d+0 - xp(1,1)
         x0(1)  = xp(1,1)
         x0(2)  = xp(2,1)
         j0 = 1
         ider = 0
         istb = 1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,gt0)

         etap(1) = eta0
         etap(2) = 0.0d+00
         vp(1) = v0
         vp(2) = 0.0d+00
         gp(1) = gt0
         gp(2) = 0.0d+00
         gt = gt0
         npt(1) = ntt
         npt(2) = 0.0d+00
         densp(2) = 0.0d+00
         zp(2) = 0.0d+00
!c          vmolp(2) = 0.0d+00
         do i = 1,nco
            xp(i,1) = x0(i)
            xp(i,2) = 0.0d+00
            np(i,1) = nt(i)
            np(i,2) = 0.0d+00
         end do
      else !iph = 2
         xp(2,1) = 1.0d+0 - xp(1,1)
         x0(1)  = xp(1,1)
         x0(2)  = xp(2,1)
         j0 = 1
         ider = 0
         istb = 1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,gt0)
         lnphi_L(1:nco) = lphi(1:nco)
!         print*,'ln phi L = ',lphi(1:nco)
         etap(1) = eta0
         vp(1) = v0
         gp(1) = gt0
         npt(1) = (nt(1) - ntt * xp(1,2))/(xp(1,1) - xp(1,2))
!c
         xp(2,2) = 1.0d+0 - xp(1,2)
         x0(1)  = xp(1,2)
         x0(2)  = xp(2,2)
         j0 = 1
         ider = 0
         istb = 1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,gt0)
         lnphi_V(1:nco) = lphi(1:nco)
!         print*,'ln phi V = ',lphi(1:nco)
         etap(2) = eta0
!c
!c        pour les flash tres tres proches des corps purs, il peut y avoir
!c        un probleme (eta_liq = eta_gaz)
!c
         if((xp(1,1).gt.0.999.and.xp(1,2).gt.0.999).or. &
     &      (xp(2,1).gt.0.999.and.xp(2,2).gt.0.999)) then !corps pur
            if(dabs(etap(2)-etap(1)).lt.0.01) then
!c
!c              phas1 a choisi 2 fois la meme phase
!c
               x0(1)  = xp(1,2)
               x0(2)  = xp(2,2)
               j0 = 1
               ider = 0
               istb = 10
               call phas1(j0,ider,istb,p00,x0,eta0,v0,gt0)
               etap(2) = eta0
            endif
         endif
         vp(2) = v0
         gp(2) = gt0
         npt(2) = ntt - npt(1)

         do j=1,iph
            mmolp(j) = 0.0d+00
            vph(j) = vp(j)
            do i=1,nco
               mmolp(j) = mmolp(j) + xp(i,j)*mmol(i)
            end do
            densp(j) = mmolp(j)/vph(j)
         end do

!c         if(etap(2).gt.etap(1)) then
         if(densp(2).gt.densp(1)) then
            xxx = npt(1)
            npt(1) = npt(2)
            npt(2) = xxx
            xxx = etap(1)
            etap(1) = etap(2)
            etap(2) = xxx
            xxx = vp(1)
            vp(1) = vp(2)
            vp(2) = xxx
            xxx = gp(1)
            gp(1) = gp(2)
            gp(2) = xxx
            xxx = mmolp(1)
            mmolp(1) = mmolp(2)
            mmolp(2) = xxx
            xxx = zp(1)
            zp(1) = zp(2)
            zp(2) = xxx
            do i=1,nco
               xxx = xp(i,1)
               xp(i,1) = xp(i,2)
               xp(i,2) = xxx
               xxx = np(i,1)
               np(i,1) = np(i,2)
               np(i,2) = xxx
            end do
         endif
      endif

       ph = iph
       if(iph.eq.0) ph = 1
       rti = 1.0d+0/r/t0
       vt = 0.0d+00
       mmolt = 0.0d+00
       do j=1,ph
          mmolp(j) = 0.0d+00
          do i=1,nco
             mmolp(j) = mmolp(j) + xp(i,j)*mmol(i)
          end do
          mmolt = mmolt + npt(j)*mmolp(j)
          zp(j) = p0*vp(j)*rti
          densp(j) = mmolp(j)/vp(j)
          vp(j) = npt(j)*vp(j)
          vt = vt + vp(j)
       end do
       denst = mmolt/vt
       mmolt = mmolt/ntt
       zt = p0*vt*rti/ntt

       do i=1,nco
          volat(i) = 0.0D+00
       end do

       if(iph.eq.2) then
          dropout = 1.00d+02*vp(1)/vt   ! % volumique de la phase liquide
          dropout = dabs(dropout)
          do i=1,nco
             volat(i) = xp(i,2)/xp(i,1)   !y/x
          end do
       endif
       fl_ok = .false.
       if(iph.eq.2) then
          if(ygtx.and.xp(1,2).gt.xp(1,1)) fl_ok = .true.
          if(.not.ygtx.and.xp(1,1).gt.xp(1,2)) fl_ok = .true.
       endif
!c
! 1483 continue

!c
!c     on teste la barotropie !
!c
      densmax = max(densp(1),densp(2))
      densmin = min(densp(1),densp(2))
      ecart_densite = 100.0d+00*(dabs(densmax-densmin)/densmin)

      if(iph.eq.2.and.(nom(1)(1:3).eq.'H2O' &
     &    .or.nom(2)(1:3).eq.'H2O')) then

         call calcul_psat(1,t0,psat1,er1,vliq,vgaz)
         if(er1.ne.0) then
            rien = 7.0d+0/3.0d+00*(fac(1)+1.0d+0)
            psat1 = pc(1)*10.0d+00**(rien*(1.0d+0-tc(1)/t0))
         endif

         call calcul_psat(2,t0,psat2,er2,vliq,vgaz)
         if(er2.ne.0) then
            rien = 7.0d+0/3.0d+00*(fac(2)+1.0d+0)
            psat2 = pc(2)*10.0d+00**(rien*(1.0d+0-tc(2)/t0))
         endif

         Ps_max = max(Psat1,Psat2)
         Ps_min = min(Psat1,Psat2)

         if(p00.gt.1.0d+0*ps_min) then
!c
!c           on a un systeme renfermant de l'eau.
!c
!c           on est tres probablement en ELL et il y a risque de barotropie :
!c           les densites des 2 phases liquides s'inversent.
!c
!c           1 cas typique est le systeme : eau + hexane. A basse pression,
!c           la phase aqueuse est plus dense que la phase organique. Vers
!c           2000 bar, la phase orga devient plus dense que la phase aqueuse.
!c
!c           Le flash classe les phases en fonction de leur densite. Il y aura
!c           donc inversion des phases (la phase aqueuse initialement classe
!c           "phase 1", devient "phase 2"). C'est genant ! DONC :
!c
!c           la phase dite "liquide" i.e. phase 1 sera TOUJOURS celle riche en
!c           compose lourd (constituant 2).
!c           la phase dite "vapeur" i.e. phase 2 sera la phase riche en compose
!c           leger (constituant 1)
!c
            on_inverse = .false.

             if(nom(2)(1:3).eq.'H2O') then !l'eau est le constituant 2 (le + lourd)
               if(p00.gt.(1.0d+0*Ps_max+0.0d+0*Ps_min)) then  !! if(p00.gt.(1.0d+0*Ps_min+2.0d+0*Ps_min))
                if(xp(1,1).gt.xp(1,2)) then
                 if(ecart_densite.lt.200.0d+00) then
                   on_inverse = .true.
                 endif
                endif
               endif
             endif
             if(nom(1)(1:3).eq.'H2O') then!l'eau est le constituant 1 (le + leger)
               if(p00.gt.(1.0d+0*Ps_max+0.0d+0*Ps_min)) then
                if(xp(1,1).lt.xp(1,2)) then
                 if(ecart_densite.lt.200.0d+00) then  !AVANT ECART_DENSITE.LT.30
                   on_inverse = .true.
                 endif
                endif
               endif

               if(max(xp(1,1),xp(1,2)).lt.0.9.and.  & !AVANT  Min(XP11,XP12).LT.0.9
     &             p00.lt.(2*Ps_min+1.5*Ps_max).and. &
     &             abs(xp(1,1)-xp(1,2)).lt.0.5) then   !AVANT P00.LT.1.5*PS_MAX H2O-C6
                   on_inverse = .false.
               endif

             endif

!C           on_inverse = .FALSE.
            if(on_inverse) then
               vp(1) = vp(1)/npt(1)
               vp(2) = vp(2)/npt(2)
               xxx = npt(1)
               npt(1) = npt(2)
               npt(2) = xxx
               xxx = etap(1)
               etap(1) = etap(2)
               etap(2) = xxx
               xxx = vp(1)
               vp(1) = vp(2)
               vp(2) = xxx
               xxx = gp(1)
               gp(1) = gp(2)
               gp(2) = xxx
               xxx = mmolp(1)
               mmolp(1) = mmolp(2)
               mmolp(2) = xxx
               xxx = zp(1)
               zp(1) = zp(2)
               zp(2) = xxx
               do i=1,nco
                  xxx = xp(i,1)
                  xp(i,1) = xp(i,2)
                  xp(i,2) = xxx
                  xxx = np(i,1)
                  np(i,1) = np(i,2)
                  np(i,2) = xxx
               end do
               ph = iph
               rti = 1.0d+0/r/t0
               vt = 0.0d+00
               mmolt = 0.0d+00
               do j=1,ph
                  mmolp(j) = 0.0d+00
                  do i=1,nco
                     mmolp(j) = mmolp(j) + xp(i,j)*mmol(i)
                  end do
                  mmolt = mmolt + npt(j)*mmolp(j)
                  zp(j) = p0*vp(j)*rti
                  densp(j) = mmolp(j)/vp(j)
                  vp(j) = npt(j)*vp(j)
                  vt = vt + vp(j)
               end do
               denst = mmolt/vt
               mmolt = mmolt/ntt
               zt = p0*vt*rti/ntt

               do i=1,nco
                  volat(i) = 0.0D+00
               end do

               if(iph.eq.2) then
                  dropout = 1.00d+02*vp(1)/vt   ! % volumique de la phase liquide
                  dropout = dabs(dropout)
                  do i=1,nco
                     volat(i) = xp(i,2)/xp(i,1)   !y/x
                  end do
               endif
               fl_ok = .false.
               if(iph.eq.2) then
                  if(ygtx.and.xp(1,2).gt.xp(1,1)) fl_ok = .true.
                  if(.not.ygtx.and.xp(1,1).gt.xp(1,2)) fl_ok = .true.
               endif
            endif
         endif
      endif
      if(imp.ge.1) then
         call prtfl
      endif
      return
      end subroutine flash

end module
!==============================================================================
! END: MODULE SECTION
!==============================================================================

!c*****************************************************************************
!   SETUP SUBROUTINES/FUNCTIONS
!c*****************************************************************************
      subroutine setup(PGLinputDir,err0)
		USE GlobConst, only:iEosOpt
        use comflash
        implicit none
        integer :: err0,i,iGEoption
        real(8) :: etac,zc
        character*255 :: fileName,PGLinputDir
!C--------------------------------------------------------------------
        err0 = 0
!C      Load model definition and parameters
!c      gE model          : 1 = van Laar, 2 = Wilson, 3 = UNIQUAC, 4 = NRTL
		iGEoption=2
		if(iEosOpt==21)iGEoption=1
        call read_option(err0)
        if(err0.ne.0) return
!C      EoS copmponent-independent parameters
        etac = ((1.d0-r1)*(1.d0-r2)**2)**(1.d0/3.d0) + &
     &         ((1.d0-r2)*(1.d0-r1)**2)**(1.d0/3.d0) + 1.d0
        etac = 1.d0 / etac

        zc = 1.d0/(3.d0 - etac*(1.d0+r1+r2))	!r1= -(1+sqrt2), r2 = sqrt2-1 => r1+r2= -1-sqrt2+sqrt2-1 = -2

        omegab = zc*etac

        omegaa = (1.d0-etac*r1) * (1.d0-etac*r2) * (2.d0-etac*(r1+r2)) &
     &          / ((1.d0-etac)*(3.d0-etac*(1.d0+r1+r2))**2)

!C      Load Gauss coefficients for Gauss_Legendre integration
        err0 = 1
        fileName = "gauss_legendre.dat"
        fileName =TRIM(PGLinputDir)//"\PrLorraineEtcFiles\" // trim(fileName)
        OPEN(20,FILE=fileName,status='old',err=10)
        err0 = 0
        READ(20,*)
        READ(20,*)
        READ(20,*)

        DO i = 1,5
           READ(20,*) z_gauss(i),w_gauss(i)
        ENDDO

        CLOSE(20)

        z_gauss(6:10) = -z_gauss(1:5)
        w_gauss(6:10) = w_gauss(1:5)

  10    continue
       if(err0.eq.1)then
           write(lun1,*)'gauss_legendre.dat file not found'
           write(errm,*)'gauss_legendre.dat file not found'
       end if
       return


      end subroutine !FLASH
!c
!c ----------------------------------------------------------------------------
!c
      subroutine setfluid(PRLinputDir,fluids,err0)
        use comflash
        implicit none
        integer :: i,j,k,id1,id2,ioErr ! JRE 20210604 - for case when BIPs not found.
        integer :: err0,fluids(nco)
        integer :: prop(n_prop),typec
        integer,parameter :: n_pts = 10
        real(8) :: rtsp,calc_c,BIP(2)
        character*234 :: fluidFile,PRLinputDir
        character(70) :: str

        do i = 1,nco
            write(str,'(i5)') fluids(i)
            fluidFile = 'CHEM' // trim(adjustl(str))//'.fld'
            fluidFile = TRIM(PRLinputDir)//"\PrLorraineEtcFiles\fluids\" //  &
     &                  trim(fluidFile)
            prop = 0
            err0 = 1
            open(22,file=fluidFile,status='old',err=10)
            err0 = 0
            read(22,*)
            do k = 1,nb_ids
               read(22,'(A70)') idComp(i,k)
            end do

            read(22,*)
            read(22,*)
!C      store fluid constants for each of the nco components
            read(22,*) mmol(i) !wmas(i)
            read(22,*) tc(i)
            read(22,*) pc(i)
            read(22,*) fac(i)
            read(22,*) vc(i)
            read(22,*)
            read(22,*) tpt(i)
            read(22,*) tpp(i)
            read(22,*) tnbp(i)

            if(choix_r1r2.eq.1)then !tc-RK parameters
                read(22,*)
                read(22,*)
                read(22,*)
                read(22,*) parL(i),parM(i),parN(i)
                read(22,*)
                read(22,*) parC(i)
                do j = 1,6
                   read(22,*)
                enddo
            elseif(choix_r1r2.eq.2)then !tc-PR parameters
                do j = 1,9
                   read(22,*)
                enddo
                read(22,*) parL(i),parM(i),parN(i)
                read(22,*)
                read(22,*) parC(i)
            end if

            read(22,*)
            read(22,*)
            !Temperature-dependent properties

            do k = 1,nb_tpty
                read(22,*)
                read(22,*)eqIDComp(i,k), &
     &                (tdep_prop(i,k,j), j=1,nb_coeff), &
     &                 tmin_prop(i,k),tmax_prop(i,k)
                if(eqIDComp(i,k).ne.888) prop(k) = 1
            end do
          close(22)

          if((tc(i) > 2000d0) .or. &
     &        (pc(i) < 8889d0 .and. pc(i) > 8887d0)  &
     &       .or. (fac(i) < 8889d0 .and. fac(i) > 8887d0)) then
              write(lun1,1502)trim(str)
              write(errm,1502)trim(str)
 1502         format(' Molecule ',(A),' cannot be modeled by a CEoS.')
              err0 = 104
              return
          end if


!C      set the type of alpha function and volume translation parameters
          if(FoncAlpha.eq.'TWU88' .or. (FoncAlpha.eq.'TWU91' .and.   &
     &         parL(i).gt.8887d0))then

             parL(i) = 0.0728d0 + 0.6693d0*fac(i) + 0.0925d0*fac(i)**2
             parM(i) = 0.8788d0 - 0.2258d0*fac(i) + 0.1695d0*fac(i)**2
             parN(i) = 2.0d0
          end if

          rtsp = r*tc(i)/pc(i) ! pc en bar, tc en K et r = 83.14 bar.cm3/mol/K
          b(i) = omegab*rtsp ! cm3/mol
          ac(i) = omegaa*rtsp*r*tc(i) ! bar.cm6/mol2
          ac(i) = ac(i)*1.d-7 ! unité du systeme international (Pa.m6/mol2)

!c         Calcul de bc(i) = b(i)/R en uSI = K/Pa
!c          r = 83.14 bar.cm3/mol/K et b en cm3/mol
          bc(i) = 1.d-5 * b(i) / r ! uSI

          if(FoncAlpha.EQ.'SOAVE') then
!c         -------------------------
!c         Fonction alpha de SOAVE :
!c         -------------------------
              if(choix_r1r2.eq.1) then
!c            EoS Soave72 :
!c           -----------
                  msoave(i) = 0.480d0 + 1.574d0*fac(i) - &
     &                                0.176d0*fac(i)**2
              elseif(choix_r1r2.eq.2)then
!c            EoS PR78 :
!c           -----------
                 if(fac(i).le.0.491d+00) then
                    msoave(i) = 0.37464d+00 + 1.54226d+00*fac(i) &
     &                    - 0.26992d+00*fac(i)**2
                 else
                     msoave(i) = 0.379642d+00 + 1.48503d+00*fac(i) &
     &                     - 0.164423d+00*fac(i)**2 &
     &                     + 0.016666d+00*fac(i)**3
                  endif
               end if

!c         -------------------------
!c         -------------------------
            endif

           if(parC(i).gt.8887d0)then !c not available from tc-models
               if(prop(ldn_idx).eq.1 .and. &
     &            (0.8d0*tc(i).gt.tmin_prop(i,ldn_idx) .and.   &
     &             0.8d0*tc(i).lt.tmax_prop(i,ldn_idx)))then !c from exp data at Tr=0.8
                  typec = 1
               else !c from correlation c=f(omega)
                  typec = 2
               end if
               parC(i) = calc_c(i,typec)
           endif
           parC(i) = parC(i) * 1.0d-06 ! cm3/mol > m3/mol
        end do
        c(1:nco) = parC(1:nco) * 1.0d06 !m3/mol ---> cm3/mol
!c      Modèle Wilson / Uniquac
       IF(choix_ge == 2) THEN ! Wilson
          mu(1:nco) = 1.d0
          sigma(1:nco) = b(1:nco) - parC(1:nco)*1.0d06 ! cm3/mol
       ELSEIF(choix_ge == 3) THEN !UNIQUAC
          mu(1:ncomax) = qi(1:ncomax) !q UNIQUAC parameter
          sigma(1:ncomax) = qi(1:ncomax) !q UNIQUAC parameter
       ENDIF

!c    Load BIP
      agE0 = 0.0d0	! JRE 20210604 initialize to zero. If nothing found, so be it.
      bgE0 = 0.0d0
      cgE0 = 0.0d0

      gE_file = TRIM(PRLinputDir)//'\PrLorraineEtcFiles\' // trim(gE_file)
      open(5001,file=(gE_file),status='old',err=10)
      do i = 1,nco-1
        do j = i+1,nco
           do
             read(5001,*,ioStat=ioErr) id1,id2,BIP(1),BIP(2)
			 if(ioErr/=0)exit	! JRE 20210604 If end of file is reached and nothing found, so be it: age0=0
             if((fluids(i).eq.id1 .or. fluids(i).eq.id2) .and. &
     &          (fluids(j).eq.id1 .or. fluids(j).eq.id2)) then
                if(fluids(i).eq.id1)then
                    agE0(i,j) = BIP(1)
                    agE0(j,i) = BIP(2)
                else
                    agE0(i,j) = BIP(2)
                    agE0(j,i) = BIP(1)
                end if
                exit
             end if
           end do
           rewind(5001)
        end do

      end do
      CLOSE(5001)

 10    continue
       if(err0.eq.1)then
           write(lun1,1501)trim(fluidFile),trim(str)
           write(errm,1501)trim(fluidFile),trim(str)
 1501      format('The file ',(A),' does not exist.', &
     &      ' Molecule ',(A),' is not present in the database.')
       end if
       return
      end subroutine
!c
!c ----------------------------------------------------------------------------
!c
     subroutine read_option(err0)
	  USE GlobConst, only:iEosOpt
      use comflash
      implicit none
      integer :: err0
!c ----------------------------------------------------------------------------
      err0 = 0
!c    MODEL DEFINITION
!c    Peng-Robinson EoS
      choix_r1r2 = 2
      r1 = -1.d0 - DSQRT(2.d0)
      r2 = -1.d0 + DSQRT(2.d0)
      choice(1) = 2		!
!	Default is tcPRvtTwuWilson
!c    Alpha function
      FoncAlpha = 'TWU91'
      choice(2) = 3	!Designates the Twu alpha
!c    Quadratic mixing rule for b: b_ij = {0.5 * [b_i^(1/s) + b_j^(1/s)]}^s * (1 - L_ij)
      choix_b = 2
      quadra_b = .true.
!c   - Value of  parameter s   (for b quadratic)
      puisb = 1.5d0
!c   - Value of  parameter L12 (for b quadratic)
      Lij(1,2) = 0.0d0
      Lij(2,1) = Lij(1,2)
!c    gE model          : 1 = van Laar, 2 = Wilson, 3 = UNIQUAC, 4 = NRTL
	  choix_ge=2
	  if(iEosOpt==21)then	! Implement the PPR78 options.
		FoncAlpha = 'SOAVE'
		choice(2) = 1
		choix_b = 1 
		quadra_b = .FALSE.
		puisb = 1
		choix_ge = 1
	  endif

      if(choix_ge.gt.4 .or. choix_ge.lt.1) THEN
         err0 = 5
         write(lun1,'(1X,A)') 'Wrong choice for gE (OPTION.D)'
         write(errm,'(1X,A)') 'Wrong choice for gE (OPTION.D)'
         close(541)
         return
      endif

      if(choix_ge.eq.1) then ! van Laar
         modele_gE = "PPR78 (Van Laar)"
      ELSEIF(choix_ge.eq.2 .or. choix_ge .eq.3) THEN ! Wilson/UNIQUAC
         modele_gE = "HV (Uniquac-Wilson)"
      ELSEIF(choix_ge.eq.4) THEN ! NRTL
         modele_gE = "HV (NRTL)"
      ENDIF

      if(choix_ge.ge.2) then ! WILSON/UNIQUAC/NRTL
         Lambda = -0.6232d0
      else
         Lambda = 0.d0
      endif

      kij_const = .false.
      predictif = .false.

!c     --- MODEL: WILSON / UNIQUAC / NRTL ---
!c     "   (2.1) SIMULATION"
      simul_PRWil = .TRUE.
      if(choix_ge.eq.2) gE_file = 'PARAM_WILSON.txt'
      if(choix_ge.eq.3) gE_file = 'PARAM_UNIQUAC.txt'
      if(choix_ge.eq.4) gE_file = 'PARAM_NRTL.txt'

      if(modele_gE == "HV (NRTL)") alphgE = 0.3d0 !By default

      ph3 = 0
      nimp = 0
      iac = 1
      isubs = 0
      cam = 0
      epseq = 1.0d-14
      eps0 = 1.0d-14
      eps1 = 1.0d-03

      if(isubs.eq.1.and.eps1.gt.1d-6) eps1 = 1d-6

      oa = 2
      nb = 2
      l0 = 0.5d0
      nlp = 1
      nls = 3
      ld(1:3) = (/0.5d0,0.01d0,1.0d0/)
      nitflash = 20
      choix_flash = 2

      return
      end subroutine
!c
!c ----------------------------------------------------------------------------
!c
!c*****************************************************************************
!   END: SETUP SUBROUTINES/FUNCTIONS
!c*****************************************************************************
!c
!c
!c*****************************************************************************
!   PURE COMPONENTS SUBROUTINES/FUNCTIONS
!c*****************************************************************************
      subroutine calcul_psat(numero,tk,psat,err0,vol_mol_liq, &
     &          vol_mol_vap)

!c Calcul de la pression de vapeur saturante d'un corps pur
!c par une EoS cubique.
!c
!c INPUT : numero, entier = 1 ou 2 pour spécifier le constituant
!c         tk (temp. specifiee) en K
!c OUTPUT :  psat (bar), err0 (variable d'erreur, = 0 si tout va bien)
!c
      use comflash
      implicit none
      real*8 apr,bpr,a_function
      real*8 tcrit,pcrit
      real*8 psat,psat0,tk,vol_mol_liq,vol_mol_vap
      real*8 lnphi_liq,lnphi_vap,lnphi
      real*8 res1,res2,res3,vvap_old,vliq_old,eps
!c     dpsat          : d Psat(1) / dT
!c     dhvap          : enthalpie de vaporisation
!c      real*8 dhvap,dpsat,h_res
      real*8 dpsat
      logical der_calc

      integer nb000,nb_iter,nb_iter_max,err0,err1,numero

      common/ab_pr/apr,bpr
      common/deriv/der_calc

      common/der_psat/dpsat

      eps = 1.d-10
      psat = 0.d0
      err0 = 0
      nb_iter_max = 150
      nb_iter = 0  ! compteur d'iterations

      bpr = b(numero) * 1.d-6 ! cm3/mol > m3/mol

      tcrit = tc(numero)        !en K
      pcrit = pc(numero)*1.0d+5 !en Pa

      if(dabs(tcrit - tk) / tcrit.lt.1.d-5) then
         psat = pcrit
         goto 31
      endif

      if(tk.gt.tcrit) then
         err0 = 1  ! T > Tc et Psat non calculable
         goto 31
      endif

!c     procedure d'initialisation psat :
!c     -------------------------------
      call initialisation_psat(numero,tk,psat,err0)
      if(err0.ne.0) goto 31

!c     Methode de Newton pour calculer psat :
!c     ------------------------------------
!c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,tk,0)
 30   continue
      vliq_old = vol_mol_liq
      vvap_old = vol_mol_vap
      call resol_eos(tk,psat,vol_mol_liq,vol_mol_vap,nb000,err1)
      if(err1.ne.0) err0 = 10
      if(err0.ne.0) goto 31

      lnphi_liq = lnphi(tk,vol_mol_liq)
      lnphi_vap = lnphi(tk,vol_mol_vap)

      psat0 = psat

      if(nb_iter.eq.nb_iter_max) then
         if(tk/tcrit.lt.0.35d0) then
            eps = 5.d-4
         else
            eps = 1.d-6
         endif

         if(res1.lt.eps .and. res2.lt.eps .and. res3.lt.eps) goto 31

         if(tk/tcrit.lt.0.5d0 .and. psat.lt.1.d0) then
            goto 31 ! pas d'erreur dans ce cas
         endif

         err0 = 2
         goto 31
      endif
      nb_iter = nb_iter + 1

      psat = psat0 - rgas * tk * (lnphi_vap - lnphi_liq) &
     &        / (vol_mol_vap - vol_mol_liq)

!c     si la formule de Newton conduit a des pression < 0 :
      if(psat.lt.1.d-16) then
         psat = psat0 * dexp(-rgas*tk/psat0*(lnphi_vap-lnphi_liq) &
     &                       / (vol_mol_vap-vol_mol_liq))
      endif

!c     remarque : a basse T et faible Ps, ces trois criteres sont
!c                identiques et res1 = res2 = res3 car :
!c                d Delta_vap ln phi = dPsat / Psat = d vGAZ / vGAZ
!c                NB : d vLIQ / vLIQ = 0
      res1 = dabs(lnphi_vap - lnphi_liq)
      res2 = dabs(psat - psat0)/dabs(psat0)

      if(nb_iter>1) then
         res3 = dabs(vol_mol_liq - vliq_old)/dabs(vliq_old) &
     &          + dabs(vol_mol_vap - vvap_old)/dabs(vvap_old)
      endif

      if(res1.lt.eps .and. res2.lt.eps .and. res3.lt.eps &
     &   .and. nb_iter>1) then
         goto 31
      else
         goto 30
      endif

 31   continue

      psat = psat * 1.d-5 ! en bar

      end subroutine calcul_psat
!c
!c ----------------------------------------------------------------------------
!c
      subroutine initialisation_psat(numero,tk,psat,err0)

!c Initialisation de la pression de vap. sat. :
!c   1) initialisation a partir d'une formule l'approximant connaissant
!c      tc, pc et omega.
!c   2) Affinage de cette valeur : on calcule psat de sorte que la
!c      resolution de l'equation d'etat conduise a une racine liquide et
!c      une racine vapeur
      use comflash
      implicit none
      real*8 tcrit,pcrit,omega0,a_function
      real*8 apr,bpr

      real*8 psat,tk,vol_mol_liq,vol_mol_vap,pinf,psup
      integer nb0,nb_iter,nb_iter_max,err0,numero

      common/ab_pr/apr,bpr

      nb_iter = 0
      nb_iter_max = 1000
      err0 = 0

      pcrit = pc(numero) * 1.d5 ! bar > Pa
      tcrit = tc(numero) ! K
      omega0 = fac(numero)

      psat = pcrit * 10.d0**(7.d0 / 3.d0 * (omega0 + 1.d0) &
     &                    * (1.d0 - tcrit / tk))

!c     La resolution de l'EoS a T=tk et P=psat amene-t-elle 3 racines ?
!c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,tk,0)

 50   continue
      call resol_eos(tk,psat,vol_mol_liq,vol_mol_vap,nb0,err0)
      if(err0.ne.0) goto 52

      if(nb_iter.eq.nb_iter_max) then
         err0 = 2
         goto 52
      endif
      nb_iter = nb_iter + 1

!c     si la racine liq a ete trouvee mais pas la racine vap :
!c     -----------------------------------------------------
      if(vol_mol_vap.gt.1.d25) then
!c        diminuer la valeur de psat jusqu'a ce que 2 racines apparaissent
         psat = 2.d0 * psat - pcrit
         nb_iter = nb_iter + 1
         goto 50
      endif

!c     si la racine vap a ete trouvee mais pas la racine liq :
!c     -----------------------------------------------------
      if(vol_mol_liq.gt.1.d25) then
!c        recherche de la valeur de psat conduisant a 2 racines par
!c        dichotomie
         pinf = psat
         psup = pcrit
         psat = 0.5d0 * (pinf + psup)
         nb_iter = nb_iter + 1

 51      continue

!c        Resolution de l'equation d'etat a psat et tk donnees
         call resol_eos(tk,psat,vol_mol_liq,vol_mol_vap,nb0,err0)
         if(err0.ne.0) goto 52

         if(nb_iter.eq.nb_iter_max) then
            err0 = 2
            goto 52
         endif
         nb_iter = nb_iter + 1

!c        si la racine liquide a ete trouvee mais pas la racine vapeur :
         if(vol_mol_vap.gt.1.d25) then
            psup = psat
            psat = 0.5d0 * (pinf + psup)
            nb_iter = nb_iter + 1
            goto 51
         endif

!c        si la racine vapeur a ete trouvee mais pas la racine liquide :
         if(vol_mol_liq.gt.1.d25) then
            pinf = psat
            psat = 0.5d0 * (pinf + psup)
            nb_iter = nb_iter + 1
            goto 51
         endif
      endif

 52   continue

      end subroutine initialisation_psat
!c
!c ----------------------------------------------------------------------------
!c
      subroutine resol_eos(tk,pres0,vol_mol_liq,vol_mol_vap,nb0,err3)

!c sous-programme resolvant a temperature tk et pression pres0 fixees,
!c l'EoS d'un corps pur. Il calcule les racines liquide et vapeur
!c de l'equation si elles existent. Si elles n'existent pas, elles
!c valent 1.d31
      use comflash
      implicit none
      real*8 apr,bpr
      real*8 pres0,tk,eta0,vol_mol_liq,vol_mol_vap,pression,pEoS
      real*8 pi,rho(3),coefb0,coefc0,coefd0,coefb,coefc,coefd
      real*8 qq,rr,delta,ss0,tt,theta,vm(3),res,cur
      real*8 eta_c
      integer nb_liq,nb_vap,nb0,err3,i
      character*3 state

      common/ab_pr/apr,bpr
      err3 = 0
!c     Compacite critique :
      eta_c = cur((1.d0-r1)*(1.d0-r2)**2)+cur((1.d0-r2)*(1.d0-r1)**2) &
     &        +1.d0
      eta_c = 1.d0 / eta_c

      if(pres0.lt.0.d0) then
         write(*,*) 'p < 0 dans RESOL_EOS corps pur -- STOP'
         stop
      endif

!c     ***********************************
!c     RESOLUTION PAR LA METHODE DE CARDAN
!c     ***********************************

      pi = 4.d0*datan(1.d0)
      nb0 = 0
      rho(:) = 0.d0

!c     Resolution de l'equation cubique en rho=1/v
      coefb0 = -(bpr*(r1+r2+1.d0) + rgas*tk/pres0)
      coefc0 = bpr**2*(r1*r2+r1+r2)+rgas*tk*bpr*(r1+r2)/pres0+apr/pres0
      coefd0 = -bpr*((r1*r2*bpr**2+r1*r2*bpr*rgas*tk/pres0+apr/pres0))
      coefb = coefc0/coefd0
      coefc = coefb0/coefd0
      coefd = 1.d0/coefd0

      qq    = (3.d0*coefc-coefb**2)/9.d0
      rr    = (9.d0*coefb*coefc-27.d0*coefd-2.d0*coefb**3.d0)/54.d0
      delta = qq**3 + rr**2
      vm(:) = 0.d0
      if(delta.gt.0.d0) then
         ss0 = cur(rr+dsqrt(delta))
         tt = cur(rr-dsqrt(delta))
         rho(1) = ss0 + tt - coefb/3.d0
         rho(2) = 0.0d0
         rho(3) = 0.0d0
         vm(1) = 1.d0/rho(1)
         nb0 = 1
      else
         theta = acos(rr/dsqrt(-qq**3))
         rho(1) = 2.d0*dsqrt(-qq)*dcos(theta/3.d0)-coefb/3.d0
         rho(2) = 2.d0*dsqrt(-qq)*dcos(theta/3.d0+2.d0*pi/3.d0)- &
     &               coefb/3.d0
         rho(3) = 2.d0*dsqrt(-qq)*dcos(theta/3.d0+4.d0*pi/3.d0)- &
     &               coefb/3.d0
         rho(:) = 1.d0/rho(:)
         vm(1) = minval(rho)
         vm(3) = maxval(rho)
         do i = 1,3
            if(rho(i).ne.vm(1) .and. rho(i).ne.vm(3)) then
               vm(2) = rho(i)
            endif
         enddo
         rho(:) = 1.d0/rho(:)
         nb0 = 2
      endif

 90   continue
      if(nb0.eq.1) then
         eta0 = bpr*rho(1)
         if(eta0.lt.eta_c) then ! racine gaz
            vol_mol_vap = vm(1)
            vol_mol_liq = 1.d31
         else ! racine liq
            vol_mol_liq = vm(1)
            vol_mol_vap = 1.d31
         endif
      else ! nb0 = 2
         if(vm(1).lt.bpr .and. vm(2).lt.bpr .and. vm(3).gt.bpr) then
            vm(1) = vm(3)
            vm(2:3) = 0.d0
            nb0 = 1
            goto 90
         elseif(vm(1).gt.bpr .and. vm(2).gt.bpr .and. vm(3).gt.bpr) then
            vol_mol_liq = vm(1)
            vol_mol_vap = vm(3)
         else
            write(*,*) 'RESOL_EOS corps pur -- Cas non envisage - STOP'
            err3 = 21
            goto 80
         endif
      endif

!c     Verification de la qualite de la resolution
!c     A-t-on : P specifiee = P(T,v) ?

      res = 0.d0
      if(vol_mol_liq.lt.1.d25) then
         res = res + dabs(pEoS(tk,vol_mol_liq) - pres0)
      endif
      if(vol_mol_vap.lt.1.d25) then
         res = res + dabs(pEoS(tk,vol_mol_vap) - pres0)
      endif

      if((vol_mol_liq.lt.1.d25 .and. res.lt.1.d-6) &
     &    .or. (vol_mol_vap.lt.1.d25 .and. res.lt.1.d-7)) then ! cardan OK
         goto 80
      else ! pb cardan, on lance Newton
         continue
      endif

!c     ***********************************
!c     RESOLUTION PAR LA METHODE DE NEWTON
!c     ***********************************

!c     Volume molaire de la phase vapeur
      eta0 = 1.d-5
      state = 'vap'

 32   continue
      vol_mol_vap = bpr / eta0

      pression = pEoS(tk,vol_mol_vap)

      do while(pression.gt.pres0)
         eta0 = eta0 / 2.d0
         goto 32
      end do

!c     calcul du nombre de racines vapeur (0 ou 1)
      call newton_resol_eos(tk,pres0,eta0,vol_mol_vap,state,err3)
      if(err3.ne.0) goto 80

      if(vol_mol_vap.gt.1.d25) then
         nb_vap = 0
      else
         nb_vap = 1
      endif

!c     Recherche du volume molaire de la phase liquide
      eta0 = 1d0 -1d-5
      state = 'liq'

 33   continue
      vol_mol_liq = bpr / eta0

      pression = pEoS(tk,vol_mol_liq)

      do while(pression.lt.pres0)
         eta0 = (eta0 + 1.d0) / 2.d0
         goto 33
      end do

      call newton_resol_eos(tk,pres0,eta0,vol_mol_liq,state,err3)
      if(err3.ne.0) goto 80

!c     calcul du nombre de racines liquides (0 ou 1)
      if(vol_mol_liq.gt.1.d25) then
         nb_liq = 0
      else
         nb_liq = 1
      endif

!c    calcul du nombre de racines de l'equation d'etat (1 ou 2)
      nb0 = nb_liq + nb_vap

 80   continue
      end subroutine resol_eos
!c
!c ----------------------------------------------------------------------------
!c
      subroutine newton_resol_eos(tk,pres0,eta0,vol_mol,state,err4)

!c sous-programme resolvant a temperature tk et pression pres0 fixees,
!c l'EoS d'un corps pur. Il calcule une seule racine de l'equation.
!c La nature de cette racine depend de la valeur d'initialisation
!c de eta0 (=b/v), la densite reduite par rapport au covolume.


!c ------- declarations des variables
      use comflash
      implicit none
      integer compteur,nb_iter,nb_iter_max,err4
      real*8 apr,bpr
      real*8 pres0,tk,eta0,vol_mol,residu,eta_c,eta1
      real*8 pression,pression_prime,pEoS,dpEoS,cur
      character*3 state

      common/ab_pr/apr,bpr
      err4 = 0
!c     Compacite critique :
      eta_c = cur((1.d0-r1)*(1.d0-r2)**2)+cur((1.d0-r2)*(1.d0-r1)**2) &
     &        +1.d0
      eta_c = 1.d0 / eta_c

!c     initialisation du compteur de passages dans une boucle
      compteur = 0
      residu = 1.d-10

!c     nombre maximal d'iterations autorisees
      nb_iter_max = 300

!c     initialisation du compteur du nb d'iterations
      nb_iter = 0

!c     initialisation de la variable eta1
      eta1 = eta0

 33   continue
      vol_mol = bpr / eta1

      if(nb_iter.eq.nb_iter_max) then
         if(residu.lt.1.d-4) then
            residu = 1.d-3
            nb_iter = 0
         else
            err4 = 2
            goto 35
         endif
      endif

      nb_iter = nb_iter + 1

      pression = pEoS(tk,vol_mol)

      if(state.eq.'liq') then
!c
!c        CAS DE LA RACINE LIQUIDE
!c        Lorsque eta1 < eta_c et que l'on cherche la racine liquide,
!c           cela veut dire qu'elle n'existe sans doute pas.
!c           On confirme ce resultat en ramenant eta1 a eta_c et on
!c           poursuit le calcul. Si la methode de Newton calcule a nouveau
!c           eta1 tel que eta1 < eta_c alors on confirme que la racine liquide
!c           n'existe pas
         if(eta1.lt.eta_c) then
               vol_mol = 1.d31
               goto 35
         endif
      else
!c
!c        CAS DE LA RACINE VAPEUR
!c        Lorsque eta1 > eta_c et que l'on cherche la racine vapeur,
!c           cela veut dire qu'elle n'existe sans doute pas.
!c           On confirme ce resultat en ramenant eta1 a eta_c et on
!c           poursuit le calcul. Si la methode de Newton calcule a nouveau
!c           eta1 tel que eta1 > eta_c alors on confirme que la racine vapeur
!c           n'existe pas
         if(eta1.gt.eta_c) then
            if(compteur.eq.0) then
               compteur = compteur + 1
               eta1 = eta_c
               goto 33
            else
               vol_mol = 1.d31
               goto 35
            endif
         endif
      endif

!c     (dP / deta) a T fixee
      pression_prime = dpEoS(tk,vol_mol)

      eta1 = eta1 + (pres0 - pression) / pression_prime
      vol_mol = bpr / eta1

      if(pression.lt.1.d-15) then
         err4 = 3
         goto 35
      endif

      if(dabs(pres0 - pression).gt.residu) then
        goto 33
      else
        goto 35
      endif

 35   continue
      end subroutine newton_resol_eos
!c
!c ----------------------------------------------------------------------------
!c
      function pEoS(tk,vol_mol)

!c calcul de la pression (en Pa) du fluide a T (K) et v (m3/mol) fixes
!c EoS de PR78
      use comflash
      implicit none
      real*8 apr,bpr
      real*8 pEoS,tk,vol_mol

      common/ab_pr/apr,bpr

      pEoS = rgas * tk / (vol_mol - bpr) - apr / &
     &       ((vol_mol-bpr*r1)*(vol_mol-bpr*r2))

      end function pEoS
!c
!c ----------------------------------------------------------------------------
!c
      function dpEoS(tk,vol_mol)

!c calcul de dP/d eta (en Pa) du fluide a T (K) fixe
!c EoS de PR78
      use comflash
      implicit none
      real*8 apr,bpr
      real*8 dpEoS,tk,vol_mol,dpdvpur
      common/ab_pr/apr,bpr

      dpdvpur = -rgas*tk/(vol_mol-bpr)**2 &
     &         + apr*(2.d0*vol_mol - (r1+r2)*bpr)/ &
     &         ((vol_mol-bpr*r1)*(vol_mol-bpr*r2))**2

      dpEoS = -vol_mol**2/bpr * dpdvpur

      end function dpEoS

!c
!c ----------------------------------------------------------------------------
!c
      function lnphi(tk,vol_mol)

!c calcul du logarithme du coefficient de fugacite du corps pur
!c a la pression pres0, la temperature tk
!c (volume molaire correspondant : vol_mol) en utilisant
!c l'equation de peng et robinson de parametres apr et bpr


!c ------- declarations des variables
      use comflash
      implicit none
      real*8 apr,bpr
      real*8 pres0,tk,vol_mol,z0,lnphi,pEoS

      common/ab_pr/apr,bpr

      pres0 = pEoS(tk,vol_mol)
      z0 = pres0 * vol_mol / (rgas * tk)

      lnphi = z0 - 1.d0 - dlog(z0 - bpr*pres0/rgas/tk) &
     &        + apr / (rgas*tk*bpr*(r1-r2)) &
     &          *dlog((vol_mol-bpr*r1)/(vol_mol-bpr*r2))

      end function lnphi
!c
!c ----------------------------------------------------------------------------
!c
      function h_res(numero,tk,vol_mol)

!c calcul de l'enthalpie residuelle du fluide (J / mol) a la temperature tk
!c (volume molaire correspondant : vol_mol) en utilisant
!c l'equation de peng et robinson de parametres apr et bpr
      use comflash
      implicit none
      integer numero
      real*8 apr,bpr
      real*8 tk,vol_mol,daprdt,z0,pres0,h_res,pEoS
      real*8 a_function

      common/ab_pr/apr,bpr

      pres0 = pEoS(tk,vol_mol)
      z0 = pres0 * vol_mol / (rgas * tk)
      daprdt = a_function(numero,tk,1) ! da / dT

      h_res = rgas * tk * (z0 - 1.d0) &
     &        + 1.d0/((r1-r2) * bpr) * (apr - tk*daprdt) &
     &          *dlog((vol_mol-bpr*r1)/(vol_mol-bpr*r2))

      end function h_res
!c
!c ----------------------------------------------------------------------------
!c
      function cur(x)

!c Definition de la fonction : x --> cur(x) = x^(1 / 3)
!c Domaine de definition     : ]-infini ; + infini[

      implicit none
      real*8 x,cur

      if(x.le.0.d0) then
         cur = -dabs(x)**(1.d0 / 3.d0)
      else
         cur = x**(1.d0 / 3.d0)
      endif

      end function cur
!c
!c ----------------------------------------------------------------------------
!c
       function a_function(numero,tk,ideriv)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c
!c Calcul du parametre a(T) dans le terme attractif de l'EoS cubique
!c des corps purs
!c numero = entier utilise pour designer le corps pur considere
!c Unite de a(T) : Pa.m6/mol2 (systeme international)
!c ideriv : indice de dérivation :
!c           - si ideriv = 0 : a_function renvoie a/(RT)
!c           - si ideriv = 1 : a_function renvoie d[a/(RT)] / dT
!c           - si ideriv = 2 : a_function renvoie a/(RT)
!c
!c ATTENTION : - la fonction a(T) est egalement definie dans la sous-routine
!c             parametres_EoS (pour le calcul de Psat et Teb des corps purs)
!c             - la fonction a(T) est egalement definie dans la sous-routine
!c             temset (pour le calcul des hM et cPM)
!c
       use comflash
       implicit none
       integer numero,ideriv
       real*8 pcrit
       real*8 a_function,tk
       real*8 alpha,rac_alpha
       real*8 delt,gamm,theta
       real*8 polTH

       if(ideriv==0) then ! calcul de la fonction a(T)
            if(FoncAlpha.EQ.'SOAVE') then
!c         -------------------------
!c         Fonction alpha de SOAVE :
!c         -------------------------
                rac_alpha = 1.d0+msoave(numero)*(1.d0- &
     &                                            DSQRT(tk/tc(numero)))
                alpha = rac_alpha**2
!c         -------------------------
!c         -------------------------
            elseif(FoncAlpha.EQ.'TWU91'.OR.FoncAlpha.EQ.'TWU88') then
!c         -------------------------
!c         Fonction alpha de TWU:
!c         -------------------------
                delt = parN(numero) * (parM(numero) - 1.d0)
                gamm = parN(numero) * parM(numero)
                alpha = (tk/tc(numero))**delt * DEXP(parL(numero)*(1.d0 &
     &                                         - (tk/tc(numero))**gamm))
!c         -------------------------
!c         -------------------------
            endif
            a_function = ac(numero)*alpha ! Pa.m6/mol2
       elseif(ideriv == 1) then ! calcul de d a(T) / dT
            if(FoncAlpha.EQ.'SOAVE') then
!c         -------------------------
!c         Fonction alpha de SOAVE :
!c         -------------------------
                rac_alpha = 1.d0+msoave(numero)*(1.d0- &
     &                                            DSQRT(tk/tc(numero)))
                pcrit = pc(numero) * 1.d5 ! bar > Pa
!c         da / dT dans les unites du systeme international
                a_function = -msoave(numero)*ac(numero)*rac_alpha &
     &                      /(dsqrt(tc(numero)*tk))
!c         -------------------------
!c         -------------------------
            elseif(FoncAlpha.EQ.'TWU91'.OR.FoncAlpha.EQ.'TWU88') then
!c         -------------------------
!c         Fonction alpha de TWU:
!c         -------------------------
                delt = parN(numero) * (parM(numero) - 1.d0)

                gamm = parN(numero) * parM(numero)
                alpha = (tk/tc(numero))**delt * DEXP(parL(numero)*(1.d0 &
     &                                         - (tk/tc(numero))**gamm))
                a_function = ac(numero) * alpha * tc(numero) / tk * &
     &          (delt - parL(numero) * gamm * (tk/tc(numero))**gamm)
                a_function = a_function / tc(numero)
            endif
       elseif(ideriv == 2) then ! calcul de d2 a(T) / dT2
            if(FoncAlpha.EQ.'SOAVE') then
!c         -------------------------
!c         Fonction alpha de SOAVE :
!c         -------------------------
                 a_function = msoave(numero)*(1.d0+msoave(numero))* &
     &                  ac(numero)/(2.d0* tk**1.5d0 * dsqrt(tc(numero)))
            elseif(FoncAlpha.EQ.'TWU91'.OR.FoncAlpha.EQ.'TWU88') then
!c         -------------------------
!c         Fonction alpha de TWU:
!c         -------------------------
                delt = parN(numero) * (parM(numero) - 1.d0)
                gamm = parN(numero) * parM(numero)
                alpha = (tk/tc(numero))**delt * DEXP(parL(numero)*(1.d0 &
     &                                         - (tk/tc(numero))**gamm))
                theta = parL(numero) * gamm * (tk/tc(numero))**gamm
                polTH = theta**2.d0 - (2.d0*delt + gamm - 1.d0)*theta + &
     &                                              delt * (delt - 1.d0)
                a_function = ac(numero) * alpha * polTH * (tc(numero) / &
     &                                                     tk)**2
                a_function = a_function / tc(numero)**2
             endif
       else
          write(lun1,*) 'ideriv non defini (FUNCTION a_function)'
!          write(lun2,*) 'ideriv non defini (FUNCTION a_function)'
          stop
       endif

       end function a_function
!c
!c ----------------------------------------------------------------------------
!c
      subroutine calcul_teb(numero,p_in,teb,err0,vol_mol_liq, &
     &                                                     vol_mol_vap)

!c Calcul de la temperature d'ebullition d'un corps pur
!c par un EoS cubique.
!c
!c INPUT : numero, entier = 1 ou 2 pour spécifier le constituant
!c         p_in (pression specifiee) en bar
!c OUTPUT :  teb (K), err0 (variable d'erreur, = 0 si tout va bien)
!c
      use comflash
      implicit none
      real*8 apr,bpr,pcrit,tcrit
      real*8 teb,teb0,p_in,vol_mol_liq,vol_mol_vap
      real*8 lnphi_liq,lnphi_vap,lnphi,h_res,h_res_liq,h_res_vap
      real*8 pres0,a_function
      real*8 res1,res2,res3,vvap_old,vliq_old,eps

      integer numero,nb000,nb_iter,nb_iter_max,err0,err1
      common/ab_pr/apr,bpr

      eps = 1.d-9
      teb = 0.d0
      err0 = 0
      nb_iter_max = 300
      nb_iter = 0  ! compteur d'iterations
      pres0 = p_in*1.0d+5  !en Pa
      bpr = b(numero) * 1.d-6 ! cm3/mol > m3/mol

      pcrit = pc(numero) * 1.d5 ! bar > Pa
      tcrit = tc(numero) ! K
      if(dabs(pcrit - pres0) / pcrit.lt.1.d-5) then
         teb = tcrit
         goto 31
      endif

      if(pres0.gt.pcrit) then
         err0 = 1  ! P > Pc et Teb non calculable
         goto 31
      endif
!c     procedure d'initialisation teb :
!c     ------------------------------
      call initialisation_teb(numero,teb,pres0,err0)
      if(err0.ne.0) goto 31

!c     Methode de Newton pour calculer teb :
!c     -----------------------------------
 30   continue
!c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,teb,0)

      vliq_old = vol_mol_liq
      vvap_old = vol_mol_vap
      call resol_eos(teb,pres0,vol_mol_liq,vol_mol_vap,nb000,err1)

      if(err1.ne.0) err0 = 10
      if(err0.ne.0) goto 31

      lnphi_liq = lnphi(teb,vol_mol_liq)
      lnphi_vap = lnphi(teb,vol_mol_vap)

      h_res_liq = h_res(numero,teb,vol_mol_liq)
      h_res_vap = h_res(numero,teb,vol_mol_vap)

      teb0 = teb

      if(nb_iter.eq.nb_iter_max) then
         eps = 1.d-6

         if(res1.lt.eps .and. res2.lt.eps .and. res3.lt.eps) goto 31
         err0 = 2
         goto 31
      endif
      nb_iter = nb_iter + 1

      teb = teb0 + rgas * teb**2 * (lnphi_vap - lnphi_liq) &
     &             / (h_res_vap - h_res_liq)

      res1 = dabs(lnphi_vap - lnphi_liq)
      res2 = dabs(teb - teb0)/dabs(teb0)
      res3 = 0d0
      if(nb_iter>1) then
         res3 = dabs(vol_mol_liq - vliq_old)/dabs(vliq_old) &
     &          + dabs(vol_mol_vap - vvap_old)/dabs(vvap_old)
      endif

      if(res1.lt.eps .and. res2.lt.eps .and. res3.lt.eps &
     &   .and. nb_iter>1) then
         goto 31
      else
         goto 30
      endif

 31   continue

      end subroutine calcul_teb
!c
!c ----------------------------------------------------------------------------
!c
      subroutine initialisation_teb(numero,teb,pres0,err0)

!c Initialisation de la temperature d'ebullition  :
!c   1) initialisation a partir d'une formule l'approximant connaissant
!c      tc, pc et omega.
!c   2) Affinage de cette valeur : on calcule teb de sorte que la
!c      resolution de l'equation d'etat conduise a une racine liquide et
!c      une racine vapeur
      use comflash
      implicit none
      real*8 tcrit,pcrit,omega0
      real*8 apr,bpr,a_function
      real*8 pres0,teb,vol_mol_liq,vol_mol_vap,tinf,tsup
      integer nb0,nb_iter,nb_iter_max,err0,numero

      common/ab_pr/apr,bpr

      pcrit = pc(numero)*1.d5 ! bar>Pa
      tcrit = tc(numero)
      omega0 = fac(numero)

      nb_iter = 0
      nb_iter_max = 1000
      err0 = 0

      teb = tcrit*(omega0 + 1.d0) / (omega0 + 1.d0 - 3.d0/7.d0 &
     &                             * dlog10(pres0 / pcrit))

!c     La resolution de l'EoS a T=teb et P=pres0 amene-t-elle 3 racines ?
 50   continue
!c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,teb,0)
      call resol_eos(teb,pres0,vol_mol_liq,vol_mol_vap,nb0,err0)
      if(err0.ne.0) goto 52

      if(nb_iter.eq.nb_iter_max) then
         err0 = 2
         goto 52
      endif
      nb_iter = nb_iter + 1

!c     si la racine vapeur a ete trouvee mais pas la racine liquide :
!c     ------------------------------------------------------------
      if(vol_mol_liq.gt.1.d25) then
!c        diminuer la valeur de teb jusqu'a ce que 3 racines apparaissent
         teb = 2d0 * teb - tcrit
         nb_iter = nb_iter + 1
         goto 50
      endif

!c     si la racine liquide a ete trouvee mais pas la racine vapeur :
!c     ------------------------------------------------------------
      if(vol_mol_vap.gt.1.d25) then
!c        Recherche de la valeur de teb conduisant a 3 racines par
!c        dichotomie

!c        Valeur de Teb initiale :
         tinf = teb
         tsup = tcrit
         teb = 0.5d0 * (tinf + tsup)
         nb_iter = nb_iter + 1

 51      continue

!c        Resolution de l'equation d'etat a pres0 et teb donnees
!c        a_function donne a(T) en Pa.m6/mol2
         apr = a_function(numero,teb,0)
         call resol_eos(teb,pres0,vol_mol_liq,vol_mol_vap,nb0,err0)
         if(err0.ne.0) goto 52

         if(nb_iter.eq.nb_iter_max) then
            err0 = 2
            goto 52
         endif
         nb_iter = nb_iter + 1

!c        si la racine vapeur a ete trouvee mais pas la racine liquide :
         if(vol_mol_liq.gt.1.d25) then
            tsup = teb
            teb = 0.5d0 * (tinf + tsup)
            nb_iter = nb_iter + 1
            goto 51
         endif

!c        si la racine liquide a ete trouvee mais pas la racine vapeur :
         if(vol_mol_vap.gt.1.d25) then
            tinf = teb
            teb = 0.5d0 * (tinf + tsup)
            nb_iter = nb_iter + 1
            goto 51
         endif
      endif

 52   continue

      end subroutine initialisation_teb
!c
!c ----------------------------------------------------------------------------
!c
      subroutine calcul_racine_stable_pur(numero,tk,pbar00,vm,err1)
!c
!c INPUT : numero (du corps pur => 1 ou 2)
!c         tk (temp. specifiee) en K, pbar00 en bar
!c OUTPUT :  vm (m3/mol), err0 (variable d'erreur, = 0 si tout va bien)
!c
      use comflash
      implicit none
      real*8 apr,bpr
      real*8 pres,tk,vm1,vm2,vm
      real*8 lnphi1,lnphi2,lnphi
      real*8 pbar00,a_function
      integer err1,nb0,numero
      common/ab_pr/apr,bpr

      pres = pbar00*1.0d+5 !en Pa
      bpr = b(numero)*1.d-6 ! cm3/mol > m3/mol

!c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,tk,0)
      err1 = 0
      call resol_eos(tk,pres,vm1,vm2,nb0,err1)
      if(err1.ne.0) return

      if(nb0.eq.1) then
         if(vm1.lt.1.d30) then
            vm = vm1
         else
            vm  = vm2
         endif
      else
         lnphi1 = lnphi(tk,vm1)
         lnphi2 = lnphi(tk,vm2)
         if(lnphi1.lt.lnphi2) then
            vm = vm1
         else
            vm = vm2
         endif
      endif

      end subroutine calcul_racine_stable_pur
!c
!c ----------------------------------------------------------------------------
!c
      subroutine ispure(x,icomp)
!c  Determine if the user has requested the properties of a pure fluid.
!c  This happens if 1) nc=1, 2)  one of the compositions in the x array is one.
!c
        use comflash
        implicit none
        integer :: icomp,i
        real(8) :: x(nco)
!c-----------------------------------------------------------------------------
        icomp = 0
        if(nco .eq.1)then
            icomp = 1
        else
            do i = 1,nco
               if (DABS(x(i)-1.d0).lt.1.d-12) then
                 icomp=i
                 RETURN
               endif
            end do
        end if

      end subroutine
!c
!c-------------------------------------------------------------------------------
!c
      function calc_c(icomp,typec)
        use comflash
        implicit none
        real(8) :: tk,calc_c,coeff(nb_coeff),dippr_y
        real(8) :: vliq_exp,vliq_eos,psat,vliq,vvap
        integer :: icomp,eqID,typec,err0

        coeff(:) = tdep_prop(icomp,ldn_idx,:)
        eqID = eqIDComp(icomp,ldn_idx)

        if(typec.eq.1)then !from correlation at Tr = 0.8
            tk = 0.8d0 * tc(icomp)
            vliq_exp = dippr_y(eqId,tc(icomp),coeff,tk)
            vliq_exp = 1.0d0 / vliq_exp * 1d03 ! cm3 mol-1

            call calcul_psat(icomp,tk,psat,err0,vliq,vvap)
            vliq_eos = vliq * 1.0d06 ! cm3 mol-1
            calc_c = vliq_eos - vliq_exp
        elseif(typec.eq.2)then
             calc_c = rgas* tc(icomp) / pc(icomp) * &
     &                (5.0d-06 + 0.0047d0 * fac(icomp)) * 10.0d0 !cm3 mol-1
        end if
      end function
!c
!c-------------------------------------------------------------------------------
!c
!c*****************************************************************************
!   END: PURE COMPONENTS SUBROUTINES/FUNCTIONS
!c*****************************************************************************
!c
!c
!c*****************************************************************************
!   MIXTURES SUBROUTINES/FUNCTIONS
!c*****************************************************************************
!c     thermodynamics block
!c
      subroutine FuPrLorraine_TP(icon,ntyp,mtyp,ierr,tKelvin,pMPa,zn,rho,zFactor, &
     &                       fg,ft,fp,fx,uRes_RT,aRes_RT,CvRes_R,CpRes_R,dpdRho_RT)	!JRE: fg ~ fugc
	  !call FuPrLorraine_TP(icon,ntyp,mtyp,iErrPRL,T,P,X,rho,Z, &	 ! JRE 20210510
	  !	&                FUGC,ft,fp,fx,uRes_RT,aRes_RT,CvRes_R,CpRes_R,dpdRho_RT)
	  use GlobConst, only:dumpUnit,form610
      use comflash
      implicit none
      integer :: icon,ntyp,mtyp,ierr
      real(8) :: tKelvin,pMPa,zn(ncomax),rho,zFactor,uRes_RT,denom,etaPass !,xFrac(ncomax)
      real(8) :: aRes_RT,hRes_RT,sRes_R,CvRes_R,CpRes_R,dpdRho_RT,sResCheck
      real(8) :: fg(ncomax),ft(ncomax),fp(ncomax),fx(ncomax,ncomax),aux(15)
	  LOGICAL FORT	! JRE20210615
!c-----------------------------------------------------------------------------
	  FORT=.TRUE.	! JRE20210615
      ntemp = 0
      nder = 1
      ierr = 0
      if ( icon.gt.2 ) nder = 2
      if ( icon.eq.2 .or. icon.gt.3 .or. icon.eq.0 ) ntemp = 2
      if ( icon.gt.4 ) ntemp = 2

      call cubgen(ntyp,mtyp,ierr,tKelvin,pMPa,zn,zFactor,fg,ft,fp,fx,aux) !JRE: pause here to enter sample calcs.
	  if(ierr .ne. 0)then
		  if(FORT)write (dumpUnit,*) ' FuPrLorraine_TP: CubGen(), iErr=',iErr
		  return	! JRE20210615. Otherwise, zFactor=0 in denom can occur.
	  endif

      !Outputs
	  denom=(zFactor * rgas * tKelvin)
	  if(denom < 1.d-11)then  ! JRE20210615
		if(FORT)write (dumpUnit,form610)' FuPrLorraine_TP: denom <= 0. Z,T,R',zFactor,tKelvin,rgas
		if(FORT)write (dumpUnit,form610) ' Check denom=',denom
	  endif
      rho = pMPa/(zFactor * rgas * tKelvin) !mol/cm3
	  etaPass=eta ! eta is USEd from comflash, needed for PsatEar

      uRes_RT = aux(3)/(rgas*tKelvin)										
      aRes_RT = aux(4)/(rgas*tKelvin)
      uRes_RT = aux(11)										
      aRes_RT = aux(12)
      sRes_R  = aux(14)
      sResCheck=uRes_RT - aRes_RT
      hRes_RT = aux(1) 

      CpRes_R = aux(6)
      CvRes_R = aux(7)
      dpdRho_RT = - aux(8) / rho**2 * 1.0d-12 !aux(10) = (dPdv)/RT (uSI)


      end subroutine !FuPrLorraine_TP
!c
!c----------------------------------------------------------
!c
      subroutine FuPrLorraine_TV(icon,ierr,tKelvin,vTotCm3,zn,pMPa,zFactor, &
     &                       fg,ft,fp,fx,uRes_RT,aRes_RT,CvRes_R,CpRes_R,dpdRho_RT)

      use comflash
      implicit none
      integer :: icon,ierr
      real(8) :: tKelvin,pMPa,zn(ncomax),vTotCm3,zFactor,uRes_RT
      real(8) :: aRes_RT,CvRes_R,CpRes_R,dpdRho_RT,rho
      real(8) :: fg(ncomax),ft(ncomax),fp(ncomax),fx(ncomax,ncomax),aux(15)
!c-----------------------------------------------------------------------------
      ntemp = 0
      nder = 1
      ierr = 0
	  zn=zn/sum(zn)
      if ( icon.gt.2 ) nder = 2
      if ( icon.eq.2 .or. icon.gt.3 .or. icon.eq.0 ) ntemp = 2
      if ( icon.gt.4 ) ntemp = 2

      call volgen(tKelvin,vTotCm3,zn,pMPa,fg,ft,fp,fx,aux,ierr)

      if(ierr.ne.0) return

      !Outputs
      rho = sum(zn(1:nco)) / vTotCm3 !rho = n / V (mol/cm3)
      zFactor = pMPa / rgas / tKelvin / rho

      uRes_RT = aux(11)
      aRes_RT = aux(12)

      CpRes_R = aux(6)/rgas
      CvRes_R = aux(7)/rgas

      dpdRho_RT = - aux(8) / rho**2 * 1.0d-12 !aux(10) = (dPdv)/RT (uSI)

      end subroutine
!c
!c-----------------------------------------------------------------------------
!c
      subroutine cubgen(ntyp,mtyp,ierr,t,p,xx,zz,fg,ft,fp,fx,aux)
!c
!c     parameters:
!c
!c     mt:     (i):      phase type desired
!c                       1 = liq, -1 = vap, 0 = min g phase
!c     t:      (i):      temperature (k)
!c     p:      (i):      pressure (MPa)
!c     xx:     (i):      mixture total composition
!c
!c     ic:     (o):      indicator for phase type; ic returns
!c                       1:   if called with mt = 1/-1 and phase is found
!c                       -1:  if called with mt = 1/-1 and phase is not found
!c                       2:   if called with mt = 0 and min g phase is liquid
!c                       -2:  if called with mt = 0 and min g phase is vapour
!c     z:      (o):      calculated compressibility factor
!c     fg:     (o):      log fugacity coefficient vector
!c     ft:     (o):      t-derivative of fg
!c     fp:     (o):      p-derivative of fg
!c     fx:     (o):      scaled composition derivative of fg
!c     aux:    (o):      various residual properties (cf. volgen for defs)
!c
      use comflash
      implicit none
      integer :: ntyp,mtyp
!c-----dummy variables
      real(8) :: t,p,zz
      real(8) :: xx(ncomax),fg(ncomax),ft(ncomax),fp(ncomax),fx(ncomax,ncomax)
      real(8) :: aux(15)
!c-----local variables
      integer :: i,j,k,n,ider,istab,icx,i3r,ierr
!      integer :: i3r
      real(8) :: v0,lnphi_c(ncomax),rho,zDeTrans
      real(8) :: x(ncomax)
      real(8) :: Cpec,Cvec,eta0,g0,g0r,hec,sec,sumr,sumx
!c-----------------------------------------------------------------------------
      if(ntyp.eq.1) j = 1
      if(ntyp.eq.-1) j = 2
      istab = 0

      if(ntyp.eq.0) then
          istab = 1
          j = 3
      end if

      if(nder.eq.1) ider = 0
      if(nder.eq.2 .or. ntemp.eq.2) ider = 1


      n = nco
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

      call regles_melange(ider,ntemp,x)

      eta0 = 8888d0
      call phasLorraine(j,ider,istab,p0,x,eta0,v0,g0,i3r,ierr)	!v0 is molar volume 

      if(ierr.ne.0) return !JRE20210615. Otherwise, Z=0 in denom.
        !ierr = 21
        !'RESOL_EOS -- Problem - eta1 and eta2 < 0'

        !ierr = 22
        !'RESOL_EOS -- Problem - eta1 and eta2 > 1'

        !ierr = 23
        !'RESOL_EOS -- Problem - vapor root (eta1) < 0'

        !ierr = 24
        !'RESOL_EOS -- Problem - liquid root (eta2) > 1'

      icx = 1
      if (v0 .gt. 3.d0*bm) icx = -1
!c
!c     phase type return indicator
!c
      if (ntyp.ne. 0) then
          mtyp = icx*ntyp
      else
          mtyp = 2*icx
      endif

      if(i3r.ge.2) mtyp = 3    !eta is equal to the real part of the imaginary root

!      if(mtyp.eq.-1 .and. istab.eq.0) then
!         eta0 = 8888d0
!         call rcub(i3r,eta0)
!         if(i3r.eq.2) mtyp = 3    !eta is equal to the real part of the imaginary root
!         eta = eta0
!		 if( eta < 0 .or. eta > 1)then !JRE 20210513
!			pause 'cubgen: eta violates 0 < eta < 1'
!			continue
!		 endif
!         call eqet
!         call fug
!         v0 = v
!      end if

      zz = z							! "z"=zFactor is computed in phasLorraine and is a member of comflash
      g0r = 0.0d0 !gRes/RT
      do i=1,nco
         g0r = g0r + x(i)*lphi(i)
      end do

      lnphi_c(1:nco) = p0 * c(1:n) / rt

      aux(8) = dpdv / t / rgas * 1.0d11 !dpdv / RT (uSI)
      aux(9) = dpdt * 1.0d-01           !dpdt      bar --> MPa
      aux(10) = dpdt / dpdv * 1.0d-06   !dvdt  (uSI)
!c
!c    are temperature derivatives required
!c
      if ( ntemp.ne.0 ) then
         call hec_calc(x,t,v0,hec)
        zDeTrans = z - cm * p0 / rt
         rho=p0/(zDeTrans*rgas*t)
!c------excess enthalpy/r stored in aux(3)
         aux(1)  = (hec - cm*p) / rgas
         aux(13) = (hec - cm*p) / rgas/t    ! hRes/RT
         sec = hec / t / rgas - g0r !sDep/R = (hDep - gDep)/RT
!c------Sdep/Rgas stored in aux(4)
         aux(2) = sec
		 aux(14)= sec +LOG(zz)       !Sres/R
         aux(15)= sec +LOG(zDeTrans)
!c------excess internal energy/RT stored in aux(11)
         aux(11) = hec / t / rgas - (zz-1) !uRes/RT = (hRes - PV)/RT = hRes/RT - Z
!c------excess Helmholtz energy/RT stored in aux(12)
         aux(12) = g0r - (zz-1) + LOG(zz) !aDep/RT = (gDep - PV)/RT = gDep/RT - Z
         aux(4) = aux(12) * rgas * t		! aux( 4)=Adep
!c
!c     residual cp required
!c
         if ( ntemp.ge.2 ) then
            call Cp_Cvec_calc(x,t,v0,Cpec,Cvec)
!c------excess heat capacity/r in aux(5)
            aux(6) = Cpec / rgas
            aux(7) = Cvec / rgas
         endif
      endif

      if ( nder.ne.0 ) then
!c
!c     log fugacity coefficients
!c
         fg(1:n) = lphi(1:n) - lnphi_c(1:n)
         if ( nder.ge.2 .or. ntemp.ne.0 ) then
!c
!c     composition derivatives required
!c
            if ( nder.ge.2 ) then
               do i = 1,n
                 do k = 1,n
                    fx(i,k) = dlphidn(i,k)
                    fx(k,i) = fx(i,k)
                 enddo
               enddo
            endif
!c
!c     temperature and pressure derivatives
!c
            if ( ntemp.ne.0 ) then
               fp(1:n) = (dlphidp(1:n) - lnphi_c(1:n) / p0)*10.0d0 !bar-1 -----> MPa-1
               ft(1:n) = dlphidt(1:n) + lnphi_c(1:n) / t
            endif
         endif
      endif
      zz = z - cm * p0 / rt
      aux(5) = (g0r - cm*p / rgas / t) !gRes/RT (J/mol)

      end !subroutine cubgen
!c
!c-----------------------------------------------------------------------------
!c		
      subroutine volgen(t,vv,x,p,fg,ft,fv,fx,aux,ierr)
	  use GlobConst, only: dumpUnit
      use comflash
      implicit none
      integer :: ierr,i,j,k,n
      integer :: ider
      real(8) :: t,vv,p,som,sumn,vv2,aux(15)
      real(8) :: x(ncomax),fg(ncomax),ft(ncomax),fv(ncomax)
      real(8) :: fx(ncomax,ncomax),x2(ncomax),lnphi_c(ncomax)
      real(8) :: dfgdv(ncomax),dfgdn(ncomax,ncomax)
      real(8) :: g0r,hec,sec,Cpec,Cvec,zz
!c
!c-------parameters of volgen (crit. point calc)
!c
!c       t     temperature  (k)                          (input)
!c       v     total volume (cm3)                        (input)
!c       x---  mixture mole numbers                     (input)
!c       p     pressure    (mpa)                        (output)
!c       fg    log fugacity coefficient vector          (output)
!c       ft    t-derivative of fg (const. vol)          (output)
!c       fv    vol-derivative of fg (const temp)        (output)
!c       fx    comp-derivative of fg (const t. and v)   (output)
!c       aux   various residual properties              (output)
!              1=Hdep [J/mol]
!              2=Sdep [J/molK]
!              3=Udep [J/mol]
!              4=Adep [J/mol]
!              5=excess heat capacity/r (?)
!              6=(Cp-Cpig)/Rgas
!              7=(Cv-Cvig)/Rgas
!              8=dpdv_RT (uSI)
!              9=dpdt bar --> MPa
!              10=dvdt  (uSI)
!              11=Ures/RT
!              12=Ares/RT
!c---------------------------------------------------
!c---  modified and corrected april 2021
!c---
!c---------------------------------------------------
!c
!c     sumn is total moles; volgen calculates in terms of overall moles
!c     and overall volume
!c
      if(nder.eq.1) ider = 0
      if(nder.eq.2 .or. ntemp.eq.2) ider = 2

      n = nco
      sumn = 0.d0
      do i = 1,n
         sumn = sumn + x(i)
      enddo

      x2(1:n) = x(1:n) / sumn

      t0 = t
      call regles_melange(ider,ntemp,x2)

!c     Eta calculation
      vv2 = vv/sumn + cm !"Translation". vv is the input value of vTotCm3
      eta = bm/vv2
      if(eta < 0 .or. eta > 1)then
        ierr = 21
        write (dumpUnit,'(a,8E12.4)') ' volgen: eta,bm,vv2',eta,bm,vv2
        continue
        return
      endif

      call eqet !Calculation of z and P
      p = peq * 1.0d-01 ! bar ---> MPa
      call fug !Calculation of ln phi
      zz = z

      lnphi_c(1:nco) = peq * c(1:n) / rt

      g0r = 0.0d0
      do i=1,nco
         g0r = g0r + x2(i)*lphi(i)
      end do

      if ( nder.ne.0 ) then
!c       fg    logarithm of fugacity coefficient  (output)
         do i = 1,n
           fg(i) = lphi(i) - lnphi_c(i)
         end do
         if ( nder.ge.2 .or. ntemp.ne.0 ) then
            ider = 2
            call dertp(ider,x2)
!c
!c  volume derivative of fugacity
!c
            fv(1:n) = dlphidv(1:n)
            dfgdv(1:n) = - fv(1:n) * bm/eta / sumn
            fv(1:n) = (fv(1:n) - lnphi_c(1:n) * dpdv/peq) / sumn
            aux(8) = dpdv / rgas / t / sumn * 1.0d11 !dpdv_RT (uSI)
            aux(9) = dpdt * 1.0d-01 !dpdt bar --> MPa
            aux(10) = dpdt / dpdv * 1.0d-06   !dvdt  (uSI)

            if ( nder.ge.2 ) then
                do i = 1,n
                  do k = i,n
                    som = 0.0d0
                    do j = 1,n
                       if(j == k) som = som + dlphidn(i,j) * (1.0d0-x2(j))
                       if(j /= k) som = som - dlphidn(i,j) * x2(j)
                    enddo
                     dfgdn(i,k) = som / sumn
                     fx(i,k) = dfgdv(i) + dfgdn(i,k) - lnphi_c(i)*dpdn(k)/peq
                     fx(k,i) = fx(i,k)
                  enddo
                enddo
            endif
            if ( ntemp.ne.0 ) then
                call hec_calc(x,t,vv2,hec)
!c--------------excess enthalpy/rt stored in aux(1)
                aux(1) = (hec - cm*p) / rgas 		 ! JRE I changed peq(bar) to p(MPa)
                sec = (hec - g0r*rgas*t) / t / rgas
!c--------------excess entropy/r stored in aux(2) = Sdep/R.
                aux(2) = sec
!c--------------excess internal energy/RT stored in aux(3,11)
                aux(3) = hec - zz * rgas * t
				aux(11)= aux(3)/rgas/t					!JRE: I added #11, 12.
!c------ -------excess Helmholtz energy/RT stored in aux(4,12)
                aux(4) = g0r - zz * rgas * t		! aux( 4)=Adep
				aux(12)= aux(4)/rgas/t !+ DLOG(zz)		! JRE: aux(12)=Adep/RT, don't subtracd ln(Z) yet. 

                ft(1:n) = dlphidt(1:n) - lnphi_c(1:n) * (dpdt/peq + 1.0d0/t)    !d ln phi(i) / dT

!c
!c              residual cp/cv required
!c
                if ( ntemp.ge.2 ) then
                    call Cp_Cvec_calc(x,t,vv2,Cpec,Cvec)
!c------------------excess heat capacity/r in aux(5)
                    aux(6) = Cpec / rgas
                    aux(7) = Cvec / rgas
                 endif
            endif
         endif
      endif
      zz = z - cm * peq / rt	!Zactual=Pbar/RT*(Vtrans-c)
      aux(5) = (g0r - cm*p / rgas / t) !gDep/RT (J/mol) = GresTrans/RT -pMPa/RT*c
      end ! subroutine volgen()
!c
!c-----------------------------------------------------------------------------
!c
      subroutine phas1(j,ider,istab,p00,x,eta0,v0,g0)
!c
!c  calcul des proprietes d'un systeme suppose monophasique, avec si
!c  si istab=1 on retient la solution la plus stable
!c  si istab=2 on etudie la stabilite par rapport a la diffusion.
!c  si istab=10 (on veut choisir sa solution)
!c
       use comflash
       implicit none
!c
       integer nchol
       parameter ( nchol=ncomax*(ncomax+3)/2 )
!c
!c*** variables de la liste d'appel
!c
       integer j,ider,istab
       real*8 p00,eta0,v0,g0
       real*8 x(ncomax)
!c
!c*** variables internes
!c
       integer i,i3r,isym,irtn
       real*8 eta1,g,eta_in
       real*8 ttt(nchol)
!c-------------------------------------------------
!c      Si perte de precision de l'ordinateur :
       if(x(1).eq.0.d0) then
          x(1) = 1.d-15
          x(2) = 1.d0 - x(1)
       endif

       if(x(2).eq.0.d0) then
          x(2) = 1.d-15
          x(1) = 1.d0 - x(2)
       endif


       if(istab.eq.10) then
          eta_in = eta0
!c
!c         en sortie, on veut 1 valeur de eta differente
!c
       endif
       n0 = ncomax
       m0 = ncomax
!c
       call regles_melange(ider,8888,x)
       if(j.eq.3) eta = eta0
       call rop(j,p00,eta0)
       call fug

       eta0 = eta
       v0 = v
       g0 = 0.0D+00
       do i=1,nco
          if(x(i).le.0.d0) then
             write(lun1,*) "Erreur x = 0 ou x < 0 (PHAS1)"
             write(lun1,*) "Systeme : ",nom(1)(1:len_trim(nom(1))), &
     &                     " - ",nom(2)(1:len_trim(nom(2)))
             write(lun1,*) "T = ",t0,"P = ",p00
             write(lun1,*) "x1 = ",x(1),"x2 = ",x(2)
!             write(lun2,*) "Erreur x = 0 ou x < 0 (PHAS1)"
!             write(lun2,*) "Systeme : ",nom(1)(1:len_trim(nom(1))), &
!     &                     " - ",nom(2)(1:len_trim(nom(2)))
!             write(lun2,*) "T = ",t0,"P = ",p00
!             write(lun2,*) "x1 = ",x(1),"x2 = ",x(2)
          endif
          g0 = g0 + x(i)*(dlog(x(i)) + lphi(i))
       end do

       if(istab.eq.0) goto 10 ! on se contente de la solution
!c                             ! trouvee en premier
       call rcub(i3r,eta1)
       if(i3r.eq.0) goto 10   ! si l'equation n'a qu'une racine
!c
       eta = eta1             ! autre solution
       call eqet
       call fug

       g = 0.0D+00
       do i=1,nco
          g = g + x(i)*(dlog(x(i)) + lphi(i))
       end do

       if(istab.eq.10) then
          if(dabs(eta1-eta_in).gt.0.001) then !on veut eta1
             eta0 = eta
             g0 = g
             v0 = v
          else !on veut eta0
             eta = eta0          ! retour a la premiere solution
             call eqet
             call fug
          endif
          istab = 1
          goto 10
       endif

       if(g.lt.g0) then       ! si la 2eme solution est plus stable
          eta0 = eta
          g0 = g
          v0 = v
       else
          eta = eta0          ! retour a la premiere solution
          call eqet
          call fug
       endif

  10   if(ider.eq.1) call dertp(ider,x)
       if(istab.ne.2) return

       isym = 1
       irtn = 1                                         ! stabilite par rapport
       call chol(nco,nco,dmudn,x,ttt,x,isym,istab,irtn) ! a la diffusion
!c                                                      ! istab = 0, stable
       return                                           ! istab = 1, instable
       end
!c
!c**********************************************************************
!c
       subroutine phasLorraine(j,ider,istab,p00,x,eta0,v0,g0,i3r,err0)
!c
!c  calcul des proprietes d'un systeme suppose monophasique, avec si
!c  si istab=1 on retient la solution la plus stable
!c  si istab=2 on etudie la stabilite par rapport a la diffusion.
!c  si istab=10 (on veut choisir sa solution)
!c
       use comflash
       implicit none
!c
       integer nchol
       parameter ( nchol=ncomax*(ncomax+3)/2 )
!c
!c*** variables de la liste d'appel
!c
       integer j,ider,istab,j_ini,err0
       real*8 p00,eta0,v0,g0
       real*8 x(ncomax)
!c
!c*** variables internes
!c
       integer i,i3r,isym,irtn
       real*8 eta1,g,eta_in,eta2
       real*8 ttt(nchol)
	   LOGICAL FORT
!c-------------------------------------------------
	   FORT=.TRUE. 	! ENABLES DEBUGGING MESSAGES WRITTEN TO CONSOLE JRE 20210517. "FORT" = "LOUD" IN FRENCH.
	   FORT=.FALSE.
       err0 = 0
       j_ini = j
!c      Si perte de precision de l'ordinateur :
       if(x(1).eq.0.d0) then
          x(1) = 1.d-15
          x(2) = 1.d0 - x(1)
       endif

       if(x(2).eq.0.d0) then
          x(2) = 1.d-15
          x(1) = 1.d0 - x(2)
       endif


       if(istab.eq.10) then
          eta_in = eta0
!c
!c         en sortie, on veut 1 valeur de eta differente
!c
       endif
       n0 = ncomax
       m0 = ncomax
!c
       call regles_melange(ider,8888,x)
       if(j.eq.3) eta = eta0

       if(j.eq.3) j = 1 !Search the liquid root, then the vapor one and search for min g0.

       eta1 = eta0
       call rcubLorraine(i3r,p00,eta1,eta2)

       if(eta2 .lt. 0.0d0)then
          if(FORT)write(*,*) 'RESOL_EOS -- Problem - eta1 and eta2 < 0'
          err0 = 21
          return
       elseif(eta1 .gt. 1.0d0)then
          if(FORT)write(*,*) 'RESOL_EOS -- Problem - eta1 and eta2 > 1'
          err0 = 22
          return
       elseif(eta1.lt.0.0d0 .and. j.eq.2) then
          if(FORT)write(*,*) 'RESOL_EOS -- Problem - vapor root (eta1) < 0'
          err0 = 23
          return
       elseif(eta2.gt.1.0d0 .and. j.eq.1) then
          if(FORT)write(*,*) 'RESOL_EOS -- Problem - liquid root (eta2) > 1'
          err0 = 24
          return
       endif

       if(j.eq.1) eta = eta2 !Liquid
       if(j.eq.2) eta = eta1 !Vapor
       eta0 = eta
       call eqet	! compute Z = 1/(1-eta)-(a/bRT)/[(1-r1*eta)(1-r2*eta)] given eta. => this eta is NOT etaVT.
       call fug

       !i3r = 1 (three real solutions)
       !i3r = 2 (one real and two imaginary solutions). The real solution is eta1
       !i3r = 3 (one real and two imaginary solutions). The real solution is eta2
       if(i3r.eq.1 .or. (j.eq.1 .and. i3r.eq.3) .or. &
     &     (j.eq.2 .and. i3r.eq.2)) then
          i3r = 1
          !Refine with a Newton step
          call rop(j,p00,eta0)
          call fug
          eta0 = eta
       end if
       v0 = bm / eta0  ! is "bm" = bmVT???
       g0 = 0.0D+00
       do i=1,nco
          if(x(i).le.0.d0 .AND. FORT) then
             write(lun1,*) "Erreur x = 0 ou x < 0 (PHAS1)"
             write(lun1,*) "Systeme : ",nom(1)(1:len_trim(nom(1))), &
     &                     " - ",nom(2)(1:len_trim(nom(2)))
             write(lun1,*) "T = ",t0,"P = ",p00
             write(lun1,*) "x1 = ",x(1),"x2 = ",x(2)
          endif
          g0 = g0 + x(i)*(dlog(x(i)) + lphi(i))
       end do

      if(istab.eq.0) goto 10 ! on se contente de la solution
!c                             ! trouvee en premier

!c     Take the remaining solution
       if(j.eq.1) eta = eta1 !Liquid
       if(j.eq.2) eta = eta2 !Vapor

       call eqet
       call fug

       g = 0.0D+00
       do i=1,nco
          g = g + x(i)*(dlog(x(i)) + lphi(i))
       end do

       if(istab.eq.10) then
          if(dabs(eta1-eta_in).gt.0.001) then !on veut eta1
             eta0 = eta
             g0 = g
             v0 = v
          else !on veut eta0
             eta = eta0          ! retour a la premiere solution
             call eqet
             call fug
          endif
          istab = 1
          goto 10
       endif

       if(g.lt.g0) then       ! si la 2eme solution est plus stable
          eta0 = eta
          g0 = g
          v0 = v
       else
          eta = eta0          ! retour a la premiere solution
          call eqet
          call fug
       endif

  10   if(ider.eq.1) call dertp(ider,x)
       j = j_ini
       if(istab.ne.2) return

       isym = 1
       irtn = 1                                         ! stabilite par rapport
       call chol(nco,nco,dmudn,x,ttt,x,isym,istab,irtn) ! a la diffusion
!c                                                      ! istab = 0, stable
       return                                           ! istab = 1, instable
       end

!c
!c**********************************************************************
!c
       subroutine regles_melange(ider,icalltempe,x)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c      Variables qui doivent être spécifiées préalablement :
!c      1) quadra_b (= .TRUE. si expression de b quadratique)
!c      2) Lij(:,:) et puisb (intervenant dans l'expression de b)
!c      3) modele_gE (chaîne de caractères indiquant le choix de la
!c                    fonction gE/RT)
!c      4) Lambda : constante universelle au denominateur du terme
!c                  d'exces de la regle de melange sur le a

!c      N.B. : en sortie de la subroutine, on a
!c         am = paramètre a dans les uSI
!c         bm = paramètre b en cm3/mol
!c         cm = paramètre c en cm3/mol

!c      ider (entier) :
!c          si ider > 0, calcul de daidx(:,:) pour permettre le calcul des
!c                       derivees de lphi dans dertp
!c
!c      icalltempe (entier) :
!c                 Si icalltempe < 8888, appel de tempe
!c                 pour le calcul des grandeurs dépendant de T
!c                 uniquement (a corps pur et ses dérivées,
!c                 éventuellement, kij & Eij et leurs dérivées)
!c
!c                 Si icalltempe = 8888 : pas d'appel de tempe
!c                 Si icalltempe = 0, 1 ou 2 : appel de tempe
!c                 avec icalltempe, l'indice de dérivation
!c                 par rapport a T

!c                 *** icalltempe doit être < 8888 si
!c                 regles_melange est appelée indépendemment
!c                 de phas1 (i.e., en dehors d'un calcul de
!c                 flash) ; c'est typiquement le cas quand
!c                 on veut calculer hM et cPM.
       use comflash
       implicit none

       integer   ider   ! ider=1 s'il faut mettre à jour les ln phi (Input)
       integer   i,j,k,l
       integer   icalltempe ! < 8888 pour mettre à jour propriétés dépendant
!c                             de T ou pour calculer prop dérivées

       real*8    x(ncomax) ! (Input)
       real*8    bxen,aapur,sumb(ncomax),unspuisb, truc
!c      Pour règles de mélange UNIQUAC-WILSON
       real*8    sompsi(ncomax),psi(ncomax),psi_sur_x(ncomax)
       real*8    lnsompsi(ncomax)
       real*8    mux(ncomax),som,som_mux,terme2,terme3
       real*8    dsompsi_T(ncomax),d2sompsi_T(ncomax),dE_dx(ncomax)
!c      Pour règles de mélange NRTL
       real*8    dsom_ExpAlphTau_dT(ncomax),som_ExpAlphTau(ncomax)
       real*8    E_ij(ncomax,ncomax),dEij_dT(ncomax,ncomax),d2Eij_dT2(ncomax,ncomax)
       real*8    dei_dT(ncomax)
!c-----------------------------------------------------------------------------
          !call tempe(icalltempe)
       if(icalltempe<8887) then
          call tempe(icalltempe)
       endif


!c      ======================================================
!c      Règle de mélange sur le paramètre c :
!c      ======================================================
!c      c linéaire est un cas particulier du c quadratique à condition
!c      de choisir cij = 0.5*(c_i + c_j)
!c      On maintient pourtant la possibilité d'utiliser la formulation linéaire
!c      directe (sans passer par la formulation quadratique) car elle
!c      est moins couteuse en temps de calcul
       cx(:) = 0.d0
       cm = 0.0d0
       do i=1,nco
          cx(i) = c(i)*x(i) ! (cm3/mol)
          cm = cm + cx(i)   ! (cm3/mol)
          dcmdx(i) = c(i)    ! d(cm)/dx(i) = c(i)
       end do

!c      ======================================================
!c      Règle de mélange sur le paramètre b :
!c      ======================================================
!c      b linéaire est un cas particulier du b quadratique à condition
!c      de choisir bij = 0.5*(b_i + b_j)
!c      On pourrait donc utiliser une formulation générale (quadratique)
!c      pour traiter le cas du b linéaire.
!c      On maintient pourtant la possibilité d'utiliser la formulation linéaire
!c      directe (sans passer par la formulation quadratique) car elle
!c      est moins couteuse en temps de calcul

       unspuisb = 1.d0/puisb
       bx(:) = 0.d0
       do i=1,nco
          bx(i) = b(i)*x(i) ! cm3/mol
       end do

       bm = 0.0d0   ! en cm3/mol
       dbmdx(:) = 0.d0

       do i=1,nco
          if(quadra_b) then     ! b quadratique (cm3/mol)
             sumb(i) = 0.d0
             do j = 1,nco
!c               -----------------------------------------
!c               Expression de bij (formulation générale)
!c               -----------------------------------------
                bij(i,j) = 0.5d0 * (b(i)**unspuisb + b(j)**unspuisb)
                bij(i,j) = bij(i,j)**puisb * (1.d0-Lij(i,j))

                sumb(i) = sumb(i) + x(j)*bij(i,j)
             enddo
             bm = bm + x(i)*sumb(i)
             dbmdx(i) = 2.d0*sumb(i) ! d(bm)/dx(i)
          else                  ! b lineaire (cm3/mol)
             bm = bm + bx(i)    ! (cm3/mol)
             dbmdx(i) = b(i)    ! d(bm)/dx(i) = b(i)
          endif
       end do

       if(quadra_b) then  ! b quadratique
!c         d(n.bm)/dn(i) = d(bm)/dx(i) - bm
          derbm(1:nco) = dbmdx(1:nco) - bm
          do i=1,nco
             do k=1,nco
                der2bm(i,k) = 2.d0*bij(i,k) - dbmdx(k) ! d derbm(i) / d x(k)
             enddo
          enddo
       else               ! b lineaire
          derbm(1:nco) = dbmdx(1:nco) ! d(n.bm)/d n(i)
!c                                    ! = d(bm)/dx(i) = b(i)
          der2bm(1:nco,1:nco) = 0.d0  ! d derbm(i) / d x(k)
       endif
       bmix = bm*1.d-6/rgas ! uSI
!c
!c      ======================================================
!c      Regle de melange sur la paramètre a
!c      ======================================================
!c
!c      Calcul de somme des x(i).a(i)/[b(i)RT]
!c      --------------------------------------
       aapur = 0.d0
       do i=1,nco
          aapur = aapur + asb(i)*x(i) ! asb = a/bRT
       enddo


       SELECT CASE(modele_gE)
!c         Règle de mélange sur le a :
!c         ~~~~~~~~~~~~~~~~~~~~~~~~~~
!c         (a/RTb) = [somme des x(i) * a(i)/RTb(i)] + E(T,x)
!c         --> E(T,x) = fonction d'exces (adimensionnelle)
!c
          CASE("PPR78 (Van Laar)")
!c            - Cas d'un kij constant
!c            - Cas (E-)PPR78 ...
!c
!c            Fonction d'exces de type Van Laar - Peneloux
!c             --> ee = E(T,x) = -0.5 * double somme des
!c                               x(i).x(j).b(i).b(j).[E(i,j)/RT]
!c                               / somme des x(k).b(k)
             ee = 0.d0
             do i=1,nco
                e(i) = 0.d0
                do k=1,nco
                   bxen = bx(k)*en(k,i)  ! bxen = x(k).b(k).E(k,i)/RT
                   e(i) = e(i) + bxen    ! e(i) = somme des b(k).x(k).E(i,k)/(RT)
                end do
!c             double somme x(i).x(k).b(i).b(k).E(i,k)/(RT)
!c             (sans dimension)
              ee = ee + e(i)*bx(i)
             end do
             ee = -0.5d0*ee/bm ! E(T,x), sans dimension

             if(icalltempe>0 .and. icalltempe<8888) then
!c               Calcul des dérivées p/r à T (pour estimation
!c               de hM / cPM par exemple)
                dE_dT = 0.d0
                do i = 1,nco
                   truc = 0.d0
                   do j = 1,nco
                      truc = truc + bx(j)*dEijsRT(i,j)
                   enddo
                   dei_dT(i) = b(i)*truc*1.0d-06  !Facteur de conversion (cm3) > (m3) : 1.d-6
                   dE_dT = dE_dT + bx(i)*truc
                enddo


!c               les bx(:) sont en cm3/mol ; bm est en cm3/mol
!c               dEijsRT est dans les uSI
!c               Facteur de conversion (cm3) > (m3) : 1.d-6
                dE_dT = -0.5d-6*dE_dT / bm


                if(icalltempe >= 2) then
                   d2E_dT2 = 0.d0
                   do i = 1,nco
                      truc = 0.d0
                      do j = 1,nco
                         truc = truc + bx(j)*d2EijsRT(i,j)
                      enddo
                      d2E_dT2 = d2E_dT2 + bx(i)*truc
                   enddo

!c                  les bx(:) sont en cm3/mol ; bm est en cm3/mol
!c                  d2EijsRT est dans les uSI
!c                  Facteur de conversion (cm3) > (m3) : 1.d-6
                   d2E_dT2 = -0.5d-6*d2E_dT2 / bm
                endif
             endif
!c
!c            Calcul de aa (alpha = am/bmRT) et de ee (fonction d'exces E)
!c            des aai(i) (= d n.alpha / dn(i))
!c            et daidx(k) (= d aai(i) / d x(k))
             aa = ee + aapur ! aapur = sum(xi*ai/(biRT)), sans dimension
             am = aa * (bm*1.d-6) * rgas*t0 ! uSI
             amix = am/(rgas**2*t0) ! uSI

!c            Calcul de d(n.alpha) / d n(i) note aa(i)
!c            --> d(n.alpha) / d n(i)
             do i=1,nco
                aai(i) = asb(i) - (ee*derbm(i) + e(i)*b(i))/bm
             end do

!c            d(alpha)/dx(i)
             dalphadx(1:nco) = asb(1:nco) - (b(1:nco)*e(1:nco) &
     &                                  + ee*dbmdx(1:nco))/bm

             if(ider.ne.0) then ! calcul de fonctions utiles pour les derivees
!c              de ln phi par rapport aux fractions mol.
!c              Calcul des derivees des aai(i) par rapport
!c              aux fractions molaires x(k)
                do i=1,nco
                   do k=1,nco
                      daidx(i,k) = -b(i)*b(k)*bm*en(i,k) &
     &                     + derbm(i)*(b(k)*e(k) + 2.d0*ee*dbmdx(k)) &
     &                     - bm*ee*der2bm(i,k) &
     &                     + dbmdx(k)*b(i)*e(i)
                      daidx(i,k) = daidx(i,k) / bm**2
                   end do
                end do
             endif

             if(icalltempe>0 .and. icalltempe<8888) then
!c               Calcul des dérivées p/r à T (pour estimation
!c               de hM / cPM par exemple)
!c               Calcul de d alpha / dT :
                dalpha_dT = 0.d0
                do i = 1,nco
                   dalpha_dT = dalpha_dT + asb(i)*x(i)*dlna_T(i)
                enddo
!c               Pour rappel aapur = somme des x(i)*asb(i)
                dalpha_dT = dalpha_dT - aapur/t0 + dE_dT

!c               at = d amix / dT = d racine(am/R²/T) / dT (uSI)
                at = dalpha_dT*bmix

                if(icalltempe >= 2) then
!c                  Calcul de d² alpha / dT² :
                   d2alpha_dT2 = 0.d0
                   do i = 1,nco
                      d2alpha_dT2 = d2alpha_dT2 + asb(i)*x(i)*d2lna_T(i)
                   enddo
                   d2alpha_dT2 = d2alpha_dT2 - 2.d0/t0* &
     &                           (dalpha_dT - dE_dT) + d2E_dT2
!c                  att = d at / dT (uSI) :
                   att = d2alpha_dT2*bmix
!c                  Calcul de d² alpha / dTdxi :
                   do i = 1,nco
                      daai_dT(i) = asb(i)*(dlna_T(i)- 1.0d0/t0) &
     &                            - (dei_dT(i) + dE_dT*derbm(i))/bm
                   enddo
                endif
             endif
          CASE("HV (Uniquac-Wilson)")
!c
!c            Fonction d'exces de type Wilson - UNIQUAC
!c             --> ee = E(T,x) = gE_res / Lambda
!c                    gE_res = Wilson ou Uniquac
!c
             ee = 0.d0
             som = 0.d0
             do i=1,nco
                som = som + x(i)*sigma(i)
                mux(i) = mu(i)*x(i)			   ! Wilson => mui = 1; UNIQUAC => mui = qi 
             end do

             som_mux = SUM(mux)
             psi_sur_x(1:nco) = sigma(1:nco)/som   ! Wilson => sigmai = bi - ci & psi(j) = x(j)*sigma(j)/sigmaMix
             psi(1:nco) = x(1:nco)*psi_sur_x(1:nco)

             do i=1,nco
                sompsi(i) = 0.d0
                do j=1,nco
                   sompsi(i) = sompsi(i) + psi(j)*expA(i,j)	! = sum(PHIj*Lamij) where PHIj = xj*bjVT/sum(xi*biVT), biVT=biPR-ciVT
                end do
             end do
             lnsompsi(1:nco) = DLOG(sompsi(1:nco)) !Wilson's eq => GE/RT= -sum( xi*ln(sum(PHIj*Lamij) )

!c            E(T,x), sans dimension
             ee = 0.d0
             DO i = 1,nco
                ee = ee - mux(i)*lnsompsi(i)	   !Wilson's eq => GE/RT= -sum( xi*ln(sum(PHIj*Lamij) )
             ENDDO

             ee = ee / Lambda						!am/bm= ee*RT + sum(xi*ai/bi)

             if(icalltempe>0 .and. icalltempe<8888) then
!c               Calcul des dérivées p/r à T (pour estimation
!c               de hM / cPM par exemple)
                dE_dT = 0.d0
                DO i = 1,nco
                   dsompsi_T(i) = 0.d0
                   DO j = 1,nco
                      dsompsi_T(i) = dsompsi_T(i) + psi(j)* &
     &                               expA(i,j)*dpargE(i,j)
                   ENDDO
                   dE_dT = dE_dT + mux(i)*dsompsi_T(i)/sompsi(i)
                ENDDO
                dE_dT = dE_dT / Lambda

                if(icalltempe >= 2) then
                   d2E_dT2 = 0.d0
                   DO i = 1,nco
                      d2sompsi_T(i) = 0.d0
                      DO j = 1,nco
                         som = 0.d0
                         DO k = 1,nco
                            som = som + psi(k)*expA(i,k)*( &
     &                          d2pargE(i,k) +dpargE(i,k)*( &
     &                          dpargE(i,j) - dpargE(i,k)))
                         ENDDO
                         d2sompsi_T(i) = d2sompsi_T(i) + psi(j)* &
     &                               expA(i,j)*som
                      ENDDO
                      d2E_dT2 = d2E_dT2 + mux(i)*d2sompsi_T(i)/ &
     &                          sompsi(i)**2
                   ENDDO
                   d2E_dT2 = d2E_dT2 / Lambda
                endif
             endif
!c
!c            Calcul de aa (alpha = am/bmRT)
!c            des aai(i) (= d n.alpha / dn(i))
!c            et daidx(k) (= d aai(i) / d x(k))
             aa = ee + aapur ! alpha = am / (bm.RT), sans dimension
             am = aa * (bm*1.d-6) * rgas*t0 ! uSI
             amix = am/(rgas**2*t0) ! uSI

!c            Calcul de d(n.alpha) / d n(i) = aa(i)
!c            --> d(n.alpha) / d n(i)
             do i=1,nco
                terme3 = 0.d0
                do l = 1,nco
                   terme3 = terme3 + mux(l)*expA(l,i)/sompsi(l)
                enddo
                terme2 = psi_sur_x(i)*(terme3 - Som_mux)
                dE_dx(i) = mu(i)*lnsompsi(i) + terme2
                if(ider.ne.0) then ! calcul de fonctions utiles
!c                  pour les derivees de ln phi par rapport aux
!c                  fractions mol.
!c                  Calcul des derivees des aai(i) par rapport
!c                  aux fractions molaires x(k)
                   do k=1,i
                      som = 0.d0
                      do l = 1,nco
                         som = som + mux(l)*expA(l,i)*expA(l,k) &
     &                               /sompsi(l)**2
                      enddo
                      daidx(i,k) = mu(i)*psi_sur_x(k)*(expA(i,k)/ &
     &                      sompsi(i) - 1.d0) &
     &                      -psi_sur_x(k)*terme2 +psi_sur_x(i)* &
     &                      (mu(k)*(expA(k,i)/sompsi(k) - 1.d0) &
     &                       -psi_sur_x(k)*(som-terme3))
                      daidx(i,k) = -daidx(i,k) / Lambda
                      if(k.NE.i) daidx(k,i) = daidx(i,k)
                   end do
                endif
             end do
             dE_dx(1:nco) = -dE_dx(1:nco)/Lambda

             do i=1,nco
                aai(i) = asb(i) + dE_dx(i)
             end do

!c            d(alpha)/dx(i) = d(n.alpha) / dn(i)
             dalphadx(1:nco) = aai(1:nco)

             if(icalltempe>0 .and. icalltempe<8888) then
!c               Calcul des dérivées p/r à T (pour estimation
!c               de hM / cPM par exemple)
!c               Calcul de d alpha / dT :
                dalpha_dT = 0.d0
                do i = 1,nco
                   dalpha_dT = dalpha_dT + asb(i)*x(i)*dlna_T(i)
                enddo
!c               Pour rappel aapur = somme des x(i)*asb(i)
                dalpha_dT = dalpha_dT - aapur/t0 + dE_dT
!c               at = d amix / dT = d racine(am/R²/T) / dT (uSI)
                at = dalpha_dT*bmix

                if(icalltempe >= 2) then
!c                  Calcul de d² alpha / dT² :
                   d2alpha_dT2 = 0.d0
                   do i = 1,nco
                      d2alpha_dT2 = d2alpha_dT2 + asb(i)*x(i)*d2lna_T(i)
                   enddo
                   d2alpha_dT2 = d2alpha_dT2 - 2.d0/t0* &
     &                           (dalpha_dT - dE_dT) + d2E_dT2
!c                  att = d at / dT (uSI) :
                   att = d2alpha_dT2*bmix
                endif
             endif
          CASE("HV (NRTL)")
!c
!c            Fonction d'exces de type NRTL
!c             --> ee = E(T,x) = gE_res / Lambda
!c                    gE_res = NRTL
!c
!c            E(T,x), sans dimension
             ee = 0.d0
             do i=1,nco
                som_ExpAlphTau(i) = 0.d0
                do k=1,nco
                   som_ExpAlphTau(i) = som_ExpAlphTau(i) + &
     &                     x(k) * expA(k,i)
                end do
                do j=1,nco
                   E_ij(i,j) = pargE(j,i)*expA(j,i)/som_ExpAlphTau(i)
                   ee = ee + x(i)*x(j)*E_ij(i,j)
                end do
             end do

             ee = ee / Lambda

             if(icalltempe>0 .and. icalltempe<8888) then
!c               Calcul des dérivées p/r à T (pour estimation
!c               de hM / cPM par exemple)
               dE_dT = 0.d0
               do i=1,nco
                  dsom_ExpAlphTau_dT(i) = 0.d0
                  do k=1,nco
                     dsom_ExpAlphTau_dT(i) = dsom_ExpAlphTau_dT(i) &
     &                         + x(k)*alphgE(k,i)* dpargE(k,i)*expA(k,i)
                  end do
                  do j=1,nco
                     dEij_dT(i,j) = expA(j,i)/som_ExpAlphTau(i)*( &
     &                         dpargE(j,i)*(1.d0-pargE(j,i)*alphgE(j,i)) &
     &                         + pargE(j,i)*dsom_ExpAlphTau_dT(i) &
     &                                               /som_ExpAlphTau(i))
                     dE_dT = dE_dT + x(i)*x(j)*dEij_dT(i,j)
                 end do
               end do
               dE_dT = dE_dT / Lambda

               if(icalltempe >= 2) then
                   d2E_dT2 = 0.d0
                   do i=1,nco
                     som = 0.d0
                     do k=1,nco
                       som = som + x(k)*alphgE(k,i)*expA(k,i)*( &
     &                             d2pargE(k,i) - alphgE(k,i) &
     &                                                  *dpargE(k,i)**2)
                     end do
                     do j=1,nco
                        terme2 = d2pargE(j,i) - &
     &                           alphgE(j,i)*(dpargE(j,i)**2 &
     &                          + pargE(j,i)*d2pargE(j,i)) &
     &                          + pargE(j,i)*som/som_ExpAlphTau(i)
                        terme3 = dsom_ExpAlphTau_dT(i)/som_ExpAlphTau(i) &
     &                          *(dpargE(j,i) + pargE(j,i)* &
     &                         dsom_ExpAlphTau_dT(i)/som_ExpAlphTau(i))
                        d2Eij_dT2(i,j) = dEij_dT(i,j)* &
     &                       (-alphgE(j,i)*dpargE(j,i) + &
     &                        dsom_ExpAlphTau_dT(i)/som_ExpAlphTau(i)) &
     &                       + expA(j,i)/som_ExpAlphTau(i)* &
     &                                               (terme2 + terme3)
                        d2E_dT2 = d2E_dT2 + x(i)*x(j)*d2Eij_dT2(i,j)
                     end do
                   end do
                   d2E_dT2 = d2E_dT2 / Lambda
               endif
             endif
!c
!c            Calcul de aa (alpha = am/bmRT)
!c            des aai(i) (= d n.alpha / dn(i))
!c            et daidx(k) (= d aai(i) / d x(k))
             aa = ee + aapur ! alpha = am / (bm.RT), sans dimension
             am = aa * (bm*1.d-6) * rgas*t0 ! uSI
             amix = am/(rgas**2*t0) ! uSI

!c            Calcul de d(n.alpha) / d n(i) = aa(i)
!c            --> d(n.alpha) / d n(i)
             do i=1,nco
                som = 0.d0
                do j = 1,nco
                   som = som + x(j)*E_ij(i,j)
                enddo
                terme2 = 0.d0
                do l = 1,nco
                   terme3 = 0.d0
                   do j = 1,nco
                      terme3 = terme3 + x(j)*E_ij(l,j)*expA(i,l)/ &
     &                                         som_ExpAlphTau(l)
                   enddo
                   terme2 = terme2 + x(l)*(E_ij(l,i) - terme3)
                enddo
                dE_dx(i) = som + terme2
                if(ider.ne.0) then ! calcul de fonctions utiles
!c                  pour les derivees de ln phi par rapport aux
!c                  fractions mol.
!c                  Calcul des derivees des aai(i) par rapport
!c                  aux fractions molaires x(k)

                   do k=1,i
                      terme2 = 0.d0
                      do l = 1,nco
                         som = 0.d0
                         do j = 1,nco
                            som = som + x(j)*E_ij(l,j)*expA(k,l) &
     &                                                    *expA(i,l)
                         enddo
                         terme2 = terme2 + x(l)/som_ExpAlphTau(l)*( &
     &                            E_ij(l,i)*expA(k,l) + &
     &                              E_ij(l,k)*expA(i,l) - &
     &                               2.d0*som/som_ExpAlphTau(l))
                      enddo

                      terme3 = 0.d0
                      do j = 1,nco
                         terme3 = terme3 + x(j)*( &
     &                            E_ij(i,j)*expA(k,i)/som_ExpAlphTau(i) &
     &                          + E_ij(k,j)*expA(i,k)/som_ExpAlphTau(k))
                      enddo

                      daidx(i,k) = E_ij(i,k)+E_ij(k,i)-terme2-terme3
                      daidx(i,k) = daidx(i,k) / Lambda
                      if(k.NE.i) daidx(k,i) = daidx(i,k)
                   end do
                endif
             end do
             dE_dx(1:nco) = dE_dx(1:nco)/Lambda

             do i=1,nco
                aai(i) = asb(i) + dE_dx(i)
             end do

!c            d(alpha)/dx(i) = d(n.alpha) / dn(i)
             dalphadx(1:nco) = aai(1:nco)
             if(icalltempe>0 .and. icalltempe<8888) then
!c               Calcul des dérivées p/r à T (pour estimation
!c               de hM / cPM par exemple)
!c               Calcul de d alpha / dT :
                dalpha_dT = 0.d0
                do i = 1,nco
                   dalpha_dT = dalpha_dT + asb(i)*x(i)*dlna_T(i)
                enddo
!c               Pour rappel aapur = somme des x(i)*asb(i)
                dalpha_dT = dalpha_dT - aapur/t0 + dE_dT
!c               at = d amix / dT = d racine(am/R²/T) / dT (uSI)
                at = dalpha_dT*bmix

                if(icalltempe >= 2) then
!c                  Calcul de d² alpha / dT² :
                   d2alpha_dT2 = 0.d0
                   do i = 1,nco
                      d2alpha_dT2 = d2alpha_dT2 + asb(i)*x(i)*d2lna_T(i)
                   enddo
                   d2alpha_dT2 = d2alpha_dT2 - 2.d0/t0* &
     &                           (dalpha_dT - dE_dT) + d2E_dT2
!c                  att = d at / dT (uSI) :
                   att = d2alpha_dT2*bmix
                endif
             endif
       END SELECT
!c
!c       read(*,*)
       end subroutine regles_melange

!c
!c**********************************************************************
!c
       subroutine tempe(ideriv)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c      toutes les fois que la valeur de la temperature est modifiee,
!c      call tempe adapte les parametres a cette nouvelle valeur.
!c      ideriv = 0 pour calculer fonctions de la T
!c      ideriv = 1 pour calculer en sus, dérivées 1res des fonctions de la T
!c      ideriv = 2 pour calculer en sus, dérivées 2ndes des fonctions de la T
!c
!c      Variables calculees par ce sous-programme :
!c      - k12, kij(i,j)
!c      - EN(i,j) = Eij/RT (Eij = coeff. d'interaction Van Laar en mol/cm3)
!c      - pvap(i) (estimation Psat par Clapeyron en bar)
!c      - delta(i) = [racine de a_i/(RT)] / b_i (uSI)
!c      - si ideriv = 1 : on calcule aussi dk12_dT
!c      - si ideriv = 2 : on calcule aussi dk12_dT et dk12_dT2
!c
       use comflash
       implicit none
       integer i,j,k,l,tnum,ideriv,mm0,nn0,intemp
       real*8 r7s3,dixsrt,t298st
       real*8 a_function,unst,racT,racRT
       real*8 alpha1(nbgr_tot_max),alpha2(nbgr_tot_max)
       real*8 dsum,dsum2,high,rtuSI
       real*8 delta(ncomax),delta1(ncomax),delta2(ncomax)
       real*8 f0,g0,h0,tk,sum_gr1,sum_gr2
       real*8 a_de_T,der1_a,der2_a, truc,som
       real*8 Eij,dEij_dT,d2Eij_dT2
	   real*8 expArg ! JRE 20210601
	   LOGICAL FORT
!c-------------------------------------------------
	   FORT=.TRUE. 	! ENABLES DEBUGGING MESSAGES WRITTEN TO CONSOLE JRE 20210517. "FORT" = "LOUD" IN FRENCH.
	   FORT=.FALSE.
!c-----------------------------------------------------------------------------
       !if(ideriv<0 .or. ideriv>2) then
       !   IF(FORT)write(lun1,*) 'Argument tempe loufoque : ideriv = ',ideriv
!          write(lun2,*) 'Argument tempe loufoque : ideriv = ',ideriv
       !   return
       !endif

       rt = r*t0 ! t0 en K et r = 83.14 bar.cm3/mol/K
       r7s3 = 7.0d0/3.0d0
       unst = 1.d0/t0
       racT = dsqrt(t0) ! uSI
       racRT = dsqrt(rgas*t0) ! uSI

       t298st = 298.15d0*unst ! uSI
       dixsrt = 10.d0/rt ! unités Péneloux

       tk = t0
!c      Toutes les grandeurs sont exprimées dans le systeme
!c      international

       do i = 1,nco
          a_de_T = a_function(i,tk,0) ! a(T)
          der1_a = a_function(i,tk,1) ! da/dT
          der2_a = a_function(i,tk,2) ! d2a / dT2
          if(ideriv>=1) dlna_T(i) = der1_a / a_de_T
          if(ideriv>=2) d2lna_T(i) = der2_a / a_de_T

!c         Calcul de ac0, ac1 et ac2 utilisées pour calculs de hM / cPM :
!c         -------------------------------------------------------------
          ac0(i) = dsqrt(a_de_T / tk) / rgas ! racine(a / R²T)
          if(ideriv>=1) then
             ac1(i) = (der1_a*tk - a_de_T)/(2.d0*ac0(i)*tk**2*rgas**2) ! d ac0(i) / dT
             if(ideriv>=2) then
                ac2(i) = (-der1_a**2 * tk**2 - 2.d0*der1_a*a_de_T*tk  & ! d ac1(i) / dT
     &              + 3.d0*a_de_T**2 + 2.d0*a_de_T*tk**2*der2_a) &
     &              / (4.d0*a_de_T*tk**3*ac0(i)*rgas**2)
             endif
          endif
       enddo
!c
!c calcul de asb et pvap pour les constituants purs
!c
       do i=1,nco
          pvap(i) = 10.0d0**(r7s3*(fac(i) + 1.0d0) &
     &              *(1.0d0 - tc(i)*unst))*pc(i)

!c         -- calcul de Racine[a(i)]/b(i) = delta(i)
          delta(i) = ac0(i) * racT / bc(i) ! uSI
!c         -- calcul de Racine[a(i)/RT] / b(i) = delta_sur_RacRT(i)
          delta_sur_RacRT(i) = delta(i) / racRT ! uSI
          delta_sur_RacRT(i) = 1.d-3 * delta_sur_RacRT(i) ! en (cm3/mol)^(-0.5)

!c         -- calcul de delta1(i) = d(delta(i)) / dT
          if(ideriv >= 1) then
             delta1(i) = (ac0(i) + 2.d0 * tk * ac1(i)) &
     &       / (2.d0*racT*bc(i)) ! uSI
          endif

!c         -- calcul de delta2(i) = d(delta1(i)) / dT
          if(ideriv >= 2) then
             delta2(i) = (tk * ac2(i) + ac1(i) - ac0(i) / (4.d0*tk)) &
     &               / (racT*bc(i)) ! uSI
          endif
          asb(i) = ac0(i)**2/bc(i) ! a/(bRT) du corps pur i [sans dimension]
       end do
!c
!c recherche pour savoir si un calcul de en a deja ete fait
!c a la temperature demandee
!c
!c      A l'issue du paragrphe ci-dessous, tnum vaut 0 si
!c      nouvelle température ; est > 0 si T déjà en mémoire
       tnum = 0 ! init
       intemp = itempc
       if(itempc.gt.tempcmax)intemp = tempcmax
       do i=1,intemp ! itempc = nb de T en mémoire
          if(dabs(t0-tmemc(i)).lt.1.d-13) tnum = i
       end do
       if(tnum.eq.0.and.itempc.lt.tempcmax) then
          itempc = itempc + 1
          tmemc(itempc) = t0
       elseif(itempc.ge.tempcmax)then
          itempc = itempc + 1
       endif

       if(modele_gE == "PPR78 (Van Laar)") then
!c        ================================================================
!c         Cas d'un gE de Van Laar :
!c                 - kij constant
!c                 - kij(T) et Eij(T) (PPR78 et compagnie)
!c        ================================================================
!c
!c        calcul des parametres EN(i,j) (coeff. d'interaction de Van Laar)
!c        et des kij :
!c
          kij(1,1) = 0.d0
          kij(2,2) = 0.d0

          if(k12_const>8887.d0) then ! si kij(T)
!c            calcul de alpha1(j) = proportion de groupes j dans la molecule 1
!c            Ce calcul est utile pour le calcul de kij(T) et de ses dérivées
!c
!c            -- calcul de sum_gr1 : nb de groupes total dans molecule 1
             sum_gr1 = 0.d0
             do i = 1,ngr1
                sum_gr1 = sum_gr1 + occurgr1(i)
             end do

!c            -- calcul de alpha1
             alpha1 = 0.d0
             do i = 1,ngr1
                alpha1(nomgr1(i)) = occurgr1(i) / sum_gr1
             enddo

!c            calcul de alpha2(j) = proportion de groupes j dans la molecule 2
!c            -- calcul de sum_gr2 : nb de groupes total dans molecule 2
             sum_gr2 = 0.d0
             do i = 1,ngr2
                sum_gr2 = sum_gr2 + occurgr2(i)
             end do

!c            -- calcul de alpha2
             alpha2 = 0.d0
             do i = 1,ngr2
                alpha2(nomgr2(i)) = occurgr2(i) / sum_gr2
             enddo
          endif

          if(tnum == 0) then ! si T pas en mémoire
!c             N.B. : k12_const = 8888 si kij(T)
!c             sinon, k12_const = valeur du kij constant.
             if(k12_const.lt.8887.d0) then ! ******** Cas d'un kij constant
                k12 = k12_const
                Eij = (delta(1) - delta(2))**2 &
     &                     + 2.d0*k12*delta(1)*delta(2)     ! uSI
                EN(1,2) = (delta_sur_RacRT(1)-delta_sur_RacRT(2))**2 &
     &                + 2.d0*k12*delta_sur_RacRT(1)*delta_sur_RacRT(2)
                EN(2,1) = EN(1,2)
                if(itempc.le.tempcmax) ENT(itempC,1,2) = EN(1,2)
                kij(1,2) = k12
                kij(2,1) = k12
                if(itempc.le.tempcmax) kijT(itempC,1,2) = kij(1,2)
             else
!c               -- calcul double somme par contributions de groupes :
!c               N.B. : les agr0 et bgr0 sont exprimés en Pa
!c                      dsum est donc en Pa
                dsum = 0.d0
                do mm0 = 1,nbgr_used
                   do nn0 = 1,nbgr_used
                      k = liste_gr(mm0)
                      l = liste_gr(nn0)
                      if(agr0(k,l).ne.0.d0) then
                         high = bgr0(k,l) / agr0(k,l) - 1.d0
                         dsum = dsum + (alpha1(k) - alpha2(k)) &
     &                           * (alpha1(l) - alpha2(l)) &
     &                           * t298st**high &
     &                           * agr0(k,l)


                      end if
                   enddo
                enddo

!c             -- calcul des fonctions f0, g0 et h0 telles que :
!c              f0 = E12
!c              g0 = [delta(1)-delta(2)]
!c              h0 = [2 * delta(1) * delta(2)]

                f0 = -0.5d0 * dsum ! E12 en Pa
                g0 = (delta(1)-delta(2))**2
                h0 = 2.d0 * delta(1) * delta(2)
                k12 = (f0 - g0) / h0
                kij(1,2) = k12
                kij(2,1) = k12
                if(itempc.le.tempcmax) kijT(itempC,1,2) = kij(1,2)

!c               Calcul des E12/RT en (cm3/mol)^{-1}
                en(1,2) = f0 * 1.d-6 * dixsrt
                en(2,1) = en(1,2)
                if(itempc.le.tempcmax) ent(itempc,1,2) = en(1,2)

             endif
          else ! temperature déjà considérée (effet mémoire)
             en(1,2) = ent(tnum,1,2)
             en(2,1) = en(1,2)

             kij(1,2) = kijt(tnum,1,2)
             kij(2,1) = kij(1,2)
             k12 = kij(1,2)
          endif ! fin boucle tnum = 0

!c         -----------------------------------------
!c         Calcul de d (E12/RT) / dT = dEijsRT(1,2)
!c         Si ideriv = 2 : on calcule de plus :
!c            d² (E12/RT) / dT² = d2EijsRT(1,2)
!c         -----------------------------------------

          if(ideriv >= 1 ) then
             rtuSI = rgas*t0
             dEijsRT(1,1) = 0.d0
             dEijsRT(2,2) = 0.d0
             dsum = 0.d0
             if(ideriv >= 2) then
                d2EijsRT(1,1) = 0.d0
                d2EijsRT(2,2) = 0.d0
                dsum2 = 0.d0
             endif
             if(k12_const>8887.d0) then
                    do mm0 = 1,nbgr_used
                       do nn0 = 1,nbgr_used
                          k = liste_gr(mm0)
                          l = liste_gr(nn0)
                          if(agr0(k,l).ne.0.d0) then
                             high = bgr0(k,l) / agr0(k,l) - 1.d0
                             truc = (alpha1(k) - alpha2(k)) &
     &                               * (alpha1(l) - alpha2(l)) &
     &                               * t298st**high &
     &                               * agr0(k,l) * (high + 1.d0)

                             dsum = dsum + truc
                             if(ideriv >= 2) then
                                dsum2 = dsum2 + truc * (high + 2.d0)
                             endif
                          end if
                       enddo
                    enddo

                    dEijsRT(1,2) = 0.5d0 * dsum / (rtuSI*t0) ! uSI
                    dEijsRT(2,1) = dEijsRT(1,2)

                    if(ideriv >= 2) then
                       d2EijsRT(1,2) = -0.5d0 * dsum2 / (rtuSI*t0**2) ! uSI
                       d2EijsRT(2,1) = d2EijsRT(1,2)
                    endif
             else
                 som = 0.d0
                 do i = 1,nco
                    som = som + delta(i)*delta1(i)
                 enddo
                 dEij_dT = 2.d0*(som-(1.0d0 - k12)*( &
     &                delta(1)*delta1(2)+delta(2)*delta1(1)))
                 dEijsRT(1,2) = (dEij_dT - Eij/t0) / rtuSI ! uSI
                 dEijsRT(2,1) = dEijsRT(1,2)

                 if(ideriv >= 2) then
                    som = 0.d0
                    do i = 1,nco
                       som = som + delta1(i)**2 + delta(i)*delta2(i)
                    enddo
                    truc = delta(1)*delta2(2) + 2.d0*delta1(1)*delta1(2) &
     &                    + delta(2)*delta2(1)
                    d2Eij_dT2 = 2.0d0*(som - (1.0d0 -k12)*truc)
                    d2EijsRT(1,2) = (d2Eij_dT2 - 2.d0*dEij_dT/t0 & ! uSI
     &                           + 2.d0*Eij/t0**2) / rtuSI
                    d2EijsRT(2,1) = d2EijsRT(1,2)
                 endif
             endif


          endif ! fin ideriv >= 1
	elseif(modele_gE == "HV (Uniquac-Wilson)".or. &
	&         modele_gE == "HV (NRTL)") then
		if(modele_gE == "HV (Uniquac-Wilson)")then
			alphgE = 1.0d00
		end if
		DO i = 1,nco
			pargE(i,i) = 0.d0
			expA(i,i) = 1.d0
			DO j = 1,nco
			  IF(i.NE.j) THEN
				  pargE(i,j) = agE0(i,j)/t0 + bgE0(i,j) + cgE0(i,j)*t0
				  expArg= -alphgE(i,j) * pargE(i,j)
				  if(expArg > 222)expArg=222 ! JRE 20210601: prevent exponential overflow. This may happen for bad guess of age0()
				  expA(i,j) = DEXP(expArg)
				  if(itempc.le.tempcmax) then
					 pargEt(itempc,i,j) = pargE(i,j)
					 expAt(itempc,i,j) = expA(i,j)
				  end if
			   ENDIF
			ENDDO
		ENDDO
		goto 861 ! JRE 20210531 Disable all this temperature management business. It randomizes the values of age0().
		if(tnum == 0) then ! si T pas en mémoire
		!c            Parametres fonction de T :
			 DO i = 1,nco
				pargE(i,i) = 0.d0
				expA(i,i) = 1.d0
				DO j = 1,nco
				  IF(i.NE.j) THEN
					  pargE(i,j) = agE0(i,j)/t0 + bgE0(i,j) &
		&                             + cgE0(i,j)*t0
				  expArg= -alphgE(i,j) * pargE(i,j)
				  if(expArg > 222)expArg=222 ! JRE 20210601: prevent exponential overflow. This may happen for bad guess of age0()
					  if(itempc.le.tempcmax) then
						 pargEt(itempc,i,j) = pargE(i,j)
						 expAt(itempc,i,j) = expA(i,j)
					  end if
				   ENDIF
				ENDDO
			 ENDDO
		  else ! temperature déjà considérée (effet mémoire)
			 DO i = 1,nco
				pargE(i,i) = 0.d0
				expA(i,i) = 1.d0
				DO j = 1,ncomax
				  IF(i.NE.j) THEN
					  pargE(i,j) = pargEt(tnum,i,j)
					  expA(i,j) = expAt(tnum,i,j)
				   ENDIF
				ENDDO
			 ENDDO
		endif ! fin boucle tnum = 0
861		continue		

          if(ideriv >= 1) then
!c            d / dT :
             DO i = 1,nco
                dpargE(i,i) = 0.d0
                DO j = 1,ncomax
                  IF(i.NE.j) THEN
                     dpargE(i,j) = -agE0(i,j)/t0**2 + cgE0(i,j)
                   ENDIF
                ENDDO
             ENDDO
          endif

          if(ideriv >= 2) then
!c            d^2 / dT^2 :
             DO i = 1,nco
                dpargE(i,i) = 0.d0
                DO j = 1,ncomax
                  IF(i.NE.j) THEN
                      d2pargE(i,j) = 2.d0*agE0(i,j)/t0**3
                   ENDIF
                ENDDO
             ENDDO
          endif
       endif ! fin choix modele_gE
!C
       end subroutine tempe
!c
!c**********************************************************************
!c
       subroutine eqet0
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c      Calcul des constantes universelles des EoS cubiques :
!c      - etac = b/vc (compacite critique), sans unite
!c      - omegab et omegaa (cstes univ. intervenant dans a et b)
!c      - b(i) = covolume du constituant i en (cm3/mol)
!c      - ac(i) en bar.cm6/mol2
!c      - msoave(i) = facteur de forme de Soave (sans unite)
       use comflash
       implicit none
!c
!c*** variables internes
       integer i
       real*8 etac,zc
       real*8 rtsp
!c-----------------------------------------------------------------------------

!c  1) calcul des parametres independants des constituants
!c
       etac = ((1.d0-r1)*(1.d0-r2)**2)**(1.d0/3.d0) + &
     &        ((1.d0-r2)*(1.d0-r1)**2)**(1.d0/3.d0) + 1.d0
       etac = 1.d0 / etac

       zc = 1.d0/(3.d0 - etac*(1.d0+r1+r2))

       omegab = zc*etac

       omegaa = (1.d0-etac*r1) * (1.d0-etac*r2) * (2.d0-etac*(r1+r2)) &
     &          / ((1.d0-etac)*(3.d0-etac*(1.d0+r1+r2))**2)

!c  2) calcul des parametres qui dependent de la nature des constituants
!c     - calcul de ac, b, msoave
!c
       do i=1,nco
          rtsp = r*tc(i)/pc(i) ! pc en bar, tc en K et r = 83.14 bar.cm3/mol/K
          b(i) = omegab*rtsp ! cm3/mol
          ac(i) = omegaa*rtsp*r*tc(i) ! bar.cm6/mol2
          ac(i) = ac(i)*1.d-7 ! unité du systeme international (Pa.m6/mol2)

!c         Calcul de bc(i) = b(i)/R en uSI = K/Pa
!c          r = 83.14 bar.cm3/mol/K et b en cm3/mol
          bc(i) = 1.d-5 * b(i) / r ! uSI

          if(FoncAlpha.EQ.'SOAVE') then
!c         -------------------------
!c         Fonction alpha de SOAVE :
!c         -------------------------
!c         EoS de PR78 :
!c         -----------
            if(fac(i).le.0.491d+00) then
               msoave(i) = 0.37464d+00 + 1.54226d+00*fac(i) &
     &                     - 0.26992d+00*fac(i)**2
            else
               msoave(i) = 0.379642d+00 + 1.48503d+00*fac(i) &
     &                     - 0.164423d+00*fac(i)**2 &
     &                     + 0.016666d+00*fac(i)**3
            endif
!c         -------------------------
!c         -------------------------
          endif

       end do

       end subroutine
!c
!c**********************************************************************
!c
       subroutine rop(j,p00,eta0)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c
!c  pilote de la resolution de l'equation d'etat
!c  call rop(j,p00) provoque la resolution, a temperature constante, de
!c  l'equation peq(eta=b/v) = p00 par la methode de newton (paragraphe
!c  ii-ii.1, page 41 de la these de A. Gramajo)
!c  - p00 est la pression en bar
!c  - la temperature a laquelle l'EoS est resolue est t0 (K)
!c  - la composition du fluide est x(:) (fractions molaires)
!c  - l'initialisation de eta depend de la valeur de j
!c    j = 1 : initialisation correspondant au liquide
!c    j = 2 : initialisation correspondant au gaz
!c    j = 3 : initialisation a partir de la valeur de eta se
!c            trouvant deja en memoire (dans com2 du fichier COMFLASH.CMN).
!c
       use comflash
	   use GlobConst, only: dumpUnit
	   use portlib
       implicit none
!c
!c*** liste d'appel
       integer j,nech
       parameter(nech=5)
       real*8 p00
	   real*8 bestEta, bestErr !JRE
!c
!c*** variables internes
       integer nit,ij,j_in,kk
       real*8 eta0
       logical echec(nech)
	   LOGICAL FORT
!c-------------------------------------------------
	   FORT=.TRUE. 	! ENABLES DEBUGGING MESSAGES WRITTEN TO CONSOLE JRE 20210517. "FORT" = "LOUD" IN FRENCH.
	   FORT=.FALSE.
!c-----------------------------------------------------------------------------
       j_in = j
       echec(1:nech) = .FALSE.
       if(eta0.lt.0.0d0 .or. eta0.gt.1.0d0) echec(1) = .true.
       ij = 0
	   bestEta=0.05D0				!JRE
	   if(j==1)bestEta=1-1.D-9		!JRE
	   bestErr = 12345				!JRE
 10    nit = 1
       goto(1,2,3), j

   1   continue
       if(.NOT.echec(1)) then
          eta = eta0
       elseif(.NOT.echec(2)) then
          eta = 0.999d0
       elseif(.NOT.echec(3)) then
          eta = 0.9999d0
       elseif(.NOT.echec(4)) then
          eta = 0.99999d0
       else
          eta = 0.999999d0
       endif
       goto 3
!c
   2   continue
       if(.NOT.echec(1)) then
          eta = eta0
       elseif(.NOT.echec(2)) then
          eta = 1.d-3
       elseif(.NOT.echec(3)) then
          eta = 1.d-4
       elseif(.NOT.echec(4)) then
          eta = 1.d-5
       else
          eta = 1.d-6
       endif

   3   continue
       if( eta < 0 .or. eta > 1)then !JRE 20210513
            write (dumpUnit,*) 'rop1: eta violates 0 < eta < 1 before eqet'
            continue
       endif
       call eqet
       if(nit.eq.1) goto 4

       if(nit.eq.50) then           ! echec de l'initialisation essayee
          if(ij.eq.2) then          ! echec definitif de la resolution
             do kk = 1,nech
                if(.NOT.echec(kk)) then
                   echec(kk) = .TRUE.
                   exit
                endif
             enddo
             if(echec(nech)) then
                j = 0                  ! echec definitif de la resolution
                eta = dabs(eta)
                z = dabs(z)
                qpseta = dabs(qpseta)
                y = dabs(y)

                IF(FORT)write(lun1,141) nom(1)(1:len_trim(nom(1))), &
     &                          nom(2)(1:len_trim(nom(2)))
                IF(FORT)write(lun1,*) 'P/bar = ',p00
                IF(FORT)write(lun1,*) 'T/K   = ',t0

 141            format(1x,'EOS NON SOLUBLE - systeme ',a,' - ',a)
!c                stop
                return   !sortie OK.
             endif
             ij = 0
             j = j_in
             goto 10
          else
             if(j.eq.3) then        ! changement d'initialisation
                j = 1
                ij = 1
             else
                j = 3 - j
                ij = 2
             endif
             goto 10
          endif
       endif
!c
       if(nit.eq.48) then  !il se peut que epseq soit un peu trop petit
          if (dabs(eta0-eta).lt.1000.0d+00*epseq) then !on converge sur eta
             if (eta.lt.0.00d+00.or.eta.gt.1.00d+00) then
                nit = 49
             else                  !eta est bon - on regarde si p est bonne.
                if(dabs(p00 - peq).ge.1.0d-08) then
                   nit = 49
                else
                   return     !   sortie favorable
                endif
             end if
          end if   
       endif

       if (dabs(eta0-eta).lt.epseq) then !on converge sur eta
          if (eta.lt.0.00d+00.or.eta.gt.1.00d+00) then
             nit = 49
          else                  !eta est bon - on regarde si p est bonne.
             if(p00<=1000.d0 .AND. dabs(p00 - peq).ge.1.0d-10) then
                nit = 49
             elseif(p00>10000.d0 .AND. dabs(p00 - peq).ge.1.0d-6) then
                nit = 49
             elseif(p00>1000.d0 .AND. p00<=10000.d0 .AND. &
     &              dabs(p00 - peq).ge.1.0d-8) then
                nit = 49
             else
                return     !   sortie favorable
             endif
          end if
       end if

   4   eta0 = eta
		if( dabs(p00-peq) < bestErr)then
			bestErr=dabs(p00-peq)
			bestEta=eta
		endif

       eta = eta + (p00 - peq)/dpdeta
       if( eta < 0 .or. eta > 1)then !JRE 20210513
           !pause 'rop2: eta violates 0 < eta < 1 after Newton step.'
           eta=bestEta/(1+0.2*random(0)) !this will just keep the iteration going until it maxes out the iterations. 
        endif
!       print*,'nit = ',nit,'eta = ',eta,'z = ',z, &
!    &     'f = ',dabs(eta0-eta)
!       if(eta.gt.1.00d0.and.j.eq.2) eta = 0.10d0
       nit = nit + 1
       goto 3  ! JRE: cycle to the beginning with new guess for eta.

       end subroutine rop
!c
!c**********************************************************************
!c
       subroutine eqet
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c
!c  utilisation de l'equation cubique generale
!c
       use comflash
       implicit none

!c      Forme de l'EoS cubique consideree :
!c      -----------------------------------
!c      z = 1/(1-eta) - alpha * eta * qpseta [avec alpha = a/(bRT)]
!c
!c
       if(eta.eq.1.0D+00) eta = 1.0d0 - 1.0d-10
!c
!c      Terme repulsif de z = v/(v-b) = 1/(1-eta) avec eta = b/v
       y = 1.0d0/(1.0d0 - eta)

!c      Fraction rationnelle en v de la partie attractive de z :
       qpseta = 1.d0 / ((1.d0 - r1*eta)*(1.d0 - r2*eta))

!c      aa = a/(RTb) = alpha
!c      z = v/(v-b) - a.v/[RT.(1-r1.v).(1-r2.v)]
!c        = 1/(1-eta) - alpha.eta/[(1-r1.eta)(1-r2.eta)]
       z = y - aa*eta*qpseta

!c      Arrangement douteux pour eviter probleme dans calcul de LN Phi
        z = dabs(z)

!c      dzdeta = d(z) / d(eta) a T et compo fixees
       dzdeta = y**2 - aa*(1.d0-r1*r2*eta**2)*qpseta**2

       rtsb = rt/bm
       peq = z*rtsb*eta
       dpdeta = (dzdeta*eta + z)*rtsb
       return
       end subroutine eqet
!c
!c**********************************************************************
!c
      subroutine fug
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c     Calcul des coeff. de fugacite par EoS cubique
!c     La subroutine EQET doit avoir ete prealablement appelee
!c     (pour le calcul de z, y, etc.)
!c
       use comflash
       implicit none
       integer i
       real*8 lphi0
!c
       v = bm/eta
       q = 1.d0/(r1-r2)*DLOG((1.d0-r2*eta)/(1.d0-r1*eta))
!c
!c      y = 1 / (1 - eta) = v / (v - bm)
       lphi0 = dlog(y/z) ! -Ln(z) - Ln[(v-bm)/v] = -lnZ - ln(1-bm*rho) = -lnZ-ln(1-B/Z)= -lnZ-ln(Z-B)+lnZ=ln(Z-B)

       do i=1,nco
          lphi(i) = lphi0 + (z-1.d0)*derbm(i)/bm -q*aai(i)
       end do
!c      print*,'phi = ',dexp(lphi(1:2))
       end
!c
!c**********************************************************************
!c
       subroutine rcub(i3r,eta1)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c
!c  connaissant une racine eta de l'equation du 3eme degre
!c  on calcule celle, eta1, qui ne correspond pas au systeme instable
!c
       use comflash
       implicit none
!c
!c*** liste d'appel
!c
       integer i3r
       real*8 eta1,coefa,coefb,coefc
!c
!c*** variables internes
!c
       real*8 pbsrt,bbb,ccc,del,eta2
!c-----------------------------------------------------------------------------

       pbsrt = peq/rtsb

!c      L'EoS cubique s'ecrit :
!c      coefa.eta^3 + coefb.eta^2 + coefc.eta -pbsrt = 0
!c
       coefa = pbsrt*r1*r2 + r1*r2 + aa ! aa = +alpha
       coefb = -pbsrt*(r1*r2 + r1 + r2) - r1 - r2 - aa
       coefc = pbsrt*(r1 + r2 + 1.d0) + 1.d0

!c      on re-ecrit l'EoS cubique sous la forme :
!c      (eta - eta connu).(coefa.eta^2 + bbb.eta + ccc) = 0
!c
       bbb = coefb + coefa*eta
       ccc = coefc + bbb*eta
       del = bbb**2 - 4.0d0*coefa*ccc
       if(del.le.0.0d0) then
          i3r = 0                        ! l'equation n'a qu'une racine
       else
          i3r = 1                        ! l'equation a trois racines
          del = dsqrt(del)
          eta1 = 0.5d0*(-bbb + del)/coefa
          eta2 = 0.5d0*(-bbb - del)/coefa
          if(dabs(eta2-eta).gt.dabs(eta1-eta)) eta1 = eta2 ! elimination
!c                        ! de la racine correspondant au systeme instable
          if(eta1.lt.0.0d0.or.eta1.gt.1.0d0) i3r = 0
       endif
       return
       end
!c
!c**********************************************************************
!c
       subroutine rcubLorraine(i3r,p00,eta1,eta2)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c
!c  connaissant une racine eta de l'equation du 3eme degre
!c  on calcule celle, eta1, qui ne correspond pas au systeme instable
!c
       use comflash
       implicit none
!c
!c*** liste d'appel
!c
       integer i3r
       real*8 eta2
       real*8 eta1,coefa,coefb,coefc
!c
!c*** variables internes
!c
       real*8 pbsrt,A0,A1,A2,QQ,RR,TEST,S1,S2
       real*8 cur,pi,ZZ(3),ADD1,ADD2,theta,p00
!c-----------------------------------------------------------------------------
       pi = 4.0d0 * datan(1.0d0)

       peq = p00
       rtsb = rt / bm
       pbsrt = peq/rtsb

!c      L'EoS cubique s'ecrit :
!c      coefa.eta^3 + coefb.eta^2 + coefc.eta -pbsrt = 0
!c
       coefa = pbsrt*r1*r2 + r1*r2 + aa ! aa = +alpha
       coefb = -pbsrt*(r1*r2 + r1 + r2) - r1 - r2 - aa
       coefc = pbsrt*(r1 + r2 + 1.d0) + 1.d0

       !this part of the code is used when phas1 did not find the desired root
       A2 = coefb / coefa / 3.0d0
       A1 = coefc / coefa
       A0 = - pbsrt / coefa

       QQ = A1 / 3.0D0 - A2 * A2

       RR = (A1 * A2 - A0) / 2.0D0 - A2 * A2 * A2

       TEST = QQ * QQ * QQ + RR * RR
       if(TEST.lt.0.0d0) then
          TEST = DSQRT(-TEST)
          QQ = 2.0D0 * cur(DSQRT(RR*RR+TEST*TEST))
          THETA = pi / 2.0d0
          IF(DABS(RR) > 1.0d-10) THETA = ATAN(TEST/RR)
          IF(THETA < 0) THETA = THETA + PI
          THETA = THETA / 3.0D0
          ADD1 = 2.0d0 * pi / 3.0d0
          ADD2 = 2.0d0 * ADD1


          ZZ(1) = QQ * COS(THETA) - A2
          ZZ(2) = QQ * COS(THETA+ADD1) - A2
          ZZ(3) = QQ * COS(THETA+ADD2) - A2

          i3r = 1
          !rop found a root eta. We choose eta1 in such way to take the remaining one
!          if(dabs(eta-zmin).lt.1.0d-12) eta1 = zmax
       else
           !Test > 0 !One real and two imaginary roots
           TEST = DSQRT(TEST)
           S1 = cur(RR+TEST)
           S2 = cur(RR-TEST)

           ZZ(1) = S1 + S2 - A2
           ZZ(2) = -0.5d0*(S1+S2) - A2 !eta1 is assumed to be the real part of the
                                       !imaginary solutions
           ZZ(3) = ZZ(2)
           i3r = 2
           if(TEST > 1.0d-12) goto 10
           !Test = 0          !Three real solutions with, at least, two identical solutions
           ZZ(1) = S1 + S2 - A2
           ZZ(2) = - S1 - A2
           ZZ(3) = ZZ(2)
           i3r = 1
!           if(dabs(eta-eta1).lt.1.0d-12) eta1 = eta2
       end if

  10   eta1 = minval(ZZ)
       eta2 = maxval(ZZ)

       if(i3r.eq.2)then
          if(dabs(eta1 - ZZ(2)).lt.1.0d-12) i3r = 3
       end if

       return
       end
!c
!c**********************************************************************
!c
       subroutine dertp(ider,x)
!c      -------------------
!c      SUBROUTINE NETTOYEE
!c      -------------------
!c  calcul des derivees de lphi en variables t,p,n  (ider=1)
!c  ou de lphi et p en variables t,v,n (ider=2)
!c
       use comflash
       implicit none

       integer ider
       integer i,k
       real*8 x(ncomax)
       real*8 dzdesz,qp,ertsb,detadp,somd,etasb,dpdvseta
       real*8 somdp,dlphidvseta,dzdt
       real*8 dzdx(ncomax),dpdx(ncomax)
       real*8 dlphidx(ncomax,ncomax)
!c-----------------------------------------------------------------------------
       dzdesz = dzdeta/z ! d ln z / d eta
       qp = eta*qpseta   ! eta / [(1-r1.eta)(1-r2.eta)]
       ertsb = eta*rtsb
       detadp = 1.0d0/dpdeta ! d eta / dP a T et compo fixees

!c
!c  derivees en variables t,eta,x
!c  =============================
!c
       do i=1,nco
!c         qpseta = 1 / [(1-r1.eta)(1-r2.eta)]
!c         y = 1 / (1-eta)
!c         derbm(i) = d(n.bm)/dn(i)
          dlphideta(i) = -dzdesz + y +derbm(i)/bm*dzdeta &
     &                   -aai(i)*qpseta

          dzdx(i) = -dalphadx(i)*qp    ! d(z)/dx(i)
          dpdx(i) = ertsb*(dzdx(i) - z*dbmdx(i)/bm)
       end do

       do i=1,nco
          do k=1,nco
             dlphidx(i,k) = - dzdx(k)/z &
     &                      + derbm(i)/bm*dzdx(k) &
     &                      - derbm(i)*dbmdx(k)*(z - 1.d0)/bm**2 &
     &                      + (z - 1.d0)/bm*der2bm(i,k) &
     &                      - daidx(i,k)*q
          end do
       end do

      if(ider.eq.2) goto 10
!c
!c  derivees en variables t,p,n      ! calcul de flash (ider = 1)
!c  =============================================================
!c
       etasb = eta/bm
       dpdvseta = -etasb*dpdeta ! 1/eta * (dP/dv)
       dpdv = eta*dpdvseta      ! dP/dv
       dpdt = ertsb*(- dalpha_dT * qp + z/t0) ! dP/dT
       detadt = eta*etasb*dpdt/dpdv ! deta / dT
       dzdt = dzdeta*detadt - dalpha_dT * qp ! dz /dT

       do i=1,nco
          dlphidt(i) = dlphideta(i)*detadt  -  & ! d ln phi(i) / dT
     &           dalpha_dT*qp*(derbm(i)/bm - 1.0d0/z) - q*daai_dT(i)

          dlphidp(i) = dlphideta(i)*detadp       ! d ln phi(i) / dP
          somd = 0.0d0
          do k=1,nco
             dlphidx(i,k) = dlphidx(i,k) - dpdx(k)*dlphidp(i)
             somd = somd + x(k)*dlphidx(i,k)
          end do
          do k=1,nco
             dlphidn(i,k) = dlphidx(i,k) - somd ! d ln phi(i) / dnk
             dmudn(i,k) = dlphidn(i,k) - 1.0d0 ! d ln fi(i) / dnk
          end do
          if(x(i).lt.1.0d-300) x(i) = 1.0d-250  !traficos
          dmudn(i,i) = dmudn(i,i) + 1.0d0/x(i)   ! d ln fi(i) / dnk
       end do
       return

!c
!c  derivees en variab     les t,v,n      ! (ider = 2)
!c  ============================================================
!c
 10    etasb = eta/bm
       dpdvseta = -etasb*dpdeta ! 1/eta * (dP/dv)
       dpdv = eta*dpdvseta      ! dP/dv
       dzdt = -dalpha_dT * qp
       dpdt = ertsb*(dzdt + z/t0) ! dP/dT

       somdp = 0.d0
       do i=1,nco
          somd = 0.d0
          dlphidvseta = -etasb*dlphideta(i) ! 1/eta * d ln phi(i) / dv
          dlphidv(i) = eta*dlphidvseta        ! d ln phi(i) / dv
          dlphidt(i) = (derbm(i)/bm - 1.0d0/z)*dzdt  - q*daai_dT(i) ! d ln phi(i) / dT
          do k=1,nco
             dlphidx(i,k) = dlphidx(i,k) - dlphidvseta*dbmdx(k)
             somd = somd + x(k)*dlphidx(i,k)
          end do
          do k=1,nco
            dlphidn(i,k) = dlphidx(i,k) - somd + eta*dlphideta(i) ! d ln phi(i) / dnk
            dmudn(i,k) = dlphidn(i,k) - 1.0d0  ! d ln f(i) / dnk
          end do
          if(x(i).lt.1.0d-300) x(i) = 1.0d-250  !traficos
          dmudn(i,i) = dmudn(i,i) + 1.0d0/x(i)   ! d ln fi(i) / dnk
          dpdx(i) = dpdx(i) - dpdvseta*dbmdx(i)
          somdp = somdp + x(i)*dpdx(i)
       end do

       do i=1,nco
          dpdn(i) = eta*dpdeta + dpdx(i) - somdp
       end do        ! dpdn(i)    = nj*(dp/dnk),

       do i=1,nco
          do k=1,nco
               dmudn(i,k) = dmudn(i,k) + dpdn(k)/peq ! d ln fi(i) / dnk
          enddo
       enddo

       return

       end
!c
!c**********************************************************************
!c
       subroutine chol(n,m,d,f,t,x,isym,itest,irtn)
!c
!c  isym = 1 la matrice  d  est symetrique, on definit  t=d
!c       = 0 la matrice  d  n'est pas symetrique, on calcule  t=dt*d
!c  irtn = 0 on resout le systeme  f=d*x (on calcule t=f ou t=dt*f)
!c       = 1 on triangule uniquement la matrice t (calcul du determinant)
!c       = 2 on inverse la matrice triangulaire t (a la sortie t est la
!c           matrice inverse)
!c
       use comflash
       implicit none
!c
!c*** variables de la liste d'appel
!c
       integer n,m,itest,irtn,isym
       real*8 f(2*nco+1),t(2*nco+1),x(2*nco+1),s
       real*8 d(n0,m0)
!c
!c*** variables internes
!c
       integer i,ii,im1,j,jj,jm1,jp1,k,k1,kk,l,mf,mm1,mp1
       real*8 t1(100)
       logical negatif
!c-------------------------------------------------
       itest = 0
       k = 0
!c
!c  - cas d'une matrice  d  symetrique  (isym=1)
!c
       if(isym.eq.1) then
          do j=1,m
             do k1=1,j
                k = k + 1
                t(k) = d(k1,j)
             end do
          end do
          do j=1,m
             k = k + 1
             t(k) = f(j)
          end do
          go to 200
       endif
!c
!c  - cas d'une matrice  d  non symetrique  (isym=0)
!c
       do j=1,m
          do k1=1,j
             k = k + 1
             t(k) = 0.0d0
             do i=1,n
                t(k) = t(k) + d(i,j)*d(i,k1)
             end do
          end do
       end do
       do j=1,m
          k = k + 1
          t(k) = 0.0d0
          do i=1,n
             t(k) = t(k) + f(i)*d(i,j)
          end do
       end do
!c
!c-------------------------------------------------
!c
!c
!c  triangulation de la matrice  t
!c
!c
 200   mp1 = m + 1
       k = 0
       do i=1,mp1
          ii = i*(i-1)/2
          do j=1,i
             if(j.eq.mp1) goto 29

             k = k + 1
             jm1 = j - 1
             s = 0.0d0
             jj = j*(j-1)/2
             if(jm1.eq.0) goto 22

             do l=1,jm1
                s = s + t(jj+l)*t(ii+l)
             end do

  22         s = t(k) - s
             if(i.eq.j) goto 25

             t(k) = s/t(jj+j)
             goto 29

  25         if(s.gt.0.0d0) goto 26

             if(s.gt.-1.d-13) then
                s = 1.d-16
             else
                if(irtn.ne.1) then
                   itest = 1
                else
                   if(k.lt.m) itest = 1
                endif
                return
             endif

  26         t(k) = dsqrt(s)

  29      end do
       end do

       if(irtn.eq.1) return

!c
!c-------------------------------------------------
!c
!c
!c  inversion de la matrice  t
!c
!c
       if(irtn.eq.2) then
          do i=1,m
             ii = i*(i-1)/2
             t1(ii+i) = 1.0d0/t(ii+i)
             im1 = i - 1
             if(im1.eq.0) go to 32
             do j=1,im1
                s = 0.0d0
                do k=j,im1
                   kk = k*(k-1)/2
                   s = s + t(ii+k)*t1(kk+j)
                end do
                t1(ii+j) = -s/t(ii+i)
             end do
  32      end do

          mf = m*(m+1)/2
          do i=1,mf
             t(i) = t1(i)
          end do
          return
       endif
!c
!c-------------------------------------------------
!c
!c
!c  resolution du systeme  f=t*x
!c
!c

       negatif = .true.
       if(t(k).gt.0.0D+00.and.t(k-m).gt.0.0D+00) negatif = .false.
       if(t(k).lt.0.0D+00.and.t(k-m).lt.0.0D+00) negatif = .false.
       if(t(k).eq.0.0D+00) then
          x(m) = 0.0D+00
       elseif(t(k-m).eq.0.0D+00) then
          x(m) = 1.0D+250
          if(negatif) x(m) = -x(m)
       else
          if((dlog(dabs(t(k)))-dlog(dabs(t(k-m)))).gt.709.0D+00) then
             x(m) = 1.0D+250
             if(negatif) x(m) = -x(m)
          else
             x(m) = t(k)/t(k-m)
          endif
       endif
       if(m.eq.1) return

       mm1 = m - 1
       do l=1,mm1
          j = m - l
          jj = j*(j+1)/2
          jp1 = j + 1
          s = 0.0d0
          do i=jp1,m
             ii = i*(i-1)/2
             s = s + t(ii+j)*x(i)
          end do
          negatif = .true.
          if((t(k-l)-s).gt.0.0D+00.and.t(jj).gt.0.0D+00) negatif=.false.
          if((t(k-l)-s).lt.0.0D+00.and.t(jj).lt.0.0D+00) negatif=.false.
          if((t(k-l)-s).eq.0.0D+00) then
             x(j) = 0.0D+00
          elseif(t(jj).eq.0.0D+00) then
             x(j) = 1.0D+250
             if(negatif) x(j) = -x(j)
          else
             if((dlog(dabs((t(k-l)-s)))-dlog(dabs(t(jj)))).gt.709.0D+00) &
     &       then
                x(j) = 1.0D+250
                if(negatif) x(j) = -x(j)
             else
                x(j) = (t(k-l)-s)/t(jj)
             endif
          endif
       end do

       return
       end
!c
!c-----------------------------------------------------------------------------
!c
      subroutine Cp_Cvec_calc(xcomp,tk,volmol,Cpec,Cvec)
!c     ENTREE : tk, temperature en K
!c              volmol, volume molaire en cm3
!c              xcomp(1:2), fractions mol. de 1 et 2
!c     SORTIE : CP,ec, capacite calo. a P cste d'ecart molaire en J/mol/K
!c
      use comflash
      implicit none

      real*8  xcomp(nco),tk,volmol
      real*8  amel,bmel,damel_T,damel_TT,vm
      real*8  Cvec,betaP,ktv,Cpec
!c-----------------------------------------------------------------------------
!c     Calcul de amix et bmix ainsi que des dérivées at et att :
      t0 = tk
      call regles_melange(0,2,xcomp)

!c     amix = a/(R².T) - unités : SI
      amel = amix * rgas**2*tk  ! a(T) dans les uSI

!C************calcul damel_T = da/dT *************
!c     at = d[a/(R.T)] / dT
!c     R.at = d(a/T) / dT = 1/T.da/dT - a/T
!c     da/dT = R.T.at + a/T = R.(T.at + a/(R.T))
!c     da/dT = R.(T.at + amix)
      damel_T = rgas**2*(tk*at + amix) ! da / dT dans les uSI

!C************calcul damel_TT = d2a / dT2 **********
!c     att = d(at) / dT
!c     R.att = d(1/T.da/dT - a/T) / dT
!c            = 1/T.da/dT - 2/T.da/dT + 2.a/T^3
!c     da/dT= R(T.att + 2/(RT).at - 2/T.a/(R.T))
!c            = R(T.att + 2/(RT).at - 2/T.amix)
      damel_TT = rgas**2*(tk*att + 2.d+0/(rgas**2*tk)*damel_T &
     &           - 2.d+0/tk*amix) ! d2a / dT2 dans les uSI

!c*********************************************
      bmel = bmix*rgas ! m3/mol
      vm = volmol*1.d-6 ! m3/mol

       betaP = rgas/(vm-bmel) - damel_T/(vm-bmel*r1)/(vm-bmel*r2)

       ktv = 1.d+0/(rgas*tk/(vm-bmel)**2 - amel*(2.d+0*vm-(r1+r2)*bmel) &
     &        /((vm-bmel*r1)*(vm-bmel*r2))**2)

       Cvec = -tk/(bmel*(r1-r2)) * damel_TT &
     &        * dlog((vm-bmel*r1)/(vm-bmel*r2))
       Cpec = Cvec - rgas + tk*ktv*betaP**2
!c       write(10,*)xcomp(1:2),Cvec,ktv,betaP,Cpec
      end subroutine Cp_Cvec_calc
!c
!c-----------------------------------------------------------------------------
!c
      subroutine hec_calc(xcomp,tk,volmol,hec)
!c     ENTREE : tk, temperature en K
!c              volmol, volume molaire en cm3
!c              xcomp(1:2), fractions mol. de 1 et 2
!c     SORTIE : hec, enthalpie d'ecart molaire en J/mol
!c
      use comflash
      implicit none
      real*8  xcomp(nco),tk,volmol,hec
      real*8  amel,bmel,damel_T,vm
      real*8  frac

!c     Calcul de amix et bmix ainsi que de la dérivée at :
      t0 = tk
      call regles_melange(0,1,xcomp)

!c     amix = a/(R².T) - unités : SI
      amel = amix * rgas**2*tk  ! a(T) dans les uSI

!c     at = d[a/(R.T)] / dT
!c     R.at = d(a/T) / dT = 1/T.da/dT - a/T
!c     da/dT = R.T.at + a/T = R.(T.at + a/(R.T))
!c     da/dT = R.(T.at + amix)
      damel_T = rgas**2*(tk*at + amix) ! uSI

!c     bmix = b/R
      bmel = bmix*rgas ! m3/mol
      vm = volmol*1.d-6 ! m3/mol

      frac = (vm-bmel*r1)/(vm-bmel*r2)
      hec  = rgas*tk*bmel/(vm-bmel) &
     &       -amel*vm / (vm-bmel*r1) / (vm-bmel*r2) &
     &       + 1.d0/(bmel*(r1-r2)) * (amel-tk*damel_T)*dlog(frac)

!c      print*,'t1 = ',rgas*tk*bmel/(vm-bmel),'t2 = ',
!c     &       -amel*vm / (vm-bmel*r1) / (vm-bmel*r2),'t3 = ',
!c     &       + 1.d0/(bmel*(r1-r2)) * (amel-tk*damel_T)*dlog(frac),
!c     &       'hec = ',hec
!c      write(10,*)xcomp(1:2),hec
      end
!c*****************************************************************************
!   END: MIXTURES SUBROUTINES/FUNCTIONS
!c*****************************************************************************
!c
!c
!c*****************************************************************************
!   AUXILIAR SUBROUTINES/FUNCTIONS
!c*****************************************************************************
!----------------------------------------------------
      function dippr_y(eqId,t_c,coeff, t)
       implicit none
       integer :: eqID
       real(8) :: t, t_c, tr, tt
       real(8) :: y,dippr_y
       real(8) :: a, b, c, d, e, f, g,coeff(7)

       a = coeff(1)
        b = coeff(2)
        c = coeff(3)
        d = coeff(4)
        e = coeff(5)
        f = coeff(6)
        g = coeff(7)


        selectcase(eqID)
            case(100)
                y = a + b * t + c * t**2 + d * t**3 + e * t**4
            case(101)
                y = dexp(a + b / t + c * dlog(t) + d * t**e)
            case(102)
                y = a * t**b / (1.0d0 + c / t + d / t**2)
            case(103)
                y = a + b * dexp(-c / t**d)
            case(104)
                y = a + b / t + c / t**3 + d / t**6 + e / t**9
            case(105)
                y = a / (b**(1.0d0 + (1.0d0 - t / c)**d))
            case(106)
                tr = t / t_c
                y = a * (1.0d0 - tr)**(b + c * tr + d * tr**2 &
     &                      + e * tr**3)
            case(107)
            y = a + b * (c / t / (sinh(c / t)))**2 + &
     &                             d * (e / t / (cosh(e / t)))**2
            case(114)
                tt = 1.0d0 - t / t_c
                y = a**2 / tt + b - 2.0d0 * a * c * tt - &
     &                      a * d * tt**2 - c**2 * tt**3 / 3.0d0 &
     &                      - c * d * tt**4 / 2.0d0 &
     &                      - d**2 * tt**5 / 5.0d0
            case(115)
                y = a + b / t + c * dlog(t) + d * t**2 &
     &                       + e / t**2
            case(116)
                tr = t / t_c
                y = a + b * (1.0d0 - tr)**0.35d0 + &
     &                     c * (1.0d0 - tr)**(2.0d0/3.0d0) + &
     &                     d * (1.0d0 - tr) + &
     &                     e * (1.0d0 - tr)**(4.0d0/3.0d0)
            case(119)
                tt = 1.0d0 - t / t_c
                y = a + b * tt**(1.0d0/3.0d0) + &
     &                     c * tt**(2.0d0/3.0d0) + &
     &                     d * tt**(5.0d0/3.0d0) + &
     &                     e * tt**(16.0d0/3.0d0) + &
     &                     f * tt**(43.0d0/3.0d0) + &
     &                     g * tt**(110.0d0/3.0d0)

            case(120)
                y = a * dexp(b*t) + c * dexp(d*t)
            case(123)
                tt = 1.0d0 - t / t_c
                y = a * (1.0d0 + b * tt**(1.0d0/3.0d0) &
     &             + c * tt**(2.0d0/3.0d0) + d * tt)
            case(124)
                tt = 1.0d0 - t / t_c
                y = a + b / tt + c * tt + d * tt**2 + &
     &               e * tt**3
            case(127)
                y = a + &
     &             b * (c/t)**2 * dexp(c/t) &
     &             /(dexp(c/t) - 1)**2 &
     &             + d * (e/t)**2 * dexp(e/t) &
     &             /(dexp(e/t) - 1)**2 &
     &             + f * (g/t)**2 * dexp(g/t) &
     &             / (dexp(g/t) - 1)**2
            case(150) !IU: K / OU: Torr from DDB
                y =  10.0d0**(a - b / (t + c ))
                y =  y * 101325d0 / 760d0 ! Torr to Pa
            case(151) !IU: K / OU: bar from DDB
                y =  10.0d0**(a - b / (t + c ))
                y =  y * 1.0d05 ! bar to Pa
            case(152) !IU: K / OU: kPa from DDB
                y =  10.0d0**(a - b / (t + c ))
                y =  y * 1.0d03 ! kPa to Pa
            case(153) !IU: K / OU: Pa from DDB
                y =  10.0d0**(a - b / (t + c ))
            case(154) !IU: K / OU: atm from DDB
                y =  10.0d0**(a - b / (t + c ))
                y =  y * 101325d0  ! atm to Pa
            case(155) !IU: °C / OU: Torr from DDB
                y =  10.0d0**(a - b / (t + c - 273.15d0))
                y =  y * 101325d0 / 760d0 ! Torr to Pa
            case(156) !IU: °C / OU: bar from DDB
                y =  10.0d0**(a - b / (t + c - 273.15d0))
                y =  y * 1.0d05 ! bar to Pa
            case(157) !IU: °C / OU: kPa from DDB
                y =  10.0d0**(a - b / (t + c - 273.15d0))
                y =  y * 1.0d03 ! kPa to Pa
            case(158) !IU: °C / OU: Pa from DDB
                y =  10.0d0**(a - b / (t + c - 273.15d0))
            case(159) !IU: °C / OU: atm from DDB
                y =  10.0d0**(a - b / (t + c - 273.15d0))
                y =  y * 101325d0  ! atm to Pa
            case default
                y = huge(1.0d0)
         end select
         dippr_y = y
      end function

!==============================================================================
!c*****************************************************************************
!   END: AUXILIAR SUBROUTINES/FUNCTIONS


!c
!c---------------------------------------------------------------------------
!c
       subroutine flash2(init)
!c
!c  flash diphasique a temperature et pression constante.
!c  methode de substitution suivie par la methode de newton,
!c  avec rattrapage d'echec par application du critere du plan tangent
!c
       use comflash
       implicit none
       integer nchol
       parameter ( nchol=ncomax*(ncomax+3)/2 )
!c
!c*** variables de la liste d'appel
!c
       integer init
!c
!c*** variables internes
!c
       integer jp,newt,istab,ider,istb,jpp
       integer initptg,js,passage,iterAPM
       integer echec,isym,irtn,itend,tend(500)
       integer i,j,k,j0,nit,nitn,iptg,nit0,mono,itest,nitn0,neg,ph

       real*8 xxx,yyy,l00,eps11,eta0,v0,eta1,vv,g,acc,som,dsom,fmax0
       real*8 s1,s2,sk,rok,rok1,rok2
       real*8 sf,sf0,sfn,sfn0,sfd,rap,rapmax,coefrap,rti

       real*8 x0(ncomax),xs(ncomax),etas(2)
       real*8 f(ncomax),ttt(nchol),dn(ncomax),d(ncomax,ncomax)
       real*8 zzz(ncomax),x(ncomax),ki(ncomax),f0(ncomax),test(ncomax)
       real*8 etap1(phmax),etap2(phmax),etap3(phmax)
       real*8 xp1(ncomax,phmax),xp2(ncomax,phmax),xp3(ncomax,phmax)
!c
!c
!c-------------------------------------------------
!c
!c
!c  initialisations preliminaires
!c

       pcmax = maxval(pc(1:nco))
       m0 = ncomax
       n0 = ncomax
       j0 = 1
       newt = 0
       l00 = l0
       acc = 1.0D+00
       nit = 0
       nit0 = 0
       nitn = 0
       nitn0 = 0
       iptg = 0
       itend = 0
       passage = 0
       iterAPM = 0
       eps11 = eps1
       sf = 0.0d0
       sfn = 0.0d0
       do j=1,500
          tend(j) = 0
       end do
       do j = 1,phmax
          etap1(j) = 0.0D+00
          etap2(j) = 0.0D+00
          etap3(j) = 0.0D+00
       end do
       do j = 1,ncomax
          x0(j) = 0.0D+00
          xs(j) = 0.0D+00
          f(j) = 0.0D+00
          dn(j) = 0.0D+00
          x(j) = 0.0D+00
          zzz(j) = 0.0D+00
          ki(j) = 0.0D+00
          f0(j) = 0.0D+00
          test(j) = 0.0D+00
       end do
       do i=1,ncomax
          do j=1,ncomax
             d(i,j) = 0.0d0
          end do
       end do
       do i=1,ncomax
          do j=1,phmax
             xp1(i,j) = 0.0d0
             xp2(i,j) = 0.0d0
             xp3(i,j) = 0.0d0
          end do
       end do
       do j=1,2
          etas(j) = 0.0d0
       end do
       do j=1,nchol
          ttt(j) = 0.0d0
       end do

       xxx = 20.0D+00*eps0
       do i=1,nco
          test(i) = xxx
       end do
!c
!c  nombres de moles et fractions molaires initiales
!c
       ntt = 0.0D+00
       do i=1,nco
          ntt = ntt + nt(i) ! ntt nombre total de moles dans le systeme
       end do

       do i=1,nco
          x0(i) = nt(i)/ntt            ! fractions molaires globales
       end do

       ider = 0
       istb = 1                                   ! calcul de l'enthalpie libre
       call phas1(j0,ider,istb,p0,x0,eta0,v0,gt0) ! gt0 du systeme global
                                                 ! suppose monphasique.
       if(nimp.ge.2) then
          write(lun1,9000) gt0
          write(lun2,9000) gt0
          write(lun1,9010)
          write(lun2,9010)
       endif
 9000  format('phase initiale supposee monophasique, gm/rt=',d17.10)
 9010  format('  nit nit0 nitn',4x,'acc',7x,'fmax',6x,'nl/ng' &
     & ,6x,'etal',6x,'etag  flash2')

       if(init.eq.-1) then
          iph = 1
          goto 500
       endif
       if(init.eq.3) newt = 1
!c
!c
!c-------------------------------------------------
!c
!c
!c*** entree a - initialisations des substitutions
!c
!c
  50   if(init.eq.0) then                 ! calcul des nombres de moles
          do i=1,nco                      ! dans chaque phase a l'initialisation.
             xxx = l00*p0/pvap(i)         ! si init = 1, on a
             np(i,2) = nt(i)/(1.d0 + xxx) ! une initialisation par nombre
             np(i,1) = np(i,2)*xxx        ! de moles et on doit disposer de
          end do                          ! valeurs pour np(i,j).
       endif

       do j=1,2
          npt(j) = 0.0D+00             ! nombres de moles totaux dans
          do i=1,nco                   ! chaque phase.
             npt(j) = npt(j) + np(i,j)
          end do
       end do

       if(init.eq.2) then               ! initialiation par valeurs des
          xxx = npt(2)/npt(1)           ! ki. on doit disposer de
          do i=1,nco                    ! de valeurs initiales pour
            ki(i) = xxx*np(i,1)/np(i,2) ! les nombres de moles des
          end do                        ! deux phases.
          nit0 = 1
          goto 160                      ! on va en b (substitution classique)
       endif
!c
!c
!c-------------------------------------------------
!c
!c
!c *** point de retour pour continuer apres toute iteration reussie
!c
!c
 100   nit = nit + 1
       nit0 = nit0 + 1
       l12 = npt(1)/npt(2)
       do j=1,2
          do i=1,nco                   ! fractions molaires dans chaque
             xp(i,j) = np(i,j)/npt(j)  ! phase
             xp3(i,j) = xp2(i,j)       ! et
             xp2(i,j) = xp1(i,j)       ! memoire
             xp1(i,j) = xp(i,j)        ! des
          end do                       ! solutions
          etap3(j) = etap2(j)          ! precedentes
          etap2(j) = etap1(j)          !
          etap1(j) = etap(j)           !
       end do
!c
!c  critere de solution triviale
!c
       mono = 0
       do i=1,nco
          if(dabs(xp(i,2)-xp(i,1)).lt.1d-4) mono = mono + 1
       end do
       if(mono.eq.nco) then            ! solution triviale.
          if(nimp.gt.1) then
             write(lun1,131)
             write(lun2,131)
 131         format(1x,'solution triviale')
          endif
          goto 400                     ! on va en e (critere du plan tangent)
       endif
!c
!c
!c-------------------------------------------------
!c
!c
!c  calcul des coefficients de fugacite pour toute methode
!c  et de leurs derivees pour la methode de newton
!c
!c
       do i=1,nco
          ki(i) = 0.0D+00
       end do
       gt = 0.0d+0
       if(newt.eq.1) then
          xxx = 20.0D+00*eps0
          do i=1,nco
             test(i) = xxx
                do k=1,nco
                   d(k,i) = 0.0D+00
                end do
          end do
       endif

       do j=1,2
          do i=1,nco
             x(i) = xp(i,j)
          end do
          if(nit0.lt.3.and.init.lt.1) then
             jp = j
          else
             eta = etap(j)
             jp = 3
          endif
          istb = 0
          call phas1(jp,newt,istb,p0,x,etap(j),vp(j),gp(j))
          gt = gt + npt(j)*gp(j)
          do i=1,nco
             ki(i) = ki(i) + lphi(i)*(2.0d0*j - 3.0d0)
          end do
          if(newt.eq.1) then
             xxx = 1.0D+00/npt(j)
             do i=1,nco
                test(i) = test(i) + dabs(dlphideta(i))*epseq
                do k=1,nco
                   d(k,i) = d(k,i) + dmudn(k,i)*xxx
                end do
             end do
          endif
       end do
       gt = gt/ntt

       fmax0 = fmaxi
       fmaxi = 0.0D+00
       itest = 0
       do i=1,nco
          f0(i) = f(i)
          f(i) = dlog(xp(i,2)/xp(i,1)) + ki(i)
          if(ki(i).gt.700.0D+00) ki(i) = 700.0D+00  !traficos JN
          ki(i) = dexp(ki(i))
          if(ki(i).lt.1.0D-14) ki(i)=1.0D-13  !traficos JN
          if(dabs(f(i)).gt.fmaxi) fmaxi = dabs(f(i))
          if(dabs(f(i)).gt.test(i)) itest = itest + 1
       end do
       if(nitn0.eq.0) then
          sf0 = sf
          sf = 0.0D+00
          sfd = 0.0D+00
          do i=1,nco
             sf = sf + f(i)*f(i)
             if(iac.ge.2) sfd = sfd + f0(i)*(f(i) - f0(i))
          end do
          if(iac.eq.1) then
             sfn0 = sfn
             sfn = dsqrt(sf)
          endif
       endif

       if(nimp.ge.2) then
          write(lun1,9050) nit,nit0,nitn,acc,fmaxi,l12,etap(1),etap(2)
          write(lun2,9050) nit,nit0,nitn,acc,fmaxi,l12,etap(1),etap(2)
       endif
 9050  format(3i5,f7.2,2d11.4,2f10.5)

       if(init.ne.3.and.newt.eq.1.and.nitn0.eq.0.and.fmaxi.gt.fmax0) &
     &    then
          newt = 0
          nit0 = nit0 + 1
       endif


       if(newt.eq.1) goto 300          ! on va en d (methode de newton).
!c
!c
!c-------------------------------------------------
!c
!c
!c  methode de substitution simple
!c
!c
       if(oa.eq.2.or.init.eq.1) goto 150 ! option a2. on passe
!c                                        ! en b (substitution classique).
!c
       if(nit0.eq.1) then
          npt(1) = 0.0D+0                  ! option a1: une iteration de
          npt(2) = 0.0D+0                  ! substitution simple.
          do i=1,nco
             xxx = l12*ki(i)
             np(i,2) = nt(i)/(1.0D+00 + xxx)
             np(i,1) = np(i,2)*xxx
             npt(1) = npt(1) + np(i,1)
             npt(2) = npt(2) + np(i,2)
          end do
          goto 100
       endif
!c
!c
!c-------------------------------------------------
!c
!c
!c  entree b ou c
!c  methode de substitution classique
!c
!c  bc - 1: determination du coefficient d'acceleration
!c
!c

 150   if(iac.eq.0.or.nit0.le.2) then
          acc = 1.0D+00
          goto 160
       endif

       if(nit0/2*2.eq.nit0.and.iac.ne.3) then
          acc = 1.0D+00
       else
          if(iac.eq.1) then
             if(sfn0.ne.sfn) then
                acc = dabs(sfn0/(sfn0 - sfn))
             else
                acc = 1.00D+00
             endif
          else
             acc = acc*dabs(sf0/sfd)
          endif
          do i=1,nco
             if(dabs(acc*f(i)).gt.4.6D+0) acc = 1.0D+00
          end do
       endif

       if(acc.ne.1.0D+00) then
          do i=1,nco
             ki(i) = xp(i,1)/xp(i,2)*dexp(acc*f(i))
          end do
       endif
!c
!c
!c  bc - 2: on verifie si le systeme a resoudre pour calculer les nombres
!c          de moles possede une solution
!c

 160   s1 = 0.0D+00
       s2 = 0.0D+00
       do i=1,nco
          zzz(i) = x0(i)*(ki(i) - 1.0D+00)
          s1 = s1 + zzz(i)
          s2 = s2 + zzz(i)/ki(i)
       end do

       if((s1.gt.0.0D+00.and.s2.gt.0.0D+00).or. &
     &    (s1.lt.0.0D+00.and.s2.lt.0.0D+00)) then !ce systeme n'a pas de solution
          if(acc.ne.1.0D+00) then
             do i=1,nco
                if(f(i).gt.700.00D+00) f(i) = 700.0D+00    !traficos JN
                ki(i) = xp(i,1)/xp(i,2)*dexp(f(i))
             end do                ! on examine si on peut le resoudre
             acc = 1.0D+00         ! en supprimant l'acceleration.
             goto 160
          else                     !*** reinitialisation b1
             if(nit0.le.nb+1.and.l00.gt.1d-5.and.l00.lt.1d+3) then
                itend = itend + 1
                if(l12.ge.l0) then
                   l00 = l00*5.0D+00
                   tend(itend) = 1
                else
                   l00 = l00*0.2D+00
                   tend(itend) = 2
                endif
                if(itend.gt.1) then
                   if(tend(itend).ne.tend(itend-1)) then
                      if(nimp.gt.1) then
                         write(lun1,132)
                         write(lun2,132)
 132                     format(1x,'substitution: impossible')
                      endif
                      goto 400  ! on va en e (critere du plan tangent)
                   endif
                endif
                nit0 = 0
                goto 50   ! on va en a avec un nouveau l00
             else
                if(nimp.gt.1) then
                   write(lun1,133)
                   write(lun2,133)
 133               format(1x,'substitution: impossible')
                endif
                goto 400  ! on va en e (critere du plan tangent)
             endif
          endif
       endif
!c
!c
!c bc - 3. cycle de sous-iterations
!c
!c

       iterAPM = 0
 200   som = 0.0D+00
       dsom = 0.0D+00
       iterAPM = iterAPM + 1
       do i=1,nco
          xxx = nt(i)*(ki(i) - 1.0D+00)/(ntt + (ki(i) - 1.0d+0)*npt(1))
          som = som + xxx
          dsom = dsom - xxx*xxx/nt(i)
       end do
       npt(1) = npt(1) - som/dsom

       if(npt(1).le.0.0D+00.or.npt(1).ge.ntt) then
          passage = passage + 1
          if(passage.gt.1) then
            passage = 0
            goto 400
          endif
                                             ! il n'est pas possible
          rok1 = 0.0D+00                      ! de resoudre le systeme a
          rok2 = 1.0D+00                      ! partir de l'initialisation
                                              ! choisie. on en cherche une
 220      rok = 0.5d+0*(rok1 + rok2)          ! autre par dichotomie.
          sk = 0.0D+00
          do i=1,nco
             sk = sk + zzz(i)/(1.0d0 + (ki(i) - 1.0d0)*rok)
          end do
          if(sk*s1.lt.0.0d0) then
             s2 = sk
             rok2 = rok
          else
             s1 = sk
             rok1 = rok
          endif
          if(dabs(rok1 - rok2).gt.1d-7) goto 220

          npt(1) = rok*ntt
          goto 200
       endif
       if(dabs(som).lt.1.d-9 .and. iterAPM .gt.50) goto 229

       if(dabs(som).gt.1.d-10) goto 200 ! critere d'arret des
!c                                       ! sous-iterations
!c
!c  bc - 4. fin du calcul des nombres de moles
!c
!c
 229   npt(2) = ntt - npt(1)
       l12 = npt(1)/npt(2)
       do i=1,nco
          if(ki(i).gt.1.0D+300) ki(i) = 1.0D+300  !traficos JN
          xxx = l12*ki(i)
          np(i,2) = nt(i)/(1.0d0 + xxx)
          np(i,1) = np(i,2)*xxx
       end do
!c
!c
!c  bc - 5. critere de sortie de la methode de substitution
!c
!c
       if(nit0.ge.100 .OR. (nit0>nitflash .AND. p0>2.d0*pcmax)) then
          if(nimp.gt.1) then
             write(lun1,134)
             write(lun2,134)
 134         format(1x,'flash2 - substitution: convergence lente')
          endif
          goto 400
       endif

       if(nit0.lt.4.or.fmaxi.gt.0.1d0) then
          goto 100
       endif
       if(iac.ne.3.and.acc.ne.1.0d0) then
          goto 100
       endif

       if(isubs.eq.1) then             ! methode de substitution
          if(fmaxi.lt.eps11) then      ! utilisee seule.
             iph = 2                   ! systeme diphasique
             goto 500                  ! sortie favorable
          endif
          goto 100                     ! autre iteration
       endif

       if(iptg.eq.1) then
          if(nit0.gt.nitflash.or.fmaxi.le.1d-8) then
             newt = 1
          else
             do i=1,nco
                x(i) = xp(i,1)
             end do
             eta1 = etap(1)
             jpp = 3
             ider = 1
             istab = 2
             call phas1(jpp,ider,istab,p0,x,eta1,vv,g)
             if(istab.eq.1) goto 100

             do i=1,nco
                x(i) = xp(i,2)
             end do
             eta1 = etap(2)
             jpp = 3
             ider = 1
             istab = 2
             call phas1(jpp,ider,istab,p0,x,eta1,vv,g)
             if(istab.eq.0) newt = 1
          endif
       else
          if(isubs.ne.1) newt = 1
          if(cam.eq.1) then
             if(gt.lt.gt0.or.fmaxi.lt.eps11) newt = 1
          else
             if(fmaxi.lt.eps11) newt = 1
          endif
       endif
       if(newt.eq.1) then
          nitn0 = 0
       end if
       goto 100                 ! autre iteration
!c
!c
!c-------------------------------------------------
!c
!c
!c  entree c: methode de newton raphson
!c
!c

 300   if(itest.eq.0) then
          iph = 2               ! systeme diphasique
          goto 500              ! sortie favorable
       endif
!c                              ! non convergence de la methode de newton
       if(nitn0.eq.nitflash) goto 400! on va en e (critere du plan tangent).
!c
       nitn = nitn + 1
       nitn0 = nitn0 + 1
       isym = 1
       irtn = 0
       call chol(nco,nco,d,f,ttt,dn,isym,echec,irtn)
       if(echec.eq.1) then
          if(nimp.gt.1) then
             write(lun1,135)
             write(lun2,135)
 135         format(1x,'newton: echec')
          endif
          goto 400             ! on va en e (critere du plan tangent)
       endif

       xxx = npt(1)
       yyy = npt(2)
       npt(1) = 0.0d0             ! calcul des nouveaux nombres de moles
       npt(2) = 0.0d0
       do i=1,nco
          np(i,1) = np(i,1) + dn(i)
          np(i,2) = np(i,2) - dn(i)
          npt(1) = npt(1) + np(i,1)
          npt(2) = npt(2) + np(i,2)
          if(np(i,1).le.0.0d0.or.np(i,2).le.0.0d0) then
             if(nimp.gt.1) then
                write(lun1,136)
                write(lun2,136)
 136            format(1x,'newton: nombres de moles negatifs')
             endif
             do k=1,i
                np(k,1) = np(k,1) - dn(k)
                np(k,2) = np(k,2) + dn(k)
             end do
             npt(1) = xxx
             npt(2) = yyy
             goto 400          ! on va en e (critere du plan tangent)
          endif
       end do
       goto 100                 ! nouvelle iteration
!c
!c
!c-------------------------------------------------
!c
!c
!c  entree e: critere du plan tangent
!c
!c
 400   if(iptg.eq.1) goto 420 ! le critere du plan tangent a deja ete
!c                             ! applique, on modifie l'initialisation
       iptg = 1
!c
!c
!c  e-1. initialisation de la phase xs (pages 69-70)
!c
!c
       if(nit0.lt.5) then ! si moins de 5 iterations ont ete effectuees
          initptg = 0     ! on se contentera des initialisations raoult
          js = 1
          if(l12.gt.1.0d0) js = 2
       else               ! sinon on prend pour xs initial la phase qui
          initptg = 1     ! 3 iterations avant, est la plus differente
          js = 1          ! de x0
          if(dabs(etap3(2)-eta0).gt.dabs(etap3(1)-eta0)) js = 2
          do i=1,nco
             xs(i) = xp3(i,js)
          end do
          etas(2) = etap3(js)
       endif
!c
!c
!c  e-2. critere du plan tangent
!c
!c
       call ptg(initptg,js,p0,x0,xs,etas)

       if(iph.eq.1) goto 500   ! sortie "systeme monophasique"

       if(etas(2).gt.eta0) then
          js = 1               ! phase "s" liquide
       else
          js = 2               ! phase "s" gaz
       endif
       rapmax = 1d+30
       do i=1,nco
          xxx = x0(i)/(xs(i) - x0(i))
          if(xxx.gt.0.0D+00.and.xxx.lt.rapmax) rapmax = xxx
       end do
       neg = 0
       if(dabs(etas(1)-etas(2)).gt.0.1d0) then
          neg = 0
          rap = 0.5d-6
       else
          neg = 1
          rap = 1.99d0*rapmax
       endif

 420   if(rap.lt.0.1d0.or.rap.gt.10.0d0) then
          coefrap = 2.0d0
       else
          coefrap = 2.0d0
       endif
       if(neg.eq.0) then
          rap = rap*coefrap
       else
          rap = rap/coefrap
       endif

       if(rap.lt.0.99d-6) then
          iph = 0
          goto 500
       endif

       if(rap.ge.rapmax) then
          rap = 0.999d0*rapmax
          neg = 1
       endif

       do i=1,nco
          x(i) = (1.0D+00 + rap)*x0(i) - rap*xs(i)
       end do

       if(neg.eq.0) then
          jpp = 3 - js
          ider = 1
          istab = 2
          call phas1(jpp,ider,istab,p0,x,eta1,vv,g)
          if(istab.ne.0) goto 420
       else
          jpp = 3 - js
          ider = 0
          istb = 1
          call phas1(jpp,ider,istb,p0,x,eta1,vv,g)
       endif

       if(etas(2).gt.eta1) then
          js = 1
       else
          js = 2
       endif
       npt(3-js) = ntt/(1.0d0 + rap)
       npt(js) = rap*npt(3-js)
       do i=1,nco
          np(i,js) = xs(i)*npt(js)
          np(i,3-js) = x(i)*npt(3-js)
       end do
       if(nimp.gt.1) then
          write(lun1,137) rap
          write(lun2,137) rap
 137      format(1x,'retour, rap=',f17.12)
       endif
       if(nimp.ge.2) then
          write(lun1,9010)
          write(lun2,9010)
       endif
       nit0 = 0
       nitn0 = 0
       newt = 0
       goto 100

 500   if (iph.le.1) then
          vp(1) = v0
          vp(2) = 0.0d+00
          gp(1) = gt0
          gp(2) = 0.0d+00
          gt = gt0
          npt(1) = ntt
          npt(2) = 0.0d+00
          densp(2) = 0.0d+00
          zp(2) = 0.0d+00
!c          vmolp(2) = 0.0d+00
          do i = 1,nco
             xp(i,1) = x0(i)
             xp(i,2) = 0.0d+00
             np(i,1) = nt(i)
             np(i,2) = 0.0d+00
          end do
       end if

       do j=1,iph
          mmolp(j) = 0.0d+00
          vph(j) = vp(j)
          do i=1,nco
             mmolp(j) = mmolp(j) + xp(i,j)*mmol(i)
          end do
          densp(j) = mmolp(j)/vph(j)
       end do

!c       if(etap(2).gt.etap(1).and.iph.gt.1) then
       if(densp(2).gt.densp(1)) then
          xxx = npt(1)
          npt(1) = npt(2)
          npt(2) = xxx
          xxx = etap(1)
          etap(1) = etap(2)
          etap(2) = xxx
          xxx = vp(1)
          vp(1) = vp(2)
          vp(2) = xxx
          xxx = gp(1)
          gp(1) = gp(2)
          gp(2) = xxx
          xxx = mmolp(1)
          mmolp(1) = mmolp(2)
          mmolp(2) = xxx
          xxx = zp(1)
          zp(1) = zp(2)
          zp(2) = xxx
          do i=1,nco
             xxx = xp(i,1)
             xp(i,1) = xp(i,2)
             xp(i,2) = xxx
             xxx = np(i,1)
             np(i,1) = np(i,2)
             np(i,2) = xxx
          end do
       endif

       ph = iph
       if(iph.eq.0) ph = 1
       rti = 1.0d+0/r/t0
       vt = 0.0d+00
       mmolt = 0.0d+00
       do j=1,ph
          mmolp(j) = 0.0d+00
          do i=1,nco
             mmolp(j) = mmolp(j) + xp(i,j)*mmol(i)
          end do
          mmolt = mmolt + npt(j)*mmolp(j)
          zp(j) = p0*vp(j)*rti
          densp(j) = mmolp(j)/vp(j)
          vp(j) = npt(j)*vp(j)
          vt = vt + vp(j)
       end do
       denst = mmolt/vt
       mmolt = mmolt/ntt
       zt = p0*vt*rti/ntt

       do i=1,nco
          volat(i) = 0.0D+00
       end do

       if(iph.eq.2) then
          dropout = 1.00d+02*vp(1)/vt   ! % volumique de la phase liquide
          dropout = dabs(dropout)
          do i=1,nco
             volat(i) = xp(i,2)/xp(i,1)   !y/x
          end do
       endif
       fl_ok = .false.
       if(iph.eq.2) then
          if(ygtx.and.xp(1,2).gt.xp(1,1)) fl_ok = .true.
          if(.not.ygtx.and.xp(1,1).gt.xp(1,2)) fl_ok = .true.
       endif

       return
       end
!c
!c---------------------------------------------------------------------------
!c
       subroutine ptg(init,jd0,p00,x0,xs,etas)
!c
!c  analyse de la stabilite de la phase x0 par le critere du plan tangent
!c  pages 68-73
!c
       use comflash
       implicit none

       integer nchol
       parameter (nchol=ncomax*(ncomax+3)/2 )
!c
!c*** variables de la liste d'appel
!c
       integer init,jd0
       real*8 p00,x0(ncomax),xs(ncomax),etas(2)
!c
!c*** variables internes
!c
       integer j0,js,newt,istab,ider
       integer echec,isym,irtn
       integer i,k,j,jd,nit,nitn,itest,neg,nimp0,num,num0

       real*8 v0,g
       real*8 xxx,s,sm1,ly,sf,sfn,sfn0,acc,testmax,accmin,gs,gsmin

       real*8 etasn(3),gsn(3),xsn(ncomax,3)
       real*8 f(ncomax),ttt(nchol),dn(ncomax),d(ncomax,ncomax)
       real*8 f0(ncomax),test0(ncomax),yy(ncomax),y0(ncomax),y1(ncomax),test(ncomax)
       real*8 yjn(ncomax)
	   LOGICAL FORT
!c-------------------------------------------------
	   FORT=.TRUE. 	! ENABLES DEBUGGING MESSAGES WRITTEN TO CONSOLE JRE 20210517. "FORT" = "LOUD" IN FRENCH.
	   FORT=.FALSE.
!c
!c      initialisation des variables internes a zero
!c
       v0 = 0.0d0
       g = 0.0d0
       xxx = 0.0d0
       s = 0.0d0
       sm1 = 0.0d0
       ly = 0.0d0
       sf = 0.0d0
       sfn = 0.0d0
       sfn0 = 0.0d0
       acc = 0.0d0
       testmax = 0.0d0
       accmin = 0.0d0
       gs = 0.0d0
       gsmin = 0.0d0
       do i=1,ncomax
          f(i) = 0.0d0
          dn(i) = 0.0d0
          f0(i) = 0.0d0
          test0(i) = 0.0d0
          yy(i) = 0.0d0
          y0(i) = 0.0d0
          y1(i) = 0.0d0
          test(i) = 0.0d0
          yjn(i) = 0.0d0
          do j=1,ncomax
             d(i,j) = 0.0d0
          end do
          do j=1,3
             xsn(i,j) = 0.0d0
          end do
       end do
       do i=1,3
          etasn(i) = 0.0d0
          gsn(i) = 0.0d0
       end do
!c
!c
!c-------------------------------------------------
!c
!c
       m0 = ncomax
       n0 = ncomax
       num = 0
       nimp0 = 0
       jd = jd0      ! parametre d'orientation de la nature des phases
       js = jd       ! pour la phase xs
       j0 = 3 - js   ! pour la phase x0
!c
!c
!c-------------------------------------------------
!c
!c
!c  proprietes de la phase donnee x0
!c
       ider = 1
       istab = 1
       call phas1(j0,ider,istab,p00,x0,etas(1),v0,g)

       do i=1,nco
          f0(i) = dlog(x0(i)) + lphi(i)
          test0(i) = 10.0D+00*eps0 + dabs(dlphideta(i)*epseq)
       end do
!c
!c
!c-------------------------------------------------
!c
!c  initialisation de la phase recherchee xs
!c
!c
       if(init.eq.1) then  ! on dispose d'une composition xs qui sert
          s = 1.0d0        ! d'initialisation
          do i=1,nco
             yy(i) = xs(i)
          end do
          goto 110
       endif

 100   s = 0.0d0           ! initialisations suivant la loi de raoult
       if(jd.eq.1) then    ! xs supposee liquide
          do i=1,nco
             yy(i) = x0(i)*p00/pvap(i)
             s = s + yy(i)
          end do
       else                ! xs supposee gaz
          do i=1,nco
             yy(i) = x0(i)*pvap(i)/p00
             s = s + yy(i)
          end do
       endif

 110   nit = 0
       nitn = 0
       newt = 0
       if(nimp.eq.3 .AND.FORT) then
          write(lun1,9000)
          write(lun2,9000)
       endif
 9000  format('  nit nitn',4x,'acc',8x,'fmax',9x,'s-1' &
     & ,6x,'eta0',6x,'etas     ptg')
!c
!c
!c-------------------------------------------------
!c
!c  debut des iterations
!c  proprietes de la phase xs
!c
!c
 120   nit = nit + 1
       do i=1,nco
          xs(i) = yy(i)/s
       end do
       if(init.eq.0) then
          j = js
       else
          j= 3
       endif
       istab = 0
       call phas1(j,newt,istab,p00,xs,etas(2),v0,gs)

       if(nimp.eq.3) then
          write(lun1,9100) nit,nitn,acc,fmaxi,sm1,etas(1),etas(2)
          write(lun2,9100) nit,nitn,acc,fmaxi,sm1,etas(1),etas(2)
       endif
9100   format(2i5,f7.2,2d12.2,2f10.5)

       if(newt.eq.1) goto 130  ! vers la methode de newton
!c
!c
!c-------------------------------------------------
!c
!c  methode de substitution
!c
!c
       fmaxi = 0.0d0
       s = 0.0d0
       sf = 0.0d0
       do i=1,nco
          ly = f0(i) - lphi(i)
          f(i) = ly - dlog(yy(i))
          sf = sf + f(i)*f(i)
          if(dabs(f(i)).gt.fmaxi) fmaxi = dabs(f(i))
          y1(i) = y0(i)
          if(ly.gt.700.0D+00) ly = 700.0D+00  !traficos JN
          yy(i) = dexp(ly)
          if(yy(i).lt.1.0D-14) yy(i)=1.0D-13  !traficos JN
          y0(i) = yy(i)
          s = s + yy(i)
       end do
       sfn0 = sfn
       sfn = dsqrt(sf)
       xxx = dabs(sfn0 - sfn)
       if(xxx.eq.0.0d+0) xxx = 1d-15
!c
!c  - calcul du coefficient d'acceleration
!c
       if(nit.lt.2.or.nit/2*2.ne.nit) then
          acc = 1.0D+00
       else
          acc = sfn0/xxx
          do i=1,nco
             if(dabs(acc*f(i)).gt.1.0D+00) acc = 1.0D+00
          end do
          if(acc.gt.100.0d0) acc = 1.0D+00
       endif

       if(acc.ne.1.0D+00) then
          s = 0.0D+00
          do i=1,nco
             yy(i) = y1(i)*dexp(acc*f(i))
             s = s + yy(i)
          end do
       endif
       sm1 = s - 1.0D+00
!c
!c  - critere de passage a la methode de newton
!c
       if(fmaxi.le.1d-3.and.acc.eq.1.and.nit.gt.4) newt = 1
!c       if(nit.gt.100) newt = 1
       if(nit.gt.100 .OR. (p0>2.d0*pcmax .AND. nit>nitflash)) newt = 1
       goto 120
!c
!c
!c-------------------------------------------------
!c
!c  methode de newton
!c
!c
 130   nitn = nitn + 1
       itest = 0
       fmaxi = 0.0D+00
       testmax = 0.0D+00
       xxx = 1.0d0/s
       do i=1,nco
          f(i) = f0(i) - dlog(yy(i)) - lphi(i)
          if(dabs(f(i)).gt.fmaxi) fmaxi = dabs(f(i))
          test(i) = test0(i) + 10.0D+0*eps0 + dabs(dlphideta(i)*epseq)
          if(test(i).ge.testmax) testmax = test(i)
          if(dabs(f(i)).gt.test(i)) itest = itest + 1
          do k=1,nco
             d(k,i) = dmudn(k,i)*xxx + xxx
          end do
       end do
       isym = 1
       irtn = 0
       call chol(nco,nco,d,f,ttt,dn,isym,echec,irtn)

       if(echec.ne.0.or.nitn.eq.20) then ! echec ou non convergence de
          sm1 = 0.0D+00                  ! la methode de newton. xs est
          goto 150                       ! presumee non stable.
       endif

 140   s = 0.0D+00
       neg = 0
       do i=1,nco
          if(yy(i).gt.1.0D+300) yy(i)=1.0D+300
          yjn(i) = yy(i)                      ! calcul des nouveaux yy, avec
          yy(i) = yy(i) + dn(i)               ! verification de leur
          if(yy(i).le.0.0D+00) neg = neg + 1  ! caractere positif.
          s = s + yy(i)
       end do
       sm1 = s - 1.0D+00

       if(neg.ne.0) then                ! diminution des accroissements
          accmin = 1d+30                ! afin d'obtenir des yy tous
          do i=1,nco                    ! positifs.
             yy(i) = yjn(i)
             if(dn(i).ne.0.0d0) acc = - yy(i)/dn(i)
             if(acc.gt.0.0d0.and.acc.le.accmin) accmin = acc
          end do
          acc = 0.5D+00*accmin
          do i=1,nco
             dn(i) = acc*dn(i)
          end do
          goto 140
       endif

       if(itest.ne.0) goto 120 ! nouvelle iteration, sinon
!c                              ! fin de la methode de newton.
!c
!c
!c-------------------------------------------------
!c
!c
       if(nimp.ge.2) then
          if(nimp0.eq.0) then
             nimp0 = 1
             IF(FORT)write(lun1,9110)
             IF(FORT)write(lun2,9110)
          endif
          IF(FORT)write(lun1,9120) init,jd,etas(1),etas(2),fmaxi,sm1,testmax
          IF(FORT)write(lun2,9120) init,jd,etas(1),etas(2),fmaxi,sm1,testmax
       endif
 9110  format(' init   jd',6x,'eta0',6x,'etas',8x,'fmax' &
     & ,9x,'s-1',5x,'testmax     ptg')
 9120  format(2i5,2f10.5,3d12.2)
!c
!c
!c-------------------------------------------------
!c
!c  conclusions sur la stabilite du systeme de composition x0
!c
!c
 150   if(sm1.gt.testmax) then
          iph = 2             ! systeme diphasique
          if(ph3.eq.0) return

          num = num + 1
          etasn(num) = etas(2)
          gsn(num) = gs
          do i=1,nco
             xsn(i,num) = xs(i)
          end do
       endif

       if(init.eq.1) then     ! si on est parti d'une initialisation xs
          init = 0            ! on essaie une premiere init. de raoult
          goto 100
       else
          if(jd.eq.jd0) then
             jd = 3 - jd0     ! puis l'autre.
             js = jd
             goto 100
          else
             if(ph3.eq.0.or.num.eq.0) then
                iph = 1       ! systeme presume monophasique stable
                return
             else
                iph = 2
                gsmin = 1d+30
                do j=1,num
                   if(gsn(j).lt.gsmin) then
                      gsmin = gsn(j)
                      num0 = j
                   endif
                end do
                etas(2) = etasn(num0)
                do i=1,nco
                   xs(i) = xsn(i,num0)
                end do
                return
             endif
          endif
       endif

       end
!c
!c---------------------------------------------------------------------------
!c
      subroutine tpd_rp(z1,t,p00)
      use comflash
      implicit none

      real*8 eta0,v0,g0
      real*8 x0(ncomax),lnphi0(ncomax),lnphi(ncomax),x1,tpd,z1,x2,z2
      real*8 p00,t,d1,d2,eps,x1old,eps_tpd,stepx,n10,n1,n2
      real*8 tempo,xxmin,xxmax,eps_tpd2
      integer j0,ider,istb,iter,niter_max,compteur,i,npoints
      logical stable,actif
      common/is_stable/stable
      real*8 tpdminimum ! modif rp tempo
      common/mapoubelle/tpdminimum ! modif rp tempo
      common/debugging/actif ! modif rp tempo

      eps_tpd  = -1.d-10
      eps_tpd2 = -1.d-5
      stable   = .true.
      j0       = 1
      ider     = 0
      istb     = 1
      t0       = t
      x0(1)    = z1
      x0(2)    = 1.d0 - z1
      z2       = 1.d0 - z1
      niter_max = 500

      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
      lnphi0(1:ncomax) = lphi(1:ncomax)
      d1 = dlog(z1) + lnphi0(1)
      d2 = dlog(z2) + lnphi0(2)

!c 1ere etape : on fait un balayage en compo du TPD
!c ================================================

      xxmin = 1.d-9
      xxmax = 1.d0 - xxmin
      stepx = (xxmax - xxmin) / 150.d0
      compteur = 0
      npoints = 151
      do i = 1,npoints
         x1 = xxmin + dfloat(i - 1)*stepx
!      do x1 = xxmin,xxmax,stepx
         compteur = compteur + 1
         x0(1) = x1
         x0(2) = 1.d0 - x1
         x2 = 1.d0 - x1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
         lnphi(1:ncomax) = lphi(1:ncomax)
         tpd = x1 * (lnphi(1)-lnphi0(1) + dlog(x1/z1)) + &
     &         x2 * (lnphi(2)-lnphi0(2) + dlog(x2/z2))

         if(tpd.lt.eps_tpd) then
!c           Si la composition nonstable est environ egale a la composition
!c           de la phase testee, on se mefie (par exemple, si
!c           z1 et x1 sont distants de 0.01 et si tdp = -2.d-10, on est
!c           proche d'un min qui devrait valoir zero mais qui numeriquement
!c           peut etre legerement < 0)
            if(DABS(x1-z1).lt.5.d-3 .AND. tpd.gt.eps_tpd2) then
!c              on estime que le test de stabilite en ce point n'est pas
!c              concluant
               goto 400 ! on ne declare pas la compo non-stable et on
!c                         poursuit l'analyse du TPD
            endif
            stable = .false.
            goto 300 ! return
         end if
 400     continue
      enddo

      goto 300 ! EXIT / FIN DES TESTS **************************
!c      if(z1<0.999d0 .and. z1>0.001d0) goto 300 ! EXIT / FIN DES TESTS **************************

!c 2eme etape : Si stable, on cherche precisement a localiser les minima
!c =====================================================================

!c     a) Initialisation a x1 proche de 0
!c     ----------------------------------

      n1 = 1.d-12
      n2 = 1.d0 - n1
      x1 = n1 / (n1 + n2)
      x2 = 1.d0 - x1
      eps = 1.d0
      iter = 1

      do while(eps.gt.1.d-12)
         x0(1) = x1
         x0(2) = x2
         x1old = x1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

         n1 = dexp(d1-lphi(1))
         n2 = dexp(d2-lphi(2))

         if(n1.lt.1.d-13 .or. n2.lt.1.d-13) goto 297 ! probleme

         tempo = n1 + n2
         x1 = n1 / (n1 + n2)
         x2 = 1.d0 - x1

         eps = dabs(x1 - x1old)
         iter = iter + 1
         if(iter.gt.niter_max) goto 297
      enddo
 297  continue

      x0(1) = x1
      x0(2) = x2
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) + &
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

!c     b) Initialisation a x1 proche de 1
!c     ----------------------------------

      n1 = 0.999d0
      n2 = 1.d0 - n1
      x1 = n1 / (n1 + n2)
      x2 = 1.d0 - x1
      eps = 1.d0
      iter = 1

      do while(eps.gt.1.d-12)
         x0(1) = x1
         x0(2) = x2
         x1old = x1

         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

         n1 = dexp(d1-lphi(1))
         n2 = dexp(d2-lphi(2))

         if(n1.lt.1.d-13 .or. n2.lt.1.d-13) goto 298 ! probleme

         tempo = n1 + n2
         x1 = n1 / (n1 + n2)
         x2 = 1.d0 - x1

         eps = dabs(x1 - x1old)
         iter = iter + 1
         if(iter.gt.niter_max) goto 298
      enddo
 298  continue

      x0(1) = x1
      x0(2) = x2
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) + &
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

!c     c) Initialisation a x1 tres proche de 1
!c     ---------------------------------------

      n1 = 0.9999999d0
      n2 = 1.d0 - n1
      x1 = n1 / (n1 + n2)
      x2 = 1.d0 - x1
      eps = 1.d0
      iter = 1

      do while(eps.gt.1.d-12)
         x0(1) = x1
         x0(2) = x2
         x1old = x1

         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

         n1 = dexp(d1-lphi(1))
         n2 = dexp(d2-lphi(2))

         if(n1.lt.1.d-13 .or. n2.lt.1.d-13) goto 399 ! probleme

         tempo = n1 + n2
         x1 = n1 / (n1 + n2)
         x2 = 1.d0 - x1

         eps = dabs(x1 - x1old)
         iter = iter + 1
         if(iter.gt.niter_max) goto 399
      enddo
 399  continue

      x0(1) = x1
      x0(2) = x2
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) + &
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

      goto 300 ! EXIT / FIN DES TESTS **************************


!c     d) Initialisation de x1 entre xp(1,1) et xp(1,2)
!c     ------------------------------------------------

      n1 = 0.5d0 * (xp(1,1) + xp(1,2))
      n2 = 1.d0 - n1
      x1 = n1 / (n1 + n2)
      x2 = 1.d0 - x1
      eps = 1.d0
      iter = 1

      do while(eps.gt.1.d-12)
         x0(1) = x1
         x0(2) = x2
         x1old = x1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

         tempo = n1 + n2
         n1 = dexp(d1-lphi(1))
         n2 = dexp(d2-lphi(2))
         if(n1.lt.1.d-13 .or. n2.lt.1.d-13) goto 299 ! probleme

         x1 = n1 / (n1 + n2)
         x2 = 1.d0 - x1

         eps = dabs(x1 - x1old)
         iter = iter + 1
         if(iter.gt.niter_max) goto 299
      enddo
 299  continue

      x0(1) = x1
      x0(2) = x2
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) + &
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif


!c     e) Initialisation de x1 a diverses compo entre xp(1,1) et xp(1,2)
!c     -----------------------------------------------------------------

      stepx = 5.d-2
      n10 = 0.d0

 987  n10 = n10 + stepx

      n1 = n10
      if(n1.ge.1.d0) goto 300

      n2 = 1.d0 - n1
      x1 = n1 / (n1 + n2)
      x2 = 1.d0 - x1
      eps = 1.d0
      iter = 1

      do while(eps.gt.1.d-12)
         x0(1) = x1
         x0(2) = x2
         x1old = x1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

         tempo = n1 + n2
         n1 = dexp(d1-lphi(1))
         n2 = dexp(d2-lphi(2))

         if(n1.lt.1.d-13 .or. n2.lt.1.d-13) goto 301 ! probleme

         x1 = n1 / (n1 + n2)
         x2 = 1.d0 - x1

         eps = dabs(x1 - x1old)
         iter = iter + 1
         if(iter.gt.niter_max) goto 301
      enddo
 301  continue

      x0(1) = x1
      x0(2) = x2
      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) + &
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

      goto 987

 300  continue
      tpdminimum = tpd ! modif rp tempo

      end
!c
!c---------------------------------------------------------------------------
!c
      subroutine search_min_tpd(z1,t,p00,xmin)

!c     Recherche des points stationnaires du TPD
!c     Si TPDmin = 0, on cherche la solution qui ne correspond pas a z1
      use comflash
      implicit none

      real*8 eta0,v0,g0,xstat_up(20),xstat_dn(20),xstat(20)
      real*8 x0(ncomax),lnphi0(ncomax),x1,tpd,z1,x2,z2
      real*8 p00,t,d1,d2,x1old,eps_tpd,stepx
      real*8 tempo,xmin,xmax,dtpd,tpd1,tpd2,dtpdold,tpdold
      real*8 tpdstat(20),tpdmin
      integer j0,ider,istb,iter,niter_max,i,k,nstat,npoints
      logical stable,croissant(20)
      common/is_stable/stable

      xstat  = 0.d0
      eps_tpd= -1.d-12
      stable = .true.
      j0     = 1
      t0     = t
      ider   = 0
      istb   = 1
      x0(1)  = z1
      x0(2)  = 1.d0 - z1
      z2     = 1.d0 - z1
      niter_max = 1000

      call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
      lnphi0(1:ncomax) = lphi(1:ncomax)
      d1 = dlog(z1) + lnphi0(1)
      d2 = dlog(z2) + lnphi0(2)

      xmin = 1.d-12
      xmax = 1.d0 - xmin
      stepx = (xmax - xmin) / 5000.d0
      dtpd = 0.d0
      tpd = 0.d0
      k = 1

!c     Resolution de l'equation (dTPD / dx1) = 0 (points stationnaires de TPD)
!c     xstat_dn(i) est une approximation par valeurs inf de la i-eme solution
!c     xstat_up(i) est une approximation par valeurs sup de la i-eme solution
      npoints = 5001
      do i = 1,npoints
!      do x1 = xmin,xmax,stepx
         x1 = xmin + dfloat(i-1) * stepx
         x2 = 1.d0 - x1
         x0(1)  = x1
         x0(2)  = 1.d0 - x1
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
         tpd1 = dlog(x1) + lphi(1) - d1
         tpd2 = dlog(x2) + lphi(2) - d2
         dtpdold = dtpd
         tpdold = tpd
         dtpd = tpd1 - tpd2
         tpd  = x1 * tpd1 + x2 * tpd2

         if(x1.gt.xmin) then
            tempo = dtpd * dtpdold
            if(tempo.lt.0.d0) then ! tangente horizontale
               xstat_up(k) = x1
               xstat_dn(k) = x1 - stepx
               if(dtpd.gt.0.d0) then
                  croissant(k)= .true.
               else
                  croissant(k)= .false.
               endif
               k = k + 1
               if(k.gt.21) then
                  print*,'Probleme search_min_tpd'
                  print*,'Modifier dimension de XSTAT'
                  stop
               end if
            endif
         endif
      enddo

      nstat = k - 1 ! nombre de points stationnaires

!c     Une fois les points stationnaires approximes, on les calcule precisement
!c     par dichotomie

      do i = 1,nstat
         xmin = xstat_dn(i)
         xmax = xstat_up(i)

         iter = 0
 400     continue
         x1old = x1
         x1 = 0.5d0 * (xmin + xmax)
         x2 = 1.d0 - x1
         x0(1)  = x1
         x0(2)  = x2
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
         tpd1 = dlog(x1) + lphi(1) - d1
         tpd2 = dlog(x2) + lphi(2) - d2
         tpd  = x1 * tpd1 + x2 * tpd2
         dtpd  = tpd1 - tpd2
         if(dtpd.gt.0.d0) then
            if(croissant(i)) then
               xmax = x1
            else
               xmin = x1
            endif
         else
            if(croissant(i)) then
               xmin = x1
            else
               xmax = x1
            endif
         end if
         iter = iter + 1

         if(dabs(x1old-x1).gt.1.d-12) then
            goto 400
         endif
         xstat(i)   = x1
         tpdstat(i) = tpd
      enddo

!c     A present, on calcule les min
      tpdmin  = 1.d31

      do i = 1,nstat
!c        On prend la solution autre que z1
         if(dabs(xstat(i)-z1).lt.1.d-5) goto 401

         if(tpdstat(i).lt.tpdmin) then
            tpdmin = tpdstat(i)
            xmin   = xstat(i)
         end if

 401     continue
      enddo

      if(tpdmin .gt.1.d-3) then
         xmin = z1
      endif
! 420  continue

      end
!c
!c---------------------------------------------------------------------------
!c
     subroutine prtfl
!c
      use comflash
      implicit none

!c*** variables internes
!c
      integer i,j,iphase
      real*8 ki(phmax)
      real*8 gorst
	LOGICAL FORT
!c-------------------------------------------------
	FORT=.TRUE. 	! ENABLES DEBUGGING MESSAGES WRITTEN TO CONSOLE JRE 20210517. "FORT" = "LOUD" IN FRENCH.
	FORT=.FALSE.
!c
!c
!c-------------------------------------------------
!c
!c
	if(iph.ne.0) then
	  iphase = iph
	else
	  iphase = 1
	endif
	do i=1,nco
	  do j=1,iphase-1
		 ki(j) = xp(i,iphase)/xp(i,j)
	  end do
	ENDDO
    if(iphase.eq.2)gorst = 83.14411d+00*288.15d+00*npt(2)/(1.01325d+00*vp(1))
	IF(.NOT.FORT)RETURN
	IF(corrFL) then
	  write(lun1,*)
	  write(lun2,*)
	  write(lun1,*) 'Flash RP'
	  write(lun2,*) 'Flash RP'
	ELSE
	  write(lun1,*)
	  write(lun2,*)
	  write(lun1,*) 'Flash Peneloux'
	  write(lun2,*) 'Flash Peneloux'
	ENDIF

   if(iph.eq.0) then
      write(lun1,9000) t0,p0
      write(lun2,9000) t0,p0
   elseif(iph.eq.1) then
      write(lun1,9002) t0,p0
      write(lun2,9002) t0,p0
   elseif(iph.eq.2) then
      write(lun1,9004) t0,p0
      write(lun2,9004) t0,p0
   endif
 9000  format(1x,'T =',f15.10,' K,   P =',f18.10,' bar,' &
     & ,3x,'monophasique (incertain)')
 9002  format(1x,'T =',f15.10,' K,   P =',f18.10,' bar,' &
     & ,3x,'monophasique')
 9004  format(1x,'T =',f15.10,' K,   P =',f18.10,' bar,' &
     & ,3x,'diphasique')

       if(iph.le.1) then
          write(lun1,9030)
          write(lun2,9030)
       else
          if(iph.eq.2) then
             write(lun1,9040)
             write(lun2,9040)
          else
             write(lun1,9050)
             write(lun2,9050)
          endif
       endif
 9030  format(1x,'constituant',8x,'moles',5x,'fr. mol.')
 9040  format(1x,'constituant',6x,'moles',9x,'fractions molaires' &
     & ,7x,'k',/,18x,'global',7x,'liquide',6x,'vapeur')
 9050  format(1x,'constituant',8x,'moles',8x,'fractions molaires',15x &
     & ,'k31',8x,'k32',/,18x,'global  liquide 1  liquide 2',5x,'vapeur')

       do i=1,nco
          write(lun1,9080) nom(i),nt(i),(xp(i,j),j=1,iphase) &
     &                     ,(ki(j),j=1,iphase-1)
          write(lun2,9080) nom(i),nt(i),(xp(i,j),j=1,iphase) &
     &                     ,(ki(j),j=1,iphase-1)
 9080     format(1x,a10,3x,f11.5,2x,f11.5,2x,f11.5,3x,g13.6)
       end do

       write(lun1,9100) ntt,(npt(j),j=1,iphase)
       write(lun2,9100) ntt,(npt(j),j=1,iphase)
       write(lun1,9120) mmolt,(mmolp(j),j=1,iphase)
       write(lun2,9120) mmolt,(mmolp(j),j=1,iphase)
       if(iphase.eq.2) then
          gorst = 83.14411d+00*288.15d+00*npt(2)/(1.01325d+00*vp(1))
          write(lun1,9141) vt,(vp(j),j=1,iphase),dropout
          write(lun2,9141) vt,(vp(j),j=1,iphase),dropout
          write(lun1,9142) vt/ntt,((vp(j)/npt(j)),j=1,iphase)
          write(lun2,9142) vt/ntt,((vp(j)/npt(j)),j=1,iphase)
       else
          write(lun1,9140) vt,(vp(j),j=1,iphase)
          write(lun2,9140) vt,(vp(j),j=1,iphase)
          write(lun1,9143) vt/ntt,((vp(j)/npt(j)),j=1,iphase)
          write(lun2,9143) vt/ntt,((vp(j)/npt(j)),j=1,iphase)
       endif
       write(lun1,9160) denst,(densp(j),j=1,iphase)
       write(lun2,9160) denst,(densp(j),j=1,iphase)
       write(lun1,9180) zt,(zp(j),j=1,iphase)
       write(lun2,9180) zt,(zp(j),j=1,iphase)
       if(iphase.eq.2) then
          write(lun1,9182) gorst
          write(lun2,9182) gorst
       endif
 9100  format(/,1x,'moles        ',f11.5,2x,f11.5,2X,f11.5)
 9120  format(1x,'M g/mol      ',f11.2,2x,f11.2,2x,f11.2)
 9140  format(1x,'Vol_tot/cm3  ',3x,g14.7,1x,g14.7)
 9141  format(1x,'Vol_tot/cm3  ',3x,g14.7,1x,g14.7,1x,g14.7, &
     &        2x,'drop out :',f6.2,' %')
 9142  format(1x,'vol_mol/cm3  ',3x,g14.7,1x,g14.7,1x,g14.7)
 9143  format(1x,'vol_mol/cm3  ',3x,g14.7,1x,g14.7)
 9160  format(1x,'densite g/cm3',f11.5,2x,f11.5,2X,f11.5)
 9180  format(1x,'z            ',f11.5,2x,f11.5,2X,f11.5)
 9182  format(1x,'g.o.r.       ',f11.2,' sm3/m3')
       return
       end
       subroutine longueur (nom,long)
       implicit none
!c
!c  sous-programme de recherche de la longueur d'un mot
!c
       integer i,long
       character*80 nom
       character*3 blanc
!c
       data blanc/'   '/
!c
       long = 0
       do i=1,80
         if(nom(i:i+4).eq.blanc) goto 10
         long = long + 1
       end do
!c
 10    return
       end
!c*****************************************************************************
!	Added for PGL6ed
!c*****************************************************************************
!c Notes for general fugacity calls:
!c       icon  level of derivatives
!c       1     p and fg only
!c       2     also t- and v-derivatives
!c       3     also n-derivativs
!c       4     all derivatives
!c     ntyp:     (i):      phase type desired
!c                       1 = liq, -1 = vap, 0 = min g phase
!c Notes for fugacity_TP calls:
!c     mtyp:     (o):      indicator for phase type; ic returns
!c                       1:   if called with mt = 1/-1 and phase is found
!c                       -1:  if called with mt = 1/-1 and phase is not found
!c                       2:   if called with mt = 0 and min g phase is liquid
!c                       -2:  if called with mt = 0 and min g phase is vapour
!c  JRE20210513: Restricting GlobConst with "only" avoids conflicts with other GlobConst variables (e.g. Rgas).
!c               FYI, SetNewEos is a function "contain"ed in GlobConst.
Subroutine GetPrLorraine(nComps,iErr)
	USE GlobConst, only:SetNewEos,PGLInputDir,LOUD,bVolCc_mol,zeroTol,etaMax,ID,iEosOpt,tKmin
	!USE PrTcParms, only: ! so we can call GetPrTc and set tKmin
	USE comflash
    integer :: i,fluids(ncomax),err0,ierr
	LOGICAL LOUDER
	LOUDER=LOUD
	LOUDER=.TRUE.
	!LOUDER=.FALSE.
	iErr=SetNewEos(iEosOpt) ! returns 0. Wipes out previous possible declarations of isTPT or isPcSaft.
	call GetPrTc(nComps,iErrGet) ! sets tKmin, as well as some other parameters that PrLorraine doesn't really use. 
	tmemc = 0.0d0	!Array that stores old temperatures in memory in order to	save some time
	itempc = 0		!Both of these need to be initialized to wipe out values stored for previous compounds.
	nco=nComps		! copies nComps from GlobConst to comflash
	if(nco > ncomax)then
		iErr=1
		if(LOUDER)write(*,*)'GetPrLorraine: nCoMax, nComps =',nCoMax,nComps
		return
	endif
    call setup(PGLinputDir,err0) ! COMPUTES OMEGAa, OMEGAb, Zc, 
	!FYI: setup calls READ_OPTION() which hard codes choix_ge=2, meaning the "HV-uniquac/wilson" model. 
	if(err0 > 0)then
		iErr=2
		if(LOUDER)write(*,*)'GetPrLorraine: Setup(), iErr=',err0
		return
	endif
	fluids(1:2)=ID(1:2)
    call setfluid(PGLinputDir,fluids,err0)
	if(err0 > 0)then
		iErr=3
		if(LOUDER)write(*,*)'GetPrLorraine: SetFluid(), iErr=',err0
		return
	endif
	etaMax=1-zeroTol
	do i=1,nComps
		bVolCc_mol(i)=b(i) !found in comflash ~line 109. Divide by 10 because Lorraine assumes Rgas=83.1447.
	enddo
	if(LOUDER)then
        write(*,*)'Pure-compound properties'
        write(*,105)'Component','Tc (K)','Pc (bar)','L','M','N','b (cm3/mol)','c (cm3/mol)'
        do i = 1,nco
            write(*,106)trim( idComp(i,1) ),tc(i),pc(i),parL(i),parM(i),parN(i),b(i),c(i)
        enddo
        !BIP
        write(*,*)
        write(*,'(A)')' BIP (K)'
        do i = 1,nco
            write(*,'(1x,A15,A,I1,A,10(X,F12.4))')trim( adjustl(idComp(i,1)) ), ' (',i,')',agE0(i,1:nco)
        enddo
        write(*,*)
	endif
105 FORMAT(T2,A,T28,A,T41,A,T54,A,T65,A,T74,A,T84,A,T96,A)
106 FORMAT(T2,A,T25,F10.4,T40,F8.4,T52,F8.4,T63,F8.4,T72,F8.4,T84,F8.3,T96,F8.3)
	return
end subroutine !GetPrLorraine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine SetParMixPrLorraine(iParm,value,iErr) ! called from UaEosTools.
	USE comflash      !Just for PrLorraine.
	DoublePrecision value
	iErr=0
	if(iParm==1)then
				 age0(1,2)=value
	elseif(iParm==2)then
				 age0(2,1)=value
	else
				 iErr=1
	endif
	return
end	subroutine !SetParMixPrLorraine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine QueryParMixPrLorraine(iParm,value,iErr) ! 
	USE comflash      !Just for PrLorraine.
	DoublePrecision value
	iErr=0
	if(iParm==1)then
				 value=age0(1,2)
	elseif(iParm==2)then
				 value=age0(2,1)
	else
				 iErr=1
	endif
	return
end	subroutine !SetParMixPrLorraine
