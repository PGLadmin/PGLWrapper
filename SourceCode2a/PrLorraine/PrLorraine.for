
c
c**********************************************************************
c
      program EoS_AE
c      -- version DETILDEE juillet 2018
c      -- ajout code ELLV novembre 2018
c
c     ce programme ajuste les parametres des regles de melange d'une
c     EoS cubique avec regles de melange de type AE (a d'exces)
c
c     NB : PR2SRK non programme dans cette version (A FAIRE)
c
      use load_bdd_mod
      implicit none
c
      include '..\source\comflash.cmn'
c
      integer ncpmax,nmax,itermax,nseriemax
      parameter (ncpmax=2*nsysmax+1,nmax=10,itermax=3000)
      parameter(nseriemax=30) ! pour azeo et LLV
c
c nsysmax = nbre max de systemes binaires dans le fichier de donnees (common)
c nsermax = nbre de series de points maximal pour un binaire donne (common)
c nptmax  = nbre max de points dans une serie donnee (pour un binaire donne) en common
c ncpmax  = nbre max de corps purs differents dans le fichier de donnees
c npmax  = nbre maxi de points calcules sur une courbe isot ou isop (common)
c nmax = nbre maxi de parametres a ajuster
c
      character*80 name,nomc1(ncpmax),nom_ok(ncpmax),directory
      character*80 nomy,mapres,texte
      character*800 phrase,phras,phrase0
      character*40 fichres
      character*12 nomfich
      character*3 nam(999)
      character*15 nomcomp
      character*60 tit
      character*1  blanc
      character*8 date_in,date_out
      character*1022 dexcri,dexfil,dex_cri,dex_fil
      character*180 creation

      integer i,ii,jj,nexp,npoin(nsysmax),j,nserie,k,nsyst
      integer kn,kk,code,ndata,it,ip,ndat(nsysmax),rien,codc1(ncpmax)
      integer cod_ok(ncpmax),compteur,itr(nsysmax),itrh(nsysmax)
      integer itrcp(nsysmax),i2,nsertot,ly,rr,rp
      integer ipr(nsysmax),er1,er2,npsat,nteb
      integer l1,ll1,l2,l3,jno7,lun,answer,m_bfgs
      integer nbg,ncod,cod(pmax),ij,ijj,ik
      integer npt_x,npt_y,npt_c,npt_tot,npt_h,npt_cp,npt_az,n_az
      integer npt_LLV,n_LLV
      integer passage,ajout,ldir,jt,jp,nl_titre,n_serie
      integer n_bu,n_ro,n_cr,n_cp,n_h
      integer n_seriet,n_seriep,compteur_serie,compt_dest,compt_desp
      integer compt_desh,compt_descp,n_seriecp,n_serieh
      integer ntemper,ntbltemp,npressi,ntblpres,num(nsermax)
      integer ntbltemph,ntbltempcp
      integer codes(nsermax),col_bulle,col_rosee,npoint(nsermax)
      integer n_param,ecrit,filtre,des_filtre,des_cri
      integer time_in,time_out,mois_in,jour_in,mois_out,jour_out
      integer njour,n_jour,n_hour,n_min,n_sec,duree_s
      integer iter,reste,imini,ia1_2,ib1_2,ia1_3,ib1_3,ia1_4,ib1_4
      integer n_passage_1,n_passage_2,n_passage_3,n_passage_4
      integer n_passage_5,n_passage_6,n_hand
      integer n_courbe,jt_init,jp_init,begin,err0
      integer group_nb(50),lcri,lfil,l_cri,l_fil,start,finis
      integer rpaz,nser_az(nseriemax),rpLLV,nser_LLV(nseriemax)
      integer rpcri,nser_cr(nseriemax) !modif tempo APM 31/10/19     
      integer nptsh,nptscp,h,nser_non_tracees,truc,ans
      integer nbphases
      integer end_courbe(5)
      integer n_couple,nb_copies,nmesh1,nmesh2
      parameter(nptsh=399,nptscp=99)

      real*8 mmc1(ncpmax),tcc1(ncpmax),pcc1(ncpmax),omegac1(ncpmax)
      real*8 mm_ok(ncpmax),tc_ok(ncpmax),pc_ok(ncpmax),om_ok(ncpmax)
      real*8 tde(nsermax),tabltemp(nsermax),pde(nsermax),TildAB
      real*8 tablpres(nsermax),oldpres(200)
      real*8 tabltemph(nsermax),tabltempcp(nsermax)
      real*8 huit,temp(nsermax),p(nsermax,nptmax)
      real*8 pres(nsermax),temptr(nsysmax,nsermax)
      real*8 prespr(nsysmax,nsermax),temptrcp(nsysmax,nsermax)
      real*8 temptrh(nsysmax,nsermax)
      real*8 t(nsermax,nptmax),diff
      real*8 xl(nsermax,nptmax),xg(nsermax,nptmax),pression
      real*8 xl2(nsermax,nptmax)
      real*8 temperature,factp(nsermax)
      real*8 factt1(nsermax),factt2(nsermax)
      real*8 sk(nbgr_tot_max),tdet,pdet
      real*8 xb(npmax),xr(npmax),pb(npmax),pr(npmax),p_min,p_max
      real*8 tb(npmax),tr(npmax),t_min,t_max
      real*8 xmin,xmax,pasx,ymin,ymax,pasy,yymax
      real*8 valeurx(50),valeury(50)
      real*8 theta(nmax),fobj
      real*8 affiche,fo(itermax),step_a,step_b,ecart_min
      real*8 init_1,init_2,init_3,init_4,init_5,init_6
      real*8 agr1_2(itermax),bgr1_2(itermax)
      real*8 agr1_3(itermax),bgr1_3(itermax)
      real*8 agr1_4(itermax),bgr1_4(itermax)
      real*8 psat0,psat1,psat2,psatura(2),xsatura(2)
      real*8 teb0,teb1,teb2,tebull(2),xebull(2)
      real*8 p_llv,p_ll,t_ll,t_llv,xI_llv,xII_llv,y1_llv,corr
      real*8 poub,pazeo(nseriemax),tazeo(nseriemax),xazeo(nseriemax)
      real*8 pLLV(nseriemax),tLLV(nseriemax)
      real*8 xLLV(nseriemax),xLLV2(nseriemax),xLLV3(nseriemax)
      real*8 prcri(nseriemax),tmcri(nseriemax),xxcri(nseriemax)
      real*8 h_min,h_max,hMexp(nsermax,nptmax),zexp(nsermax,nptmax)
      real*8 cp_min,cp_max,cpMexp(nsermax,nptmax),xcomp(2)
      real*8 presh(nsermax),patm,hMcalc(nsermax,nptsh+2,2)
      real*8 prescp(nsermax),cpMcalc(nsermax,nptscp+2,2)
      real*8 hMel,cpMel,Cpec,x1L,x1G,hmL,hmG,tempo,x1gaz,hmgaz
      real*8 fobj_opt,theta_opt(nmax),gg(nmax),gg_opt,tempo_const
      real*8 x,nflash(pmax),p00,aamin,aamax,bbmin,bbmax
      real*8 cppur(pmax),cpgaz(pmax),val_a
      
      
      parameter(patm=101325.d-5)
      logical posit,onephase(nptsh+2)
      logical het_gauche,het_droite,llv_gauche,llv_droite
      logical is_ELL0,ligne_triph
      logical azeo_a_prendre(nseriemax),exist_Paz
      logical LLV_a_prendre(nseriemax),exist_PLLV
      logical cri_a_prendre(nseriemax),exist_Pcri
      logical is_EGG,change_echelle,onlyFobj,onlyGraph
      logical tilded,maillage,maillage_PRWil

      common/eq_liqliq/het_gauche,het_droite,llv_gauche,llv_droite,
     &                 ligne_triph,p_llv,p_ll,xI_llv,xII_llv,y1_llv,
     &                 is_ELL0,t_ll,t_llv
      common/eq_gasgas/is_EGG

      data blanc/' '/
      option_debug = .false. 
c
c    Lecture de la BDD 
      call load_bdd()
      OPEN(10,FILE='verif_eos.txt')
      OPEN(220,file='hM_cPM.txt')
c    le programme est ici limite a 999 series de points differentes isot/isop
c    (nam(i) definit le nombre de dessins isot ou isop)
c
      call timer(time_in)  !temps en centisecondes ecoule depuis minuit
      time_in = int(time_in/100.0d+00) !temps en secondes ecoule depuis minuit
      call date(date_in) !date sous la forme "MM/JJ/AA" (mois/jour/annee)
      open(3000,status='scratch')
      write(3000,*) date_in(1:2) !le mois
      write(3000,*) date_in(4:5) !le jour
      rewind(3000)
      read(3000,*) mois_in
      read(3000,*) jour_in
      close(3000)

c     Variable nam utilisee pour generer noms de fichiers de dessin 
      do i = 1,999
         write(nam(i),'(i3)') i
      enddo
      nam = ADJUSTL(nam)
      do i=1,1022
         dexcri(i:i) = blanc
         dexfil(i:i) = blanc
         dex_cri(i:i) = blanc
         dex_fil(i:i) = blanc
      end do
      fl_ok = .true.  !initialisation
      ygtx = .true.   !initialisation
      isGPED = .false. !initialisation
      nombre_flash = 0
      appel_FO = 0
      choix_eos = -1  !initialisation
      choix_ge = -1  !initialisation
      p_max_exp = 0.0d0
      T_min_exp = 1.0d99

c
c     Paramètres pour les règles de mélange :
      choix_b = -1  !initialisation
      quadra_b = .FALSE. !initialisation
      Lij(1:nco,1:nco) = 0.d0 !initialisation
      agE0(:,:) = 0.d0
      bgE0(:,:) = 0.d0
      cgE0(:,:) = 0.d0

c     Conversion Akl / Bkl tildés -> Akl / Bkl détildés
      TildAB = 1.5d0 + DSQRT(2.d0)
      
      nom(:) = ''
      tc(:) = 0.d0
      pc(:) = 0.d0
      fac(:) = 0.d0
      zra(:) = 0.d0
      qi(:) = 0.d0
c
c     ouverture et lecture du fichier option.d qui va piloter le programme
c
      open(541,file='option.d',status='old')

      do i = 1,6
         read(541,*) texte ! commentaires
      enddo
      
c     DEFINITION DU MODELE      
      read(541,*) texte,choix_eos
      if(choix_eos.gt.4 .or. choix_eos.lt.1) THEN
         write(lun1,'(1X,A)') 'Choix EoS a revoir (OPTION.D)'
         write(lun2,'(1X,A)') 'Choix EoS a revoir (OPTION.D)'
         stop
      endif
      IF(choix_eos.eq.1) THEN ! PPR78
         modele_gE = "PPR78 (Van Laar)" ! initialisation
      ELSEIF(choix_eos.le.3) THEN ! PR-Wilson/UNIQUAC
         if(choix_eos == 2) choix_ge = 1 !Wilson
         if(choix_eos == 3) choix_ge = 2 !UNIQUAC
         choix_eos = 2
         modele_gE = "HV (Uniquac-Wilson)" ! initialisation
      ELSEIF(choix_eos.eq.4) THEN ! PR-NRTL
         choix_ge = 3
         modele_gE = "HV (NRTL)" ! initialisation
c         alphgE = 0.6d0
      ENDIF
      read(541,*)
      read(541,*)
      read(541,*) texte,FoncAlpha
      IF(FoncAlpha.NE.'SOAVE' .AND. FoncAlpha.NE.'TWU91'
     &    .AND. FoncAlpha.NE.'TWU88') THEN
         WRITE(lun1,'(1X,A)') "Fonction alpha non geree (OPTION.D)"
         WRITE(lun2,'(1X,A)') "Fonction alpha non geree (OPTION.D)"
         stop
      ENDIF
      read(541,*)
        
      read(541,*) texte,choix_b
      read(541,*) 
      read(541,*) 
      if(choix_b.gt.2 .or. choix_b.lt.1) THEN
         write(lun1,'(1X,A)') 'Choix b a revoir (OPTION.D)'
         write(lun2,'(1X,A)') 'Choix b a revoir (OPTION.D)'
         stop
      endif
      if(choix_b.eq.1) then
         quadra_b = .FALSE. ! RM sur b lineaire
         puisb = 1.d0
         read(541,*) texte ! L12
         read(541,*) texte ! puisb
      else
         quadra_b = .TRUE. ! RM sur b quadratique
         read(541,*) texte,Lij(1,2)
         read(541,*) texte,puisb
         Lij(2,1) = Lij(1,2)
      endif
      read(541,*)
      read(541,*)
      if(choix_eos == 1) then ! PPR78
         read(541,*) texte,tabingr_file
      else
         tabingr_file = ""
         nbgr_tot = nbgr_tot_max
         read(541,*)
      endif

      if(choix_eos.ge.2) then ! PR-WILSON/UNIQUAC/NRTL
         read(541,*) texte,Lambda
      else
         Lambda = 0.d0
         read(541,*)
      endif

c     PARAMETRES POUR LE CALCUL DES EQUILIBRES
c     ET LE TRACE DES GRAPHES
      do i = 1,4
         read(541,*) ! lignes de commentaires
      enddo
      read(541,*) texte,answer  !(type de flash : simple ou renforce)
      if(answer.eq.1 .or. answer.eq.2) then
         choix_flash = answer
      else
         write(lun1,*) 'Choix incorrect du type de flash (OPTION.D)'
         write(lun2,*) 'Choix incorrect du type de flash (OPTION.D)'
         stop
      endif
      read(541,*)  !ligne de commentaire
      read(541,*) texte,answer  !(variable voir_ecart)
      voir_ecart = .false.
      if(answer.eq.1) voir_ecart = .true.
      if(voir_ecart) then
         read(541,*) texte,ecmax0
      else
         ecmax0 = 1.0d31
         read(541,*)
      endif
      read(541,*)  !ligne de commentaire
      read(541,*)  !"Affichage en cours de calcul/ajustement ... "
      read(541,*) texte,answer
      voir_bug = .false.
      onlyFobj = .false.
      onlyGraph = .false.      
      if(answer.eq.1 .or. answer.eq.2 .or. answer.eq.3
     &   .or. answer.eq.4)voir_bug = .true.
      if(answer.eq.2) option_debug = .true.
      if(answer.eq.3) onlyFobj = .true.
      if(answer.eq.4) onlyGraph = .true.
      read(541,*)  !ligne de commentaire
      read(541,*) texte,yrejet
      read(541,*) texte,ecmax_y
      read(541,*) texte,xrejet
      read(541,*) texte,ecmax_x
      read(541,*)  !ligne de commentaire
      read(541,*)  !'gestion des diagrammes'
      read(541,*)  !'isotherme'
      read(541,*) texte,pmaxim
      read(541,*) texte,pminim
      read(541,*)  !'isobare'
      read(541,*) texte,tmaxim
      read(541,*) texte,tminim
      read(541,*) texte,dxmax_global
      read(541,*)  !ligne de commentaire

      pmax_opt = pmaxim     !on sauve les valeurs de option.d
      pmin_opt = pminim
      tmax_opt = tmaxim
      tmin_opt = tminim

      do i = 1,10
         read(541,*)  !commentaire
      enddo
c
c     MODELE E-PPR78
c     **************
c
c     Calcul avec valeur constante de k12 (kij_const=TRUE)
c     Calcul predictif avec le modele PPR78 (predictif=TRUE)
c     Ajustement automatique des parametres du modele PPR78 (fitting=TRUE)
c     SUPPRIME DE CETTE VERSION : ajustement manuel des parametres du modele PPR78 (fit_hand=TRUE)
c
      kij_const = .false.
      predictif = .false.
      fitting = .false.
      simul_PRWil = .false.
      fitting_PRWil = .false.
	  tilded = .true.
	  maillage = .false.  
      maillage_PRWil = .false.	  
c     "   (1.1) SIMULATION - k12 = constante"       
      read(541,*) texte,answer  !(variable kij_const)
      if(answer.eq.1) kij_const = .true.
      read(541,*) texte,k12_const
      if(.not.kij_const) k12_const = 8888.d0

c     "   (1.2) SIMULATION PREDICTIVE E-PPR78 "
      g1 = 0
      g2 = 0
      g3 = 0
      g4 = 0
      read(541,*) texte,answer  !(variable predictif)
c      if(answer.eq.1 .and. .not.kij_const) then
      if((answer.eq.1 .or.answer.eq.3) .and. .not.kij_const) then   !modif tempo APM
c         ans = answer
         if (answer.eq.3) tilded = .false.
         predictif = .true.
         read(541,*) texte,n_param
         read(541,*) !Si 2 param, definir les groupes g1 et g2
         read(541,*) !Si 4 param, definir les groupes g1, g2, g3
         read(541,*) !Si 6 param, definir les groupes g1, g2, g3, g4
         if(predictif .AND. n_param.eq.2) then
            read(541,*) texte,g1,g2
            read(541,*) 
         elseif(predictif .AND. n_param.eq.4) then
            read(541,*) texte,g1,g2,g3
            read(541,*) 
         elseif(predictif .AND. n_param.eq.6) then
            read(541,*) texte,g1,g2,g3,g4
            read(541,*) 
         else
            read(541,*) 
            read(541,*) 
         endif
      else
         do i = 1,6
            read(541,*) 
         enddo
      endif
         
c      "   (1.3) AJUSTEMENT AUTOMATIQUE E-PPR78 (BFGS) "
      read(541,*) texte,answer  !(variable fitting)
c      if(answer.eq.1 .and. .not.kij_const.and..not.predictif) then
      if((answer.eq.1 .or. answer.eq.3) .and. .not.kij_const
     &                                   .and..not.predictif) then    !modif tempo APM
         if (answer.eq.3) maillage = .true.
         fitting = .true.
         read(541,*) texte,n_param
         m_bfgs = n_param
         if(fitting .AND. n_param.NE.2 .AND. n_param.NE.4
     &      .AND. n_param.NE.6) THEN
            write(lun1,*) 'Ajustement auto E-PPR78'
            write(lun2,*) 'Ajustement auto E-PPR78'
            write(lun1,*) 'PROBLEME n_param incorrect - STOP'
            write(lun2,*) 'PROBLEME n_param incorrect - STOP'
            stop
         endif
         if(n_param.eq.2) then
            read(541,*) texte,g1,g2
            read(541,*) 
         elseif(n_param.eq.4) then
            read(541,*) texte,g1,g2,g3
            read(541,*) 
         elseif(n_param.eq.6) then
            read(541,*) texte,g1,g2,g3,g4
            read(541,*) 
         else
            read(541,*) 
            read(541,*) 
         endif
      else
         do i = 1,3
            read(541,*) 
         enddo
      endif
c      "   (1.4) AJUSTEMENT AUTOMATIQUE PAR SYSTEME BINARE (E-)PPR78 (BFGS) "
      read(541,*) texte,answer  !(variable fitting)
      if(answer.eq.1 .and. .not.kij_const.and..not.predictif
     &                                    .and..not.fitting) then
c        fitting = .true.
         borneinf(:) = 8888.d0
         bornesup(:) = 8888.d0
         fitbin  = .true.
         n_param = 2
         m_bfgs = n_param
         g1 = 1
         g2 = 2
c    No des couples de valeurs intiales aleatoires (A_kl,B_kl)
         read(541,*) texte,n_couple
         read(541,*) texte,borneinf(1),texte,bornesup(1)
         read(541,*) texte,borneinf(2),texte,bornesup(2)
         do i = 1,6
            if(borneinf(i)<8887.d0 .or. borneinf(i)>8889.d0) then
c              on dé-tilde :
               borneinf(i) = borneinf(i)/TildAB
            endif

            if(bornesup(i)<8887.d0 .or. bornesup(i)>8889.d0) then
c              on dé-tilde :
               bornesup(i) = bornesup(i)/TildAB
            endif
         enddo         
      else
         do i = 1,3
            read(541,*) 
         enddo
      endif
      
      do i = 1,7
         read(541,*) 
      enddo
c     lecture du fichier PARFLASH qui contient des parametres
c     de l'EOS
c     Si choix_eos = 1, lecture des agr (Akl) et bgr (Bkl) du modele PPR78
c     Cet appel est fait avant la lecture des agr et bgr du fichier
c     d'option de sorte que les params du fichier d'option
c     ont le dernier mot. 
      call lecflash()

      if(fitting .OR. predictif) then
c       if(fitting) then
c        FITTING : le vecteur theta(:) contient :
c        - tous les parametres Akl a ajuster
c        - puis tous les parametres Bkl a ajuster
c        PREDICTIF : les coefficients agr et bgr doivent
c        etre renseignes ; les theta sont inutiles
         if(n_param == 2) then
            read(541,*) texte,theta(1)
            read(541,*) texte,theta(2)
c           on dé-tilde les Akl et Bkl
            if(tilded)then
               theta(1:2) = theta(1:2)/TildAB !modif tempo APM
            end if   
            if(predictif) then
               agr(g1,g2) = theta(1)  ! en MPa
               bgr(g1,g2) = theta(2)  ! en MPa
               
               agr(g2,g1) = agr(g1,g2)
               bgr(g2,g1) = bgr(g1,g2)
            endif
         elseif(n_param == 4) then
            read(541,*) texte,theta(1)
            read(541,*) texte,theta(3)
            read(541,*) texte,theta(2)
            read(541,*) texte,theta(4)
c           on dé-tilde les Akl et Bkl
            if(tilded)then
               theta(1:4) = theta(1:4)/TildAB !modif tempo APM
            end if             
            if(predictif) then
               agr(g1,g2) = theta(1)  ! en MPa
               bgr(g1,g2) = theta(3)  ! en MPa
               agr(g1,g3) = theta(2)  ! en MPa
               bgr(g1,g3) = theta(4)  ! en MPa

               agr(g2,g1) = agr(g1,g2)
               agr(g3,g1) = agr(g1,g3)
               bgr(g2,g1) = bgr(g1,g2)
               bgr(g3,g1) = bgr(g1,g3)
            endif
         elseif(n_param == 6) then
            read(541,*) texte,theta(1)
            read(541,*) texte,theta(4)
            read(541,*) texte,theta(2)
            read(541,*) texte,theta(5)
            read(541,*) texte,theta(3)
            read(541,*) texte,theta(6)
c           on dé-tilde les Akl et Bkl
            if(tilded)then
               theta(1:6) = theta(1:6)/TildAB !modif tempo APM
            end if 
            if(predictif) then
               agr(g1,g2) = theta(1) ! en MPa
               bgr(g1,g2) = theta(4) ! en MPa
               agr(g1,g3) = theta(2) ! en MPa
               bgr(g1,g3) = theta(5) ! en MPa
               agr(g1,g4) = theta(3) ! en MPa
               bgr(g1,g4) = theta(6) ! en MPa

               agr(g2,g1) = agr(g1,g2)
               agr(g3,g1) = agr(g1,g3)
               agr(g4,g1) = agr(g1,g4)
               bgr(g2,g1) = bgr(g1,g2)
               bgr(g3,g1) = bgr(g1,g3)
               bgr(g4,g1) = bgr(g1,g4)
            endif
         endif
         do i = 1,6-n_param
            read(541,*) 
         enddo
      else
         do i = 1,6
            read(541,*) 
         enddo
      endif
c
c     MODELE PR-WILSON/UNIQUAC/NRTL
c     ********************************
c
      do i = 1,8
         read(541,*) 
      enddo
c     "   (2.1) SIMULATION"
      read(541,*) texte,answer
      if(answer.eq.1 .and. .not.kij_const.and..not.predictif
     &   .and. .not.fitting) then
         simul_PRWil = .TRUE.
      endif
      read(541,*) !Les valeurs initiales des BIPs doivent etre fournies ci-dessous.
c     (2.2) AJUSTEMENT AUTOMATIQUE - BFGS
      read(541,*) texte,answer     
      read(541,*) texte,n_param,texte,val_a
      if((answer.eq.1 .or. answer.eq.3)
     &    .and. .not.kij_const.and..not.predictif
     &   .and. .not.fitting .and. .not.simul_PRWil) then
	     
         if (answer.eq.3) maillage_PRWil = .true.
         fitting_PRWil = .TRUE.
         m_bfgs = n_param
      endif
      
      if(modele_gE == "HV (NRTL)") alphgE = val_a       
      
      read(541,*) !Les valeurs initiales des BIPs doivent etre fournies ci-dessous.
      read(541,*) !Les valeurs des parametres non ajustes doivent aussi etre fournies ci-dessous.
      read(541,*) !Les bornes min et max du domaine d'ajustement sont a renseigner en section 3.1
      
      do i = 1,5
         read(541,*)
      enddo
      if(simul_PRWil .or. fitting_PRWil) then
         do i = 1,6
            read(541,*) texte,theta(i)
         enddo
         agE0(1,2) = theta(1)
         agE0(2,1) = theta(2)
         bgE0(1,2) = theta(3)
         bgE0(2,1) = theta(4)
         cgE0(1,2) = theta(5)
         cgE0(2,1) = theta(6)
      else
         do i = 1,6
            read(541,*) 
         enddo
      endif

      if(.not.kij_const .and. .not.predictif.and.
     &   .not.fitting.and..not.fitbin.and.
     &   .not.simul_PRWil.and..not.fitting_PRWil) then
         write(lun1,145)
         write(lun2,145)
 145     format(/,1x,'Mettre OUI=1 a l''une des 5 propositions '
     &             'du fichier option.d')
         stop
      endif

c     PARAMETRES DES AJUSTEMENTS
c     **************************
c     (3.1) BORNES DES DOMAINES DE RECHERCHE DES PARAMETRES
      if(fitting) then 
         do i = 1,8
            read(541,*) 
         enddo
         if(n_param == 2) then
            read(541,*) texte,borneinf(1),texte,bornesup(1)
            read(541,*) texte,borneinf(2),texte,bornesup(2)
         elseif(n_param == 4) then
            read(541,*) texte,borneinf(1),texte,bornesup(1)
            read(541,*) texte,borneinf(3),texte,bornesup(3)
            read(541,*) texte,borneinf(2),texte,bornesup(2)
            read(541,*) texte,borneinf(4),texte,bornesup(4)         
         elseif(n_param == 6) then
            read(541,*) texte,borneinf(1),texte,bornesup(1)
            read(541,*) texte,borneinf(4),texte,bornesup(4)
            read(541,*) texte,borneinf(2),texte,bornesup(2)
            read(541,*) texte,borneinf(5),texte,bornesup(5)
            read(541,*) texte,borneinf(3),texte,bornesup(3)
            read(541,*) texte,borneinf(6),texte,bornesup(6)   
         endif
         do i = 1,6-n_param
            read(541,*) 
         enddo
         do i = 1,n_param
             if(borneinf(i).gt.bornesup(i))then
                write(lun1,*)'Problème de bornes - Regardez fichier '//
     &                        'option.d'
                write(lun2,*)'Problème de bornes - Regardez fichier '//
     &                        'option.d'
                stop
             end if
         end do
         do i = 1,6
            if(borneinf(i)<8887.d0 .or. borneinf(i)>8889.d0) then
c              on dé-tilde :
               borneinf(i) = borneinf(i)/TildAB
            endif

            if(bornesup(i)<8887.d0 .or. bornesup(i)>8889.d0) then
c              on dé-tilde :
               bornesup(i) = bornesup(i)/TildAB
            endif
         enddo
      endif

      if(fitting_PRWil) then 
         do i = 1,8
            read(541,*)
         enddo
         do i = 1,6
            read(541,*) texte,borneinf(i),texte,bornesup(i)
         enddo
      endif

      if(.NOT.fitting_PRWil .AND. .NOT.fitting) then 
         do i = 1,14
            read(541,*) 
         enddo
      endif

c   (3.2) DEFINITION DE LA FONCTION OBJECTIF
      poids_fob(1:10) = 0.d0
      do i = 1,15
         read(541,*) texte
         write(*,*) texte
      enddo
      read(541,*) texte,choix_h,texte,choix_cp
      do i = 1,13
         read(541,*) texte
         write(*,*) texte
      enddo      
      read(541,*) texte,choix_fob
      if(choix_fob<1 .or. choix_fob>3) then
         write(lun1,*) 'Choix Fob = ',choix_fob
         write(lun2,*) 'Choix Fob = ',choix_fob
         write(lun1,*) 'Choix Fob incorrect - stop'
         write(lun2,*) 'Choix Fob incorrect - stop'
         stop
      endif
      if(choix_fob==3) then
         read(541,*) 
         do i = 1,6
c           poids_fob(1) -> w_elv2
c           poids_fob(2) -> w_crit
c           poids_fob(3) -> w_az
c           poids_fob(4) -> w_LLV
c           poids_fob(5) -> w_cp
c           poids_fob(6) -> w_h
            read(541,*) texte,poids_fob(i)
         enddo
         if(SUM(DABS(poids_fob(1:6)))<1.d-5) then
            write(lun1,*) 'Fob : definir des poids > 0'
            write(lun2,*) 'Fob : definir des poids > 0'
            stop
         endif
      endif
      close(541)

      lg1 = 2
      if(g1.LT.10) lg1 = 1
      lg2 = 2
      if(g2.LT.10) lg2 = 1
      lg3 = 2
      if(g3.LT.10) lg3 = 1
      lg4 = 2
      if(g4.LT.10) lg4 = 1
      open(2619,status='scratch')
      write(2619,2111) g1,g2,g3,g4
 2111 format(4(i2,/))
      rewind(2619)
      read(2619,2113) tg1
      read(2619,2113) tg2
      read(2619,2113) tg3
      read(2619,2113) tg4
 2113 format(a)
      close(2619)
      if(lg1.eq.1) tg1(1:1)=tg1(2:2)
      if(lg2.eq.1) tg2(1:1)=tg2(2:2)
      if(lg3.eq.1) tg3(1:1)=tg3(2:2)
      if(lg4.eq.1) tg4(1:1)=tg4(2:2)

      lun = 50
      huit = 8888.0d+00
      itempc = 0
c
c     ouverture des fichiers de travail
c
      call system ('cls')
 2140 write(lun1,214)
 214  format(//,15x,'nom du fichier de donnees experimentales ? ',$)
c      read(5,215) fichier
       fichier = 'data.d'
      l1 = len_trim(fichier)
      if(fichier(l1-1:l1).ne.'.d') fichier = fichier(1:l1)//'.d'

 215  format(a)

      call longueur(fichier,l1)
      open(40,file=fichier(1:l1),status='old',iostat=err0)  !fichier de donnees
      if(err0.ne.0) then
         write(*,*) 'Nom de fichier incorrect : ',fichier
         goto 2140
      endif

      l2 = l1
      do i=1,l1
         if(fichier(i:i).eq.'.') then
            l2 = i - 1
            goto 216
         endif
      end do
 216  continue
      fichres = fichier(1:l2)//'.res' !une copie d'ecran - contient les details du calcul.
      l2 = l2 + 4
      open(lun2,file=fichres(1:l2))

      if(kij_const) then    !Calcul avec valeur fixee de k12
         predictif = .false.
         fitting = .false.
         simul_PRWil = .false.
         fitting_PRWil = .false.
      endif
c
c     initialisation de variables generales
c
      nimp = 0
      nco  = 2   
c
c     lecture du fichier de donnees pour determiner :
c     le nombre de series de points (nserie) tous systemes confondus
c     le nombre de series de points (nser) par systeme
c     et le nombre de donnees experimentales (nexp) tous systemes confondus
c
      nserie = 0
      nexp = 0
      npt_x = 0
      npt_y = 0
      npt_c = 0
      npt_az = 0
      npt_LLV = 0
      npt_h = 0
      npt_cp = 0
      read(40,*) nsyst  !nombre de systemes binaires differents
      if(nsyst.gt.nsysmax) then
         write(lun1,5431) nsyst,nsysmax
 5431    format(1x,'le fichier de donnees contient ',i3,' systemes',/,
     &   1x,'augmenter la valeur de nsysmax qui vaut:',i3)
         stop
      endif
      do i=1,nsyst
         compteur = 0
         nser(i) = 0
         read(40,*)   !ligne d'etoiles
         read(40,*) npoin(i)  !nbre de points experimentaux pour le systeme
         ndat(i) = 0
         read(40,*)  !code + nom
         read(40,*)  !code + nom

         group_nb(:) = 0
         read(40,*) codc1(i),nomc1(i),nbg,(group_nb(kk),poub,kk=1,nbg)
         if(choix_eos.eq.1) then ! PPR78/E-PPR78
            if(MAXVAL(group_nb).gt.nbgr_tot_max) then
               print*,group_nb
               write(lun1,*) 'Augmenter nbgr_tot_max dans COMFLASH'
               write(lun2,*) 'Augmenter nbgr_tot_max dans COMFLASH'
               STOP
            endif
   
            if(MAXVAL(group_nb).gt.nbgr_tot) then
               write(lun1,*) 'Groupes manquants dans TABINGR PPR/EPPR'
               write(lun2,*) 'Groupes manquants dans TABINGR PPR/EPPR'
               STOP
            endif
         endif
         read(40,*) mmc1(i),tcc1(i),pcc1(i),omegac1(i)

         group_nb(:) = 0
         read(40,*) codc1(i+nsyst),nomc1(i+nsyst),nbg,
     &              (group_nb(kk),poub,kk=1,nbg)
         if(choix_eos.eq.1) then ! PPR78/E-PPR78
            if(MAXVAL(group_nb).gt.nbgr_tot_max) then
               write(lun1,*) 'Augmenter nbgr_tot_max dans COMFLASH'
               write(lun2,*) 'Augmenter nbgr_tot_max dans COMFLASH'
               STOP
            endif
   
            if(MAXVAL(group_nb).gt.nbgr_tot) then
               write(lun1,*) 'Groupes manquants dans TABINGR PPR/EPPR'
               write(lun2,*) 'Groupes manquants dans TABINGR PPR/EPPR'
               STOP
            endif
         endif
         read(40,*) mmc1(i+nsyst),tcc1(i+nsyst),pcc1(i+nsyst),
     &              omegac1(i+nsyst)
         read(40,*)  !etoiles
c
c        test de la volatilite
c
         if(tcc1(i).lt.tcc1(i+nsyst)) then
            tdet = 0.95*tcc1(i)
         else
            tdet = 0.95*tcc1(i+nsyst)
         endif

         psat0 = 7.0d+0/3.0d+00*(omegac1(i)+1.0d+0)
         psat1 = pcc1(i)*10.0d+00**(psat0*(1.0d+0-tcc1(i)/tdet))
         psat0 = 7.0d+0/3.0d+00*(omegac1(i+nsyst)+1.0d+0)
         psat2 = pcc1(i+nsyst)*1.d+1**(psat0*(1.d+0-tcc1(i+nsyst)/tdet))
         if(psat1.lt.psat2) then
            write(lun1,2000) nomc1(i),nomc1(i+nsyst)
            write(lun2,2000) nomc1(i),nomc1(i+nsyst)
 2000       format(1X,'Pour le systeme : ',A5,' + ',A5,'inverser SVP ',
     &             'l''ordre des constituants - MERCI -')
         endif

 220     read(40,*)       !reference bibliographique
         read(40,*,err=221) nom(1),nom(2),code,ndata
         if(ndata.gt.nptmax) then
            write(lun1,9657) ndata,nptmax
 9657       format(1x,'une serie contient ',i3,' points',/,
     &      1x,'alors que nptmax = ',i3)
            stop
         endif
         if(code.eq.1.or.code.eq.11.or.code.eq.2.or.code.eq.200.or.code.
     &eq.21.or.code.eq.210.or.code.eq.3.or.code.eq.300.or.code.eq.4.or.c
     &ode.eq.41.or.code.eq.42.or.code.eq.5.or.code.eq.51.or.code.eq.6.or
     &.code.eq.600.or.code.eq.61.or.code.eq.610.or.code.eq.7.or.code.eq.
     &700.or.code.eq.71.or.code.eq.710.or.code.eq.65.or.code.eq.66
     &.or.code.eq.67.or.code.eq.55.or.code.eq.56.or.code.eq.57
     &.or.code.eq.8.or.code.eq.9) then
            continue
         else
            write(lun1,224) i,compteur
 224        format(1x,'erreur pour le systeme numero : ',i3,/,
     &       1x,'verifiez le nombre de points de la serie : ',i3)
            stop            
         endif
         ndat(i) = ndat(i) + ndata   !nombre de points pour un systeme donne
         read(40,*)      ! t ou p + facteur(s) de conversion
         if(code.eq.1.or.code.eq.11.or.code.eq.2.or.code.eq.200.or.code.
     &eq.21.or.code.eq.210.or.code.eq.5.or.code.eq.51.or.code.eq.6.or.co
     &de.eq.600.or.code.eq.61.or.code.eq.610) then !x-data
            npt_x = npt_x + ndata
         endif
         if(code.eq.1.or.code.eq.11.or.code.eq.3.or.code.eq.300.or.code.
     &eq.5.or.code.eq.51.or.code.eq.7.or.code.eq.700.or.code.eq.71.or.co
     &de.eq.710) then !y-data
            npt_y = npt_y + ndata
         endif
         if(code.eq.4.or.code.eq.41.or.code.eq.42) then
            npt_c = npt_c + ndata
         endif
         if(code.eq.8) then
            npt_az = npt_az + ndata
         endif
         if(code.eq.9) then
            npt_LLV = npt_LLV + ndata
         endif
         if(code.eq.55.or.code.eq.56.or.code.eq.57) then
            npt_h = npt_h + ndata
         endif
         if(code.eq.65.or.code.eq.66.or.code.eq.67) then
            npt_cp = npt_cp + ndata
         endif
         nserie = nserie + 1
         nexp = nexp + ndata
         do j=1,ndata
            read(40,*)
         enddo
         if(ndat(i).ne.npoin(i)) then
            compteur = compteur + 1
            goto 220
         else
            nser(i) = compteur + 1
            if(nser(i).gt.nsermax) then
               write(lun1,5540) i,nser(i),nsermax
 5540          format(1x,'le systeme numero:',i3,' contient ',i3,
     &         ' series de points mais nsermax = ',i3)
               stop
            endif
         endif
      end do
      goto 223
 221  write(lun1,222) i,npoin(i),ndat(i)
 222  format(1x,'erreur pour le systeme numero : ',i3,/,
     &       1x,'le nombre de donnees declarees est : ',i4,/,
     &       1x,'le programme en trouve : ',i4)
      stop
 223  continue
c
c     verification de l'homogeneite des mm,tc,pc,omega
c     determination du nombre de corps purs differents
c

      ii = 0
      do i=1,2*nsyst-1
         passage = 0
         do j=i+1,2*nsyst
            if(codc1(i).eq.codc1(j)) then
               if(mmc1(i).ne.mmc1(j)) then
                  write(lun1,320)
                  write(lun2,320)
 320              format(1x,'mm non homogenes')
                  l1 = len_trim(nomc1(i))
                  l2 = len_trim(nomc1(j))
                  write(lun1,*) nomc1(i)(1:l1),' Mw = ',mmc1(i)
                  write(lun1,*) nomc1(j)(1:l2),' Mw = ',mmc1(j)
                  write(lun2,*) nomc1(i)(1:l1),' Mw = ',mmc1(i)
                  write(lun2,*) nomc1(j)(1:l2),' Mw = ',mmc1(j)
c                  stop
               endif
               if(tcc1(i).ne.tcc1(j)) then
                  write(lun1,321)
                  write(lun2,321)
 321              format(1x,'tc non homogenes')
                  l1 = len_trim(nomc1(i))
                  l2 = len_trim(nomc1(j))
                  write(lun1,*) nomc1(i)(1:l1),' Tc = ',tcc1(i)
                  write(lun1,*) nomc1(j)(1:l2),' Tc = ',tcc1(j)
                  write(lun2,*) nomc1(i)(1:l1),' Tc = ',tcc1(i)
                  write(lun2,*) nomc1(j)(1:l2),' Tc = ',tcc1(j)
c                  stop
               endif
               if(pcc1(i).ne.pcc1(j)) then
                  write(lun1,322)
                  write(lun2,322)
 322              format(1x,'pc non homogenes')
                  l1 = len_trim(nomc1(i))
                  l2 = len_trim(nomc1(j))
                  write(lun1,*) nomc1(i)(1:l1),' Pc = ',pcc1(i)
                  write(lun1,*) nomc1(j)(1:l2),' Pc = ',pcc1(j)
                  write(lun2,*) nomc1(i)(1:l1),' Pc = ',pcc1(i)
                  write(lun2,*) nomc1(j)(1:l2),' Pc = ',pcc1(j)
c                  stop
               endif
               if(omegac1(i).ne.omegac1(j)) then
                  write(lun1,322)
                  write(lun2,322)
 323              format(1x,'omega non homogenes')
                  l1 = len_trim(nomc1(i))
                  l2 = len_trim(nomc1(j))
                  write(lun1,*) nomc1(i)(1:l1),' om = ',omegac1(i)
                  write(lun1,*) nomc1(j)(1:l2),' om = ',omegac1(j)
                  write(lun2,*) nomc1(i)(1:l1),' om = ',omegac1(i)
                  write(lun2,*) nomc1(j)(1:l2),' om = ',omegac1(j)
c                  stop
               endif
               if(passage.eq.0) then
                  passage = 1
                  if(ii.ne.0) then
                     ajout = 1
                     do jj=1,ii
                        if(codc1(i).eq.cod_ok(jj)) then !on ne le rajoute pas
                           ajout = 0
                        endif
                     end do
                     if(ajout.ne.0) then
                        ii = ii + 1
                        cod_ok(ii) = codc1(i)
                        mm_ok(ii) = mmc1(i)
                        nom_ok(ii) = nomc1(i)
                        tc_ok(ii) = tcc1(i)
                        pc_ok(ii) = pcc1(i)
                        om_ok(ii) = omegac1(i)
                     endif
                  else
                     ii = ii + 1
                     cod_ok(ii) = codc1(i)
                     mm_ok(ii) = mmc1(i)
                     nom_ok(ii) = nomc1(i)
                     tc_ok(ii) = tcc1(i)
                     pc_ok(ii) = pcc1(i)
                     om_ok(ii) = omegac1(i)
                  endif
               endif
            endif
            if(j.eq.(2*nsyst)) then
               if(passage.eq.0) then
                  passage = 1
                  if(ii.ne.0) then
                     ajout = 1
                     do jj=1,ii
                        if(codc1(i).eq.cod_ok(jj)) then !on ne le rajoute pas
                           ajout = 0
                        endif
                     end do
                     if(ajout.ne.0) then
                        ii = ii + 1
                        cod_ok(ii) = codc1(i)
                        mm_ok(ii) = mmc1(i)
                        nom_ok(ii) = nomc1(i)
                        tc_ok(ii) = tcc1(i)
                        pc_ok(ii) = pcc1(i)
                        om_ok(ii) = omegac1(i)
                     endif
                  else
                     ii = ii + 1
                     cod_ok(ii) = codc1(i)
                     mm_ok(ii) = mmc1(i)
                     nom_ok(ii) = nomc1(i)
                     tc_ok(ii) = tcc1(i)
                     pc_ok(ii) = pcc1(i)
                     om_ok(ii) = omegac1(i)
                  endif
               endif
            endif
         end do
      end do
      i=2*nsyst !on teste le dernier code.
      ajout = 1
      do jj=1,ii
         if(codc1(i).eq.cod_ok(jj)) then !on ne le rajoute pas
            ajout = 0
         endif
      end do
      if(ajout.ne.0) then
         ii = ii + 1
         cod_ok(ii) = codc1(i)
         mm_ok(ii) = mmc1(i)
         nom_ok(ii) = nomc1(i)
         tc_ok(ii) = tcc1(i)
         pc_ok(ii) = pcc1(i)
         om_ok(ii) = omegac1(i)
      endif
c
c     on a donc ii corps purs differents
c
      write(lun2,600)
 600  format(/,1x,'proprietes physico-chimiques des differents corps'//
     & ' purs :',/,1x,'  nom  ',1x,'mw(g/mol)',1x,' tc/K  ',2x,'pc/bar '
     & ,1x,' omega ')
 601  format(1x,a7,2x,f7.3,2x,f7.3,1x,f7.3,1x,f7.4)
c
c     avant d'afficher les proprietes des corps purs, on les trie par masse
c     molaire croissante
c
c      call trie_data(mm_ok,ii,nom_ok,tc_ok,pc_ok,om_ok)

      call trie_data2(mm_ok,ii,nom_ok,tc_ok,pc_ok,om_ok,cod_ok) !modif APM
c      goto 2811
      write(lun1,*)''
      write(lun2,*)''
      write(lun1,*)'VERIFICATION DE DONNEES'
      write(lun2,*)'VERIFICATION DE DONNEES'      
      do i=1,ii
         j = cod_ok(i)
         tempo_const = bdd%molecule(j)%cstvalue(mw_idx)%val
         if(mm_ok(i).ne.tempo_const) then
            write(lun1,*)cod_ok(i),' masse mol. AVANT / APRES : ',
     &               mm_ok(i),tempo_const
            write(lun2,*)cod_ok(i),' masse mol. AVANT / APRES : ',
     &               mm_ok(i),tempo_const   
         endif                      
         mm_ok(i) = bdd%molecule(j)%cstvalue(mw_idx)%val
         
         tempo_const = bdd%molecule(j)%cstvalue(tc_idx)%val
         if(tc_ok(i).ne.tempo_const) then
             write(lun1,*)cod_ok(i),' Tc AVANT / APRES : ',
     &              tc_ok(i),tempo_const
             write(lun2,*)cod_ok(i),' Tc AVANT / APRES : ',
     &              tc_ok(i),tempo_const    
         endif                      
         tc(i) = bdd%molecule(j)%cstvalue(tc_idx)%val
         
         tempo_const = bdd%molecule(j)%cstvalue(pc_idx)%val
         if(pc_ok(i).ne.tempo_const) then
            write(lun1,*)cod_ok(i),' Pc AVANT / APRES : ',
     &               pc_ok(i),tempo_const
            write(lun2,*)cod_ok(i),' Pc AVANT / APRES : ',
     &               pc_ok(i),tempo_const      
         endif             
         pc(i) = bdd%molecule(j)%cstvalue(pc_idx)%val         

         tempo_const = bdd%molecule(j)%cstvalue(acen_idx)%val
         if(om_ok(i).ne.tempo_const) then
            write(lun1,*)cod_ok(i),' OMEGA AVANT / APRES : ',
     &               om_ok(i),tempo_const
            write(lun2,*)cod_ok(i),' OMEGA AVANT / APRES : ',
     &               om_ok(i),tempo_const      
         endif            
         om_ok(i) = bdd%molecule(j)%cstvalue(acen_idx)%val
            
         write(lun2,601) nom_ok(i),mm_ok(i),tc_ok(i),pc_ok(i),om_ok(i)
      end do
c      read(*,*)
c 2811 continue     
      close(40)  !fermeture du fichier de donnees

      npt_tot = npt_x + npt_y + npt_c + npt_h + npt_cp
     &          + npt_az + npt_LLV
      write(lun1,1000) nserie,npt_x,npt_y,npt_c,npt_az,npt_LLV,
     &                 npt_h,npt_cp,npt_tot
      write(lun2,1000) nserie,npt_x,npt_y,npt_c,npt_az,npt_LLV,
     &                 npt_h,npt_cp,npt_tot
 1000 format(/,1x,'le fichier de donnees contient : ',i4,
     &      ' series de points',/,
     &       1x,i5,' points de bulle',/,
     &       1x,i5,' points de rosee',/,
     &       1x,i5,' points critiques',/,
     &       1x,i5,' azeotropes',/,
     &       1x,i5,' donnees LLV',/,
     &       1x,i5,' donnees hM',/,
     &       1x,i5,' donnees cPM',/,
     &       1x,'-----',/,
     &       1x,i5,' points experimentaux')
     
      if(npt_h.gt.0) then ! au moins une donnee hM
c                           => on lit le fichier des Cp gaz parfait
c                              pour pouvoir calculer les
c                              DELTA T = DELTA hM / cP
         open(78,file='..\source\cpGP.dat')
         indice_cp(:) = 0 ! init

         do j = 1,ncompmaxcp
c           ! tabcpgp(j,1) contient le code du constituant sur la jieme
c           !              ligne de cpGP.dat
c           ! tabcpgp(j,2:6) contient les 5 coeff. de Cp sur la jieme
c           !                ligne de cpGP.dat
            read(78,*,END=54) tabcpgp(j,1:6)
c           vecteur qui associe le code du constituant au numero de la ligne
c           correspondante dans cpGP.dat
c           EXEMPLE :
c           indice_cp(cod(1)) renvoie le nø de ligne du tableau tabcpgp
c           qui contient les coeff du CP gaz parfait du constituant
c           associe au code cod(1).
            if(int(tabcpgp(j,1)).gt.ncompmaxcp2) then
               write(lun1,*) 'Augmenter ncompmaxcp2 dans COMFLASH '//
     &                       'ou (mieux) diminuer le code constituant'
               write(lun2,*) 'Augmenter ncompmaxcp2 dans COMFLASH '//
     &                       'ou (mieux) diminuer le code constituant'
               stop
            endif
            indice_cp(int(tabcpgp(j,1))) = j

            if(j.eq.ncompmaxcp) then
               write(lun1,*) 'Augmenter ncompmaxcp dans COMFLASH'
               write(lun2,*) 'Augmenter ncompmaxcp dans COMFLASH'
               stop
            endif
         enddo
 54      continue

         close(78)
      endif
c
      write(lun1,*) 
      call func_min(n_param,theta,fobj,ecrit)
      
      call date(date_out)
      open(3000,status='scratch')
      write(3000,*) date_out(1:2) !le mois
      write(3000,*) date_out(4:5) !le jour
      rewind(3000)
      read(3000,*) mois_out
      read(3000,*) jour_out
      close(3000)
      call timer(time_out)
      time_out = int(time_out/100.0d+00) !secondes
      if(date_out.eq.date_in) then       ! calcul effectue sur 1 seule journee
         njour = 0
      else                               !calcul etale sur plusieurs jours
         if(mois_out.eq.mois_in) then    !pas de changement de mois.
            njour = jour_out - jour_in
         else
            if(mois_in.eq.1.or.mois_in.eq.3.or.mois_in.eq.5.or.
     &         mois_in.eq.7.or.mois_in.eq.8.or.mois_in.eq.10.or.
     &         mois_in.eq.12) then !mois de 31 jours
               njour = jour_out - jour_in + 31
             elseif(mois_in.eq.2) then !fevrier
                njour = jour_out - jour_in + 28
             else
                njour = jour_out - jour_in + 30
             endif
          endif
      endif

      duree_s = time_out - time_in + njour*24.0d0*3600.0d0

      n_jour = int(duree_s/(3600.0D+0*24.0d0))
      n_hour = int((duree_s - (n_jour*3600.0D+0*24.0d0))/3600.0D+0)
      n_min = int((duree_s - (n_jour*3600.0D0*24.0d0)-(n_hour*3600.0D0))
     &        /60.0D+0)
      n_sec = duree_s - (n_jour*3600.0D+0*24.0d0)- (n_hour*3600.0D+0)
     &        - (n_min*60.0d0)
      write(lun1,8880) n_jour,n_hour,n_min,n_sec
      write(lun2,8880) n_jour,n_hour,n_min,n_sec
 8880 FORMAT(/,1X,70('*'),/,2X,'LE PROGRAMME A TOURNE : ',I2,' jours ',
     &I2,' heures ',
     &I2,' minutes et ',I2,' secondes',/,1X,70('*'),/)
     
      close(220)

      end
c
c---------------------------------------------------------------------------
c
      subroutine func_min(n,theta,f,ecrit)
      use load_bdd_mod !modif tempo APM
      implicit none
      include '..\source\comflash.cmn'
c
c     Calcule la valeur de la fonction objectif f a minimiser
c     n = nombre de parametres. ARGUMENTS d'entree
c     theta(i) = les parametres i=1 a n. ARGUMENTS d'entree
c     f = valeur de la fonction objectif. ARGUMENT de sortie
c     ecrit = argument d'entree (si non nul, on ecrit les resultats i.e. FO)
c
      integer numf,n,no,i,j,k,ii,nexp,nsyst,ecrit,ecrit1
      integer kn,kk,code
      integer l1,l2,lun,jn,jj
      integer nbg,ncod,cod(pmax)
      integer npt_x,npt_y,npt_c,npt_tot,npt_out,dew_rej,bub_rej
      integer npt_h,npt_cp,npt_az,npt_LLV,nelv
      integer jt,jp,er1,er2,iii,rp,npt_deltaT
      integer npt_pn
      integer npoint(nsermax),nbphases
      integer i1,i2

      real*8 theta(n),f
      real*8 nflash(pmax),temp(nsermax),p(nsermax,nptmax)
      real*8 pres(nsermax)
      real*8 t(nsermax,nptmax)
      real*8 xl(nsermax,nptmax),xg(nsermax,nptmax),pression
      real*8 xl2(nsermax,nptmax)
      real*8 hM(nsermax,nptmax),cPM(nsermax,nptmax),xM(nsermax,nptmax)
      real*8 temperature,factp(nsermax)
      real*8 factt1(nsermax),factt2(nsermax)
      real*8 ecartx,ecarty,ecartpc,ecartpca,ecartxc
      real*8 fo_x,fo_y,fo_pc,fo_pca,fo_xc,fo_h,fo_cp
      real*8 fo_pLLV,fo_xLLV,fo_pLLVa,fo_xLLVa
      real*8 ecartxLLV,ecartxLLV2,ecartxLLV3
      real*8 ecartxLLVa,ecartxLLVa2,ecartxLLVa3
      real*8 ecartpLLV,ecartpLLVa
      real*8 sk(nbgr_tot_max),tdet,press
      real*8 ecartxa,ecartya,ecartxca,fo_xa,fo_ya,fo_xca,pun,pze
      real*8 ecarth,ecartcp,fo_ha,fo_cpa,ecartha,ecartcpa
      real*8 ecarth2,ecartcp2,fo_h1,fo_cp1,fo_h2,fo_cp2
      real*8 xb(npmax),xr(npmax),pb(npmax),pr(npmax)
      real*8 tb(npmax),tr(npmax)
      real*8 nume,nume2,nume3,deltaT_hM,xcomp(2),cpmcalc,hmcalc,Cpec
      real*8 pcrit,xcrit,pc_min,pc_max,xc_max,xc_min
      real*8 ecartxaz,ecartpaz
      real*8 ecartxaza,ecartpaza
      real*8 fo_xaz,fo_paz
      real*8 fo_xaza,fo_paza,felv
      real*8 psat1,psat2,patm,corr,Cpreal
      real*8 cpgaz(2),coeffCp(5),Cpmel,deltaT
      real*8 x1L,x1G,hmL,hmG
      real*8 felv2,fcrit,faz,fcp,fh,fLLV
      real*8 ecartpn,ecartpna,fo_pn,fo_pna,ppmax
      real*8 untiers,tempo_const
      real*8 parWilson(200,2)
c      real*8 pLLV(nseriemax),tLLV(nseriemax) ! données exp
c      real*8 xLLV(nseriemax),xLLV2(nseriemax),xLLV3(nseriemax) ! données exp

c     test vthermo, thermo subroutines
      integer ierr
      real*8 dp_dt,dp_dv,fg(pmax),ft(pmax),fv(pmax)
      real*8 fx(pmax,pmax),fp(pmax)

      real*8  zz(pmax),zmix,fug(pmax),fugt(pmax),fugp(pmax)
      real*8  fugx(pmax,pmax),aux(13),vtest


      
      parameter(patm=101325.d-5) ! bar

      character*80 name
      character*15 nomcomp,nomcomp1,nomcomp2
      character*60 tit
      character*7 bul_ros
      character*8 reject
      character*180 biblio
      character*28 type_fh,type_fcp
      character*3  type_pn 

      logical bug,is_EGG
      integer codeCP
      common/eq_gasgas/is_EGG
      common/compCP/codeCP
      common/ectpen/ppmax
      
c      open(200,file='calc_pres.txt')
      
      nombre_flash = 0
      eval_fo = .true.
      bul_ros = 'undefin'
      reject = 'indefini'
      bug = voir_bug
      if(ecrit.eq.0) then
         ecrit1 = 0             !on n'ecrit ni la FO ni les ecarts > ecmax0
         voir_bug = .false.     !on n'indique pas les erreurs des diagrammes
      endif
      if(ecrit.ne.0) then !on ecrira la FO
         if(voir_ecart) then
            ecrit1 = 1 !on ecrit les ecarts > ecmax0
         else
            ecrit1 = 0 !on n'ecrit PAS les ecarts > ecmax0
         endif
      endif

      appel_FO = appel_FO + 1
      nimp = 0
      pun = 1.d0 - 1d-8 !presque 1 ; autrefois : 0.99999d+00
      pze = 1.d-8       !presque zero ; autrefois : 0.00001d+00

      numf = 40

      call longueur(fichier,l1)
      open(numf,file=fichier(1:l1),status='old') !variable fichier dans common

      
      if(fitting .OR. fitbin) then  !Ajustement automatique des parametres du modele PPR78
         if(n.eq.2) then
            agr(g1,g2) = theta(1)
            agr(g2,g1) = agr(g1,g2)
            bgr(g1,g2) = theta(2)
            bgr(g2,g1) = bgr(g1,g2)
         elseif(n.eq.4) then
            agr(g1,g2) = theta(1)
            agr(g2,g1) = agr(g1,g2)
            agr(g1,g3) = theta(2)
            agr(g3,g1) = agr(g1,g3)
            bgr(g1,g2) = theta(3)
            bgr(g2,g1) = bgr(g1,g2)
            bgr(g1,g3) = theta(4)
            bgr(g3,g1) = bgr(g1,g3)
         elseif(n.eq.6) then
            agr(g1,g2) = theta(1)
            agr(g2,g1) = agr(g1,g2)
            agr(g1,g3) = theta(2)
            agr(g3,g1) = agr(g1,g3)
            agr(g1,g4) = theta(3)
            agr(g4,g1) = agr(g1,g4)
            bgr(g1,g2) = theta(4)
            bgr(g2,g1) = bgr(g1,g2)
            bgr(g1,g3) = theta(5)
            bgr(g3,g1) = bgr(g1,g3)
            bgr(g1,g4) = theta(6)
            bgr(g4,g1) = bgr(g1,g4)
         endif
      endif

      if(fitting_PRWil) then  !Ajustement automatique des parametres du modele PR-WILSON
         if(n.eq.2) then
            agE0(1,2) = theta(1)
            agE0(2,1) = theta(2)
         elseif(n.eq.4) then
            agE0(1,2) = theta(1)
            agE0(2,1) = theta(2)
            bgE0(1,2) = theta(3)
            bgE0(2,1) = theta(4)
         elseif(n.eq.6) then
            agE0(1,2) = theta(1)
            agE0(2,1) = theta(2)
            bgE0(1,2) = theta(3)
            bgE0(2,1) = theta(4)
            cgE0(1,2) = theta(5)
            cgE0(2,1) = theta(6)
         endif
      endif

      read(numf,*)nsyst  !nbre de systemes binaires
      read(numf,*)       !etoiles


      do k=1,nsyst
        tmemc = 0.0d0
         read(numf,*) nexp
         read(numf,*) cod(1),nom(1)  !comp
         read(numf,*) cod(2),nom(2)  !comp              
                 
         name = nom(1)
         call longueur(name,l1)
         name = nom(2)
         call longueur(name,l2)
c         write(lun1,1313),k,nom(1)(1:l1),nom(2)(1:l2)
 1313    format(1x,'systeme number : ',i2,' (',a,'-',a,')')

         do i=1,nco
            do j=1,nbgr_tot
               sk(j) = 0.0d+00
            end do
            read(numf,*)ncod,tit,nbg,(j,sk(j),kk=1,nbg)  !hp.d

c           On relit la ligne qu'on vient de lire et l'on stocke
c           les informations dans d'autres variables utilisees par
c           la subroutine TEMPE (pour le calcul des hM / CPM)
            backspace(numf)
            if(i.eq.1) then
               read(numf,*) codcomp1,comp1,ngr1,(nomgr1(kk),
     &                      occurgr1(kk),kk=1,ngr1)
            else
               read(numf,*) codcomp2,comp2,ngr2,(nomgr2(kk),
     &                      occurgr2(kk),kk=1,ngr2)     
            endif

            if(ncod.ne.cod(i)) then
               write(lun1,*) 'probleme de code'
               write(lun2,*) 'probleme de code'
               stop
            endif
            read(numf,*) mmol(i),tc(i),pc(i),fac(i),
     &      zra(i)   !hp.d                 
c
            ss(i) = 0.0d+00
            do j=1,nbgr_tot
               ss(i) = ss(i) + sk(j)   !sk(j) lus dans hp.d
            end do
            do j=1,nbgr_tot
               if(ss(i).eq.0) then
                  write(lun1,*) 'ss = 0 - PROBLEME !'
                  write(lun2,*) 'ss = 0 - PROBLEME !'
                  stop
               endif
               srg(i,j) = sk(j)/ss(i)
            end do
         end do

c        -- TABINGR light (portion de prog utilisee par TEMPE)
         liste_gr = 0
         nbgr_used = 0
         do i = 1,ngr1
            liste_gr(i) = nomgr1(i)
         enddo
         kk = 0
         nbgr_used = ngr1
         do i = 1,ngr2
            do j = 1,ngr1
               if(nomgr2(i).eq.liste_gr(j)) then
                  goto 3333
               endif
            enddo
            kk = kk + 1
            nbgr_used = nbgr_used + 1
            liste_gr(ngr1+kk) = nomgr2(i)
 3333       continue
         enddo
         agr0 = 0.d0
         bgr0 = 0.d0
         agr0(1:nbgr_tot_max,1:nbgr_tot_max)
     &              = agr(1:nbgr_tot_max,1:nbgr_tot_max)
         bgr0(1:nbgr_tot_max,1:nbgr_tot_max)
     &              = bgr(1:nbgr_tot_max,1:nbgr_tot_max)
         corr = 1.d6 ! MPa -> Pa
         do i=1,nbgr_tot_max
            do j=1,nbgr_tot_max
               agr0(i,j) = agr0(i,j) * corr ! conversion en uSI
               bgr0(i,j) = bgr0(i,j) * corr ! conversion en uSI
            end do
         end do
c        -- Fin TABINGR light
c
c       Modif tempo APM 
        do i=1,nco
          j = cod(i)                     
          mmol(i) = bdd%molecule(j)%cstvalue(mw_idx)%val                     
          tc(i) = bdd%molecule(j)%cstvalue(tc_idx)%val          
          pc(i) = bdd%molecule(j)%cstvalue(pc_idx)%val           
          fac(i) = bdd%molecule(j)%cstvalue(acen_idx)%val          
          zra(i) = bdd%molecule(j)%cstvalue(zr_idx)%val
          if(FoncAlpha.EQ.'TWU91')then
             parL(i) = bdd%molecule(j)%parTwu(PR_idx,1)
             parM(i) = bdd%molecule(j)%parTwu(PR_idx,2)            
             parN(i) = bdd%molecule(j)%parTwu(PR_idx,3)                  
          elseif(FoncAlpha.EQ.'TWU88')then
             parL(i) = 0.0728d0 + 0.6693d0*fac(i) + 0.0925d0*fac(i)**2
             parM(i) = 0.8788d0 - 0.2258d0*fac(i) + 0.1695d0*fac(i)**2            
             parN(i) = 2.0d0            
          end if                      
          qi(i) = bdd%molecule(j)%cstvalue(qi_idx)%val
          parC(i) = bdd%molecule(j)%trans(PR_idx)
        end do
c       ----------------
         call eqet0   ! anciennement : CALL NEWFLU
         itempc = 0   ! on change de systeme => on vide la memoire des EN
         t0 = 0.0d+00

         read(numf,*)        !serie d'etoiles
c
c        nser(k) est le nombre de serie de points pour le binaire k
c
         do i=1,nser(k)     !boucle sur les differentes series de points

            read(numf,*) biblio   !reference biblio
            read(numf,*) nomcomp,nomcomp,code,npoint(i)

c            read(numf,*)    !reference biblio
c            read(numf,*) nomcomp,nomcomp,code,npoint(i)            
c            print*,'i=',i,' code=',code

            if(code.eq.1.or.code.eq.11.or.code.eq.2.or.
     &         code.eq.200.or.code.eq.21.or.code.eq.210.
     &         or.code.eq.3.or.code.eq.300.or.code.eq.4.or.
     &         code.eq.41.or.code.eq.42.or.code.eq.8.or.code.eq.9) then   !donnees isothermes
               read(numf,*) temp(i),factp(i)  !t/k et conversion de pression

               l1 = len_trim(nom(1))
               l2 = len_trim(nom(2))
               if(choix_eos==1) then
                  g1 = nomgr1(1)
                  g2 = nomgr2(1)
                  write(lun1,912) nom(1)(1:l1),nom(2)(1:l2),i,
     &                  nser(k),temp(i),agr(g1,g2),bgr(g1,g2)
                  write(lun2,912) nom(1)(1:l1),nom(2)(1:l2),i,
     &                  nser(k),temp(i),agr(g1,g2),bgr(g1,g2)
               elseif(choix_eos==2 .or. choix_eos==4) then
                  write(lun1,1912) nom(1)(1:l1),nom(2)(1:l2),i,
     &                  nser(k),temp(i),agE0(1,2),agE0(2,1)
                  write(lun2,1912) nom(1)(1:l1),nom(2)(1:l2),i,
     &                  nser(k),temp(i),agE0(1,2),agE0(2,1)
               endif
 912           format(1x,a,'-',a,t15,'(ELV)',t21,i3,'/',i3,'   T/K ='
     &               ,f7.2,' Ag1g2 =',f21.16,' Bg1g2 =',f21.16)
 1912          format(1x,a,'-',a,t15,'(ELV)',t21,i3,'/',i3,'   T/K ='
     &               ,f7.2,' A12 =',f25.16,' A21 =',f25.16) !modif tempo APM 

               write(*,*)
               do j=1,npoint(i)
                  nombre_flash = 0
 1001                format(1x,'corps pur decele systeme',i3,
     &               ' serie',i3,' /',i3)
                  if(code.eq.21) then !p,x isotherme
                     read(numf,*) p(i,j),xl(i,j)  !p,x1
                     if(code.eq.21) ygtx = .true.
                     if(code.eq.210) ygtx = .false.
                     if(xl(i,j).lt.pze.or.xl(i,j).gt.pun)then
                        write(lun1,1001) k,i,nser(k)
                        write(lun2,1001) k,i,nser(k)
                        stop
                     endif
                     xg(i,j) = 1.0d+00
                  endif
                  p(i,j) = p(i,j)*factp(i)
                  temperature = temp(i)
                  pression = p(i,j)
                  if(code.ne.4.and.code.ne.41.and.code.ne.42
     &               .and.code.ne.8.and.code.ne.9) then  !pas un point critique ou un azeo
                     nflash(1) = xl(i,j)
                     nflash(2) = 1.0d+00 - nflash(1)
                     
                     write(*,*) 'thermo call'
                    pression = pression * 0.1d0 !bar to MPa 
                    call thermo(1,0,1,temperature,pression,zmix,
     &                        nflash,fg,ft,fp,fx,aux)
     
     
                    write(*,*) 'Temperature (K) = ',temperature
                    write(*,*) 'Pressure (MPa)  = ',pression
                    write(*,*) 'Composition 1   = ',nflash(1)
                    write(*,*) 'z               = ',zmix
                    
                    write(*,*)
                    write(*,*) 'vthermo call'
                    vtest = 555d0 !total volume/r    (k*mol/mpa)
                    call vthermo(1,temperature,vtest,pression,dp_dt,
     &                         dp_dv,nflash,fg,ft,fp,fx,ierr)                    

                    write(*,*) 'Temperature (K) = ',temperature
                    write(*,*) 'Pressure (MPa)  = ',pression
                    write(*,*) 'Composition 1   = ',nflash(1)
                  else
                     write(*,*) 'Probleme de code - FIN'
                     STOP
                  endif
               end do !boucle sur le nombre de points
            else
               write(lun1,*) 'code non identifie'
               write(lun2,*) 'code non identifie'
               stop
            endif
            
         end do   !boucle sur le nombre de series de points
         read(numf,*)  !etoiles
      end do   !boucle sur le nombre de systemes binaires

      close(200)
      end
c
c--------------------------------------------------------------------
c
      subroutine trie(ntemper,tde,ntbltemp,tabltemp,nsermax,typev)
      implicit none
c
c     cherche le nombre d'elements differents d'un tableau (vecteur)
c
      integer numt,fff,filtre,nsermax
      integer ntemper,ntbltemp,typev
      real*8 ttt
      real*8 tde(nsermax),tabltemp(nsermax)
      logical newtemp
c
c     ntemper = nombre de temperatures en entree (certaines sont identiques)
c     tde(1 a ntemper) : tableau des temperatures en entree
c
c     ntbltemp = nombre de temperatures differentes (en sortie)
c     tabltemp(1 a ntbltemp) : tableau des temperatures en sortie
            
      ntbltemp = 0
      do numt = 1,ntemper
         newtemp = .false.
         ttt = tde(numt)
         if (ntbltemp.gt.0) then
            fff = filtre(ttt,tabltemp,ntbltemp,nsermax,typev)
            if (fff.eq.0) newtemp = .true.
         else
            newtemp = .true.
         end if
         if (newtemp) then
            ntbltemp = ntbltemp + 1
            tabltemp(ntbltemp) = ttt
         end if
      end do

      return
      end
c
c--------------------------------------------------------------------
c
       function filtre(nombre,liste_nombre,nbelem,nsermax,typev)
       implicit none
c
c       si filtre different de zero "nombre" (reel) est un element du tableau
c       liste_nombre qui renferme nbelem valeurs.
c
       integer filtre,nsermax,typev
       real*8 nombre,liste_nombre(nsermax)
       integer i,nbelem
       logical is_equal

       do i=nbelem+1,nsermax
          liste_nombre(i) = 0.0d+00
       end do

       i = 1
c       do while (nombre.ne.liste_nombre(i).and.i.le.nbelem)
       do while (.not.is_equal(typev,nombre,liste_nombre(i))
     &                                      .and.i.le.nbelem)
          i = i + 1
       end do
c       if (nombre.eq.liste_nombre(i)) then
       if(is_equal(typev,nombre,liste_nombre(i)))then
           filtre = i
       else 
           filtre = 0
       end if

       return 
       end
c
c--------------------------------------------------------------------
c
       function is_equal(typev,xx,yy)
       implicit none
       logical is_equal
       integer typev
       real*8 xx,yy
       
       if(typev.eq.1)then
         is_equal = (dabs(xx - yy).lt.1.0d-02)
       elseif(typev.eq.2)then
         if(xx.gt.1.0d00)then
           is_equal = (dabs(xx - yy).lt.1.0d-01)
         else
           is_equal = (dabs(xx - yy).lt.1.0d-02)
         end if
       else
         is_equal = (xx.eq.yy) 
       end if
       end function
c
c-----------------------------------------------------------------------
c
c
      subroutine trie_data(vecteur,nombre_d_elements,v2,v3,v4,v5)
      implicit none
c
c     tri des elements d'un tableau par algorithme de la bulle non optimise
c                     (range par ordre croissant)
c     on trie selon la variable vecteur, les autres tableaux v2,v3,v4,v5
c     suivent.
c
      integer i,j,nombre_d_elements
      real*8 vecteur(*),v3(*),v4(*),v5(*)
      character*80 v2(*)

      do i = 1,nombre_d_elements-1
         do j = 1,nombre_d_elements-i
            if (vecteur(j).gt.vecteur(j+1)) then
               call echange_reel(vecteur(j),vecteur(j+1))
               call echange_chaine(v2(j),v2(j+1))
               call echange_reel(v3(j),v3(j+1))
               call echange_reel(v4(j),v4(j+1))
               call echange_reel(v5(j),v5(j+1))
            end if
         end do
      end do

      return
      end
c
c-----------------------------------------------------------------------
c
c
      subroutine trie_data2(vecteur,nombre_d_elements,v2,v3,v4,v5,v6)
      implicit none
c
c     tri des elements d'un tableau par algorithme de la bulle non optimise
c                     (range par ordre croissant)
c     on trie selon la variable vecteur, les autres tableaux v2,v3,v4,v5
c     suivent.
c
      integer i,j,nombre_d_elements,v6(*)
      real*8 vecteur(*),v3(*),v4(*),v5(*)
      character*80 v2(*)

      do i = 1,nombre_d_elements-1
         do j = 1,nombre_d_elements-i
            if (vecteur(j).gt.vecteur(j+1)) then
               call echange_reel(vecteur(j),vecteur(j+1))
               call echange_chaine(v2(j),v2(j+1))
               call echange_reel(v3(j),v3(j+1))
               call echange_reel(v4(j),v4(j+1))
               call echange_reel(v5(j),v5(j+1))
               call echange_entier(v6(j),v6(j+1))
            end if
         end do
      end do

      return
      end
c
c-----------------------------------------------------------------------

c
c
      subroutine trie_bulle(vecteur,nombre_d_elements)
      implicit none
c
c     tri des elements d'un tableau par algorithme de la bulle non optimise
c                     (range par ordre croissant)
c
      integer i,j,nombre_d_elements
      real*8 vecteur(*)

      do i = 1,nombre_d_elements-1
         do j = 1,nombre_d_elements-i
            if (vecteur(j).gt.vecteur(j+1)) then
               call echange_reel(vecteur(j),vecteur(j+1))
            end if
         end do
      end do

      return
      end
c
c-----------------------------------------------------------------------
c
       subroutine echange_reel(a,b)
       real*8 a,b
       real*8 x
       x = a
       a = b
       b = x
       return
       end
c
c***************************************************************************
c
       subroutine echange_chaine(a,b)
       character*80 a,b
       character*80 x
       x = a
       a = b
       b = x
       return
       end
c
c***************************************************************************
c
       subroutine echange_entier(a,b)
       integer a,b
       integer x
       x = a
       a = b
       b = x
       return
       end
c
c**********************************************************************
c
       subroutine phas1(j,ider,istab,p00,x,eta0,v0,g0)
c
c  calcul des proprietes d'un systeme suppose monophasique, avec si
c  si istab=1 on retient la solution la plus stable
c  si istab=2 on etudie la stabilite par rapport a la diffusion.
c  si istab=10 (on veut choisir sa solution)
c
       implicit none
c
       include '..\source\comflash.cmn'
c
       integer nchol
       parameter ( nchol=pmax*(pmax+3)/2 )
c
c*** variables de la liste d'appel
c
       integer j,ider,istab
       real*8 p00,eta0,v0,g0
       real*8 x(pmax)
c
c*** variables internes
c
       integer i,i3r,isym,irtn
       real*8 eta1,g,eta_in
       real*8 ttt(nchol)
c
c
c-------------------------------------------------
c
c
c      Si perte de precision de l'ordinateur :
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
c
c         en sortie, on veut 1 valeur de eta differente
c
       endif
       n0 = pmax
       m0 = pmax
c
       call regles_melange(ider,8888,x)   
       if(j.eq.3) eta = eta0
       call rop(j,p00)       
       call fug      
       eta0 = eta
       v0 = v
       g0 = 0.0D+00
       do i=1,nco
          if(x(i).le.0.d0) then
             write(lun1,*) "Erreur x = 0 ou x < 0 (PHAS1)"
             write(lun1,*) "Systeme : ",nom(1)(1:len_trim(nom(1))),
     &                     " - ",nom(2)(1:len_trim(nom(2)))
             write(lun1,*) "T = ",t0,"P = ",p00
             write(lun1,*) "x1 = ",x(1),"x2 = ",x(2)
             write(lun2,*) "Erreur x = 0 ou x < 0 (PHAS1)"
             write(lun2,*) "Systeme : ",nom(1)(1:len_trim(nom(1))),
     &                     " - ",nom(2)(1:len_trim(nom(2)))
             write(lun2,*) "T = ",t0,"P = ",p00
             write(lun2,*) "x1 = ",x(1),"x2 = ",x(2)
          endif
          g0 = g0 + x(i)*(dlog(x(i)) + lphi(i))
       end do
       if(istab.eq.0) goto 10 ! on se contente de la solution
c                             ! trouvee en premier
       call rcub(i3r,eta1)
       if(i3r.eq.0) goto 10   ! si l'equation n'a qu'une racine
c
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
c
  10   if(ider.eq.1) call dertp(ider,x)
       if(istab.ne.2) return
c
       isym = 1
       irtn = 1                                         ! stabilite par rapport
       call chol(nco,nco,dmudn,x,ttt,x,isym,istab,irtn) ! a la diffusion
c                                                       ! istab = 0, stable
       return                                           ! istab = 1, instable
       end

c
c**********************************************************************
c
       subroutine regles_melange(ider,icalltempe,x)
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c      Variables qui doivent être spécifiées préalablement :
c      1) quadra_b (= .TRUE. si expression de b quadratique)
c      2) Lij(:,:) et puisb (intervenant dans l'expression de b)
c      3) modele_gE (chaîne de caractères indiquant le choix de la
c                    fonction gE/RT)
c      4) Lambda : constante universelle au denominateur du terme
c                  d'exces de la regle de melange sur le a
c
c      N.B. : en sortie de la subroutine, on a
c         am = paramètre a dans les uSI
c         bm = paramètre b en cm3/mol 
c
c      ider (entier) :
c          si ider > 0, calcul de daidx(:,:) pour permettre le calcul des
c                       derivees de lphi dans dertp
c
c      icalltempe (entier) :
c                 Si icalltempe < 8888, appel de tempe
c                 pour le calcul des grandeurs dépendant de T
c                 uniquement (a corps pur et ses dérivées,
c                 éventuellement, kij & Eij et leurs dérivées)
c
c                 Si icalltempe = 8888 : pas d'appel de tempe
c                 Si icalltempe = 0, 1 ou 2 : appel de tempe
c                 avec icalltempe, l'indice de dérivation 
c                 par rapport a T
c
c                 *** icalltempe doit être < 8888 si
c                 regles_melange est appelée indépendemment
c                 de phas1 (i.e., en dehors d'un calcul de
c                 flash) ; c'est typiquement le cas quand
c                 on veut calculer hM et cPM.
c
       implicit none
c
       include '..\source\comflash.cmn'
       integer   ider   ! ider=1 s'il faut mettre à jour les ln phi (Input)
       integer   i,j,k,l
       integer   icalltempe ! < 8888 pour mettre à jour propriétés dépendant
c                             de T ou pour calculer prop dérivées       

       real*8    x(pmax) ! (Input)
       real*8    bxen,aapur,sumb(pmax),unspuisb, truc,pro
c      Pour règles de mélange UNIQUAC-WILSON
       real*8    sompsi(pmax),psi(pmax),psi_sur_x(pmax)
       real*8    lnsompsi(pmax)
       real*8    mux(pmax),som,som_mux,terme2,terme3
       real*8    dsompsi_T(pmax),d2sompsi_T(pmax),dE_dx(pmax)
c      Pour règles de mélange NRTL
       real*8    dsom_ExpAlphTau_dT(pmax),som_ExpAlphTau(pmax)
       real*8    E_ij(pmax,pmax),dEij_dT(pmax,pmax),d2Eij_dT2(pmax,pmax)
       real*8    dEijT_sur_Eij(pmax,pmax),dei_dT(pmax),dasb_dT(pmax)
       real*8    dEijsRT1,dEijsRT2,dEijsRT3,d2EijsRT1,d2EijsRT2
       real*8    d2EijsRT3,d2EijsRT4
c

       if(icalltempe<8887) then
          call tempe(icalltempe)
       endif

c      ======================================================
c      Règle de mélange sur le paramètre b :
c      ======================================================
c      b linéaire est un cas particulier du b quadratique à condition
c      de choisir bij = 0.5*(b_i + b_j)
c      On pourrait donc utiliser une formulation générale (quadratique)
c      pour traiter le cas du b linéaire. 
c      On maintient pourtant la possibilité d'utiliser la formulation linéaire
c      directe (sans passer par la formulation quadratique) car elle
c      est moins couteuse en temps de calcul
c
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
c               -----------------------------------------
c               Expression de bij (formulation générale)
c               -----------------------------------------
                bij(i,j) = 0.5d0 * (b(i)**unspuisb + b(j)**unspuisb)
                bij(i,j) = bij(i,j)**puisb * (1.d0-Lij(i,j))
c
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
c         d(n.bm)/dn(i) = d(bm)/dx(i) - bm
          derbm(1:nco) = dbmdx(1:nco) - bm
          do i=1,nco
             do k=1,nco
                der2bm(i,k) = 2.d0*bij(i,k) - dbmdx(k) ! d derbm(i) / d x(k) 
             enddo
          enddo
       else               ! b lineaire
          derbm(1:nco) = dbmdx(1:nco) ! d(n.bm)/d n(i)
c                                     ! = d(bm)/dx(i) = b(i)
          der2bm(1:nco,1:nco) = 0.d0  ! d derbm(i) / d x(k) 
       endif
       bmix = bm*1.d-6/rgas ! uSI
c
c      ======================================================
c      Regle de melange sur la paramètre a
c      ======================================================
c
c      Calcul de somme des x(i).a(i)/[b(i)RT]
c      --------------------------------------
       aapur = 0.d0
       do i=1,nco
          aapur = aapur + asb(i)*x(i)
       enddo
c
c
       SELECT CASE(modele_gE)
c         Règle de mélange sur le a :          
c         ~~~~~~~~~~~~~~~~~~~~~~~~~~
c         (a/RTb) = [somme des x(i) * a(i)/RTb(i)] + E(T,x)
c         --> E(T,x) = fonction d'exces (adimensionnelle) 
c
          CASE("PPR78 (Van Laar)")
c            - Cas d'un kij constant             
c            - Cas (E-)PPR78 ...
c
c            Fonction d'exces de type Van Laar - Peneloux
c             --> ee = E(T,x) = -0.5 * double somme des
c                               x(i).x(j).b(i).b(j).[E(i,j)/RT]
c                               / somme des x(k).b(k)
             ee = 0.d0
             do i=1,nco
                e(i) = 0.d0
                do k=1,nco
                   bxen = bx(k)*en(k,i)  ! bxen = x(k).b(k).E(k,i)/RT
                   e(i) = e(i) + bxen    ! e(i) = somme des b(k).x(k).E(i,k)/(RT)
                end do
c             double somme x(i).x(k).b(i).b(k).E(i,k)/(RT) 
c             (sans dimension)
              ee = ee + e(i)*bx(i)
             end do
             ee = -0.5d0*ee/bm ! E(T,x), sans dimension
             
             if(icalltempe>0 .and. icalltempe<8888) then
c               Calcul des dérivées p/r à T (pour estimation
c               de hM / cPM par exemple) 
                dE_dT = 0.d0
                do i = 1,nco
                   truc = 0.d0
                   do j = 1,nco
                      truc = truc + bx(j)*dEijsRT(i,j)
                   enddo
                   dei_dT(i) = b(i)*truc*1.0d-06  !Facteur de conversion (cm3) > (m3) : 1.d-6
                   dE_dT = dE_dT + bx(i)*truc
                enddo
                
                            
c               les bx(:) sont en cm3/mol ; bm est en cm3/mol
c               dEijsRT est dans les uSI
c               Facteur de conversion (cm3) > (m3) : 1.d-6
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
                                      
c                  les bx(:) sont en cm3/mol ; bm est en cm3/mol
c                  d2EijsRT est dans les uSI
c                  Facteur de conversion (cm3) > (m3) : 1.d-6
                   d2E_dT2 = -0.5d-6*d2E_dT2 / bm  
                endif
             endif
c  
c            Calcul de aa (alpha = am/bmRT) et de ee (fonction d'exces E)
c            des aai(i) (= d n.alpha / dn(i))
c            et daidx(k) (= d aai(i) / d x(k))
             aa = ee + aapur ! alpha = am / (bm.RT), sans dimension
             am = aa * (bm*1.d-6) * rgas*t0 ! uSI
             amix = am/(rgas**2*t0) ! uSI
          
c            Calcul de d(n.alpha) / d n(i) note aa(i)
c            --> d(n.alpha) / d n(i)  
             do i=1,nco
                aai(i) = asb(i) - (ee*derbm(i) + e(i)*b(i))/bm                
             end do

c            d(alpha)/dx(i) 
             dalphadx(1:nco) = asb(1:nco) - (b(1:nco)*e(1:nco)
     &                                  + ee*dbmdx(1:nco))/bm
     
             if(ider.ne.0) then ! calcul de fonctions utiles pour les derivees
c              de ln phi par rapport aux fractions mol. 
c              Calcul des derivees des aai(i) par rapport
c              aux fractions molaires x(k)
                do i=1,nco
                   do k=1,nco
                      daidx(i,k) = -b(i)*b(k)*bm*en(i,k)
     &                     + derbm(i)*(b(k)*e(k) + 2.d0*ee*dbmdx(k))
     &                     - bm*ee*der2bm(i,k)
     &                     + dbmdx(k)*b(i)*e(i)
                      daidx(i,k) = daidx(i,k) / bm**2
                   end do
                end do
             endif

             if(icalltempe>0 .and. icalltempe<8888) then
c               Calcul des dérivées p/r à T (pour estimation
c               de hM / cPM par exemple)
c               Calcul de d alpha / dT :
                dalpha_dT = 0.d0
                do i = 1,nco
                   dalpha_dT = dalpha_dT + asb(i)*x(i)*dlna_T(i)
                enddo
c               Pour rappel aapur = somme des x(i)*asb(i)
                dalpha_dT = dalpha_dT - aapur/t0 + dE_dT

c               at = d amix / dT = d racine(am/R²/T) / dT (uSI)
                at = dalpha_dT*bmix

                if(icalltempe >= 2) then
c                  Calcul de d² alpha / dT² :
                   d2alpha_dT2 = 0.d0
                   do i = 1,nco
                      d2alpha_dT2 = d2alpha_dT2 + asb(i)*x(i)*d2lna_T(i)
                   enddo
                   d2alpha_dT2 = d2alpha_dT2 - 2.d0/t0*
     &                           (dalpha_dT - dE_dT) + d2E_dT2
c                  att = d at / dT (uSI) :
                   att = d2alpha_dT2*bmix
c                  Calcul de d² alpha / dTdxi :
                   do i = 1,nco
                      daai_dT(i) = asb(i)*(dlna_T(i)- 1.0d0/t0)
     &                            - (dei_dT(i) + dE_dT*derbm(i))/bm
                   enddo                                     
                endif
             endif
          CASE("HV (Uniquac-Wilson)")
c
c            Fonction d'exces de type Wilson - UNIQUAC
c             --> ee = E(T,x) = gE_res / Lambda
c                    gE_res = Wilson ou Uniquac
c
             ee = 0.d0
             som = 0.d0
             do i=1,nco
                som = som + x(i)*sigma(i)
                mux(i) = mu(i)*x(i)
             end do
             som_mux = SUM(mux)
             psi_sur_x(1:nco) = sigma(1:nco)/som
             psi(1:nco) = x(1:nco)*psi_sur_x(1:nco)

             do i=1,nco
                sompsi(i) = 0.d0
                do j=1,nco                
                   sompsi(i) = sompsi(i) + psi(j)*expA(i,j)
                end do
             end do
             lnsompsi(1:nco) = DLOG(sompsi(1:nco))
             
c            E(T,x), sans dimension
             ee = 0.d0
             DO i = 1,nco
                ee = ee - mux(i)*lnsompsi(i)
             ENDDO

             ee = ee / Lambda

             if(icalltempe>0 .and. icalltempe<8888) then
c               Calcul des dérivées p/r à T (pour estimation
c               de hM / cPM par exemple)
                dE_dT = 0.d0
                DO i = 1,nco
                   dsompsi_T(i) = 0.d0
                   DO j = 1,nco
                      dsompsi_T(i) = dsompsi_T(i) + psi(j)*
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
                            som = som + psi(k)*expA(i,k)*(
     &                          d2pargE(i,k) +dpargE(i,k)*(
     &                          dpargE(i,j) - dpargE(i,k)))
                         ENDDO
                         d2sompsi_T(i) = d2sompsi_T(i) + psi(j)*
     &                               expA(i,j)*som
                      ENDDO
                      d2E_dT2 = d2E_dT2 + mux(i)*d2sompsi_T(i)/
     &                          sompsi(i)**2
                   ENDDO
                   d2E_dT2 = d2E_dT2 / Lambda
                endif
             endif
c  
c            Calcul de aa (alpha = am/bmRT)
c            des aai(i) (= d n.alpha / dn(i))
c            et daidx(k) (= d aai(i) / d x(k))
             aa = ee + aapur ! alpha = am / (bm.RT), sans dimension
             am = aa * (bm*1.d-6) * rgas*t0 ! uSI
             amix = am/(rgas**2*t0) ! uSI
          
c            Calcul de d(n.alpha) / d n(i) = aa(i)
c            --> d(n.alpha) / d n(i)  
             do i=1,nco
                terme3 = 0.d0
                do l = 1,nco
                   terme3 = terme3 + mux(l)*expA(l,i)/sompsi(l)
                enddo
                terme2 = psi_sur_x(i)*(terme3 - Som_mux)
                dE_dx(i) = mu(i)*lnsompsi(i) + terme2
                if(ider.ne.0) then ! calcul de fonctions utiles 
c                  pour les derivees de ln phi par rapport aux 
c                  fractions mol. 
c                  Calcul des derivees des aai(i) par rapport
c                  aux fractions molaires x(k)
                   do k=1,i
                      som = 0.d0
                      do l = 1,nco
                         som = som + mux(l)*expA(l,i)*expA(l,k)
     &                               /sompsi(l)**2
                      enddo
                      daidx(i,k) = mu(i)*psi_sur_x(k)*(expA(i,k)/
     &                      sompsi(i) - 1.d0)
     &                      -psi_sur_x(k)*terme2 +psi_sur_x(i)*
     &                      (mu(k)*(expA(k,i)/sompsi(k) - 1.d0)
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

c            d(alpha)/dx(i) = d(n.alpha) / dn(i)
             dalphadx(1:nco) = aai(1:nco)
     
             if(icalltempe>0 .and. icalltempe<8888) then
c               Calcul des dérivées p/r à T (pour estimation
c               de hM / cPM par exemple)
c               Calcul de d alpha / dT :
                dalpha_dT = 0.d0
                do i = 1,nco
                   dalpha_dT = dalpha_dT + asb(i)*x(i)*dlna_T(i)
                enddo
c               Pour rappel aapur = somme des x(i)*asb(i)
                dalpha_dT = dalpha_dT - aapur/t0 + dE_dT
c               at = d amix / dT = d racine(am/R²/T) / dT (uSI)
                at = dalpha_dT*bmix

                if(icalltempe >= 2) then
c                  Calcul de d² alpha / dT² :
                   d2alpha_dT2 = 0.d0
                   do i = 1,nco
                      d2alpha_dT2 = d2alpha_dT2 + asb(i)*x(i)*d2lna_T(i)
                   enddo
                   d2alpha_dT2 = d2alpha_dT2 - 2.d0/t0*
     &                           (dalpha_dT - dE_dT) + d2E_dT2
c                  att = d at / dT (uSI) :
                   att = d2alpha_dT2*bmix
                endif
             endif
          CASE("HV (NRTL)")
c
c            Fonction d'exces de type NRTL
c             --> ee = E(T,x) = gE_res / Lambda
c                    gE_res = NRTL
c            
c            E(T,x), sans dimension
             ee = 0.d0
             do i=1,nco
                som_ExpAlphTau(i) = 0.d0
                do k=1,nco                
                   som_ExpAlphTau(i) = som_ExpAlphTau(i) +
     &                     x(k) * expA(k,i)                          
                end do
                do j=1,nco                
                   E_ij(i,j) = pargE(j,i)*expA(j,i)/som_ExpAlphTau(i)
                   ee = ee + x(i)*x(j)*E_ij(i,j)
                end do                
             end do
 
             ee = ee / Lambda
                              
             if(icalltempe>0 .and. icalltempe<8888) then
c               Calcul des dérivées p/r à T (pour estimation
c               de hM / cPM par exemple)
               dE_dT = 0.d0
               do i=1,nco
                  dsom_ExpAlphTau_dT(i) = 0.d0
                  do k=1,nco                
                     dsom_ExpAlphTau_dT(i) = dsom_ExpAlphTau_dT(i)
     &                         + x(k)*alphgE(k,i)* dpargE(k,i)*expA(k,i)                 
                  end do
                  do j=1,nco
                     dEij_dT(i,j) = expA(j,i)/som_ExpAlphTau(i)*(
     &                         dpargE(j,i)*(1.d0-pargE(j,i)*alphgE(j,i))  
     &                         + pargE(j,i)*dsom_ExpAlphTau_dT(i)
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
                       som = som + x(k)*alphgE(k,i)*expA(k,i)*(
     &                             d2pargE(k,i) - alphgE(k,i)      
     &                                                  *dpargE(k,i)**2)                
                     end do
                     do j=1,nco
                        terme2 = d2pargE(j,i)
     &                          - alphgE(j,i)*(dpargE(j,i)**2
     &                          + pargE(j,i)*d2pargE(j,i)) 
     &                          + pargE(j,i)*som/som_ExpAlphTau(i)                 
                        terme3 = dsom_ExpAlphTau_dT(i)/som_ExpAlphTau(i)
     &                          *(dpargE(j,i) + pargE(j,i)*
     &                         dsom_ExpAlphTau_dT(i)/som_ExpAlphTau(i))
                        d2Eij_dT2(i,j) = dEij_dT(i,j)*
     &                       (-alphgE(j,i)*dpargE(j,i) + 
     &                        dsom_ExpAlphTau_dT(i)/som_ExpAlphTau(i))
     &                       + expA(j,i)/som_ExpAlphTau(i)*
     &                                               (terme2 + terme3)    
                        d2E_dT2 = d2E_dT2 + x(i)*x(j)*d2Eij_dT2(i,j)
                     end do                
                   end do
                   d2E_dT2 = d2E_dT2 / Lambda
               endif              
             endif
c  
c            Calcul de aa (alpha = am/bmRT)
c            des aai(i) (= d n.alpha / dn(i))
c            et daidx(k) (= d aai(i) / d x(k))
             aa = ee + aapur ! alpha = am / (bm.RT), sans dimension
             am = aa * (bm*1.d-6) * rgas*t0 ! uSI
             amix = am/(rgas**2*t0) ! uSI
          
c            Calcul de d(n.alpha) / d n(i) = aa(i)
c            --> d(n.alpha) / d n(i)  
             do i=1,nco
                som = 0.d0 
                do j = 1,nco
                   som = som + x(j)*E_ij(i,j)
                enddo               
                terme2 = 0.d0
                do l = 1,nco
                   terme3 = 0.d0
                   do j = 1,nco
                      terme3 = terme3 + x(j)*E_ij(l,j)*expA(i,l)/
     &                                         som_ExpAlphTau(l) 
                   enddo
                   terme2 = terme2 + x(l)*(E_ij(l,i) - terme3)
                enddo
                dE_dx(i) = som + terme2
                if(ider.ne.0) then ! calcul de fonctions utiles 
c                  pour les derivees de ln phi par rapport aux 
c                  fractions mol. 
c                  Calcul des derivees des aai(i) par rapport
c                  aux fractions molaires x(k)

                   do k=1,i
                      terme2 = 0.d0
                      do l = 1,nco
                         som = 0.d0
                         do j = 1,nco
                            som = som + x(j)*E_ij(l,j)*expA(k,l)
     &                                                    *expA(i,l)
                         enddo
                         terme2 = terme2 + x(l)/som_ExpAlphTau(l)*(
     &                            E_ij(l,i)*expA(k,l) +
     &                              E_ij(l,k)*expA(i,l) -
     &                               2.d0*som/som_ExpAlphTau(l))
                      enddo

                      terme3 = 0.d0
                      do j = 1,nco
                         terme3 = terme3 + x(j)*( 
     &                            E_ij(i,j)*expA(k,i)/som_ExpAlphTau(i) 
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

c            d(alpha)/dx(i) = d(n.alpha) / dn(i)
             dalphadx(1:nco) = aai(1:nco)
             if(icalltempe>0 .and. icalltempe<8888) then
c               Calcul des dérivées p/r à T (pour estimation
c               de hM / cPM par exemple)
c               Calcul de d alpha / dT :
                dalpha_dT = 0.d0
                do i = 1,nco
                   dalpha_dT = dalpha_dT + asb(i)*x(i)*dlna_T(i)
                enddo
c               Pour rappel aapur = somme des x(i)*asb(i)
                dalpha_dT = dalpha_dT - aapur/t0 + dE_dT
c               at = d amix / dT = d racine(am/R²/T) / dT (uSI)
                at = dalpha_dT*bmix

                if(icalltempe >= 2) then
c                  Calcul de d² alpha / dT² :
                   d2alpha_dT2 = 0.d0
                   do i = 1,nco
                      d2alpha_dT2 = d2alpha_dT2 + asb(i)*x(i)*d2lna_T(i)
                   enddo
                   d2alpha_dT2 = d2alpha_dT2 - 2.d0/t0*
     &                           (dalpha_dT - dE_dT) + d2E_dT2
c                  att = d at / dT (uSI) :
                   att = d2alpha_dT2*bmix
                endif
             endif             
       END SELECT
c
c       read(*,*)
       end subroutine regles_melange
c
c
c**********************************************************************
c
       subroutine tempe(ideriv)
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c      toutes les fois que la valeur de la temperature est modifiee,
c      call tempe adapte les parametres a cette nouvelle valeur.
c      ideriv = 0 pour calculer fonctions de la T
c      ideriv = 1 pour calculer en sus, dérivées 1res des fonctions de la T
c      ideriv = 2 pour calculer en sus, dérivées 2ndes des fonctions de la T
c
c      Variables calculees par ce sous-programme :
c      - k12, kij(i,j)
c      - EN(i,j) = Eij/RT (Eij = coeff. d'interaction Van Laar en mol/cm3)
c      - pvap(i) (estimation Psat par Clapeyron en bar)
c      - delta(i) = [racine de a_i/(RT)] / b_i (uSI)
c      - si ideriv = 1 : on calcule aussi dk12_dT
c      - si ideriv = 2 : on calcule aussi dk12_dT et dk12_dT2
c
       implicit none
c
       include '..\source\comflash.cmn'

       integer i,j,k,l,m,mm,tnum,ideriv,mm0,nn0,intemp
       real*8 r7s3,cof,dixsrt,enij,fct,t298st
       real*8 a_function,unst,cofk,racT,racRT
       real*8 alpha1(nbgr_tot_max),alpha2(nbgr_tot_max)
       real*8 dsum,dsum2,high,rtuSI
       real*8 delta(pmax),delta1(pmax),delta2(pmax)
       real*8 f0,g0,h0,tk,sum_gr1,sum_gr2
       real*8 a_de_T,der1_a,der2_a, truc,som
       real*8 Eij,dEij_dT,d2Eij_dT2
c
       if(ideriv<0 .or. ideriv>2) then
          write(lun1,*) 'Argument tempe loufoque : ideriv = ',ideriv
          write(lun2,*) 'Argument tempe loufoque : ideriv = ',ideriv
          stop
       endif

       rt = r*t0 ! t0 en K et r = 83.14 bar.cm3/mol/K
       r7s3 = 7.0d0/3.0d0
       unst = 1.d0/t0
       racT = dsqrt(t0) ! uSI
       racRT = dsqrt(rgas*t0) ! uSI
c
       t298st = 298.15d0*unst ! uSI
       dixsrt = 10.d0/rt ! unités Péneloux

       tk = t0
c      Toutes les grandeurs sont exprimées dans le systeme
c      international
       do i = 1,2
          a_de_T = a_function(i,tk,0) ! a(T)
          der1_a = a_function(i,tk,1) ! da/dT
          der2_a = a_function(i,tk,2) ! d2a / dT2

          if(ideriv>=1) dlna_T(i) = der1_a / a_de_T
          if(ideriv>=2) d2lna_T(i) = der2_a / a_de_T

c         Calcul de ac0, ac1 et ac2 utilisées pour calculs de hM / cPM :
c         -------------------------------------------------------------
          ac0(i) = dsqrt(a_de_T / tk) / rgas ! racine(a / R²T)
          if(ideriv>=1) then
             ac1(i) = (der1_a*tk - a_de_T)/(2.d0*ac0(i)*tk**2*rgas**2) ! d ac0(i) / dT
             if(ideriv>=2) then
                ac2(i) = (-der1_a**2 * tk**2 - 2.d0*der1_a*a_de_T*tk ! d ac1(i) / dT
     &              + 3.d0*a_de_T**2 + 2.d0*a_de_T*tk**2*der2_a)
     &              / (4.d0*a_de_T*tk**3*ac0(i)*rgas**2)
             endif
          endif
       enddo
c
c calcul de asb et pvap pour les constituants purs
c
       do i=1,nco          
          pvap(i) = 10.0d0**(r7s3*(fac(i) + 1.0d0)
     &              *(1.0d0 - tc(i)*unst))*pc(i)

c         -- calcul de Racine[a(i)]/b(i) = delta(i)
          delta(i) = ac0(i) * racT / bc(i) ! uSI
c         -- calcul de Racine[a(i)/RT] / b(i) = delta_sur_RacRT(i)         
          delta_sur_RacRT(i) = delta(i) / racRT ! uSI
          delta_sur_RacRT(i) = 1.d-3 * delta_sur_RacRT(i) ! en (cm3/mol)^(-0.5)
          
c         -- calcul de delta1(i) = d(delta(i)) / dT
          if(ideriv >= 1) then 
             delta1(i) = (ac0(i) + 2.d0 * tk * ac1(i))
     &       / (2.d0*racT*bc(i)) ! uSI
          endif

c         -- calcul de delta2(i) = d(delta1(i)) / dT
          if(ideriv >= 2) then 
             delta2(i) = (tk * ac2(i) + ac1(i) - ac0(i) / (4.d0*tk))
     &               / (racT*bc(i)) ! uSI
          endif
          asb(i) = ac0(i)**2/bc(i) ! a/(bRT) du corps pur i [sans dimension]
       end do
c
c recherche pour savoir si un calcul de en a deja ete fait
c a la temperature demandee
c
c      A l'issue du paragrphe ci-dessous, tnum vaut 0 si 
c      nouvelle température ; est > 0 si T déjà en mémoire
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
c        ================================================================
c         Cas d'un gE de Van Laar :
c                 - kij constant 
c                 - kij(T) et Eij(T) (PPR78 et compagnie)
c        ================================================================
c
c        calcul des parametres EN(i,j) (coeff. d'interaction de Van Laar)
c        et des kij :
c
          kij(1,1) = 0.d0
          kij(2,2) = 0.d0

          if(k12_const>8887.d0) then ! si kij(T)
c            calcul de alpha1(j) = proportion de groupes j dans la molecule 1
c            Ce calcul est utile pour le calcul de kij(T) et de ses dérivées
c
c            -- calcul de sum_gr1 : nb de groupes total dans molecule 1
             sum_gr1 = 0.d0
             do i = 1,ngr1
                sum_gr1 = sum_gr1 + occurgr1(i)
             end do

c            -- calcul de alpha1
             alpha1 = 0.d0
             do i = 1,ngr1
                alpha1(nomgr1(i)) = occurgr1(i) / sum_gr1
             enddo
  
c            calcul de alpha2(j) = proportion de groupes j dans la molecule 2
c            -- calcul de sum_gr2 : nb de groupes total dans molecule 2
             sum_gr2 = 0.d0
             do i = 1,ngr2
                sum_gr2 = sum_gr2 + occurgr2(i)
             end do
  
c            -- calcul de alpha2
             alpha2 = 0.d0
             do i = 1,ngr2
                alpha2(nomgr2(i)) = occurgr2(i) / sum_gr2
             enddo
          endif
  
          if(tnum == 0) then ! si T pas en mémoire
c             N.B. : k12_const = 8888 si kij(T)
c             sinon, k12_const = valeur du kij constant.
             if(k12_const.lt.8887.d0) then ! ******** Cas d'un kij constant
                k12 = k12_const
                Eij = (delta(1) - delta(2))**2
     &                     + 2.d0*k12*delta(1)*delta(2)     ! uSI
                EN(1,2) = (delta_sur_RacRT(1)-delta_sur_RacRT(2))**2
     &                + 2.d0*k12*delta_sur_RacRT(1)*delta_sur_RacRT(2)
                EN(2,1) = EN(1,2)
                if(itempc.le.tempcmax) ENT(itempC,1,2) = EN(1,2)
                kij(1,2) = k12
                kij(2,1) = k12  
                if(itempc.le.tempcmax) kijT(itempC,1,2) = kij(1,2)
             else  
c               -- calcul double somme par contributions de groupes :
c               N.B. : les agr0 et bgr0 sont exprimés en Pa
c                      dsum est donc en Pa
                dsum = 0.d0
                do mm0 = 1,nbgr_used
                   do nn0 = 1,nbgr_used
                      k = liste_gr(mm0)
                      l = liste_gr(nn0)
                      if(agr0(k,l).ne.0.d0) then
                         high = bgr0(k,l) / agr0(k,l) - 1.d0
                         dsum = dsum + (alpha1(k) - alpha2(k))
     &                           * (alpha1(l) - alpha2(l))
     &                           * t298st**high
     &                           * agr0(k,l)
                         
                         
                      end if
                   enddo
                enddo
   
c             -- calcul des fonctions f0, g0 et h0 telles que :
c              f0 = E12
c              g0 = [delta(1)-delta(2)]
c              h0 = [2 * delta(1) * delta(2)]
               
                f0 = -0.5d0 * dsum ! E12 en Pa
                g0 = (delta(1)-delta(2))**2
                h0 = 2.d0 * delta(1) * delta(2)
                k12 = (f0 - g0) / h0
                kij(1,2) = k12
                kij(2,1) = k12
                if(itempc.le.tempcmax) kijT(itempC,1,2) = kij(1,2)

c               Calcul des E12/RT en (cm3/mol)^{-1}
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

c         -----------------------------------------
c         Calcul de d (E12/RT) / dT = dEijsRT(1,2)
c         Si ideriv = 2 : on calcule de plus :
c            d² (E12/RT) / dT² = d2EijsRT(1,2)
c         -----------------------------------------

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
                             truc = (alpha1(k) - alpha2(k))
     &                               * (alpha1(l) - alpha2(l))
     &                               * t298st**high
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
                 dEij_dT = 2.d0*(som-(1.0d0 - k12)*(
     &                delta(1)*delta1(2)+delta(2)*delta1(1)))
                 dEijsRT(1,2) = (dEij_dT - Eij/t0) / rtuSI ! uSI
                 dEijsRT(2,1) = dEijsRT(1,2)
                  
                 if(ideriv >= 2) then
                    som = 0.d0
                    do i = 1,nco
                       som = som + delta1(i)**2 + delta(i)*delta2(i) 
                    enddo
                    truc = delta(1)*delta2(2) + 2.d0*delta1(1)*delta1(2)
     &                    + delta(2)*delta2(1)               
                    d2Eij_dT2 = 2.0d0*(som - (1.0d0 -k12)*truc)
                    d2EijsRT(1,2) = (d2Eij_dT2 - 2.d0*dEij_dT/t0 ! uSI
     &                           + 2.d0*Eij/t0**2) / rtuSI 
                    d2EijsRT(2,1) = d2EijsRT(1,2)             
                 endif                 
             endif

        
          endif ! fin ideriv >= 1
       elseif(modele_gE == "HV (Uniquac-Wilson)".or.
     &         modele_gE == "HV (NRTL)") then
              
          if(modele_gE == "HV (Uniquac-Wilson)")then
            alphgE = 1.0d00
          end if
          if(tnum == 0) then ! si T pas en mémoire            
c            Parametres fonction de T :
             DO i = 1,pmax
                pargE(i,i) = 0.d0
                expA(i,i) = 1.d0
                DO j = 1,pmax
                  IF(i.NE.j) THEN
                      pargE(i,j) = agE0(i,j)/t0 + bgE0(i,j)
     &                             + cgE0(i,j)*t0
                      expA(i,j) = DEXP(-alphgE(i,j) * pargE(i,j))
                      if(itempc.le.tempcmax) then
                         pargEt(itempc,i,j) = pargE(i,j) 
                         expAt(itempc,i,j) = expA(i,j)
                      end if
                   ENDIF
                ENDDO
             ENDDO
          else ! temperature déjà considérée (effet mémoire)
             DO i = 1,pmax
                pargE(i,i) = 0.d0
                expA(i,i) = 1.d0
                DO j = 1,pmax
                  IF(i.NE.j) THEN
                      pargE(i,j) = pargEt(tnum,i,j)
                      expA(i,j) = expAt(tnum,i,j)
                   ENDIF
                ENDDO
             ENDDO
        endif ! fin boucle tnum = 0

          if(ideriv >= 1) then
c            d / dT :
             DO i = 1,pmax
                dpargE(i,i) = 0.d0
                DO j = 1,pmax
                  IF(i.NE.j) THEN
                     dpargE(i,j) = -agE0(i,j)/t0**2 + cgE0(i,j)
                   ENDIF
                ENDDO
             ENDDO
          endif

          if(ideriv >= 2) then
c            d^2 / dT^2 :
             DO i = 1,pmax
                dpargE(i,i) = 0.d0
                DO j = 1,pmax
                  IF(i.NE.j) THEN
                      d2pargE(i,j) = 2.d0*agE0(i,j)/t0**3
                   ENDIF
                ENDDO
             ENDDO
          endif
       endif ! fin choix modele_gE
C
       end subroutine tempe
c
c**********************************************************************
c
       subroutine eqet0
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c      Calcul des constantes universelles des EoS cubiques :
c      - etac = b/vc (compacite critique), sans unite
c      - omegab et omegaa (cstes univ. intervenant dans a et b) 
c      - b(i) = covolume du constituant i en (cm3/mol)
c      - ac(i) en bar.cm6/mol2
c      - msoave(i) = facteur de forme de Soave (sans unite)

       implicit none
c
       include '..\source\comflash.cmn'
c
c*** variables internes
       integer i
       real*8 etac,zc
       real*8 rtsp
c
c
c---------------------------------------------------
c
c  1) calcul des parametres independants des constituants
c
       etac = ((1.d0-r1)*(1.d0-r2)**2)**(1.d0/3.d0) +
     &        ((1.d0-r2)*(1.d0-r1)**2)**(1.d0/3.d0) + 1.d0
       etac = 1.d0 / etac

       zc = 1.d0/(3.d0 - etac*(1.d0+r1+r2))

       omegab = zc*etac

       omegaa = (1.d0-etac*r1) * (1.d0-etac*r2) * (2.d0-etac*(r1+r2))
     &          / ((1.d0-etac)*(3.d0-etac*(1.d0+r1+r2))**2)
c
c
c  2) calcul des parametres qui dependent de la nature des constituants
c
c     - calcul de ac, b, msoave
c
       do i=1,nco
          rtsp = r*tc(i)/pc(i) ! pc en bar, tc en K et r = 83.14 bar.cm3/mol/K
          b(i) = omegab*rtsp ! cm3/mol
          ac(i) = omegaa*rtsp*r*tc(i) ! bar.cm6/mol2
          ac(i) = ac(i)*1.d-7 ! unité du systeme international (Pa.m6/mol2)

c         Calcul de bc(i) = b(i)/R en uSI = K/Pa
c          r = 83.14 bar.cm3/mol/K et b en cm3/mol      
          bc(i) = 1.d-5 * b(i) / r ! uSI
          
          if(FoncAlpha.EQ.'SOAVE') then
c         -------------------------
c         Fonction alpha de SOAVE :
c         -------------------------
c         EoS de PR78 :
c         -----------
            if(fac(i).le.0.491d+00) then
               msoave(i) = 0.37464d+00 + 1.54226d+00*fac(i)
     &                     - 0.26992d+00*fac(i)**2
            else
               msoave(i) = 0.379642d+00 + 1.48503d+00*fac(i)
     &                     - 0.164423d+00*fac(i)**2
     &                     + 0.016666d+00*fac(i)**3
            endif
c         -------------------------
c         -------------------------              
          endif 
         
       end do
c
c      Modèle Wilson / Uniquac
       IF(choix_ge == 1) THEN ! Wilson        
          mu(1:pmax) = 1.d0 
c          sigma(1:pmax) = b(1:pmax) ! cm3/mol modif tempo APM 03/03/20
          sigma(1:pmax) = b(1:pmax) - parC(1:pmax) ! cm3/mol           
       ELSEIF(choix_ge == 2) THEN !UNIQUAC
          mu(1:pmax) = qi(1:pmax) !q UNIQUAC parameter
          sigma(1:pmax) = qi(1:pmax) !q UNIQUAC parameter        
       ENDIF

       end subroutine
c
c
c**********************************************************************
c
       function a_function(numero,tk,ideriv)
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c
c Calcul du parametre a(T) dans le terme attractif de l'EoS cubique
c des corps purs
c numero = entier utilise pour designer le corps pur considere
c Unite de a(T) : Pa.m6/mol2 (systeme international)
c ideriv : indice de dérivation :
c           - si ideriv = 0 : a_function renvoie a/(RT)
c           - si ideriv = 1 : a_function renvoie d[a/(RT)] / dT
c           - si ideriv = 2 : a_function renvoie a/(RT)
c
c ATTENTION : - la fonction a(T) est egalement definie dans la sous-routine
c             parametres_EoS (pour le calcul de Psat et Teb des corps purs)
c             - la fonction a(T) est egalement definie dans la sous-routine
c             temset (pour le calcul des hM et cPM)
c
       implicit none
c
       include '..\source\comflash.cmn'
c
       integer numero,ideriv
       integer pcrit
       real*8 a_function,tk
       real*8 alpha,rac_alpha
       real*8 delt,gamm,theta
       real*8 polTH

       if(ideriv==0) then ! calcul de la fonction a(T)
            if(FoncAlpha.EQ.'SOAVE') then
c         -------------------------
c         Fonction alpha de SOAVE :
c         -------------------------
                rac_alpha = 1.d0+msoave(numero)*(1.d0-
     &                                            DSQRT(tk/tc(numero)))
                alpha = rac_alpha**2
c         -------------------------
c         -------------------------              
            elseif(FoncAlpha.EQ.'TWU91'.OR.FoncAlpha.EQ.'TWU88') then 
c         modif tempo APM
c         -------------------------
c         Fonction alpha de TWU:
c         -------------------------
                delt = parN(numero) * (parM(numero) - 1.d0)
                gamm = parN(numero) * parM(numero)
                alpha = (tk/tc(numero))**delt * DEXP(parL(numero)*(1.d0
     &                                         - (tk/tc(numero))**gamm))  
c         -------------------------
c         -------------------------             
            endif
            a_function = ac(numero)*alpha ! Pa.m6/mol2            
       elseif(ideriv == 1) then ! calcul de d a(T) / dT
            if(FoncAlpha.EQ.'SOAVE') then
c         -------------------------
c         Fonction alpha de SOAVE :
c         -------------------------       
                rac_alpha = 1.d0+msoave(numero)*(1.d0-
     &                                            DSQRT(tk/tc(numero)))
                pcrit = pc(numero) * 1.d5 ! bar > Pa
c         da / dT dans les unites du systeme international          
                a_function = -msoave(numero)*ac(numero)*rac_alpha
     &                      /(dsqrt(tc(numero)*tk))
c         -------------------------
c         -------------------------              
            elseif(FoncAlpha.EQ.'TWU91'.OR.FoncAlpha.EQ.'TWU88') then 
c         modif tempo APM
c         -------------------------
c         Fonction alpha de TWU:
c         -------------------------     
                delt = parN(numero) * (parM(numero) - 1.d0)
                gamm = parN(numero) * parM(numero)
                alpha = (tk/tc(numero))**delt * DEXP(parL(numero)*(1.d0
     &                                         - (tk/tc(numero))**gamm))     
                a_function = ac(numero) * alpha * tc(numero) / tk *
     &          (delt - parL(numero) * gamm * (tk/tc(numero))**gamm)
                a_function = a_function / tc(numero)               
            endif     
       elseif(ideriv == 2) then ! calcul de d2 a(T) / dT2
            if(FoncAlpha.EQ.'SOAVE') then
c         -------------------------
c         Fonction alpha de SOAVE :
c         -------------------------        
                 a_function = msoave(numero)*(1.d0+msoave(numero))*
     &                  ac(numero)/(2.d0* tk**1.5d0 * dsqrt(tc(numero)))          
            elseif(FoncAlpha.EQ.'TWU91'.OR.FoncAlpha.EQ.'TWU88') then 
c         modif tempo APM
c         -------------------------
c         Fonction alpha de TWU:
c         -------------------------         
                delt = parN(numero) * (parM(numero) - 1.d0)
                gamm = parN(numero) * parM(numero)
                alpha = (tk/tc(numero))**delt * DEXP(parL(numero)*(1.d0
     &                                         - (tk/tc(numero))**gamm))               
                theta = parL(numero) * gamm * (tk/tc(numero))**gamm
                polTH = theta**2.d0 - (2.d0*delt + gamm - 1.d0)*theta +
     &                                              delt * (delt - 1.d0)           
                a_function = ac(numero) * alpha * polTH * (tc(numero) /
     &                                                     tk)**2
                a_function = a_function / tc(numero)**2
             endif            
       else
          write(lun1,*) 'ideriv non defini (FUNCTION a_function)'
          write(lun2,*) 'ideriv non defini (FUNCTION a_function)'
          stop
       endif
       
       end function a_function
c
c**********************************************************************
c
       subroutine rop(j,p00)
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c
c  pilote de la resolution de l'equation d'etat
c  call rop(j,p00) provoque la resolution, a temperature constante, de
c  l'equation peq(eta=b/v) = p00 par la methode de newton (paragraphe
c  ii-ii.1, page 41 de la these de A. Gramajo)
c  - p00 est la pression en bar
c  - la temperature a laquelle l'EoS est resolue est t0 (K)
c  - la composition du fluide est x(:) (fractions molaires)
c  - l'initialisation de eta depend de la valeur de j
c    j = 1 : initialisation correspondant au liquide
c    j = 2 : initialisation correspondant au gaz
c    j = 3 : initialisation a partir de la valeur de eta se
c            trouvant deja en memoire (dans com2 du fichier COMFLASH.CMN).
c
       implicit none
c
       include '..\source\comflash.cmn'
c
c*** liste d'appel
       integer j,nech
       parameter(nech=4)
       real*8 p00
c
c*** variables internes
       integer nit,ij,j_in,kk
       real*8 eta0,P00_in,P00_inn
       logical echec(nech)
c
c
c---------------------------------------------------
c
c
       j_in = j
       echec(1:nech) = .FALSE.
       ij = 0
c
 10    nit = 1
       goto(1,2,3), j
c
   1   continue
       if(.NOT.echec(1)) then
          eta = 0.999d0
       elseif(.NOT.echec(2)) then
          eta = 0.9999d0
       elseif(.NOT.echec(3)) then
          eta = 0.99999d0
       else
          eta = 0.999999d0
       endif

       goto 3
c
   2   continue
       if(.NOT.echec(1)) then
          eta = 1.d-3
       elseif(.NOT.echec(2)) then
          eta = 1.d-4
       elseif(.NOT.echec(3)) then
          eta = 1.d-5
       else
          eta = 1.d-6
       endif
c
   3   continue
       call eqet
       if(nit.eq.1) goto 4
c
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

                write(lun1,141) nom(1)(1:len_trim(nom(1))),
     +                          nom(2)(1:len_trim(nom(2)))
                write(lun1,*) 'P/bar = ',p00
                write(lun1,*) 'T/K   = ',t0
                write(lun1,*) 'kij = ',k12
                write(lun2,141) nom(1)(1:len_trim(nom(1))),
     +                          nom(2)(1:len_trim(nom(2)))
                write(lun2,*) 'P/bar = ',p00
                write(lun2,*) 'T/K   = ',t0
                write(lun2,*) 'kij = ',k12
                
                
 141            format(1x,'EOS NON SOLUBLE - systeme ',a,' - ',a)
c                stop
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
c
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
             elseif(p00>1000.d0 .AND. p00<=10000.d0 .AND.
     &              dabs(p00 - peq).ge.1.0d-8) then
                nit = 49
             else
                return     !   sortie favorable
             endif
          end if
       end if
c
   4   eta0 = eta
       eta = eta + (p00 - peq)/dpdeta
c       print*,'nit = ',nit,'eta = ',eta,'z = ',z,
c     &     'f = ',dabs(eta0-eta)       
c      if(eta.gt.1.00d0.and.j.eq.2) eta = 0.10d0
       nit = nit + 1
       goto 3
c
       end subroutine rop
c
c**********************************************************************
c
       subroutine eqet
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c
c  utilisation de l'equation cubique generale
c
       implicit none
c
       include '..\source\comflash.cmn'
c
c
c      Forme de l'EoS cubique consideree :
c      -----------------------------------
c      z = 1/(1-eta) - alpha * eta * qpseta [avec alpha = a/(bRT)]
c
c
       if(eta.eq.1.0D+00) eta = 1.0d0 - 1.0d-10
c
c      Terme repulsif de z = v/(v-b) = 1/(1-eta) avec eta = b/v
       y = 1.0d0/(1.0d0 - eta)

c      Fraction rationnelle en v de la partie attractive de z :
       qpseta = 1.d0 / ((1.d0 - r1*eta)*(1.d0 - r2*eta))

c      aa = a/(RTb) = alpha
c      z = v/(v-b) - a.v/[RT.(1-r1.v).(1-r2.v)]
c        = 1/(1-eta) - alpha.eta/[(1-r1.eta)(1-r2.eta)]
       z = y - aa*eta*qpseta

c      Arrangement douteux pour eviter probleme dans calcul de LN Phi
c       z = dabs(z)
       if(.not.isGPED) z = dabs(z)

c      dzdeta = d(z) / d(eta) a T et compo fixees
       dzdeta = y**2 - aa*(1.d0-r1*r2*eta**2)*qpseta**2

       rtsb = rt/bm
       peq = z*rtsb*eta
       dpdeta = (dzdeta*eta + z)*rtsb
       return
       end
c
c**********************************************************************
c
      subroutine fug
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c     Calcul des coeff. de fugacite par EoS cubique
c     La subroutine EQET doit avoir ete prealablement appelee
c     (pour le calcul de z, y, etc.)
c
       implicit none
       include '..\source\comflash.cmn'
       integer i
       real*8 lphi0
c
       v = bm/eta
       q = 1.d0/(r1-r2)*DLOG((1.d0-r2*eta)/(1.d0-r1*eta))
c
c      y = 1 / (1 - eta) = v / (v - bm)
      if(.not.isGPED) then
          lphi0 = dlog(y/z) ! -Ln(z) - Ln[(v-bm)/v]
       else
          lphi0 = dlog(y) + dlog(rtsb*eta) ! Ln(RT/v) - Ln[(v-bm)/v]
       end if
       do i=1,nco
          lphi(i) = lphi0 + (z-1.d0)*derbm(i)/bm
     &              -q*aai(i)   
       end do
c      print*,'phi = ',dexp(lphi(1:2))
       end
c
c**********************************************************************
c
       subroutine rcub(i3r,eta1)
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c
c  connaissant une racine eta de l'equation du 3eme degre 
c  on calcule celle, eta1, qui ne correspond pas au systeme instable
c
       implicit none
c
       include '..\source\comflash.cmn'
c
c*** liste d'appel
c
       integer i3r
       real*8 eta1,coefa,coefb,coefc
c
c*** variables internes
c
       real*8 pbsrt,bbb,ccc,del,eta2
c
c
c--------------------------------------------------
c
c
       pbsrt = peq/rtsb

c      L'EoS cubique s'ecrit : 
c      coefa.eta^3 + coefb.eta^2 + coefc.eta -pbsrt = 0
c
       coefa = pbsrt*r1*r2 + r1*r2 + aa ! aa = +alpha
       coefb = -pbsrt*(r1*r2 + r1 + r2) - r1 - r2 - aa
       coefc = pbsrt*(r1 + r2 + 1.d0) + 1.d0

c      on re-ecrit l'EoS cubique sous la forme :
c      (eta - eta connu).(coefa.eta^2 + bbb.eta + ccc) = 0
c
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
c                        ! de la racine correspondant au systeme instable
          if(eta1.lt.0.0d0.or.eta1.gt.1.0d0) i3r = 0
       endif
       return
       end
c
c**********************************************************************
c
       subroutine dertp(ider,x)
c      -------------------
c      SUBROUTINE NETTOYEE
c      -------------------
c  calcul des derivees de lphi en variables t,p,n  (ider=1)
c  ou de lphi et p en variables t,v,n (ider=2)
c
       implicit none
c
       include '..\source\comflash.cmn'
c
       integer ider
       integer i,k
       real*8 x(pmax)
       real*8 dzdesz,qp,ertsb,detadp,somd,etasb,dpdvseta
       real*8 somdp,dlphidvseta,dzdt,dvdt
       real*8 dzdx(pmax),dpdx(pmax)
       real*8 dlphidx(pmax,pmax)

c
       dzdesz = dzdeta/z ! d ln z / d eta
       qp = eta*qpseta   ! eta / [(1-r1.eta)(1-r2.eta)]
       ertsb = eta*rtsb
       detadp = 1.0d0/dpdeta ! d eta / dP a T et compo fixees
c
c
c  derivees en variables t,eta,x
c  =============================
c
       do i=1,nco
c         qpseta = 1 / [(1-r1.eta)(1-r2.eta)]
c         y = 1 / (1-eta)
c         derbm(i) = d(n.bm)/dn(i)
          dlphideta(i) = -dzdesz + y +derbm(i)/bm*dzdeta
     &                   -aai(i)*qpseta

          dzdx(i) = -dalphadx(i)*qp    ! d(z)/dx(i)
          dpdx(i) = ertsb*(dzdx(i) - z*dbmdx(i)/bm)
       end do

c
       do i=1,nco
          do k=1,nco
             dlphidx(i,k) = - dzdx(k)/z
     &                      + derbm(i)/bm*dzdx(k)
     &                      - derbm(i)*dbmdx(k)*(z - 1.d0)/bm**2
     &                      + (z - 1.d0)/bm*der2bm(i,k)
     &                      - daidx(i,k)*q
          end do
       end do

c
      if(ider.eq.2) goto 10
c
c  derivees en variables t,p,n      ! calcul de flash (ider = 1)
c  =============================================================
c
       etasb = eta/bm
       dpdvseta = -etasb*dpdeta ! 1/eta * (dP/dv)
       dpdv = eta*dpdvseta      ! dP/dv
       dpdt = ertsb*(- dalpha_dT * qp + z/t0) ! dP/dT
       detadt = eta*etasb*dpdt/dpdv ! deta / dT
       dzdt = dzdeta*detadt - dalpha_dT * qp ! dz /dT

       do i=1,nco
          dlphidt(i) = dlphideta(i)*detadt  -   ! d ln phi(i) / dT
     &           dalpha_dT*qp*(derbm(i)/bm - 1.0d0/z) - q*daai_dT(i)

          dlphidp(i) = dlphideta(i)*detadp
          somd = 0.0d0
          do k=1,nco
             dlphidx(i,k) = dlphidx(i,k) - dpdx(k)*dlphidp(i)
             somd = somd + x(k)*dlphidx(i,k)
          end do
          do k=1,nco
             dmudn(i,k) = dlphidx(i,k) - somd - 1.0d0
          end do
          if(x(i).lt.1.0d-300) x(i) = 1.0d-250  !traficos
          dmudn(i,i) = dmudn(i,i) + 1.0d0/x(i)
c
       end do        ! dmudn(i,k) = nj*dik (i.7, page 17), nj etant le
c                    ! nombre de moles total dans la phase j consideree
c       print*,'dlphidt = ',dlphidt(1:2)
c       read(*,*)
       return
c
c
c  derivees en variab     les t,v,n      ! (ider = 2)
c  ============================================================
c
 10    etasb = eta/bm
       dpdvseta = -etasb*dpdeta ! 1/eta * (dP/dv)
       dpdv = eta*dpdvseta      ! dP/dv
       dzdt = -dalpha_dT * qp
       dpdt = ertsb*(dzdt + z/t0) ! dP/dT       
c
       somdp = 0.d0
       do i=1,nco
          somd = 0.d0
          dlphidvseta = -etasb*dlphideta(i) ! 1/eta * d ln phi(i) / dv
          dmudv(i) = eta*dlphidvseta        ! d ln phi(i) / dv
          dmudt(i) = (derbm(i)/bm - 1.0d0/z)*dzdt  - q*daai_dT(i) ! d ln phi(i) / dT
          do k=1,nco
             dlphidx(i,k) = dlphidx(i,k) - dlphidvseta*dbmdx(k)
             somd = somd + x(k)*dlphidx(i,k)
          end do
          do k=1,nco
            dmudn(i,k) = dlphidx(i,k) - somd -1.0d0 + eta*dlphideta(i)           
          end do
c          dmudn(i,i) = dmudn(i,i) + 1.0d0/x(i)
          dpdx(i) = dpdx(i) - dpdvseta*dbmdx(i)
          somdp = somdp + x(i)*dpdx(i)
       end do
c 
       do i=1,nco
          dpdn(i) = eta*dpdeta + dpdx(i) - somdp
c
       end do        ! dpdn(i)    = nj*(dp/dnk),
       do i=1,nco
          do k=1,nco
               dmudn(i,k) = dmudn(i,k) + dpdn(k)/peq     
          enddo
       enddo    
       return        ! dmudn(i,k) = nj*(dmui/dnk),  nj etant le
c                    ! nombre de moles total dans la phase j consideree

       end
c
c**********************************************************************
c
       subroutine chol(n,m,d,f,t,x,isym,itest,irtn)
c
c  isym = 1 la matrice  d  est symetrique, on definit  t=d
c       = 0 la matrice  d  n'est pas symetrique, on calcule  t=dt*d
c  irtn = 0 on resout le systeme  f=d*x (on calcule t=f ou t=dt*f)
c       = 1 on triangule uniquement la matrice t (calcul du determinant)
c       = 2 on inverse la matrice triangulaire t (a la sortie t est la
c           matrice inverse)
c
c
       implicit none
c
       include '..\source\comflash.cmn'
c
c*** variables de la liste d'appel
c
       integer n,m,itest,irtn,isym
       real*8 f(1),t(1),x(1),s
       real*8 d(n0,m0)
c
c*** variables internes
c
       integer i,ii,im1,j,jj,jm1,jp1,k,k1,kk,l,mf,mm1,mp1
       real*8 t1(100)
       logical negatif
c
c
c-------------------------------------------------
c
c
       itest = 0
       k = 0
c
c  - cas d'une matrice  d  symetrique  (isym=1)
c
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
c
c  - cas d'une matrice  d  non symetrique  (isym=0)
c
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
c
c
c-------------------------------------------------
c
c
c  triangulation de la matrice  t
c
c
 200   mp1 = m + 1
       k = 0
       do i=1,mp1
          ii = i*(i-1)/2
          do j=1,i
             if(j.eq.mp1) goto 29
c
             k = k + 1
             jm1 = j - 1
             s = 0.0d0
             jj = j*(j-1)/2
             if(jm1.eq.0) goto 22
c
             do l=1,jm1
                s = s + t(jj+l)*t(ii+l)
             end do
c
  22         s = t(k) - s
             if(i.eq.j) goto 25
c
             t(k) = s/t(jj+j)
             goto 29
c
  25         if(s.gt.0.0d0) goto 26
c
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
c
  26         t(k) = dsqrt(s)
c
  29      end do
       end do
c
       if(irtn.eq.1) return
c
c
c-------------------------------------------------
c
c
c  inversion de la matrice  t
c
c
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
c
          mf = m*(m+1)/2
          do i=1,mf
             t(i) = t1(i)
          end do
          return
       endif
c
c
c-------------------------------------------------
c
c
c  resolution du systeme  f=t*x
c
c

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
c
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
             if((dlog(dabs((t(k-l)-s)))-dlog(dabs(t(jj)))).gt.709.0D+00)
     &       then
                x(j) = 1.0D+250
                if(negatif) x(j) = -x(j)
             else
                x(j) = (t(k-l)-s)/t(jj)
             endif
          endif
       end do
c
       return
       end
c
c**********************************************************************
c
       subroutine lecflash
c
c      lecture du fichier parflash qui determine l'eos et pilote le flash
c      Si le modele est (E-)PPR78 : lecture du fichier "tabingr.*" qui donne 
c      les Akl et Bkl de PPR78 ou E-PPR78
c      choix_eos sert a savoir dans quel fichier on lit
c
       implicit none
c
       include '..\source\comflash.cmn'
c
c*** variables internes:
c
       character*60 fich,texte
       real*8 coeff
       integer i,j
       integer nsyst
       integer id1,id2
       real*8 parAA,parBB       
c
c      fichier ouvert
c
       open(2,file='..\source\parflash')
c
       read(2,*) 
       read(2,*) texte,ph3
       read(2,*) texte,nimp
c
       read(2,*) texte,choix_r1r2
c
       if(choix_r1r2.eq.1) then ! VDW
          r1 = 0.d0
          r2 = 0.d0
          write(lun1,*) "ATTENTION : Fichier TABINGR "//
     &      "non compatible avec valeurs de r1 & r2 selectionnees"
          write(lun2,*) "ATTENTION : Fichier TABINGR "//
     &      "non compatible avec valeurs de r1 & r2 selectionnees"
          pause
       elseif(choix_r1r2.eq.2) then ! SRK
          r1 = -1.d0
          r2 = 0.d0
          write(lun1,*) "ATTENTION : Fichier TABINGR "//
     &      "non compatible avec valeurs de r1 & r2 selectionnees"
          write(lun2,*) "ATTENTION : Fichier TABINGR "//
     &      "non compatible avec valeurs de r1 & r2 selectionnees"
          pause
       elseif(choix_r1r2.eq.3) then ! Peng-Robinson
          r1 = -1.d0 - DSQRT(2.d0)
          r2 = -1.d0 + DSQRT(2.d0)
       else
          write(lun1,*) "Valeur de choix_r1r2 non geree "
          write(lun1,*) "(fichier PARFLASH)"
          write(lun2,*) "Valeur de choix_r1r2 non geree "
          write(lun2,*) "(fichier PARFLASH)"
          stop
       endif

       read(2,*) texte,iac
       read(2,*) texte,isubs
       read(2,*) texte,cam
       read(2,*) texte,epseq
       read(2,*) texte,eps0
       read(2,*) texte,eps1
c
       if(isubs.eq.1.and.eps1.gt.1d-6) eps1 = 1d-6
c
       read(2,*) texte,oa
       read(2,*) texte,nb
       read(2,*) texte,l0
       read(2,*) texte,nlp,nls
       read(2,*) texte,(ld(i),i=1,3)
       read(2,*) texte,nitflash
c
       close(2)
       
       if(choix_eos==1) then ! PPR78/E-PPR78
c
c      lecture du fichier "tabingr.*" qui donne les Akl et Bkl de
c      PPR78/E-PPR78
c
         fich = '..\source\'//tabingr_file(1:LEN_TRIM(tabingr_file))
         open(14,file=fich,status='old')
         read(14,*) nbgr_tot
         if(nbgr_tot.gt.nbgr_tot_max) then
            write(lun1,*) 'Augmenter nbgr_tot_max dans COMFLASH'
            write(lun1,*) 'nbgr_tot = ',nbgr_tot,
     &                    'nbgr_tot_max = ',nbgr_tot_max
            write(lun2,*) 'Augmenter nbgr_tot_max dans COMFLASH'
            write(lun2,*) 'nbgr_tot = ',nbgr_tot,
     &                    'nbgr_tot_max = ',nbgr_tot_max
            STOP
         endif

c        Pour de-tilder :            
         coeff = 1.5d0 + DSQRT(2.d0)

         do i=1,nbgr_tot
            read(14,*) (agr(i,j),j=1,i) ! en MPa
            do j=1,i
c              On de-tilde Akl :
               agr(i,j) = agr(i,j)/coeff ! en MPa
               agr(j,i) = agr(i,j) ! en MPa
            end do
         end do
         do i=1,nbgr_tot
            read(14,*) (bgr(i,j),j=1,i) ! en MPa
            do j=1,i
c              On de-tilde Bkl :
               bgr(i,j) = bgr(i,j)/coeff ! en MPa
               bgr(j,i) = bgr(i,j) ! en MPa
            end do
         end do
         close(14)

         if(tabingr_file(1:LEN_TRIM(tabingr_file))=='TABINGR.YLG')then
         
          open(14,file='database_names.d',status='old')
          read(14,*) nsyst
          read(14,*) texte
          do i = 1,nsyst
              read(14,*) id1,id2,parAA,parBB
                agr(id1,id2) = parAA
                bgr(id1,id2) = parBB

                agr(id2,id1) = agr(id1,id2)
                bgr(id2,id1) = bgr(id1,id2)              
          end do
          read(14,*) texte
          close(14)         
         
         end if
       endif
       return
       end
c
c---------------------------------------------------------------------------
c
      subroutine tpd_rp(z1,t,p00)
      implicit none
      include '..\source\comflash.cmn'

      real*8 eta0,v0,g0
      real*8 x0(pmax),lnphi0(pmax),lnphi(pmax),x1,tpd,z1,x2,z2
      real*8 p00,t,d1,d2,eps,x1old,eps_tpd,y1,stepx,n10,n20,n1,n2
      real*8 tempo,xxmin,xxmax,eps_tpd2
      integer j0,ider,istb,iter,niter_max,compteur
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
      lnphi0(1:pmax) = lphi(1:pmax)
      d1 = dlog(z1) + lnphi0(1)
      d2 = dlog(z2) + lnphi0(2)
      
c 1ere etape : on fait un balayage en compo du TPD
c ================================================

      xxmin = 1.d-9
      xxmax = 1.d0 - xxmin
      stepx = (xxmax - xxmin) / 150.d0
      compteur = 0
      do x1 = xxmin,xxmax,stepx
         compteur = compteur + 1
         x0(1) = x1
         x0(2) = 1.d0 - x1
         x2 = 1.d0 - x1 
         call phas1(j0,ider,istb,p00,x0,eta0,v0,g0)
         lnphi(1:pmax) = lphi(1:pmax)
         tpd = x1 * (lnphi(1)-lnphi0(1) + dlog(x1/z1)) +
     &         x2 * (lnphi(2)-lnphi0(2) + dlog(x2/z2))
         
         if(tpd.lt.eps_tpd) then
c           Si la composition nonstable est environ egale a la composition
c           de la phase testee, on se mefie (par exemple, si 
c           z1 et x1 sont distants de 0.01 et si tdp = -2.d-10, on est
c           proche d'un min qui devrait valoir zero mais qui numeriquement
c           peut etre legerement < 0)
            if(DABS(x1-z1).lt.5.d-3 .AND. tpd.gt.eps_tpd2) then
c              on estime que le test de stabilite en ce point n'est pas
c              concluant
               goto 400 ! on ne declare pas la compo non-stable et on
c                         poursuit l'analyse du TPD
            endif
            stable = .false. 
            goto 300 ! return
         end if
 400     continue
      enddo

      goto 300 ! EXIT / FIN DES TESTS **************************
c      if(z1<0.999d0 .and. z1>0.001d0) goto 300 ! EXIT / FIN DES TESTS **************************

c 2eme etape : Si stable, on cherche precisement a localiser les minima
c =====================================================================

c     a) Initialisation a x1 proche de 0
c     ----------------------------------

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

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) +
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

c     b) Initialisation a x1 proche de 1
c     ----------------------------------

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

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) +
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

c     c) Initialisation a x1 tres proche de 1
c     ---------------------------------------

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

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) +
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

      goto 300 ! EXIT / FIN DES TESTS **************************


c     d) Initialisation de x1 entre xp(1,1) et xp(1,2)
c     ------------------------------------------------

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

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) +
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif


c     e) Initialisation de x1 a diverses compo entre xp(1,1) et xp(1,2)
c     -----------------------------------------------------------------

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

      tpd = x1 * (lphi(1)-lnphi0(1) + dlog(x1/z1)) +
     &      x2 * (lphi(2)-lnphi0(2) + dlog(x2/z2))
      if(tpd.lt.eps_tpd .and. dabs(x1 - z1).gt.1.d-8) then
         stable = .false.
         goto 300
      endif

      goto 987

 300  continue
      tpdminimum = tpd ! modif rp tempo

      end
c
c---------------------------------------------------------------------------
c
      subroutine search_min_tpd(z1,t,p00,xmin)

c     Recherche des points stationnaires du TPD
c     Si TPDmin = 0, on cherche la solution qui ne correspond pas a z1

      implicit none
      include '..\source\comflash.cmn'

      real*8 eta0,v0,g0,xstat_up(20),xstat_dn(20),xstat(20)
      real*8 x0(pmax),lnphi0(pmax),lnphi(pmax),x1,tpd,z1,x2,z2
      real*8 p00,t,d1,d2,eps,x1old,eps_tpd,y1,stepx,n10,n20,n1,n2
      real*8 tempo,xmin,xmax,dtpd,tpd1,tpd2,dtpdold,tpdold,res
      real*8 tpdstat(20),tpdmin,ln_fi(pmax),tpd01,tpd02
      integer j0,ider,istb,iter,niter_max,err,i,k,nstat
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
      lnphi0(1:pmax) = lphi(1:pmax)
      d1 = dlog(z1) + lnphi0(1)
      d2 = dlog(z2) + lnphi0(2)

      xmin = 1.d-12
      xmax = 1.d0 - xmin
      stepx = (xmax - xmin) / 5000.d0
      dtpd = 0.d0
      tpd = 0.d0
      k = 1

c     Resolution de l'equation (dTPD / dx1) = 0 (points stationnaires de TPD)
c     xstat_dn(i) est une approximation par valeurs inf de la i-eme solution
c     xstat_up(i) est une approximation par valeurs sup de la i-eme solution

      do x1 = xmin,xmax,stepx
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

c     Une fois les points stationnaires approximes, on les calcule precisement 
c     par dichotomie

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

c     A present, on calcule les min
      tpdmin  = 1.d31

      do i = 1,nstat
c        On prend la solution autre que z1
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
 420  continue

      end
c
c---------------------------------------------------------------------------
c
      subroutine calcul_teb(numero,p_in,teb,err0)

c Calcul de la temperature d'ebullition d'un corps pur
c par un EoS cubique.
c
c INPUT : numero, entier = 1 ou 2 pour spécifier le constituant
c         p_in (pression specifiee) en bar
c OUTPUT :  teb (K), err0 (variable d'erreur, = 0 si tout va bien)
c
      implicit none
      include '..\source\comflash.cmn'
      real*8 apr,bpr,pcrit,tcrit
      real*8 teb,teb0,p_in,vol_mol_liq,vol_mol_vap
      real*8 lnphi_liq,lnphi_vap,lnphi,h_res,h_res_liq,h_res_vap
      real*8 pres0,pliq,pEoS,a_function
      real*8 res1,res2,res3,cur,vvap_old,vliq_old,eps

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
c     procedure d'initialisation teb :
c     ------------------------------
      call initialisation_teb(numero,teb,pres0,err0)
      if(err0.ne.0) goto 31

c     Methode de Newton pour calculer teb :
c     -----------------------------------
 30   continue
c     a_function donne a(T) en Pa.m6/mol2
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

      teb = teb0 + rgas * teb**2 * (lnphi_vap - lnphi_liq) 
     &             / (h_res_vap - h_res_liq)
     
      res1 = dabs(lnphi_vap - lnphi_liq)
      res2 = dabs(teb - teb0)/dabs(teb0)
      res3 = 0d0
      if(nb_iter>1) then
         res3 = dabs(vol_mol_liq - vliq_old)/dabs(vliq_old)
     &          + dabs(vol_mol_vap - vvap_old)/dabs(vvap_old)
      endif

      if(res1.lt.eps .and. res2.lt.eps .and. res3.lt.eps
     &   .and. nb_iter>1) then
         goto 31
      else
         goto 30
      endif

 31   continue

      end subroutine calcul_teb
c
c ------------------------------------------------------------------
c
      subroutine initialisation_teb(numero,teb,pres0,err0)

c Initialisation de la temperature d'ebullition  : 
c   1) initialisation a partir d'une formule l'approximant connaissant 
c      tc, pc et omega.
c   2) Affinage de cette valeur : on calcule teb de sorte que la 
c      resolution de l'equation d'etat conduise a une racine liquide et 
c      une racine vapeur

      implicit none
      include '..\source\comflash.cmn'
      real*8 tcrit,pcrit,omega0
      real*8 apr,bpr,a_function
      real*8 pres0,teb,vol_mol_liq,vol_mol_vap,tinf,tsup
      integer nb0,nb_iter,nb_iter_max,err1,err0,numero

      common/ab_pr/apr,bpr

      pcrit = pc(numero)*1.d5 ! bar>Pa
      tcrit = tc(numero)
      omega0 = fac(numero)

      nb_iter = 0
      nb_iter_max = 1000

      teb = tcrit*(omega0 + 1.d0) / (omega0 + 1.d0 - 3.d0/7.d0 
     &                             * dlog10(pres0 / pcrit))

c     La resolution de l'EoS a T=teb et P=pres0 amene-t-elle 3 racines ?
 50   continue
c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,teb,0)
      call resol_eos(teb,pres0,vol_mol_liq,vol_mol_vap,nb0,err0)
      if(err0.ne.0) goto 52

      if(nb_iter.eq.nb_iter_max) then
         err0 = 2
         goto 52
      endif
      nb_iter = nb_iter + 1

c     si la racine vapeur a ete trouvee mais pas la racine liquide :
c     ------------------------------------------------------------
      if(vol_mol_liq.gt.1.d25) then
c        diminuer la valeur de teb jusqu'a ce que 3 racines apparaissent
         teb = 2d0 * teb - tcrit
         nb_iter = nb_iter + 1
         goto 50
      endif

c     si la racine liquide a ete trouvee mais pas la racine vapeur :
c     ------------------------------------------------------------
      if(vol_mol_vap.gt.1.d25) then
c        Recherche de la valeur de teb conduisant a 3 racines par
c        dichotomie

c        Valeur de Teb initiale :
         tinf = teb
         tsup = tcrit
         teb = 0.5d0 * (tinf + tsup)
         nb_iter = nb_iter + 1

 51      continue

c        Resolution de l'equation d'etat a pres0 et teb donnees
c        a_function donne a(T) en Pa.m6/mol2
         apr = a_function(numero,teb,0)
         call resol_eos(teb,pres0,vol_mol_liq,vol_mol_vap,nb0,err0)
         if(err0.ne.0) goto 52

         if(nb_iter.eq.nb_iter_max) then
            err0 = 2
            goto 52
         endif
         nb_iter = nb_iter + 1

c        si la racine vapeur a ete trouvee mais pas la racine liquide :
         if(vol_mol_liq.gt.1.d25) then
            tsup = teb
            teb = 0.5d0 * (tinf + tsup)
            nb_iter = nb_iter + 1
            goto 51
         endif

c        si la racine liquide a ete trouvee mais pas la racine vapeur :
         if(vol_mol_vap.gt.1.d25) then
            tinf = teb
            teb = 0.5d0 * (tinf + tsup)
            nb_iter = nb_iter + 1
            goto 51
         endif
      endif

 52   continue

      end subroutine initialisation_teb
c
c---------------------------------------------------------------------------
c
      subroutine calcul_psat(numero,tk,psat,err0)

c Calcul de la pression de vapeur saturante d'un corps pur
c par une EoS cubique.
c
c INPUT : numero, entier = 1 ou 2 pour spécifier le constituant
c         tk (temp. specifiee) en K
c OUTPUT :  psat (bar), err0 (variable d'erreur, = 0 si tout va bien)
c
      implicit none
      include '..\source\comflash.cmn'
      real*8 apr,bpr,a_function
      real*8 tcrit,pcrit
      real*8 psat,psat0,tk,vol_mol_liq,vol_mol_vap
      real*8 lnphi_liq,lnphi_vap,lnphi
      real*8 res1,res2,res3,cur,vvap_old,vliq_old,eps
c     dpsat          : d Psat(1) / dT
c     dhvap          : enthalpie de vaporisation
      real*8 dhvap,dpsat,h_res
      logical der_calc

      integer nb000,nb_iter,nb_iter_max,err0,err1,numero

      common/ab_pr/apr,bpr
      common/deriv/der_calc

      common/der_psat/dpsat      
      
      eps = 1.d-9
      psat = 0.d0
      err0 = 0
      nb_iter_max = 150
      nb_iter = 0  ! compteur d'iterations

      bpr = b(numero) * 1.d-6 ! cm3/mol > m3/mol
      
      tcrit = tc(numero)        !en K
      pcrit = pc(numero)*1.0d+5 !en Pa
      if(isGPED) pcrit = pcrit*10d0
      if(dabs(tcrit - tk) / tcrit.lt.1.d-5) then
         psat = pcrit
         goto 31
      endif

      if(tk.gt.tcrit) then
         err0 = 1  ! T > Tc et Psat non calculable
         goto 31
      endif

c     procedure d'initialisation psat :
c     -------------------------------
      call initialisation_psat(numero,tk,psat,err0)
      if(err0.ne.0) goto 31

c     Methode de Newton pour calculer psat :
c     ------------------------------------
c     a_function donne a(T) en Pa.m6/mol2
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

      psat = psat0 - rgas * tk * (lnphi_vap - lnphi_liq) 
     &        / (vol_mol_vap - vol_mol_liq)

c     si la formule de Newton conduit a des pression < 0 :
      if(psat.lt.1.d-16) then
         psat = psat0 * dexp(-rgas*tk/psat0*(lnphi_vap-lnphi_liq) 
     &                       / (vol_mol_vap-vol_mol_liq))
      endif

c     remarque : a basse T et faible Ps, ces trois criteres sont
c                identiques et res1 = res2 = res3 car :
c                d Delta_vap ln phi = dPsat / Psat = d vGAZ / vGAZ
c                NB : d vLIQ / vLIQ = 0
      res1 = dabs(lnphi_vap - lnphi_liq)
      res2 = dabs(psat - psat0)/dabs(psat0)

      if(nb_iter>1) then
         res3 = dabs(vol_mol_liq - vliq_old)/dabs(vliq_old)
     &          + dabs(vol_mol_vap - vvap_old)/dabs(vvap_old)
      endif

      if(res1.lt.eps .and. res2.lt.eps .and. res3.lt.eps
     &   .and. nb_iter>1) then
         goto 31
      else
         goto 30
      endif

 31   continue
c    Code rajouté pour crit.for (modif APM 11/09/20)  
      if(der_calc .and. err0.eq.0) then
         if(vol_mol_vap.gt.1.d6 .or. vol_mol_liq.gt.1.d6) then
            err0 = 3
            goto 1477
         end if

         dhvap = h_res(numero,tk,vol_mol_vap)
     &                        - h_res(numero,tk,vol_mol_liq)
         dpsat = dhvap / (tk * (vol_mol_vap - vol_mol_liq))
      end if
 1477 continue
c     Fin  (modif APM 11/09/20)

      psat = psat * 1.d-5 ! en bar
      if(isGPED) psat = psat *1.0d-1 !bar---> MPa

      end subroutine calcul_psat
c
c ------------------------------------------------------------------
c
      subroutine initialisation_psat(numero,tk,psat,err0)

c Initialisation de la pression de vap. sat. : 
c   1) initialisation a partir d'une formule l'approximant connaissant 
c      tc, pc et omega.
c   2) Affinage de cette valeur : on calcule psat de sorte que la 
c      resolution de l'equation d'etat conduise a une racine liquide et 
c      une racine vapeur

      implicit none
      include '..\source\comflash.cmn'
      real*8 tcrit,pcrit,omega0,a_function
      real*8 apr,bpr

      real*8 psat,tk,vol_mol_liq,vol_mol_vap,pinf,psup
      integer nb0,nb_iter,nb_iter_max,err1,err0,numero

      common/ab_pr/apr,bpr
      
      nb_iter = 0
      nb_iter_max = 1000

      pcrit = pc(numero) * 1.d5 ! bar > Pa
      tcrit = tc(numero) ! K
      omega0 = fac(numero)
      
      psat = pcrit * 10.d0**(7.d0 / 3.d0 * (omega0 + 1.d0)
     &                    * (1.d0 - tcrit / tk))

c     La resolution de l'EoS a T=tk et P=psat amene-t-elle 3 racines ?
c     a_function donne a(T) en Pa.m6/mol2
      apr = a_function(numero,tk,0)

 50   continue
      call resol_eos(tk,psat,vol_mol_liq,vol_mol_vap,nb0,err0)
      if(err0.ne.0) goto 52

      if(nb_iter.eq.nb_iter_max) then
         err0 = 2
         goto 52
      endif
      nb_iter = nb_iter + 1

c     si la racine liq a ete trouvee mais pas la racine vap :
c     -----------------------------------------------------
      if(vol_mol_vap.gt.1.d25) then
c        diminuer la valeur de psat jusqu'a ce que 2 racines apparaissent
         psat = 2.d0 * psat - pcrit
         nb_iter = nb_iter + 1
         goto 50
      endif

c     si la racine vap a ete trouvee mais pas la racine liq :
c     -----------------------------------------------------
      if(vol_mol_liq.gt.1.d25) then
c        recherche de la valeur de psat conduisant a 2 racines par
c        dichotomie
         pinf = psat
         psup = pcrit
         psat = 0.5d0 * (pinf + psup)
         nb_iter = nb_iter + 1

 51      continue

c        Resolution de l'equation d'etat a psat et tk donnees
         call resol_eos(tk,psat,vol_mol_liq,vol_mol_vap,nb0,err0)
         if(err0.ne.0) goto 52

         if(nb_iter.eq.nb_iter_max) then
            err0 = 2
            goto 52
         endif
         nb_iter = nb_iter + 1

c        si la racine liquide a ete trouvee mais pas la racine vapeur :
         if(vol_mol_vap.gt.1.d25) then
            psup = psat
            psat = 0.5d0 * (pinf + psup)
            nb_iter = nb_iter + 1
            goto 51
         endif

c        si la racine vapeur a ete trouvee mais pas la racine liquide :
         if(vol_mol_liq.gt.1.d25) then
            pinf = psat
            psat = 0.5d0 * (pinf + psup)
            nb_iter = nb_iter + 1
            goto 51
         endif
      endif

 52   continue

      end subroutine initialisation_psat
c
c ------------------------------------------------------------------
c
      subroutine resol_eos(tk,pres0,vol_mol_liq,vol_mol_vap,nb0,err3)

c sous-programme resolvant a temperature tk et pression pres0 fixees, 
c l'EoS d'un corps pur. Il calcule les racines liquide et vapeur
c de l'equation si elles existent. Si elles n'existent pas, elles
c valent 1.d31

      implicit none
      include '..\source\comflash.cmn'

      real*8 apr,bpr
      real*8 pres0,tk,eta0,vol_mol_liq,vol_mol_vap,pression,pEoS
      real*8 pi,rho(3),coefb0,coefc0,coefd0,coefb,coefc,coefd
      real*8 qq,rr,delta,ss0,tt,theta,vm(3),res,cur
      real*8 eta_c
      integer nb_liq,nb_vap,nb0,err3,i
      character*3 state

      common/ab_pr/apr,bpr

c     Compacite critique :      
      eta_c = cur((1.d0-r1)*(1.d0-r2)**2)+cur((1.d0-r2)*(1.d0-r1)**2)
     &        +1.d0
      eta_c = 1.d0 / eta_c

      if(pres0.lt.0.d0) then
         write(*,*) 'p < 0 dans RESOL_EOS corps pur -- STOP'
         stop
      endif

c     ***********************************
c     RESOLUTION PAR LA METHODE DE CARDAN
c     ***********************************

      pi = 4.d0*datan(1.d0)
      nb0 = 0   
      rho(:) = 0.d0
        
c     Resolution de l'equation cubique en rho=1/v
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
         rho(2) = 2.d0*dsqrt(-qq)*dcos(theta/3.d0+2.d0*pi/3.d0)-
     &               coefb/3.d0
         rho(3) = 2.d0*dsqrt(-qq)*dcos(theta/3.d0+4.d0*pi/3.d0)-
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
            stop
         endif
      endif

c     Verification de la qualite de la resolution 
c     A-t-on : P specifiee = P(T,v) ?

      res = 0.d0
      if(vol_mol_liq.lt.1.d25) then
         res = res + dabs(pEoS(tk,vol_mol_liq) - pres0)
      endif
      if(vol_mol_vap.lt.1.d25) then
         res = res + dabs(pEoS(tk,vol_mol_vap) - pres0)
      endif

      if(res.lt.1.d-7) then ! cardan OK
         goto 80
      else ! pb cardan, on lance Newton
         continue
      endif

c     ***********************************
c     RESOLUTION PAR LA METHODE DE NEWTON
c     ***********************************

c     Volume molaire de la phase vapeur
      eta0 = 1.d-5
      state = 'vap'

 32   continue
      vol_mol_vap = bpr / eta0

      pression = pEoS(tk,vol_mol_vap)

      do while(pression.gt.pres0)
         eta0 = eta0 / 2.d0
         goto 32
      end do

c     calcul du nombre de racines vapeur (0 ou 1)
      call newton_resol_eos(tk,pres0,eta0,vol_mol_vap,state,err3)
      if(err3.ne.0) goto 80

      if(vol_mol_vap.gt.1.d25) then
         nb_vap = 0
      else
         nb_vap = 1
      endif

c     Recherche du volume molaire de la phase liquide
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

c     calcul du nombre de racines liquides (0 ou 1)
      if(vol_mol_liq.gt.1.d25) then
         nb_liq = 0
      else
         nb_liq = 1
      endif

c    calcul du nombre de racines de l'equation d'etat (1 ou 2)
      nb0 = nb_liq + nb_vap

 80   continue
      end subroutine resol_eos
c
c ------------------------------------------------------------------
c
      subroutine newton_resol_eos(tk,pres0,eta0,vol_mol,state,err4)

c sous-programme resolvant a temperature tk et pression pres0 fixees, 
c l'EoS d'un corps pur. Il calcule une seule racine de l'equation.
c La nature de cette racine depend de la valeur d'initialisation 
c de eta0 (=b/v), la densite reduite par rapport au covolume. 


c ------- declarations des variables
      implicit none
      include '..\source\comflash.cmn'
      integer compteur,nb_iter,nb_iter_max,err4
      real*8 apr,bpr
      real*8 pres0,tk,eta0,vol_mol,residu,eta_c,eta1
      real*8 pression,pression_prime,pEoS,dpEoS,cur
      character*3 state
      character*1 choix

      common/ab_pr/apr,bpr

c     Compacite critique :      
      eta_c = cur((1.d0-r1)*(1.d0-r2)**2)+cur((1.d0-r2)*(1.d0-r1)**2)
     &        +1.d0
      eta_c = 1.d0 / eta_c
      
c     initialisation du compteur de passages dans une boucle
      compteur = 0
      residu = 1.d-10

c     nombre maximal d'iterations autorisees
      nb_iter_max = 300

c     initialisation du compteur du nb d'iterations
      nb_iter = 0

c     initialisation de la variable eta1
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
c
c        CAS DE LA RACINE LIQUIDE
c        Lorsque eta1 < eta_c et que l'on cherche la racine liquide,
c           cela veut dire qu'elle n'existe sans doute pas.
c           On confirme ce resultat en ramenant eta1 a eta_c et on
c           poursuit le calcul. Si la methode de Newton calcule a nouveau
c           eta1 tel que eta1 < eta_c alors on confirme que la racine liquide
c           n'existe pas
         if(eta1.lt.eta_c) then
               vol_mol = 1.d31
               goto 35
         endif
      else
c
c        CAS DE LA RACINE VAPEUR
c        Lorsque eta1 > eta_c et que l'on cherche la racine vapeur,
c           cela veut dire qu'elle n'existe sans doute pas.
c           On confirme ce resultat en ramenant eta1 a eta_c et on
c           poursuit le calcul. Si la methode de Newton calcule a nouveau
c           eta1 tel que eta1 > eta_c alors on confirme que la racine vapeur
c           n'existe pas
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

c     (dP / deta) a T fixee
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
c
c ------------------------------------------------------------------
c
      function lnphi(tk,vol_mol)

c calcul du logarithme du coefficient de fugacite du corps pur
c a la pression pres0, la temperature tk 
c (volume molaire correspondant : vol_mol) en utilisant 
c l'equation de peng et robinson de parametres apr et bpr


c ------- declarations des variables

      implicit none
      include '..\source\comflash.cmn'
      real*8 apr,bpr
      real*8 pres0,tk,vol_mol,z0,lnphi,pEoS

      common/ab_pr/apr,bpr

      pres0 = pEoS(tk,vol_mol)
      z0 = pres0 * vol_mol / (rgas * tk)

      lnphi = z0 - 1.d0 - dlog(z0 - bpr*pres0/rgas/tk)
     &        + apr / (rgas*tk*bpr*(r1-r2))
     &          *dlog((vol_mol-bpr*r1)/(vol_mol-bpr*r2))

      end function lnphi
c
c ------------------------------------------------------------------
c
      function h_res(numero,tk,vol_mol)

c calcul de l'enthalpie residuelle du fluide (J / mol) a la temperature tk 
c (volume molaire correspondant : vol_mol) en utilisant 
c l'equation de peng et robinson de parametres apr et bpr

      implicit none
      include '..\source\comflash.cmn'
      integer numero
      real*8 apr,bpr
      real*8 tk,vol_mol,daprdt,z0,pres0,h_res,pEoS
      real*8 a_function

      common/ab_pr/apr,bpr

      pres0 = pEoS(tk,vol_mol)
      z0 = pres0 * vol_mol / (rgas * tk)
      daprdt = a_function(numero,tk,1) ! da / dT

      h_res = rgas * tk * (z0 - 1.d0) 
     &        + 1.d0/((r1-r2) * bpr) * (apr - tk*daprdt)
     &          *dlog((vol_mol-bpr*r1)/(vol_mol-bpr*r2))

      end function h_res
c
c ------------------------------------------------------------------
c
      function pEoS(tk,vol_mol)

c calcul de la pression (en Pa) du fluide a T (K) et v (m3/mol) fixes
c EoS de PR78

      implicit none
      include '..\source\comflash.cmn'
      real*8 apr,bpr
      real*8 pEoS,tk,vol_mol

      common/ab_pr/apr,bpr
     
      pEoS = rgas * tk / (vol_mol - bpr) - apr / 
     &       ((vol_mol-bpr*r1)*(vol_mol-bpr*r2))

      end function pEoS

c ------------------------------------------------------------------
c
      function dpEoS(tk,vol_mol)

c calcul de dP/d eta (en Pa) du fluide a T (K) fixe
c EoS de PR78

      implicit none
      include '..\source\comflash.cmn'
      real*8 apr,bpr
      real*8 dpEoS,tk,vol_mol,dpdvpur
      common/ab_pr/apr,bpr

      dpdvpur = -rgas*tk/(vol_mol-bpr)**2
     &         + apr*(2.d0*vol_mol - (r1+r2)*bpr)/
     &         ((vol_mol-bpr*r1)*(vol_mol-bpr*r2))**2

      dpEoS = -vol_mol**2/bpr * dpdvpur

      end function dpEoS
c
c----------------------------------------------------------

      subroutine Cpec_calc(xcomp,tk,volmol,Cpec)
c     ENTREE : tk, temperature en K
c              volmol, volume molaire en cm3
c              xcomp(1:2), fractions mol. de 1 et 2
c     SORTIE : CP,ec, capacite calo. a P cste d'ecart molaire en J/mol
c
      implicit none
      include '..\source\comflash.cmn'

      integer i,j
      real*8  xcomp(2),tk,volmol,hec,aec,sec
      real*8  amel,bmel,damel_T,damel_TT,vm
      real*8  Cvec,betaP,ktv,Cpec
      real*8  frac

c     Calcul de amix et bmix ainsi que des dérivées at et att :
      t0 = tk
      call regles_melange(0,2,xcomp)
c      call anew(2,xcomp,tk,amix,at,att,bmix)

c     amix = a/(R².T) - unités : SI
      amel = amix * rgas**2*tk  ! a(T) dans les uSI

C************calcul damel_T = da/dT *************
c     at = d[a/(R.T)] / dT
c     R.at = d(a/T) / dT = 1/T.da/dT - a/T
c     da/dT = R.T.at + a/T = R.(T.at + a/(R.T))
c     da/dT = R.(T.at + amix)
      damel_T = rgas**2*(tk*at + amix) ! da / dT dans les uSI

C************calcul damel_TT = d2a / dT2 **********
c     att = d(at) / dT
c     R.att = d(1/T.da/dT - a/T) / dT
c            = 1/T.da/dT - 2/T.da/dT + 2.a/T^3
c     da/dT= R(T.att + 2/(RT).at - 2/T.a/(R.T))
c            = R(T.att + 2/(RT).at - 2/T.amix)
      damel_TT = rgas**2*(tk*att + 2.d+0/(rgas**2*tk)*damel_T
     &           - 2.d+0/tk*amix) ! d2a / dT2 dans les uSI

c*********************************************
      bmel = bmix*rgas ! m3/mol
      vm = volmol*1.d-6 ! m3/mol

       betaP = rgas/(vm-bmel) - damel_T/(vm-bmel*r1)/(vm-bmel*r2)

       ktv = 1.d+0/(rgas*tk/(vm-bmel)**2 - amel*(2.d+0*vm-(r1+r2)*bmel)
     &        /((vm-bmel*r1)*(vm-bmel*r2))**2)

       Cvec = -tk/(bmel*(r1-r2)) * damel_TT
     &        * dlog((vm-bmel*r1)/(vm-bmel*r2))
       Cpec = Cvec - rgas + tk*ktv*betaP**2
c       write(10,*)xcomp(1:2),Cvec,ktv,betaP,Cpec
      end subroutine Cpec_calc

c----------------------------------------------------------

      subroutine hec_calc(xcomp,tk,volmol,hec)
c     ENTREE : tk, temperature en K
c              volmol, volume molaire en cm3
c              xcomp(1:2), fractions mol. de 1 et 2
c     SORTIE : hec, enthalpie d'ecart molaire en J/mol
c
      implicit none
      include '..\source\comflash.cmn'

      integer i,j
      real*8  xcomp(2),tk,volmol,hec,aec,sec
      real*8  amel,bmel,damel_T,vm
      real*8  frac

c     Calcul de amix et bmix ainsi que de la dérivée at :
      t0 = tk
      call regles_melange(0,1,xcomp)
      
c     amix = a/(R².T) - unités : SI
      amel = amix * rgas**2*tk  ! a(T) dans les uSI

c     at = d[a/(R.T)] / dT
c     R.at = d(a/T) / dT = 1/T.da/dT - a/T
c     da/dT = R.T.at + a/T = R.(T.at + a/(R.T))
c     da/dT = R.(T.at + amix)
      damel_T = rgas**2*(tk*at + amix) ! uSI

c     bmix = b/R
      bmel = bmix*rgas ! m3/mol
      vm = volmol*1.d-6 ! m3/mol

      frac = (vm-bmel*r1)/(vm-bmel*r2)
      hec  = rgas*tk*bmel/(vm-bmel)
     &       -amel*vm / (vm-bmel*r1) / (vm-bmel*r2)
     &       + 1.d0/(bmel*(r1-r2)) * (amel-tk*damel_T)*dlog(frac)
     
c      print*,'t1 = ',rgas*tk*bmel/(vm-bmel),'t2 = ',
c     &       -amel*vm / (vm-bmel*r1) / (vm-bmel*r2),'t3 = ',
c     &       + 1.d0/(bmel*(r1-r2)) * (amel-tk*damel_T)*dlog(frac),
c     &       'hec = ',hec
c      write(10,*)xcomp(1:2),hec
      end
c----------------------------------------------------------
      subroutine cp_idgaz(tk,coeffCp,cpGP)

c     ENTREE : tk, temperature en K
c              coeffCp : cinq param pour cp du gaz parfait
c                        issus de la DIPPR (en J/kmol/K)
c
c     SORTIE : cpGP, CP du gaz parfait en J/mol/K
c
      implicit none

      real*8  tk,coeffCp(5),cpGP,CsurT,EsurT

      CsurT = coeffCp(3)/tk
      EsurT = coeffCp(5)/tk

      cpGP  = coeffCp(1) + coeffCp(2)*(CsurT / DSINH(CsurT))**2
     &        + coeffCp(4)*(EsurT / DCOSH(EsurT))**2 ! en J/kmol/K
      cpGP = cpGP/1000.d0 ! en J/mol/K

      end

c----------------------------------------------------------
c
      subroutine calcul_racine_stable_pur(numero,tk,pbar00,vm,err1)
c
c INPUT : numero (du corps pur => 1 ou 2)
c         tk (temp. specifiee) en K, pbar00 en bar
c OUTPUT :  vm (m3/mol), err0 (variable d'erreur, = 0 si tout va bien)
c
      implicit none
      include '..\source\comflash.cmn'
      real*8 apr,bpr
      real*8 pres,tk,vm1,vm2,vm,cur
      real*8 lnphi1,lnphi2,lnphi,r11,r22
      real*8 pbar00,a_function
      integer err1,nb0,numero
      common/ab_pr/apr,bpr

      pres = pbar00*1.0d+5 !en Pa
      bpr = b(numero)*1.d-6 ! cm3/mol > m3/mol

c     a_function donne a(T) en Pa.m6/mol2
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
c
c======================================================================================
c
      subroutine interpolation3(xvar,yvar,xnew,ynew)
c	xvar(1:4) = [x1 ... x4]
c	yvar(1:4) = [y1 ... y4]
c     coeff = coeffs du polynome de degre 3
c
      implicit none
      integer nspec,i,j,ipiv(4)
      real*8  xvar(4),yvar(4)
      real*8  vect(4),mat(4,4)
      real*8  coeff(4),xnew,ynew


      mat = 0.d0

      do i = 1,4
         do j = 1,4
		   mat(i,j) = xvar(i)**(j-1)
         enddo
      enddo


c     Resolution de l'equation vect = coeff * mat  (inconnue : coeff)

      call lu(4,4,ipiv,mat)
      coeff(1:4) = yvar
      call back(4,4,ipiv,mat,coeff) 
      ynew = 0.d0
      do i = 1,4
         ynew = ynew + coeff(i)*xnew**(i-1)
	 end do
      end
c
c----------------------------------------------------------
c
      subroutine lu(nd,n,indx,a)
c
c     plain fortran version
c
c   lu-decomposition of matrix a by crout's method
c
c   nd:   row dimension of a in calling program
c   n:    actual size of a
c   indx: pivot vector for row interchange during factorization
c   a:    matrix to factorize; on exit: factorized matrix
c
      implicit real*8(a-h,o-z)
      dimension a(nd,*) , indx(*)
c
c     set indx-vector
c
      do 200 i = 1 , n
         indx(i) = i
         if ( i.gt.1 ) then
            do 20 k = i , n
               do 10 j = 1 , i - 1
                  a(k,i) = a(k,i) - a(j,i)*a(k,j)
 10            continue
 20         continue
         endif
c
c     pivot
c
         xmax = abs(a(i,i))
         ipiv = i
         do 50 k = i + 1 , n
            y = abs(a(k,i))
            if ( y.gt.xmax ) then
               ipiv = k
               xmax = y
            endif
 50      continue
         if ( ipiv.ne.i ) then
            indx(i) = ipiv
            do 60 k = 1 , n
               x = a(i,k)
               a(i,k) = a(ipiv,k)
               a(ipiv,k) = x
 60         continue
         endif
         fact = 1.d0/a(i,i)
         do 100 k = i + 1 , n
            a(k,i) = a(k,i)*fact
 100     continue
         do 120 k = i + 1 , n
             do 110 j = 1 , i - 1
                a(i,k) = a(i,k) - a(j,k)*a(i,j)
 110         continue
 120     continue
 200  continue
      return
      end
c----------------------------------------------------------
      subroutine back(nd,n,indx,a,v)
c
c     back-substitution for crout-based lu-routine
c
c     plain fortran version
c
c     nd:       row dimension of a
c     n:        actual size of a
c     indx:     pivot vector calculated in lu
c     a:        factorized matrix
c     v:        rhs-vector, on exit soln. to linear eqns.
c
      implicit real*8(a-h,o-z)
      dimension a(nd,*) , v(*) , indx(*)
c
c     reorder
c
      do 100 i = 1 , n - 1
         ipiv = indx(i)
         if ( ipiv.ne.i ) then
            x = v(i)
            v(i) = v(ipiv)
            v(ipiv) = x
         endif
 100  continue
c
c     l-invers
c
      do 200 i = 2 , n
         do 150 j = 1 , i - 1
            v(i) = v(i) - v(j)*a(i,j)
 150     continue
 200  continue
c
c     u-invers
c
      do 300 i = n , 1 , -1
         if ( i.ne.n ) then
             do 220 j = i + 1 , n
                v(i) = v(i) - v(j)*a(i,j)
  220        continue
         endif
         v(i) = v(i)/a(i,i)
 300  continue
      end
c=====================================================================      
      function cur(x)

c **********
c
c Definition de la fonction : x --> cur(x) = x^(1 / 3)
c Domaine de definition     : ]-infini ; + infini[
c
c **********

      implicit none
      real*8 x,cur
      if(x.le.0.d0) then
         cur = -dabs(x)**(1.d0 / 3.d0)
      else
         cur = x**(1.d0 / 3.d0)
      endif
      end  
c
c**********************************************************************
c
       subroutine longueur (nom,long)
       implicit none
c
c  sous-programme de recherche de la longueur d'un mot
c
       integer i,long
       character*80 nom
       character*3 blanc
c
       data blanc/'   '/
c
       long = 0
       do i=1,80
         if(nom(i:i+4).eq.blanc) goto 10
         long = long + 1
       end do
c
 10    return
       end      
c
c----------------------------------------------------------------------
c 
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      
      