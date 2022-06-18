c -----------------------------------------------------------------------

      subroutine plot_des(tmin,tmax,pmin,pmax,file_in,file_out,
     &                    n,code_plot,abscisse,ordonnee,titre)

c -----------------------------------------------------------------------

c input :
c    file_in  : nom du fichier texte (80 caract. max) 
c               qui contient les s‚ries de donn‚es
c               ces s‚ries sont s‚par‚es par une ligne du couple
c               (8888.0 8888.0)
c    file_out : nom du fichier texte (80 caract. max) sans extension
c               qui peut ˆtre converti en hgl
c    n        : nombre de s‚ries (maximum 50)


c output
c    un fichier outfile.des que l'on peut lire avec ecran.exe et
c    convertir en hgl

c fonctionnement
c    le programme place tous les couples (x,y) dans la mˆme s‚rie de points
c    lorsqu'il rencontre le couple (8888.0 8888.0), il passe … une autre
c    s‚rie de points.
c    le fichier form‚ contient une d‚claration des ‚chelles en bonne et due
c    forme et peut ˆtre visualis‚ … l'aide de l'executable 'ecran'

      implicit none

      integer n,i,l1,error,lun,ij
      integer nser_ps1,nser_az1,k,long
      parameter(lun=101)
      integer*1 code_plot(50,3)

      real*8 x,y,xmax,ymax,xmin,ymin,y_min,y_max,x_min,x_max,yymax
      real*8 pasy,pasx,valeurx(50),valeury(50),huit,tmin,tmax
      real*8 pmin,pmax
 
      character*80 file_out,file_in,abscisse,ordonnee,titre
      character*50 phrase

      logical pres_temp,temp_comp,default_title
      logical affich_comment

      common/graph/   pres_temp,temp_comp
      common/plot/    nser_ps1,default_title

      huit = 8888.d0
      affich_comment = .FALSE.

      call longueur(file_in,l1)

      open(100,file=file_in(1:l1),status='old') ! fichier de donn‚es


      if(pres_temp) then
         ordonnee = 'P/bar'
         pmin = 10.d0 * pmin
         pmax = 10.d0 * pmax
      end if

c ouverture du fichier de sortie :
      call longueur(file_out,l1)
      do i = 1,l1
         if(file_out(i:i).eq.'.') then
            l1 = i - 1
            goto 200
         end if
      enddo

 200  file_out = file_out(1:l1)//'.des'

      open(lun,file=file_out)

c recherche de l'‚chelle

      y_max = pmax
      ymax  = y_max
      y_min = pmin
      ymin  = y_min


      if((y_max - y_min).gt.15000.0d+00) then
         pasy = 5000.0d+00
      elseif((y_max - y_min).gt.10000.0d+00) then
         pasy = 2000.0d+00
      elseif((y_max - y_min).gt.5000.0d+00) then
         pasy = 1000.0d+00
      elseif((y_max - y_min).gt.3000.0d+00) then
         pasy = 500.0d+00
      elseif((y_max - y_min).gt.1000.0d+00) then
         pasy = 250.0d+00
      elseif((y_max - y_min).gt.700.0d+00) then
         pasy = 200.0d+00
      elseif ((y_max - y_min).gt.600.0d+00) then
         pasy = 150.0d+00
      elseif ((y_max - y_min).gt.200.0d+00) then
         pasy = 100.0d+00
      elseif ((y_max - y_min).gt.100.0d+00) then
         pasy = 20.0d+00
      else if ((y_max - y_min).gt.50.0d+00) then
         pasy = 10.00d+00
      else if ((y_max - y_min).gt.20.0d+00) then
         pasy = 5.00d+00
      else if ((y_max - y_min).gt.7.5d+00) then
         pasy = 2.50d+00
      else if ((y_max - y_min).gt.5.0d+00) then
         pasy = 2.0d+00
      else if ((y_max - y_min).gt.3.0d+00) then
         pasy = 1.00d+00
      else if ((y_max - y_min).gt.1.0d+00) then
         pasy = 0.50d+00
      else if ((y_max - y_min).gt.0.2d+00) then
         pasy = 0.100d+00
      else if ((y_max - y_min).gt.0.02d+00) then
         pasy = 0.0100d+00
      else
         pasy = 0.00100d+00
      endif
c      if(ymin.lt.0.0d+00) ymin = 0.0d+00

c      if(ymax.ge.0) then
c         ymax = pasy*dint(0.999999 + y_max/pasy)
c      else
c         ymax = pasy*(dint(0.999999 + y_max/pasy) + 1.d0)
c      end if
c      if(ymin.ge.0) then
c         ymin = pasy*dint(y_min/pasy + 0.000001)
c      else
c         ymin = pasy*(dint(y_min/pasy + 0.000001) - 1.d0)
c      end if

      x_max = tmax
      xmax  = x_max
      x_min = tmin
      xmin  = x_min
      
      if((x_max - x_min).gt.700.0d+00) then
         pasx = 200.0d+00
c      elseif ((x_max - x_min).gt.300.0d+00) then
c         pasx = 150.0d+00
      elseif ((x_max - x_min).gt.200.0d+00) then
         pasx = 100.0d+00
      elseif ((x_max - x_min).gt.100.0d+00) then
         pasx = 20.0d+00
      else if ((x_max - x_min).gt.50.0d+00) then
         pasx = 10.00d+00
      else if ((x_max - x_min).gt.20.0d+00) then
         pasx = 5.00d+00
      else if ((x_max - x_min).gt.5.0d+00) then
         pasx = 2.50d+00
      else if ((x_max - x_min).gt.1.0d+00) then
         pasx = 0.50d+00
      else if ((x_max - x_min).gt.0.2d+00) then
         pasx = 0.100d+00
      else if ((x_max - x_min).gt.0.02d+00) then
         pasx = 0.0100d+00
      else
         pasx = 0.00100d+00
      endif

c      if(xmax.ge.0.d0) then
c         xmax = pasx*dint(0.999999 + x_max/pasx)
c      else
c         xmax = pasx*(dint(0.999999 + x_max/pasx) + 1.d0)
c      endif
c      if(xmin.ge.0.d0) then
c         xmin = pasx*dint(x_min/pasx + 0.000001)
c      else
c         xmin = pasx*(dint(x_min/pasx + 0.000001) - 1.d0)
c      endif


      if(dabs(xmin-0.d0).lt.1.d-10 .and. dabs(xmax-1.d0).lt.1.d-10) then
         write(lun,61) -0.02d0,1.02d0,xmin,pasx,
     &                 ymin,ymax,ymin,pasy
      else
         write(lun,61) xmin,xmax,xmin,pasx,
     &                 ymin,ymax,ymin,pasy
      end if
 61   format('''  ''  1',/,8(f11.4,1x))

      pasx = 2.0d+00*pasx
      pasy = 2.0d+00*pasy
      do ij = 1,2*int((xmax-xmin)/pasx+1.000001)
          valeurx(ij) = xmin + pasx*(ij-1)
      end do
      do ij = 2*int((xmax-xmin)/pasx+1.000001),2,-2
          valeurx(ij) = valeurx(ij/2)
          valeurx(ij-1) = valeurx(ij)
      end do

      do ij = 1,2*int((ymax-ymin)/pasy+1.000001)
          valeury(ij) = ymin + pasy*(ij-1)
      end do
      do ij = 2*int((ymax-ymin)/pasy+1.000001),2,-2
          valeury(ij) = valeury(ij/2)
          valeury(ij-1) = valeury(ij)
      end do

      if (xmax.lt.0.01.or.pasx.lt.0.005d+0) then
         write(lun,7772) int((xmax-xmin)/pasx+1.000001),
     &   (valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
      elseif (xmax.lt.1.1d0) then
         write(lun,70) int((xmax-xmin)/pasx+1.000001),
     &   (valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
      elseif (xmax.lt.10.0) then
         write(lun,7771) int((xmax-xmin)/pasx+1.000001),
     &   (valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
      elseif (xmax.lt.100.00000) then
         if(pasx.lt.0.05d+00) then
            write(lun,7773) int((xmax-xmin)/pasx+1.000001),
     &     (valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
         else
            write(lun,71) int((xmax-xmin)/pasx+1.000001),
     &      (valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
         endif
      elseif (xmax.lt.1000.00000) then
         write(lun,771) int((xmax-xmin)/pasx+1.000001),
     &  (valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
      else
        write(lun,6661) int((xmax-xmin)/pasx+1.000001),
     &	(valeurx(ij),ij = 1,2*int((xmax-xmin)/pasx+1.000001))
      end if
	

      if(ymin.lt.0.d0) then
         if(dabs(ymax).gt.ymin) then
            yymax = dabs(ymax)
         else
            yymax = dabs(ymin)
         end if
         if (yymax.lt.0.01.or.pasy.lt.0.005d+0) then
            write(lun,7782) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         elseif (yymax.lt.1.1d0) then
            write(lun,7782) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         elseif (yymax.lt.10.0) then
            write(lun,7781) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         elseif (yymax.lt.100.00000) then
            if(pasy.lt.0.05d+00) then
               write(lun,7783) int((ymax-ymin)/pasy+1.000001),
     &         (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
            else
               write(lun,91) int((ymax-ymin)/pasy+1.000001),
     &         (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
            endif
         elseif (yymax.lt.1000.00000) then
            write(lun,781) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         else
            write(lun,6671) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         end if

         goto 666
      endif

      if (ymax.lt.0.01.or.pasy.lt.0.005d+0) then
         write(lun,7772) int((ymax-ymin)/pasy+1.000001),
     &   (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
      elseif (ymax.lt.1.1d0) then
         write(lun,7772) int((ymax-ymin)/pasy+1.000001),
     &   (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
      elseif (ymax.lt.10.0) then
         write(lun,7771) int((ymax-ymin)/pasy+1.000001),
     &   (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
      elseif (ymax.lt.100.00000) then
         if(pasy.lt.0.05d+00) then
            write(lun,7773) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         else
            write(lun,71) int((ymax-ymin)/pasy+1.000001),
     &      (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
         endif
      elseif (ymax.lt.1000.00000) then
         write(lun,771) int((ymax-ymin)/pasy+1.000001),
     &   (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
      elseif (ymax.lt.10000.00000) then
         write(lun,6661) int((ymax-ymin)/pasy+1.000001),
     &   (valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
      else
        write(lun,6662) int((ymax-ymin)/pasy+1.000001),
     &	(valeury(ij),ij = 1,2*int((ymax-ymin)/pasy+1.000001))
      end if

 70         format(i2,2x,50(f3.1,1x,'''',f3.1,'''',1x,'3',2x))
 71         format(i2,2x,50(f4.1,1x,'''',f4.1,'''',1x,'4',2x))
 771        format(i2,2x,50(f5.1,1x,'''',f5.1,'''',1x,'5',2x))
 7771       format(i2,2x,50(f4.2,1x,'''',f4.2,'''',1x,'4',2x))
 7772       format(i2,2x,50(f5.3,1x,'''',f5.3,'''',1x,'5',2x))
 7773       format(i2,2x,50(f5.2,1x,'''',f5.2,'''',1x,'5',2x))
 6661       format(i2,2x,50(f6.1,1x,'''',f6.1,'''',1x,'6',2x))
 6662       format(i2,2x,50(f7.1,1x,'''',f7.1,'''',1x,'7',2x))

 90         format(i2,2x,50(f5.1,1x,'''',f5.1,'''',1x,'5',2x))
 91         format(i2,2x,50(f6.1,1x,'''',f6.1,'''',1x,'6',2x))
 781        format(i2,2x,50(f7.1,1x,'''',f7.1,'''',1x,'7',2x))
 7781       format(i2,2x,50(f6.2,1x,'''',f6.2,'''',1x,'6',2x))
 7782       format(i2,2x,50(f7.3,1x,'''',f7.3,'''',1x,'7',2x))
 7783       format(i2,2x,50(f7.2,1x,'''',f7.2,'''',1x,'7',2x))
 6671       format(i2,2x,50(f8.1,1x,'''',f8.1,'''',1x,'8',2x))
 6672       format(i2,2x,50(f9.1,1x,'''',f9.1,'''',1x,'9',2x))


 666  continue

      call longueur(abscisse,l1)
      write(lun,'(a6,a,a,i2)') '1 0 ''',abscisse(1:l1),'''',l1 
      call longueur(ordonnee,l1)
      write(lun,'(a6,a,a,i2)') '1 0 ''',ordonnee(1:l1),'''',l1 
      write(lun,'(a)') '2 1 0 '
      if(default_title) then
         write(lun,'(a)') '6'
      else
         write(lun,'(a)') '1'
      end if
      l1 = len_trim(titre)
      write(lun,'(a1,a,a1,2x,i3)') '''',titre(1:l1),'''',l1


      if(.not.default_title) goto 999

      titre =
     & '----------- legende [c. = courbe | p. = point] -----------'
      long = len_trim(titre)
      write(lun,'(3a,i)') '''',titre(1:long),'''',long

      long = 58

      titre =
     & '   c. rouge = lpc          |   c. violette = lpc liq-liq  '
      write(lun,'(3a,i)') '''',titre(1:long),'''',long

      titre =
     & '   c. verte = az‚o         |   c. blanche  = psat         '
      write(lun,'(3a,i)') '''',titre(1:long),'''',long

      titre =
     & '   c. bleu fonc‚e, jaune, bleu clair    = llv          '
      write(lun,'(3a,i)') '''',titre(1:long),'''',long

      titre =
     & '  o = cep | + = p. crit pur | * = p. de bancroft '//
     & '| <> = p. exp.'
      write(lun,'(3a,i)') '''',titre(1:long),'''',long

 999  continue

c series de points :
      rewind(100)

      do i = 1,n
         k = 1
         write(lun,'(3(i1,1x))') code_plot(i,1),code_plot(i,2),
     &                           code_plot(i,3)

         do
 80         continue
            IF(affich_comment) THEN
               read(100,*,iostat=error) x,y,phrase
            ELSE
               read(100,*,iostat=error) x,y
            ENDIF

            if(error.ne.0) goto 81 ! end of file
            if(dabs(x-8888.d0).lt.1.d-5.and.
     &         dabs(y-8888.d0).lt.1.d-5) then
               write(lun,82) huit,huit
               goto 81
            end if

            if(pres_temp) then
               IF(affich_comment) THEN
                  write(lun,82) x,y*10.d0,phrase
               ELSE
                  write(lun,82) x,y*10.d0
               ENDIF
            else
               IF(affich_comment) THEN
                  write(lun,82) x,y,phrase
               ELSE
                  write(lun,82) x,y
               ENDIF
            end if

            k = k + 1
         end do
 81      continue
      enddo

      if(pres_temp) then
         ordonnee = 'P/bar'
         pmin = pmin / 10.d0
         pmax = pmax / 10.d0
      end if




 82   format(1x,f18.12,5x,f18.12,1x,a)

      close(101)

      close(100)

      end 

c
c -----------------------------------------------------------------------
c
       subroutine longueur (nom,long)
       implicit none
c
c  sous-programme de recherche de la longueur d'un mot
c
       integer*4 i,long
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
c
c**********************************************************************
c
       subroutine trax(chemin,lch,fic)
c
c      chemin = la directory o— se trouve le fichier … transformer au format
c               hpgl (exemple : chemin = dessin\1-2)
c      lch = la longueur de la chaine de caractere chemin
c      fic = nom du fichier … transformer (ex: t001.des)
c      lfic = longueur de la chaine de caractšre fic
c
       character*50 titre
       character*80 chemin,acol
       character*12 fic,fic1,fic2
       integer nbgradumax
       parameter (nbgradumax = 20)
       integer*4 unitres,usource
       integer*4 lch,lfic
       character*50 tx(nbgradumax),ty(nbgradumax),tcx(nbgradumax)
       character*50 tcy(nbgradumax),tc(nbgradumax)
       integer ntx(nbgradumax),nty(nbgradumax),ncx(nbgradumax)
       integer ncy(nbgradumax),icx(nbgradumax),icy(nbgradumax)
       integer nt,ic(nbgradumax),nc(nbgradumax)
       integer wb,wh,wd,wg,x(3000),y(3000)
       integer*4 big,nligne
       real*4 xt(nbgradumax),yt(nbgradumax)
       logical comment
       common /ctex/ nxt,xt,ntx,nyt,yt,nty,nxc,icx,ncx,nyc,icy,ncy
     &              ,tcy,tx,ty,tcx
       common /zf/ unitres
       common /opt/ icar,igrad,ie,ixb,ixg,iyb,iyg,ihx,ihy
       common /cgraph/ wg,wd,wb,wh,xg,xd,x0,px,yb,yh,y0,py,ax,bx,ay,by
       parameter (usource = 18)

       call system('md '//chemin(1:lch)//'\hgl')

c       lfic = 8
       lfic = len_trim(fic)
       unitres = 19
       big = 0
       icar = 70   !  diminue la taille des caractšres et des points
       igrad = 120 !  augmente la longueur des graduations et des flšches
       ie = 140    !  augmente l'espace compt‚ pour chaque caractšre
       ixb = 480   ! d‚place vers le bas les chiffres sur l'axe des x
       ixg = 20    ! d‚place vers la gauche les chiffres sur l'axe des x
       iyb = 100   ! d‚place vers le bas les chiffres sur l'axe des y
       iyg = 160   ! d‚place vers la gauche les chiffres sur l'axe des y
       ihx = 0     ! d‚place vers le haut le texte de l'axe des x
       ihy = 40    ! d‚place vers le haut le texte de l'axe des y
c
       dim = 110
       fic1 = fic(1:lfic)
       open(usource,file = chemin(1:lch)//'\'//fic1
     &             ,status = 'old')

       fic2 = fic(1:lfic-4)//'.hgl'
       open(unitres,file = chemin(1:lch)//'\hgl\'//fic2)

       if (dim.eq.0) then
          print*, ' gauche, droite, bas, haut (de 0 … 10000) ?'
          read*, wg,wd,wb,wh
c
c       pour ‚tirer le dessin vers la droite, il faut augmenter wd
c       exemple : wg = 1000 - valeur minimale pour laisser place aux l‚gendes
c                 wd = 10000
c                 wb = 1000 - valeur minimale pour les l‚gendes
c                 wh = 5000
c
c
c       pour ‚tirer le dessin vers le haut, il faut augmenter wh
c       exemple : wg = 1000 - valeur minimale pour laisser place aux l‚gendes
c                 wd = 5000
c                 wb = 1000 - valeur minimale pour les l‚gendes
c                 wh = 10000
c

       else
          wg = 1000.
          wd = icar*dim
          wb = 1000.
          wh = wd
       endif

       read(usource,*)
       i9 = 1
       read(usource,*) xg,xd,x0,px,yb,yh,y0,py
       if(i9.eq.0) goto 3016
       read(usource,*) nxt,(xt(i),tx(i),ntx(i),i=1,nxt)
       read(usource,*) nyt,(yt(i),ty(i),nty(i),i=1,nyt)
       read(usource,*) nxc,(icx(i),tcx(i),ncx(i),i=1,nxc)
       read(usource,*) nyc,(icy(i),tcy(i),ncy(i),i=1,nyc)
       iu = 0.6*dim
       read(usource,*)
       read(usource,*) nligne
       do i = 1,nligne
          read(usource,*) 
       end do

 3016  ax = (wd - wg)/(xd - xg)
       bx = wg - ax*xg
       ay = (wh - wb)/(yh - yb)
       by = wb - ay*yb
       x0 = ax*x0 + bx
       y0 = ay*y0 + by
       px = ax*px
       py = ay*py

       do 3020 i = 1,nxt
          xt(i) = ax*xt(i) + bx
 3020  continue
       do 3030 i = 1,nyt
          yt(i) = ay*yt(i) + by
 3030  continue

       if (i9.ne.0) call axe

c       call texte(0,0,1,1,1,' ')
c       call texte(0,8500,1,1,50,titre)
       
 110   read(usource,*,end = 140) ilp,icol,icar
c
c       rajout des couleurs
c
       if(icol.ne.0) then
           open(3,file='provi')
           write(3,9020) icol
           rewind(3)
           read(3,*) acol
           write(unitres,9040) acol(1:1)
           close(3)
           call system('del provi')
       endif
 9020  format('''',i1,'''')
 9040  format('sp',a1,';')
c
c       fin du rajout des couleurs
c
c
 111   continue
c
       if (ilp.eq.2) then
          read(usource,*) xx,yy,nt,(ic(i),tc(i),nc(i),i = 1,nt)
          if (icar.eq.0) then
             nx = ax*xx + bx
             ny = ay*yy + by
          else
             nx = wg + 0.8*xx*(wd - wg)
             ny = wb + 0.8*yy*(wh - wb)
          end if
          call texte(nx,ny,nt,ic,nc,tc)
          goto 110
       end if
c
       i = 1
 120   read(usource,*) xx,yy
       if(xx.eq.8888.and.yy.eq.8888) goto 130
       if(i.gt.1000) then
          backspace(usource)
          big = 1   !on coupe la s‚rie de points
          goto 130
       endif
       x(i) = int(xx*ax + bx) 
       y(i) = int(yy*ay + by)
       i = i + 1
       goto 120
 130   n = i - 1
       if (icol.eq.0) goto 110
c
       if (ilp.eq.0) then
          if(n.gt.1) call trait(n,x,y,icar)
c          call texte(0,0,1,1,1,' ')
       else if (ilp.eq.1) then
          np = 0
          do 3040 i = 1,n
             if (x(i).ge.wg-1.and.x(i).le.wd+1.and.
     &          y(i).ge.wb-1.and.y(i).le.wh+1) then
                np = np + 1
                x(np) = x(i)
                y(np) = y(i)
             endif
 3040     continue  
          call point(np,x,y,icar)
c          call texte(0,0,1,1,1,' ')
       end if
       if(big.eq.1) then
          big = 0
          goto 111
       endif
       goto 110
c
 140   close(usource)
       close(unitres)

       comment = .false.
       if (comment) then
          write(*,57) chemin(1:lch)//'\hgl\'//fic(1:lfic-4)//'.hgl'
       end if
 57    format(1x,'fichier graphique ',a,' cree')

       end
c
c----------------------------------------------------------------------
c
       subroutine axe
       integer*4 unitres
       integer nbgradumax
       parameter (nbgradumax = 20)
       common /zf/ unitres
       character*50 tx(nbgradumax),ty(nbgradumax),tcx(nbgradumax)
       character*50 tcy(nbgradumax),tc(nbgradumax)
       character a*6,aj*6,aa*2000
       integer ntx(nbgradumax),nty(nbgradumax),ncx(nbgradumax)
       integer ncy(nbgradumax),icx(nbgradumax),icy(nbgradumax)
       integer wg,wd,wb,wh,x(3000),y(3000)
       integer nx,ny,nt,ic(nbgradumax),nc(nbgradumax)
       real*4 xt(nbgradumax),yt(nbgradumax)
       common /ctex/ nxt,xt,ntx,nyt,yt,nty,nxc,icx,ncx,nyc,icy,ncy
     &              ,tcy,tx,ty,tcx
       common /cgraph/ wg,wd,wb,wh,xg,xd,x0,px,yb,yh,y0,py,ax,bx,ay,by
       common /opt/ icar,igrad,ie,ixb,ixg,iyb,iyg,ihx,ihy
       aa = 'ip0,0,10000,10000;sr1,2;'
       write(unitres,9010) aa(1:24)
 9010  format(a)
c
       i = 1
       x(i) = wg
       y(i) = wb
       if(int(x0).gt.wg+1) then
          i = i + 1
          x(i) = x0
          y(i) = wb
       endif
110    if(x(i).le.wd+1) then
          x(i+1) = x(i)
          y(i+1) = wb - igrad
          x(i+2) = x(i)
          y(i+2) = wb
          x(i+3) = x(i) + px
          y(i+3) = wb
          i = i + 3
          goto 110
       else
          i = i - 1
       endif
       i = i + 1
       x(i) = wd
       y(i) = wb
       if(int(y0).gt.wb+1) then
          i = i + 1
          x(i) = wd
          y(i) = y0
       endif
120    if(y(i).le.wh+1) then
          y(i+1) = y(i)
          x(i+1) = wd - igrad
          y(i+2) = y(i)
          x(i+2) = wd
          y(i+3) = y(i) + py
          x(i+3) = wd
          i = i + 3
          goto 120
       else
          i = i - 1
       endif
       i = i + 1
       x(i) = wd
       y(i) = wh
       call ligne(i,x,y)
c
       i = 1
       x(i) = wg
       y(i) = wb
       if(int(y0).gt.wb+1) then
          i = i + 1
          x(i) = wg
          y(i) = y0
       endif
130    if(y(i).le.wh+1) then
          y(i+1) = y(i)
          x(i+1) = wg - igrad
          y(i+2) = y(i)
          x(i+2) = wg
          y(i+3) = y(i) + py
          x(i+3) = wg
          i = i + 3
          goto 130
       else
          i = i - 1
       endif
       i = i + 1
       x(i) = wg
       y(i) = wh
       if(int(x0).gt.wg+1) then
          i = i + 1
          x(i) = x0
          y(i) = wh
       endif
  140  if(x(i).le.wd+1) then
          x(i+1) = x(i)
          y(i+1) = wh - igrad
          x(i+2) = x(i)
          y(i+2) = wh
          x(i+3) = x(i) + px
          y(i+3) = wh
          i = i + 3
         goto 140
       else
          i = i - 1
       endif
       i = i + 1
       x(i) = wd
       y(i) = wh
       call ligne(i,x,y)
c
       if (nxt.eq.0) goto 150
       j = wb - ixb
       call convia(j,aj,naj)
       do 3010 i = 1,nxt
          j = xt(i) - ixg - ntx(i)*ie/2
          call convia(j,a,na)
          aa = 'pu;pa'//a(1:na)//','//aj(1:naj)
          naa = 6 + na + naj
          aa = aa(1:naa)//';lb'//tx(i)(1:ntx(i))//char(3)
          naa = naa + 4 + ntx(i)
          write(unitres,9010) aa(1:naa)
 3010  continue
c
 150   if (nyt.eq.0) goto 160
       do 3020 i=1,nyt
          aa = 'pu;pa'
          j = wg - igrad - iyg - nty(i)*ie 
          call convia(j,a,na)
          aa = 'pu;pa'//a(1:na)//','
          naa = 6 + na
          j = yt(i) - iyb
          call convia(j,a,na)
          aa = aa(1:naa)//a(1:na)
          naa = naa + na
          aa = aa(1:naa)//';lb'//ty(i)(1:nty(i))//char(3)
          naa = naa + 4 + nty(i)
          write(unitres,9010) aa(1:naa)
 3020  continue
c
 160   nl = 0
       do 3030 i = 1,nxc
          nl = nl + ncx(i)
 3030  continue  
       x(1) = wd + 2*igrad
       y(1) = wb
       x(2) = x(1) + (nl + 2)*ie
       y(2) = wb
       x(3) = x(2) 
       y(3) = wb - igrad
       x(4) = x(3) + igrad
       y(4) = wb
       x(5) = x(3)
       y(5) = wb + igrad
       x(6) = x(2)
       y(6) = wb
       call ligne(6,x,y)
c
       x(1) = wd + 2*igrad
       y(1) = wh - 5*igrad
       x(2) = x(1)
       y(2) = wh - igrad
       x(3) = x(1) - igrad
       y(3) = y(2) 
       x(4) = x(1) 
       y(4) = wh
       x(5) = x(1) + igrad
       y(5) = y(2)
       x(6) = x(2)
       y(6) = y(2)
       call ligne(6,x,y)
c
       if(nxc.ne.0) then
          nx = wd + 2*igrad
          ny = wb + 2*igrad + ihx
          nt = nxc
          do 3040 i=1,nt
             ic(i) = icx(i)
             nc(i) = ncx(i)
             tc(i) = tcx(i)
 3040     continue
          call texte(nx,ny,nt,ic,nc,tc)
c
          nx = wd + 3*igrad
          ny = wh - 5*igrad + ihy
          nt = nyc
          do 3050 i = 1,nt
             ic(i) = icy(i)
             nc(i) = ncy(i)
             tc(i) = tcy(i)
 3050     continue
          call texte(nx,ny,nt,ic,nc,tc)
       endif
c
       return
       end
c
c--------------------------------------------------------------------
c
       subroutine texte(nx,ny,nt,ic,nc,tc)
       integer nbgradumax
       parameter (nbgradumax = 20)
       integer*4 unitres
       common /zf/ unitres
       character tc(nbgradumax)*50,a*6,aa*2000
       integer nx,ny,nt,ic(nbgradumax),nc(nbgradumax)
       iu = 60

       call convia(nx,a,na)
       aa = 'pu;pa'//a(1:na)//','
       naa = 6 + na
       call convia(ny,a,na)
       aa = aa(1:naa)//a(1:na)//';'
       naa = naa + na + 1
       do 3010 i=1,nt
          if(i.eq.1) then
             j = 3*iu*ic(i)/2
          else
             j = 3*iu*(ic(i) - ic(i-1))/2
          endif
          if(j.ge.0) then    
             call convia(j,a,na)
             aa=aa(1:naa)//'pu;pr0,'//a(1:na)//';lb'
             naa = naa + 10 + na
          else
             call convia(-j,a,na)
             aa=aa(1:naa)//'pu;pr0,-'//a(1:na)//';lb'
             naa = naa + 11 + na
          endif
          aa = aa(1:naa)//tc(i)(1:nc(i))//char(3)
          naa = naa + nc(i) + 1
 3010  continue
       write(unitres,9010) aa(1:naa)
 9010  format(a)
       return
       end
c
c--------------------------------------------------------------------
c
       subroutine ligne(n,x,y)
       integer*4 unitres
       common /zf/ unitres
       character a*6,aa*20000
       integer n,na,naa,i,x(3000),y(3000)
       aa(1:5) = 'pu;pa'
       naa = 5
       do 3010 i=1,n
          if(i.eq.2) then
            aa = aa(1:naa-1)//';pd;pa'
            naa = naa + 5
          endif
          call convia(x(i),a,na)
          aa = aa(1:naa)//a(1:na)//','
          naa = naa + na + 1
          call convia(y(i),a,na)
          aa = aa(1:naa)//a(1:na)//','
          naa = naa + na + 1
 3010  continue
       aa = aa(1:naa-1)//';'
       write(unitres,9010) aa(1:naa)
 9010  format(a)
       return
       end
c
c--------------------------------------------------------------------
c
       subroutine point(n,x,y,icar)
       integer*4 unitres
       common /zf/ unitres
       integer n,icar,x(3000),y(3000),z(3000),t(3000),zmin,tmin,tmax
       pi = 3.14159
       iu = 60
       iu2 = 3*iu/2
c
       if(icar.eq.1) then
          do 3010 i=1,n
             z(1) = x(i) + iu2
             t(1) = y(i)
             z(2) = x(i) - iu2
             t(2) = y(i)
             call ligne(2,z,t)
             z(1) = x(i)
             t(1) = y(i) + iu2
             z(2) = x(i)
             t(2) = y(i) - iu2
             call ligne(2,z,t)
 3010     continue
       else if(icar.eq.2) then
          do 3020 i=1,n
             z(1) = x(i) + iu
             t(1) = y(i) - iu
             z(2) = x(i) - iu
             t(2) = y(i) + iu
             call ligne(2,z,t)
             z(1) = x(i) - iu
             t(1) = y(i) - iu
             z(2) = x(i) + iu
             t(2) = y(i) + iu
             call ligne(2,z,t)
 3020     continue
       else if(icar.eq.3) then
          do 3030 i=1,n
             z(1) = x(i) + iu
             t(1) = y(i) + iu
             z(2) = x(i) - iu
             t(2) = t(1)
             z(3) = z(2)
             t(3) = y(i) - iu
             z(4) = z(1)
             t(4) = t(3)
             z(5) = z(1)
             t(5) = t(1)
             call ligne(5,z,t)
 3030     continue
       else if(icar.eq.4) then
          do 3040 i=1,n
             z(1) = x(i) + iu
             t(1) = y(i) 
             z(2) = x(i) 
             t(2) = y(i) + iu2
             z(3) = x(i) - iu
             t(3) = y(i) 
             z(4) = x(i) 
             t(4) = y(i) - iu2
             z(5) = z(1)
             t(5) = t(1)
             call ligne(5,z,t)
 3040     continue
       else if(icar.eq.5) then
          do 3050 i=1,n
             z(1) = x(i) + iu2
             t(1) = y(i)
             z(2) = x(i) - iu2
             t(2) = y(i)
             call ligne(2,z,t)
             z(1) = x(i)
             t(1) = y(i) + iu2
             z(2) = x(i)
             t(2) = y(i) - iu2
             call ligne(2,z,t)
             z(1) = x(i) + iu
             t(1) = y(i) - iu
             z(2) = x(i) - iu
             t(2) = y(i) + iu
             call ligne(2,z,t)
             z(1) = x(i) - iu
             t(1) = y(i) - iu
             z(2) = x(i) + iu
             t(2) = y(i) + iu
             call ligne(2,z,t)
 3050     continue
       else if(icar.eq.6) then
          do 3060 i=1,n
             do 3055 j=1,40
                alpha = 0.05*j*pi
                z(j) = x(i) + cos(alpha)*iu
                t(j) = y(i) + sin(alpha)*iu
 3055        continue
             call ligne(40,z,t)
 3060     continue
       else if(icar.eq.31) then
          do 3070 i=1,n
             zmin = x(i) - iu
             tmax = y(i) + iu
             tmin = y(i) - iu
             jj = 2
             do 3065 j=1,2*iu+1,4
                z(jj-1) = zmin + j - 1
                t(jj-1) = tmin
                z(jj) = z(jj-1)
                t(jj) = tmax
                jj = jj + 2
 3065        continue
             jj = jj - 2     
             call ligne(jj,z,t)
 3070     continue
       else if(icar.eq.41) then
          do 3080 i=1,n
             zmin = x(i) - iu
             tmin = y(i) - iu2
             jj = 2
             do 3075 j=1,iu+1,4
                z(jj-1) = zmin + j - 1
                t(jj-1) = y(i) + 1.5*(j - 1)
                z(jj) = x(i) + j - 1
                t(jj) = tmin + 1.5*(j - 1)
                jj = jj + 2
 3075        continue
             jj = jj - 2
             call ligne(jj,z,t)
 3080     continue
       else if(icar.eq.61) then
          do 3090 i=1,n
             jj = 2
             do 3085 j=1,40
                alpha = 0.05*j*pi
                z(jj-1) = x(i)
                t(jj-1) = y(i)
                z(jj) = x(i) + cos(alpha)*iu
                t(jj) = y(i) + sin(alpha)*iu
                jj = jj + 2
 3085        continue
             jj = jj - 2
             call ligne(jj,z,t)
 3090     continue
       endif                        
       return
       end
c
c--------------------------------------------------------------------
c
       subroutine trait(n,x,y,icar)
       integer*4 unitres
       common /zf/ unitres
       integer x(3000),y(3000),z(3000),t(3000)
       integer wg,wd,wb,wh
       real*4 xr(3000),yr(3000)
       common /cgraph/ wg,wd,wb,wh,xg,xd,x0,px,yb,yh,y0,py,ax,bx,ay,by
c
       do 3005 i=1,n
           xr(i) = real(x(i))
           yr(i) = real(y(i))
 3005  continue
       iu = 60
       if(n.gt.2) goto 110
c
       if(icar.eq.0) then
          call ligne(2,x,y)
          return
       endif
c     
       dx = xr(2) - xr(1)
       dy = yr(2) - yr(1)
       dist = sqrt(dx*dx + dy*dy)
       nseg = dist/((icar + 1)*iu) + 1
       dx0 = dx/(nseg)
       dy0 = dy/(nseg)
       rap = real(icar)/real(icar + 1)
       do 3010 j=1,nseg
          z(1) = xr(1) + (j - 1)*dx0
          t(1) = yr(1) + (j - 1)*dy0
          z(2) = z(1) + rap*dx0
          t(2) = t(1) + rap*dy0
          call ligne(2,z,t)
 3010  continue
       z(1) = z(2)
       t(1) = t(2)
       z(2) = x(n)
       t(2) = y(n)
       call ligne(2,z,t)
       return

 110   io = 1
 120   np = 0
       ig = 0
       do 3020 i=io,n
          if(int(xr(i)).ge.wg-1.and.int(xr(i)).le.wd+1.and.
     &       int(yr(i)).ge.wb-1.and.int(yr(i)).le.wh+1) then
             if(i.ne.io.and.ig.eq.0) then
                ig = 1
                np = np + 1
                if(int(xr(i-1)).lt.wg-1) then
                   x(np) = wg
                   cy = (yr(i) - yr(i-1))/(xr(i) - xr(i-1))
                   y(np) = yr(i-1) + cy*(wg - xr(i-1))
                else if(int(xr(i-1)).gt.wd+1) then 
                   x(np) = wd
                   cy = (yr(i) - yr(i-1))/(xr(i) - xr(i-1))
                   y(np) = yr(i-1) + cy*(wd - xr(i-1))
                else if(int(yr(i-1)).lt.wb-1) then
                   cx = (xr(i) - xr(i-1))/(yr(i) - yr(i-1))
                   x(np) = xr(i-1) + cx*(wb - yr(i-1))
                   y(np) = wb
                else if(int(yr(i-1)).gt.wh+1) then
                   cx = (xr(i) - xr(i-1))/(yr(i) - yr(i-1))
                   x(np) = xr(i-1) + cx*(wh - yr(i-1))
                   y(np) = wh
                endif
             endif
             ig = 1                 
             np = np + 1
             x(np) = xr(i)
             y(np) = yr(i)
          else if(ig.eq.1) then
             np = np + 1
             if(int(xr(i)).lt.wg-1) then
                x(np) = wg
                cy = (yr(i) - yr(i-1))/(xr(i) - xr(i-1))
                y(np) = yr(i) + cy*(wg - xr(i))
             else if(int(xr(i)).gt.wd+1) then 
                x(np) = wd
                cy = (yr(i) - yr(i-1))/(xr(i) - xr(i-1))
                y(np) = yr(i) + cy*(wd - xr(i))
             else if(int(yr(i)).lt.wb-1) then
                cx = (xr(i) - xr(i-1))/(yr(i) - yr(i-1))
                x(np) = xr(i) + cx*(wb - yr(i))
                y(np) = wb
             else if(int(yr(i)).gt.wh+1) then
                cx = (xr(i) - xr(i-1))/(yr(i) - yr(i-1))
                x(np) = xr(i) + cx*(wh - yr(i))
                y(np) = wh
             endif
             goto 130
          endif
 3020  continue
       if(ig.eq.0) return
c
 130   if(icar.eq.0) then
          call ligne(np,x,y)
       else
          ig1 = 1
          i1 = 1
          do 3030 i=1,np
             dx = real(x(i) - x(i1))
             dy = real(y(i) - y(i1))
             dist = sqrt(dx*dx + dy*dy)
             if(ig1.eq.1) then
                np1 = 1 + i - i1
                z(np1) = x(i)
                t(np1) = y(i)
                if(int(dist).ge.icar*iu) then
                   call ligne(np1,z,t)
                   ig1 = 0
                   i1 = i + 1
                endif
             else if(ig1.eq.0) then
                if(int(dist).ge.iu) then
                   ig1 = 1
                   i1 = i + 1
                endif
             endif
 3030     continue
       endif
       io = i + 1
       goto 120      
c
       end
c
c--------------------------------------------------------------------
c
       subroutine convia(i,a,na)
       integer i,i1,ip,j,k,ii
       character a(5)*1,b*5
       logical negatif
       data b /'     '/
c
       do 3010 j=1,5
          a(j) = b(j:j)
 3010  continue
c
       if(i.ge.0) then
          if(i.le.9) then
             a(1) = char(i+48)
             na = 1
             return
          endif
       endif
       negatif = .false.

       if(i.lt.0) then
          negatif = .true.
          i = -i
       endif

c
       k = i
       ip = 5
       na = 6
       do 3020 j=1,5
c
 110   i1 = k/10**ip
       if(i1.eq.0) then
          if(j.eq.1) then
             na = ip
             ip = ip - 1
             goto 110
          else
             a(j) = char(48)
             goto 120
          endif
       endif
c
       a(j) = char(i1+48)
       k = k - i1*10**ip
c
 120   ip = ip - 1
       if(ip.eq.0) then
          a(j+1) = char(k+48)
          if(negatif) then
             na = na + 1
             do ii = na,2,-1
               a(ii) = a(ii-1)
             end do
             a(1) = "-"
          endif
          return
       endif
c
 3020  continue
       end
c
c--------------------------------------------------------------------
c
       subroutine convai(a,i)
       integer i,j,k
       character a*30
c
       do 3010 k=1,29
          if (a(k+1:k+1).eq.' ') goto 110
 3010  continue
c
 110   i = 0
       do 3020 j=1,k
          i = i + (ichar(a(j:j)) - 48)*10**(k-j)
 3020  continue
       return
       end
c
c----------------------------------------------------------------------

       module inf_nan_detection

c!!     inf_nan_detection module 
c!!     copyright(c) 2003, lahey computer systems, inc.
c!!     copies of this source code, or standalone compiled files 
c!!     derived from this source may not be sold without permission
c!!     from lahey computers systems. all or part of this module may be 
c!!     freely incorporated into executable programs which are offered
c!!     for sale. otherwise, distribution of all or part of this file is
c!!     permitted, provided this copyright notice and header are included.
c
c!!     this module exposes four elemental functions:
c!!
c!!     isnan(x)    - test for a "not a number" value
c!!
c!!     isinf(x)    - test for either a positive or negative "infinite" value
c!!
c!!     isposinf(x) - test for a positive "infinite" value
c!!
c!!     isneginf(x) - test for a negative "infinite" value
c!!
c!!     each function accepts a single or double precision real argument, and
c!!     returns a true or false value to indicate the presence of the value 
c!!     being tested for. if the argument is array valued, the function returns
c!!     a conformable logical array, suitable for use with the any function, or
c!!     as a logical mask.
c!!
c!!     each function operates by transferring the bit pattern from a real 
c!!     variable to an integer container. unless testing for + or - infinity,
c!!     the sign bit is cleared to zero. the value is exclusive ored with
c!!     the value being tested for. the integer result of the ieor function is
c!!     converted to a logical result by comparing it to zero.
c!!

      implicit none

      private

      public :: isnan, isinf, isposinf, isneginf

c    ! kind numbers for single and double precision integer containers
      integer, parameter :: single = selected_int_kind(precision(1.e0))
      integer, parameter :: double = selected_int_kind(precision(1.d0))

c    ! single precision ieee values
      integer(single), parameter :: snan    = z"7fc00000"
      integer(single), parameter :: sposinf = z"7f800000"
      integer(single), parameter :: sneginf = z"ff800000"

c    ! double precision ieee values
      integer(double), parameter :: dnan    = z"7ff8000000000000"
      integer(double), parameter :: dposinf = z"7ff0000000000000"
      integer(double), parameter :: dneginf = z"fff0000000000000"

c    ! locatation of single and double precision sign bit (intel)
c    ! subtract one because bit numbering starts at zero
      integer, parameter :: spsb = bit_size(snan) - 1
      integer, parameter :: dpsb = bit_size(dnan) - 1
    
      interface isnan
         module procedure sisnan
         module procedure disnan
      end interface   

      interface isinf
         module procedure sisinf
         module procedure disinf
      end interface   
   
      interface isposinf
         module procedure sisposinf
         module procedure disposinf
      end interface   
   
      interface isneginf
         module procedure sisneginf
         module procedure disneginf
      end interface   
   
      contains    

c  ! single precision test for nan
      elemental function sisnan(x) result(res)
         real(kind(1.e0)), intent(in) :: x
         logical :: res
         res = ieor(ibclr(transfer(x,snan),spsb), snan) == 0
      end function  

c  ! double precision test for nan
      elemental function disnan(d) result(res)
         real(kind(1.d0)), intent(in) :: d
         logical :: res
         res = ieor(ibclr(transfer(d,dnan),dpsb), dnan) == 0
      end function  
  
c  ! single precision test for inf
      elemental function sisinf(x) result(res)
         real(kind(1.e0)), intent(in) :: x
         logical :: res
         res = ieor(ibclr(transfer(x,sposinf),spsb), sposinf) == 0
      end function  

c  ! double precision test for inf
      elemental function disinf(d) result(res)
         real(kind(1.d0)), intent(in) :: d
         logical :: res
         res = ieor(ibclr(transfer(d,dposinf),dpsb), dposinf) == 0
      end function  
  
c  ! single precision test for +inf
      elemental function sisposinf(x) result(res)
         real(kind(1.e0)), intent(in) :: x
         logical :: res
         res = ieor(transfer(x,sposinf), sposinf) == 0
      end function  

c  ! double precision test for +inf
      elemental function disposinf(d) result(res)
         real(kind(1.d0)), intent(in) :: d
         logical :: res
         res = ieor(transfer(d,dposinf), dposinf) == 0
      end function  
  
c  ! single precision test for -inf
      elemental function sisneginf(x) result(res)
         real(kind(1.e0)), intent(in) :: x
         logical :: res
         res = ieor(transfer(x,sneginf), sneginf) == 0
      end function  

c  ! double precision test for -inf
      elemental function disneginf(d) result(res)
         real(kind(1.d0)), intent(in) :: d
         logical :: res
         res = ieor(transfer(d,dneginf), dneginf) == 0
      end function  
  
      end module

c----------------------------------------------------------------------


