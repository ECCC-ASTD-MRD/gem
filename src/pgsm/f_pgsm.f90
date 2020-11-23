!     
!**S/P CALCUL  MOYENNE ZONALE OU MERIDIONALE D UN CHAMP
!     
      subroutine calcul(entre, sortie, ni, nj, poids, ccoupe, cigtyp)
      implicit none
!
!AUTEUR 
!    P. SARRAZIN  DORVAL QUEBEC AVRIL 85 DRPN 
!REVISION 4.0.2
!   CONVERSION DES VARIABLES HOLLERITH EN CARACTERE
!   Y. CHARTIER AOUT 90 DRPN DORVAL QUEBEC. 
!
!LANGAGE RATFOR
!
!OBJET(CALCUL)
!            CALCUL MOYENNE ZONALE OU MERIDIONALE 
!            LA MOYENNE ZONALE(EST-OUEST) CONTIENT NJ POINTS
!            LA MOYENNE MERIDIONALE(NORD-SUD) CONTIENT NI POINTS
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN   ENTRE  -CHAMP(NI,NJ) LUT PAR COUPZM ROUTINE
!  IN   SORTIE -CHAMP(NI) MERIDIONAL  CHAMP(NJ) ZONALE
!               CONTENANT MOYENNE ZONALE OU MERIDIONALE
!  IN   NI     -DIMENSION DU CHAMP ENTRE EST-OUEST
!  IN   NJ     -DIMENSION DU CHAMP ENTRE NORD-SUD 
!  IN   POIDS  -UTILISE POUR LE CALCUL DE LA MOYENNE MERIDIONALE D UN 
!               CHAMP GAUSSIEN MAX NJ
!  IN   JCOUP  -COUPE ZONALE JCOUP="ZON"
!               COUPE MERIDIONALE JCOUP="MER"
!  IN   IGTYP  -TYPE DE GRILLE SI IGTYP="B" ELIMINER DERNIERE LONG
!
!APPEL
!         -VIA ROUTINE COUPZM 
!          CALL CALCUL(ENTRE, SORTIE, NI, NJ, POIDS, JCOUP, IGTYP)
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
      integer i,ni,j,nj,nii
!
      real entre(ni,nj),sortie(1),poids(1),sum,sump
      character ccoupe*8, cigtyp*1
!
!     SI COUPE MERIDIONALE
!     
      if (ccoupe.eq.'MER ') then
!
         do i=1,ni
            sum=0.0
            sump=0.0
            do j=1,nj
               sum=sum + entre(i,j)*poids(j)
               sump = sump + poids(j)
               sortie(i)=sum/sump
            enddo
         enddo
!
!   COUPE ZONALE
!
      else
         nii=ni
!     
!     SI IGTYP='B' ON ELIMINE LA DERNIERE LONGITUDE 
!     
         if (cigtyp.eq.'B') then
            nii=ni-1
         endif
!     
         do j=1,nj
            sum=0.0
            do i=1,nii
               sum=sum + entre(i,j)
!     
               sortie(j)=sum/nii
            enddo
         enddo
      endif
      return
      end
!
!**   S/P CHAMP, IDENTIFICATION DU CHAMP APPELER ROUTINE APPROPRIEE
!
   subroutine champ(nom, ipr1, ipr2, ipr3, ipr4, ipr5, ipr6, ipr7, ipr8, ipr9, ipr10, &
      ipr11,ipr12,ipr13,ipr14,ipr15,ipr16,ipr17, ipr18, ipr19, ipr20, &
      ipr21, ipr22, ipr23, ipr24, ipr25, ipr26, ipr27, ipr28, ipr29, ipr30)
      implicit none
!
!AUTEUR
!   P. SARRAZIN JANVIER 82 DRPN DORVAL P.Q. CANADA
!REVISION 4.0.2
!         CONVERSION DES VARIABLES HOLLERITH EN CARACTERE
!         MODIFS SUR LES TESTS TOUCHANT LE TABLEAU "PAIRE"
!         Y. CHARTIER AOUT 90 DRPN DORVAL QUEBEC.
!
!LANGAGE RATFOR
!
!OBJET(CHAMP)
!         POINT D'ENTREE APPELE PAR LA DIRECTIVE
!         CHAMP(NOM,  PARM1,PARM2....PARMX)
!         EX.  CHAMP(Z,1000,850,700,500)
!         LE S/P CHAMP APPELLE LE SOUS PROGRAMME APPROPRIE
!         POUR LE TYPE DE CHAMP DEMANDE
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!
!ARGUMENTS
!  IN    NOM         NOM DU CHAMP (DIRECTIVES DE L'USAGER)
!  IN    IPR1-IPR30  DESCRIPTEURS SUPLEMENTAIRES (HEURES - NIVEAUX)
!                    DEMANDEES PAR L'USAGER
!
!IMPLICITES
!
!     MODULES
   external fstcvt, argdims, grille2,epaisur,macpcp,uvectur,scalair
!
!APPEL
!          VIA DIRECTIVE
!          CHAMP(NOM,IPR1......IPR30)
!          MAXIMUM DE 30 IPR
!
!MESSAGES
!          TYPE DE GRILLE NON DEFINI DEFAUT=GRILLE P.S.NORD(2805)
!          DIRECTIVE HEURE PAS NECESSAIRE POUR MAC(ST)
!          DIRECTIVE HEURE PAS NECESSAIRE POUR PRECIP
!
!-----------------------------------------------------------------------
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer champpr,nchamp,nchmp,npar
      common / champs/ nchamp, nchmp, npar,champpr(31)
!     
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
   integer heures,nhur,nheure
   common / heures/ nhur,nheure,heures(40)
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
      integer ip1style, dateform
      common /styles/ ip1style, dateform
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!----------------------------------------------------------------------
!
!
!
   integer ihr,ihrs,nparm
   integer ipr1(2),ipr2(2),ipr3(2),ipr4(2),ipr5(2),ipr6(2),ipr7(2)
   integer ipr8(2),ipr9(2),ipr10(2),ipr11(2),ipr12(2),ipr13(2),ipr14(2),ipr15(2)
   integer ipr16(2),ipr17(2),ipr18(2),ipr19(2),ipr20(2),ipr21(2),ipr22(2),ipr23(2)
   integer ipr24(2),ipr25(2),ipr26(2),ipr27(2),ipr28(2),ipr29(2),ipr30(2)
   integer i,np,nw,trouve
   integer fstcvt,argdims
   integer kinds(30)
   integer nom(2)
   character*8 cnom, string
   real p
!
!
   if (inputmod == SEQUENTIEL) then
	   print *,'***************************************************'
	   print *,'* ON NE PEUT UTILISER LA DIRECTIVE "CHAMP"        *'
	   print *,'* AVEC UN FICHIER D ENTREE SEQUENTIEL             *'
	   print *,'*                                                 *'
	   print *,'* UTILISEZ PLUTOT LA DIRECTIVE "CHAMP_SEQ"        *'
	   print *,'***************************************************'
      return
   endif
   nchamp = min0(31,nchamp)
   champpr(30) = ipr30(1)
   champpr(29) = ipr29(1)
   champpr(28) = ipr28(1)
   champpr(27) = ipr27(1)
   champpr(26) = ipr26(1)
   champpr(25) = ipr25(1)
   champpr(24) = ipr24(1)
   champpr(23) = ipr23(1)
   champpr(22) = ipr22(1)
   champpr(21) = ipr21(1)
   champpr(20) = ipr20(1)
   champpr(19) = ipr19(1)
   champpr(18) = ipr18(1)
   champpr(17) = ipr17(1)
   champpr(16) = ipr16(1)
   champpr(15) = ipr15(1)
   champpr(14) = ipr14(1)
   champpr(13) = ipr13(1)
   champpr(12) = ipr12(1)
   champpr(11) = ipr11(1)
   champpr(10) = ipr10(1)
   champpr(9)  = ipr9(1)
   champpr(8)  = ipr8(1)
   champpr(7)  = ipr7(1)
   champpr(6)  = ipr6(1)
   champpr(5)  = ipr5(1)
   champpr(4)  = ipr4(1)
   champpr(3)  = ipr3(1)
   champpr(2)  = ipr2(1)
   if (nchamp >= 2) then
      champpr(1)  = ipr1(1)
   else
      champpr(1)  = -1
   endif
!
   kinds(30) = ipr30(2)
   kinds(29) = ipr29(2)
   kinds(28) = ipr28(2)
   kinds(27) = ipr27(2)
   kinds(26) = ipr26(2)
   kinds(25) = ipr25(2)
   kinds(24) = ipr24(2)
   kinds(23) = ipr23(2)
   kinds(22) = ipr22(2)
   kinds(21) = ipr21(2)
   kinds(20) = ipr20(2)
   kinds(19) = ipr19(2)
   kinds(18) = ipr18(2)
   kinds(17) = ipr17(2)
   kinds(16) = ipr16(2)
   kinds(15) = ipr15(2)
   kinds(14) = ipr14(2)
   kinds(13) = ipr13(2)
   kinds(12) = ipr12(2)
   kinds(11) = ipr11(2)
   kinds(10) = ipr10(2)
   kinds(9)  = ipr9(2)
   kinds(8)  = ipr8(2)
   kinds(7)  = ipr7(2)
   kinds(6)  = ipr6(2)
   kinds(5)  = ipr5(2)
   kinds(4)  = ipr4(2)
   kinds(3)  = ipr3(2)
   kinds(2)  = ipr2(2)
   kinds(1)  = ipr1(2)
!
   nw = min(argdims(1), 2)
   if (nw == 1) then
      if (nom(1) == -1) then
         cnom = ' '
      else
         write(cnom, '(A4)') nom(1)
      endif
   else
      write(cnom,'(3A4)') (nom(i), i=1,nw)
   endif
   call low2up(cnom,cnom)
!
   nchmp = nchamp
!
!
   nparm = max0(1,nchmp - 1)
   do i=1,nparm
      if (argdims(i+1) > 1) then
         p = transfer(champpr(i), p)
         call convip_plus(champpr(i), p, -1*kinds(i)-1000, ip1style, string, &
                          .false.)
         endif
   enddo
   if (.not.associated(tmplat).and.cgrtyp  /=  '*') then
      if (message) then
         write(6,*)'GRILLE NON DEFINIE ..GRILLE DE DEFAUT P.S.(2805)'
      endif
      ngr=8
      call grille2(3,51,55,26.,28.,381000.,350.,1)
   endif
!
!   VERIFIER SI DIRECTIVE HEURE EXISTE OBLIGATOIRE AVEC CHAMP
!
   if (nhur == 0) then
      if (cnom /= 'DFPR'.or.cnom /= 'DFST') then
         print *, '  ON DOIT DEFINIR DIRECTIVE HEURE '
         return
      endif
   endif
!
!
   if (npair>40) npair=40
!
   do ihrs = 1, nhur
      ihr = heures(ihrs)
!
!   CALCUL DES VECTEURS OU DE LA VITESSE DU VENT
!
      trouve=0
      do np=1,npair
         if (cnom == paire(np)(1:8)) trouve=np
      enddo
!
!  SI ON A TROUVE ON VA A L'INTERPOLATION
!
      if (trouve  /=  0) then
         vvent  = .false.
         wdvent = .false.
         if (paire(trouve)(17:20) /= '??  ') then
            vvent= .true.
         endif
         if (paire(trouve)(21:24) == 'WD  ') then
            wdvent= .true.
            vvent = .true.
         endif
         if (cgrtyp == '*') then
            write (6,*) 'GRILLE(AUCUNE) NE FONCTIONNE QUE POUR LES VARIABLES'
            write (6,*) 'PCP, EPAIS, DFST ET NUAG'
         else
            call uvectur(paire(trouve)(9:12), paire(trouve)(13:16),               paire(trouve)(17:20),ihr,nparm,champpr)
         endif
!
!     CALCUL LA DIFFERENCE ENTRE DEUX "GZ"
!
!
      else if (cnom == 'DFGZ') then
         call epaisur(ihr, nparm, champpr)
!
!     CALCUL LA DIFFERENCE ENTRE DEUX "ST ACCUMULATEUR D'AJUSTEMENT"
!
!
      else if (cnom == 'DFST') then
         if (ihrs == 1)   then
            call macpcp('ST  ', nparm, champpr)
            if (message) then
               if (nhur>1)                  write(6,*)' HEURE PAS NECESSAIRE  (ST)'
            endif
         endif
!
!
!     CALCUL LA DIFFERENCE ENTRE DEUX "PR" PRECIPITATION
!
      else if (cnom == 'DFPR') then
         if (ihrs == 1) then
            call macpcp('PR  ', nparm, champpr)
            if (message) then
               if (nhur>1) then
                  write(6,*)                     'DIRECTIVE HEURE PAS NECESSAIRE POUR PRECIP'
               endif
            endif
         endif
!
!     INTERPOLATION DES NUAGES BAS,MOYEN,HAUT
!
      else if (cnom == 'NUAG') then
         call scalair('NB  ', ihr, 1, champpr)
         call scalair('NM  ', ihr, 1, champpr)
         call scalair('NH  ', ihr, 1, champpr)
!     AUTRE NOM  (GZ,TT,DD,WW,QQ,ES,DZ,ST,PR........)
!
      else
         if (cgrtyp == '*') then
            write (6,*) 'GRILLE(AUCUNE) NE FONCTIONNE QUE POUR LES VARIABLES'
            write (6,*) 'PCP, EPAIS, DFST ET NUAG'
         else
            call scalair(cnom, ihr, nparm, champpr)
         endif
      endif
   enddo
   return
   end
!
!**s/p champ_seq  Miroir de la directive champ pour fichiers sequentiels
!
      subroutine champ_seq (listn,listip1,waitOrGo)
      implicit none
      integer listn(*),listip1(*),waitOrGo
      external ecritur,fstrwd,pgsmlir,fstprm,symetri,fstsel,fstsui,pgsmluk
      external loupneg,loupsou,fstopc,argdims,pgsmabt,imprims,grille2
      external imprime,messags,fstcvt
      external liraxez
      integer  fstinf,pgsmlir,fstprm,fstopc,fstcvt,fstsel,fstsui,fstrwd,pgsmluk
      integer ezqkdef, ezsint, ezdefset
!
!auteur  Yves Chartier drpn Dorval Quebec Avril 1996
!revision
!
!langage fortran
!
!objet(champ_seq)
!
!arguments
!  in    listn    liste de nomvar
!  in    listip1  liste de niveau
!  in    waitOrGo commutateur d'accumulation de directives
!
!
!implicites
!
!messages
!          'mauvais appel a champdif il devrait y avoir 3 arguments'
!
!
!modules  fstinf,memoir,fstprm,pgsmlir,symetry,rgscint,ecritur,argdims
!         imprime,loupsou,pgsmabt
!
!appel     via directive champ_seq(listn, listip1, waitOrGo)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer champpr,nchamp,nchmp,npar
      common / champs/ nchamp, nchmp, npar,champpr(31)
!     
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
      integer nmaxlist1,nmaxlist2,wait,go
      parameter (nmaxlist1=16,nmaxlist2=16)
      parameter (wait=0,go=1)
      character*2 listnom(nmaxlist1,nmaxlist2)
      integer listniv(nmaxlist1,nmaxlist2)
      integer ntitems,nitems1(nmaxlist1),nitems2(nmaxlist2)
      common /seq/  ntitems,nitems1,nitems2
      common /cseq/ listnom
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
   integer heures,nhur,nheure
   common / heures/ nhur,nheure,heures(40)
!
      integer ip1style, dateform
      common /styles/ ip1style, dateform
!
      character*12 etiket
      character*4 nomvar
      character*2 typvar
      character*1 cigtyp
      integer ig1,ig2,ig3,ig4,irec,iunit
      integer num1,num2,num3,nloop,deet
      integer ip1,ip2,ip3,i,j,k, date
      integer ni,nj,nk,nbits,datyp,swa, lng, dltf, ubc, extra1, extra2, extra3
      integer argdims
      logical symetri,sym,heureok,ip3ok,processed
      character*8 string
      real p
      real fbidon
      iunit = 1
!
      if (npar.ne. 3) then
         if (message) then
            write(6,*) 'DIRECTIVE CHAMP_SEQ IL DEVRAIT Y AVOIR 3 ARGUMENTS (CHAMP_SEQ)'
         endif
         return
      endif
      if (.not.associated(tmplat)) then
         if (message) then
            write(6,*)'GRILLE NON DEFINIE ..GRILLE P.S.(2805)'
         endif
         call grille2(3,51,55,26.,28.,381000.,350.,1)
      endif
!
!   trouver nombre d'arguments dans une liste (ip1,ip2,ip3)
!
      ntitems = ntitems + 1
      if (ntitems.gt.nmaxlist1) then
         print *,'*******************************************************'
         print *,'* LA LIMITE DE 16 DIRECTIVES CHAMP_SEQ A ETE DEPASSEE *'
         print *,'*******************************************************'
         call pgsmabt
      endif
      if (argdims(1).gt.nmaxlist2) then
         print *,'*******************************************************'
         print *,'* LA LIMITE DE 16 NOMS DE VARIABLES A ETE DEPASSEE    *'
         print *,'*******************************************************'
         call pgsmabt
      endif
      if (argdims(2).gt.nmaxlist2) then
         print *,'*******************************************************'
         print *,'* LA LIMITE DE 16 NIVEAUX VERTICAUX A ETE DEPASSEE    *'
         print *,'*******************************************************'
         call pgsmabt
      endif
      nitems1(ntitems) = argdims(1)
      nitems2(ntitems) = argdims(2)
      do i=1,argdims(1)
         write(listnom(ntitems,i),'(A2)') listn(i)
      enddo
      do i=1,argdims(2)
         listniv(ntitems,i) = listip1(i)
      enddo
      if (listniv(ntitems, 1) > 1000000 .and. listniv(ntitems,2)  < 0) then
         do i=1,argdims(2),2
            p = transfer(listniv(ntitems,i), p)
            call convip_plus(listniv(ntitems, i/2+1), p, &
                             -1*listniv(ntitems, (i+1))-1000, &
                             ip1style, string, .false.)
         enddo
         nitems2(ntitems) = argdims(2)/2
      endif
      if (waitOrGo.eq.WAIT) then
         return
      endif
      ier =fstrwd(lnkdiun(1))
      irec=fstsel(1,ni,nj,nk,-1,'        ',-1,-1,-1,' ','  ')
 200  irec = fstsui(1,ni,nj,nk)
      if (irec.ge.0) then
         processed = .false.
         ier = fstprm(irec, date,deet,npas,ni, nj, nk,          nbits,datyp,         ip1,ip2,ip3,typvar,nomvar,etiket,         cigt&
     &yp,ig1,ig2,ig3,ig4,         swa, lng, dltf, ubc, extra1, extra2, extra3)
!         print *,nomvar,typvar,ip1,ip2,ip3,etiket,date
         heureok = .false.
         if (heures(1).eq.-1) then
            heureok=.true.
         else
            do k=1,nhur
               if (ip2.eq.heures(k)) then
                  heureok = .true.
               endif
            enddo
         endif
 100     if (heureok.and..not.processed) then
            do i=1,ntitems
               if (.not.processed) then
                  do j=1,nitems1(i)
                     if (listnom(i,j).eq.nomvar.or.listnom(i,j).eq.' '.and..not.processed) then
                        do k=1,nitems2(i)
                           if (listniv(i,k).eq.ip1.or.listniv(i,k).eq.-1.and..not.processed) then
                              allocate(tmpif1(ni,nj))
                              allocate(tmpif2(li,lj))
                              ier=pgsmluk(tmpif1,irec,ni,nj,nk,nomvar,cigtyp)
!
                              if (nk .gt. 1) then
                                 write(6,*)'***********************************************'
                                 write(6,*)'         PGSM N ACCEPTE PAS UN          '
                                 write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (CHMPDIF)'
                                 write(6,*)'***********************************************'
                                 call pgsmabt
                              endif
                              gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit, fbidon, fbidon)
                              ier = ezdefset(gdout, gdin)
                              ier = ezsint(tmpif2, tmpif1)
                              call ecritur(tmpif2,npack,date,deet,npas,                              li,lj,1,ip1,ip2,ip3,          &
     &                    typvar,nomvar,etiket,cgrtyp,lg1,lg2,lg3,lg4)
!
                              deallocate(tmpif2)
                              deallocate(tmpif1)
                              processed=.true.
                           endif
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
         goto 200
      endif
!**  l'interpolation est terminee - On a passé a travers le fichier
      do i=1,ntitems
         do j=1,nitems2(i)
            listnom(i,j) = '  '
            listniv(i,j) = -1
         enddo
         nitems2(i)=0
      enddo
      ntitems=0
      return
      end
  subroutine chk_extrap(gdout, gdin, li, lj, ni, nj)
      implicit none
  integer gdout, gdin, li, lj, ni, nj
  real, allocatable, dimension(:,:) ::  lat_dest, lon_dest, xdest, ydest
  integer ier,i,j
  integer gdll
  external gdll
  logical outside
  real r_ni, r_nj
  allocate(lat_dest(li,lj), lon_dest(li,lj), xdest(li,lj), ydest(li,lj))
!  print *, 'chk_extrap', li, lj, ni, nj
  ier = gdll(gdout, lat_dest, lon_dest)
!  print *, lat_dest
!  print *, lon_dest
  call gdxyfll(gdin, xdest, ydest, lat_dest, lon_dest, li*lj)
  outside = .false.
  r_ni = real(ni)
  r_nj = real(nj)
!  print *, xdest
!  print *, ydest
  do j=1,lj
     do i=1,li
	if (xdest(i,j)<0.5.or.xdest(i,j)>(r_ni+0.5)) outside = .true.
	if (ydest(i,j)<0.5.or.ydest(i,j)>(r_nj+0.5)) outside = .true.
     enddo
  enddo
  if (outside) then
     print *, ' (CHK_EXTRAP) LA GRILLE DE DESTINATION CONTIENT '
     print *, '              DES POINTS HORS DE LA GRILLE SOURCE'
     print *, ' (CHK_EXTRAP) TERMINAISON ABRUPTE...'
     stop
  endif
  deallocate(lat_dest, lon_dest, xdest, ydest)
  return
  end
      subroutine chk_hy(lu_in, lu_out)
      implicit none
      integer lu_in, lu_out
      external fstinl,fstprm, fstecr, fstluk, fstinf
      integer fstinl, fstprm, fstecr, fstluk, fstinf
      integer liste(256)
      logical rewrit
      character *12 cetiket
      character *4 cnomvar,cnomx
      character *2 ctypvar
      character *1 cigtyp
      integer dateo,datev,i
      integer deet,ig1,ig2,ig3,ig4
      integer irec,irec_out,ip1,ip2,ip3,ni,nj,nk,nrecs,ier,npas
      integer nbits,cdatyp,cswa,clng,cdltf,cubc,extra1,extra2,extra3
      real hydata(4096)
      irec=fstinl(lu_in,ni,nj,nk,-1,'            ',-1,-1,-1,'  ','HY  ', liste, nrecs, 256)
      rewrit = .true.
      do i=1,nrecs
         irec = liste(i)
         ier=fstprm(irec, dateo,deet,npas,ni, nj, nk, nbits,cdatyp,         ip1,ip2,ip3,ctypvar,cnomvar,cetiket,         cigtyp,ig1&
     &,ig2, ig3, ig4, cswa, clng, cdltf, cubc,          datev, extra2, extra3)
         irec_out=fstinf(lu_out,ni,nj,nk,dateo,cetiket,ip1,ip2,ip3,ctypvar,cnomvar)
	 if (irec_out < 0) then
	    ier=fstluk(hydata, irec, ni, nj, nk)
	    ier = fstecr(hydata,hydata,-nbits,lu_out,dateo,deet,npas,            ni,nj,nk,ip1,ip2,ip3,ctypvar,cnomvar,cetiket,     &
     &       cigtyp,ig1,ig2,ig3,ig4,cdatyp,rewrit )
	 endif
      enddo
      return
      end
      subroutine chk_toctoc(lu_in, lu_out)
      implicit none
      integer lu_in, lu_out
      external fstinl,fstprm, fstecr, fstluk, fstinf, fstopl
      integer fstinl, fstprm, fstecr, fstluk, fstinf, fstopl
      integer liste(256)
      logical rewrit
      character *12 cetiket
      character *4 cnomvar,cnomx
      character *2 ctypvar
      character *1 cigtyp
      integer dateo,datev,i
      integer deet,ig1,ig2,ig3,ig4
      integer irec,irec_out,ip1,ip2,ip3,ni,nj,nk,nrecs,ier,npas
      integer nbits,cdatyp,cswa,clng,cdltf,cubc,extra1,extra2,extra3
      real hydata(4096)
      ier = fstopl('REDUCTION32',.false.,.false.)
      irec=fstinl(lu_in,ni,nj,nk,-1,'            ',-1,-1,-1,'  ','!!  ', liste, nrecs, 256)
      rewrit = .true.
      do i=1,nrecs
         irec = liste(i)
         ier=fstprm(irec, dateo,deet,npas,ni, nj, nk, nbits,cdatyp,         ip1,ip2,ip3,ctypvar,cnomvar,cetiket,         cigtyp,ig1&
     &,ig2, ig3, ig4, cswa, clng, cdltf, cubc,          datev, extra2, extra3)
         irec_out=fstinf(lu_out,ni,nj,nk,dateo,cetiket,ip1,ip2,ip3,ctypvar,cnomvar)
	      if (irec_out < 0) then
	         ier=fstluk(hydata, irec, ni, nj, nk)
	         ier = fstecr(hydata,hydata,-nbits,lu_out,dateo,deet,npas,            ni,nj,nk,ip1,ip2,ip3,ctypvar,cnomvar,cetiket,&
     &            cigtyp,ig1,ig2,ig3,ig4,cdatyp,rewrit )
	      endif
      enddo
      return
      end
      integer function chkenrpos(luin, luout, ip1, ip2, ip3)
      implicit none
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      real, dimension (:), allocatable :: ax, ay
      integer luin, luout, ip1, ip2, ip3
      integer lip1, lip2, lip3
      integer fstinf, fstluk, fstecr, fstprm
      external fstinf, fstluk, fstecr, fstprm
      integer ni1, nj1, nk1, ni2, nj2, nk2
      integer ier,ier1, ier2, yy_key, ax_key, ay_key
      logical yinyang_grid
      character*2 typvarx, typvary, grtyp, grref
      character*12 etikx, etiky
      character*4 nomx, nomy
      character*24 chaine
      integer  cdatyp
      integer dateo, deet, npas, nbits,datyp, ig1, ig2, ig3,      ig1ref, ig2ref, ig3ref, ig4ref,       swa, lng, dltf, ubc, extra1&
     &, extra2, extra3, npak
      yinyang_grid = .false.
      chkenrpos = -1
      if (mode.eq.1) then
         ax_key = fstinf(luout, ni1, nj1, nk1, -1, '            ', ip1, ip2, ip3, '  ', '>>  ')
         ay_key = fstinf(luout, ni1, nj1, nk1, -1, '            ', ip1, ip2, ip3, '  ', '^^  ')
         if (ax_key.ge.0.and.ay_key.ge.0) then
            chkenrpos = 0
         endif
         yy_key = fstinf(luout, ni1, nj1, nk1, -1, '            ', ip1, ip2, ip3, '  ', '^>  ')
         if (yy_key.ge.0) then
            chkenrpos = 1
         endif
         if (chkenrpos == 1 .and. ax_key >= 0) then
            if (yy_key > ax_key.or.yy_key > ay_key) then
               chkenrpos = 0
            endif
         endif
         if (chkenrpos >= 0) then
            return
         endif
      endif
      if (mode.eq.5) then
         chkenrpos = 0
         return
      endif
      ax_key = fstinf(luin, ni1, nj1, nk1, -1, '            ', ip1, ip2, ip3, '  ', '>>  ')
      ay_key = fstinf(luin, ni1, nj1, nk1, -1, '            ', ip1, ip2, ip3, '  ', '^^  ')
      yy_key = fstinf(luin, ni1, nj1, nk1, -1, '            ', ip1, ip2, ip3, '  ', '^>  ')
      if (ax_key.lt.0.or.ay_key.lt.0) then
         chkenrpos = -1
      else
         chkenrpos = 0
      endif
      if (yy_key.ge.0) then
         chkenrpos = 1
         yinyang_grid = .true.
      endif
      if (chkenrpos == 1 .and. ax_key >= 0) then
         if (yy_key > ax_key.or.yy_key > ay_key) then
            chkenrpos = 0
            yinyang_grid = .false.
         endif
      endif
      if (chkenrpos == -1) then
         return
      endif
      if (yinyang_grid) then
         ier = fstprm(yy_key, dateo, deet, npas, ni1, nj1, nk1, nbits,      datyp, lip1, lip2, lip3, typvarx, nomx, etikx,      grr&
     &ef, ig1ref, ig2ref, ig3ref, ig4ref,       swa, lng, dltf, ubc, extra1, extra2, extra3)
         allocate(ax(ni1*nj1*nk1))
         ier = fstluk(ax, yy_key, ni1, nj1, nk1)
         call ecritur(ax,-nbits,dateo,deet,npas,ni1,nj1,nk1,      lip1,lip2,lip3,       typvarx,nomx,etikx,grref,ig1ref,ig2ref,ig3r&
     &ef,ig4ref)
         deallocate(ax)
      else
         ier = fstprm(ax_key, dateo, deet, npas, ni1, nj1, nk1, nbits,      datyp, lip1, lip2, lip3, typvarx, nomx, etikx,      grr&
     &ef, ig1ref, ig2ref, ig3ref, ig4ref,       swa, lng, dltf, ubc, extra1, extra2, extra3)
         ier = fstprm(ay_key, dateo, deet, npas, ni2, nj2, nk2, nbits,      datyp, lip1, lip2, lip3, typvary, nomy, etiky,      grr&
     &ef, ig1ref, ig2ref, ig3ref, ig4ref,       swa, lng, dltf, ubc, extra1, extra2, extra3)
         allocate(ax(ni1*nj1*nk1))
         allocate(ay(ni2*nj2*nk2))
         ier = fstluk(ax, ax_key, ni1, nj1, nk1)
         ier = fstluk(ay, ay_key, ni2, nj2, nk2)
         call ecritur(ax,-nbits,dateo,deet,npas,ni1,nj1,nk1,      lip1,lip2,lip3,       typvarx,nomx,etikx,grref,ig1ref,ig2ref,ig3r&
     &ef,ig4ref)
         call ecritur(ay,-nbits,dateo,deet,npas,ni2,nj2,nk2,      lip1,lip2,lip3,       typvary,nomy,etiky,grref,ig1ref,ig2ref,ig3r&
     &ef,ig4ref)
         deallocate(ax)
         deallocate(ay)
      endif
      return
      end
   subroutine chk_userdate(datev)
   implicit none
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   integer datev
   if (userdate == -1) then
      datev = -1
      return
   endif
   if (userdate .eq. 1) then
      if (date3.eq.-1) then
         datev=date2
      else
         call newdate(datev, date2, date3, 3)
      endif
   else if (userdate == 0) then
      date2 = -1
      date3 = -1
      datev = -1
   else
      print *, '(chk_userdate) suspicious date value :', userdate
   endif
   return
   end subroutine chk_userdate
!**s/p chmpdif  interpole difference entre deux champs
!
      subroutine chmpdif (noment,nomsrt,ip1tab,ip2tab,ip3tab,ip1s, ip2s, ip3s)
      implicit none
      external ecritur,fstinf,pgsmlir,memoir,fstprm,symetri
      external loupneg,loupsou,fstopc,argdims,pgsmabt,imprims,grille2
      external imprime,messags,fstcvt
      external liraxez
      integer  fstinf,pgsmlir,fstprm,fstopc,fstcvt
      integer ezqkdef, ezsint, ezdefset
      integer ip1s, ip2s, ip3s
      integer incip1
!
!auteur  p.sarrazin juillet 86  drpn dorval p.q. canada
!revision
!    4.0.2 conversion des variables hollerith en caracteres
!          y. chartier aout 90 drpn dorval quebec
!    5.2   Support des grilles sources Z
!    5.7.7 npas et deet prennent la valeur du premier champ
!
!langage ratfor
!
!objet(chmpdif)
!          extraire la difference entre deux champs par rapport
!          a ip1, ip2, ou ip3 determine par l'usager:
!          chmpdif ("gz","dz",[500,1000],6,0) liste sur ip1
!          chmpdif ("pr","pr",0,[0,12],0) liste sur ip2
!          chmpdif ("tz","zt",0,12,[1,2,3,4] liste sur ip3
!          on peut changer le nom du resultat sur le fichier
!          de sortie apres interpolation de la difference
!          avec routine fstinf on extrait le record necessaire pour
!          routine fstprm qui identifie les parametres utilises
!          par routine lire
!          on reserve la memoire pour les deux champs de travail
!          routine ecritur identifit la sorte de fichier utilise
!          pour ecrire
!
!librairies
!         -source  armnsrc,drpn
!         -objet   pgsmlib,id=armnpjs.
!
!arguments
!  in    noment  nom du champ sur fichier d'entre
!  out   nomsrt  nom du champ sur fichier de sorti defaut=-1
!                nomsrt=noment
!  in    ip1tab  peut etre une liste nombre pair ou le niveau du
!                champ d'entre
!  in    ip2tab  peut etre une liste nombre pair ou l'heure du
!                champ d'entre
!  in    ip3tab  peut etre une liste nombre pair valeur determinee
!                par l'usager sur le champ d'entre
!
!
!implicites
!
!messages
!          'mauvais appel a champdif il devrait y avoir 5 arguments'
!          'record n existe pas sur fichier d entre (chmpdif)'
!
!modules  fstinf,memoir,fstprm,pgsmlir,symetry,rgscint,ecritur,argdims
!         imprime,loupsou,pgsmabt
!
!appel     via directive chmpdif (noment,nomsrt,ip1tab,ip2tab,ip3tab)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer champpr,nchamp,nchmp,npar
      common / champs/ nchamp, nchmp, npar,champpr(31)
!     
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
      integer ip1style, dateform
      common /styles/ ip1style, dateform
!
      character*12 cetiket
      character*4 cnoment, cnomsrt
      character*2 ctypvar
      character*1 cigtyp
      integer ig1,ig2,ig3,ig4,irec1,irec2
      integer jp1,jp2,jp3,jp01,jp02,jp03,jp11,jp12,jp13,ni,nj,nk,nn
      integer lesips(3),jp(3)
      integer ip1tab(40),ip2tab(40),ip3tab(40),noment,nomsrt
      integer lcl_ip1tab(40)
      integer num1,num2,num3,nloop,dat1,dat2,deet1,deet2,npas1,npas2
      integer datsrt,deetsrt,npassrt,i,j,k,ii,jj,kk,iloop,n, datev
      integer ni1,ni2,nj1,nj2,nk1,nk2,cnbits,cdatyp,cswa, clng, cdltf, cubc, extra1, extra2, extra3
      integer argdims
      logical symetri,sym
      integer iunit, chkenrpos
      real fbidon, p
      character*8 string
      integer npts
      real, dimension(:),allocatable :: lclif1, lclif2
!
      iunit = 1
      nk = 1
      if (npar.lt. 5) then
         if (message) then
            write(6,*)            'DIRECTIVE CHMPDIF IL DEVRAIT Y AVOIR AU MOINS 5 ARGUMENTS (CHMPDIF)'
         endif
         return
      endif
      if (.not.associated(tmplat).and.cgrtyp.ne.'*') then
         if (message) then
            write(6,*)'GRILLE NON DEFINIE ..GRILLE DE DEFAUT P.S.(2805)'
         endif
         ngr=8
         call grille2(3,51,55,26.,28.,381000.,350.,1)
      endif
!
!   trouver nombre d'arguments dans une liste (ip1,ip2,ip3)
!
      num1=argdims(3)
      num2=argdims(4)
      num3=argdims(5)
!
      nloop=0
!
      if (num1.gt.1) nloop=num1
      if (num2.gt.1) nloop=num2
      if (num3.gt.1) nloop=num3
!
      if (nloop.eq.0) then
         write(6,*)         ' AUCUNE LISTE [IP1], [IP2], [IP3] DIRECTIVE CHMPDIF'
         return
      endif
!
!     Ajustement pour IP1s reels
      if (argdims(3) > 1) then
         if (ip1tab(1) > 1000000 .and. ip1tab(2) < 0) then
            do i=1,num1,2
               p = transfer(ip1tab(i), p)
               call convip_plus(lcl_ip1tab(i/2+1), p, &
                                -1*ip1tab(i+1)-1000, &
                                ip1style, string, .false.)
            enddo
            num1 = num1 / 2
            nloop = nloop / 2
         else
            do i=1,num1
               lcl_ip1tab(i) = ip1tab(i)
            enddo
         endif
      else
         do i=1,num1
            lcl_ip1tab(i) = ip1tab(i)
         enddo
      endif
!
!     verifier si il ,y a plus d'une liste
!
      if (num1.gt.1.and.num2.gt.1) then
         write(6,*)' IP1 ET IP2 CONTIENNENT UNE LISTE DE VARIABLES'
         write(6,*)' VERIFIER DIRECTIVE CHMPDIF'
         return
      endif
      if (num1.gt.1.and.num3.gt.1) then
         write(6,*)' IP1 ET IP3 CONTIENNENT UNE LISTE DE VARIABLES'
         write(6,*)' VERIFIER DIRECTIVE CHMPDIF'
         return
      endif
      if (num2.gt.1.and.num3.gt.1) then
         write(6,*)' IP2 ET IP3 CONTIENNENT UNE LISTE DE VARIABLES'
         write(6,*)' VERIFIER DIRECTIVE CHMPDIF'
         return
      endif
!     execution de chaque paire dans la liste
!
      do iloop=1,nloop,2
         if (num1.gt.1) then
            i=iloop
            ii=iloop+1
            j=1
            jj=1
            k=1
            kk=1
         endif
!
         if (num2.gt.1) then
            j=iloop
            jj=iloop+1
            i=1
            ii=1
            k=1
            kk=1
         endif
!
         if (num3.gt.1) then
            k=iloop
            kk=iloop+1
            j=1
            jj=1
            i=1
            ii=1
         endif
!
         call chk_userdate(datev)
!
!     identifier le numero de chaque record avec fstinf
!
!     modification de hollerith a caractere
!
         write(cnoment,'(A4)') noment
         write(cnomsrt,'(A4)') nomsrt
         if (etikent(1) .ne. -1) then
            write(cetiket,'(3A4)') (etikent(n), n=1,nwetike)
         else
            cetiket = '        '
         endif
         if (typeent .ne. -1) then
            write(ctypvar, '(A2)') typeent
         else
            ctypvar = '  '
         endif
         irec1=fstinf(1,ni1,nj1,nk1,datev,cetiket,lcl_ip1tab(i),ip2tab(j),ip3tab(k),ctypvar,cnoment)
         irec2=fstinf(1,ni2,nj2,nk2,datev,cetiket,lcl_ip1tab(ii),ip2tab(jj),ip3tab(kk),ctypvar,cnoment)
         if (irec2 .lt. 0 .or.irec1 .lt. 0) then
          write(6,*)'RECORD N EXISTE PAS SUR FICHIER D ENTRE (CHMPDIF)'
          write(6,*)' VERIFIER NOM,IP1,IP2,IP3 SUR DIRECTIVE CHMPDIF'
         return
         endif
!
!
!
         if (nk2 .gt. 1 .or. nk1.gt.1  ) then
            write(6,*)'***********************************************'
            write(6,*)'         PGSM N ACCEPTE PAS UN          '
            write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (CHMPDIF)'
            write(6,*)'***********************************************'
            call pgsmabt
         endif
!
!     verifier dimension des deux champs d'entre
!
         if (ni1.ne.ni2.or.nj1.ne.nj2.or.nk1.ne.nk2) then
            write(6,*)' DIMENSION DES DEUX CHAMPS DE CHMPDIF DIFFERENT'
            write(6,*)' VERIFIER FICHIER D ENTREE NI, NJ, NK'
            return
         endif
!
!
!
!     identifier parametres pour champ 1
!
         ier = fstprm( irec1, dat1,deet1,npas1,         ni, nj, nk, cnbits,cdatyp,         jp01,jp02, jp03,ctypvar,cnoment,cetiket,&
     &         cigtyp, ig1,ig2,ig3,ig4,          cswa, clng, cdltf, cubc, extra1, extra2, extra3)
         if (ier .lt. 0) then
            write(6,*)' IER = FSTPRM NEGATIF VOIR CHMPDIF'
         endif
!
!     verifier si grille gaussienne ni doit etre pair
!
         if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  then
            call messags(ni)
         endif
!
!
!     lire champ no 1
!
         allocate(lclif1(ni*nj))
         if (.not.message) ier = fstopc('TOLRNC','DEBUGS',.true.)
         call chk_userdate(datev)
!
         irec1=pgsmlir(lclif1,1,ni,nj,nk,datev,cetiket,lcl_ip1tab(i),         ip2tab(j), ip3tab(k),ctypvar,cnoment,cigtyp)
!
         if (printen)  call imprime(cnoment,lclif1,ni,nj)
!
!     identifier parametres pour champ 2
!
         ier = fstprm( irec2, dat2,deet2,npas2,ni, nj, nk,          cnbits,cdatyp,         jp11,jp12,jp13,ctypvar,cnoment,cetiket, &
     &        cigtyp,ig1,ig2,ig3,ig4,         cswa, clng, cdltf, cubc, extra1, extra2, extra3)
         if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR CHMPDIF'
!
!
!     verifier si grille gaussienne ni doit etre pair
!
         if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  call messags(ni)
!
!
!     si les deux variables identiques on la transferre
!     dans le fichier sorti
!
         datsrt=dat2
!         if (dat1.eq.dat2) datsrt=dat1
         deetsrt=deet2
!         if (deet1.eq.deet2) deetsrt=deet1
         npassrt=npas2
!         if (npas1.eq.npas2) npassrt=npas1
!
!     lire champ 2
!
         npts = li*lj
         if (npts < ni*nj) then
            npts = ni*nj
         endif
         allocate(lclif2(npts))
         if (.not.message) ier = fstopc('TOLRNC','DEBUGS',.true.)
         call chk_userdate(datev)
!
         irec2=pgsmlir(lclif2,1,ni,nj,nk,datev,cetiket,lcl_ip1tab(ii),         ip2tab(jj), ip3tab(kk),ctypvar,cnoment,cigtyp)
!
!
         if (printen)  call imprime(cnoment,lclif2,ni,nj)
!
!     difference entre les deux champs
!
         nn = ni*nj
         call loupsou(lclif1,lclif2,nn)
!
!     interpolation horizontale
!
         if (cgrtyp.eq.'*') then
            ier = chkenrpos(1,2,ig1,ig2,ig3)
         else
!     #  variable symetrique oui=.true.
            gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
            ier = ezdefset(gdout, gdin)
            ier = ezsint(lclif2, lclif1)
         endif
!
!     ecrire le ip1,ip2,ip3 correspondant aux definitions
!
         if (num1.gt.1) then
            jp(1)=lcl_ip1tab(i)
            jp(2)=lcl_ip1tab(ii)
            jp(3)=ip2tab(j)
         endif
!
         if (num2.gt.1) then
            jp(1)=lcl_ip1tab(i)
            jp(2)=ip2tab(j)
            jp(3)=ip2tab(jj)
         endif
!
         if (num3.gt.1) then
            jp(1)=lcl_ip1tab(i)
            jp(2)=ip3tab(k)
            jp(3)=ip3tab(kk)
         endif
         if (cnomsrt.eq.'    ') then
            cnomsrt=cnoment
         endif
         lesips(1) = ip1s
         lesips(2) = ip2s
         lesips(3) = ip3s
         do i=1,3
            if (lesips(i).eq.65001) then
               jp(i) = jp01
            endif
            if (lesips(i).eq.65002) then
               jp(i) = jp11
            endif
            if (lesips(i).eq.65003) then
               jp(i) = jp02
            endif
            if (lesips(i).eq.65004) then
               jp(i) = jp12
            endif
            if (lesips(i).eq.65005) then
               jp(i) = jp03
            endif
            if (lesips(i).eq.65006) then
               jp(i) = jp13
            endif
         enddo
!
!
!     ecrire sur fichier standard,ms,sequentiel
!
!
      if (cgrtyp.eq.'*') then
         call ecritur(lclif1,npack,datsrt,deetsrt,npassrt,ni,nj,nk,         jp(1),jp(2),jp(3),         ctypvar,cnomsrt,cetiket,cigt&
     &yp,ig1,ig2,ig3,ig4)
      else
         call ecritur(lclif2,npack,datsrt,deetsrt,npassrt,         li,lj,1,jp(1),jp(2),jp(3),         ctypvar,cnomsrt,cetiket,cgrty&
     &p,lg1,lg2,lg3,lg4)
      endif
!
!     remetre espace des champs de travail
!
      deallocate(lclif2)
      deallocate(lclif1)
!
      enddo
!
      return
      end
!
!**S/P   COMME   LIRE UN CHAMP DANS ACCUMULATEUR
!
      subroutine comme(iunit, nom, type, idat, niv, ihr, ip3, etiqet)
      implicit none
      external fstinf,pgsmlir,memoir,fstprm,pgsmabt,imprime
      external fstopc,messags,fstcvt
      integer fstinf,pgsmlir,fstprm,fstopc,fstcvt
!
!AUTEUR Y. CHARTIER
!
!LANGAGE FORTRAN 77
!
!OBJET(COMME)
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN   NOM     NOM DU CHAMP LCAR(GZ),"TT"......
!   IN   TYPE    TYPE DE CHAMP "P"=PREVISION  "A" ANALYSE
!   IN   NIV     NIVEAU DU CHAMP
!   IN   IHR     HEURE DU CHAMP
!   IN   IP3     LIBRE(USAGER) COMPTEUR POUR MOYENE UTILISER PAR ECRITS
!   IN   ETIQET  ETIQUETTE 10 CARACTERES
!
!IMPLICITES
!MESSAGES
!         RECORD N EXISTE PAS SUR FICHIER (FSTINF DANS COMME)
!         RECORD N EXISTE PAS (PGSMLIR DANS ROUTINE COMME)
!
!MODULES  FSTINF,PGSMABT,FSTPRM,MEMOIR,PGSMLIR
!
!APPEL     VIA DIRECTIVE
!         LIREE(NOM, TYPE, IDAT, NIV, IHR, IP3, ETIQUET)
!         LIRES(NOM, TYPE, IDAT, NIV, IHR, IP3, ETIQUET)
!
! -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
      integer blancs
      data blancs /4H    /
!
!
      character *12 cetiqet
      character *4 cnomvar
      character *2 ctypvar
      character *1 cigtyp
      character *4 cbidon
      integer etiqet(3),idat,ihr, ip3,irec1,iunit,niv,nom,num1,type
      integer cnbits,cdatyp,extra1,extra2,extra3,cubc,cdltf,clng,cswa
      integer iopc, bidon,i
      integer ezqkdef, ezgxprm,gdll, argdims, letiket(3)
      external  ezqkdef, ezgxprm, gdll, argdims
!
!
!     MODIFICATION DE HOLLERITH A CARACTERE
!
      bidon = 0
      cnomvar = '    '
      ctypvar = '  '
      cetiqet = '            '
      cigtyp  = ' '
      letiket(1) = etiqet(1)
      letiket(2) = blancs
      letiket(3) = blancs
      if (argdims(8).gt.1) then
         letiket(2) = etiqet(2)
      endif
      if (argdims(8).gt.2) then
         letiket(3) = etiqet(3)
      endif
 100  ier = fstcvt(    nom,   type, letiket,    bidon,              cnomvar,ctypvar,cetiqet,cbidon,.true.)
      if (etiqet(1) .ne. -1) then
         write(cetiqet,'(3A4)') (etiqet(i), i=1,argdims(9))
      else
         cetiqet = '            '
      endif
      if (cnomvar=='    '.and.ctypvar.eq.'  '.and.cetiqet=='            '.and. niv == -1 .and.ihr==-1.and.ip3==-1.and.idat==-1) the&
     &n
         print *, '****************************************************'
         print *, '* Parametres de selection trop vagues... Sorry...  *'
         print *, '****************************************************'
         call pgsmabt
      endif
!      print *, 'COMME', iunit
!      print *, idat, cetiqet, niv,ihr,ip3,ctypvar,cnomvar
      irec1=fstinf(iunit,nni,nnj,nnk,idat,cetiqet,niv,ihr,ip3,      ctypvar,cnomvar)
      if (irec1 .lt. 0)   then
         write(6,*)         'RECORD N EXISTE PAS (routine COMME)'
         call pgsmabt
      endif
!
      if (nnk.gt.1)   then
         write(6,*)'*************************************************'
         write(6,*)'         PGSM N ACCEPTE PAS UN          '
         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (LIREE-LIRES)'
         write(6,*)'*************************************************'
         call pgsmabt
      endif
!
!
!  #  clef pour directive pluse,moinse,ecrits....
!
!
      ier = fstprm(irec1,idatt,ideet,npas,nni,nnj,nnk, cnbits,cdatyp, &
         jpp1,jpp2,jpp3,ctypvar,cnomvar,cetiqet,cigtyp,igg1,igg2,igg3,      igg4,cswa, clng, cdltf, cubc,&
         extra1, extra2, extra3)
      if (ier .lt. 0) then
         write(6,*)'RECORD N EXISTE PAS (comme)'
         call pgsmabt
      endif
!    ALLOCATION DE LA MEMOIRE
!
      if (cigtyp.ne.'Z'.and.cigtyp.ne.'Y') then
         gdout = ezqkdef(nni,nnj, cigtyp,igg1,igg2,igg3,igg4,iunit)
         ier = ezgxprm(gdout,li,lj,cgrtyp,         lg1,lg2,lg3,lg4,cgtypxy,ig1ref,ig2ref,ig3ref,ig4ref)
         allocate(tmplon(li,lj))
         allocate(tmplat(li,lj))
      else
         print *, 'avant gritp12', igg1, igg2, igg3
	 if (iunit == 1) then
            call gritp12(7,igg1,igg2,igg3)
         else
            call gritp12(8,igg1,igg2,igg3)
         endif
      endif
!
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
      return
!
      end
!
!**S/P CONLALO   CALCUL LAT LONG DE CHAQUE PT D'UNE GRILLE TYPE "Y" OU "Z"
!
      subroutine conlalo(lat,lon,ni,nj,grtyp,grtypxy,ig1,ig2,ig3,ig4)
      implicit none
      external conlal2
      integer lat,lon,ni,nj,grtyp,grtypxy,ig1,ig2,ig3,ig4
      character*1 cgrtyp, cgtypxy
      write(cgrtyp    , '(A1)') grtyp
      write(cgtypxy, '(A1)') grtypxy
      write(6,101) cgrtyp, cgtypxy
 101  format(' CONLALO:','CGRTYP: ',a1, 'CGTYPXY: ', a1)
      call conlal2(lat,lon,ni,nj,cgrtyp,cgtypxy,ig1,ig2,ig3,ig4)
      return
      end
      subroutine conlal2(lat,lon,ni,nj,cgrtyp,cgtypxy,ig1,ig2,ig3,ig4)
!     
      implicit none
!     
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(CONLALO)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE DE TYPE "Y" OU "Z"
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
      external cigaxg,llfxy,pgsmabt,messags
!     
      integer ni,nj,ig1,ig2,ig3,ig4,i,j,hem
      real lat(ni,nj),lon(ni,nj),dlat,lat0,dlon,lon0
      real pii,pjj,d60,dgrw,buflat,buflon,dla,dlo
      real xlat1,xlon1,xlat2,xlon2
!     
      character*1 cgrtyp, cgtypxy
!     
      if (cgrtyp.eq.'Z') then
!
!   ATTENTION BOUCLE SUIVANTE NJ PERMET D'AVOIR DES NI>NJ SANS PROBLEME
!
         do j=nj,1,-1
            lat(1,j)=lat(j,1)
         enddo
!     
         do i=1,ni
            do j=1,nj
               lat(i,j)=lat(1,j)
               lon(i,j)=lon(i,1)
            enddo
         enddo
      endif
!     
      hem=1
      if (cgtypxy.eq.'S') hem=2
!     
      if (cgtypxy.eq.'N'.or.cgtypxy.eq.'S')  then
         call cigaxg(cgtypxy,pii,pjj,d60,dgrw,ig1,ig2,ig3,ig4)
!     
         do i=1,ni
            do j=1,nj
               buflat=lat(i,j) - pjj
               buflon=lon(i,j) - pii
               call llfxy(dla,dlo,buflon,buflat,d60,dgrw,hem)
               if (dlo.le.0.0) dlo=dlo + 360.0
               lat(i,j)=dla
               lon(i,j)=dlo
            enddo
         enddo
      else if (cgtypxy.eq.'L') then
         call cigaxg(cgtypxy,lat0,lon0,dlat,dlon,ig1,ig2,ig3,ig4)
!     
         do i=1,ni
            do j=1,nj
               lat(i,j) = lat(i,j)*dlat + lat0
               lon(i,j) = lon(i,j)*dlon + lon0
            enddo
         enddo
      else
         write(6,*)' TYPE DE GRILLE PAS "N","S","L" '
         write(6,*) 'DIRECTIVE GRILLE(TAPE1/TAPE2.....OUCH???'
         call pgsmabt
      endif
!     
      return
      end
      subroutine conlale(lat,lon,latg,long,ni,nj,                    cgrtyp,cgtypxy,ig1,ig2,ig3,ig4)
!     
      implicit none
!
!AUTEUR   - Y. Chartier DRPN Dorval Avril 94
!
!LANGAGE - RATFOR
!
!OBJET(CONLALO)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE DE TYPE "Y" OU "Z"
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
      external cigaxg,llfxy,pgsmabt,messags
!
      integer ni,nj,ig1,ig2,ig3,ig4,i,j,hem
      real lat(ni,nj),lon(ni,nj),latg(ni,nj),long(ni,nj),dlat,      lat0,dlon,lon0
      real pii,pjj,d60,dgrw,buflat,buflon,dla,dlo
      real xlat1,xlon1,xlat2,xlon2
!     
      character*1 cgrtyp, cgtypxy
!
      if (cgrtyp.eq.'Z') then
         do j=nj,1,-1
            lat(1,j)=lat(j,1)
         enddo
         do i=1,ni
            do j=1,nj
               latg(i,j)=lat(1,j)
               long(i,j)=lon(i,1)
            enddo
         enddo
!
         call cigaxg(cgtypxy,xlat1,xlon1,xlat2,xlon2,ig1,ig2,ig3,ig4)
         call ez_gfllfxy(lon,lat,long,latg,ni*nj,xlat1,xlon1,xlat2,xlon2)
      else
         write(6,*)' TYPE DE GRILLE PAS "E"'
         write(6,*) 'DIRECTIVE GRILLE(TAPE1/TAPE2.....OUCH???'
         call pgsmabt
      endif
!     
      return
      end
!     
!**   S/P  CONVER, PLMNMOD ECART AU CHAMP ET MULTIPLIER PAR FACTEUR
!     
!     AUTEUR P.SARRAZIN MAI 82 DRPN DORVAL QUEBEC CANADA
!     
      subroutine conver(z, ni, nj, cnom)
      implicit none
!     
!LANGAGE RATFOR
!
!OBJET(CONVER)
!          AUGMENTE UN CHAMP D UNE VALEUR UNIFORME (ECARTS) ET
!          MULTIPLIER PAR UN FACTEUR APPROPRIE ELIMINE LES VALEURS TROP PETITES 
!          OU TROP GRANDES PREDETERMINEES
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN-OUT Z   - CHAMP(NI,NJ) QUI SERA MODIFIE
!  IN     NI  - NOMBRE DE PTS DANS LA DIRECTION EST-OUET
!  IN     NJ  - NOMBRE DE PTS DANS LA DIRECTION NORD-SUD
!  IN     NOM - NOM DE LA VARIABLE DU CHAMP A MODIFIER
!
!APPEL   VIA CALL
!        CALL CONVER(CHAMP,NI,NJ,NOM)  APPELLE DANS ECRITUR 
!
!  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
!
      integer ncon,nomb
      real ecarts,facts,bass,hauts
      common/tabls/ncon,nomb,ecarts(256),facts(256),bass(256),hauts(256)
      character*4 nomss(256)
      common /ctabls/ nomss
!
!
!
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
!
      integer ni,nj,i,k,j
      real z(ni,nj)
      character*4 cnom
!     
      if ( nomb .eq. 0 ) return
      k = 0
      do i = 1,nomb
         if (cnom .eq. nomss(i)) then
            k = i
            goto 10
         endif
      enddo
 10   continue
      if ( k .eq. 0 ) return
!     
      do j=1,nj
         do i=1,ni
            z(i,j) = amax1(bass(k),amin1(hauts(k),(z(i,j) +             ecarts(k))*facts(k)))
         enddo
      enddo
!     
      return
      end
!
!**   s/p convs, batir une table avec noms,ecart,facteur,bas,haut
!
      subroutine convs(nom,  ecart,  facteur, bas, haut)
      implicit none
!
!auteur p. sarrazin mai 82 drpn dorval quebec canada
!
!     langage ratfor
!
!     objet(convs)
!     la directive convs assigne a chaque table la valeur appropriee
!     les tables augmentent a chaque appel convs
!
!arguments
!     in    nom    nom du champ que l on veut modifier
!     in    ecart  valeur plmnmodr au champ
!     in    facteur valeur utiliser pour multiplication
!     champ(i,j)=(champ(i,j) + ecart)*facteur
!     in    bas    valeur < bas auront la valeur de bas
!     in    haut   valeur du champ > haut auront la valeur de haut
!
!implicites
!     appel   via directive
!     conv(nom,ecart,facteur,bas,haut)
!
!     messages
!     plus de 40 changements d'echelle routine convs
!
!   - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer ncon,nomb
      real ecarts,facts,bass,hauts
      common/tabls/ncon,nomb,ecarts(256),facts(256),bass(256),hauts(256)
      character*4 nomss(256)
      common /ctabls/ nomss
!
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
!
!     trouver si le nom existe
!     oui- remplacer
!     non- plmnmodr
!
      external fstcvt
      integer  nom,i,ier, fstcvt
      character*4 cnom, cnoma
      real bas,ecart,facteur,haut
!
      write(cnom,'(A4)') nom
      cnoma = cnom
      i = 1
 10   if (i.le.nomb) then
         if (cnoma.ne.nomss(i)) then
            i = i + 1
            goto 10
         endif
      endif
!
!     definir nomb danger si plus grand que 40
!
      nomb = max0(nomb,i)
      if (nomb.lt.256) then
         hauts(i)=1.e+30
         bass(i)= -hauts(i)
         ecarts(i) = ecart
         facts(i) = facteur
         nomss(i) = cnoma
         if (ncon.ge.4) bass(i)=bas
         if (ncon.eq.5) hauts(i)=haut
      else
         if (message) then
            write(6,*)'PLUS DE 256 DIRECTIVES "CONVS" --- RESTE IGNORE'
         endif
      endif
      return
      end
      subroutine coord(lescoords,mode)
      implicit none
      real lescoords(*)
      integer mode
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      external argdims
      integer argdims
      integer i,localncoords
      localncoords = argdims(1)
      if (mod(localncoords,2).ne.0) then
         print *, 'Malheureux(se)! Nombre de coordonnees impaires!'
         stop
      endif
      if (mode.eq.0) then
         ncoords = 0
      endif
      do i=1,localncoords,2
         if (i.lt.nmaxcoords) then
            coordll(ncoords+i/2+1,1) = lescoords(i)
            coordll(ncoords+i/2+1,2) = lescoords(i+1)
         else
            print *, '(COORD) TROP DE POINTS! MAX=', nmaxcoords, ' !!!'
         endif
      enddo
      ncoords = ncoords+localncoords/2
      return
      end
!
!**s/p coupe, calcul coupe zonale\meridionale
!
      subroutine pgcoupe(nom,lcoupe,ipr1,ipr2,ipr3,ipr4,ipr5,ipr6,      ipr7,ipr8,ipr9,ipr10,ipr11,ipr12,ipr13,ipr14,ipr15,ipr16,  &
     &    ipr17,ipr18,ipr19,ipr20,ipr21,ipr22,ipr23,ipr24,ipr25,      ipr26,ipr27,ipr28,ipr29,ipr30)
      implicit none
      external coupzm,messags,fstcvt,pgsmabt
      integer  fstcvt
!     
!     auteur p. sarrazin avril 85 drpn dorval p.q. canada
!     
!     revision 
!     4.0.2 - conversion en caracteres de toutes les variables
!     de type hollerith
!     y. chartier- dorval quebec juillet 90 drpn.
!
!     
!     langage ratfor
!     
!     objet(coupe)
!     lire un champ sur une grille "g","l","b","c","a"  et calcul une
!     moyenne zonale est-ouest ou meridionale nord-sud pour chaque niveau
!     max(30) de l'usager 
!     
!librairies
!     -source  armnsrc,drpn
!     -objet   pgsmlib,id=armnpjs.
!     
!     arguments
!     in    nom        nom du champ requis.....z,tt,es......
!     in    lcoupe     lcoupe=lcar(zon)  coupe zonale est-ouest
!     lcoupe=lcar(mer)  coupe meridionale nord-sud
!     in    ipr1-ipr30 niveau de l'usager optionel
!     
!     implicites
!     
!     
!modules
!     coupzm
!     
!     appel
!     via directive
!     moyent(nom,lcoup,ipr..........)
!     moysrt(nom,lcoup,ipr..........)
!     maximum de 30 ipr
!     
!     messages 
!     pas assez d'arguments directive moyent/moysrt
!
!     
!     implicites
!     
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!     
!     
   integer heures,nhur,nheure
   common / heures/ nhur,nheure,heures(40)
!
!     
!     
      integer nivospr,nmoy,nmo
      common / nivos/ nmoy, nmo, nivospr(31)
!
!
!
!     
!     
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
!     
!     
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
!     
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
!
!
      integer champpr,nchamp,nchmp,npar
      common / champs/ nchamp, nchmp, npar,champpr(31)
!     
!
!
!
!
      character cnomvar*4,ctypvar*2,cigtyp*1,cetiket*12
      common /cfldinf/ cnomvar, ctypvar, cetiket, cigtyp
!
!
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
!
      integer nom,lcoupe,ipr1,ipr2,ipr3,ipr4,ipr5,ipr6,ipr7,ipr8
      integer ipr9,ipr10,ipr11,ipr12,ipr13,ipr14,ipr15,ipr16,ipr17
      integer ipr18,ipr19,ipr20,ipr21,ipr22,ipr23,ipr24,ipr25
      integer ipr26,ipr27,ipr28,ipr29,ipr30
      integer iunit,nparm,i
!     
      character*8 cjcoup
!     
!     initialiser nivospr
!     
      do i=1,30
         nivospr(i)=-1
      enddo
!     
!     
      iunit=1
 1000 nmoy = min0(31,nmoy)
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,        21,22,23,24,25,26,27,28,29,30,31) nmoy
 31   nivospr(30) = ipr30
 30   nivospr(29) = ipr29
 29   nivospr(28) = ipr28
 28   nivospr(27) = ipr27
 27   nivospr(26) = ipr26
 26   nivospr(25) = ipr25
 25   nivospr(24) = ipr24
 24   nivospr(23) = ipr23
 23   nivospr(22) = ipr22
 22   nivospr(21) = ipr21
 21   nivospr(20) = ipr20
 20   nivospr(19) = ipr19
 19   nivospr(18) = ipr18
 18   nivospr(17) = ipr17
 17   nivospr(16) = ipr16
 16   nivospr(15) = ipr15
 15   nivospr(14) = ipr14
 14   nivospr(13) = ipr13
 13   nivospr(12) = ipr12
 12   nivospr(11) = ipr11
 11   nivospr(10) = ipr10
 10   nivospr(9) = ipr9
 9    nivospr(8) = ipr8
 8    nivospr(7) = ipr7
 7    nivospr(6) = ipr6
 6    nivospr(5) = ipr5
 5    nivospr(4) = ipr4
 4    nivospr(3) = ipr3
 3    nivospr(2) = ipr2
 2    nivospr(1) = ipr1
!     
!   1 nomvar = nom
 1    continue
!
!     sauve la valeur de nmoy a cause de readlx
      nmo = nmoy
!     
      nparm = max0(1,nmo-2)
!     
!     
      write(cnomvar,'(a2)') nom
!     jcoup=lcoupe
      if (lcoupe.eq.1) then
         cjcoup = 'ZON'
      endif
      if (lcoupe.eq.2) then
         cjcoup = 'MER'
      endif
      call coupzm(iunit,cnomvar,cjcoup)
      return
!     
!     directive moysrt lire sur fichier de sorti
!     
      entry moysrt(nom,lcoupe,ipr1,ipr2,ipr3,ipr4,ipr5,ipr6,ipr7,ipr8,      ipr9,ipr10,ipr11,ipr12,ipr13,ipr14,ipr15,ipr16,      ip&
     &r17,ipr18,ipr19,ipr20,ipr21,ipr22,ipr23,ipr24,ipr25,      ipr26,ipr27,ipr28,ipr29,ipr30)
      iunit=2
      go to 1000
!     
      end
!
!**   s/p coupzm  coupe zonale ou meridionale d un champ
!
      subroutine coupzm(iunit, cnom, cjcoup)
      implicit none
      external calcul,ecritur,gauss,fstinf,pgsmlir,memoir,fstprm,fstcvt,          pgsmabt,imprime,loupmir,louptra,loupin1,fstopc,me&
     &ssags
      integer  fstinf,pgsmlir,fstprm,fstopc,fstcvt
!
!auteur  p. sarrazin  dorval quebec avril 85 drpn
!
!revision
!     4.0.2 - conversion en caracteres de toutes les variables
!             de type hollerith
!             y. chartier- dorval quebec juillet 90 drpn.
!
!langage ratfor
!
!objet(coupzm)
!            faire une coupe zonale ou meridionale sur un champ dont
!            le type de grille est "g"-"a"-"l"-"b"-"c"
!            calcul pour une coupe meridionale sur un champ gaussien
!            "g" est fait a partir de poids calcule par gaussg
!            la moyenne zonale(est-ouest) contient nj points
!            la moyenne meridionale(nord-sud) contient ni points
!
!librairies
!         - source pgsm/un=armnsrc
!         -  objet pgsmlib,id=armnpjs sur xmp
!
!arguments
!  in   iunit  numero du fichier a lire
!  in   cnom    nom du champ 2 caracteres gz.tt.dd.......
!  in   cjcoup  coupe zonale='zon' meridionale='mer'
!
!appel
!         -via routine coupe
!         call coupzm(iunit, cnom, cjcoup)
!
!messages
!          record n'existe pas verifier directive moyent/moysrt
!          iunit=  niveaux=   heure=   nom=
!          mauvais type de grille directive moyent/moysrt
!          grtyp=
!          doit-etre "g"-"a"-"l"-"c"-"b"
!
!modules pgsmabt,rfl,fstinf,fstprm,pgsmlir,ecritur
!
!-----------------------------------------------------------------
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
!
!
      integer dat,deet
!
!
!
!
   integer heures,nhur,nheure
   common / heures/ nhur,nheure,heures(40)
!
!
!
      integer nivospr,nmoy,nmo
      common / nivos/ nmoy, nmo, nivospr(31)
!
!
!
!
!
      character cnomvar*4,ctypvar*2,cigtyp*1,cetiket*12
      common /cfldinf/ cnomvar, ctypvar, cetiket, cigtyp
!
!
!
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      character cnom*3, cjcoup*8
      real coa(500),w(500),sia(500),rad(500),poids(500),sinm2(500),      sinm1(500),sin2(500),champ(1000)
      integer i, iunit, datev
      integer ihr,iheur,iprs,npres,irec,ni,nj,nk
      integer jp1,jp2,jp3,ig1,ig2,ig3,ig4
      integer num,ilath,j,cnbits,cdatyp,iopc,      cswa, clng, cdltf, cubc, extra1, extra2, extra3
      character*8 cdummy
      integer      dummy
      integer un
      un = 1
      cnomvar = cnom
!
!
!     loop des heures
!
      do ihr=1,nhur
         iheur=heures(ihr)
      enddo
!
      npres=max0(1,nmo-2)
      do iprs = 1,npres
!
!     identifier numero du record
!
         call chk_userdate(datev)
!
!     modification de hollerith a caractere
!
         if (etikent(1) .ne. -1) then
            write(cetiket,'(3A4)') (etikent(i), i=1,nwetike)
         else
            cetiket = '            '
         endif
         if (typeent .ne. -1) then
            write(ctypvar, '(A2)') typeent
         else
            ctypvar = '  '
         endif
         cigtyp = ' '
         irec = fstinf(iunit,ni,nj,nk,datev,cetiket,nivospr(iprs),          iheur,ip3ent,ctypvar,cnomvar)
!
         if (irec .lt. 0) then
            write(6,*)            'RECORD N EXISTE PAS VERIFIER DIRECTIVE MOYENT/MOYSRT'
            write(6,*)            ' IUNIT=',iunit,' NIVEAU=',nivospr(iprs),            ' HEURE=',iheur
            write(6,600) cnom
            return
         endif
         if (nk.gt.1) then
            write(6,*)'**********************************************'
            write(6,*)'         PGSM N ACCEPTE PAS UN          '
            write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (MOYENT-MOYSRT)'
            write(6,*)'***********************************************'
            call pgsmabt
         endif
!
!
!     identifier parametres si type g-a-b-l-c
!
!
         ier = fstprm( irec, dat,deet,npas,ni, nj, nk, cnbits,cdatyp,         jp1,jp2, jp3,ctypvar,cnomvar,cetiket,cigtyp,         &
     & ig1,ig2,ig3,ig4,          cswa, clng, cdltf, cubc, extra1, extra2, extra3)
         if (ier .lt. 0) then
            write(6,*)' IER = FSTPRM NEGATIF VOIR CHMPDIF'
         endif
!
!     verifier si grille gaussienne ni doit etre pair
!
         if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  then
            call messags(ni)
         endif
!
!
!     nombre de longitude max=1000
!
         if (ni.gt.1000) then
            write(6,*)'  TROP DE LONGITUDES CHANGER DIMENSION DANS'
            write(6,*)' COUPZM CHAMP'
            call pgsmabt
         endif
!
!     verifier si dimension nj .gt. 500
!
         if (ig1.eq.0.and.nj.gt.500) then
            write(6,*)' CHAMP D ENTRE GLOBALE POUR MOYENT-MOYSRT'
            write(6,*)' DIMENSION NJ .GT.500 STOP'
            write(6,*)' MODIFIER ROUTINE COUPZM DANS PGSM'
            call pgsmabt
         endif
!
!     verifier type de grille
!
         if (cigtyp.ne.'G'.and.cigtyp.ne.'A'.and.cigtyp.ne.'L'         .and.cigtyp.ne.'B'.and.cigtyp.ne.'C') then
            write(6,*)' MAUVAIS TYPE DE GRILLE DIRECTIVE MOYENT/MOYSRT'
            write(6,601) cigtyp
            write(6,*)' DOIT-ETRE G - L - B - A - C  (MOYENT/MOYSRT) '
            return
         endif
!
!     lire champ
!
         allocate(tmpif1(ni,nj))
!
         if (.not.message) then
            iopc= fstopc('TOLRNC','DEBUGS',.true.)
         endif
         call chk_userdate(datev)
!
         num=pgsmlir(tmpif1,iunit,ni,nj,nk,datev,cetiket,jp1,jp2,         jp3,ctypvar,cnomvar,cigtyp)
        if (printen)  then
           call imprime(cnom,tmpif1,ni,nj)
        endif
        if (num .lt. 0)  then
           write(6,*)' CHAMP N EXISTE PAS LIRE DANS MOYENT/MOYSRT'
           call pgsmabt
        endif
!
!     initialiser les poids pour grille gaussienne meridionale
!
        if (cjcoup.eq.'MER'.and.cigtyp.eq.'G') then
           ilath=nj
           if (ig1.eq.0) then
              ilath=nj/2
           endif
!
           call gauss(ilath,coa,poids,sia,rad,w,sinm1,sinm2,sin2)
!
!     sauve les poids dans coa pour renverser le champ si necessaire
!
           do j=1,ilath
              coa(j)=poids(j)
           enddo
!
!     si gaussienne globale  transfer hemis nord dans hemis sud (miroir)
!
           if (ig1.eq.0)  then
              call loupmir(poids(ilath),poids(ilath),ilath)
           endif
!
!  hemisphere sud  orientation nord-sud renverse
!
           if ((ig1.eq.2.and.ig2.eq.1) .or.(ig1.eq.1.and.ig2.eq.0))then
              call louptra(poids(ilath),coa,ilath)
           else
              call loupin1(poids(1),poids(1),nj)
           endif
        endif
!     calcul moyenne zonale ou meridional
!
!
        call calcul(tmpif1,champ,ni,nj,poids,cjcoup,cigtyp)
!
!     si type de grille 'b' on reduit ni=ni-1
!
        if (cigtyp.eq.'B') then
           ni=ni-1
        endif
!
!     ecrit champ zonal
!
        if (cjcoup.eq.'ZON') then
           call ecritur(champ,npack,dat,deet,npas,un,nj,nk,jp1,jp2,           jp3,ctypvar,cnomvar,cetiket,cigtyp,ig1,ig2,ig3,ig4)
!
!     ecrit champ  meridional
!
        else
           call ecritur(champ,npack,dat,deet,npas,ni,un,un,jp1,jp2,           jp3,ctypvar,cnomvar,cetiket,cigtyp,ig1,ig2,ig3,ig4)
        endif
        deallocate(tmpif1)
      enddo
!
 600  format(' NOM= ',a2)
 601  format(' IGTYP= ',a1)
!
      return
      end
!
!**S/P ECRITS   ECRIRE SUR FICHIER STANDARD, MS, SEQUENTIEL
!
      subroutine ecrits(nom,npac,idat,ip1,ip2,ip3,type,      etiqet,igtyp,imprim,ig1srt,ig2srt,ig3srt,ig4srt)
      implicit none
      external conver,fstecr,fclos,memoir,pgsmabt,imprims,fstopc,messags,fstcvt,putfld
      integer fstecr,fstopc,fstcvt,iopc
!
!AUTEUR  P.SARRAZIN  AOUT 82  DRPN DORVAL P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(ECRITS)
!          ECRIRE SUR FICHIER STANDARD AVEC ROUTINE ECRIRE
!          ECRIRE SUR FICHIER MS AVEC PUTFLD
!          ECRIRE SUR FICHIER SEQUENTIEL AVEC PUTFLD
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN    NOM     NOM DU CHAMP 2 CARACTERES
!   IN    NPAC    COMPACTION DU DATA DANS CHAMP
!   IN    IDAT    DATE DU CHAMP (CMC STAMP)
!   IN    IP1     NORMALEMENT NIVEAU DU CHAMP
!   IN    IP2     HEURE DU CHAMP
!   IN    IP3     LIBRE
!   IN    TYPV    TYPE DU CHAMP 1 CARACTERE
!   IN    ETIK    ETIQUETTE 1 MOT CDC (USAGER)
!   IN    IGTYP   TYPE DE GRILLE 1 CARACTERE
!   IN    IMPRIM   IMPRIMME LES ELEMENTS DES FICHIERS
!                  D'ENTRE OU DE SORTI
!
!IMPLICITES
!MESSAGES
!          MAUVAISE DIRECTIVE NUMERO=0 FICHIER MS
!          FICHIER INCONNU ROUTINE ECRITUR
!          FIN DU FICHIER CA DEBORDE
!
!MODULES  FSTECR,PUTFLD
!
!APPEL   VIA DIRECTIVE
!        ECRITS(NOM,NPACK,IDAT,IP1,IP2,IP3,TYPE ETIQET,IGTYP,IMPRIM)
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      integer ncon,nomb
      real ecarts,facts,bass,hauts
      common/tabls/ncon,nomb,ecarts(256),facts(256),bass(256),hauts(256)
      character*4 nomss(256)
      common /ctabls/ nomss
!
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer noenrg,numero,nbrow,numdel,istart
   common/enrege/ noenrg,numero,nbrow,numdel,istart
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer blancs
      data blancs /4H    /
      integer ip1style, dateform
      common /styles/ ip1style, dateform
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
      character *12 cetiqet
      character *4 cnomvar
      character *2 ctypvar
      character *1 cigtyp
      integer i
      integer nom,npac,idat,ip1(2),ip2,ip3,igtyp,imprim,npkc
      integer iun,istamp,etiqet(*),type,cdatyp
      integer ig1srt,ig2srt,ig3srt,ig4srt,ig1s,ig2s,ig3s,ig4s
      logical rewrit
      integer letiket(3)
      integer argdims
      external argdims
      integer lip1
      real p
      character*8 string
!-----------------------------------------------------------------
!
      cnomvar = '    '
      ctypvar = '  '
      cetiqet = '            '
      cigtyp  = ' '
      letiket(1) = etiqet(1)
      letiket(2) = blancs
      letiket(3) = blancs
      if (argdims(8).gt.1) then
         letiket(2) = etiqet(2)
      endif
      if (argdims(8).gt.2) then
         letiket(3) = etiqet(3)
      endif
      lip1 = ip1(1)
      if (argdims(4) > 1) then
         p = transfer(ip1(1), p)
         call convip_plus(lip1, p, -1*ip1(2)-1000, 2, string, .false.)
      endif
      ier = fstcvt(      nom,    type,  letiket,  igtyp,       cnomvar, ctypvar, cetiqet, cigtyp, .true.)
      print *, cnomvar, '--', ctypvar, '--' , cetiqet, '--', cigtyp
      if (nom.eq.-1)        cnomvar=cnumv
      if (type.eq.-1)       ctypvar=ctypv
      if (etiqet(1).eq.-1)  cetiqet=cetik
      if (igtyp.eq.-1)      cigtyp=cigty
      npkc=npac
      if (npac.eq.-1)  npkc=-16
      if (idat.eq.-1)  idat=idatt
      if (lip1.eq.-1)  lip1=jpp1
      if (ip2.eq.-1)  ip2=jpp2
      if (ip3.eq.-1)  ip3=jpp3
      if (ip3.eq.   4095)  ip3=icnt
!
      if (necrt.lt.11) then
         ig4s=igg4
         ig3s=igg3
         ig2s=igg2
         ig1s=igg1
      endif
      if (necrt.eq.11) then
         ig4s=igg4
         ig3s=igg3
         ig2s=igg2
         ig1s=ig1srt
      endif
      if (necrt.eq.12)  then
         ig4s=igg4
         ig3s=igg3
         ig2s=ig2srt
         ig1s=ig1srt
      endif
      if (necrt.eq.13) then
         ig4s=igg4
         ig3s=ig3srt
         ig2s=ig2srt
         ig1s=ig1srt
      endif
      if (necrt.eq.14) then
         ig4s=ig4srt
         ig3s=ig3srt
         ig2s=ig2srt
         ig1s=ig1srt
      endif
!
      if (necrt.gt.9) then
         write(6,660) cnumv,idatt,jpp1,jpp2,jpp3,ctypv,cetik,cigty
 660     format('*   ENTRE     ',a2,2x,i10,3x,i5,3x,i2,         3x,i3,4x,a1,4x,a10,3x,a1)
         write(6,670) cnomvar,idat,lip1,ip2,ip3,ctypvar,cetiqet,cigtyp
 670     format('*   SORTIE    ',a2,2x,i10,3x,i5,3x,i2,         3x,i3,4x,a1,4x,a10,3x,a1)
!
      endif
!
      if (ichck.eq.0)  then
         write(6,*)         'DIRECTIVE LIREN OU LIRSR DOIT-ETRE APPELE AVANT ECRITS'
         call pgsmabt
      endif
!
!     SI LE NOM EXISTE DANS LA TABLE BATIE PAR L USAGER ALORS
!     LE CHAMP EST MODIFIE
!     EX: ACUMULA(NNI,NNJ)= AMAX1(BAS,AMIN1(HAUT,(ACUMULA(NNI,NNJ) +
!     ECART)*FACT))
!
      call conver(tmpif0, nni, nnj, cnomvar)
!
      if (printsr) then
         call imprims(cnomvar,tmpif0,nni,nnj)
      endif
!
      iun=lnkdiun(idx_ozsrt)
      if (mode.eq.1) then
         if (compression_level.eq.0) then
            cdatyp = 1
         else
            if (npac <= -16) then
            cdatyp = 134
            else if (npac == -32) then
              cdatyp = 133
            else
              cdatyp = 1
            endif
!          call armn_compress_setlevel(compression_level)
        endif
         if (iwrit.eq.+1) then
            if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
            rewrit = .true.
!
            ier = fstecr(tmpif0,tmpif0,npkc,iun,idat,ideet,npas,            nni,nnj,nnk,lip1,ip2,ip3,ctypvar,cnomvar,cetiqet,cigtyp&
     &,            ig1s,ig2s,ig3s,ig4s,cdatyp,rewrit )
         else
            if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
            rewrit = .false.
!
            ier = fstecr(tmpif0,tmpif0,npkc,iun,idat,ideet,npas,            nni,nnj,nnk,lip1,ip2,ip3,ctypvar,cnomvar,cetiqet,cigtyp&
     &,            ig1s,ig2s,ig3s,ig4s,cdatyp,rewrit)
         endif
!
      else if (mode.eq.3) then
         if (valid) then
            istamp=idat
         else
            istamp=0
         endif
         call putfld(tmpif0,tmpif0,iun,0,iwrit,nni,nnj,         nbrow,npkc,istamp)
         if (message) then
            write(6,610)ctypvar,cnomvar,lip1,ip2,ip3,nni,nnj,iun
         endif
      else
         if (message) then
            write(6,*)'FICHIER INCONNU ROUTINE ECRITS'
         endif
      endif
!
!
      deallocate(tmpif0)
!
 600  format(2x,' ENREG.ECRIT ',2(a2,'- '),3(i5,'- '),      'TAILLE ',2(i5,'- '),      'FICHIER MS ',i4,'   REC=',i4)
 610  format(2x,' ENREG.ECRIT ',2(a2,'- '),3(i5,'- '),'TAILLE ',      2(i5,'- '), 'FICHIER SEQUENTIEL',i4)
!
      return
      end
!**S/P ECRITUR   ECRIRE SUR FICHIER STANDARD, MS, SEQUENTIEL
!
      subroutine ecritur(fld,npac,idat,deet,npas,ni,nj,nk,ip1,ip2,ip3,ctypvar,cnomvar,cetiket,cgtyp,llg1,llg2,llg3,llg4)
      implicit none
      external conver,fstecr,fclos,memoir,pgsmabt,      imprims,fstopc,messags,fstcvt,putfld
      integer fstopc,fstcvt,fstecr
!
!AUTEUR  P.SARRAZIN  JANVIER 82  DRPN DORVAL P.Q. CANADA
!REVISION 4.0.2
!   MODIF. VARIABLES TYPE HOLLERITH EN CARACTERE
!   Y. CHARTIER -AOUT 90- DRPN DORVAL QUEBEC.
!Revision 5.0.0
!   Elimination appels a fstcvt pour etiksrt
!   Y. Chartier -mai  91- drpn Dorval Quebec
!
!LANGAGE RATFOR
!
!OBJET(ECRITUR)
!          ECRIRE SUR FICHIER STANDARD AVEC ROUTINE FSTECR
!          ECRIRE SUR FICHIER MS AVEC PUTFLD
!          ECRIRE SUR FICHIER SEQUENTIEL AVEC PUTFLD
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN    fld   fld(NI,NJ,NK) A ECRIRE
!   IN    NPAC    COMPACTION DU DATA DANS fld
!   IN    IDAT    DATE DU fld (CMC STAMP)
!   IN    DEET    PAS DE TEMPS ENTIER SECONDES
!   IN    NPAS    NUMERO DU PAS DE TEMPS
!   IN    NI      1 ER DIMENSION DU fld
!   IN    NJ      2 IEM DIMENSION DU fld
!   IN    NK      3 IEME DIMENSION DU fld
!   IN    IP1     NORMALEMENT NIVEAU DU fld
!   IN    IP2     HEURE DU fld
!   IN    IP3ENT     LIBRE
!   IN    TYPEENT   TYPE DU fld 1 CARACTERE
!   IN    NOM     NOM DU fld 2 CARACTERES
!   IN    ETIKE   ETIQUETTE 1 MOT CDC (USAGER)
!   IN    GRTYPE  TYPE DE GRILLE 1 CARACTERE
!   IN    LLG1    DESCRIPTEUR DE GRILLE
!   IN    LLG2    PS - LLG1 POSITION J DU POLE
!   IN    LLG3         LLG2 POSITION I DU POLE
!   IN    LLG4         LLG3 DGRW*100
!                      LLG4 D60 HETOMETRE(0-36000)
!                 LAT-LON  LLG1- DLAT*100
!                          LLG2- DLON*100
!                          LLG3- (90-LAT)*100 0=POLE SUD
!                          LLG4- (LON*100) (0-36000) LONGITUDE COIN
!                 GAUSSIEN  LLG1= 1 HEMISPHERE NORD
!                           LLG1= 2 HEMISPHERE SUD
!                           LLG1= 3 GLOBALE
!
!MESSAGES
!         FIN DU FICHIER CA DEBORDE
!         MAUVAISE DIRECTIVE NUMERO=0 FICHIER MS
!         FICHIER INCONNU ROUTINE ECRITUR
!
!MODULES  FSTECR,PUTFLD
!
!APPEL
!         CALL ECRITUR(fld,NPAC,IDAT,DEET,NPAS,NI,NJ,NK,IP1,IP2,IP3,
!                      TYPEV,NOM,ETIKE,GRTYPE,LLG1,LLG2,LLG3,LLG4)
!
!- -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer noenrg,numero,nbrow,numdel,istart
   common/enrege/ noenrg,numero,nbrow,numdel,istart
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer qposition, qitems(16), qnitems
   character*1 qcsepar
   character*16 qcform
   common /idents/  qposition, qitems, qnitems
   common /cidents/ qcsepar,qcform
      integer fltwgtlng(2),fltntimes(2),fltverb(2),fltlist(9,2)
      logical fltoggle(2)
      common /qfilter/ fltoggle,fltwgtlng,fltntimes,fltverb,fltlist
      integer ip1style, dateform
      common /styles/ ip1style, dateform
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
!
!-----------------------------------------------------------------
!
      character *24 chaine
      character *12 cetiket, cetksrt
      character *4 cnomvar
      character *2 ctypvar
      character *1 cgtyp
      character*12 letiksrt
      character*4 lnomvar
      character*2 ltypsrt
      character*72 form1,form2
      integer i, npac,idat,idatv,npas,ni,nj,nk,ip1,ip2,ip3,deet
      real fld(ni,nj,nk)
      real dummy(2)
      integer llg1,llg2,llg3,llg4,iun,istamp,ip3o,ip2o
      integer cdatyp,iopc,ier,gdout, datev
      logical rewrit
      integer gdll,ezgetgdout
      external gdll,ezgetgdout
      integer local_npac
      real*8 delta_t
      if (etiksrt(1) .eq. -1) then
         cetksrt = cetiket
      else
         write(cetksrt,'(3A4)') (etiksrt(i), i=1,nwetiks)
      endif
      if (typesrt .eq. -1) then
         ltypsrt = ctypvar
      else
         write(ltypsrt, '(A2)') typesrt
      endif
!
      if (ip3srt.ne.-1) then
         ip3o=ip3srt
      else
         ip3o=ip3
      endif
!
      if (ip2srt.ne.-1) then
         ip2o=ip2srt
      else
         ip2o=ip2
      endif
!      ltypsrt = ctypsrt
      lnomvar = cnomvar
      letiksrt = cetksrt
      if(npac == 1023) then
        local_npac = npack_orig
      else
        local_npac = npac
      endif
!      print *, npac, npack_orig, local_npac
!
!
!
!
!     SI LE NOM EXISTE DANS LA TABLE BATIT PAR L USAGER ALORS
!     LE fld EST MODIFIE  EX: fld(NI,NJ) = (fld(NI,NJ)+ECART)*FACTEUR
!
      call conver(fld, ni, nj, cnomvar)
!
!
!     filtrage du fld de sortie si le fld n'est pas un stream latlon
!
      if (fltoggle(2)) then
        if (cgtyp.eq.'Y') then
          write (6, *) '(ECRITUR) Impossible de filtrer des flds sur grille Y'
        else
          write (6, *) ' fld FILTRE A L''ECRITURE'
!          call statfld4 (fld,cnomvar,0,'AVANFFLT',ni,nj,1,ni,nj,1,0,0
!     &         ,0)
          call filtre (fld, NI, NJ, 0, fltntimes(2), fltlist(1,2), fltwgtlng(2))
!          call statfld4 (fld,cnomvar,1,'APRESFLT',ni,nj,1,ni,nj,1,0,0
!     &         ,0)
        endif
      endif
!
!
!
      if (printsr)  then
         call imprims(cnomvar,fld,ni,nj)
      endif
!
      iun=lnkdiun(idx_ozsrt)
      if (mode.eq.1) then
!     IWRIT=+1  SORTI(STD,500,R)
         if (compression_level.eq.0) then
            cdatyp = 1
         else
            if (local_npac < 0 .and. local_npac >= -16) then
             cdatyp = 134
            else if (local_npac == -32) then
               cdatyp = 133
            else
               cdatyp = 1
            endif
!          call armn_compress_setlevel(compression_level)
         endif
         if (iwrit.eq.+1) then
            if (.not.message) then
               iopc= fstopc('TOLRNC','DEBUGS',.true.)
            endif
            rewrit = .true.
!
            ier = fstecr(fld,dummy,local_npac,iun,idat,deet,npas,            ni,nj,nk,ip1,ip2o,ip3o,ltypsrt,cnomvar,cetksrt,       &
     &     cgtyp,llg1,llg2,llg3,llg4,cdatyp,rewrit )
         else
            if (.not.message) then
               iopc= fstopc('TOLRNC','DEBUGS',.true.)
            endif
            rewrit = .false.
            ier = fstecr(fld,dummy,local_npac,iun,idat,deet,npas,            ni,nj,nk,ip1,ip2o,ip3o,ltypsrt,cnomvar,cetksrt,       &
     &     cgtyp,llg1,llg2,llg3,llg4,cdatyp,rewrit )
         endif
!
      else
         if (mode.eq.2) then
            write(6,*)            'LES FICHIERS DE TYPE "MS" NE SONT PLUS SUPPORTES'
         else
         endif
!
         if (mode.eq.3.or.mode.eq.4) then
            if (mode.eq.4) then
               cdatyp = 1
               write (chaine, 10) ltypsrt,lnomvar,letiksrt,cgtyp
 10            format(a2,2x,a4,a12,a1,3x)
               write (iun) npac, idat, deet, npas, ni, nj, nk,                ip1, ip2o, ip3o, llg1, llg2, llg3, llg4, cdatyp,     &
     &          chaine
            endif
            write(iun) fld
            if (message) then
               write(6,610)ltypsrt,cnomvar,ip1,ip2o,ip3o,ni,nj,iun
            endif
         else if (mode.eq.5) then
            if (valid) then
               call chk_userdate(datev)
               if (datev .ne. -1) then
                  istamp = datev
               else
                  istamp = idat
               endif
            else
               istamp=0
            endif
            delta_t = deet*npas/3600.0
            call incdatr(idatv,idat,delta_t)
            gdout = ezgetgdout()
            if (gdout .lt.0) gdout = 0
            if (cnomvar(1:2).ne.'LA') then
              if (associated(tmplat)) then
                deallocate(tmplat)
                nullify(tmplat)
                allocate(tmplat(ni,nj))
              endif
            endif
            if (cnomvar(1:2).ne.'LO') then
              if (associated(tmplon)) then
                deallocate(tmplon)
                nullify(tmplon)
                allocate(tmplon(ni,nj))
              endif
            endif
            ier = gdll(gdout, tmplat,tmplon)
            call pgsmwr(2,fld,ni,nj,nk,qcform,qposition,qitems,qcsepar,cnomvar,ctypvar,cetiket,            idat,idatv,dateform,ip1,&
     &ip2,ip3,tmplat,tmplon)
!
         else
            if (message) then
               write(6,*)'FICHIER INCONNU ROUTINE ECRITUR'
            endif
         endif
      endif
!
!
 600  format(2x,' ENREG.ECRIT ',2(a2,'- '),3(i5,'- '),      'TAILLE ',2(i5,'- '),'FICHIER MS ',i4,'   REC=',i4)
 610  format(2x,' ENREG.ECRIT ',2(a2,'- '),3(i5,'- '),      'TAILLE ',2(i5,'- '),'FICHIER SEQUENTIEL',i4)
      return
      end
      subroutine iecritur(fld,npac,idat,deet,npas,ni,nj,nk,ip1,ip2,ip3,ctypvar,cnomvar,cetiket,cgtyp,llg1,llg2,llg3,llg4)
      implicit none
      external conver,fstecr,fclos,memoir,pgsmabt,      imprims,fstopc,messags,fstcvt,putfld
      integer fstopc,fstcvt,fstecr
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer noenrg,numero,nbrow,numdel,istart
   common/enrege/ noenrg,numero,nbrow,numdel,istart
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer qposition, qitems(16), qnitems
   character*1 qcsepar
   character*16 qcform
   common /idents/  qposition, qitems, qnitems
   common /cidents/ qcsepar,qcform
      integer fltwgtlng(2),fltntimes(2),fltverb(2),fltlist(9,2)
      logical fltoggle(2)
      common /qfilter/ fltoggle,fltwgtlng,fltntimes,fltverb,fltlist
      integer ip1style, dateform
      common /styles/ ip1style, dateform
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
!
!-----------------------------------------------------------------
!
      character *24 chaine
      character *12 cetiket, cetksrt
      character *4 cnomvar
      character *2 ctypvar
      character *1 cgtyp
      character*12 letiksrt
      character*4 lnomvar
      character*2 ltypsrt
      character*72 form1,form2
      integer i, npac,idat,idatv,npas,ni,nj,nk,ip1,ip2,ip3,deet, datev
      integer fld(ni,nj,nk)
      real dummy(2)
      integer llg1,llg2,llg3,llg4,iun,istamp,ip3o,ip2o
      integer cdatyp,iopc,ier,gdout,local_npac
      logical rewrit
      integer gdll,ezgetgdout
      external gdll,ezgetgdout
      real*8 delta_t
      if (etiksrt(1) .eq. -1) then
         cetksrt = cetiket
      else
         write(cetksrt,'(3A4)') (etiksrt(i), i=1,nwetiks)
      endif
      if (typesrt .eq. -1) then
         ltypsrt = ctypvar
      else
         write(ltypsrt, '(A2)') typesrt
      endif
!
      if (ip3srt.ne.-1) then
         ip3o=ip3srt
      else
         ip3o=ip3
      endif
!
      if (ip2srt.ne.-1) then
         ip2o=ip2srt
      else
         ip2o=ip2
      endif
!      ltypsrt = ctypsrt
      lnomvar = cnomvar
      letiksrt = cetksrt
!
      if(npac == 1023) then
        local_npac = npack_orig
      else
        local_npac = npac
      endif
      iun=lnkdiun(idx_ozsrt)
      if (mode.eq.1) then
!     IWRIT=+1  SORTI(STD,500,R)
         if (compression_level.eq.0) then
            cdatyp = 2
         else
            if (local_npac < 0 .and. local_npac >= -16) then
               cdatyp = 130
            else
               cdatyp = 2
            endif
         endif
         if (iwrit.eq.+1) then
            if (.not.message) then
               iopc= fstopc('TOLRNC','DEBUGS',.true.)
            endif
            rewrit = .true.
!
            ier = fstecr(fld,dummy,local_npac,iun,idat,deet,npas,ni,nj,nk,ip1,ip2o,ip3o,ltypsrt,cnomvar,cetksrt, cgtyp,llg1,llg2,ll&
     &g3,llg4,cdatyp,rewrit )
         else
            if (.not.message) then
               iopc= fstopc('TOLRNC','DEBUGS',.true.)
            endif
            rewrit = .false.
            ier = fstecr(fld,dummy,local_npac,iun,idat,deet,npas, &
               ni,nj,nk,ip1,ip2o,ip3o,ltypsrt,cnomvar,cetksrt, cgtyp,llg1,llg2,llg3,llg4,cdatyp,rewrit )
         endif
!
      else
         if (mode.eq.2) then
            write(6,*)            'LES FICHIERS DE TYPE "MS" NE SONT PLUS SUPPORTES'
         else if (mode.eq.3.or.mode.eq.4) then
            if (mode.eq.4) then
               cdatyp = 1
               write (chaine, 10) ltypsrt,lnomvar,letiksrt,cgtyp
 10            format(a2,2x,a4,a12,a1,3x)
               write (iun) npac, idat, deet, npas, ni, nj, nk,  ip1, ip2o, ip3o, &
                  llg1, llg2, llg3, llg4, cdatyp, chaine
            endif
            write(iun) fld
            if (message) then
               write(6,610)ltypsrt,cnomvar,ip1,ip2o,ip3o,ni,nj,iun
            endif
         else if (mode.eq.5) then
            if (valid) then
               call chk_userdate(datev)
               if (datev .ne. -1) then
                  istamp = datev
               else
                  istamp = idat
               endif
            else
               istamp=0
            endif
            delta_t = deet*npas/3600.0
            call incdatr(idatv,idat,delta_t)
         else
            if (message) then
               write(6,*)'FICHIER INCONNU ROUTINE ECRITUR'
            endif
         endif
      endif
!
!
 600  format(2x,' ENREG.ECRIT ',2(a2,'- '),3(i5,'- '),      'TAILLE ',2(i5,'- '),'FICHIER MS ',i4,'   REC=',i4)
 610  format(2x,' ENREG.ECRIT ',2(a2,'- '),3(i5,'- '),      'TAILLE ',2(i5,'- '),'FICHIER SEQUENTIEL',i4)
      return
      end
!
!**s/p epaisur  difference entre 2 champs de hauteur
!
      subroutine epaisur(iheur, npar, niveau)
      implicit none
      external ecritur,fstinf,pgsmlir,memoir,fstprm,pgsmabt,      fstcvt,symetri, imprime,loupsou,fstopc,messags,      liraxez
      integer fstinf,pgsmlir,fstprm,fstopc,fstcvt
      integer ezsint, ezqkdef, ezdefset
!
!auteur  p.sarrazin janvier 82  drpn dorval p.q. canada
!revision 4.0.2
!
!   conversion des variables hollerith en caracteres
!   y. chartier -aout 90- drpn dorval quebec.
!
!langage ratfor
!
!objet(epaisur)
!          lire sur fichier standard 2 champs de hauteur
!          prendre la difference entre les 2 champs ecrire le
!          resultat sur fichier approprie(standard,ms,seq)
!
!librairies
!         -source  armnsrc,drpn
!         -objet   pgsmlib,id=armnpjs.
!
!arguments
!   in    iheur   heure des champs (gz)
!   in    npar    nombre de niveaux ( 2)
!   in    niveau  table(2) contenant les 2 niveaux
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!messages
!     record n existe pas sur fichier d entre (epaisur)
!
!modules  fstinf,fstprm,pgsmlir,rgscint,ecritur
!
!appel   via champ
!        call epaisur(iheur,npar,niveau)
!
! - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer dat,deet
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
!
   character *12 cetiket
   character *4 cnomvar
   character *2 ctypvar
   character*1 cigtyp
   real fbidon
   integer iunit
   integer i, niveau(2), iheur, npar, ni, nj, nk, jp1, jp2, jp3, ig1, iopc
   integer ig2,ig3,ig4,irec1,irec2,num1,num2,nn,cnbits,cdatyp,cswa,clng
   integer cdltf,cubc,extra1,extra2,extra3, datev
   integer chkenrpos
      logical sym,symetri
      iunit = 1
!
!     heure ou iheur dans cette routine ne peut-etre -1 heure(tout) pas valide
!
      if (iheur.eq.-1) then
         write(6,*)         'HEURE NE PEUT-ETRE -1(TOUT/ALL) AVEC DIRECTIVE EPAIS'
         call pgsmabt
      endif
!
!  identifier le numero de chaque record avec fstinf
!
      call chk_userdate(datev)
!
!  modification de hollerith a caractere
!
      if (etikent(1) .ne. -1) then
         write(cetiket,'(3A4)') (etikent(i), i=1,nwetike)
      else
         cetiket = '        '
      endif
      if (typeent .ne. -1) then
         write(ctypvar, '(A2)') typeent
      else
         ctypvar = '  '
      endif
      cigtyp = ' '
!
      irec1=fstinf(1,ni,nj,nk,datev,cetiket,niveau(1),iheur,ip3ent,      ctypvar,'GZ  ')
      irec2=fstinf(1,ni,nj,nk,datev,cetiket,niveau(2),iheur,ip3ent,      ctypvar,'GZ  ')
      if (irec2 .lt. 0 .or. irec1 .lt. 0) then
         write(6,*)'RECORD N EXISTE PAS SUR FICHIER D ENTRE (EPAISEUR)'
         return
      endif
      if (nk.gt.1) then
         write(6,*)'**************************************************'
         write(6,*)'         PGSM N ACCEPTE PAS UN          '
         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (EPAISUR)'
         write(6,*)'**************************************************'
         call pgsmabt
      endif
!
!     identifier parametres pour champ 1
!
      ier = fstprm( irec1, dat,deet,npas,ni, nj, nk, cnbits,cdatyp,      jp1,jp2, jp3,ctypvar,cnomvar,cetiket,cigtyp,      ig1,ig2,&
     &ig3,ig4,cswa, clng, cdltf, cubc,       extra1, extra2, extra3)
      if (ier .lt. 0) then
         write(6,*)' IER = FSTPRM NEGATIF VOIR EPAISUR'
      endif
!
!     verifier si grille gaussienne ni doit etre pair
!
      if (cigtyp.eq.'G'.and.mod(ni,2).ne.0) then
         call messags(ni)
      endif
!
!     lire champ no 1
!
      allocate(tmpif1(ni,nj))
      if (.not.message) then
         iopc= fstopc('TOLRNC','DEBUGS',.true.)
      endif
      call chk_userdate(datev)
!
      num1 = pgsmlir(tmpif1,1,ni,nj,nk,datev,cetiket,jp1,jp2,jp3,ctypvar, 'GZ  ', cigtyp)
!
      if (printen) then
         call imprime(cnomvar,tmpif1,ni,nj)
      endif
!
!     identifier parametres pour champ 2
!
      ier = fstprm( irec2, dat,deet,npas,ni, nj, nk, cnbits,cdatyp,      jp1,jp2, jp3,ctypvar,cnomvar,cetiket,cigtyp,       ig1,ig2&
     &,ig3,ig4, cswa,clng,cdltf,cubc,extra1,extra2,extra3)
      if (ier .lt. 0) then
         write(6,*)' IER = FSTPRM NEGATIF VOIR EPAISUR'
      endif
!
!     verifier si grille gaussienne ni doit etre pair
!
      if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  then
         call messags(ni)
      endif
!
!     lire champ 2
!
      allocate(tmpif2(max0(li,ni),max0(nj,lj)))
      if (.not.message)  then
         iopc= fstopc('TOLRNC','DEBUGS',.true.)
      endif
      call chk_userdate(datev)
      num2 = pgsmlir(tmpif2,1,ni,nj,nk,datev,cetiket,jp1,jp2,jp3,  ctypvar, 'GZ  ', cigtyp)
      if (printen)  call imprime(cnomvar,tmpif1,ni,nj)
!
!  difference entre les deux champs
!
      nn = ni*nj
      call loupsou(tmpif1,tmpif2,nn)
!
!     interpolation horizontale
!
      if (cigtyp == 'A' .or. cigtyp == 'B' .or. cigtyp == 'G') then
         if (ig1 /= 0) sym = symetri(cnomvar)
      endif
      if (cgrtyp.eq.'*') then
         ier = chkenrpos(1,2,ig1,ig2,ig3)
         call ecritur(tmpif1,npack,dat,deet,npas,ni,nj,nk,         niveau(1),niveau(2),iheur,         ctypvar,'DZ  ',cetiket,cigtyp&
     &,ig1,ig2,ig3,ig4)
      else
         gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
         ier = ezdefset(gdout, gdin)
         ier = ezsint(tmpif2, tmpif1)
!
!
         jp1 = niveau(1)
         jp2 = niveau(2)
         jp3 = iheur
!
!     ecrire sur fichier standard,ms,sequentiel
!
         call ecritur(tmpif2,npack,dat,deet,npas,li,lj,nk,jp1,jp2,jp3,         ctypvar,'DZ  ',cetiket,cgrtyp,lg1,lg2,lg3,lg4)
      endif
      deallocate(tmpif1)
      deallocate(tmpif2)
!
      return
      end
   integer function fst_get_mask_key(mask_key, fld_key, mask_flags, iun) result(status)
   implicit none
   integer mask_key, fld_key, mask_flags, iun, ier
   integer fstprm, fstinf
   external fstprm, fstinf
   character(len=4)  :: fld_nomvar, mask_nomvar
   character(len=2)  :: fld_typvar, mask_typvar
   character(len=12) :: fld_etiket, mask_etiket
   character(len=1)  :: fld_grtyp,  mask_grtyp
   integer :: fld_dateo, fld_deet, fld_npas, fld_ni, fld_nj, fld_nk, fld_nbits, fld_datyp
   integer :: fld_ip1, fld_ip2, fld_ip3, fld_ig1, fld_ig2, fld_ig3, fld_ig4
   integer :: fld_lng, fld_dltf, fld_ubc, fld_swa, fld_datev, fld_extra2, fld_extra3
   integer :: mask_dateo, mask_deet, mask_npas, mask_ni, mask_nj, mask_nk, mask_nbits, mask_datyp
   integer :: mask_ip1, mask_ip2, mask_ip3, mask_ig1, mask_ig2, mask_ig3, mask_ig4
   integer :: mask_lng, mask_dltf, mask_ubc, mask_swa, mask_datev, mask_extra2, mask_extra3
   integer :: allones
   integer ip_allones
   integer :: sorte = 3
   integer :: mode  = 1
   character(len=32) :: ip_string
   logical :: flag = .false.
!    allones = int(.not.ishft(-1,28))
  allones = -1
   ier = fstprm(fld_key, fld_dateo, fld_deet, fld_npas, fld_ni, fld_nj, fld_nk, &
            fld_nbits, fld_datyp, fld_ip1, fld_ip2, fld_ip3, fld_typvar, fld_nomvar, fld_etiket, &
            fld_grtyp, fld_ig1, fld_ig2, fld_ig3, fld_ig4, fld_swa, fld_lng, fld_dltf, fld_ubc, &
            fld_datev, fld_extra2, fld_extra3)
   if (fld_typvar(2:2) /= '@') then
      print *, ' (fst_get_mask_key) This is not a masked field', fld_nomvar, fld_typvar
      status = -1
      return
   endif
   if (fld_typvar(1:1) == '@') then
      print *, ' (fst_get_mask_key) This is a mask field', fld_nomvar, fld_typvar
      status = -1
      return
   endif
   mask_datev  = fld_datev
   mask_ip1    = fld_ip1
   mask_ip2    = fld_ip2
   mask_ip3    = fld_ip3
   mask_nomvar = fld_nomvar
   mask_typvar = '@@'
   mask_key = fstinf(iun, mask_ni,mask_nj, mask_nk, mask_datev, mask_etiket, mask_ip1, &
      mask_ip2, mask_ip3, mask_typvar, mask_nomvar)
   if (mask_key >= 0) then
      ier = fstprm(mask_key, mask_dateo, mask_deet, mask_npas, mask_ni, mask_nj, mask_nk, &
               mask_nbits, mask_datyp, mask_ip1, mask_ip2, mask_ip3, mask_typvar, mask_nomvar, mask_etiket, &
               mask_grtyp, mask_ig1, mask_ig2, mask_ig3, mask_ig4, mask_swa, mask_lng, mask_dltf, mask_ubc, &
               mask_datev, mask_extra2, mask_extra3)
      if (mask_ni == fld_ni .and. mask_nj == fld_nj .and. mask_nk == fld_nk .and. &
          mask_grtyp == fld_grtyp .and. mask_ig1 == fld_ig1 .and. mask_ig2 == fld_ig2 .and. &
          mask_ig3 == fld_ig3 .and. mask_ig4 == fld_ig4) then
      status = 0
      return
      endif
   endif
   mask_ip1    = allones
   mask_key = fstinf(iun, mask_ni,mask_nj, mask_nk, mask_datev, mask_etiket, mask_ip1, &
      mask_ip2, mask_ip3, mask_typvar, mask_nomvar)
   if (mask_key >= 0) then
      ier = fstprm(mask_key, mask_dateo, mask_deet, mask_npas, mask_ni, mask_nj, mask_nk, &
               mask_nbits, mask_datyp, mask_ip1, mask_ip2, mask_ip3, mask_typvar, mask_nomvar, mask_etiket, &
               mask_grtyp, mask_ig1, mask_ig2, mask_ig3, mask_ig4, mask_swa, mask_lng, mask_dltf, mask_ubc, &
               mask_datev, mask_extra2, mask_extra3)
      if (mask_ni == fld_ni .and. mask_nj == fld_nj .and. mask_nk == fld_nk .and. &
          mask_grtyp == fld_grtyp .and. mask_ig1 == fld_ig1 .and. mask_ig2 == fld_ig2 .and. &
          mask_ig3 == fld_ig3 .and. mask_ig4 == fld_ig4) then
      status = 0
      return
      endif
   endif
   mask_nomvar = '@@@@'
   mask_ip1 = fld_ip1
   mask_key = fstinf(iun, mask_ni,mask_nj, mask_nk, mask_datev, mask_etiket, mask_ip1, &
      mask_ip2, mask_ip3, mask_typvar, mask_nomvar)
   if (mask_key >= 0) then
      ier = fstprm(mask_key, mask_dateo, mask_deet, mask_npas, mask_ni, mask_nj, mask_nk, &
               mask_nbits, mask_datyp, mask_ip1, mask_ip2, mask_ip3, mask_typvar, mask_nomvar, mask_etiket, &
               mask_grtyp, mask_ig1, mask_ig2, mask_ig3, mask_ig4, mask_swa, mask_lng, mask_dltf, mask_ubc, &
               mask_datev, mask_extra2, mask_extra3)
      if (mask_ni == fld_ni .and. mask_nj == fld_nj .and. mask_nk == fld_nk .and. &
          mask_grtyp == fld_grtyp .and. mask_ig1 == fld_ig1 .and. mask_ig2 == fld_ig2 .and. &
          mask_ig3 == fld_ig3 .and. mask_ig4 == fld_ig4)  then
      status = 0
      return
      endif
   endif
   mask_ip1 = allones
   mask_key = fstinf(iun, mask_ni,mask_nj, mask_nk, mask_datev, mask_etiket, mask_ip1, &
      mask_ip2, mask_ip3, mask_typvar, mask_nomvar)
   if (mask_key >= 0) then
      ier = fstprm(mask_key, mask_dateo, mask_deet, mask_npas, mask_ni, mask_nj, mask_nk, &
               mask_nbits, mask_datyp, mask_ip1, mask_ip2, mask_ip3, mask_typvar, mask_nomvar, mask_etiket, &
               mask_grtyp, mask_ig1, mask_ig2, mask_ig3, mask_ig4, mask_swa, mask_lng, mask_dltf, mask_ubc, &
               mask_datev, mask_extra2, mask_extra3)
      if (mask_ni == fld_ni .and. mask_nj == fld_nj .and. mask_nk == fld_nk .and. &
          mask_grtyp == fld_grtyp .and. mask_ig1 == fld_ig1 .and. mask_ig2 == fld_ig2 .and. &
          mask_ig3 == fld_ig3 .and. mask_ig4 == fld_ig4)  then
      status = 0
      return
      endif
   endif
   print *, ' (fst_get_mask_key) Associated mask not found'
   status = -1
   return
   end function fst_get_mask_key
!====================================================================================================
!  Reserve pour provisions futures
!    mask_ip3 = allones
!
!    mask_key = fstinf(iun, mask_ni,mask_nj, mask_nk, mask_datev, mask_etiket, mask_ip1, &
!       mask_ip2, mask_ip3, mask_typvar, mask_nomvar)
!
!    if (mask_key >= 0) then
!       if (mask_ni == fld_ni .and. mask_nj == fld_nj .and. mask_nk == fld_nk .and. &
!           mask_grtyp == fld_grtyp .and. mask_ig1 == fld_ig1 .and. mask_ig2 == fld_ig2 .and. &
!           mask_ig3 == fld_ig3 .and. mask_ig4 == fld_ig4)  then
!       status = 0
!       return
!       endif
!    endif
!
!    mask_ip2 = allones
!
!    mask_key = fstinf(iun, mask_ni,mask_nj, mask_nk, mask_datev, mask_etiket, mask_ip1, &
!       mask_ip2, mask_ip3, mask_typvar, mask_nomvar)
!
!    if (mask_key >= 0) then
!       if (mask_ni == fld_ni .and. mask_nj == fld_nj .and. mask_nk == fld_nk .and. &
!           mask_grtyp == fld_grtyp .and. mask_ig1 == fld_ig1 .and. mask_ig2 == fld_ig2 .and. &
!           mask_ig3 == fld_ig3 .and. mask_ig4 == fld_ig4)  then
!       status = 0
!       return
!       endif
!    endif
      subroutine fillcoord(lat,lon)
      implicit none
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      real lat(*),lon(*)
      integer i
      do i=1,ncoords
         lat(i) = coordll(i,1)
         lon(i) = coordll(i,2)
         enddo
      return
      end
!
!**S/P GRIGAUS   CALCUL LAT LONG DE CHAQUE PT D'UNE GRILLE GAUSSIENNE
!
      subroutine grigaus(nni,nnj,nhem)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRIGAUS)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE GAUSSIENNE LONGITUDE EQUIDISTANTE
!          LATITUDE GAUSSIENNE
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
!
      external memoir,pgsmabt,grgg,messags
      external ezqkdef, gdll
      integer ezqkdef, gdll
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
!
      integer nni,nnj,nhem,ier,nroot
!
!   RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
      allocate(tmplat(nni,nnj))
      allocate(tmplon(nni,nnj))
!
      lg1=nhem
!
      li=nni
      lj=nnj
!
      if (lg1.lt.0 .or.lg1.gt.2) then
         write(6,*)'LG1 DOIT ETRE  GLOBAL,NORD,SUD   GRILLE(GAUSS.....'
         call pgsmabt
      endif
!
      cgrtyp='G'
      lg2=0
      lg3=0
      lg4=0
!
      if (lg1.eq.0) then
         nroot=nnj
      else
         nroot=nnj*2
      endif
!     BUFL(IROOT)  -TABLE DE RACINES DES POLYNOMES DE LEGENDRE
!     UTILISER PAR GRIGAUS POUR LE CALCUL DES LATITUDES
!     LONGITUDES DANS GRGG ROUTINE
!
      allocate(tmproot(nroot))
      tmproot(1)= 100.0
      gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
      ier = gdll(gdout, tmplat, tmplon)
!
      deallocate(tmproot)
!
      return
      end
!
!**S/P GRIGEF   CALCUL LAT LONG DE CHAQUE PT D'UNE GRILLE "E"
!
      subroutine grigef(it,nni,nnj,xlat1,xlon1,xlat2,xlon2)
      implicit none
!
!AUTEUR   - y. chartier april 94
!
!LANGAGE - RATFOR
!
!     OBJET(GRISTDB)
!     CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!     DE LA GRILLE DE SORTIE STANDARD 'e' LAT ET LONG EQUIDISTANT
!     LONGITUDE ZERO ET 360 PRESENT.
!
!
!------------------------------------------------------
      external memoir,pgsmabt,grll,lastcol,messags
      external ezqkdef, gdll
      integer ezqkdef, gdll
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
      integer nni,nnj,it,ier
      real xlat1,xlon1,xlat2,xlon2
      real xla0,xlo0,dlat,dlon,valeur
!
      li=nni
      lj=nnj
!
!
!     RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
      allocate(tmplat(nni,nnj))
      allocate(tmplon(nni,nnj))
      allocate(tmplatg(nni,nnj))
      allocate(tmplong(nni,nnj))
!
      cgrtyp='E'
      call cxgaig(cgrtyp,lg1,lg2,lg3,lg4,xlat1,xlon1,xlat2,xlon2)
      gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
      ier = gdll(gdout, tmplat, tmplon)
      return
      end
!
!**   S/P GRIGRIB  CALCUL LATITUDE LONGITUDE DE CHAQUE PT D'UNE GRILLE GRIB
!
      subroutine grigrib(ig1,ig2,ig3,ig4)
      implicit none
!
!     AUTEUR   -  Y. CHARTIER DRPN DORVAL MAI 1996
!
!
!     LANGAGE - RATFOR
!
!     OBJET(GRIGRIB)
!     CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!     DE LA GRILLE DE SORTIE POLAIRE STEREOGRAPHIQUE
!
!
!     LIBRAIRIES
!     -SOURCE  ARMNSRC,DRPN
!     -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
      external memoir,pgsmabt,grps,cigaxg,cxgaig,messags
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
      character*1 gtyout
      real xg(20)
      integer nni,nnj,ihm,hem,ier,npts
      integer ig1,ig2,ig3,ig4
      integer tmpig1, tmpig2, tmpig3, tmpig4
      external ezqkdef, gdll
      integer ezqkdef, gdll
      real, dimension(:,:), allocatable :: x, y
!
!     RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
      cgrtyp = '!'
      call igaxg95(gtyout,xg,15,cgrtyp,ig1,ig2,ig3,ig4)
      if (gtyout.eq.'H') then
         nni = nint(xg(8))
         nnj = nint(xg(9))
      endif
      allocate(tmplat(nni,nnj))
      allocate(tmplon(nni,nnj))
      allocate(tmplon(nni,nnj))
      allocate(x(nni,nnj))
      allocate(y(nni,nnj))
!
!
      li=nni
      lj=nnj
      lg1 = ig1
      lg2 = ig2
      lg3 = ig3
      lg4 = ig4
!
!
      npts = nni*nnj
      call ez_llflamb(tmplat,tmplon,x,y,npts,cgrtyp,ig1,ig2,ig3,ig4)
!
!      call cxgaig('L', tmpig1, tmpig2, tmpig3, tmpig4,0.,0.,1.0,1.0)
!      gdout = ezgdef(npts,1,'Y','L', tmpig1, tmpig2, tmpig3, tmpig4,
!     $        tmplon,tmplat)
      gdout = ezqkdef(nni,nnj,cgrtyp,ig1,ig2,ig3,ig4,1)
      ier = gdll(gdout, tmplat, tmplon)
      deallocate(x)
      deallocate(y)
      return
      end
!
!**S/P GRILLE   DETERMINE LA SORTE DE GRILLE DEMANDE PAR USAGER
!
      subroutine grille2(it,p1,p2,p3,p4,p5,p6,p7,p8)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 82 DRPN DORVAL P.Q. CANADA
!           MODIFIER JANVIER 87 P.SARRAZIN DORVAL P.Q. CANADA
!           MODIFIER MAI 87 P. SARRAZIN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRILLE)
!          VERIFI  LE NOMBRE D ARGUMENTS LA VALEUR DU PREMIER ARGUMENT
!          DETERMINE LE TYPE DE GRILLE DE SORTIE, RETOURNE LA MEMOIRE
!          POUR IXLAT INITIALISE A ZERO DANS LE MAIN PGSM.
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN   IT    - 1= GRILLE STANDARD
!                2= GRILLE LAT-LON
!                3= GRILLE P.S. POLAIRE STEREOGRAPHIQUE
!                4= TAPE4 FICHIER CONTENANT LATITUDES LONGITUDES
!                 OU COORDONNEES EST-OUEST OU NORD-SUD
!                5- STDB GRILLE STANDARD 'B'
!                6= GRILLE GAUSSIENNE
!                7= GRILLE TAPE1 1 REC LAT,Y  1 REC LON,X
!                8= GRILLE TAPE2 ECRIT 1 REC LAT,Y 1 REC LON,X SUR TAPE2
!                9= GRILLE GEF
!               10= GRILLE GRIB
!               11= COORDONNEES LOCALES
!
!MESSAGES
!         DIRECTIVE SORTIE DOIT-ETRE APPELE AVANT DIRECTIVE GRILLE
!         GRILLE INCONUE(GRILLE)
!
!MODULES PGSMABT,GRILSTD,GRLALON,GRILLPS,GRILTP4,GRIGAUS,GRISTDB
!MODULES GRITP12,MEMOIR,LLFXY
!------------------------------------------------------
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
      external pgsmabt,grilstd,grlalon,grillps,griltp4,grigaus,gristdb,messags
      external gritp12,memoir
!     external verlalo, chklalo
!
      integer ier,it,p1,p2,p3,p4,p5,p6,p7,p8
!
!
      if (iset.eq.-2) then
         write(6,*)       'DIRECTIVE SORTIE DOIT-ETRE APPELEE AVANT DIRECTIVE GRILLE'
         call pgsmabt
      endif
!
      cgtypxy='L'
!
      if (associated(tmplat)) deallocate(tmplat)
      if (associated(tmplon)) deallocate(tmplon)
      if (associated(tmplatg)) deallocate(tmplatg)
      if (associated(tmplong)) deallocate(tmplong)
!
      if (it.lt.0 .or. it.gt.14) then
         write(6,*)'GRILLE INCONNUE (GRILLE)'
         call pgsmabt
      endif
!
      if (it.eq.gr_a) then
!
!
!     CALCUL GRILLE STD
!     ----------
!
         if (ngr.eq.4) then
            call grilstd(p1,p2,p3)
!
!     NI =P1  NOMBRE DE POINTS EST-OUEST
!     NJ =P2  NOMBRE DE POINTS NORD SUD
!     LG1=P3  0=GLOBAL;  1=H. NORD;  2=H. SUD
!
         else
            write(6,*) ' MAUVAIS APPEL GRILLE(STD,NI,NJ,NORD/SUD/GLOBAL)'
            call pgsmabt
         endif
!...................................................................
!
!     CALCUL GRILLE LATLON
!     -------------
      elseif (it.eq.gr_latlon) then
         if (ngr.eq.7) then
            call grlalon(p1,p2,p3,p4,p5,p6)
!
!     NI=P1      NOMBRE DE POINTS EST-OUEST
!     NJ=P2      NOMBRE DE POINTS NORD-SUD
!     XLAT0=P3   1 IERE LAT EN BAS A GAUCHE
!     XLON0=P4   1 IERE LONG EN BAS A GAUCHE
!     DLAT=P5    ESPACEMENT ENTRE CHAQUE LATITUDE
!     DLON=P6    ESPACEMENT ENTRE CHAQUE LONGITUDE
!
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(LATLON,NI,NJ,XLAT0,XLON0,DLAT,DLON)'
            call pgsmabt
         endif
!
!...................................................................
      elseif (it.eq.gr_ps) then
!     CALCUL GRILLE PS
!     ---------
!
         if (ngr.eq.7) then
            call grillps(p1,p2,p3,p4,p5,p6,1)
         else if (ngr.eq.8) then
            call grillps(p1,p2,p3,p4,p5,p6,p7)
!
!     NI=P1.........NOMBRE DE POINTS DANS DIRECTION EST-OUEST
!     NJ=P2.........NOMBRE DE POINTS DANS DIRECTION NORD-SUD
!     PI=P3.........POSITION DU POLE DIRECTION EST-OUEST (GRID POINT)
!     PJ=P4.........POSITION DU POLE DIRECTION NORD SUD (GRID POINT)
!     D60=P5........DISTANCE ENTRE 2 GRID POINTS EN METRES
!     DGRW=P6.......ORIENTATION DE LA GRILLE PAR RAPPORT A GREENWICH
!     NORD/SUD=P7...HEMISPHERE NORD/SUD
!
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(PS,NI,NJ,PI,PJ,D60,DGRW,NORD/SUD)'
            call pgsmabt
         endif
!
!...................................................................
!
!   CALCUL GRILLE TAPE4
!          ------------
      elseif (it.eq.gr_tape4) then
         if (ngr.eq.6) then
            call griltp4(p1,p2,p3,p4,p5)
         else if (ngr.eq.3) then
            call griltp4(p1,p2,-1,-1,-1)
!
!     NI=P1.....NOMBRE DE POINTS EST-OUEST
!     NJ=P2.....NOMBRE DE POINTS NORD-SUD
!     IP1=P3....VALEUR DE IP1 TRANSFER DANS IG1 POUR TAPE2
!     IP2=P4....VALEUR DE IP2 TRANSFER DANS IG2 POUR TAPE2
!     IP3=P5....VALEUR DE IP3 TRANSFER DANS IG3 POUR TAPE2
!
         else
            write(6,*) 'MAUVAIS APPEL  GRILLE(TAPE4,NI,NJ [,IP1,IP2,IP3] )'
         write(6,*)' CETTE DIRECTIVE CONTIENT 3 OU 6 ARGUMENTS'
         call pgsmabt
      endif
!
      elseif (it.eq.gr_g) then
!...................................................................
!
!     CALCUL GRILLE GAUSSIENNE
!     -----------------
!
         if (ngr.eq.4) then
            if (mod(p1,2).ne.0) then
               write(6,*)' ON NE PEUT PRODUIRE UN CHAMP GAUSSIEN'
               write(6,*)'  AVEC UN NOMBRE DE LONGITUDES IMPAIRS'
               call pgsmabt
            endif
            call grigaus(p1,p2,p3)
!
!     NI =P1  NOMBRE DE POINTS EST-OUEST
!     NJ =P2  NOMBRE DE POINTS NORD SUD
!     LG1=P3  0=GLOBAL;  1=H. NORD;  2=H. SUD
!
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(GAUSS,NI,NJ,NORD/SUD/GLOBAL)'
            call pgsmabt
         endif
      elseif (it.eq.gr_b) then
!......................................................................
!
!     CALCUL GRILLE STDB
!     -----------
!
         if (ngr.eq.4) then
            call gristdb(p1,p2,p3)
!
!     NI =P1  NOMBRE DE POINTS EST-OUEST
!     NJ =P2  NOMBRE DE POINTS NORD SUD
!     LG1=P3  0=GLOBAL;  1=H. NORD;  2=H. SUD
!
         else
            write(6,*)' MAUVAIS APPEL GRILLE(STDB,NI,NJ,NORD/SUD/GLOBAL) '
            call pgsmabt
         endif
      elseif (it.eq.gr_tape1.or.it.eq.gr_tape2.or.it.eq.gr_stations) then
!...................................................................
!
!     CALCUL GRILLE TAPE1/TAPE2 LAT-LON OU X-Y
!
!     SI IT=7 TAPE1 ENTRE
!     SI IT=8 TAPE2 SORTI
!
         if (ngr.eq.4) then
            call gritp12(it,p1,p2,p3)
!
!     IT=7......LIRE TAPE 1 ECRIT SUR TAPE2
!     IT=8......LIRE TAPE 2 ECRIT SUR TAPE2
!     IP1=P1....IDENTIFICATION DU RECORD VALEUR MAX=2047
!     IP2=P2....IDENTIFICATION DU RECORD VALEUR MAX=2047
!     IP3=P3....IDENTIFICATION DU RECORD VALEUR MAX=2047
!
         else
            write(6,*)'MAUVAIS APPEL GRILLE(TAPE1/TAPE2,IP1,IP2,IP3)'
            call pgsmabt
         endif
!
      elseif (it.eq.gr_comme) then
         if (ngr.eq.9) then
            call comme(p1,p2,p3,p4,p5,p6,p7,p8)
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(COMME,FENTREE/FSORTIE,NOMVAR,TYPVAR,DATEV,IP1,IP2,IP3,ETIKET)'
            call pgsmabt
         endif
      elseif (it.eq.gr_grib) then
         if (ngr.eq.5) then
            call grigrib(p1,p2,p3,p4)
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(GRIB,IG1,iG2,IG3,IG4)'
            call pgsmabt
         endif
      elseif (it.eq.15) then
         if (ngr.eq.7) then
            call grigef(it,p1,p2,p3,p4,p5,p6)
!
!     IT=7......LIRE TAPE 1 ECRIT SUR TAPE2
!     IT=8......LIRE TAPE 2 ECRIT SUR TAPE2
!     IP1=P1....IDENTIFICATION DU RECORD VALEUR MAX=2047
!     IP2=P2....IDENTIFICATION DU RECORD VALEUR MAX=2047
!     IP3=P3....IDENTIFICATION DU RECORD VALEUR MAX=2047
!
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(E,NI,NJ,XLAT1,XLON1,XLAT2,XLON2)'
            call pgsmabt
         endif
      elseif (it.eq.gr_stereo) then
         if(ngr.eq.7) then
            call gristereo(p1,p2,p3,p4,p5,p6)
         else
            write(6,*)            'MAUVAIS APPEL GRILLE(STEREO,NI,NJ,D60,DGRW,CLAT,CLON)'
            call pgsmabt
         endif
      endif
      return
      end
!
!**S/P GRILLPS  CALCUL LATITUDE LONGITUDE DE CHAQUE PT D'UNE GRILLE P.S.
!
      subroutine grillps(nni,nnj,pi,pj,d60,dgrw,hem)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRILLPS)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE POLAIRE STEREOGRAPHIQUE
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
       external memoir,pgsmabt,grps,cigaxg,cxgaig,messags
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
       external ezqkdef, gdll
       integer ezqkdef, gdll
       integer nni,nnj,ihm,hem,ier
       real pi,pj,d60,dgrw,pp1,pp2,pp3,pp4
!
!   RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
       allocate(tmplat(nni,nnj))
       allocate(tmplon(nni,nnj))
!
!
       li=nni
       lj=nnj
!
!
       ihm=hem
       if (ihm.eq.1) then
          cgrtyp='N'
       else
          cgrtyp='S'
       endif
       call cxgaig(cgrtyp,lg1,lg2,lg3,lg4,pi,pj,d60,dgrw)
       call cigaxg(cgrtyp,pp1,pp2,pp3,pp4,lg1,lg2,lg3,lg4)
      if (ihm.lt.1.or.ihm.gt.2) then
         write(6,*)'GRILLE P.S. CODE D HEMISPHERE DOIT-ETRE NORD OU SUD'
         call pgsmabt
      endif
!
      gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
      ier = gdll(gdout, tmplat, tmplon)
!
!
!
      return
      end
!
!**S/P GRILSTD   CALCUL LATITUDE LONGITUDE DE CHAQUE PT D'UNE GRILLE STD
!
      subroutine grilstd(nni,nnj,hem)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRILSTD)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE STANDARD INTERVAL REGULIER MAIS DECALE
!          1/2 POINT DU POLE ET DE L'EQUATEUR
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
       external memoir,pgsmabt,grll,messags
       external ezqkdef, gdll
       integer ezqkdef, gdll
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
       integer nni,nnj,hem,ier
       real xla0,xlo0,dlat,dlon
!
       li=nni
       lj=nnj
!
!
!   RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
       allocate(tmplat(nni,nnj))
       allocate(tmplon(nni,nnj))
!
!
       lg1=hem
       cgrtyp='A'
       lg2=0
       lg3=0
       lg4=0
       gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
       ier = gdll(gdout, tmplat, tmplon)
!
!
       if (lg1.ne.0.and.lg1.ne.1.and.lg1.ne.2) then
          write(6,*)'GRILLE(STD........ DOIT ETRE NORD/SUD/GLOBAL'
          call pgsmabt
       endif
!
       return
       end
!     
!**S/P GRILTP4   LIRE LAT LONG DE CHAQUE PT D'UNE GRILLE TYPE "X" OU "Y"
!
      subroutine griltp4(nni,nnj,ip1,ip2,ip3)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRILTP4)
!          LIRE LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE TAPE4 LAT LON ET ECRIRE LES DEUX
!          RECORDS LAT LONG SUR FICHIER STANDARD TYPE "X" SI
!          IP1,IP2,IP3 SONT DEFINIS TYPE "Y" ET SERONT TRANSFERRE DANS
!          IG1,IG2,IG3   IG4=0
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
      integer nni,nnj,ip1,ip2,ip3
      write(6,*)  'CE TYPE DE GRILLE N EST PAS SUPPORTE DANS CETTE VERSION DE PGSM'
      call pgsmabt
      return
      end
!
!**S/P GRISTDB   CALCUL LAT LONG DE CHAQUE PT D'UNE GRILLE STD "B"
!
      subroutine gristdb(nni,nnj,hem)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRISTDB)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE STANDARD 'B' LAT ET LONG EQUIDISTANT
!          LONGITUDE ZERO ET 360 PRESENT.
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
      external memoir,pgsmabt,grll,lastcol,messags
      external ezqkdef, gdll
      integer ezqkdef, gdll
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
      integer nni,nnj,hem,ier
      real xla0,xlo0,dlat,dlon,valeur
!
      li=nni
      lj=nnj
!
!
!   RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
      allocate(tmplat(nni,nnj))
      allocate(tmplon(nni,nnj))
!
      lg1=hem
      cgrtyp='B'
      lg2=0
      lg3=0
      lg4=0
!
      if (lg1.lt.0 .or. lg1.gt.2) then
         write(6,*)'GRILLE(STDB.....DOIT-ETRE  GLOBAL,NORD,SUD '
         call pgsmabt
      endif
!
      gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
      ier = gdll(gdout, tmplat, tmplon)
!
!
      return
      end
!
!**S/P GRISTEREO  CALCUL LATITUDE LONGITUDE DE CHAQUE PT D'UNE GRILLE P.S.
!
      subroutine gristereo(nni,nnj,d60,dgrw,clat,clon)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRISTEREO)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE STEREOGRAPHIQUE GENERALISEE
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
       external memoir,pgsmabt,grps,cigaxg,cxgaig,messags
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
       external ezqkdef, gdll
       integer ezqkdef, gdll
       integer nni,nnj,ihm,hem,ier
       real pi,pj,d60,dgrw,pp1,pp2,pp3,pp4,clat,clon
!
!   RESERVER MEMOIR POUR LATITUDE ET LONGITUDE
!
       allocate(tmplat(nni,nnj))
       allocate(tmplon(nni,nnj))
!
!
       li=nni
       lj=nnj
!
!
       cgrtyp='T'
       call cxgaig(cgrtyp,lg1,lg2,lg3,lg4,d60,dgrw,clat,clon)
       call cigaxg(cgrtyp,pp1,pp2,pp3,pp4,lg1,lg2,lg3,lg4)
      gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
      ier = gdll(gdout, tmplat, tmplon)
!
!
!
      return
      end
!
!**S/P GRITP12   LIRE LAT LONG DE CHAQUE PT D'UNE GRILLE TYPE "Y" OU "Z"
!
      subroutine gritp12(it,ip1,ip2,ip3)
      implicit none
      integer it,ip1,ip2,ip3
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRITP12)
!          LIRE UN REC DE LAT ET UN REC DE LONG
!          POUR LA GRILLE DE TYPE "Y" OU "Z".
!          Y = LISTE DE LAT-LON(NI,NJ) OU X-Y(NI,NJ)
!          Z = COLONNE DE LAT(NJ) OU RANGEE DE LONG(NI)
!              OU COLONNE DE Y(NJ)   RANGEE DE X(NI)
!          COORDONNEE X,Y SONT POLAIRE STEREOGRAPHIQUE
!
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
!
      external pgsmabt,memoir,fstprm,cigaxg,messags,fstcvt
      external imprime,ecritur,fstopc,fstinf,fstlir,conlal2
      integer fstprm,fstopc,fstinf,fstlir,fstcvt
      integer ezqkdef, ezgdef,ezgfstp, ezgxprm, gdgaxes, gdll, chkenrpos
      external ezqkdef, ezgdef,ezgfstp, ezgxprm, gdgaxes, gdll, chkenrpos
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer ig1la,ig2la,ig3la,ig4la,nilo,njlo,nklo
      integer ig1lo,ig2lo,ig3lo,ig4lo
      common /tp12ig/ ig1la,ig2la,ig3la,ig4la,nilo,njlo,nklo,ig1lo,ig2lo,ig3lo,ig4lo
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
      character *12 etikx
      character *4 nomx
      character *2 ctpvrla,ctpvrlo,cgtyplo,cgtypla
      character*2 grref
      integer i
      integer swa, lng, dltf, ubc, extra1, extra2, extra3
      integer iunit,irecla,ireclo,irecyy,ni,nj,nk
      integer dateo, deet, npas
      integer numla,numlo,ier,iopc
      integer ig1, ig2, ig3, ig4, ig1ref, ig2ref, ig3ref, ig4ref
      integer ip1x, ip2x, ip3x
      character*4 nomvarx, nomvary, nomu
      character*12  etiketx, etikety, etiku
      character*2   typvarx, typvary, typvaru
      integer nix, njx, niu, niy, njy, nju, nkx, nky, nku
      integer dateox, deetx,npasx, nbitsx, datypx
      integer dateou, deetu,npasu, nbitsu, datypu
!
      real pidum,pjdum,d60dum
      integer lip1, lip2, lip3
      logical grille_z, grille_u
!
!
      if (it.eq.gr_tape2) then
         iunit=lnkdiun(idx_ozsrt)
      else
         iunit=lnkdiun(1)
      endif
!
!  VERIFICATION DES PARAMETRES DES 2 RECORDS
!  SI LU SUR UNIT 1 ON ECRIT LES 2 RECORDS SUR UNIT 2
!
      nix = -1
      njx = -1
      nkx = -1
      niy = -1
      njy = -1
      nky = -1
      niu = -1
      nju = -1
      nku = -1
      grille_z = .false.
      grille_u = .false.
      if (it.ne.gr_stations) then
!         if (it.ne.gr_tape2) then
            ier = chkenrpos(lnkdiun(1),lnkdiun(idx_ozsrt),ip1,ip2,ip3)
            if (ier.lt.0) then
               print *, '<gritp12> enregistrements positionnels absents!'
               print *, '          impossible de continuer...'
!               call exit(13)
               call pgsmabt
            elseif (ier == 0) then
	            grille_z = .true.
	         else
	            grille_u = .true.
!	         endif
         endif
         if (mode.eq.1) then
	         if (grille_z) then
               ireclo = fstinf(iunit, nix, njx, nkx, -1, '            ', ip1, ip2, ip3, '  ', '>>  ')
               irecla = fstinf(iunit, niy, njy, nky, -1, '            ', ip1, ip2, ip3, '  ', '^^  ')
               cgrtyp = 'Z'
               gdout = ezqkdef(nix,njy,cgrtyp,ip1,ip2,ip3,0,iunit)
            endif
	         if (grille_u) then
               irecyy = fstinf(iunit, niu, nju, nku, -1, '            ', ip1, ip2, ip3, '  ', '^>  ')
               cgrtyp = 'U'
	            ier = fstprm(irecyy, dateou, deetu, npasu, niu, nju, nku, nbitsu, datypu, lip1, lip2, lip3, typvaru, nomu,etiku&
     &, grref, ig1ref, ig2ref, ig3ref, ig4ref, swa, lng, dltf, ubc, extra1, extra2, extra3)
	            gdout = ezqkdef(niu,nju,cgrtyp,lip1,lip2,lip3,0,iunit)
	            ier = ezgxprm(gdout, nix, njy, cgrtyp, ig1, ig2, ig3, ig4, grref, ig1ref, ig2ref, ig3ref, ig4ref)
            endif
	         li = nix
            lj = njy
            lg1 = ig1
            lg2 = ig2
            lg3 = ig3
            lg4 = ig4
         else
            ireclo = fstinf(iunit, nix, njx, nkx, -1, '            ', ip1, ip2, ip3, '  ', '>>  ')
            irecla = fstinf(iunit, niy, njy, nky, -1, '            ', ip1, ip2, ip3, '  ', '^^  ')
            ier = fstprm(ireclo, dateox, deetx, npasx, nix, njx, nkx, nbitsx, datypx, lip1, lip2, lip3, typvarx, nomx, etikx, grref&
     &, ig1ref, ig2ref, ig3ref, ig4ref, swa, lng, dltf, ubc, extra1, extra2, extra3)
            li = nix
            lj = njy
            gdout = ezqkdef(li,lj,'Z',lip1,lip2,lip3,0,1)
         endif
         allocate(tmplon(li,lj))
         allocate(tmplat(li,lj))
         ier = gdll(gdout, tmplat, tmplon)
      else
         nix = ncoords
         niy = ncoords
         njx = 1
         njy = 1
         cgrtyp = 'Y'
         cgtypxy= 'L'
         nomvarx = '>>  '
         nomvary = '^^  '
         etikety=  'NORDSUD     '
         etiketx = 'ESTOUEST    '
         typvarx = 'C '
         typvary = 'C '
         ip1x = ip1
         ip2x = ip2
         ip3x = ip3
         deetx= 0
         npasx= 0
         dateox = 017901000
         li = ncoords
         lj = 1
         allocate(tmplon(li,lj))
         allocate(tmplat(li,lj))
         call cxgaig('L',ig1la,ig2la,ig3la,ig4la,0.,0.,1.0,1.0)
         call cxgaig('L',ig1lo,ig2lo,ig3lo,ig4lo,0.,0.,1.0,1.0)
         npack = -32
          do i=1,li*lj
             print *, i, dateox,typvarx,li,lj
          enddo
      endif
!
!     INITIALISATION DE DGRWXY POUR UTILISATION DANS  ROUTINE VDAUV
!
      if (cgtypxy.eq.'E') then
         allocate(tmplong(li,lj))
         allocate(tmplatg(li,lj))
         allocate(tmplon(li,lj))
         allocate(tmplat(li,lj))
      endif
!
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
!
      if (it.eq.gr_stations) then
         call fillcoord(tmplat,tmplon)
         gdout = ezgdef(nix,njx,cgrtyp,cgtypxy,ig1la,ig2la,         ig3la,ig4la,tmplon,tmplat)
         nj = 1
         nk = 1
         if (mode.eq.1) then
            call ecritur(tmplon,npack,dateox,deetx,npasx,ncoords,nj,nk,            ip1x,ip2x,ip3x,            typvarx,nomvarx,etike&
     &tx,cgtypxy,ig1lo,ig2lo,ig3lo,ig4lo)
            call ecritur(tmplat,npack,dateox,deetx,npasx,ncoords,nj,nk,            ip1x,ip2x,ip3x,            typvary,nomvary,etike&
     &ty,cgtypxy,ig1la,ig2la,ig3la,ig4la)
            lg1 = ip1x
            lg2 = ip2x
            lg3 = ip3x
            lg4 = 0
         endif
      else
!         ier = gdgaxes(gdout, tmplon, tmplat)
          ier = gdll(gdout, tmplat, tmplon)
      endif
      if (printen)  call imprime(nomvary,tmplat,niy,njy)
      if (printen)  call imprime(nomvarx,tmplon,nix,njx)
!
!     CALCUL LATITUDES LONGITUDES  DU TYPE "Z" OU "Y"
!
      if (cgrtyp.eq.'Z') then
          ier = ezgxprm(gdout,li,lj,cgrtyp,lg1,lg2,lg3,lg4,cgtypxy,ig1ref,ig2ref,ig3ref,ig4ref)
          ier = gdll(gdout,tmplat,tmplon)
      endif
      if (printen) write(6,*)' IMPRIME LAT LON APRES CALL GDLL'
      if (printen)  call imprime(nomvarx,tmplat,niy,njy)
      if (printen)  call imprime(nomvary,tmplon,nix,njx)
!
      return
      end
!
!**S/P GRLALON   CALCUL LATITUDE LONGITUDE DE CHAQUE PT D'UNE GRILLE LATLON
!
      subroutine grlalon(nni,nnj,p1,p2,p3,p4)
      implicit none
!
!AUTEUR   - P. SARRAZIN JANVIER 87 DRPN DORVAL P.Q. CANADA
!
!LANGAGE - RATFOR
!
!OBJET(GRLALON)
!          CALCULER LA LATITUDE ET LA LONGITUDE DE TOUS LES POINTS
!          DE LA GRILLE DE SORTIE LATLON INTERVAL REGULIER TYPE "L"
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!------------------------------------------------------
       external memoir,pgsmabt,grll,cigaxg,cxgaig,messags
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
!
!
      character(len=512) :: defo(990)
      character(len=512) :: listl(990), form
      character(len=512) :: lfn(990)
      common/carac/defo,listl,form,lfn
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
       external ezqkdef, gdll
       integer ezqkdef, gdll
       integer nni,nnj,ier
       real p1,p2,p3,p4,pp1,pp2,pp3,pp4
!
       li=nni
       lj=nnj
!
!   RESERVER MEMOIR POUR LATITUDE ET LONGITUDE ET CALL ECRIRE (IWRK)
!
       allocate(tmplat(nni,nnj))
       allocate(tmplon(nni,nnj))
!
!
!   RECALCUL  LAT,LONG,DELTA LAT,DELTA LONG
!   MEME VALEUR A L'ENTRE COMME A LA SORTIE
!
!
       cgrtyp='L'
       call cxgaig(cgrtyp,lg1,lg2,lg3,lg4,p1,p2,p3,p4)
       call cigaxg(cgrtyp,pp1,pp2,pp3,pp4,lg1,lg2,lg3,lg4)
!
       gdout = ezqkdef(li,lj,cgrtyp,lg1,lg2,lg3,lg4,0)
       ier = gdll(gdout, tmplat, tmplon)
       return
       end
!
!**S/P HEURE INITIALISER TABLE HEURE
!
      subroutine heure(ih1,ih2,ih3,ih4,ih5,ih6,ih7,ih8,ih9,ih10,      ih11,ih12,ih13,ih14,ih15,ih16,ih17,ih18,ih19,ih20,      ih21,&
     &ih22,ih23,ih24,ih25,ih26,ih27,ih28,ih29,       ih30,ih31,ih32,ih33,ih34,ih35,ih36,ih37,ih38,ih39,ih40)
      implicit none
      external pgsmabt,messags
!
!AUTEUR P. SARRAZIN JANVIER 82 DRPN DORVAL P.Q. CANADA
!
!REVISION 
!        P. SARRAZIN JAN 85 POUR AUGMENTER DE 20 A 40 HEURES
!
!LANGAGE RATFOR
!
!OBJET(HEURE)
!         EXTRAIRE LES HEURES DEMANDER PAR L'USAGER ECRIRE DANS
!         LA TABLE HEURE
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN     IH1.....IH40  HEURE DEMANDER PAR L'USAGER (READLX)
!
!MESSAGES 
!         MAUVAISE DIRECTIVE HEURE NHEURE=
!
!MODULES
!         PGSMABT
!
!APPEL   VIA DIRECTIVE
!         HEURE(IH1,IH2,IH3.....................IH40)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   integer heures,nhur,nheure
   common / heures/ nhur,nheure,heures(40)
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
!
      integer ih1,ih2,ih3,ih4,ih5,ih6,ih7,ih8,ih9,ih10,ih11,ih12,ih13
      integer ih14,ih15,ih16,ih17,ih18,ih19,ih20,ih21,ih22,ih23,ih24
      integer ih25,ih26,ih27,ih28,ih29,ih30,ih31,ih32,ih33,ih34,ih35
      integer ih36,ih37,ih38,ih39,ih40
!     
!     
      nheure = min0(40,nheure)
!     
!     
!     
      go to (1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,      23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40)n&
     &heure
 40   heures(40) = ih40
 39   heures(39) = ih39
 38   heures(38) = ih38
 37   heures(37) = ih37
 36   heures(36) = ih36
 35   heures(35) = ih35
 34   heures(34) = ih34
 33   heures(33) = ih33
 32   heures(32) = ih32
 31   heures(31) = ih31
 30   heures(30) = ih30
 29   heures(29) = ih29
 28   heures(28) = ih28
 27   heures(27) = ih27
 26   heures(26) = ih26
 25   heures(25) = ih25
 24   heures(24) = ih24
 23   heures(23) = ih23
 22   heures(22) = ih22
 21   heures(21) = ih21
 20   heures(20) = ih20
 19   heures(19) = ih19
 18   heures(18) = ih18
 17   heures(17) = ih17
 16   heures(16) = ih16
 15   heures(15) = ih15
 14   heures(14) = ih14
 13   heures(13) = ih13
 12   heures(12) = ih12
 11   heures(11) = ih11
 10   heures(10) = ih10
 9    heures(9)  = ih9
 8    heures(8)  = ih8
 7    heures(7)  = ih7
 6    heures(6)  = ih6
 5    heures(5)  = ih5
 4    heures(4)  = ih4
 3    heures(3)  = ih3
 2    heures(2)  = ih2
 1    heures(1)  = ih1
!     
!     
!     #  change de location a cause de readlx
      nhur = nheure
      return
      end
!
!**S/P  IMPRIME CHAMP LUT SUR FICHIER D ENTRE OU DE SORTI
!
      subroutine imprime(cnom,champ,ni,nj)
      implicit none
!
!AUTEUR P. SARRAZIN JUIN 85 DRPN DORVAL P.Q. CANADA
!
!REVISION 4.0.2
!   MODIF. ARGUMENT NOM (ENTIER -> CHARACTER*2)
!   Y. CHARTIER DRPN DORVAL QUEBEC
!
!LANGAGE RATFOR
!
!OBJET(IMPRIME)
!        IMPRIME AVEC LA DIRECTIVE PRINTEN RECORD LUT SUR FICHIER D ENTRE
!        OU IMPRIME AVEC LA DIRECTIVE PRINTSR RECORD QUE L ON VA ECRIRE
!        L USAGER CONTROL LE NOMBRE DE LOCATIONS A IMPRIMER 
!        FENETRE DU CHAMP A IMPRIMER DEFINIT PAR L'USAGER
!        DANS LA DIRECTIVE PRINTEN/PRINTSR MODIFIE LE COMMON
!        LIRES OU ECRIRES PRINTEN=OUI,NIS,NJS,NIF,NJF,NINC,NJNC
!        NIS = POINT DE DEPART DANS LA DIRECTION I (EST-OUEST)
!        NJS = POINT DE DEPART DANS LA DIRECTION J (NORD-SUD)
!        NIF = DERNIER POINT DANS LA DIRECTION I (EST-OUEST)
!        NJF = DERNIER POINT DANS LA DIRECTION J (NORD-SUD) 
!        NINC= INTERVAL DANS LA DIRECTION I
!        NJNC= INTERVAL DANS LA DIRECTION J
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN     NOM   -NOM DU CHAMP QUE L ON VEUT IMPRIMER 2 CARACTERES
!  IN     CHAMP -CONTIENT LE CHAMP QUE L ON VEUT IMPRIMER
!  IN     NI    -DIMENSION DU CHAMP EST=OUEST
!  IN     NJ    -DIMENSION DU CHAMP NORD-SUD
!
!
!MESSAGES 
!
!MODULES
!         PGSMABT
!
!APPEL   VIA DIRECTIVE
!         PRINTEN(OUI,NIS,NJS,NIF,NJF,NINC,NJNC)
!         PRINTSR(OUI,NIS,NJS,NIF,NJF,NINC,NJNC)
!
! - -  - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
!
      integer i,j,ni,nj,njff,niff,niifs,njjfs
      character*4 cnom
      real champ(ni,nj)
      niff=nif
      njff=njf
      if (nif.gt.ni)  niff=ni
      if (njf.gt.nj)  njff=nj
!     
      write(6,600) cnom,nis,njs,niff,njff,ninc,njnc
 600  format(' PRINT CHAMP(LU) NOM=',a2,'  NIS=',i3,'  NJS=',i3,      '  NIFF=',i3,      '  NJFF=',i3,'  NINC=',i3,'  NJNC=',i3)
 620  format(' PRINT CHAMP(ECRIT) NOM=',a2,'  NIS=',i3,'  NJS=',i3,      '  NIFF=',i3,      '  NJFF=',i3,'  NINC=',i3,'  NJNC=',i3)
!     
      do j=njs,njff,njnc
         write(6,630) j
         write(6,610) (champ(i,j),i=nis,niff,ninc)
      enddo
!     
      if (niff.lt.nis)  then
         write(6,*)    ' NIS.lt.NIF DIRECTIVE PRINTEN=OUI,NIS,NJS,NIF,NJF,NINC,NJNC'
      endif
      if (njff.lt.njs)  then
         write(6,*)      ' NJS.lt.NJF DIRECTIVE PRINTEN=OUI,NIS,NJS,NIF,NJF,NINC,NJNC'
      endif
 610  format(1h ,10e13.5)
 630  format('  RANGEE NO ',i3)
      return
      entry imprims(cnom,champ,ni,nj)
!     
      niifs=niif
      njjfs=njjf
      if (niif.gt.ni)  niifs=ni
      if (njjf.gt.nj)  njjfs=nj
!     
      write(6,620) cnom,niis,njjs,niifs,njjfs,niinc,njjnc
!     
      do j=njjs,njjfs,njjnc
         write(6,630) j
         write(6,610) (champ(i,j),i=niis,niifs,niinc)
      enddo
!     
      if (niifs.lt.niis)  then
         write(6,*)      ' NIS.lt.NIF DIRECTIVE PRINTSR=OUI,NIS,NJS,NIF,NJF,NINC,NJNC'
      endif
      if (njjfs.lt.njjs)  then
         write(6,*)      ' NJS.lt.NJF DIRECTIVE PRINTSR=OUI,NIS,NJS,NIF,NJF,NINC,NJNC'
      endif
      return
      end
      subroutine initid
      implicit none
   integer qposition, qitems(16), qnitems
   character*1 qcsepar
   character*16 qcform
   common /idents/  qposition, qitems, qnitems
   common /cidents/ qcsepar,qcform
      integer i
      qposition = 0
      qcsepar = 'T'
      qcform = 'f12.5'
      qnitems = 0
      do i=1,16
         qitems(i) = 0
      enddo
      return
      end
      subroutine initseq
      implicit none
      integer nmaxlist1,nmaxlist2,wait,go
      parameter (nmaxlist1=16,nmaxlist2=16)
      parameter (wait=0,go=1)
      character*2 listnom(nmaxlist1,nmaxlist2)
      integer listniv(nmaxlist1,nmaxlist2)
      integer ntitems,nitems1(nmaxlist1),nitems2(nmaxlist2)
      common /seq/  ntitems,nitems1,nitems2
      common /cseq/ listnom
!     
!     Initialisation des listes utilisees par champ_seq
!     
      integer  i,j
      do i=1,nmaxlist1
         do j=1,nmaxlist2
            listnom(i,j) = '    '
            listniv(i,j) = -1
         enddo
         nitems1(i)=0
         nitems2(i)=0
      enddo
      ntitems = 0
      return
      end
!**   S/P ITROUVE VERIFIER DANS LA LISTE(NOMBRE) SI IVARIA EXISTE
!     
      integer function itrouve(liste,nombre,ivaria)
!AUTEUR P. SARRAZIN RPN DORVAL FEV 81
!
!LANGAGE RATFOR
!
!OBJET(ITROUVE)
!         VERIFIER SI IVARIA EXISTE DANS LISTE SI OUI ITROUVE >0
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN     LISTE =TABLE DE DONNEES
!  IN     NOMBRE=NOMBRE DE DONNEES DANS LISTE
!  IN     IVARIA=ITEM A VERIFIER DANS LISTE
!
! --------------------------------------------------------------------
!
      implicit none
!
!
      integer liste(1),nombre,ivaria,ntr
      itrouve=0
!  #  defense contre index de zero 
      if (nombre.le.0)  return
      do ntr=1,nombre
         if (ivaria.eq.liste(ntr)) itrouve=ntr
      enddo
      return
      end
!
      subroutine lastcol(sortie,valeur,istart,ifini,incre)
      implicit none
!
!AUTEUR P. SARRAZIN DORVAL QUEBEC CANADA (DRPN)
!
!OBJET(LASTCOL)
!         LASTCOL - LA DERNIERE COLONNE DU CHAMP PREND LA VALEUR(VALEUR)
!         LOUPNEG - ELIMINE VALEUR NEGATIVE DANS LE CHAMP SI VALEUR=0.0
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   OUT   SORTIE - RESULTAT DE L'OPERATION SUR CHAMP(NOMBRE)
!   IN    VALEUR - VALEUR DE LA DERNIERE COLONNE DU CHAMP
!                  OU VALEUR MINIMAL DU CHAMP A GARDE
!   IN    ISTART - PREMIERE INDEX DU CHAMP
!   IN    IFINI  - DERNIERE INDICE DU CHAMP
!   IN    INCRE  - NOMBRE DE PTS AUGMENTATION DE L'INDICE ISTART
!
!                  NOMBRE EST (MAXIMUM) LA MOITIE DU NOMBRE DES POINTS
!
!APPEL
!     - VIA GRILLE,MACPCP
!     - CALL LASTCOL(SORTIE,VALEUR,ISTART,IFINI,INCRE)
!     - CALL LOUPNEG(SORTIE,VALEUR,ISTART,IFINI,INCRE)
!
!--------------------------------------------------------------------------
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
!
      real sortie(1),valeur
      integer istart,ifini,incre,i
!     
!         INITIALISE LA DERNIERE COLONNE DU CHAMP 
!
      do i=istart,ifini,incre
         sortie(i)=valeur
      enddo
      return
!
!----------------------------------------------------------------------
!
      entry loupneg(sortie,valeur,istart,ifini,incre)
!
!              ELIMINE VALEUR NEGATIVE
!
      do i=istart,ifini,incre
         sortie(i)=amax1(sortie(i),valeur)
      enddo
!     
      return
      end
      integer function legvar(x)
      legvar=0
      return
      end
      subroutine liraxez(iun, ni, nj, nk, ig1, ig2, ig3, ig4)
      implicit none
      integer iun, ni, nj, nk, ig1, ig2, ig3, ig4
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
      data niz, njz, nkz 	/0, 0, 0/
      data ig1ref, ig2ref, ig3ref, ig4ref /-1, -1, -1, -1/
      data ig1z,   ig2z,   ig3z,   ig4z   /-1, -1, -1, -1/
      data grref          /'Z'/
      integer  fstinf, fstprm, fstluk
      external fstinf, fstprm, fstluk, pgsmabt,messags
      character *12 citiky,citikx
      character *4 cnmvar1,cnmvar2
      character *2 ctpvry,ctpvrx
      character *1 cgtypx,cgtypy
      integer irecy,irecx
      integer idatt,idett,npas,niy,njy,nky,jjp1,jjp2,jjp3
      integer ig1y,ig2y,ig3y,ig4y,nix,njx,nkx
      integer ig1x,ig2x,ig3x,ig4x
      integer numy,numx,ier,iopc
      integer cnbits,cdatyp,extra3,extra2,extra1,cubc,cdltf,clng,cswa
!     
      if (ni.eq.niz.and.nj.eq.njz.and.nk.eq.nkz.and.ig1.eq.ig1z.and.      ig2.eq.ig2z.and.ig3.eq.ig3z.and.ig4.eq.ig4z) then
         return
      endif
      if (associated(axex)) then
	 deallocate(axex)
      endif
      if (associated(axey)) then
	 deallocate(axey)
      endif
      irecx = fstinf(iun, nix, njx, nkx, -1,' ',       ig1, ig2, ig3,'  ','>>  ')
      irecy = fstinf(iun, niy, njy, nky, -1,' ',       ig1 ,ig2, ig3,'  ','^^  ')
      if (irecy .lt. 0 .or.irecx .lt. 0) then
         print *,          'FSTINF ROUTINE LIRAXEZ RECORD ^^ OU >> MANQUANT...'
         call pgsmabt
      endif
!     
      if (nky .gt. 1) then
         write(6,*)'************************************************'
         write(6,*)'         PGSM N''ACCEPTE PAS UN          '
         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 (RECORD ^^ OU >>)'
         write(6,*)'************************************************'
         call pgsmabt
      endif
      ier = fstprm(irecy,idatt,idett,npas,niy,njy,nky, cnbits,cdatyp,      jjp1,jjp2,jjp3,ctpvry,cnmvar1,citiky,cgtypy,      ig1y,i&
     &g2y,ig3y,ig4y, cswa, clng, cdltf, cubc,      extra1, extra2, extra3)
      if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR LIRAXEZ'
!     
!     verifier si grille gaussienne ni doit etre pair
!     
      if (cgtypy.eq.'G'.and.mod(niy,2).ne.0)  call messags(niy)
!     
      ier = fstprm(irecx,idatt,idett,npas,nix,njx,nkx, cnbits,cdatyp,      jjp1,jjp2,jjp3,ctpvrx,cnmvar2,citikx,cgtypx,      ig1x,i&
     &g2x,ig3x,ig4x, cswa, clng, cdltf, cubc,      extra1, extra2, extra3)
      if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR LIRAXEZ'
!     
!     verifier si grille gaussienne ni doit etre pair
!     
      if (cgtypx.eq.'G'.and.mod(nix,2).ne.0)  call messags(nix)
!     
!     
      if (cgtypy.ne.cgtypx .or. ig1y.ne.ig1x .or. ig2y.ne.ig2x .or.       ig3y.ne.ig3x .or. ig4y .ne.ig4x) then
         write(6,*)' VERIFIER RECORD "^^" OU ">>" '
         write(6,*)' GRTYP,IG1,IG2,IG3,IG4 SONT DIFFERENTS'
         call pgsmabt
      endif
!     lecture des axes
      allocate(axex(nix))
      allocate(axey(niy))
      ier = fstluk(axex, irecx, nix, njx, nkx)
      if (ier .lt. 0) then
         print *, 'PROBLEME LECTURE ENREGISTREMENT >> '
         stop
      endif
      ier = fstluk(axey, irecy, niy, njy, nky)
      if (ier .lt. 0) then
         print *, 'PROBLEME LECTURE ENREGISTREMENT ^^ '
         stop
      endif
      niz = nix
      njz = njy
      nkz = nkx
      ig1z = ig1
      ig2z = ig2
      ig3z = ig3
      ig4z = ig4
      grref = cgtypx
      ig1ref = ig1x
      ig2ref = ig2x
      ig3ref = ig3x
      ig4ref = ig4x
      return
      end
!
!**S/P   LIREN   LIRE UN CHAMP DANS ACCUMULATEUR
!
      subroutine liren(nom, type, idat, niv, ihr, ip3, etiqet)
      implicit none
      external fstinf,pgsmlir,memoir,fstprm,pgsmabt,imprime
      external fstopc,messags,fstcvt
      integer fstinf,pgsmlir,fstprm,fstopc,fstcvt
!
!AUTEUR P. SARRAZIN  AOUT 82 DRPN DORVAL P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(LIREN)
!         LIRE UN CHAMP SUR FICHIER 1 OU 2 ET SAUVE DANS UN ACCUMULATEUR
!         POUR ETRE UTILISER PAR LES DIRECTIVES PLUSE-PLUSS
!         MOINSE-MOINSS-PFOIS-MOYENE-RACINE-MODUL2E-MODUL2S
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN   NOM     NOM DU CHAMP LCAR(GZ),"TT"......
!   IN   TYPE    TYPE DE CHAMP "P"=PREVISION  "A" ANALYSE
!   IN   NIV     NIVEAU DU CHAMP
!   IN   IHR     HEURE DU CHAMP
!   IN   IP3     LIBRE(USAGER) COMPTEUR POUR MOYENE UTILISER PAR ECRITS
!   IN   ETIQET  ETIQUETTE 10 CARACTERES
!
!IMPLICITES
!MESSAGES
!         RECORD N EXISTE PAS SUR FICHIER (FSTINF DANS LIREN)
!         RECORD N EXISTE PAS (PGSMLIR DANS ROUTINE LIREN)
!
!MODULES  FSTINF,PGSMABT,FSTPRM,MEMOIR,PGSMLIR
!
!APPEL     VIA DIRECTIVE
!         LIREE(NOM, TYPE, IDAT, NIV, IHR, IP3, ETIQUET)
!         LIRES(NOM, TYPE, IDAT, NIV, IHR, IP3, ETIQUET)
!
! -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
      integer blancs
      data blancs /4H    /
      integer ip1style, dateform
      common /styles/ ip1style, dateform
!
!
      character *12 cetiqet
      character *4 cnomvar
      character *2 ctypvar
      character*1 cigtyp
      integer etiqet(3),idat,ihr,iip3,ip3,irec1,iunit,niv(2),nom,num1,type
      integer cnbits,cdatyp,extra1,extra2,extra3,cubc,cdltf,clng,cswa
      integer iopc
      integer argdims, letiket(3)
      external argdims
      integer lniv
      real p
      character*8 string
!
      iunit=1
      iip3=ip3
      if (ip3.eq.4095)iip3=-1
!
!     MODIFICATION DE HOLLERITH A CARACTERE
!
 100  cnomvar = '    '
      ctypvar = '  '
      cetiqet = '            '
      cigtyp  = ' '
      letiket(1) = etiqet(1)
      letiket(2) = blancs
      letiket(3) = blancs
      if (argdims(7).gt.1) then
         letiket(2) = etiqet(2)
      endif
      if (argdims(7).gt.2) then
         letiket(3) = etiqet(3)
      endif
      lniv = niv(1)
      if (argdims(4) > 1) then
         p = transfer(niv(1), p)
         call convip_plus(lniv, p, -1*niv(2)-1000, ip1style, string, .false.)
      endif
      ier = fstcvt(    nom,   type, letiket,    -1,      cnomvar,ctypvar,cetiqet,cigtyp,     .true.)
      irec1=fstinf(iunit,nni,nnj,nnk,idat,cetiqet,lniv,ihr,iip3,      ctypvar,cnomvar)
      if (irec1 .lt. 0)   then
         write(6,*)         'RECORD N EXISTE PAS SUR FICHIER (FSTINF LIREE-LIRES)'
         call pgsmabt
      endif
!
!      if (nnk.gt.1)   then
!         write(6,*)'*************************************************'
!         write(6,*)'         PGSM N ACCEPTE PAS UN          '
!         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (LIREE-LIRES)'
!         write(6,*)'*************************************************'
!         call pgsmabt
!      endif
!
!
!  #  clef pour directive pluse,moinse,ecrits....
      ichck=1
!
!
      ier = fstprm(irec1,idatt,ideet,npas,nni,nnj,nnk, cnbits,cdatyp,      jpp1,jpp2,jpp3,ctypvar,cnomvar,cetiqet,cigtyp,igg1,igg2,&
     &igg3,      igg4,cswa, clng, cdltf, cubc, extra1, extra2, extra3)
      if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR LIREN'
      cnumv = cnomvar
      ctypv = ctypvar
      cetik = cetiqet
      cigty = cigtyp
!
!       VERIFIER SI GRILLE GAUSSIENNE NI DOIT ETRE PAIR
!
      if (cigtyp.eq.'G'.and.mod(nni,2).ne.0)  call messags(nni)
!
!
!    ALLOCATION DE LA MEMOIRE
!
      allocate(tmpif0(nni,nnj))
!
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
      num1 =pgsmlir(tmpif0,iunit,nni,nnj,nnk,idat,cetiqet,jpp1,jpp2,      jpp3,ctypvar,cnomvar,cigtyp)
!
      if (num1 .lt. 0) then
         write(6,*)'RECORD N EXISTE PAS (LIRE LIREE-LIRES)'
         call pgsmabt
      endif
      if (printen)  call imprime(cnomvar,tmpif0,nni,nnj)
!
!     SI COMTEUR .NE. 4095  ICNT=1
!
      icnt = 1
      if (iunit.eq.1.and.ip3.eq.   4095)  icnt = jpp3
      if (iunit.eq.2.and.ip3.eq.   4095)  icnt = jpp3
!
!
      return
!
      entry lirsr(nom, type, idat, niv, ihr, ip3, etiqet)
!
      iunit = 2
      iip3=ip3
      if (ip3.eq.4095) iip3=-1
      go to 100
      end
!
      subroutine lopascm(sortent,entre,fact,nombre)
      implicit none
!
!AUTEUR P. SARRAZIN DORVAL QUEBEC CANADA (DRPN)
!
!OBJET(LOPASCM)
!         LOPASCM - ADDITIONNE,SOUSTRAIT,MULTIPLI DEUX CHAMPS OU SOMME
!                   DEUX CHAMPS DONT CHAQUE PT EST AU CARRE 
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
! IN-OUT  SORTENT - RESULTAT DE L'OPERATION ET ENTRE DU CHAMP
!   IN    ENTRE   - DEUXIEME CHAMP POUR OPERATION 
!   IN    FACT    - OPERATEUR 1=ADDITIONNE , -1=SOUSTRAIT 3=MULTIPLI
!                             2=ADDITIONNE CHAQUE POINT AU CARRE
!   IN    NOMBRE  - NOMBRE DE POINTS DANS LES CHAMPS SORTENT/ENTRE
!
!APPEL
!     - VIA PLMNMOD 
!      CALL LOPASCM(SORTENT,ENTRE,FACT,NOMBRE)
!
!----------------------------------------------------------------------
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
!
!
      integer nombre,i,fact
      real sortent(nombre),entre(nombre)
!     
!       SOUSTRAIT SI FACT=-1  ADDITIONNE SI FACT=1
!
      if (abs(fact).eq.1) then
         do i=1,nombre
            sortent(i)=sortent(i) + entre(i)*fact
         enddo
!     
!     ADDITIONNE CHAQUE POINT DES DEUX CHAMPS AU CARRE
!     
      else if (fact.eq.2)  then
         do i=1,nombre
            sortent(i)=sortent(i)**2 + entre(i)**2
         enddo
!     
!     MULTIPLIT CHAQUE POINT DES DEUX CHAMPS
!     
      else if (fact.eq.3) then
         do i=1,nombre
            sortent(i)=sortent(i)*entre(i)
         enddo
      else  if (fact.eq.4) then
         do i=1,nombre
            if (entre(i).eq.0.0) then
               print *, 'LOPASCM - Un des elements du tableau a une valeur de 0.0'
               print *, 'Division impossible - sortie forcee'
               call qqexit(13)
            endif
         enddo
         do i=1,nombre
            sortent(i)=sortent(i)/entre(i)
         enddo
      else
         write(6,*)' ERREUR DANS ROUTINE LOPASCM (FACT..?)'
      endif
!     
      return
      end
!
      subroutine loupmir(sortie,entre,nombre)
      implicit none
!
!AUTEUR P. SARRAZIN DORVAL QUEBEC CANADA (DRPN)
!
!OBJET(LOUPMIR)
!         LOUPMIR - TRANSFER LA MOITIE DU CHAMP DEFINI DANS L'AUTRE
!                   MOITIE RENVERSER COMME UN MIROIR
!
!         LOUPTRA - TRANSFER UN CHAMP DANS UN AUTRE
!
!         LOUPIN1 - INITIALISE UN CHAMP AVEC LA VALEUR 1
!
!         LOUPSOU - SOUSTRAIRE DEUX CHAMPS
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   OUT   SORTIE - RESULTAT DE L'OPERATION SUR CHAMP(NOMBRE)
!   IN    ENTRE  - CHAMP QUI PEUT SERVIR AU CALCUL
!   IN    NOMBRE - NOMBRE DE POINTS DANS LE CHAMP EXCEPTE LOUPTRA OU
!                  NOMBRE EST (MAXIMUM) LA MOITIE DU NOMBRE DES POINTS
!
!APPEL
!     - VIA GRILLE,COUPZM,LIREN,EPAISUR 
!     - CALL LOUPMIR(SORTIE,ENTRE,NOMBRE)
!     - CALL LOUPTRA(SORTIE,ENTRE,NOMBRE)
!     - CALL LOUPIN1(SORTIE,ENTRE,NOMBRE)
!     - CALL LOUPSOU(SORTIE,ENTRE,NOMBRE)
!
!--------------------------------------------------------------------------
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
!
!
      real sortie(1),entre(1)
      integer nombre,nombpl1,i
!     
!            TRANSFER LA MOITIE DU CHAMP DANS L'AUTRE
!
!
      nombpl1=nombre + 1
!
      do i=1,nombre
         sortie(i + nombre) = entre(nombpl1 - i)
      enddo
!     
      return
!     
!------------------------------------------------------------
!
      entry louptra(sortie,entre,nombre)
!
!
!              TRANSFER CHAMP 
!
      do i=1,nombre
         sortie(i)=entre(i)
      enddo
!     
      return
!-------------------------------------------------------------
!
       entry loupin1(sortie,entre,nombre)
!
!
!           INITIALISE LE CHAMP SORTIE AVEC 1.0
!
       do i=1,nombre
          sortie(i)=1.0
       enddo
!     
       return
!-----------------------------------------------------------
!
       entry loupsou(sortie,entre,nombre)
!
!
!           SOUSTRAIT DEUX CHAMPS
!
       do i=1,nombre
          sortie(i)=entre(i) - sortie(i)
       enddo
!     
       return
       end
!
!**S/P    CHAQUE PT D'UN CHAMP LU EST MI AU CARRE DANS ACCUMULATEUR
!
      subroutine lrsmde(nom, type, idat, niv, ihr, ip3, etiqet)
      implicit none
!
!AUTEUR   P. SARRAZIN  DORVAL QUEBEC CANADA (DRPN)
!
!LANGAGE RATFOR
!
!OBJET(LRSMDE)
!         LIRE UN CHAMP SUR FICHIER D'ENTRE OU DE SORTI ET SAUVE DANS
!         L'ACCUMULATEUR CHAQUE PT AU CARRE ET LES DIRECTIVES SUIVANTES
!         LIRMODE OU LIRMODS AJOUTERONT CHAQUE CHAMP(PT AU CARRE) A
!         L'ACCUMULATEUR LA SOMME DES CHAMPS EST GARDER DANS L'ACCUMULATEUR
!         ET PEUT ETRE SAUVE PAR LA DIRECTIVE ECRITS.
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN   NOM     NOM DU CHAMP "GZ","TT"...LCAR(GZ)
!   IN   TYPE    TYPE DE CHAMP  "P"=PREVISION   "A"=ANALYSE
!   IN   NIV     NIVEAU DU CHAMP 500MB....
!   IN   IHR     HEURE DU CHAMP  (IP2)
!   IN   IP3     LIBRE A L'USAGER ET COMPTEUR POUR MOYENNE UTILISER PAR ECRITS
!   IN   ETIQET  ETIQETTE 10 CARACTERES
!
!
!MESSAGES
!         RECORD N'EXISTE PAS SUR FICHIER (FSTINF DANS LRSMDE-LRSMDS)
!         RECORD N'EXISTE PAS (LIRE DANS ROUTINE LRSMDE-LRSMDS)
!
!MODULES
!        FSTINF,PGSMABT,FSTPRM,MEMOIR,PGSMLIR
!
!APPEL   VIA DIRECTIVE
!        LIRMDE(NOM, TYPE, IDAT, NIV, IHR, IP3, ETIQET)
!        LIRMDS(NOM, TYPE, IDAT, NIV, IHR, IP3, ETIQET)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer blancs
      data blancs /4H    /
      integer ip1style, dateform
      common /styles/ ip1style, dateform
!
      external fstinf,pgsmlir,memoir,fstprm,pgsmabt,imprime, fstopc,messags,fstcvt
      integer fstinf,pgsmlir,fstprm,fstopc,fstcvt
!
      character *12 cetiket
      character *4 cnomvar
      character *2 ctypvar
      character *1 cigtyp
      integer etiqet(3),idat,ihr,iip3,ip3,irec1,iunit,niv(2),nom,num1,type
      integer inomb,i,j,iopc,lniv
      integer cnbits,cdatyp,cswa, clng,cdltf,cubc,extra1,extra2,extra3
      integer argdims, letiket(3)
      external argdims
      real p
      character*8 string
!
      iunit=1
      iip3=ip3
      if (ip3.eq.4095) iip3=-1
!
!     MODIFICATION DE HOLLERITH A CARACTERE
!
      cnomvar = '    '
      ctypvar = '  '
      cetiket = '            '
      cigtyp  = ' '
      letiket(1) = etiqet(1)
      letiket(2) = blancs
      letiket(3) = blancs
      if (argdims(7).gt.1) then
         letiket(2) = etiqet(2)
      endif
      if (argdims(7).gt.2) then
         letiket(3) = etiqet(3)
      endif
      if (argdims(4) > 1) then
         lniv = niv(1)
         p = transfer(niv(1), p)
         call convip_plus(lniv, p, -1*niv(2)-1000, ip1style, string, .false.)
      endif
 100  ier = fstcvt(    nom,   type, letiket ,    -1, cnomvar,ctypvar, cetiket, cigtyp,     .true.)
      irec1=fstinf(iunit,nni,nnj,nnk,idat,cetiket,lniv,ihr,iip3, ctypvar,cnomvar)
      if (irec1 .lt. 0)   then
         write(6,*)' RECORD N EXISTE PAS (FSTINF LIRMDE-LIRMDS)'
         call pgsmabt
      endif
!      if (nnk.gt.1)   then
!         write(6,*)'***************************************************'
!         write(6,*)'         PGSM N ACCEPTE PAS UN          '
!         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (LIRMDE-LIRMDS)'
!         write(6,*)'***************************************************'
!         call pgsmabt
!      endif
!
      ichck=1
!
!
      ier = fstprm( irec1,idatt,ideet,npas,nni,nnj,nnk,cnbits,cdatyp,      jpp1,jpp2,jpp3,ctypvar,cnomvar,cetiket,cigtyp,igg1,igg2,&
     &      igg3,igg4,cswa, clng, cdltf, cubc, extra1, extra2, extra3)
      if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR LIRMDE'
!
!     MODIFICATION DE CARACTERE A HOLLERITH
!
      cnumv = cnomvar
      ctypv = ctypvar
      cetik = cetiket
      cigty = cigtyp
!
!     VERIFIER SI GRILLE GAUSSIENNE NI DOIT ETRE PAIR
!
      if (cigtyp.eq.'G'.and.mod(nni,2).ne.0)  call messags(nni)
!
!
      inomb=nni*nnj*nnk
!
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
!
!
!      if (if9.eq.1)  then
!         num1=pgsmlir(tmpif0,iunit,nni,nnj,nnk,idat,cetiket,jpp1,
!     $        jpp2,jpp3,ctypvar,cnomvar,cigtyp)
!         if (num1 .lt.0) then
!            write(6,*) 'RECORD N EXISTE PAS LIRE (LIRMDE-LIRMDS)'
!            call pgsmabt
!         endif
!
!         if (printen)  call imprime(cnomvar,tmpift,nni,nnj)
!         do i=1,inomb
!            tmpif0(i)=(tmpift(i))**2 + tmpif0(i)
!         enddo
!
!         if (printsr)  then
!            call imprime(cnomvar,tmpif0,nni,nnj)
!         endif
!         icnt = icnt + 1
!      endif
!
!      if (if9.eq.0) then
!
!     ALLOCATION DE LA MEMOIRE
!
      allocate(tmpif0(nni,nnj))
!
!     SI COMPTEUR .NE. 4095 ICNT=1
!
         icnt=1
         if (iunit.eq.1.and.ip3.eq.4095)  icnt=jpp3
!
         if (iunit.eq.2.and.ip3.eq.4095)  icnt=jpp3
!
         num1 = pgsmlir(tmpif0,iunit,nni,nnj,nnk,idat,cetiket,jpp1,          jpp2,jpp3,ctypvar,cnomvar,cigtyp)
         if (num1 .lt. 0)  then
            write(6,*) 'RECORD N EXISTE PAS LIRE (LIRSMDE-LIRSMDS)'
            call pgsmabt
         endif
         if (printen)  call imprime(cnomvar,tmpif0,nni,nnj)
!
!         if9=1
         do j=1,nnj
            do i=1,nni
               tmpif0(i,j)=tmpif0(i,j)*tmpif0(i,j)
            enddo
         enddo
!      endif
!
!
      return
!
      entry lrsmds(nom, type, idat, niv, ihr, ip3, etiqet)
!
      iunit=2
      iip3=ip3
      if (ip3.eq.4095) iip3=-1
!
      go to 100
!
      end
!
!**s/p macpcp interpole ajustement convectif ou precipitation
!
      subroutine macpcp(cnom,npar,itime)
      implicit none
!
!auteur  p.sarrazin fevrier 82  drpn dorval p.q. canada
!
!revision 4.0.2
!     conversion variables type hollerith -> type caractere
!   y. chartier -aout 90- drpn dorval quebec
!
!revision 5.2
!   support des grilles sources de type Z
!
!langage ratfor
!
!objet(macpcp)
!          extraire la difference entre deux champs dont les
!          heures  sont differents
!          avec routine fstinf on extrait le record necessaire pour
!          routine fstprm qui identifit les parametres utilises
!          par routine lire
!          on reserve la memoire pour les deux champs de travail
!          routine ecritur identifit la sorte de fichier utilisee
!          pour ecrire
!
!librairies
!         -source  armnsrc,drpn
!         -objet   pgsmlib,id=armnpjs.
!
!arguments
!  in    nom    nom du champ
!  in    npar    nombre de locations utilisees dans itime
!  in    itime   table contenant 2 heures ou niveaux
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!messages
!          'mauvais appel a champ il devrait y avoir 3 arguments'
!           record n'existe pas sur fichier d'entre (macpcp)
!
!modules
      external ecritur,fstinf,pgsmlir,memoir,fstprm,symetri
      external loupneg,loupsou,fstopc,pgsmabt,imprime,messags,fstcvt
      external liraxez
      integer fstprm,fstinf,pgsmlir,fstopc,fstcvt
!
!appel     via champ
!         call macpcp(nom,npar,itime)
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer dat,deet
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
!
      character *12 cetike,cetiket
      character *4 cnomvar, cnom
      character *1 cigtyp
      character *2 ctypvar
      real valeur
      real fbidon
      real, dimension(:), allocatable :: lclif1, lclif2
      integer i, itime(2),ig1,ig2,ig3,ig4,ip1,irec1,irec2
      integer jp1,jp2,jp3,lilj,ni,nj,nk,nn,npar,num1,num2,iopc
      integer cdatyp,cnbits
      integer cswa, clng, cdltf, cubc, extra1, extra2, extra3
      integer ezqkdef, ezdefset, ezsint, chkenrpos, datev
      integer iunit
      logical symetri,sym
      iunit = 1
      nk = 1
!
      if (npar.ne.2) then
         if (message) then
            write(6,*)'MAUVAIS APPEL DOIT AVOIR 3 ARGUMENTS  (MACPCP)'
         endif
         return
      endif
!
!     identifier le numero de chaque record avec fstinf
!
!     # doit etre egal a zero dans fichier d'ENTRE
      ip1=0
      call chk_userdate(datev)
!
!     modification de hollerith a caractere
!
      cnomvar = cnom
      if (etikent(1) .ne. -1) then
         write(cetiket,'(3A4)') (etikent(i), i=1,nwetike)
      else
         cetiket = '            '
      endif
      if (typeent .ne. -1) then
         write(ctypvar, '(A2)') typeent
      else
         ctypvar = '  '
      endif
      irec1=fstinf(1,ni,nj,nk,datev,cetiket,ip1,      itime(1),ip3ent,ctypvar,cnomvar)
      irec2=fstinf(1,ni,nj,nk,datev,cetiket,ip1,      itime(2),ip3ent,ctypvar,cnomvar)
!     #  record n'EXISTE PAS
      if (irec2 .lt. 0 .or.irec1 .lt. 0) then
         write(6,*)'RECORD N EXISTE PAS SUR FICHIER D ENTRE (MACPCP)'
         write(6,*)' VERIFIER IP2-IP3ENT  IP1=0 SUR FICHIER D ENTRE'
         return
      endif
      if (nk.gt.1) then
         write(6,*)'************************************************'
         write(6,*)'         PGSM N ACCEPTE PAS UN          '
         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (MACPCP)'
         write(6,*)'************************************************'
         call pgsmabt
      endif
!
!
!     identifier parametres pour champ 1
!
      ier = fstprm( irec1, dat,deet,npas,ni, nj, nk, cnbits,cdatyp,      jp1,jp2, jp3,ctypvar,cnomvar,cetike,cigtyp,       ig1,ig2,&
     &ig3,ig4,      cswa, clng, cdltf, cubc, extra1, extra2, extra3)
      if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR MACPCP'
!
!     verifier si grille gaussienne ni doit etre pair
!
      if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  call messags(ni)
!
!
!     lire champ no 1
!
      allocate(lclif1(ni*nj))
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
!      if (userdate .eq. oui .or. userdate .eq. non) date=date2
      num1 = pgsmlir(lclif1,1,ni,nj,nk,datev,cetiket,jp1,jp2,jp3,      ctypvar,cnomvar,cigtyp)
      if (printen)  call imprime(cnomvar,lclif1,ni,nj)
!
!     identifier parametres pour champ 2
!
      ier = fstprm(irec2, dat,deet,npas,ni, nj, nk, cnbits,cdatyp,      jp1,jp2, jp3,ctypvar,cnomvar,cetike,cigtyp,      ig1,ig2,ig&
     &3,ig4,      cswa, clng, cdltf, cubc, extra1, extra2, extra3)
      if (ier .lt. 0) write(6,*)' IER = FSTPRM NEGATIF VOIR MACPCP'
!
!     verifier si grille gaussienne ni doit etre pair
!
      if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  call messags(ni)
!
!     lire champ 2
!
      allocate(lclif2(max0(li*lj,ni*nj)))
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
!      if (date .eq. 0 .or. date .eq. 1) date=date2
      num2 = pgsmlir(lclif2,1,ni,nj,nk,datev,cetiket,jp1,jp2,jp3,      ctypvar,cnomvar,cigtyp)
      if (printen)  call imprime(cnomvar,lclif2,ni,nj)
!
!     difference entre les deux champs
!
      nn = ni*nj
      call loupsou(lclif1,lclif2,nn)
!
!
!     interpolation horizontale
!
      if (cgrtyp.eq.'*') then
         ier = chkenrpos(1,2,ig1,ig2,ig3)
      else
!     #  variable symetrique oui=.true.
         if (cigtyp == 'A' .or. cigtyp == 'B' .or. cigtyp == 'G') then
            if (ig1 /= 0) sym = symetri(cnomvar)
         endif
         gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
         ier = ezdefset(gdout, gdin)
         ier = ezsint(lclif2, lclif1)
      endif
!
!     #   jp1 - contient heure du premier champ
!     #   jp2 - contient heure du deuxieme champ
      jp1 = itime(1)
      jp2 = itime(2)
      jp3 = 0
!
!     deet et npas contiennent les dernieres valeurs lues dans le dernier record
!
!     eliminer toutes les valeurs du champ negative precip
!     et acumulateur d'ajustement ne peuvent etre negatif
!
      lilj=li*lj
      valeur=0.0
      call loupneg(lclif2,valeur,1,lilj,1)
!
!     ecrire sur fichier standard,ms,sequentiel
!
      if (cgrtyp.eq.'*') then
         call ecritur(lclif1,npack,dat,deet,npas,ni,nj,nk,jp1,jp2,jp3,         ctypvar,cnomvar,cetike,cigtyp,ig1,ig2,ig3,ig4)
      else
         call ecritur(lclif2,npack,dat,deet,npas,li,lj,nk,jp1,jp2,jp3,         ctypvar,cnomvar,cetike,cgrtyp,lg1,lg2,lg3,lg4)
      endif
!
!     remetre espace des champs de travail
!
      deallocate(lclif2)
      deallocate(lclif1)
!
      return
      end
!
!**S/P MESSAGS  IMPRIME MESSAGE SUR UNE PAGE COMPLETE
      subroutine messags(ni)
      implicit none
!
!LANGAGE RATFOR
!
!OBJET(MESSAGS)
!          IMPRIME UN MESSAGE A CHAQUE FOIS UNE GRILLE GAUSSIENNE
!          N'A PAS UN NOMBRE DE LONGITUDE PAIRES
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN     NI  - NOMBRE DE LONGITUDES DANS LE RECORD DU FICHIER D'ENTRE
!     
      integer ni
      write(6,600)
 600  format(1h1)
      write(6,*)'**************************************************'
      write(6,*)'*                                                '
      write(6,*)'*                    ATTENTION                   '
      write(6,*)'*                                                '
      write(6,*)'*          NOMBRE DE LONGITUDES                  '
      write(6,*)'*          DOIT-ETRE PAIR                        '
      write(6,*)'*          POUR UNE GRILLE GAUSSIENNE   # LONG=',ni
      write(6,*)'*                                                '
      write(6,*)'*          GARBAGE IN   GARBAGE OUT     OUCH  ?? '
      write(6,*)'*                                                '
      write(6,*)'**************************************************'
      write(6,600)
!     
      return
      end
!
!**S/P CMETSYM   MISE A JOUR DES TABLES DE SYMETRIE
!
      subroutine cmetsym(cnom,sym)
      implicit none
      external pgsmabt,messags
!
!AUTEUR  P.SARRAZIN  FEVRIER 82  DRPN DORVAL  P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(CMETSYM)
!          DEFINIR NOM D UN CHAMP AVEC .TRUE. OU .FALSE.
!          VRAI=SYMETRIQUE  FAUX ANTISYMETRIQUE
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN   NOM    NOM DU CHAMP
!  IN   SYM    VALEUR  .TRUE. / .FALSE. 
!
!IMPLICITES
!MESSAGES 
!         PLUS DE PLACE DANS LES TABLES DE SYMETRIE
!         MAUVAISE DIRECTIVE METSYM DOIT AVOIR 2 ARGUMENTS
!
!
!APPEL    - VIA MAIN PGSM
!         CALL CMETSYM(CNOM,.TRUE./.FALSE.)
!
!
!MODULES  PGSMABT
!
!---------------------------------------------------------------------
!
      integer nnoms,maxnoms,nsym,nsymm
      logical ssym(256)
      common /symnom/ nnoms,maxnoms,nsym,nsymm,ssym
      character*4 noms(256)
      common /csymnom/ noms
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
!
      logical sym
      character*4 cnom
!     
!     sauve nsym dans nsymm a cause readlx qui remet a zero      nsymm = nsym 
!
      if (nnoms.lt.maxnoms) then
         nnoms = nnoms + 1
         ssym(nnoms) = sym
         noms(nnoms) = cnom
      else
         if (message) then
            write(6,*)            'PLUS DE PLACE DANS LES TABLES DE SYMETRIE(METSYM)'
         endif
      endif
      return
      end
!**   S/P METSYM   MISE A JOUR DES TABLES DE SYMETRIE
!     
      subroutine metsym(nom,sym)
      implicit none
      external cmetsym
!
!AUTEUR  P.SARRAZIN  FEVRIER 82  DRPN DORVAL  P.Q. CANADA
!
!REVISION 4.0.2
!   SEPARATION DE METSYM EN 2 PARTIES, L'UNE AVEC LE NOM STORE DANS UN ENTIER,
!      L'AUTRE AVEC LE NOM STORE DANS UNE CHAINE DE CARACTERES
!   Y. CHARTIER DRPN DORVAL QUEBEC
!LANGAGE RATFOR
!
!OBJET(METSYM)
!          INTERFACE A LA ROUTINE "CMETSYM"
!
      integer nom
      logical sym
      character*2 cnom
      write(cnom, '(A4)') nom
      write(6, *) 'METSYM: NOM - ',nom,'CNOM - ', cnom
      call cmetsym(cnom, sym)
      return
      end
!
!**   S/P EXTRAIRE L'EXPONENTIEL DE CHAQUE POINT D UN CHAMP DANS ACCUMULATEUR
!     ET MULTIPLIER LE RESULTAT PAR LE FACT
      subroutine operat(fact,ecart,divi)
      implicit none
      external pgsmabt,messags
!
!AUTEUR P.SARRAZIN JUIN 83 DRPN DORVAL P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(OPERAT)
!         EXTRAIRE L'EXPONENTIEL OU LOGARITHME D UN CHAMP DANS
!         L'ACCUMULATEUR DEJA LUT PAR LA DIRECTIVE LIREE OU LIRES
!         ET MULTIPLIER LE RESULTAT PAR FACT
!         CALCUL LA MOYENNE DE CHAQUE POINT DU CHAMP DANS ACCUMULATEUR
!         PFOIS ADDITIONNER OU SOUSTRAIRE UNE CONSTANTE DU CHAMP DANS
!         L'ACCUMULATEUR ET MULTIPLIER OU DIVISER LE RESULTAT
!         EXTRAIRE LA RACINE CARRE DE CHAQUE POINT DU CHAMP DANS ACCUMULATEUR
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN    FACT    -FACTEUR DE CONVERSION
!  IN    ECART    - UTILISER PAR PFOIS CONSTANTE AJOUTER A CHAQUE POINT
!  IN    DIVI     - UTILISER PAR PFOIS DIVISE CHAQUE POINT DU CHAMP
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!MESSAGES
!        'DIRECTIVE LIREE OU LIRES DOIT ETRE APPELE '
!        'VERIFIER EXPON-PFOIS-MOYENE-RACINE-ALOGN'
!
!MODULES
!         PGSMABT
!
!APPEL    VIA DIRECTIVE
!          PFOIS(FACT,ECART,DIVI) - EXPON(FACT)- MOYENE(FACT) - ALOGN(FACT)
!          RACINE(FACT)
!
!---------------------------------------------------------------
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
!
! -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      integer i,j, it
      real divi,ec,ecart,fact
!
!     VERIFIER SI DIRECTIVE LIREE OU LIRES A ETE APPELE
!
!     pfois(+ecart*fact/divi)
      it=1
!
!   A CAUSE DE LA DOCUMENTATION ECART DEVIENT FACT POUR PFOIS
!
      ec=fact
      fact=ecart
      ecart=ec
!     erreur faut appeler liree ou lires
 1000 if (ichck.eq.0)   then
         write(6,*)' LIREE LIRES DOIT ETRE APPELE AVANT '
         write(6,*)' EXPON-PFOIS-MOYENE-RACINE-ALOGN'
         call pgsmabt
      endif
!
!     IT=2 EXTRAIT L'EXPONENTIEL DE CHAQUE POINT ET MULTIPLIT PAR FACT
!     IT=0 EXTRAIT LE LOGARITHME DE CHAQUE POINT ET MULTIPLIT PAR FACT
!     IT=1 PFOIS AJOUTE ECART  MULTIPLIT PAR FACT ET DIVISE PAR DIVI
!
!     IT=3 MOYENE CAL CUL LA MOYENNE DE CHAQUE POINT (ICNT)
!     IT=4 EXTRAIRE LA RACINE CARRE DE CHAQUE PT DU CHAMP
!
      if (it.eq.1) then
         if (message) write(6,*)' PFOIS(ECART,FACTEUR,DIVISEUR)'
      endif
!     $(  # exponentiel
      if (it.eq.2) then
         if (message) write(6,*)' EXPON(FACTEUR)'
      endif
!     $(  # moyenne
      if (it.eq.3) then
         if (message) write(6,*)' MOYENE(FACTEUR)'
      endif
!     $(  # racine
      if (it.eq.4) then
         if (message) write(6,*)' RACINE(FACTEUR)'
      endif
!     $(  # logarithme
      if (it.eq.0)  then
         if (message) write(6,*)' ALOGN(FACTEUR)'
      endif
!
      if (it.eq.0) then
         do j=1,nnj
            do i=1,nni
               tmpif0(i,j)= (alog(tmpif0(i,j))*fact)
            enddo
         enddo
      endif
       if (it.eq.1) then
         do j=1,nnj
            do i=1,nni
             tmpif0(i,j)= (tmpif0(i,j)+ecart)*fact/divi
            enddo
         enddo
       endif
       if (it.eq.2) then
         do j=1,nnj
            do i=1,nni
             tmpif0(i,j)= (exp(tmpif0(i,j))*fact)
            enddo
         enddo
       endif
       if (it.eq.3) then
         do j=1,nnj
            do i=1,nni
             tmpif0(i,j)=(tmpif0(i,j)/icnt)*fact
            enddo
         enddo
       endif
       if (it.eq.4) then
         do j=1,nnj
            do i=1,nni
             tmpif0(i,j)= (sqrt(tmpif0(i,j))*fact)
            enddo
         enddo
       endif
       if (it.eq.5) then
         do j=1,nnj
            do i=1,nni
             tmpif0(i,j)= (abs(tmpif0(i,j))*fact)
            enddo
         enddo
       endif
       if (it.eq.6) then
         do j=1,nnj
            do i=1,nni
             tmpif0(i,j)= (tmpif0(i,j)*(tmpif0(i,j))*fact)
            enddo
         enddo
       endif
       return
!
       entry alogn(fact)
!
!    IT=0 ON CALCUL LE LOGARITE DE CHAQUE PT
!
       it=0
       go to 1000
!
!
!     EXPONENTIEL DE CHAQUE POINT DU CHAMP
!
       entry expon(fact)
       it=2
       go to 1000
!
!
!     PRENDRE LA MOYENNE DE CHAQUE PTS DES CHAMPS ACCUMULES
!     DANS L'ACCUMULATEUR ET MULTIPLIER PAR FACT
!
       entry moyene(fact)
       it=3
!
!     VERIFIER SI COMPTEUR EST PLUS GRAND QUE 1
!
!     # erreur
       if (icnt.le.1) then
          if (message) then
             write(6,*)'ON DIVISE PAR UNE VALEUR < OU = A 1 ICNT=',icnt
             if (message) then
                write(6,*)           '   ******  ATTENTION A LA DIRECTIVE MOYENE  ******'
             endif
          endif
       endif
!
       go to 1000
!
!     PRENDRE LA RACINE CARRE DE CHAQUE POINT DANS LE CHAMP ACCUMULA
!     ET MULTIPLIER CHAQUE POINT PAR FACT
!
       entry racine(fact)
       it=4
       go to 1000
      entry carre(fact)
       it=6
       go to 1000
      entry absolu(fact)
       it=5
       go to 1000
       end
!
!**S/P ECRIRE SUR TAPE2 1 REC LATITUDES ET 1 REC LONGITUDES 
!
      subroutine outlalo(ip1,ip2,ip3)
      implicit none
      integer ip1,ip2,ip3
!     
!AUTEUR P. SARRAZIN AOUT 84 DRPN DORVAL P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(OUTLALO)
!        EXTRAIRE DANS LA MEMOIRE LES LATITUDES ET LES LONGITUDES
!        ET CALCUL TOUS LES PARAMETRES NECESSAIRES POUR CALL ECRIRE
!        SUR UN FICHIER STANDARD TAPE2
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN    IP1     OPTIONEL POUR USAGER DEFAUT=0
!  IN    IP2     OPTIONEL POUR USAGER DEFAUT=0
!  IN    IP3     OPTIONEL POUR USAGER DEFAUT=0
!
!IMPLICITES
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!   DEFINITION DES MACROS DE PGSM
! 
!*********************************************************************
!
      external pgsmabt,ecritur,grille2,messags
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      character(len=512) :: defo(990)
      character(len=512) :: listl(990), form
      character(len=512) :: lfn(990)
      common/carac/defo,listl,form,lfn
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      integer jjp1,jjp2,jjp3,jp1,jp2,jp3,k,ier
      integer gdll
!
!
!
      if (.not.associated(tmplat)) then
         if (message)write(6,*)'GRILLE NON DEFINIE ..GRILLE P.S.(2805)'
         ngr=8
         call grille2(3,51,55,26.,28.,381000.,350.,1)
      endif
!
!
      if (nlalo.gt.3) then
         write(6,*) ' PLUS DE 3 ARGUMENTS DANS OUTLALO'
         write(6,*) ' SEULEMENT LES 3 PREMIERS SERONT UTILISES'
      endif
!
!    INITIALISER ARGUMENTS POUR ECRITUR 
!
!  lat lon 2 dimension seulement
      k=1
!
      jp1=ip1
      jp2=ip2
      jp3=ip3
!
!     ip2,ip3=0
      if (nlalo.eq.1) then
         jp2=0
         jp3=0
      endif
!
      if (nlalo.eq.2) jp3=0
!
!
!      if (mode.ne.1) then
!         if (message) then
!            write(6,*)' OUTLALO NE PEUT ECRIRE SUR FICHIER NON STD'
!         endif
!         return
!      endif
!     
      jjp1=MIN(32767,MAX(jp1,0))
      jjp2=MIN(32767,MAX(jp2,0))
      jjp3=MIN(4095,MAX(jp3,0))
!     
      ier = gdll(gdout, tmplat, tmplon)
      call ecritur(tmplat,npack,jdate,0,0,li,lj,k,jjp1,jjp2,jjp3,      'C ','LA  ','LATITUDES   ',cgrtyp,lg1,lg2,lg3,lg4)
!     
!     INITIALISER POUR LONGITUDES
!     
      call ecritur(tmplon,npack,jdate,0,0,li,lj,k,jjp1,jjp2,jjp3,      'C ','LO  ','LONGITUDES  ',cgrtyp,lg1,lg2,lg3,lg4)
!     
!     
      return
      end
!
!**S/P PAIRVCT  REMPLACE OU AJOUTE NOM AU DICTIONNAIRE COMMON/PAIR/...
!
      subroutine pairvct(nomusag, varuu, varvv, varmodule, vardir)
      implicit none
!
!AUTEUR P. SARRAZIN DORVAL QUE CANADA FEV 87
!
!REVISION 4.0.2
!   CONVERSION DES VARIABLES HOOLERITH EN CARACTERE
!REVISION 5.6.1
!   INCLUSION DE LA VARIABLE WD - DIRECTION DU VENT Y.Chartier - Aout 1996  
!
!LANGAGE RATFOR
!
!OBJET(PAIRVCT)
!          REMPLACE OU AJOUTE DANS LA TABLE PAIRE DU COMMON/PAIR/..
!          POUR REFERENCE PAR L'USAGER QUI PERMET CERTAINES INTERPOLATIONS
!          DE VARIABLES PAIRES. 2 SETS DE VARIABLES PAIRES INITIALISE 
!          DANS PGSM UU,VV  US,VS.
!
!LIBRAIRIES
!
!          - SOURCE  PGSM ID=ARMNSRC     MFA
!          - OBJET PGSMLIB,ID=ARMNPJS    XMP
!
!ARGUMENTS
!
!  IN    NOMUSAG  NOM DE L'USAGER DONNE PAR LA DIRECTIVE PAIRES
!  IN    VARUU  NOM DE 2 CARACTERES DE LA PREMIERE VARIABLE PAIRE
!  IN    VARUU   NOM DE 2 CARACTERES DE LA DEUXIEME VARIABLE PAIRE
!  IN    VARMODULE  SI VARMODULE .NE.0 NOM DE 2 CARACTERES IDENTIFIANT
!                 LE CHAMP DE SORTIE VITESSE DU VENT  EX:"UV"
!                 SI VARMODULE.EQ.0 INTERPOLATION DE 2 CHAMPS AVEC ORIENTATION
!                 GEOGRAPHIQUE NOM DU PREMIER CHAMP POUR LA SORTIE=VARUU
!                 NOM DU DEUXIEME CHAMP=VARVV
!
!APPEL
!         - VIA DIRECTIVE PAIRES(NOMUSAG,VARUU,VARVV,VARMODULE)
!
!MESSAGE
!         - 'VERIFIER NOMBRE D ARGUMENTS DIRECTIVE PAIRES(3 OU 4 ARGS)'
!           'PAIRES DEJA INITIALISE'
!           'PAIRES("VENT","UU","VV","UV")
!           'PAIRES("UV","UU","VV","0") 
!           'PAIRES("VENTUVS","US","VS","UV")
!           'PAIRES("UVS","UU","VV","0")
!
!IMPLICITES
!MODULES
!
!
!- - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
!
!
      integer nomusag(2),varuu,varvv,varmodule,vardir
      character*8 cnomusr
      character*4 cvaruu, cvarvv, ccontrl, cvarwd
      integer i, nw
      integer argdims
      external fstcvt
      integer  fstcvt
      if (npairuv.lt.3.or.npairuv.gt.5) then
         write(6,*)'  VERIFIER ARGUMENTS DIRECTIVE PAIRES(3 OU 4 ARGS)'
         write(6,*)' PAIRES DEJA INITIALISEES'
         write(6,*)' PAIRES("VENT","UU","VV","UV")'
         write(6,*)' PAIRES("UV","UU","VV","0")'
         write(6,*)' PAIRES("VENTUVS","UU","VV","UV")'
         write(6,*)' PAIRES("VENTUVS","UU","VV","0")'
         return
      endif
!
!   VERIFI SI NOM EXISTE DANS LA TABLE SI OUI ON REMPLACE
!   SI NON ON AJOUTE SI LA TABLE N'EST PAS PLEINE 
!
      nw = min(argdims(1), 2)
      write (cnomusr, 100) (nomusag(i), i=1,nw)
 100  format(2a4)
      write (cvaruu, 200) varuu
      write (cvarvv, 200) varvv
      write (ccontrl, 200) varmodule
      write (cvarwd, 200) vardir
 200  format(a4)
      if (varmodule.eq.0) ccontrl = '??'
      if (vardir.eq.0) cvarwd = '??'
      write (6, *) 'PAIRES: ',cnomusr, cvaruu, cvarvv, ccontrl, cvarwd
      call pairvc2(cnomusr, cvaruu, cvarvv, ccontrl, cvarwd)
      return
      end
!**S/P PAIRVC2  REMPLACE OU AJOUTE NOM AU DICTIONNAIRE COMMON/PAIR/...
      subroutine pairvc2(cnomusr,cvaruu,cvarvv,ccontrl,cvarwd)
!
!AUTEUR P. SARRAZIN DORVAL QUE CANADA FEV 87
!
!REVISION 4.0.2
!   CONVERSION DES VARIABLES HOLLERITH EN CARACTERE
!
!
!LANGAGE RATFOR
!
!     OBJET(PAIRVCT)
!          REMPLACE OU AJOUTE DANS LA TABLE PAIRE DU COMMON/PAIR/..
!          POUR REFERENCE PAR L'USAGER QUI PERMET CERTAINES INTERPOLATIONS
!          DE VARIABLES PAIRES. 2 SETS DE VARIABLES PAIRES INITIALISE 
!          DANS PGSM UU,VV  US,VS.
!
!LIBRAIRIES
!
!          - SOURCE  PGSM ID=ARMNSRC     MFA
!          - OBJET PGSMLIB,ID=ARMNPJS    XMP
!
!ARGUMENTS
!
!  IN    NOMUSAG  NOM DE L'USAGER DONNE PAR LA DIRECTIVE PAIRES
!  IN    VARUU  NOM DE 2 CARACTERES DE LA PREMIERE VARIABLE PAIRE
!  IN    VARUU   NOM DE 2 CARACTERES DE LA DEUXIEME VARIABLE PAIRE
!  IN    VARMODULE  SI VARMODULE .NE.0 NOM DE 2 CARACTERES IDENTIFIANT
!                 LE CHAMP DE SORTIE VITESSE DU VENT  EX:"UV"
!                 SI VARMODULE.EQ.0 INTERPOLATION DE 2 CHAMPS AVEC ORIENTATION
!                 GEOGRAPHIQUE NOM DU PREMIER CHAMP POUR LA SORTIE=VARUU
!                 NOM DU DEUXIEME CHAMP=VARVV
!
!APPEL
!         - VIA DIRECTIVE PAIRES(NOMUSAG,VARUU,VARVV,VARMODULE)
!
!MESSAGE
!         - 'VERIFIER NOMBRE D ARGUMENTS DIRECTIVE PAIRES(3 OU 4 ARGS)'
!           'PAIRES DEJA INITIALISE'
!           'PAIRES("VENT","UU","VV","UV")
!           'PAIRES("UV","UU","VV","0") 
!           'PAIRES("VENTUVS","US","VS","UV")
!           'PAIRES("UVS","UU","VV","0")
!
!MODULES
!
!
!-- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
implicit none
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
!
      character cnomusr*8, cvaruu*4 , cvarvv*4 , ccontrl*4 , cvarwd*4
      integer np
      logical remplac
!     
!   VERIFI SI NOM EXISTE DANS LA TABLE SI OUI ON REMPLACE
!   SI NON ON AJOUTE SI LA TABLE N'EST PAS PLEINE 
!
      remplac=.false.
      do np=1,npair
         if (cnomusr.eq.paire(np)(1:8)) then
            paire(np)( 1: 8) = cnomusr
            paire(np)( 9:12) = cvaruu
            paire(np)(13:16) = cvarvv
            paire(np)(17:20) = ccontrl
            paire(np)(21:24) = cvarwd
            remplac=.true.
            write (6, *) 'PAIRE(NP): ', paire(np)
         endif
      enddo
!     
      if (remplac)  go to 1000
!
!   SI ON N'A PAS REMPLACE DANS LA TABLE ON AJOUTE
!
      npair = npair + 1
      if (npair.gt.40) then
         write(6,666) npair,40
 666     format(1x,' TROP DE PAIRES DANS LA TABLE NPAIR=',i5,         /'   NPAIRMX=',i5)
         return
      endif
!
      paire(np)( 1: 8) = cnomusr
      paire(np)( 9:12) = cvaruu
      paire(np)(13:16) = cvarvv
      paire(np)(17:20) = ccontrl
      paire(np)(21:24) = cvarwd
      write (6, *) 'PAIRE(NP): ', paire(np)
      return
!
!
 1000 write(6,*)'  2 VARIABLES PAIRES REMPLACEES '
!     
      return
      end
!**programme pgsm
!	    programme general de sortie des modeles
!               programme utilitaire d'interpolation horizontale
!               interpoler  ajustement convectif,precip,epaisseurs
!               calcul du vent sqrt(u**2 + v**2)
!               interpolation des vecteurs u et v horizontalement
!               interpolation des 3 niveaux de nuages
!               interpolation horizontale de variables scalaires
!
!auteur    - p.sarrazin octobre 1980 drpn dorval p.q.  canada
!
!revision
!        4.0.1  - conversion au fichier standard 89  p. sarrazin
!                 modification avril 90 p. sarrazin dorval canada drpn
!
!        4.0.2  - leger nettoyage du code
!               - conversion des variables contenant des informations
!                 alphanumeriques de type hollerith a caractere
!               - elimination des macros "hcar" et "lcar"
!               - conversion appels a "lexins" par "qlxins"
!               - conversion pour cyber-910
!                 y. chartier -juillet-aout 90- drpn dorval quebec
!
!        5.0    - utilisation des nouveaux interpolateurs
!                 y. chartier - mai 1991 - drpn dorval quebec
!
!        5.1    - utilisation de fichiers d'entree lies avec "fstlnk"
!                 y. chartier - mai 1991 - drpn dorval quebec
!
!        5.2    - support des grilles source de type z
!        5.3    - conversion de RATFOR a FORTRAN
!                 Y. Chartier - aout 1995
!        5.4    - Optimisation des interpolateurs pour l'extension
!                 selective des grilles.
!        5.5    - Interpolation a partir de fichiers d'entree SQI
!                 Support des grilles Lambert
!                 Introduction de la librairie C gctpc
!        5.6    - Introduction des directives COORD et GRILLE(STATIONS)
!                 sortie(ASCII)
!        5.7    - Support des fichiers standards 98
!        6.0    - Introduction de ezscint comme interpolateur principal
!        6.8    - Support des grilles diese
!        6.9    - Support des grilles T (stereographiques generalisees)
!        7.8.2  - Reload avec librmn_015.1
!        7.8.3  - D. Bouhmemhem, Fev 2015, Reload avec librmn_015.2
!        7.8.4  - M. Valin, Avril 2015,fixed data statements, include version.inc
!        7.8.5  - M. Lepine, Nov 2015, verification des codes de retour de fnom
!
!langage   - fortran
!
!objet(pgsm)
!          interpolateur horizontal qui permet de faire des
!          interpolations (cubique,lineaire,voisin) d'une grille a
!          une autre ou d'une grille a un point
!          permet de faire des operations sur deux champs
!          (gros calculateur de poche)
!          fichier d'entre doit-etre format standard random
!          sorti standard random ou sequentiel - seq ms - random ms
!
!librairies
!		  - rmnxlib.a
!fichiers
!         - tape1  - fichier d'entree  (standard)
!         - tape2  - fichier de sortie standard..direct(writms)...sequentiel
!         - tape3  - fichier de records positionels ('^^','>>')
!         - tape5  - fichier d'entree(directives)
!         - tape6  - fichier de sortie sur imprimante
!
!----------------------------------------------------------------------------
	subroutine pgsm
      implicit none
!
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
      character(len=512) :: defo(990)
      character(len=512) :: listl(990), form
      character(len=512) :: lfn(990)
      common/carac/defo,listl,form,lfn
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
   integer heures,nhur,nheure
   common / heures/ nhur,nheure,heures(40)
!
      integer nnoms,maxnoms,nsym,nsymm
      logical ssym(256)
      common /symnom/ nnoms,maxnoms,nsym,nsymm,ssym
      character*4 noms(256)
      common /csymnom/ noms
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer ncon,nomb
      real ecarts,facts,bass,hauts
      common/tabls/ncon,nomb,ecarts(256),facts(256),bass(256),hauts(256)
      character*4 nomss(256)
      common /ctabls/ nomss
!
!
!
!
      integer champpr,nchamp,nchmp,npar
      common / champs/ nchamp, nchmp, npar,champpr(31)
!     
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   integer noenrg,numero,nbrow,numdel,istart
   common/enrege/ noenrg,numero,nbrow,numdel,istart
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nivospr,nmoy,nmo
      common / nivos/ nmoy, nmo, nivospr(31)
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
      integer nmaxlist1,nmaxlist2,wait,go
      parameter (nmaxlist1=16,nmaxlist2=16)
      parameter (wait=0,go=1)
      character*2 listnom(nmaxlist1,nmaxlist2)
      integer listniv(nmaxlist1,nmaxlist2)
      integer ntitems,nitems1(nmaxlist1),nitems2(nmaxlist2)
      common /seq/  ntitems,nitems1,nitems2
      common /cseq/ listnom
      integer ip1style, dateform
      common /styles/ ip1style, dateform
	character *8 qlxcon(128),qlxlcon(4)
	integer      qlxval(128)
	integer      qlxlval(4)
	integer ezsetopt
	external ezsetopt, heure, champ, sorti, grille2, metsym, cmetsym, convs
	external    qqqintx, setxtrap, liren, lirsr, plmnmod, pluss
	external moinse, moinss, ecrits,moyene, operat, modul2e, modul2s
	external expon, racine,alogn, absolu, carre, outlalo, foise, foiss, divisee, divises, pgcoupe
	external moysrt, imprims,chmpdif, pairvct, messags, champ_seq,qqqecho
	external qqqform,qqqident,coord,qqqfilt
!
	external ccard,fnom,exdb,qlxins,qlxinx,readlx,fstfrm,fstvoi
	external fstnbr,fstunl,fstouv
	external fclos,exfin,lrsmde,lrsmds,fstopc,fstopl,qlxopt
!
	integer exdb,exfin,fnom,fstfrm,fstvoi,fstnbr,fstopc,fstopl, fstouv
	integer i,iopc,ipose,kend,nequiv,npex,nsetin,nsetex,nlirmds,nlirmde
	real dum
        integer, parameter :: str_A=transfer("A   ",1)
        integer, parameter :: str_P=transfer("P   ",1)
        integer, parameter :: str_GZ=transfer("GZ  ",1)
        integer, parameter :: str_TT=transfer("TT  ",1)
        integer, parameter :: str_QQ=transfer("QQ  ",1)
        integer, parameter :: str_QR=transfer("QR  ",1)
        integer, parameter :: str_DD=transfer("DD  ",1)
        integer, parameter :: str_PP=transfer("PP  ",1)
        integer, parameter :: str_CC=transfer("CC  ",1)
        integer, parameter :: str_WW=transfer("WW  ",1)
        integer, parameter :: str_ES=transfer("ES  ",1)
        integer, parameter :: str_DFGZ=transfer("DFGZ",1)
        integer, parameter :: str_DFST=transfer("DFST",1)
        integer, parameter :: str_DFPR=transfer("DFPR",1)
        integer, parameter :: str_UV=transfer("UV  ",1)
        integer, parameter :: str_VENT=transfer("VENT",1)
        integer, parameter :: str_NUAG=transfer("NUAG",1)
        integer, parameter :: str_F2=transfer("F2  ",1)
        integer, parameter :: str_PN=transfer("PN  ",1)
        integer, parameter :: str_P0=transfer("P0  ",1)
        integer, parameter :: str_TS=transfer("TS  ",1)
        integer, parameter :: str_TM=transfer("TM  ",1)
        integer, parameter :: str_MT=transfer("MT  ",1)
        integer, parameter :: str_WDUV=transfer("WDUV",1)
        character(len=16) :: PGSM_VERSION
!
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	data listl/981*'IMENT:','OZSRT:','ISLL:','I.',    'L.',     'DATE.','MSGLVL.','ISENT:','IMPOS:','V'/
	data defo /981*'SCRAP', 'TAPE2', 'TAPE4','$INPUT','$OUTPUT','OPRUN','INFORMS','ISENT_SCRAP','IMPOS_SCRAP','OUI'/
	data lfn  /981*'SCRAP', 'TAPE2', 'TAPE4','$INPUT','$OUTPUT','NON',  'INFORMS','ISENT_SCRAP','IMPOS_SCRAP','NON'/
        data form/'(A8)'/
!
	data nheure,  heures, nnoms,  npack,  nhur, nomb, ichck         /0,  40*-2, 0,      -16,    1,    0,    0/
!
	data nomss /256*'  '/
	data ecarts,     facts,     pose,     ixlat, ixlon       /256*0.0, 256*1.0, .false., 0, 0 /
!
   	data nchamp,  ngr,  nsort,   nchmp,   icnt, nlalo         /  1,    0,     0,        1,      0,     0 /
!
	data valid, voire,   voirs, message,seldat       /.false.,.false., .false.,.true.,.false.  /
!
	data numero,  numdel,  iset,  nbrow,  ip4        / 1,           1,    -2,      0,    0    /
!
	data paire(1) /  'VENT    UU  VV  UV      ' /
	data paire(2) /  'UV      UU  VV  ??      ' /
	data paire(3) /  'VENTUVS US  VS  UV      ' /
	data paire(4) /  'UVS     US  VS  ??      ' /
	data paire(5) /  'WDUV    UU  VV  UV  WD  ' /
	data paire(6) /  'WDUD    UD  VD  UV  WD  ' /
	data paire(7) /  '!#@$!#@$>>  ^^  >>  ^^  ' /
!
	data unefois,once,vvent/.false.,.false.,.false./
!
	data cnomqq, cnomqr, cnommt /'QQ', 'QR', 'MT'/
!
	data printen,printsr,mtdone/.false.,.false.,.false./
!
 	data nis,njs,nif,njf,ninc,njnc,if9/1,1,1000,1000,10,10,0/
!
  	data niis,njjs,niif,njjf,niinc,njjnc/1,1,1000,1000,10,10/
!
  	data if7,if8,npairuv,npair/0,0,4,7/
!
  	data clatmin,clatmax,clonmin,clonmax,ncoords/-90.0, +90.0, 0.0, 360.0, 0/
	data qlxcon( 1) /'ZON'     /  qlxval( 1) /      1 /
	data qlxcon( 2) /'MER'     /  qlxval( 2) /      2 /
	data qlxcon( 3) /'TOUT'    /  qlxval( 3) /     -1 /
	data qlxcon( 4) /'ALL'     /  qlxval( 4) /     -1 /
	data qlxcon( 5) /'COMTEUR' /  qlxval( 5) /   4095 /
	data qlxcon( 6) /'IMPRIM'  /  qlxval( 6) / 999999 /
	data qlxcon( 7) /'STD'     /  qlxval( 7) /   gr_a /
	data qlxcon( 8) /'MS'      /  qlxval( 8) /      2 /
	data qlxcon( 9) /'SEQ'     /  qlxval( 9) /      3 /
	data qlxcon(10) /'R'       /  qlxval(10) /      1 /
	data qlxcon(11) /'A'       /  qlxval(11) /     -1 /
	data qlxcon(12) /'NORD'    /  qlxval(12) /      1 /
	data qlxcon(13) /'SUD'     /  qlxval(13) /      2 /
	data qlxcon(14) /'GLOBAL'  /  qlxval(14) /      0 /
	data qlxcon(15) /'LATLON'  /  qlxval(15) / gr_latlon /
	data qlxcon(16) /'PS'      /  qlxval(16) / gr_ps   /
	data qlxcon(17) /'TAPE4'   /  qlxval(17) / gr_tape4 /
	data qlxcon(18) /'GAUSS'   /  qlxval(18) / gr_g /
	data qlxcon(19) /'STDB'    /  qlxval(19) / gr_b /
	data qlxcon(20) /'TAPE1'   /  qlxval(20) / gr_tape1 /
	data qlxcon(21) /'TAPE2'   /  qlxval(21) / gr_tape2 /
	data qlxcon(22) /'XYLIS'   /  qlxval(22) / gr_xylis /
	data qlxcon(23) /'XYDIR'   /  qlxval(23) / gr_xydir /
	data qlxcon(24) /'LLDIR'   /  qlxval(24) / gr_lldir /
	data qlxcon(25) /'LLLIST'  /  qlxval(25) / gr_lllist/
	data qlxcon(26) /'ANAL'    /  qlxval(26) / str_A/
	data qlxcon(27) /'PREV'    /  qlxval(27) / str_P   /
	data qlxcon(28) /'Z'       /  qlxval(28) / str_GZ  /
	data qlxcon(29) /'T'       /  qlxval(29) / str_TT  /
	data qlxcon(30) /'Q'       /  qlxval(30) / str_QQ  /
	data qlxcon(31) /'QR'      /  qlxval(31) / str_QR  /
	data qlxcon(32) /'D'       /  qlxval(32) / str_DD  /
	data qlxcon(33) /'PP'      /  qlxval(33) / str_PP  /
	data qlxcon(34) /'CC'      /  qlxval(34) / str_CC  /
	data qlxcon(35) /'W'       /  qlxval(35) / str_WW  /
	data qlxcon(36) /'ES'      /  qlxval(36) / str_ES  /
	data qlxcon(37) /'EPAIS'   /  qlxval(37) / str_DFGZ/
	data qlxcon(38) /'MAC'     /  qlxval(38) / str_DFST/
	data qlxcon(39) /'PCP'     /  qlxval(39) / str_DFPR/
	data qlxcon(40) /'UV'      /  qlxval(40) / str_UV  /
	data qlxcon(41) /'VENT'    /  qlxval(41) / str_VENT/
	data qlxcon(42) /'NUAGES'  /  qlxval(42) / str_NUAG/
	data qlxcon(43) /'ECM'     /  qlxval(43) / str_F2  /
	data qlxcon(44) /'PNM'     /  qlxval(44) / str_PN  /
	data qlxcon(45) /'PSURF'   /  qlxval(45) / str_P0  /
	data qlxcon(46) /'TSRF'    /  qlxval(46) / str_TS  /
	data qlxcon(47) /'TMER'    /  qlxval(47) / str_TM  /
	data qlxcon(48) /'MT'      /  qlxval(48) / str_MT  /
	data qlxcon(49) /'VOISIN'  /  qlxval(49) /   100 /
	data qlxcon(50) /'LINEAIR' /  qlxval(50) /     1 /
	data qlxcon(51) /'CUBIQUE' /  qlxval(51) /     3 /
	data qlxcon(52) /'ABORT'   /  qlxval(52) /    13 /
	data qlxcon(53) /'MINIMUM' /  qlxval(53) /     5 /
	data qlxcon(54) /'MAXIMUM' /  qlxval(54) /     4 /
	data qlxcon(55) /'GEF'     /  qlxval(55) / gr_gem /
	data qlxcon(56) /'GEM'     /  qlxval(56) / gr_gef /
	data qlxcon(57) /'WAIT'    /  qlxval(57) /     0 /
	data qlxcon(58) /'GO'      /  qlxval(58) /     1 /
	data qlxcon(59) /'GRIB'    /  qlxval(59) / gr_grib/
	data qlxcon(60) /'FORMATEE'/  qlxval(60) /     5 /
	data qlxcon(61) /'STATIONS'/  qlxval(61) / gr_stations /
	data qlxcon(62) /'ADD'     /  qlxval(62) /     1 /
	data qlxcon(63) /'RESET'   /  qlxval(63) /     0 /
	data qlxcon(64) /'EST'     /  qlxval(64) /     3 /
	data qlxcon(65) /'OUEST'   /  qlxval(65) /     4 /
	data qlxcon(66) /'NONE'    /  qlxval(66) /     5 /
	data qlxcon(67) /'NOMVAR'  /  qlxval(67) /     1 /
	data qlxcon(68) /'TYPVAR'  /  qlxval(68) /     2 /
	data qlxcon(69) /'ETIKET'  /  qlxval(69) /     3 /
	data qlxcon(70) /'IP01'    /  qlxval(70) /     4 /
	data qlxcon(71) /'IP02'    /  qlxval(71) /     5 /
	data qlxcon(72) /'IP03'    /  qlxval(72) /     6 /
	data qlxcon(73) /'DATEO'   /  qlxval(73) /     7 /
	data qlxcon(74) /'DATEV'   /  qlxval(74) /     8 /
	data qlxcon(75) /'LAT'     /  qlxval(75) /    12 /
	data qlxcon(76) /'LON'     /  qlxval(76) /    13 /
	data qlxcon(77) /'NI'      /  qlxval(77) /     9 /
	data qlxcon(78) /'NJ'      /  qlxval(78) /    10 /
	data qlxcon(79) /'NK'      /  qlxval(79) /    11 /
	data qlxcon(80) /'WDUV'    /  qlxval(80) /str_WDUV/
	data qlxcon(81) /'ON'      /  qlxval(81) /     1 /
	data qlxcon(82) /'OFF'     /  qlxval(82) /     0 /
	data qlxcon(83) /'VERBOSE' /  qlxval(83) /     1 /
	data qlxcon(84) /'LECTURE' /  qlxval(84) /     1 /
	data qlxcon(85) /'ECRITURE'/  qlxval(85) /     2 /
	data qlxcon(86) /'SEQWPRM' /  qlxval(86) /     4 /
	data qlxcon(87) /'AUCUNE'  /  qlxval(87) /     0 /
	data qlxcon(88) /'COMME'   /  qlxval(88) / gr_comme /
	data qlxcon(89) /'LIKE'    /  qlxval(89) /    12 /
	data qlxcon(90) /'IP1A'    /  qlxval(90) / 65001 /
	data qlxcon(91) /'IP1B'    /  qlxval(91) / 65002 /
	data qlxcon(92) /'IP2A'    /  qlxval(92) / 65003 /
	data qlxcon(93) /'IP2B'    /  qlxval(93) / 65004 /
	data qlxcon(94) /'IP3A'    /  qlxval(94) / 65005 /
	data qlxcon(95) /'IP3B'    /  qlxval(95) / 65006 /
	data qlxcon(96) /'FENTREE' /  qlxval(96) /     1 /
	data qlxcon(97) /'FSORTIE' /  qlxval(97) /     2 /
	data qlxcon(98) /'LOCAL'   /  qlxval(98) /     1 /
	data qlxcon(99) /'IP1'     /  qlxval(99) /     4 /
	data qlxcon(100)/'IP3'     /  qlxval(100)/     6 /
	data qlxcon(101)/'STEREO'  /  qlxval(101)/ gr_stereo /
	data qlxcon(102)/'IPUN'    /  qlxval(102)/     4 /
	data qlxcon(103)/'IPDEUX'  /  qlxval(103)/     5 /
	data qlxcon(104)/'IPTROIS' /  qlxval(104)/     6 /
	data qlxcon(105)/'IPONE'   /  qlxval(105)/     4 /
	data qlxcon(106)/'IPTWO'   /  qlxval(106)/     5 /
	data qlxcon(107)/'IPTHREE' /  qlxval(107)/     6 /
!               KIND =0, p est en hauteur (m) par rapport au niveau de la mer
!               KIND =1, p est en sigma (0.0 -> 1.0)
!               KIND =2, p est en pression (mb)
!               KIND =3, p est un code arbitraire
!               KIND =4, p est en hauteur (M) par rapport au niveau du sol
!               KIND =5, p est en coordonnee hybride
!               KIND =6, p est en coordonnee theta
!               KIND =15, rererve (entiers)
!               KIND =21, p est en GalChen
	data qlxcon(108)/'METERS'  /  qlxval(108)/ -1000 /
	data qlxcon(109)/'SIGMA'   /  qlxval(109)/ -1001 /
	data qlxcon(110)/'MBAR'    /  qlxval(110)/ -1002 /
	data qlxcon(111)/'OTHER'   /  qlxval(111)/ -1003 /
	data qlxcon(112)/'METERAGL'/  qlxval(112)/ -1004 /
	data qlxcon(113)/'HYBRID'  /  qlxval(113)/ -1005 /
	data qlxcon(114)/'THETA'   /  qlxval(114)/ -1006 /
	data qlxcon(115)/'GALCHEN' /  qlxval(115)/ -1021 /
	data qlxcon(116)/'OLDSTYLE'/  qlxval(116)/ 3 /
	data qlxcon(117)/'NEWSTYLE'/  qlxval(117)/ 2 /
	data qlxcon(118)/'STAMP'   /  qlxval(118)/ 0 /
	data qlxcon(119)/'YMDHMS'  /  qlxval(119)/ 1 /
	data qlxcon(120)/'ISO8601' /  qlxval(120)/ 2 /
        data qlxcon(121)/'FAST'    /  qlxval(121)/ 1 /
        data qlxcon(122)/'BEST'    /  qlxval(122)/ 2 /
        data qlxcon(123)/'MOYENNE' /  qlxval(123)/ 4 /
        data qlxcon(124)/'GRIDAVG' /  qlxval(124)/ 4 /
        data qlxcon(125)/'SPHRAVG' /  qlxval(125)/ 5 /
        data qlxcon(126)/'EXCLUDE' /  qlxval(126)/ 31/
        data qlxcon(127)/'ORIGIN'  /  qlxval(127)/ 1023 /
        data qlxcon(128)/'RESV128' /  qlxval(128)/ 0 /
	data(qlxlcon(i),i=1,2)/'OUI', 'NON'/
	data(qlxlval(i),i=1,2)/1,0/
!        integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
   data idx_ozsrt  /982/  idx_isll  /983/  idx_i      /984/ idx_l /985/ idx_date /986/
   data idx_msglvl /987/  idx_isent /988/  idx_impos  /989/ idx_v /990/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!           listl=position  iment(tape1 standard),isll(tape4 sequentiel)
!                 ozsrt(tape2 - standard - seq file - random ms)
!           defo=liste des defauts pour iment,isll,ozsrt,i,l
!           lfn=liste que l usager propose pour remplacer
!           6=nombre de lfn
!           nequiv=nombre d'equivalence output de ccard
!
!
      include 'version.inc'
	nequiv=-1
	lnkdiun = 0
   lnkdiun(1) = 1
   lnkdiun(idx_ozsrt) = 2
	call ccard(listl,defo,lfn,990,nequiv)
	ier = fnom(5,lfn(idx_i),'SEQ',0)
	ier = fnom(6,lfn(idx_l),'SEQ',0)
!
!          imprime boite debut du programme
!
	jdate= exdb(' PGSM  ',PGSM_VERSION,  lfn(idx_date))
 
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	call qqqfilt(1,0,0,0)
	call qqqfilt(2,0,0,0)
        call chk_tmpdir
	if (lfn(1)(1:5).ne.'SCRAP'.and.lfn(idx_isent)(1:11).ne.'ISENT_SCRAP') then
	  print *,'***************************************************'
	  print *,'* ON NE PEUT MELANGER LES FICHIERS D ENTREE       *'
	  print *,'* SEQUENTIELS ET RANDOM                           *'
	  print *,'***************************************************'
	  jdate= exfin('  PGSM  ', 'ABORT', 'NON')
	  call qqexit(13)
	endif
	if (lfn(idx_isll).ne.'TAPE4') then
	    iun_isll=0
	    ier = fnom(iun_isll,lfn(idx_isll)(1:5),'FMT+SEQ+R/O',0)
	    isll_input = 1
	endif
	if (lfn(idx_isent)(1:11).ne.'ISENT_SCRAP') then
	  inputmod = SEQUENTIEL
	else
	  inputmod = RANDOM
	endif
	if (lfn(idx_impos)(1:11).ne.'IMPOS_SCRAP') then
	  ier = fnom(lnkdiun(idx_impos),lfn(idx_impos),'RND+OLD+R/O',0)
	  ier = fstouv(lnkdiun(idx_impos),'RND')
	endif
	if (lfn(idx_v).ne.'NON') then
	  ier = ezsetopt('verbose', 'yes')
	endif
	if (inputmod.eq.RANDOM) then
	  niun = 1
 100	  if (lfn(niun) .ne.'SCRAP') then
	    niun = niun+1
	    goto 100
	  endif
	  niun = niun - 1
	  if (niun .lt. 1) then
	    print *,'***************************************************'
	    print *,'* AUCUN FICHIER D''ENTREE DONNE EN ARGUMENT !!!'
	    print *,'***************************************************'
	    jdate= exfin('  PGSM  ', 'ABORT', 'NON')
	    call qqexit(13)
	  endif
	  do i=1, niun
	    ier = fnom(lnkdiun(i),lfn(i),'STD+RND+OLD+R/O+REMOTE',0)
	    if (ier .lt. 0) then
	      print *,'************************************************'
	      print *, '* PROBLEME D''OUVERTURE AVEC LE FICHIER ',lfn(i)
	      print *,'************************************************'
	      jdate= exfin('  PGSM  ', 'ABORT', 'NON')
	      call qqexit(13)
	    endif
	  enddo
	else
	  niun = 1
	  ier = fnom(lnkdiun(1),lfn(idx_isent),'STD+SEQ+OLD+R/O+REMOTE',0)
	  if (ier .lt. 0) then
	    print *,'************************************************'
	    print *, '* PROBLEME D''OUVERTURE AVEC LE FICHIER ',lfn(idx_isent)
	    print *,'************************************************'
	    jdate= exfin('  PGSM  ', 'ABORT', 'NON')
	    call qqexit(13)
	  endif
	endif
!
!
!
	mtype =   32
	maxnoms = 256
	call initseq
!
!
!  initialise les dictionnaires
!
	call qlxopt ('CARMOT', 4)
	call qlxinx (sorti,'SORTIE', nsort,0103,2)
!                        3 appels reconnus  1=sortie(std,noenrg) noenrg>=2
!                                           2=sortie(ms,noenrg,jwrit)
!                                           3=sortie(seq)
!
	call qlxinx (heure,'HEURE',nheure, 0140,2)
	call qlxinx (heure, 'IP2',nheure, 0140,2)
!                     2 appels  1=heure(00,12,24,25.....max20) minimum 1
!                               2=champ(mac,00,06) minimum 2 pour
!                               2=champ(pcp,00,06) minimum 2 pour
!                               accumulateur d"ajustement ou precipitation
!
	call qlxinx (qqqintx,'SETINTX',nsetin, 0101,2)
!
!                        appel - setintx(voisin) avec le plus proche
!                                setintx(lineair) interpolation lineaire
!                                setintx(cubique) interpolation cubique(defaut)
!
	call qlxinx (setxtrap,'EXTRAP',nsetex, 0101,2)
!
!                        appel - setintx(voisin) avec le plus proche
!                                setintx(lineair) interpolation lineaire
!                                setintx(cubique) interpolation cubique(defaut)
!
	call qlxinx (champ,'CHAMP',nchamp, 0131,2)
!                     appel - champ(z,niveau)  niveau=1000,850,.......
!                           - champ(t,niveau)  niveau=1000,850,.......
!                           - champ(q,niveau)  niveau=1000,850,.......
!                           - champ(d,niveau)  niveau=1000,850,.......
!                           - champ(w,niveau)  niveau=1000,850,.......
!                           - champ(es,niveau)  niveau=1000,850,.......
!                           - champ(uv,niveau)  niveau=1000,850,.......
!                           - champ(uvs)  pas de niveau vent de surface
!                           - champs(ventuvs) voir directive paires(.....
!                           - champ(vent,niveau) niveau=1000,850,.......
!                           - champ(nuage)  nuage bas,moyen,haut
!                                   rec 1=bas  rec 2= moyen  rec 3=haut
!                           - champ(ecm)  epaisseur de la couche limite
!                           - champ(pnm)  pression au niveau de la mer
!                           - champ(psurf)  pression a la surface
!                           - champ(ts)  temperature a la surface
!                           - champ(epais,niveau1,niveau2) niveau2 - niveau1
!                           - champ(mac,heure1,heure2)  heure2 - heure1
!                           - champ(pcp,heure1,heure2)  heure2 - heure1
!
	call qlxinx (chmpdif,'CHMPDIF',npar,0508,2)
!
!                 appel - chmpdif (noment,nomsrt,ip1tab,ip2tab,ip3tab)
!                    ex:  chmpdif ("gz","dz",[1000,500],12,0)
!                         z500mb - z1000mb  a 12hr
!                         fichier de sorti aura ip1=1000, ip2=500,ip3=12
!                    ex:  chmpdif ("gz","dz",1000,[6,12,18,24],0)
!                         z1000mb 6hr - z1000mb  a 12hr
!                         fichier de sorti aura ip1=1000, ip2=6, ip3=12
!                    ex:  chmpdif ("gz","dz",1000,6,[1,2,3,4])
!                         z1000mb 6hr ip3=1 - z1000mb  6hr ip3=2
!                         fichier de sorti aura ip1=1000, ip2=1, ip3=2
!
	call qlxinx (champ_seq,'CHAMPSEQ',npar,0303,2)
!
!                 appel - champseq(['GZ','TT','UU'],[1000,850,500],WAIT)
!                 appel - champ_seq(' ',[1000,850,500],WAIT)
!                 appel - champ_seq(['GZ','TT','UU'],-1,GO)
	call qlxinx (convs, 'CONV',ncon, 0305,2)
!
!                 appel - conv(nom, ecart, facteur, bas, haut) directive
!                         conv("ts", -273.16, 1.0,-280.0, -250.0)
!                         routine conver dans ecriture soustrait
!                         273.16 au champ et multiplit par 1.0
!                         enleve toutes les valeurs plus petites que -280
!                         enleve toutes les valeurs plus grandes que -250
!                         avant d ecrire le champ
	call qlxinx (grille2,'GRILLE',  ngr, 0109,2)
!          8 appels a grille    1=grille(std,nni,nnj,lg1)
!                                 std=standard lat lon
!                                 nni=nombre de pts est-ouest
!                                 nnj=nombre de pts nord-sud
!                                 lg1=0  global
!                                    =1  hem nord
!                                    =2  hem sud
!                               2=grille(latlon,nni,nnj,lat0,lon0,dlat,dlon)
!                                 latlon=grille lat lon
!                                 nni= nombre de pts est-ouest
!                                 nnj= nombre de pts nord-sud
!                                 lat0=premiere lat du coin degree
!                                 lon0=premiere lon du coin degree
!                                 dlat=espacement entre latitude  (degree)
!                                 dlon=espacement entre longitude (degree)
!
!                               3=grille(ps,nni,nnj,pi,pj,d60,dgrw)
!                                 ps  =polaire stereographique
!                                 nni =nombre pts est-ouest (dir i)
!                                 nnj =nombre de pts nord-sud (dir j)
!                                 pi  =position du pole nord(pi=26)
!                                 pj  = position du pole nord(pj=28)
!                                 d60 =distance en metres entre les pts
!                                      a 60 degrees nord (latitude)
!                                 drgw=angle entre l"axe x et greewich
!
!                               4=grille(tape4,nni,nnj,ip1,ip2,ip3)
!                                 tape4=fichier contenant nni*nnj(lat-lon)
!                                 nni  =nombre de pts est-ouest
!                                 nnj  =nombre de pts nord-sud
!                                 ip1  =definit par usager
!                                 ip2  =definit par usager
!                                 ip3  =definit par usager
!
!                               5=grille(stdb,nni,nnj,hem)
!                                 stdb=standard b
!                                 nni  =nombre de pts est-ouest
!                                 nnj  =nombre de pts nord-sud
!                                 hem  =hemisphere 0=global
!                                                  1=nord
!                                                  2=sud
!
!                               6=grille(gauss,nni,nnj,hem)
!                                 gauss=grille gaussienne lat-lon
!                                 nni  =nombre de pts est-ouest
!                                 nnj  =nombre de pts nord-sud
!                                 hem  =hemisphere 0=global
!                                                  1=nord
!                                                  2=sud
!
!                               7=grille(tape1,ip1,ip2,ip3,ip4,nord/sud)
!                                 tape1=lit sur fichier 1 lat-lon ou xy
!                                 ip1=valeur 0-32767
!                                 ip2=valeur 0-32767
!                                 ip3=valeur 0-4095
!                                 ip4=valeur "xydir" ou "llist"
!                                    =valeur "lldir" ou "xylis"
!
!                               8=grille(tape2,ip1,ip2,ip3,ip4,nord/sud)
!                                 tape2 lit sur fichier 2 lat-lon ou xy
!                                 ip1=valeur 0-32767
!                                 ip2=valeur 0-32767
!                                 ip3=valeur 0-4095
!                                 ip4=valeur "xydir" ou "llist"
!                                    =valeur "lldir" ou "xylis"
!
	call qlxinx (lrsmde,'LIRMODE',nlirmde,0708,2)
	call qlxinx (lrsmds,'LIRMODS',nlirmds,0708,2)
!
!                  lrsmde(nomvar,typvar,date,niveau,heure,ip3,etiquet)
!                  lrsmds(nomvar,typvar,date,niveau,heure,ip3,etiquet)
!
	call qlxinx (metsym,'METSYM',  nsym, 0202,2)
!                               metsym(z,oui)
!                               z  =geopotentiel "gz"
!                               oui=symetrique
!
	call qlxinx (outlalo,'OUTLALO', nlalo, 0108,2)
!     outlalo(ip1,ip2,ip3,nomlat,nomlon,grtyp,etiklat,etiklon)
!             ip1=valeur 0-32767
!             ip2=valeur 0-32767
!             ip3=valeur 0-4095
!             nomlat=nom du champ de latitude 2 car
!             nomlon=nom du champ de longitude 2 car
!             grtyp=type de grille
!             etiklat=nom de l'etiquette latitude
!             etiklon=nom de l'etiquette longitude
!
	call qlxinx (pairvct, "PAIRES",npairuv, 0305,2)
!      ex: paires("uv","uu","vv",0) vecteur "uu","vv" geographique
!                                     niveau donne par champ
!      ex: paires("ventuvs","us","vs","uv") vitesse du vent a la surface
!      ex: paires("uvs","us","vs",0) vecteurs du vent a la surface
!
   call qlxinx (pgcoupe,'MOYENT', nmoy, 0232,2)
   call qlxinx (moysrt,'MOYSRT', nmoy, 0232,2)
   call qlxinx (liren,'LIREE', nlire, 0708,2)
	call qlxinx (lirsr,'LIRES', nlire, 0708,2)
	call qlxinx(plmnmod,'PLUSE', najou, 0707,2)
	call qlxinx (pluss,'PLUSS', najou, 0707,2)
	call qlxinx (foise,'FOISE', multp, 0707,2)
	call qlxinx (foiss,'FOISS', multp, 0707,2)
	call qlxinx (divisee,'DIVE', multp, 0707,2)
	call qlxinx (divises,'DIVS', multp, 0707,2)
	call qlxinx (moinse,'MOINSE', nenle, 0707,2)
	call qlxinx (moinss,'MOINSS', nenle, 0707,2)
	call qlxinx (moyene,'MOYENE', nmoys, 0101,2)
   call qlxinx (ecrits,'ECRITS',  necrt, 0814,2)
   call qlxinx (modul2e,'MODUL2E', nmod, 0707,2)
   call qlxinx (modul2s,'MODUL2S', nmod, 0707,2)
   call qlxinx (racine,'RACINE', nraci, 0101,2)
   call qlxinx (operat, 'PFOIS', npfo, 0303,2)
   call qlxinx (expon, 'EXPON', npex, 0101,2)
   call qlxinx (alogn, 'ALOGN', npex, 0101,2)
   call qlxinx (absolu,  'ABSOLU', npex, 0101,2)
   call qlxinx (carre, 'CARRE', npex, 0101,2)
!
	call qlxinx (qqqecho, 'ECHO',   dum, 0101, 2)
	call qlxinx (qqqident,'IDENT',  npar, 0103, 2)
	call qlxinx (qqqform, 'FORMAT', dum, 0101, 2)
	call qlxinx (coord,   'COORD',  dum, 0202, 2)
	call qlxinx (qqqfilt, 'FILTRE', dum, 0204, 2)
!
  	call qlxins (npack,  'COMPAC',  dum, 1, 1)
	call qlxins (message,'MESSAGE', dum, 1, 1)
	call qlxins (numdel, 'DELTA',   dum, 1, 1)
	call qlxins (typeent,'TYPEENT', dum, 1, 1)
	call qlxins (typesrt,'TYPESRT', dum, 1, 1)
	call qlxins ( voire, 'VOIRENT', dum, 1, 1)
	call qlxins ( voirs, 'VOIRSRT', dum, 1, 1)
   call qlxins ( pose,  'PAUSE',   dum, 1, 1)
	call qlxins ( userdate,  'DATE',    dum, 3, 1)
	call qlxins (seldat, 'OPDAT',   dum, 1, 1)
	call qlxins (printen,'PRINTEN', dum,7, 1)
	call qlxins (printsr,'PRINTSR', dum,7, 1)
	call qlxins (etikent,'ETIKENT', nwetike, 3, 1)
   call qlxins (masque, 'MASQUE',  dum, 1, 1)
	call qlxins (etiksrt,'ETIKSRT', nwetiks, 3, 1)
	call qlxins (numero, 'ENREG',   dum, 1, 1)
	call qlxins (ip2srt, 'IP2SRT',  dum, 1, 1)
	call qlxins (ip3ent, 'IP3ENT',  dum, 1, 1)
	call qlxins (ip3srt, 'IP3SRT',  dum, 1, 1)
	call qlxins (unefois,'UNEFOIS',  dum, 1, 1)
	call qlxins (once,   'ONCE',  dum, 1, 1)
	call qlxins (diese,  'DIESE',dum,1,1)
	call qlxins (ip1style, 'IP1STYLE', dum, 1, 1)
	call qlxins (dateform, 'DATEFORM', dum, 1, 1)
   call qlxins (compression_level, 'COMPRESS', dum, 1, 1)
	do i=1,127
	   call qlxins(qlxval(i), qlxcon(i), dum, 1, 0)
	enddo
	do i=1,2
	   call qlxins(qlxlval(i), qlxlcon(i), dum, 1, 0)
	enddo
!
!   defaut pour lire fichier d'entre
!
	typeent = -1
	etikent(1) = -1
	etikent(2) = -1
	ip3ent = -1
	userdate  = -1
	date2 = -1
	date3 = -1
	diese = 1
	ip1style = 2
	dateform = 1
!
!   defaut pour fichier de sorti
!
	ip3srt= -1
	ip2srt=-1
	etiksrt(1) = -1
	etiksrt(2) = -1
	etiksrt(3) = -1
	typesrt= -1
   compression_level = 0
   masque = 0
!
!
!    initialiser avec .true. champ symetrique
!
	nsym = 2
	call cmetsym('GZ',.true.)
	call cmetsym('TT',.true.)
	call cmetsym('DD',.true.)
	call cmetsym('WW',.true.)
	call cmetsym('ES',.true.)
	call cmetsym('F2',.true.)
	call cmetsym('PN',.true.)
	call cmetsym('PS',.true.)
	call cmetsym('TS',.true.)
!
	call cmetsym('QQ',.false.)
!
!
!    directives de l'usager
!
!
!    initialisation parametres de sortie pour fichier formate
	call initid
	iopc= fstopc('MSGLVL',lfn(idx_msglvl),.false.)
	ier = fstopl('REDUCTION32',.true.,.false.)
	ipose= 0
	call readlx(5,kend,ipose)
!
!   initialise variable de printsr
!
	if (associated(tmplat)) deallocate(tmplat)
	if (associated(tmplon)) deallocate(tmplon)
	if (mode.eq.1) then
      call chk_hy(lnkdiun(1),lnkdiun(idx_ozsrt))
      call chk_toctoc(lnkdiun(1),lnkdiun(idx_ozsrt))
   endif
	iopc= fstopc('MSGLVL','INFORMS',.false.)
	do i=1,niun
	   ier = fstfrm(lnkdiun(i))
	   call fclos(lnkdiun(i))
	enddo
!	call fstunl
	if (mode.eq.1) then
	   if (voirs)  then
	      if (message) then
	         ier = fstvoi(lnkdiun(idx_ozsrt), 'RND')
	      endif
	   endif
      ier = fstfrm(lnkdiun(idx_ozsrt))
	   call fclos(lnkdiun(idx_ozsrt))
	else
!
!    fermer fichier sequentiel
!
     	   if (mode.eq.3)  then
	      call fclos(lnkdiun(idx_ozsrt))
	   endif
     	   if (mode.eq.4)  then
	      call pgsmcf(lnkdiun(idx_ozsrt))
	   endif
	endif
!
!
!
!  fermer fichier 4 dans grille
!
!
!  imprime boite avec le temps d execution du pgm  pgsm
!
	if (ipose.gt.0) then
	   jdate= exfin('  PGSM  ', 'ABORT', 'NON')
	   call qqexit(13)
	else
	   jdate= exfin('  PGSM  ', 'OK', 'NON')
	endif
!
!	 stop
	 end
!**s/p pgsmabt  sortie pas trop brutale en cas d'erreur
!
      subroutine pgsmabt
      implicit none
      external abort,fclos,fstfrm,messags,exfin
      integer exfin,fstfrm
!
!auteur   p. sarrazin  rpn mars 1983
!
!langage  ratfor
!
!objet(pgsmabt)
!         faire une sortie sans dommage pour les fichiers
!         ouverts en cas d'erreur fatale dans pgsm
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
!
      integer ier
      jdate= exfin('  PGSM  ', 'ABORT' , 'NON')
!
      ier = fstfrm(lnkdiun(1))
      if (mode.eq.1)ier = fstfrm(lnkdiun(idx_ozsrt))
      call qqexit(13)
      end
!**S/P IPGSMLIC
!     
      integer function ipgsmlic(ifld,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar,ig1,ig2,ig3,ig4,grtyp)
      implicit none
      integer iun,ni,nj,nk,ip1,ip2,ip3,datev,ig1,ig2,ig3,ig4
      integer ifld(ni,nj,nk)
      character*12 etiket
      character*4 nomvar
      character*2 typvar
      character*1 grtyp
      external fstlic
      integer fstlic
      integer ier,i
      ier = fstlic(ifld,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar,ig1,ig2,ig3,ig4,grtyp)
      if (ier.lt.0) then
        ipgsmlic=ier
        return
      endif
      call prefiltre(ifld,ni,nj,nomvar,grtyp)
      ipgsmlic=ier
      return
      end
!**S/P PGSMLIC
!     
      integer function pgsmlic(fld,iun,ni,nj,nk,datev,etiket,      ip1,ip2,ip3,typvar,nomvar,ig1,ig2,ig3,ig4,grtyp)
      implicit none
      integer iun,ni,nj,nk,ip1,ip2,ip3,datev,ig1,ig2,ig3,ig4
      real fld(ni,nj,nk)
      character*12 etiket
      character*4 nomvar
      character*2 typvar
      character*1 grtyp
      external fstlic
      integer fstlic
      integer ier,i
      ier = fstlic(fld,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar,ig1,ig2,ig3,ig4,grtyp)
      if (ier.lt.0) then
        pgsmlic=ier
        return
      endif
      call prefiltre(fld,ni,nj,nomvar,grtyp)
      pgsmlic=ier
      return
      end
!**S/P IPGSMLIR
!     
      integer function ipgsmlir(ifld,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar,grtyp)
      implicit none
      integer iun,ni,nj,nk,ip1,ip2,ip3,datev,ig1,ig2,ig3,ig4
      integer ifld(ni,nj,nk)
      character*12 etiket
      character*4 nomvar
      character*2 typvar
      character*1 grtyp
      external fstlir
      integer fstlir
      integer ier,i
      ier = fstlir(ifld,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar)
      if (ier.lt.0) then
        ipgsmlir=ier
        return
      endif
      call prefiltre(ifld,ni,nj,nomvar,grtyp)
      ipgsmlir=ier
      return
      end
!**S/P PGSMLIR
!     
      integer function pgsmlir(fld,iun,ni,nj,nk,datev,etiket,      ip1,ip2,ip3,typvar,nomvar,grtyp)
      implicit none
      integer iun,ni,nj,nk,ip1,ip2,ip3,datev,ig1,ig2,ig3,ig4
      real fld(ni,nj,nk)
      character*12 etiket
      character*4 nomvar
      character*2 typvar
      character*1 grtyp
      external fstlir
      integer fstlir
      integer ier,i
      ier = fstlir(fld,iun,ni,nj,nk,datev,etiket,ip1,ip2,ip3,typvar,nomvar)
      if (ier.lt.0) then
        pgsmlir=ier
        return
      endif
      call prefiltre(fld,ni,nj,nomvar,grtyp)
      pgsmlir=ier
      return
      end
!**S/P IPGSMLUK
!
      integer function ipgsmluk(ifld,key,ni,nj,nk,nomvar,grtyp)
      implicit none
      integer key,ni,nj,nk
      integer, dimension(ni,nj) :: ifld
      character*4 nomvar
      character*1 grtyp
      external fstluk
      integer fstluk
      integer ier,i,j
      ier = fstluk(ifld,key,ni,nj,nk)
!          do j=1,nj
!            do i=1,ni
!              print *, i, j, ifld(i,j)
!            enddo
!         enddo
         if (ier.lt.0) then
        ipgsmluk=ier
        return
      endif
!      call prefiltre(fld,ni,nj,nomvar,grtyp)
      ipgsmluk=ier
      return
      end
!**S/P PGSMLUK
!
      integer function pgsmluk(fld,key,ni,nj,nk,nomvar,grtyp)
      implicit none
      integer key,ni,nj,nk
      real, dimension(ni,nj) :: fld
      character*4 nomvar
      character*1 grtyp
      external fstluk
      integer fstluk
      integer ier,i,j
      ier = fstluk(fld,key,ni,nj,nk)
!          do j=1,nj
!            do i=1,ni
!              print *, i, j, fld(i,j)
!            enddo
!         enddo
         if (ier.lt.0) then
        pgsmluk=ier
        return
      endif
      call prefiltre(fld,ni,nj,nomvar,grtyp)
      pgsmluk=ier
      return
      end
!
!**S/P  ADDITIONNE SOUSTRAIT MULTIPLIT MODULE 2 CHAMPS
!
      subroutine plmnmod(nom,type,idat,niv,ihr,ip3,etiqet)
      implicit none
!
      external fstinf,fstsui,fstprm,pgsmabt,imprime
      external lopascm,messags,memoir,pgsmluk,fstcvt
      integer fstinf,fstsui,fstprm,pgsmluk,fstcvt
!
!AUTEUR P. SARRAZIN AOUT 82 DRPN DORVAL P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(PLMNMOD)
!         LIRE UN CHAMP SUR FICHIER 1 OU 2 DE MEME NATURE ET DIMENSIONS
!         CELUI DANS L ACCUMULATEUR ET QUE L ON AJOUTE , SOUSTRAIT , MULTIPLIT
!         OU FAIT LA SOMME DE CHAQUE POINT DES DEUX CHAMPS AU CARRE.
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!ARGUMENTS
!   IN   NOM    NOM DU CHAMP  LCAR(GZ),"TT"....
!   IN   TYPE   TYPE DE CHAMP P=PREVISION  A=ANALYSE
!   IN   NIV    NIVEAU DU CHAMP
!   IN   IDAT   DATE DU CHAMP CMC STAMP
!   IN   IHR    HEURE DU CHAMP
!   IN   IP3    LIBRE (USAGER)
!   IN   ETIQET ETIQUETTE DU CHAMP 10 CARACTERES
!
!- - -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!MESSAGES
!         RECORD N EXISTE PAS SUR FICHIER
!         DIRECTIVE LIREE OU LIRES DOIT-ETRE APPELE AVANT
!                 'NI   ENTRE =',NIE,'NI LIREE\LIRES=',NNI,
!                 'NJ   ENTRE =',NJE,'NJ LIREE\LIRES=',NNJ,
!                 'NK   ENTRE =',NKE,'NK LIREE\LIRES=',NNK,
!                 'DIMENSION DU CHAMP   MAUVAISE'
!                 TYPE DE GRILLE DIFFERENT FATAL CHAMP=
!                 MAUVAISE HEMISPHERE CHAMP ...DOIT-ETRE=
!                 ERREUR  2 CHAMPS DIFFERENTS
!                 GRILLE INCONU DIRECTIVE  PLUS-MOIN-MODULE
!
!APPEL VIA DIRECTIVE
!       PLUSE(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET) FICHIER D'ENTRE
!       PLUSS(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET)  FICHIER DE SORTIE
!       MOINSE(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET) FICHIER D'ENTRE
!       MOINSS(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET)  FICHIER DE SORTIE
!       MODUL2E(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET) FICHIER D'ENTRE
!       MODUL2S(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET)  FICHIER DE SORTIE
!       FOISE(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET) FICHIER D'ENTRE
!       FOISS(NOM,TYPE,IDAT,NIV,IHR,IP3,ETIQET)  FICHIER DE SORTIE
!
!MODULES  FSTINF,FSTSUI,PGSMABT,FSTPRM,MEMOIR
!
!----------------------------------------------------------------------
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer blancs
      data blancs /4H    /
      integer ip1style, dateform
      common /styles/ ip1style, dateform
!
      character*12 cetike,cetiket
      character*4 cnomvar, cnumve
      character*2 ctypvar,ctypve
      character*1 cigtyp,cigtye
!
      integer type,etiqet(3),nom,aa,letiqet(3)
      integer idat,idate,ideete,if1,ig1e,ig2e,ig3e,ig4e
      integer ihr,ip3,irec,itot,iunit,jp1e,jp2e,jp3e,nie
      integer niv(2),nje,nke,npase
      integer cnbits,cdatyp,cswa,clng,cdltf,cubc,extra1,extra2,extra3
      integer argdims,letiket(3)
      external argdims
      integer lniv
      real p
      character*8 string
!
!    AA  MULTIPLICATEUR POUR AJOUTER OU SOUSTRAIRE
!
      aa=1
      iunit=1
!
!    VERIFIER SI DIRECTIVE LIREN OU LIRSR A ETE APPELE
!
 100  if (ichck.eq.0)  then
!     erreur faut appeler liren ou lirsr
         write(6,*)'DIRECTIVE LIREE OU LIRES DOIT-ETRE APPELE AVANT '
         call pgsmabt
      endif
!
!   MODIFICATION DE HOLLERITH A CARACTERE
!
      cnomvar = '    '
      ctypvar = '  '
      cetiket = '            '
      cigtyp  = ' '
      letiket(1) = etiqet(1)
      letiket(2) = blancs
      letiket(3) = blancs
      if (argdims(7).gt.1) then
         letiket(2) = etiqet(2)
      endif
      if (argdims(7).gt.2) then
         letiket(3) = etiqet(3)
      endif
      lniv = niv(1)
      if (argdims(4) > 1) then
         p = transfer(niv(1), p)
         call convip_plus(lniv, p, -1*niv(2)-1000, ip1style, string, .false.)
      endif
      ier = fstcvt(    nom,    type,  letiket,     -1,      cnomvar, ctypvar, cetiket, cigtyp,     .true.)
!
      irec=fstinf(iunit,nie,nje,nke,idat,cetiket,lniv,ihr,ip3,      ctypvar, cnomvar)
      if (irec.lt.0)  then
!     arret record n'EXISTE PAS
         write(6,*)         ' VERIFIER PLUSE/S - MOINSE/S - MODUL2E/S - FOIS(E\S)'
         write(6,*)' FSTINF A PAS RECONNU RECORD MAL DEFINI'
         call pgsmabt
      endif
 10   if (irec.gt.-1) then
!
!
         ier = fstprm(irec,idate,ideete,npase,nie,nje,nke,          cnbits,cdatyp,         jp1e,jp2e,jp3e,ctypve,cnumve,cetike,cigt&
     &ye,         ig1e,ig2e,ig3e,ig4e,         cswa, clng, cdltf, cubc, extra1, extra2, extra3)
         if (ier.lt.0) then
            write(6,*)' IER = FSTPRM NEGATIF VOIR PLUSE/MOIN.....'
         endif
!
!     VERIFIER SI GRILLE GAUSSIENNE NI DOIT ETRE PAIR
!
         if (cigtye.eq.'G'.and.mod(nie,2).ne.0)  call messags(nie)
!
!
         if (nie.ne.nni.or.nje.ne.nnj.or.nke.ne.nnk) then
            write(6,600)nie,nni
            write(6,610)nje,nnj
            write(6,620)nke,nnk
            write(6,*)'DIMENSION DU CHAMP MAUVAISE VERIFIER   '
            write(6,*)'PLUS(E\S) - MOINS(E\S) - MODUL2(E\S) - FOIS(E\S)'
            call pgsmabt
         endif
!
         if (cigty.ne.cigtye) then
            write(6,660)cigtye,cigty
            write(6,*)            ' VERIFIER PLUS(E\S)-MOINS(E\S)-MODUL2(E\S)-FOIS(E\S)'
            call pgsmabt
         endif
!
         if (cigty.eq.'G'.or.cigty.eq.'A'.or.cigty.eq.'B') then
            if (ig1e.ne.igg1) then
               write(6,*)'MAUVAISE HEMISPHERE CHAMP DE =',ig1e
               write(6,*)'DOIT-ETRE=',igg1
               write(6,*)             ' VERIFIER PLUS(E\S)-MOINS(E\S)-MODUL2(E\S)-FOIS(E\S)'
               call pgsmabt
            endif
         else
            if (cigty.eq.'N'.or.cigty.eq.'S'.or.cigty.eq.'L') then
               if (ig1e.ne.igg1.or.ig2e.ne.igg2.or.               ig3e.ne.igg3.or.ig4e.ne.igg4) then
                  write(6,*)'ERREUR  2 CHAMPS DIFFERENTS'
                  write(6,*)            ' VERIFIER PLUS(E\S)-MOINS(E\S)-MODUL2(E\S)-FOIS(E\S)'
                  call pgsmabt
               endif
            endif
         endif
!
!
!    ALLOCATION DE LA MEMOIRE
!
         allocate(tmpif1(nni,nnj))
!
         ier = pgsmluk(tmpif1, irec, nni, nnj, nnk,cnomvar,cigty)
         if (ier.lt.0)  then
            write(6,*)' IER=PGSMLUK(..... NEGATIF VOIR PLUS/MOIN.....'
            return
         endif
         if (printen)  call imprime(cnumve,tmpif1,nni,nnj)
!
!
!     AJOUTE 1 AU COMPTEUR ICNT DANS COMMON ACCUM INITIALISER
!     A 1 DANS MAIN PROGRAM
!
         icnt = icnt + 1
!
!
!     ADDITIONNE-SOUSTRAIT-MODULE-MULTIPLIT CHAQUE PTS DES DEUX CHAMPS
!
         itot=nni*nnj*nnk
         call lopascm(tmpif0,tmpif1,aa,itot)
!
!
         deallocate(tmpif1)
!
         if (aa.ne.1.or.unefois.or.once) goto 11
         irec = fstsui(iunit,nie,nje,nke)
!
         goto 10
      endif
 11   continue
!
      return
      entry pluss(nom,type,idat,niv,ihr,ip3,etiqet)
!
!   AA=MULTIPLICATEUR POUR AJOUTER
!
      aa=1
      iunit = 2
      go to 100
!
      entry moinse(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=MULTIPLICATEUR POUR  SOUSTRAIRE
!
      aa=-1
      iunit = 1
      go to 100
!
      entry moinss(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=MULTIPLICATEUR POUR  SOUSTRAIRE
!
      aa=-1
      iunit = 2
      go to 100
!
      entry modul2e(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     2   AA=2 ADDITIONNER LES DEUX CHAMPS AU CARRE
!
      aa=2
      iunit = 1
      go to 100
!
      entry modul2s(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=2 ADDITIONNER LES DEUX CHAMPS AU CARRE
!
      aa=2
      iunit = 2
      go to 100
!
      entry foise(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=3 MULTIPLIER CHAQUE PT DES DEUX CHAMPS
!
      aa=3
      iunit = 1
      go to 100
!
      entry foiss(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=3 MULTIPLIER CHAQUE PT DES DEUX CHAMPS
!
      aa=3
      iunit = 2
      go to 100
!
      entry divisee(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=4 MULTIPLIER CHAQUE PT DES DEUX CHAMPS
!
      aa=4
      iunit = 1
      go to 100
!
      entry divises(nom,type,idat,niv,ihr,ip3,etiqet)
!
!     AA=4 MULTIPLIER CHAQUE PT DES DEUX CHAMPS
!
      aa=4
      iunit = 2
      go to 100
!
 600  format(2x,'NI  ENTRE =',i10,'NI ACCUMULATEUR=',i10)
 610  format(2x,'NJ  ENTRE =',i10,'NJ ACCUMULATEUR=',i10)
 620  format(2x,'NK  ENTRE =',i10,'NK ACCUMULATEUR=',i10)
 660  format(2x,'MAUVAISE GRILLE ENTRE=',a1,'ACCUMULATEUR=',a1)
!
      end
      subroutine prefiltre(fld,ni,nj,nomvar,grtyp)
      implicit none
      integer ni,nj
      real fld(ni,nj)
      character*4 nomvar
      character*1 grtyp
      integer fltwgtlng(2),fltntimes(2),fltverb(2),fltlist(9,2)
      logical fltoggle(2)
      common /qfilter/ fltoggle,fltwgtlng,fltntimes,fltverb,fltlist
      if (fltoggle(1)) then
        if (grtyp.eq.'Y') then
          write (6, *)' (PREFILTRE) Impossible de filtrer des champs sur grille Y'
        else
          write (6, *) ' CHAMP FILTRE A LA LECTURE'
!          call statfld4 (fld,nomvar,0,'AVANFFLT',ni,nj,1,ni,nj,1,0,0,0)
          call filtre (fld, NI, NJ, 0, fltntimes(1), fltlist(1,1), fltwgtlng(1))
!          call statfld4 (fld,nomvar,1,'APRESFLT',ni,nj,1,ni,nj,1,0,0,0)
        endif
      endif
      return
      end
      subroutine putfld(fld, buf, iun,ibidon,iwrit,ni,nj,nbrow,npkc,istamp)
      implicit none
      integer iun, ibidon, iwrit, ni,nj,nbrow, npkc, istamp
      real fld(ni, nj)
      real buf(ni, nj)
      call putfld1(fld, ni*nj, iun)
      return
      end
      subroutine putfld1(fld, ni, iun)
      implicit none
      integer ni, iun
      real fld(ni)
      integer i
      write(iun) (fld(i), i=1,ni)
      return
      end
!     
!**S/P QAAQR   CALCUL TOURBILLON RELATIF
!
      subroutine qaaqr(qaqr, li, lj, xlat)
      implicit none
!
!AUTEUR  P. SARRAZIN DORVAL QUEBEC JUIN 83 DRPN
!
!LANGAGE RATFOR
!
!OBJET(QAAQR)
!            CALCUL LE CORIOLIS PARAMETER POUR CHAQUE POINT DE LA GRILLE
!            SOUSTRAIRE CE CHAMP DU CHAMP DU TOURBILLON ABSOLU ON GENERE
!            UN TOURBILONN RELATIF
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN OUT QAQR  CHAMP CONTENANT TOURBILLON ABSOLU 
!  IN     LI    NOMBRE DE POINTS SUR UNE RANGEE DU CHAMP QAQR
!  IN     LJ    NOMBRE DE POINTS DANS UNE COLONNE DU CHAMP QAQR
!  IN     XLAT  LATITUDE POUR CHAQUE POINT DU CHAMP QAQR
!
!APPEL
!         -VIA ROUTINE SCALAIR
!         CALL QAAQR(QAQR, LI, LJ, XLAT)
!
!MESSAGES 
!          -AUCUN
!
!IMPLICITES
!         - AUCUN
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
!
!
      integer i,j,li,lj
      real xlat(li,lj),qaqr(li,lj),degarad,omega2
!
!
!  rotation de la terre 7.292*1.e-5
      degarad = 3.1415926535/180.
      omega2= 2*7.292*1.e-5
!
!
      do j= 1,lj
         do i = 1,li
            qaqr(i,j) = qaqr(i,j) - omega2*sin(xlat(i,j)*degarad)
         enddo
      enddo
!
      return
      end
      subroutine qqqecho(chaine)
      implicit none
      external argdims
      integer argdims
      integer chaine(20)
      integer i,j,longueur,iun
      character*80 message
      character*16 form
      longueur = argdims(1)
      message(1:80) = ' '
      do i=1,longueur
         j = 4*(i-1)+1
         write(message(j:j+3),'(a4)') chaine(i)
      enddo
      iun = 2
      call pgsmecho(iun, message,longueur*4)
      return
      end
      subroutine qqqfilt(inout,weightlst,ntimes,verbose)
      implicit none
      integer inout,weightlst(*),ntimes,verbose
      integer fltwgtlng(2),fltntimes(2),fltverb(2),fltlist(9,2)
      logical fltoggle(2)
      common /qfilter/ fltoggle,fltwgtlng,fltntimes,fltverb,fltlist
      integer i,j,l,sum
      integer nb_elem,lng_liste,istart, iend
      external argdims
      integer argdims
      data(fltoggle(i),i= 1,2)  /.false.,.false./
      data(fltntimes(i), i=1,2) /0,0/
      data(fltverb(i), i=1,2)   /0,0/
      if (inout.ne.1.and.inout.ne.2) then
        write (6,*) '(QQQFILT) On doit mentionner LECTURE ou ECRITURE'
        return
      endif
      if (argdims(2).gt.9) then
        write (6,*) '(QQQFILT) Liste de poids trop longue - Maximum=9'
        return
      endif
      fltoggle(inout)=.true.
      if (argdims(2).le.1) then
         fltoggle(inout)=.false.
         return
      endif
      fltwgtlng(inout)=argdims(2)
      fltntimes(inout)=ntimes
      fltverb(inout)=verbose
      if (mod(fltwgtlng(inout),2).ne.1) then
        write (6,*) '(QQQFILT) Liste de poids doit etre impaire'
        return
      endif
      do i=1, fltwgtlng(inout)
        fltlist(i,inout)=weightlst(i)
      enddo
      return
      end
      subroutine qqqform(theform)
      implicit none
      integer theform(4)
   integer qposition, qitems(16), qnitems
   character*1 qcsepar
   character*16 qcform
   common /idents/  qposition, qitems, qnitems
   common /cidents/ qcsepar,qcform
      external argdims
      integer argdims
      integer i,j,longueur
      character*80 chaine
      longueur = argdims(1)
      qcform(1:16) = ' '
      do i=1,longueur
         j = 4*(i-1)+1
         write(qcform(j:j+3),'(a4)') theform(i)
      enddo
      return
      end
      subroutine qqqident(position, separateur, items)
      implicit none
      integer position, separateur, items(16)
      integer argdims
      external argdims
      integer champpr,nchamp,nchmp,npar
      common / champs/ nchamp, nchmp, npar,champpr(31)
!     
!
!
   integer qposition, qitems(16), qnitems
   character*1 qcsepar
   character*16 qcform
   common /idents/  qposition, qitems, qnitems
   common /cidents/ qcsepar,qcform
      integer i
!      print *, 'qqqident', position, separateur
!      do i=1,16
!        print *, i, items(i)
!      enddo
      qposition = position
      if (qposition.eq.5) then
         qnitems = 0
      endif
      if (separateur.eq.-1) then
         qcsepar = 'T'
      else
         write(qcsepar,'(a1)') separateur
      endif
      do i=1,16
         qitems(i) = 0
      enddo
      if (-1.eq.items(1)) then
         qnitems = 11
         do i=1,qnitems
            qitems(i) = i
         enddo
      else
         do i=1,argdims(3)
            qitems(i) = items(i)
         enddo
      endif
      return
      end
!     ***************************************************************
!     *                     A S S E M B L E                         *
!     * Object :                                                    *
!     *         To assemble data field                              *
!     *                                                             *
!     * Arguments :                                                 *
!     *            IN     ni    : 1st dimension of field ZOUT       *
!     *            IN     nj    : 2nd dimension of field ZOUT       *
!     *            IN     nrows : 3rd dimension of field ZOUT       *
!     *            IN     slab  : data to assemble                  *
!     *            IN     nX    : dimension of hxpos                *
!     *            IN     hxpos : indicator of position in the grid *
!     *                                                             *
!     *            OUT    ZOUT  : field to return (assembled)       *
!     *                                                             *
!     ***************************************************************
      subroutine assemble(ZOUT,ni,nj,nrows,slab,nX,hxpos)
      implicit none
      integer nj, ni, nX, nrows
      real ZOUT(ni * nj, nrows)
      integer hxpos(nX)
      real slab(nX,nrows)
      integer I,k
      do k=1, nrows
         do I=1, nX
            ZOUT(hxpos(I), k) = slab(I,k)
         enddo
      enddo
      return
      end
!     ***************************************************************
!     *                     W R T S T D F                           *
!     * Object :                                                    *
!     *         To write standard file (FSTD)                       *
!     *                                                             *
!     * Arguments :                                                 *
!     *            IN    ZOUT   : data field to read                *
!     *            IN    iun    : unit number of the file           *
!     *            IN    dateo  : origin date of the field          *
!     *            IN    deet   : time step length in seconds       *
!     *            IN    npas   : time step number                  *
!     *            IN    ni     : 1st dimension of field            *
!     *            IN    nj     : 2nd dimension of field            *
!     *            IN    nrows  : 3rd dimension of field            *
!     *            IN    ip1    : descriptor 1 (1 to 32767)         *
!     *            IN    ip2    : descriptor 2 (1 to 32767)         *
!     *            IN    ip3    : descriptor 3 (1 to 32767)         *
!     *            IN    typvar : field type                        *
!     *            IN    nomvar : name of field                     *
!     *            IN    etiket : 9 caracter stamp                  *
!     *            IN    grtyp  : grid type                         *
!     *            IN    ig1    : grid descriptor 1 (0 to 2047)     *
!     *            IN    ig2    : grid descriptor 2 (0 to 2047)     *
!     *            IN    ig3    : grid descriptor 3 (0 to 65535)    *
!     *            IN    ig4    : grid descriptor 4 (0 to 65535)    *
!     *            IN    datyp  : type of data field                *
!     *            IN    Nextra : number of extra parameters        *
!     *                           (Nextra >= 0)                     *
!     *            IN    xtra   : field of optionnal variable       *
!     *                                (absent IF Nextra = 0)       *
!     *                                                             *
!     ***************************************************************
!       subroutine wrtstdf (ZOUT,iun, dateo, deet, npas, ni, nj,nxgrid,            nygrid,  nrows, ip1, ip2, ip3, typvar,nomvar, etiket,            grtyp, ig1, ig2, ig3, ig4, datyp,Nextra,xtra, nbits,            iflt,list,L, S)
!       implicit none
!       integer fstecr
!       integer ni, nj, nrows, k, Nextra, i, j
!       integer nxgrid, nygrid
!       real ZOUT(ni , nj, nrows), work(1), xtra(nrows, Nextra)
!       integer ip1(nrows), ip2(nrows), ip3(nrows),npak, nbits(nrows)
!       integer ig1, ig2, ig3, ig4
!       integer iun, datyp(nrows)
!       integer npas, deet, dateo
!       character *4 nomvar(nrows)
!       character *4 typvar(nrows)
!       character *4 grtyp
!       character *12 etiket
!       integer ierr, ier, S
!       integer L, iflt(nrows)
!       integer list(L)
!       integer sum
!       real, dimension (:,:), allocatable :: fact, temp
!
!       allocate(fact(L, (L+1)/2))
!       do k=1, (L+1)/2
!          do I=1, L
!             fact(I,k) = 0
!          enddo
!       enddo
!       do k=1,(L + 1)/2
!          sum = 0
!          do I=k, L - k + 1
!             sum = sum + list(I)
!          enddo
!          do I=k,L - k + 1
!             FACT(I,k) = float(list(I)) / float(sum)
!          enddo
!       enddo
!
!       allocate(temp(nxgrid, nj))
!       do 300 k=1, nrows
!          npak = -nbits(k)
!          if (nxgrid .eq. (ni + 1)) then
!             do j=1, nj
!                do i=1, ni
!                   Temp(i,j) = ZOUT(i,j,k)
!                enddo
!             enddo
!             do j=1,nj
!                Temp(nxgrid,j) = Temp(1,j)
!             enddo
!             if ((iflt(k) .GT. 0) .and. (L .gt. 1)) then
!                call filtre (Temp,nxgrid,nj,nrows,iflt(k),FACT,list,L)
!             endif
!             ierr = fstecr(Temp, work, npak, iun, dateo, deet,             npas, nxgrid, nj, 1, ip1(k), ip2(k), ip3(k),             typvar(k)(1:1), nomvar(k)(1:2), etiket(1:8),            grtyp(1:1), ig1, ig2, ig3, ig4, datyp(k),            .false.)
!
!          else
!             if ((iflt(k) .GT. 0) .and. (L .gt. 1)) then
!                call filtre (ZOUT(1,1,k),nxgrid,nj,nrows,iflt(k),FACT               ,list,L)
!             endif
!             ierr = fstecr(ZOUT(1,1,k), work, npak, iun, dateo, deet,             npas, ni, nj, 1, ip1(k), ip2(k), ip3(k),             typvar(k)(1:1), nomvar(k)(1:2), etiket(1:8),            grtyp(1:1), ig1, ig2, ig3, ig4, datyp(k),            .false.)
!          endif
!  300  continue
!
!       deallocate(temp)
!       deallocate(fact)
!
!       if (Nextra .ne. 0) then
!          ierr = fstecr(xtra ,WORK,npak, iun, 20002020,         1, 1, nrows, Nextra,1, 0,0,S,'|',         '||' ,'||||*||||','x',0,0,         0,0,1, .false.)
!       endif
!
!       return
!       end
!     ***************************************************************
!     *                       W S T D F X Y                         *
!     * Object :                                                    *
!     *         To write record ('>>' and '^^') in standard file    *
!     *                                                             *
!     * Arguments :                                                 *
!     *            IN    xpos   : field to write (dim : ni)         *
!     *            IN    ypos   : filed to write (dim : nj)         *
!     *            IN    iun    : unit number of the file           *
!     *            IN    datoe  : date of origine of the field      *
!     *            IN    deet   : time step lenght in seconds       *
!     *            IN    npas   : time step number                  *
!     *            IN    ni     : dimension of xpos                 *
!     *            IN    nj     : dimension of ypos                 *
!     *            IN    ip1    : descriptor 1                      *
!     *            IN    ip2    : descriptor 2                      *
!     *            IN    ip3    : descriptor 3                      *
!     *            IN    etiket : 9 caracter stamp                  *
!     *            IN    grtyp_ : grid type for ">>" and "^^"       *
!     *            IN    ig1_   : grid descriptor 1 of ">>" and "^^"*
!     *            IN    ig2_   : grid descriptor 2 of ">>" and "^^"*
!     *            IN    ig3_   : grid descriptor 3 of ">>" and "^^"*
!     *            IN    ig4_   : grid descriptor 4 of ">>" and "^^"*
!     *                                                             *
!     ***************************************************************
      subroutine wstdfxy (xpos, ypos, iun, dateo, deet, npas, ni, nj,  	                  ip1, ip2, ip3, etiket, grtyp_, ig1_,     &
     &                ig2_, ig3_, ig4_)
      implicit none
      integer fstecr
      integer ni, nj
      real xpos(ni), ypos(nj), work(1)
      integer ip1, ip2, ip3
      integer ig1_, ig2_, ig3_, ig4_
      integer datyp, npak, npas, deet, dateo
      integer i
      character *4 grtyp_
      character *12 etiket
      integer ierr, iun
      npak = -24
      datyp = 1
      ierr = fstecr(xpos, work, npak, iun, dateo, deet,                   npas, ni, 1, 1, ip1, ip2, ip3,                   'X ', '>&
     &>  ', etiket,                  grtyp_(1:1), ig1_, ig2_, ig3_, ig4_, datyp,                  .false.)
      ierr = fstecr(ypos, work, npak, iun, dateo, deet,                   npas, 1, nj, 1, ip1, ip2, ip3,                   'X ', '^&
     &^  ', etiket,                  grtyp_(1:1), ig1_, ig2_, ig3_, ig4_, datyp,                  .false.)
      return
      end
!**s/p scalair  interpolation horizontale d un champ
!               defini par l usager
!
      subroutine scalair(cnom, iheur, nniv, niveaux)
      implicit none
      external ecritur, pgsmluk, ipgsmluk, fstinf, fstsui, memoir, fstprm, qaaqr, fstcvt,     fstsel, symetri, imprime, itrouve, me&
     &ssags, pgsmabt
      external cvtrfi,cvtifr
      external liraxez
      integer  pgsmluk, ipgsmluk, fstinf, fstsui, fstprm, fstcvt, fstsel, fstinl
!      integer diesinf, dieslir, diesisincache, res
      integer ezgdef_fmem, ezqkdef, ezsint, ezdefset, ezgetopt
      logical skip
!
!auteur  p. sarrazin  dorval quebec fevrier 82 drpn
!revision 4.0.2
!   conversion des variables hollerith en caractere
!   y. chartier -aout 90- drpn dorval quebec
!revision 5.0.0
!   utilisation de la cuvee 91 des interpolateurs
!   elimination des appels a 'fstcvt' pour convertir etikent
!   y. chartier -mai 1991- drpn dorval quebec
!
!revision 5.2.1
!   utilisation de fstinl au lieu de fstinf-fstsui
!   l'utilisation des deux dernieres avec les grilles Z
!   causait un probleme lorsqu'on lisait les axes
!   y. chartier -oct. 92- drpn dorval quebec
!langage ratfor
!
!objet(scalair)
!              interpolation horizontale des scalaires gz, tt, dd, ww, qq, es.
!              d une grille a une autre pour nniv niveaux
!              ecrire resultat sur fichier standard, ms, sequentiel
!
!librairies
!         -source  armnsrc, drpn
!         -objet   pgsmlib, id=armnpjs.
!
!arguments
!  in   nom    nom du champ 2 caracteres gz.tt.dd.......
!  in   iheur  heure de la variable sur fichier d entre
!  in   nniv   nombre de niveaux
!  in   niveaux table contenant nniv niveaux
!
!appel
!         -via routine champ
!         call scalair(nom, iheur, nniv, niveaux)
!
!messages
!         mauvaise directive champ (scalair)
!         aucun niveau specifie dans directive champ
!         record n'existe pas sur fichier d'entre (scalair)
!         aucune interpolation horizontale champ
!
!
!modules pgsmabt, rfl, fstinf, fstprm, pgsmlir, rgscint, ecritur
!
!-----------------------------------------------------------------
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   character *12 cetiket
   character *4 cnom, cnomvar, cnomx
   character *2 ctypvar
   character *1 cigtyp
   character*8 extrap_option
   integer nniv, dateo, datev, i, nunv, itrouve, key, ii, j
   integer niveaux(nniv), deet, ig1, ig2, ig3, ig4, iheur
   integer iprs, irec, iunit, jp1, jp2, jp3, ne, ni, nj, nk, total_keys, nrecs
   integer cnbits, cdatyp, cswa, clng, cdltf, cubc, extra1, extra2, extra3
   logical sym, symetri
   real fbidon
   real, dimension(:), allocatable :: ax, ay
   integer, dimension(:,:), allocatable :: ifld_in, ifld_out
   real, dimension(:,:), pointer :: tmpout
!
   integer, dimension(:), allocatable :: liste
   logical, dimension(:), allocatable :: done_fields
   nunv=0
   iunit=lnkdiun(1)
   call pgsm_get_nfstkeys(total_keys)
!   print *, 'scalair: total_keys', total_keys
   allocate(liste(total_keys))
   allocate(done_fields(total_keys))
   done_fields = .false.
   cnomx = cnom
   if (cnom.eq.cnomqr) cnom=cnomqq
!
!
   do iprs = 1, nniv
!
!  conversion de l etiquette d'entree en caracteres
!
      if (etikent(1)  /=  -1) then
         write(cetiket, '(3A4)') (etikent(i), i=1, nwetike)
      else
         cetiket = '            '
      endif
      if (typeent  /=  -1) then
         write(ctypvar, '(A2)') typeent
      else
         ctypvar = '  '
      endif
!     identifier numero du record
!
      if (cnom.eq.cnommt) then
         if (.not.mtdone) then
            irec = fstinl(iunit, ni, nj, nk, -1, '            ', -1, -1, -1, '  ', cnom,           liste, nrecs, total_keys)
            if (nrecs .eq. 0) then
               if (message) then
                  write(6, *)                'RECORD FICHIER DE MONTAGNE PAS LA (SCALAIR)'
               endif
            go to 5000
            endif
            if (nk.gt.1) then
               write(6, *)'******************************************'
               write(6, *)'       PGSM N ACCEPTE PAS UN          '
               write(6, *)' CHAMP DE 3 DIMENSIONS NK>1 ??               (MT DANS SCALAIR)'
               write(6, *)'*****************************************'
               call pgsmabt
            endif
 5000       mtdone=.true.
         endif
      else
         call chk_userdate(datev)
         irec=fstinl(iunit, ni, nj, nk, datev, cetiket, niveaux(iprs), iheur, ip3ent, ctypvar, cnom, liste, nrecs, total_keys)
         if (nrecs .eq. 0) then
            if (message) then
               write(6, *)              'RECORD N EXISTE PAS SUR FICHIER D ENTRE (SCALAIR)'
            endif
         goto 22000
         endif
         if (nk.gt.1) then
            write(6, *)            '************************************************'
            write(6, *)            '         PGSM N ACCEPTE PAS UN          '
            write(6, *)            ' CHAMP DE 3 DIMENSIONS NK>1 ?? (CHAMP SCALAIR)'
            write(6, *)            '************************************************'
            call pgsmabt
         endif
      endif
!
!
      do ii=1, nrecs
!
         irec = liste(ii)
!  On verifie d'abord que le champ n'a pas ete traite, car les masques associes sont traites ailleurs.
         if (done_fields(ii)) cycle
         ier=fstprm(irec, dateo, deet, npas, ni, nj, nk, cnbits, cdatyp, jp1, jp2, jp3, ctypvar, cnomvar, cetiket, cigtyp, &
            ig1, ig2, ig3, ig4, cswa, clng, cdltf, cubc, extra1, extra2, extra3)
         npack_orig = -cnbits
!         print *, 'npack_orig=', npack_orig, 'npack=', npack
         if (ier .lt. 0) write(6, *)          ' IER = FSTPRM NEGATIF VOIR SCALAIR'
         if (cnomvar(1:2).eq.'>>'.or.cnomvar(1:2).eq.'^^'.or.cnomvar(1:2).eq.'HY'.or.cnomvar.eq.'^>'.or.cnomvar.eq.'!!') then
            cycle
         endif
         if (ctypvar(2:2).eq.'@'.and.masque == 1) then
            if (ctypvar(1:1) /= '@') then
              call scalair_msk(irec, liste, done_fields,total_keys)
            endif
            cycle
         endif
	 extrap_option = '        '
	 ier = ezgetopt('extrap_degree', extrap_option)
	 if (extrap_option(1:5) == 'abort') then
            gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
            call chk_extrap(gdout, gdin, li, lj, ni, nj)
	 endif
!
!   si la directive de champ emploi tout  champ(tout, tout) on doit faire
!     attention pour les vecteurs u-v car l'interpolation serait scalaire
!     on verifie et les interpolations pour les vecteurs ne sont pas faites
!
         do i=1, npair
            if (cnom.eq.'    ') then
               if ((paire(i)(9:10).eq.cnomvar(1:2)).or.(paire(i)(13:14).eq.cnomvar(1:2)))  then
                  nunv=nunv + 1
                  goto 99999
               endif
            endif
         enddo
!
!     lire champ
         skip = .false.
	 allocate(tmpif1(ni,nj))
!
         if (cdatyp .eq. 2 .or. cdatyp .eq. 4.or.cdatyp.eq.130.or.cdatyp.eq.132) then
            allocate(ifld_in(ni, nj))
            allocate(ifld_out(li, lj))
	    key = ipgsmluk(ifld_in, irec, ni, nj, nk, cnomvar, cigtyp)
!             print *, 'CDATYP=', cdatyp
            call cvtrfi(tmpif1, ifld_in, ni, nj)
         else
	    key = pgsmluk(tmpif1, irec, ni, nj, nk, cnomvar, cigtyp)
!            print *, 'CDATYP=', cdatyp
         endif
         if (printen)  call imprime(cnomvar, tmpif1, ni, nj)
         if (cigtyp == 'A' .or. cigtyp == 'B' .or. cigtyp == 'G') then
            if (ig1 /= 0) sym = symetri(cnomvar)
         endif
!
!     on ne fait pas d'interpolation si igtyp=grtyp  ig1=lg1  ig2=lg2
!     ig3=lg3  ig4=lg4
!
         if (.not.skip) then
            if (cigtyp /= cgrtyp.or. &
               ig1 /= lg1.or.ig2 /= lg2.or.ig3 /= lg3.or.ig4 /= lg4.or. &
               li /= ni.or.lj /= nj) then
!
!     interpolation
!
               allocate(tmpif2(li,lj))
               gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
               ier = ezdefset(gdout, gdin)
               ier = ezsint(tmpif2, tmpif1)
               tmpout => tmpif2
            else
               tmpout => tmpif1
               if (message) write(6, 660) cnomvar
 660  format(2x, 'AUCUNE INTERPOLATION HORIZONTALE CHAMP=', a4)
            endif
!
!
!     ecrire sur fichier approprie(std, ms, seq)
!
            if (cnomx.eq.cnomqr) then
               call qaaqr(tmpif2, li, lj, tmplat)
               cnomvar=cnomqr
            endif
!            print *, 'Avant iecritur : cdatyp', cdatyp
         if (cdatyp .eq. 2 .or. cdatyp .eq. 4.or.cdatyp.eq.130.or.cdatyp.eq.132) then
!            print *, 'Avant iecritur : cdatyp', cdatyp
            call cvtifr(ifld_out,tmpout,li,lj)
            call iecritur(ifld_out, npack, dateo, deet, npas, li, lj, nk, jp1, jp2, jp3, ctypvar, cnomvar, cetiket, cgrtyp, lg1, lg&
     &2, lg3, lg4)
            deallocate(ifld_in, ifld_out)
         else
            call ecritur(tmpout, npack, dateo, deet, npas, li, lj, nk, jp1, jp2, jp3, ctypvar, cnomvar, cetiket, cgrtyp, lg1, lg2, &
     &lg3, lg4)
            if (associated(tmpif2)) then
               deallocate(tmpif2)
            endif
            deallocate(tmpif1)
         endif
         endif
99999    continue
!
         if (unefois) goto 23000
!
         enddo
22000 enddo
23000 continue
!
      if (nunv.gt.0) then
         write(6, 666)
         write(6, 668)
         write(6, 669)
      endif
   deallocate(liste)
 666  format(' AUCUNE INTERPOLATION SUR VARIABLE PAIRE CHAMP(TOUT, TOUT)')
 668  format(' ON DOIT UTILISER LE NOM DE LA VARIABLE EX: CHAMP(UU, TOUT)')
 669  format(' ATTENTION L INTERPOLATION DES VECTEURS SERA SCALAIRE (!!!)')
   return
   end
   subroutine cvtrfi(rfld, ifld, ni, nj)
   implicit none
   integer ni, nj
   real rfld(ni, nj)
   integer ifld(ni, nj)
   integer i,j
   do j=1, nj
      do i=1, ni
         rfld(i, j) = real(ifld(i, j))
      enddo
   enddo
   return
   end
   subroutine cvtifr(ifld, rfld, ni, nj)
   implicit none
   integer ni, nj
   integer ifld(ni, nj)
   real rfld(ni, nj)
   integer i,j
   do j=1, nj
      do i=1, ni
         ifld(i, j) = nint(rfld(i, j))
      enddo
   enddo
   return
   end
!**s/p scalair  interpolation horizontale d un champ
!               defini par l usager
      subroutine scalair_msk(key, liste, done_liste, len_liste)
      implicit none
      integer :: key,len_liste
      integer, dimension(len_liste) :: liste
      logical, dimension(len_liste) :: done_liste
      external ecritur,pgsmluk,fstinf,fstsui,memoir,fstprm,qaaqr,fstcvt, &
         fstsel,symetri,imprime,itrouve,messags,pgsmabt
      external cvtifr
      external liraxez
      integer  pgsmluk, fstinf, fstsui, fstprm, fstcvt, fstsel, fstinl, fstluk
      integer ezgdef_fmem, ezqkdef, ezsint, ezsint_mdm, ezdefset, fst_get_mask_key, key_masq
      logical skip
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   character *12 cetiket
   character *4 cnom,cnomvar,cnomx
   character *2 ctypvar
   character *1 cigtyp
   integer nniv, i, nunv, itrouve, ii
   integer niveaux(512), deet, ig1, ig2, ig3, ig4, iheur
   integer dateo, datev, nbits, datyp, ip1, ip2, ip3
   integer iprs, irec, iunit, ne, ni, nj, nk, total_keys, nrecs
   integer cnbits, cdatyp, cswa, clng, cdltf, cubc, extra1, extra2, extra3
   logical, save :: unefoys = .true.
   integer dateo_masq, deet_masq, npas_masq, ni_masq, nj_masq, nk_masq
   integer nbits_masq, datyp_masq, ip1_masq, ip2_masq, ip3_masq
   integer ig1_masq, ig2_masq, ig3_masq, ig4_masq, cswa_masq, clng_masq, cdltf_masq, cubc_masq
   integer datev_masq, extra2_masq, extra3_masq
   logical sym, symetri
   character(len=4)  :: cnomvar_masq
   character(len=2)  :: ctypvar_masq
   character(len=12) :: cetiket_masq
   character(len=1)  :: cigtyp_masq
   real fbidon
   real, dimension(:,:), allocatable, target ::  fld, fld_out
   integer, dimension(:,:), allocatable, target :: masq, masq_out, masq_zones
   real, dimension(:,:), pointer :: tmpout
   integer, dimension(:,:), pointer :: tmpmsk
   logical masque_present, masque_done
   nunv=0
   iunit=1
   masque_present = .false.
   masque_done = .false.
   ctypvar_masq = '@@'
   ier = fstprm(key, dateo, deet, npas, ni, nj, nk, nbits, datyp, ip1, ip2, ip3, ctypvar, cnomvar, cetiket, &
            cigtyp, ig1, ig2, ig3, ig4, cswa, clng, cdltf, cubc, datev, extra2, extra3)
   allocate(fld(ni,nj), fld_out(li,lj))
   ier = fst_get_mask_key(key_masq, key, 0, iunit)
   if (key_masq >= 0) then
      masque_present = .true.
      allocate(masq(ni,nj), masq_out(li,lj))
   else
      masque_present = .false.
   endif
   key = pgsmluk(fld, key, ni,nj,nk,cnomvar,cigtyp)
   if (cdatyp  ==  2 .or. cdatyp  ==  4) then
      call cvtifr(fld, fld, ni, nj)
   endif
   if (printen)  call imprime(cnomvar,tmpif1,ni,nj)
   if (ig1 /= 0) sym = symetri(cnomvar)
!
!  on ne fait pas d'interpolation si igtyp=grtyp  ig1=lg1  ig2=lg2
!  ig3=lg3  ig4=lg4
!
   if (.not.masque_present) then
      if (cigtyp /= cgrtyp.or.cigtyp == 'Z'.or.ig1 /= lg1.or.ig2 /= lg2.or.ig3 /= lg3.or.ig4 /= lg4.or.li /= ni.or.lj /= nj) then
         gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
         ier = ezdefset(gdout, gdin)
         ier = ezsint(fld_out, fld)
         tmpout => fld_out
      else
         tmpout => fld
         if (message) write(6,660) cnom
 660           format(2x,'AUCUNE INTERPOLATION HORIZONTALE CHAMP=',a2)
      endif
!     ecrire sur fichier approprie(std,ms,seq)
      if (cnomx == cnomqr) then
         call qaaqr(tmpif2,li,lj,tmplat)
         cnomvar=cnomqr
      endif
      call ecritur(tmpout,npack,dateo,deet,npas,li,lj,nk,ip1,ip2, ip3, ctypvar, cnomvar, cetiket, cgrtyp, lg1, lg2, lg3, lg4)
      deallocate(fld_out)
   else
      if (cigtyp /= cgrtyp.or.cigtyp == 'Z'.or.ig1 /= lg1.or.ig2 /= lg2.or.ig3 /= lg3.or.ig4 /= lg4.or.li /= ni.or.lj /= nj) then
         gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
         ier = ezdefset(gdout, gdin)
         ier = fstprm(key_masq, dateo_masq, deet_masq, npas_masq, ni_masq, nj_masq, nk_masq, &
                  nbits_masq, datyp_masq, ip1_masq, ip2_masq, ip3_masq, &
                  ctypvar_masq, cnomvar_masq, cetiket_masq, &
                  cigtyp_masq, ig1_masq, ig2_masq, ig3_masq, ig4_masq, cswa_masq, clng_masq, &
                  cdltf_masq, cubc_masq, datev_masq, extra2_masq, extra3_masq)
         ier = fstluk(masq, key_masq, ni,nj,nk)
         ier = ezsint_mdm(fld_out, masq_out, fld, masq)
         tmpout => fld_out
         tmpmsk => masq_out
      else
         tmpout => fld
         tmpmsk => masq
         if (message) write(6,662) cnom
 662           format(2x,'AUCUNE INTERPOLATION HORIZONTALE CHAMP=',a2)
      endif
!     ecrire sur fichier approprie(std,ms,seq)
      if (cnomx == cnomqr) then
         call qaaqr(tmpif2,li,lj,tmplat)
         cnomvar=cnomqr
      endif
      call ecritur(tmpout,npack,dateo,deet,npas,li,lj,nk, ip1, ip2, ip3, &
         ctypvar, cnomvar, cetiket, cgrtyp, lg1, lg2, lg3, lg4)
         masque_done = .false.
         do i=1,len_liste
            if (liste(i) == key_masq) then
               if (done_liste(i)) then
                  masque_done = .true.
               endif
               exit
            endif
         enddo
      if (masque_present.and..not.masque_done) then
         if (unefoys) then
            allocate(masq_zones(li,lj))
            call ezget_mask_zones(masq_zones, masq)
            nk = 1
            call iecritur(masq_zones,-16,dateo_masq,deet_masq,npas_masq,li,lj, nk,&
               ip1_masq,ip2_masq, ip3_masq, '@Z', cnomvar_masq, cetiket_masq, &
               cgrtyp, lg1, lg2, lg3, lg4)
            unefoys = .false.
            deallocate(masq_zones)
         endif
         call iecritur(tmpmsk,-nbits_masq,dateo_masq,deet_masq,npas_masq,li,lj, nk,&
            ip1_masq,ip2_masq, ip3_masq, ctypvar_masq, cnomvar_masq, cetiket_masq, &
            cgrtyp, lg1, lg2, lg3, lg4)
         do i=1,len_liste
            if (liste(i) == key_masq) then
               done_liste(i) = .true.
               exit
            endif
         enddo
         if (allocated(masq_out)) then
            deallocate(masq_out)
         endif
      endif
      deallocate(fld_out)
   endif
   deallocate(fld)
!
   if (nunv > 0) then
      write(6,666)
 666     format(' AUCUNE INTERPOLATION SUR VARIABLE PAIRE CHAMP(TOUT,TOUT)')
      write(6,668)
 668     format(' ON DOIT UTILISER LE NOM DE LA VARIABLE EX: CHAMP(UU,TOUT)')
      write(6,669)
 669  format(' ATTENTION L INTERPOLATION DES VECTEURS SERA SCALAIRE (!!!)')
!
   endif
   return
   end
   subroutine qqqintx(ordre)
   implicit none
   integer ordre
!
!auteur  Y.Chartier Dec 91
!    Utilise le nouveau dispatcher de fscint.f pour ajuster
!    le degre d'interpolation
!
!
!  e      ordre     ordre de l'interpolation 0,1, ou 3
!
!implicites
!*
character*8 op
integer ezsetopt
integer ier
   select case (ordre)
   case (100)
      ier = ezsetopt('interp_degree', 'nearest')
   case (1)
      ier = ezsetopt('interp_degree', 'linear')
   case (3)
      ier = ezsetopt('interp_degree', 'cubic')
   case (4)
      ier = ezsetopt('interp_degree', 'average')
   case (5)
      ier = ezsetopt('interp_degree', 'sph_average')
   case default
      print *,  ' <QQQINTX>: MAUVAISE VALEUR - VALEUR DEVRAIT ETRE 0, 1 OU 3'
      print *,  ' <QQQINTX>: ORDRE D''INTERPOLATION INITIALISE A 3'
      ier = ezsetopt('interp_degree','cubic')
   end select
   return
   end
      subroutine setxtrap(val)
      implicit none
      integer val
      integer voisin, abort, valeur, maximum, minimum
      parameter (voisin  =   100)
      parameter (maximum =   4)
      parameter (minimum =   5)
      parameter (valeur  =   6)
      parameter (abort   =  13)
      integer n, i, j
      integer ival
      real    rval
      integer ier, ezsetval, ezsetopt
      character*8 op,v
      equivalence (ival, rval)
      ival = val
      op = 'EXTRAP'
      if (val .ne. voisin .and. val .ne. minimum .and. val .ne. maximum       .and. val .ne. abort .and. val .ne. 1) then
         v = 'VALEUR'
         ier = ezsetval('extrap_value', rval)
         ier = ezsetopt('extrap_degree', 'value')
      else
         if (val .eq. 100) then
            ier = ezsetopt('extrap_degree', 'NEAREST')
            v = '0'
         else if (val .eq. 1) then
            ier =  ezsetopt('extrap_degree', 'LINEAR')
            v = '1'
         else if (val .eq. 3)  then
            ier =  ezsetopt('extrap_degree', 'CUBIC')
            v = '3'
         else if (val .eq. minimum) then
            v = 'MINIMUM'
            ier =  ezsetopt('extrap_degree', v)
         else if (val .eq. maximum) then
            v = 'MAXIMUM'
            ier =  ezsetopt('extrap_degree', v)
         else
            v = 'ABORT'
            ier =  ezsetopt('extrap_degree', v)
         endif
      endif
      return
      end
!
!**S/P SORTI   IDENTIFICATION DU FICHIER..OUVRIR FICHIER...RESERVER MEMOIRE
!
      subroutine sorti( modx, norecs, jwrit)
      implicit none
      external fnom,fstnbr,fsteof,fstouv,pgsmabt,fstvoi,fstapp,messags
      external fstabt, fstlnk, exfin, pgsmof
      integer  fnom,fstnbr,fsteof,fstouv,fstvoi,fstapp,exfin,pgsmof
!
!AUTEUR P. SARRAZIN DEC 81 DRPN DORVAL P.Q. CANADA
!
!LANGAGE RATFOR
!
!OBJET(SORTI)
!          VERIFIER LA VALEUR DE MODX  RESERVER MEMOIRE, OUVRIR FICHIER
!            MODX=1 FICHIER STANDARD
!                 2 FICHIER ACCES DIRECT (READMS)
!                 3 FICHIER SEQUENTIEL
!                 4 FICHIER SEQUENTIEL AVEC PARAMETRES DE FSTECR
!                 5 Fichier sequentiel ascii (sortie(formatee))
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!  IN      MODX      1=FICHIER STANDARD
!                    2=FICHIER DIRECT (READMS)
!                    3=FICHIER SEQUENTIEL
!
!  IN     NOENRG     NOMBRE D'ENREGISTREMENT DANS FICHIER
!  IN      JWRIT     -1=REECRIRE SUR FICHIER OU ECRIRE A LA FIN(MS) SORTI(MS,500,A)
!                    -1=ECRIRE SUR FICHIER UN RECORD SANS DETRUIRE UN
!                       RECORD PAREIL   SORTI(STD,500,A)
!                    +1=REECRIRE SUR FICHIER FATAL SI RECORD PAS LA.   SORTI(MS,500,R)
!                    +1=REMPLACE UN RECORD SI DEJA EXISTANT DETRUIT   SORTI(STD,500,R)
!
!IMPLICITES
!
!MESSAGES
!         DEUXIEME APPEL A LA DIRECTIVE SORTIE APPEL IGNORE
!         MAUVAIS APPEL A DIRECTIVE SORTIE FICHIER STD
!         MAUVAISE DIRECTIVE (SORTIE) FICHIER MS
!         DIRECTIVE ENREG=0, INITIALISER A 1
!         MAUVAIS APPEL A SORTIE FICHIER SEQ
!         TYPE DE FICHIER INCONNU
!         MAUVAISE DIRECTIVE MODE DIFFERENT DE (STD,MS,SEQ)
!
!
!----------------------------------------------------------------------
!
!
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      character(len=512) :: defo(990)
      character(len=512) :: listl(990), form
      character(len=512) :: lfn(990)
      common/carac/defo,listl,form,lfn
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
   integer noenrg,numero,nbrow,numdel,istart
   common/enrege/ noenrg,numero,nbrow,numdel,istart
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
      integer i,jwrit,modx,norecs
      common/relu/dejalue
      logical dejalue
      data dejalue/.false./
!
!   DATE DE LA RUN OPERATIONNELLE UTILISEE PAR L'USAGER DIRECTIVE OPDAT=OUI
!
      if (lfn(idx_date).ne. 'NON') then
         if (seldat) userdate = jdate
         if (seldat) write(6,*)' DATE ORIGINE = DATE VALIDE DATE=',userdate
      endif
!
      if (dejalue)  then
         if (message) then
            write(6,*)            'DEUXIEME APPEL A LA DIRECTIVE SORTIE APPEL IGNORE'
            return
         endif
      endif
      dejalue=.true.
!
      mode = modx
      noenrg = norecs
!
!     NBR1 = MAX0(2,FSTNBR(1))
!
!      ier = fstnbr(1)
!      if (ier .lt. 0) call fstabt
!
!   LE FICHIER D ENTRE NE PEUT ETRE  FICHIER STANDARD SEQUENTIEL
!
      iset = 0
!
!  RESERVER MEMOIR POUR FICHIER D ENTRE
!
!
      if (mode.eq.1)  then
         if (noenrg.eq.1) then
            ier=fnom(lnkdiun(idx_ozsrt),lfn(idx_ozsrt),'STD+SEQ+FTN',0)
         else if (noenrg.eq.0) then
            ier=fnom(lnkdiun(idx_ozsrt),lfn(idx_ozsrt),'STD+SEQ+REMOTE',0)
         else
            ier = fnom(lnkdiun(idx_ozsrt),lfn(idx_ozsrt),'STD+RND+REMOTE',0)
            if (nsort.eq.2) then
               iwrit=+1
            else if (nsort.eq.3) then
               iwrit=jwrit
            else
               write(6,*)               'MAUVAIS APPEL A DIRECTIVE SORTIE FICHIER STD'
               call pgsmabt
            endif
         endif
         if (ier .lt. 0) then
	   write(6,*)  'PROBLEME A L''OUVERTURE DU FICHIER DE SORTIE'
	   call pgsmabt
	 endif
      else if (mode.eq.2) then
         write(6,*)         ' LES FICHIERS "MS" NE SONT PAS SUPPORTES DANS CETTE'
         write(6,*)         ' VERSION DE PGSM'
         call pgsmabt
      else if (mode.eq.3.or.mode.eq.4) then
         ier = fnom(lnkdiun(idx_ozsrt),lfn(idx_ozsrt),'SEQ+FTN+UNF',0)
         if (nsort.ne.1)  then
            write(6,*)'MAUVAIS APPEL A SORTIE FICHIER SEQ'
            call pgsmabt
         endif
      else if (mode.eq.5) then
         if (jwrit.eq.-1) then
            ier = pgsmof(lnkdiun(idx_ozsrt),lfn(idx_ozsrt))
!            ier = fnom(2,lfn(41),'SEQ+FMT+APPEND',0)
         else
            ier = pgsmof(lnkdiun(idx_ozsrt),lfn(idx_ozsrt))
            ier = fnom(lnkdiun(idx_ozsrt),lfn(idx_ozsrt),'SEQ+FMT+R/W',0)
         endif
      else
         if (message) write(6,*)'TYPE DE FICHIER INCONNU'
         if (message) write(6,*)         ' MAUVAISE DIRECTIVE MODE DIFFERENT DE (STD,MS,SEQ)'
         return
      endif
      if (ier .lt. 0) then
        write(6,*)  'PROBLEME A L''OUVERTURE DU FICHIER DE SORTIE'
	call pgsmabt
      endif
!
!
      if (mode.eq.2) then
         write(6,*)   'LES FICHIERS "MS" NE SONT PLUS SUPPORTES SUR LES CYBER-910-920'
         call pgsmabt
      endif
      if (mode.eq.1)  then
         if (noenrg.eq.1)  then
            ier = fstouv(lnkdiun(idx_ozsrt), 'SEQ+FTN')
         else if (noenrg.eq.0) then
            ier = fstouv(lnkdiun(idx_ozsrt), 'SEQ')
            if (jwrit.eq.-1) then
               ier = fstapp(lnkdiun(idx_ozsrt),' ')
               ier = fsteof(lnkdiun(idx_ozsrt))
               print *, 'fsteof retourne', ier
            endif
         else
            ier = fstouv(lnkdiun(idx_ozsrt), 'RND')
         endif
         if (ier .lt. 0) then
           write(6,*)  'PROBLEME A L''OUVERTURE DU FICHIER DE SORTIE'
	   call pgsmabt
         endif
      endif
!
!   OUVRIR FICHIER D'ENTREE STANDARD
      if (inputmod.eq.SEQUENTIEL) then
         ier = fstouv(lnkdiun(1), 'SEQ')
      else
      ntotal_keys = 0
      do i=1,niun
         ier = fstouv(lnkdiun(i), 'RND+R/O+OLD')
         ntotal_keys = ntotal_keys + fstnbr(lnkdiun(i))
!         print *, 'ntotal_keys : ', ntotal_keys
         if (ier .lt. 0) then
            print *, '************************************************'
            print *,             '* LE FICHIER #',lfn(i), 'N''EST PAS STANDARD RANDOM'
            print *, '*************************************************'
            jdate= exfin('  PGSM  ', 'ABORT' , 'NON')
            call qqexit(13)
         endif
      enddo
      call fstlnk(lnkdiun,niun)
      endif
      if (voire) then
         if (message) then
            do i=1,niun
               ier = fstvoi(lnkdiun(i), 'RND')
            enddo
         endif
      endif
!
      return
      end
      subroutine pgsm_get_nfstkeys(nkeys)
      implicit none
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
      integer nkeys
      nkeys = ntotal_keys
      return
      end
!     ***************************************************************
!     *                        F I L T R E                          *
!     * Object :                                                    *
!     *         To filter data.                                     *
!     *                                                             *
!     * Arguments :                                                 *
!     *            IN /  ni    : x dimension of data                *
!     *            IN /  nj    : y dimension of data                *
!     *            IN /  Npass : nombre de passes pour le filtrage  * 
!     *            IN /  list  : list des nombres de filtre        *
!     *            IN /  L     : dimension de la list              *
!     *         IN/OUT/  slab  : les donnees a filtrer              *
!     *                                                             *
!     ***************************************************************
      subroutine filtre (slab, NI, NJ, nrows, Npass, list, L)
      implicit none
      integer NI, NJ, nrows
      integer l,list(L)
      real slab(NI ,NJ)
      real facteur(-4:4,5)
      real temp
      integer k,I,J, ier
      integer nb_elm
      integer Npass, pass
      integer nb_elem, lng_list, istart, iend
      real sum
      real result1(ni), result2(nj)
      nb_elem = (l+1)/2
      istart = -nb_elem + 1
      iend = nb_elem -1
      do j=1, nb_elem
         do I=istart,iend
           facteur(i,j) = 0.0
         enddo
      enddo
      do j=1, nb_elem-1
         sum = 0.0
         do i=-j,j
            sum = sum + list(I+nb_elem)
         enddo
         do i=-j,j
            facteur(i,nb_elem-j) = 1.0*(list(i+nb_elem)) / sum
         enddo
      enddo
      do pass=1, Npass
         do J=1, NJ
            do I=2, NI-1
               temp = 0
               nb_elm = min(I-1,NI-I,L/2)
               do k = -nb_elm, nb_elm
                  temp = temp + slab(I+k,J) *                  facteur(k,(L/2+1)-nb_elm)
               enddo
               result1(I) = temp
            enddo
            do I=2, NI-1
               slab(I,J) = result1(I)
            enddo
         enddo
         do I=1, NI
            do J=2, NJ-1
               temp=0
               nb_elm = min(J-1,NJ-J,L/2)
               do k = -nb_elm, nb_elm
                  temp = temp + slab(I,J+k) *                   facteur(k,(L/2+1)-nb_elm)
               enddo
               result2(J) = temp
            enddo
            do J=2, NJ-1
               slab(I,J) = result2(J)
            enddo
         enddo
      enddo
      return
      end
!
!**FONCTION SYMETRI FUNCTION QUI RECONNAIT SI LA VARIABLE EST SYMETRIQUE
!
      logical function symetri(cnom)
      implicit none
      external cmetsym
!
!AUTEUR  P.SARRAZIN  FEVRIER  DRPN  DORVAL  P.Q.  CANADA
!
!REVISION 4.0.2
!   MODIFICATION ARGUMENT D'ENTREE "NOM"
!      DE "INTEGER" A "CHARACTER*2"
!   Y. CHARTIER DRPN DORVAL QUEBEC
!LANGAGE RATFOR
!
!OBJET(SYMETRI)
!         VERIFIER SI UNE VARIABLE EST SYMETRIQUE .TRUE.=SYMETRIQUE
!         MESSAGE SI VARIABLE N EST PAS RECONNUE DEFAULT SYMETRIQUE
!          SI UNE VARIABLE N EST PAS RECONNUE ELLE EST CONSIDEREE
!          COMME SYMETRIQUE
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN    NOM    NOM DE LA VARIABLE
!
!
!APPEL    VIA MACPCP,EPAISUR
!         SYMETRI(NOMBRE)
!
!MESSAGES
!         LA SYMETRIE DE LA VARIABLE  EST INCONNUE
!         ON LA SUPPOSE SYMETRIQUE
!
!MODULES
!
! - - - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - -*
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
      integer nnoms,maxnoms,nsym,nsymm
      logical ssym(256)
      common /symnom/ nnoms,maxnoms,nsym,nsymm,ssym
      character*4 noms(256)
      common /csymnom/ noms
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
      integer  i
      character*4 cnom
!
      symetri = .true.
!
      do i = 1,nnoms
         if (cnom.eq.noms(i)) then
            symetri = ssym(i)
            return
         endif
      enddo
!
      if (message) then
         if (cgrtyp == 'A'.or.cgrtyp == 'B'.or.cgrtyp =='G') then
            write(6,101) cnom
         endif
      endif
 101  format(//2x,' **** LA SYMETRIE DE LA VARIABLE- ',       a4,' - EST INCONNUE',      ' ON LA SUPPOSE SYMETRIQUE *****'//)
!
      call cmetsym(cnom, .true.)
      return
      end
      subroutine testseq
      implicit none
      integer nigem, njgem, nigauss, njgauss, iun
      parameter (nigem = 160)
      parameter (njgem = 1)
      integer npac,idat, ideet, npas, ni, nj, nk, ip1, ip2o, ip3o, llg1, llg2, llg3, llg4, idatyp
      integer ier, fnom
      external fnom
      real champ(nigem,njgem)
      character*24 chaine
      character*4 nomvar, cbidon2
      character*2 typvar, bidon1, grtyp, grref, cbidon1
      character*12 etiket, cbidon8
      iun = 0
      ier = fnom(iun,'bofseq','SEQ+FTN+UNF',0)
      rewind(iun)
 10   read (iun,err=13) npac,idat, ideet, npas, ni, nj, nk,                ip1, ip2o, ip3o, llg1, llg2, llg3, llg4, idatyp,        &
     &       chaine
      print *,  npac,idat, ideet, npas, ni, nj, nk,                ip1, ip2o, ip3o, llg1, llg2, llg3, llg4, idatyp,               c&
     &haine
      read (iun) champ
      print *, champ
      goto 10
 13   continue
      stop
      end
!
!**s/p uvectur interpolation des vecteurs u-v (horizontalement)
!
   subroutine uvectur (cnom1, cnom2, cnom3,iheur, npar, itabuv)
      implicit none
   external ecritur,cigaxg,fstinf,fstinl,ipgsmlic,pgsmlic,ipgsmlir,pgsmlir,memoir, fstcvt,fstprm,pgsmabt,imprime,vdauv, incdat,fsto&
     &pc,messags
   external liraxez,  cxgaig
   integer fstinf,fstinl,ipgsmlic,pgsmlic,ipgsmlir,pgsmlir,fstprm,fstopc,fstcvt
   integer ezqkdef, ezwdint, ezuvint, ezdefset, ezsint
   real fbidon
!
!auteur  p.sarrazin fevrier 82 drpn  dorval p.q. canada
!revision 4.0.2
!   conversion des variables hollerith en caracteres
!   y. chartier -aout 90- drpn dorval quebec
!
!langage ratfor
!
!objet(uvectur)
!         interpolation des vecteurs u et v
!
!librairies
!         -source  armnsrc,drpn
!         -objet   pgsmlib,id=armnpjs.
!
!arguments
!  in    nom1   nom du premier vecteur  ex:"uu","us"...
!  in    nom2   nom du deuxieme vecteur ex:"vv","vs"...
!  in    nom3   nom du champ a ecrire apres interpolation du vent ex:"uv"
!  in    iheur   heure de la variable
!  in    npar    nombre de locations dnas itabuv
!  in    itabuv table contenant les noms (niveau)
!
!appel
!         -via routine champ
!         call uvectur(iheur, npar, itabuv)
!
!
!modules  fstinf,pgsmabt,memoir,fstprm,pgsmlir,cuvint,cigaxg,cspauv,ecritur
!
!messages
!         mauvaise directive champ (uvectur)
!         record n'existe pas sur fichier d'entre (u,v) (uvectur)
!         aucune interpolation horizontale u v
!
!----------------------------------------------------------------------
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
      integer      lnkdiun(990),niun,inputmod
      integer      sequentiel,random, ntotal_keys
      integer idx_ozsrt, idx_isll, idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v
      integer iun_isll, isll_input
      parameter (sequentiel=1,random=2)
      common /lnkflds/ lnkdiun, niun,inputmod, ntotal_keys, idx_ozsrt, idx_isll, &
             idx_i, idx_l, idx_date, idx_msglvl, idx_isent, idx_impos, idx_v, &
             iun_isll, isll_input
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   integer noenrg,numero,nbrow,numdel,istart
   common/enrege/ noenrg,numero,nbrow,numdel,istart
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      integer ichck,nlire,najou,nenle,nmoys,necrt,nmod, nraci,npfo,multp
      common/chck/ ichck,nlire,najou,nenle,nmoys, necrt,nmod,nraci,npfo,multp
!
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
      integer ig1la,ig2la,ig3la,ig4la,nilo,njlo,nklo
      integer ig1lo,ig2lo,ig3lo,ig4lo
      common /tp12ig/ ig1la,ig2la,ig3la,ig4la,nilo,njlo,nklo,ig1lo,ig2lo,ig3lo,ig4lo
!
!---------------------------------------------------------------
   character*12 cetiket,cetike
   character*4 cnomvar,cnom1,cnom2, cnom3
   character*1 cigtyp
   character*2 ctypvar
   integer i,j
   integer jp1,jp2,jp3,liljt,ni,nj,nk,npar
   integer, dimension(:), allocatable :: listniv
   real, allocatable, dimension(:,:) :: latgdin, longdin
   integer itabuv(npar) ,deet,ig1,ig2,ig3,ig4
   integer iheur,ilop,iprs,irec_uu,irec_vv,iopc
   integer numu,numv,infon,dat,datdv
   integer cnbits,cdatyp,cswa,clng,cdltf,cubc,extra1,extra2,extra3,total_keys
   integer nom2
   integer tableau(1)
   integer ig1zz, ig2zz, ig3zz, ig4zz, ezsetgdout, gdwdfuv, gdll, npts
   real    d60dum, pidum, pjdum
   real xlat1,xlon1,xlat2,xlon2
   real dgtord,dumfld,xg1,xg2,xg3,xg4,datev
   integer iunit
   logical ssw
   integer entier_ou_reel
   real*8 delta_t
   real, dimension(:,:), pointer :: uuout,vvout
   iunit = lnkdiun(1)
   nk = 1
   call pgsm_get_nfstkeys(total_keys)
   allocate(listniv(total_keys))
   call chk_userdate(datev)
!
!
   do iprs = 1,npar
!
!     trouver record pour u,v  ou us,vs .....
!
!
!     modification de hollerith a caractere
!
      if (etikent(1) .ne. -1) then
         write(cetiket,'(3A4)') (etikent(i), i=1,nwetike)
      else
         cetiket = '            '
      endif
      if (typeent .ne. -1) then
         write(ctypvar, '(A2)') typeent
      else
         ctypvar = '  '
      endif
      ier = fstinl(iunit,ni,nj,nk,datev,cetiket,itabuv(iprs),iheur,ip3ent,ctypvar,cnom1,listniv,infon,total_keys)
      if (ier .lt. 0 .or. infon.eq.0) then
         write(6,610) cnom1
 610     format(' AUCUN RECORD SUR FICHIER (FSTINL-UVECTUR) NOM=',a2)
         cycle
      endif
      if (nk.gt.1) then
         write(6,*)'***********************************************'
         write(6,*)'         PGSM N ACCEPTE PAS UN          '
         write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (UVECTUR)'
         write(6,*)'***********************************************'
         call pgsmabt
      endif
!
!
      do ilop=1,infon
!
	 entier_ou_reel = 1
         irec_uu=listniv(ilop)
!
!     identifier parametres champ nom1
!
         cetike = '            '
         ier = fstprm( irec_uu, dat,deet,npas,ni, nj, nk, cnbits,cdatyp,jp1,jp2, jp3,ctypvar, cnomvar,cetike,cigtyp, ig1,ig2,ig3,ig&
     &4, cswa, clng, cdltf, cubc, extra1, extra2, extra3)
         npack_orig = -cnbits
         if (ier .lt. 0) then
            write(6,*)' IER = FSTPRM NEGATIF VOIR UVECTUR'
         endif
!
!     verifier si grille gaussienne ni doit etre pair
!
 675     format(' ITYP=',a1,'   CIGTYP DE FSTPRM= ',a1)
         if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  call messags(ni)
!
!     calcul la date pour le record de la variable nom2
!
         delta_t = deet*npas/3600.0
         call incdatr(datdv,dat,delta_t)
         irec_vv = fstinf(iunit,ni,nj,nk,datdv,cetike,jp1,jp2,jp3,ctypvar,cnom2)
         if (irec_vv .lt. 0) then
            write(6,610) nom2
            call pgsmabt
         endif
         if (nk.gt.1) then
          write(6,*)'********************************************'
          write(6,*)'         PGSM N ACCEPTE PAS UN          '
          write(6,*)' CHAMP DE 3 DIMENSIONS NK>1 ?? (UVECTUR)'
          write(6,*)'********************************************'
          call pgsmabt
       endif
!!  Switch pour champs masques
      if (masque == 1) then
         if (ctypvar(1:1) == '@') then
            cycle 
         else if (ctypvar(2:2) == '@') then
            call uvecteur_masque(irec_uu, irec_vv)
            cycle
         endif
      endif
!     allouer memoire
      if (cdatyp == 2 .or. cdatyp == 130 .or. cdatyp == 4 .or. cdatyp == 132) then
	allocate(itmpif1(ni,nj))
	allocate(itmpif2(ni,nj))
	allocate(itmpif3(li,lj))
	allocate(itmpif4(li,lj))
	entier_ou_reel = 2
      endif
      allocate(tmpif1(ni,nj))
      allocate(tmpif2(ni,nj))
      allocate(tmpif3(li,lj))
      allocate(tmpif4(li,lj))
!
!     lire champ nom1
!
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
      if (entier_ou_reel == 2) then
	numu =ipgsmlir(itmpif1,1,ni,nj,nk,datdv,cetike,jp1,jp2,jp3,ctypvar,cnom1,cigtyp)
      else
	numu = pgsmlir(tmpif1,1,ni,nj,nk,datdv,cetike,jp1,jp2,jp3,ctypvar,cnom1,cigtyp)
      endif
      if (printen)  call imprime(cnom1,tmpif1,ni,nj)
      if (.not.message) iopc= fstopc('TOLRNC','DEBUGS',.true.)
      if (entier_ou_reel == 2) then
	numv = ipgsmlic(itmpif2,1,ni,nj,nk,datdv,cetike,jp1,jp2,jp3,ctypvar,cnom2,ig1,ig2,ig3,ig4,cigtyp)
      else
	numv = pgsmlic(tmpif2,1,ni,nj,nk,datdv,cetike,jp1,jp2,jp3,ctypvar,cnom2,ig1,ig2,ig3,ig4,cigtyp)
      endif
      if (cigtyp.eq.'G'.and.mod(ni,2).ne.0)  call messags(ni)
!
      if (printen)  call imprime(cnom2,tmpif2,ni,nj)
!****************************************************************
!     si vvent=.true. on calcule la vitesse du vent
      if (entier_ou_reel == 2) then
	 call cvtrfi(tmpif1, itmpif1, ni, nj)
	 call cvtrfi(tmpif2, itmpif2, ni, nj)
      endif
      gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
      if (vvent) then
         ssw=.false.
         if (gdout == gdin) then
            if (ctypvar == '@@') then
                uuout => tmpif1
                vvout => tmpif2
            else
              if (wdvent) then
                ssw = .true.
                allocate(latgdin(ni,nj),longdin(ni,nj))
                ier = gdll(gdin, latgdin, longdin)
                npts = ni * nj
                ier = gdwdfuv(gdin, tmpif3, tmpif4, tmpif1, tmpif2, latgdin, longdin, npts)
                cgrtyp = cigtyp
                lg1 = ig1
                lg2 = ig2
                lg3 = ig3
                lg4 = ig4
                uuout => tmpif3
                vvout => tmpif4
              else
                uuout => tmpif1
                vvout => tmpif2
                uuout = sqrt(uuout*uuout+vvout*vvout)
              endif
            endif
         else
            ier = ezdefset(gdout, gdin)
            if (ctypvar == '@@') then
               ier = ezsint(tmpif3, tmpif1)
               ier = ezsint(tmpif4, tmpif2)
            else
               ier = ezwdint(tmpif3, tmpif4, tmpif1, tmpif2)
            endif
            uuout => tmpif3
            vvout => tmpif4
         endif
	 if (entier_ou_reel == 2) then
	    call cvtifr(itmpif3, uuout, li, lj)
	    call  iecritur(itmpif3, npack, dat, deet, npas, li, lj, nk, jp1, jp2, jp3, &
	      ctypvar, cnom3, cetike, cgrtyp, lg1, lg2, lg3, lg4)
	 else
	    call  ecritur(uuout, npack, dat, deet, npas, li, lj, nk, jp1, jp2, jp3, &
	    ctypvar, cnom3, cetike, cgrtyp, lg1, lg2, lg3, lg4)
         endif
      if (wdvent) then
         do j=1,lj
         do i=1,li
            if (tmpif4(i,j).lt.0.0) then
               tmpif4(i,j) = tmpif4(i,j) + 360.0
            endif
         enddo
         enddo
         vvout => tmpif4
         if (entier_ou_reel == 2) then
	    call cvtifr(itmpif4, vvout, li, lj)
	    call iecritur(itmpif4, npack, dat, deet, npas, li, lj, nk, jp1, jp2, jp3,&
		ctypvar,'WD  ',cetike,cgrtyp,lg1,lg2,lg3,lg4)
          else
	    call ecritur(vvout, npack, dat, deet, npas, li, lj, nk, jp1, jp2, jp3,&
		ctypvar,'WD  ',cetike,cgrtyp,lg1,lg2,lg3,lg4)
         endif
      endif
!
!****************************************************************
!
   else
!
!     on ne fait pas d'interpolation si igtyp=grtyp  ig1=lg1  ig2=lg2
!     ig3=lg3  ig4=lg4
!
            if (cigtyp.ne.cgrtyp.or.ig1.ne.lg1.or.ig2.ne.lg2.or.ig3.ne.lg3.or.ig4.ne.lg4.or.li.ne.ni.or.lj.ne.nj) then
!
!     interpolation u,v vecteur a vitesse et direction du vent
!
!     si ssw = vrai interpoler vitesse et direction
!     faux interpoler seulement vitesse
!
               ssw = .true.
               ier = ezdefset(gdout, gdin)
!
!              cas special pour typvar = @@
               if (ctypvar == '@@') then
                  ier = ezsint(tmpif3, tmpif1)
                  ier = ezsint(tmpif4, tmpif2)
               else
                  ier = ezuvint(tmpif3, tmpif4, tmpif1, tmpif2)
               endif
               uuout => tmpif3
               vvout => tmpif4
!
!     apres interpolation horizontale passer de vitesse et direction
!     aux composantes u et v
!
!     si type de grille "x",    u-v interpolation n\s - e\o
!
            else
               deallocate(tmpif3)
               deallocate(tmpif4)
               uuout => tmpif1
               vvout => tmpif2
               if (message) then
                 write(6,*)'AUCUNE INTERPOLATION HORIZONTALE '
              endif
           endif
!
!     ecrire vecteur u
!
         if (entier_ou_reel == 2) then
	    call cvtifr(itmpif3, uuout, li, lj)
	    call iecritur(itmpif3, npack,dat,deet,npas,li,lj,nk,jp1,jp2,jp3,ctypvar,cnom1,cetike,cgrtyp,lg1,lg2,lg3,lg4)
	 else
           call ecritur(uuout,npack,dat,deet,npas,li,lj,nk,jp1,jp2,jp3,ctypvar,cnom1,cetike,cgrtyp,lg1,lg2,lg3,lg4)
         endif
!
!
!     ecrire vecteur v
!
         if (entier_ou_reel == 2) then
	    call cvtifr(itmpif4, vvout, li, lj)
	    call iecritur(itmpif4, npack,dat,deet,npas,li,lj,nk,jp1,jp2,jp3,ctypvar,cnom2,cetike,cgrtyp,lg1,lg2,lg3,lg4)
	 else
            call ecritur(vvout,npack,dat,deet,npas,li,lj,nk,jp1,jp2,jp3,ctypvar,cnom2,cetike,cgrtyp,lg1,lg2,lg3,lg4)
         endif
!
!     fin du calcul des composantes
!
!
         endif
         if (associated(tmpif3)) deallocate(tmpif3)
         if (associated(tmpif4)) deallocate(tmpif4)
         if (associated(tmpif1)) deallocate(tmpif1)
         if (associated(tmpif2)) deallocate(tmpif2)
         if (allocated(latgdin)) deallocate(latgdin)
         if (allocated(longdin)) deallocate(longdin)
         if (entier_ou_reel == 2) then
	    if (associated(itmpif1)) deallocate(itmpif1)
	    if (associated(itmpif2)) deallocate(itmpif2)
	    if (associated(itmpif3)) deallocate(itmpif3)
	    if (associated(itmpif4)) deallocate(itmpif4)
         endif
99999 continue
      enddo
!
!
!     reinitialiser clef de controle
!
   enddo    
      vvent=.false.
!
!
      deallocate(listniv)
      return
      end
!
!**s/p uvectur interpolation des vecteurs u-v (horizontalement)
      subroutine uvecteur_masque(key_uu, key_vv)
      implicit none
      integer :: key_uu, key_vv
      external ecritur,pgsmluk,fstinf,fstsui,memoir,fstprm,qaaqr,fstcvt, &
         fstsel,symetri,imprime,itrouve,messags,pgsmabt
      external cvtr2i
      external liraxez
      integer  pgsmluk, fstinf, fstsui, fstprm, fstcvt, fstsel, fstinl, fstluk
      integer ezgdef_fmem, ezqkdef, ezuvint, ezuvint_mdm, ezwdint, ezdefset, fst_get_mask_key, key_mask
      logical skip
!
      integer nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas
      integer igg1,igg2,igg3,igg4,if0,icnt,ier
      logical unefois,once
      common /accum/ nni,nnj,nnk,idatt,jpp1,jpp2,jpp3,ideet,npas,igg1,igg2,igg3,igg4,if0,icnt,unefois,once
      character*12 cetik
      character*2 ctypv
      character*1 cigty
      character*4 cnumv
      common /caccum/ cetik, ctypv, cigty, cnumv
!
!
      logical valid,vvent,wdvent,seldat
      integer userdate,date2,date3,dateval,etikent(3),nwetike
      common /dates/ valid,vvent,wdvent,seldat,userdate,date2,date3,dateval,etikent,nwetike
      character*4 cnomqq,cnomqr,cnommt
      common /cdates/ cnomqq, cnomqr, cnommt
!
!
!
   character cdumnom*4, cdumtyp*2, cdumetk*12, cdumgtp*1
   common /cdummys/ cdumnom, cdumtyp, cdumetk, cdumgtp
   integer dumnom, dumtyp, dumetk, dumgtp
   common /dummys/  dumnom, dumtyp, dumetk, dumgtp
!
!
!
   integer typesrt,etiksrt(3),ip3srt, niis,njjs,niif,njjf,niinc,njjnc,ip2srt,nwetiks,compression_level
   logical printsr
   common /ecrires/ ip3srt,typesrt,etiksrt,nwetiks,printsr,niis, njjs,niif,njjf,niinc,njjnc,ip2srt, compression_level
!
!
!
      integer niz, njz, nkz
      integer ig1z, ig2z, ig3z, ig4z
      integer ig1ref, ig2ref, ig3ref, ig4ref
      real, pointer, dimension(:) :: axex, axey
      common /gdz/ niz, njz, nkz, ig1z, ig2z, ig3z, ig4z, ig1ref, ig2ref, ig3ref, ig4ref, axex, axey
      character*1 grref
      common /gdzchar/ grref
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
   integer iwrit,iset,iroot, indxs
   integer if1,if2,if3,if4,if5,if6,if7,if8,if9,ift,mode
   integer nsort,nggr,nlalo,jdate,ip4
   common /indptr/iwrit,iroot,indxs,iset,if1,if2,if3,if4,if9,ift,mode,nsort,nggr,nlalo,jdate,ip4
!
!
!
      integer typeent,ip3ent,mtype,nis,njs,nif,njf,ninc,njnc,diese
      logical printen,mtdone
      common /lires/ ip3ent,typeent,mtype,printen,nis,njs,nif,njf,ninc,njnc,diese,mtdone
!
!
!
   real, dimension(:,:), pointer :: tmpif0, tmpif1, tmpif2, tmpif3, tmpif4
   integer, dimension(:,:), pointer :: itmpif0, itmpif1, itmpif2, itmpif3, itmpif4
   real, dimension(:,:), pointer :: tmpift, tmplat, tmplon
   real, dimension(:,:), pointer :: tmplatg, tmplong
   real, dimension(:), pointer :: tmproot
   common /llccmm/ tmpif0,tmpif1,tmpif2,tmpif3,tmpif4, &
                   itmpif0,itmpif1,itmpif2,itmpif3,itmpif4, &
                   tmpift,tmproot,&
                   tmplat,tmplon, tmplatg,tmplong
!
!
!
      integer npack,npack_orig
      common/packin/npack,npack_orig
!
!
!
      integer npair,npairuv
      common/pair/npair,npairuv
      character*24 paire(40)
      common/cpair/paire
!
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
   character(len=12):: cetiket
   character(len=4) :: cnom_uu, cnom_vv
   character(len=2) :: ctypvar
   character(len=1) :: cigtyp
   integer i, j, nunv, itrouve, ii, key_mask_out, lk
   integer deet, ig1, ig2, ig3, ig4, iheur
   integer dateo, datev, nbits, datyp, ip1, ip2, ip3, iun_out
   integer iprs, irec, iunit, ne, ni, nj, nk, total_keys, nrecs
   integer cnbits, cdatyp, cswa, clng, cdltf, cubc, extra1, extra2, extra3
   integer dateo_mask, deet_mask, npas_mask, ni_mask, nj_mask, nk_mask
   integer nbits_mask, datyp_mask, ip1_mask, ip2_mask, ip3_mask
   integer ig1_mask, ig2_mask, ig3_mask, ig4_mask, cswa_mask, clng_mask, cdltf_mask, cubc_mask
   integer datev_mask, extra2_mask, extra3_mask
   logical sym, symetri
   character(len=4)  :: cnomvar_mask
   character(len=2)  :: ctypvar_mask
   character(len=12) :: cetiket_mask
   character(len=1)  :: cigtyp_mask
   real fbidon
   real, dimension(:,:), allocatable, target :: tmp_uuin, tmp_uuout, tmp_vvin, tmp_vvout
   integer, dimension(:,:), allocatable, target :: masque_in, masque_out
   real, dimension(:,:), pointer :: uu_in, vv_in, uu_out, vv_out
   integer, dimension(:,:), pointer :: tmpmsk
   logical masque_present, masque_done, write_masque, ssw
   nunv=0
   iunit=1
   masque_present = .false.
   masque_done = .false.
   ctypvar_mask = '@@'
   ier = fstprm(key_uu, dateo, deet, npas, ni, nj, nk, nbits, datyp, ip1, ip2, ip3, ctypvar, cnom_uu, cetiket, &
            cigtyp, ig1, ig2, ig3, ig4, cswa, clng, cdltf, cubc, datev, extra2, extra3)
   ier = fstprm(key_vv, dateo, deet, npas, ni, nj, nk, nbits, datyp, ip1, ip2, ip3, ctypvar, cnom_vv, cetiket, &
            cigtyp, ig1, ig2, ig3, ig4, cswa, clng, cdltf, cubc, datev, extra2, extra3)
   allocate(tmp_uuin(ni,nj), tmp_vvin(ni,nj), tmp_uuout(li,lj), tmp_vvout(li,lj))
   ier = fst_get_mask_key(key_mask, key_uu, 0, iunit)
   if (key_mask >= 0) then
      masque_present = .true.
      allocate(masque_in(ni,nj), masque_out(li,lj))
   else
      masque_present = .false.
   endif
   ier = pgsmluk(tmp_uuin, key_uu, ni,nj,nk,cnom_uu,cigtyp)
   ier = pgsmluk(tmp_vvin, key_vv, ni,nj,nk,cnom_vv,cigtyp)
   if (printen)  call imprime(cnom_uu,tmp_uuin,ni,nj)
   if (ig1 /= 0) sym = symetri(cnom_uu)
!  interpolation ordinaire sans masque
   if (.not.masque_present) then
!     si vvent=.true. on calcule la vitesse du vent
      gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
      ier = ezdefset(gdout, gdin)
      if (vvent) then
         ssw=.false.
         if (wdvent) then
            ssw = .true.
         endif
         ier = ezwdint(tmp_uuout, tmp_vvout, tmp_uuin, tmp_vvout)
         uu_out => tmp_uuout
         call  ecritur(uu_out,npack,dateo,deet,npas,li,lj,nk,ip1,ip2,ip3,ctypvar,'UV  ',cetiket,cgrtyp,lg1,lg2,lg3,lg4)
         if (wdvent) then
            do j=1,lj
               do i=1,li
                  if (tmp_vvout(i,j).lt.0.0) then
                     tmp_vvout(i,j) = tmp_vvout(i,j) + 360.0
                  endif
               enddo
            enddo
            vv_out => tmp_vvout
            call ecritur(vv_out,npack,dateo,deet,npas,li,lj,nk,ip1,ip2,ip3,ctypvar,'WD  ',cetiket,cgrtyp,lg1,lg2,lg3,lg4)
         endif
      else
         if (cigtyp.ne.cgrtyp.or.ig1.ne.lg1.or.ig2.ne.lg2.or.ig3.ne.lg3.or.ig4.ne.lg4.or.li.ne.ni.or.lj.ne.nj) then
            ssw = .true.
            ier = ezuvint(tmp_uuout, tmp_vvout, tmp_uuin, tmp_vvout)
            uu_out => tmp_uuout
            vv_out => tmp_vvout
         else
            deallocate(tmp_uuout)
            deallocate(tmp_vvout)
            uu_out => tmp_uuin
            vv_out => tmp_vvout
            if (message) then
              write(6,*)'AUCUNE INTERPOLATION HORIZONTALE '
            endif
         endif
         call ecritur(uu_out,npack,dateo,deet,npas,li,lj,nk,ip1,ip2,ip3,ctypvar,cnom_uu,cetiket,cgrtyp,lg1,lg2,lg3,lg4)
         call ecritur(vv_out,npack,dateo,deet,npas,li,lj,nk,ip1,ip2,ip3,ctypvar,cnom_vv,cetiket,cgrtyp,lg1,lg2,lg3,lg4)
      endif
   else  
!
      if (cigtyp /= cgrtyp.or.cigtyp == 'Z'.or.ig1 /= lg1.or.ig2 /= lg2.or.ig3 /= lg3.or.ig4 /= lg4.or.li /= ni.or.lj /= nj) then
         gdin = ezqkdef(ni, nj, cigtyp, ig1, ig2, ig3, ig4, iunit)
         ier = ezdefset(gdout, gdin)
         ier = fstprm(key_mask, dateo_mask, deet_mask, npas_mask, ni_mask, nj_mask, nk_mask, &
                  nbits_mask, datyp_mask, ip1_mask, ip2_mask, ip3_mask, &
                  ctypvar_mask, cnomvar_mask, cetiket_mask, &
                  cigtyp_mask, ig1_mask, ig2_mask, ig3_mask, ig4_mask, cswa_mask, clng_mask, &
                  cdltf_mask, cubc_mask, datev_mask, extra2_mask, extra3_mask)
         ier = fstluk(masque_in, key_mask, ni,nj,nk)
         ier = ezuvint_mdm(tmp_uuout, tmp_vvout, masque_out, tmp_uuin, tmp_vvin, masque_in)
         uu_out => tmp_uuout
         vv_out => tmp_vvout
         tmpmsk => masque_out
      else
         uu_out => tmp_uuin
         vv_out => tmp_vvin
         tmpmsk => masque_in
         if (message) write(6,662) cnom_uu
 662           format(2x,'AUCUNE INTERPOLATION HORIZONTALE CHAMP=',a4)
      endif
      call ecritur(uu_out,npack,dateo,deet,npas,li,lj,nk, ip1, ip2, ip3, &
         ctypvar, cnom_uu, cetiket, cgrtyp, lg1, lg2, lg3, lg4)
      call ecritur(vv_out,npack,dateo,deet,npas,li,lj,nk, ip1, ip2, ip3, &
         ctypvar, cnom_vv, cetiket, cgrtyp, lg1, lg2, lg3, lg4)
      
      iun_out = 2
      key_mask_out = fstinf(iun_out, ni_mask, nj_mask, nk_mask, datev_mask, cetiket_mask, ip1_mask, ip2_mask, ip3_mask, &
         ctypvar_mask, cnomvar_mask)
      write_masque = .true.
      if (key_mask_out >= 0) then
         write_masque = .false.
      endif
      if (masque_present.and.write_masque) then
         call iecritur(tmpmsk,-nbits_mask,dateo_mask,deet_mask,npas_mask,li,lj, nk,&
            ip1_mask,ip2_mask, ip3_mask, ctypvar_mask, cnomvar_mask, cetiket_mask, &
            cgrtyp, lg1, lg2, lg3, lg4)
         if (allocated(masque_out)) then
            deallocate(masque_out, masque_in)
         endif
      endif
      deallocate(tmp_uuout, tmp_uuin, tmp_vvout, tmp_vvin)
   endif
!
   if (nunv > 0) then
      write(6,666)
 666     format(' AUCUNE INTERPOLATION SUR VARIABLE PAIRE CHAMP(TOUT,TOUT)')
      write(6,668)
 668     format(' ON DOIT UTILISER LE NOM DE LA VARIABLE EX: CHAMP(UU,TOUT)')
      write(6,669)
 669  format(' ATTENTION L INTERPOLATION DES VECTEURS SERA SCALAIRE (!!!)')
!
   endif
   return
   end
!
!**S/P VDAUV  U ET V DE DIRECTION ET VITESSE DU VENT
!
      subroutine vdauv(srtentu,srtentv,clong,dgtord,nombre)
      implicit none
      integer nombre
      real srtentu(nombre),srtentv(nombre),clong(nombre),dgtord
!
!AUTEUR P. SARRAZIN DORVAL QUEBEC CANADA (DRPN)
!
!OBJET(VDAUV)
!         AVEC LA DIRECTION ET LA VITESSE DU VENT CALCUL LES VECTEURS 
!         U ET V. UTILISATION DE LA LONGITUDE
!
!ARGUMENTS
! IN-OUT   SRTENTU - ENTRE DIRECTION   SORTI VECTEUR U
! IN-OUT   SRTENTV - ENTRE VITESSE  SORTI VECTEUR V
!   IN     CLONG   - CHAMP DE LONGITUDES EN DEGRE 0-360
!   IN     DGTORD  - FACTEUR DE CONVERSION DE DEGREE A RADIAN
!   IN     NOMBRE  - NOMBRE DE POINTS DANS LES DEUX CHAMPS
!
!APPEL
!     - VIA UVECTUR 
!     - CALL VDAUV(SRTENTU,SRTENTV,CLONG,DGTORD,NOMBRE)
!
   integer li,lj,lg1,lg2,lg3,lg4,ixlat,ixlon,ixlatg,ixlong,ngr,ncoords
   real clatmin,clatmax,clonmin,clonmax,dgrwxy
   integer gdin, gdout, masque
   integer gr_a, gr_b, gr_g, gr_ps, gr_tape4, gr_latlon, gr_tape1, gr_tape2, gr_xylis, gr_xydir, gr_lldir, gr_lllist, &
      gr_grib, gr_stereo, gr_comme, gr_stations, gr_gem, gr_gef
   integer nmaxcoords
   parameter (gr_a         = 1)
   parameter (gr_latlon    = 2)
   parameter (gr_ps        = 3)
   parameter (gr_tape4     = 4)
   parameter (gr_g         = 5)
   parameter (gr_b         = 6)
   parameter (gr_tape1     = 7)
   parameter (gr_tape2     = 8)
   parameter (gr_xylis     = 9)
   parameter (gr_xydir    = 10)
   parameter (gr_lldir    = 11)
   parameter (gr_lllist   = 12)
   parameter (gr_grib     = 10)
   parameter (gr_stereo   = 13)
   parameter (gr_comme    = 12)
   parameter (gr_stations = 14)
   parameter (gr_gem      = 15)
   parameter (gr_gef      = 15)
   parameter (nmaxcoords  = 100000)
   real coordll(nmaxcoords,2)
   common /grilles / li, lj, lg1, lg2, lg3, lg4, gdin, gdout, ixlat, ixlon,ngr, clatmin,clatmax,clonmin,clonmax, &
      dgrwxy,ixlatg,ixlong,ncoords,coordll,masque
   character*1 cgrtyp, cgtypxy
   common /cgrille/ cgrtyp, cgtypxy
!
!
!
!
!
!----------------------------------------------------------------------
!
      external pgsmabt,messags
!
      real angle,u,v
      integer i
!
!    SI LE TYPE DE GRILLE  "L"
!
      if (cgtypxy.eq.'L') then
         do i=1,nombre
            angle=dgtord*(srtentv(i) - clong(i))
            u=srtentu(i)*sin(angle)
            v=-srtentu(i)*cos(angle)
            srtentu(i)=u
            srtentv(i)=v
         enddo
!     
!     SI LE TYPE DE GRILLE GRTYPXY= "N" HEM NORD
!     
      else if (cgtypxy.eq.'N') then
         do i=1,nombre
            angle=dgtord*(dgrwxy + srtentv(i))
            u=srtentu(i)*cos(angle)
            v=srtentu(i)*sin(angle)
            srtentu(i)=u
            srtentv(i)=v
         enddo
      else
         write(6,*) ' TYPE DE GRILLE PAS "L" OU "N" '
         write(6,*)          ' DANS ROUTINE VDAUV PAS DE CODE VALID POUR TYPE "S"'
         call pgsmabt
      endif
      return
      end
!
!**S/P  VERLALO  VERIFI LONGITUDE ET LATITUDE
!     
      subroutine verlalo(clat,clon,nombre)
      implicit none
      integer nombre
      real clat(nombre),clon(nombre)
!
!AUTEUR P. SARRAZIN DORVAL QUEBEC CANADA (DRPN)
!
!OBJET(VERLALO)
!         VERLALO - VERIFIER LES LATITUDES ET LONGITUDES DU CHAMP
!
!
!LIBRAIRIES
!         -SOURCE  ARMNSRC,DRPN
!         -OBJET   PGSMLIB,ID=ARMNPJS.
!
!ARGUMENTS
!   IN    CLAT   - CHAMP DE LATITUDE
!   IN    CLON   - CHAMP DE LONGITUDES
!   IN    NOMBRE - NOMBRE DE POINTS DANS LES CHAMPS CLAT/CLON
!
!-------------------------------------------------------------------- 
      external pgsmabt,messags
!
!
      logical voire,voirs,pose,message
      common/voir/ voire,voirs, pose, message
!
!
!
!
!
      integer i
!
!
!   VERIFIER SI LATITUDES ET LONGITUDES SONT A L'INTERIEUR
!   DES LIMITES
!
       do i=1,nombre
          if (clon(i).lt.0.0) clon(i)=clon(i)+360.0
          if (clon(i).ge.360.0) clon(i)=clon(i) - 360.0
!
          if (clat(i).lt.-90.005.or.clat(i).gt.90.005) then
             write(6,*)' ROUTINE VERLALO MAUVAISE LATITUDE=',clat(i)
          endif
          clat(i)=max(-90.0,min(90.0,clat(i)))
       enddo
       return
       end
