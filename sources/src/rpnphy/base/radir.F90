!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------

subroutine radir9(taux, fd, fm, fnuage, t, q, sh, ts, ps, nn, &
     nk, oz, nkmax, ni2, nl,ni1, nt, mx, mxx, no3, &
     ncx, nco2,g1,g2,g3,th2o,tro3,yg3,bcn,dbcn,bo3,dbo3, &
     to3, uu, f, ku, opnua, &
     xin, v, su, itemp, xpro, s, suo3, del, au, az, &
     nhaut, nmoy, nbas, &
     uco2j, tco2j, wkco2j,trmin,tmem, &
     flsh,flss,fldel,flnk,flnn)
   use, intrinsic :: iso_fortran_env, only: REAL64
   use tdpack_const, only: CPD, GRAV, STEFAN
   use phy_options
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   !     variables a une dimension de travail
   integer nn,nkmax,ni2,nl,ni1,nt,mx,i1,i2,m,j,jind, &
        mxx,no3,ncx,nk,ip,i,l,it
   integer iwing,icen,nco2
   !     variables dimensionnees
   real to,t1,t2,ap,po,elsa,epsil,corfac,eo3,eh2o, &
        borne,z,aa,bb,xnu,bz,a,tx,wc,zx, &
        r,y,wc2

   real tro3(mxx),to3(no3),g1(mxx,nt),g2(mxx,nt), &
        g3(mxx,nt),th2o(mxx,nco2),bcn(nt,nco2),dbcn(nt,nco2), &
        uu(mxx),yg3(mxx,ncx,nco2),bo3(nt),dbo3(nt)

   real trnuage,diftlim
   real sh(ni2,nk)
   real taumax,secday
   real voigth2o,voigto3
   integer flnk, flnn

   !     variables dependant de ni2 ou ni1
   real taux(ni2,flnk),fnuage(ni2,nk),t(ni1,nk+1),q(ni1,nk), &
        fm(ni2,flnn),fd(ni2,flnn),oz(ni2,nk), &
        ts(nl),ps(nl), f(ni2,nkmax,nkmax)
   integer ku(ni2,nkmax,nkmax) ,itemp(ni2,nn)
   real xin(ni2,flnn), v(ni2,flnn), su(ni2,nn), &
        xpro(ni2,nn), s(ni2,nn), trmin(nl), tmem(nl), &
        suo3(ni2,nn), del(ni2,nk) ,au(nl), az(nl), &
        nhaut(nl), nmoy(nl), nbas(nl)
   real uco2j(ni2,nkmax,2),tco2j(ni2,nkmax,nkmax), &
        wkco2j(ni2,nkmax,3)
 
   logical opnua

   real flsh(ni2,flnk), flss(ni2,flnn), fldel(ni2,flnk)
 
   !@Author l.garand (june 1989)
   !@Object
   !          to calculate infra-red radiation rate of cooling and flux
   !          according to garand jas 1983
   !@Arguments
   !          - output -
   ! taux     rate of cooling/warming in k/s
   ! fd       downward flux at the surface in w/m**2
   ! fm       upward flux at the surface in w/m**2
   !          - input -
   ! fnuage   cloud amount in each layer (0.0 - 1.0) times emissivity
   ! t        temperature in kelvin
   ! q        specific humidity (kg/kg)
   ! sh       sigma levels for t, q, oz
   ! ts       surface temperature in kelvins
   ! ps       surface pressure (n/m**2)
   ! nn       number of flux levels (nk+1)
   ! nk       number of layers
   ! oz       o3 mixing ratio in kg/kg
   ! nkmax    maximum number of levels (set in main program)
   ! ni2      maximum number of profiles to process
   ! nl       actual number of profiles to process
   ! ni1      1st dimension of t and q
   ! nt       dimension of table for temperature
   ! mx       number of values of u by decade in the table
   ! mxx      dimension of u in the table
   ! no3      dimension of ozone table: g1, g2, g3, th2o, bcn, dbcn, uu,
   !          yg3, bo3, dbo3, yg3o3, tro3, to2, uu: table precalculated
   !          in fg123
   ! ncx      dimension of co2 for yg3
   ! nco2     (parameter nco2); (pointer from subroutine pntg123)
   ! g1       5+1; (pointer from subroutine pntg123)
   ! g2       g1 + mxx*nt (pointer)
   ! g3       g2 + mxx*nt (pointer)
   ! th2o     g3 + mxx*nt (pointer)
   ! tro3     th2o + mxx*nco2 (pointer)
   ! yg3      tro3 + mxx (pointer)
   ! bcn      yg3 + nco2*mxx*ncx (pointer)
   ! dbcn     bcn + nco2*nt (pointer)
   ! bo3      dbcn+ nt*nco2; (pointer)
   ! dbo3     bo3 + nt (pointer)
   ! to3      dbo3 + nt (pointer)
   ! uu       to3 + no3 (pointer)
   ! f        work field: (2 parts)
   !          upper triangle contains cloud transmission
   !          lower triangle contains h2o-co2 transmission overlap
   ! ku       work field: (2 parts)
   !          upper triangle contains ozone table indices
   !          lower triangle contains h2o table indices
   ! opnua    .true. for cloud transmissivity to be .5 to .95
   !          (normal mode of operation, but unavailable in actual code)
   !          .false. for cloud transmissivity to be 100% (clear sky)
   !          (for research purposes)
   ! xin      work field
   ! v        work field
   ! su       work field
   ! itemp    work field
   ! xpro     work field
   ! suo3     work field  2d for local sigma
   ! au       work field
   ! reduc    .true. to use interpolation
   !          .false. means we are working on full levels
   ! radfix   .true. to use simple curve fits instead of full
   !          physics, from the lower stratosphere upwards
   ! radfltr  .true. to apply smoothing of net fluxes
   ! s        sigma flux levels (may be reduced)
   ! flsh     full sigma levels
   ! flss     full sigma flux levels
   ! del      thickness between flux levels (may be reduced)
   ! fldel    full ...
   ! flnk     full number of levels excluding ground
   ! flnn     full number of levels including ground
   !          - output -
   ! uco2j    amount of co2 in each layer in kg/m**2
   ! tco2j    transmission of co2
   !          - input -
   ! wkco2j   work field
   ! trmin    work field
   ! tmem     work field
   !          - output -
   ! az       work field; also used as output for zonal diagnostics
   !          of total clouds (variable nt)
   ! nhaut    integrated high-level clouds
   ! nmoy         "       mid-  "     "
   ! nbas         "       low-  "     "
   !@notes
   !          su and fm   can share the same space
   !          xpro and fd  "    "    "    "    "
   ! uco2     amount of co2 in each layer in kg/m**2
   ! tco2     precalculated transmission of co2 according to
   !          exp(-tco2*ps) which gives the matrix of co2 transmission.
   !          uco2 and tco2 are calculated in co2info.  upper triangle
   !          for co2 central band.  lower triangle for wings
   !          (average of the two wings)
   !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 
   !     limites de temperature
   parameter (t1=180., t2=320., to=260.)

   !     dependence en temperature de l'absorbant h2o
   parameter (ap=0.021)

   !     autres constantes physiques
   parameter ( po=101325.)
   parameter (elsa=1.66)


   parameter (epsil=1.e-12)
   parameter(corfac=600.)
   !     corfac modifie la temperature au premier niveau de flux
   !     pour le flux descendant pour tenir compte du toit trop bas
   !     en principe a modifier pour une discredisation choisie

   !     Beware of eo3 and eh2o! If values are changed, code
   !     must be modified as well because exponentiations have
   !     been eliminated
   parameter (eo3=1.00)
   parameter (eh2o=1.00)
   parameter (voigth2o=30., voigto3=400.)
   !     exposants pour calcul des masses d'ozone et d'h2o
   parameter (diftlim=8.)
   !     diftlim limite ts a + ou - diftlim de t(nn) dans radir
   !     ajoutant une  robustesse dans le cas de fortes erreurs de ts
   real tsair,ddd
   integer iv

   real stefinv



   integer idco2,ido3,idwat,ncen,nwng
   integer i_tempo
   real xlco2,xlo3,xlwat,aco2wng,aco2cen
   real rido3mx,ridwatmx
   real ln10inv, pogravinv

   real, pointer, dimension(:) :: auzt

   integer,    dimension(ni2)      :: ih, ib
   real,target,dimension(nl*2)     :: azzt
   real,       dimension(nl)       :: ps2pograv, z1d, z1d4, ts4
   real,       dimension(ni2,flnn) :: flfm, flfd
   real,       dimension(nl,nn)    :: sudelq, sudeloz
   !      REAL, dimension(ni2,flnn) :: flfm
   !      REAL, dimension(ni2,flnn) :: flfd
   !      INTEGER, dimension(ni2     ) :: ih
   !      INTEGER, dimension(ni2     ) :: ib
   !      real, dimension(nl*2    ) :: azzt
   !      real, dimension(nl      ) :: ps2pograv
   !      real, dimension(nl ,nn  ) :: sudelq
   !      real, dimension(nl ,nn  ) :: sudeloz
   !      real, dimension(nl      ) :: z1d
   !      real, dimension(nl      ) :: z1d4
   !      real, dimension(nl      ) :: ts4

   external intchamps, vco2inf2

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



!     azzt et auzt doivent etre contigus en memoire
!     (voir appel a vslog apres la boucle 111)
      auzt(1:nl) => azzt(nl+1:2*nl)

!     normalisation par stefan
      stefinv = 1./stefan

      ln10inv=1.0/log(10.)
      pogravinv=1.0/po/grav
      icen=2
      iwing=1
!     bande centrale =2 et aile=1 pour l'indice des tableaux co2
      trnuage=0.
      if(.not. opnua) trnuage=1.
!     transmissivite des nuages sera 100% (=clair) si opnua est false
!     typical clouds have 97% emissivity
!     clouds colder than 220k or above 100mb have 0.50 emissivity

!     matrix of co2 absorber amount starts at 10**-idco2
! similarly for ido3,idwat
! log10 of absorber amount starts at xlco2,xl03,xlwat
! caution idco2,ido3,idwat therefore needs to be changed if tables
! are changed in there beginning decade
! below are values for original tables (IRTAB4)
      idco2=5
      ido3=7
      idwat=6
      if(mxx.eq.501)then
! below are values for extended tables
      idco2=9
      ido3= 10
      idwat= 7
                     endif
      xlco2=-float(idco2) +1.e-6
      xlo3= -float(ido3) + 1.e-6
      xlwat=-float(idwat)+ 1.e-6
      borne=uu(mxx)-1.e-6
      rido3mx =(ido3*mx +1 +0.5)
      ridwatmx=(idwat*mx +1 +0.5)


      call vco2inf2(uco2j,tco2j,nl,nn,nk,nl,ni1,nkmax,sh,t,ps, &
           wkco2j(1,1,1),wkco2j(1,1,2),wkco2j(1,1,3),qco2)

      do 100 l=1,nl
!        Beware of eo3 and eh2o! If values are changed,
!        the following code must be modified.
!        suo3(l,1)=(s(l,1)+voigto3/ps(l))**eo3
         suo3(l,1)=s(l,1)+voigto3/ps(l)
         suo3(l,nn)=1.+voigto3/ps(l)
!        z=(s(l,1)+voigto3/ps(l))**eh2o
         z=s(l,1)+voigto3/ps(l)
         v(l,1)=t(l,1) - (sh(l,1)-s(l,1))/(sh(l,2)-sh(l,1))* &
                         (t(l,2)-t(l,1))


           if (TS_FLXIR) then

             tsair=t(l,nk+1)
             v(l,nn)=tsair
             su(l,1)= exp(ap*(v(l,1)-to)) *z
             su(l,nn)= exp(ap*(tsair-to))

           else

             v(l,nn)=ts(l)
!   discontinuite ts et tair non tenue en compte
!   mais ts limite par diftlim
!     v(l,nn)=amax1(ts(l),t(l,nk)-diftlim)
!     v(l,nn)=amin1(ts(l),t(l,nk)+diftlim)
             su(l,1)= exp(ap*(v(l,1)-to)) *z
             su(l,nn)= exp(ap*(ts(l)-to))

           endif



 100  continue

      do 20 i=2,nk
         do 120 l=1,nl
!        Beware of eo3 and eh2o! If values are changed,
!        the following code must be modified.
!           z=(s(l,i)+voigth2o/ps(l))**eh2o
!           suo3(l,i)=(s(l,i)+voigto3/ps(l))**eo3
            z=s(l,i)+voigth2o/ps(l)
            suo3(l,i)=s(l,i)+voigto3/ps(l)

            v(l,i)=(t(l,i)+t(l,i-1))*0.5
            su(l,i)= exp(ap*(v(l,i)-to)) *z
            sudelq (l,i)=(su(l,i)+su(l,i-1))*del(l,i-1)*0.5 *q(l,i-1)
            sudeloz(l,i)=(suo3(l,i)+suo3(l,i-1))*del(l,i-1)*0.5 *oz(l,i-1)
 120     continue

 20   continue
      do l=1,nl
         sudelq (l,nn)=(su(l,nn)+su(l,nn-1))*del(l,nn-1)*0.5 *q(l,nn-1)
         sudeloz(l,nn)=(suo3(l,nn)+suo3(l,nn-1))*del(l,nn-1)*0.5 *oz(l,nn-1)
         ps2pograv(l)=ps(l)*ps(l)*pogravinv*elsa
      enddo

      aa=(nt-1)/(t2-t1)
      bb=1.-aa*t1

      do 1 i=1,nn

         do 101 l=1,nl
            ku(l,i,i)=1
!           it=max0(1,nint(aa*v(l,i)+bb))
            i_tempo=aa*v(l,i)+bb+0.5
            itemp(l,i)=min(max(1,i_tempo),nt)
            f(l,i,i)=1.
            tmem(l)=1.
            trmin(l)=1.
            au(l)=0.
            az(l)=0.
!     maintenant v(l,i) remplace par v(l,i)**4
!           v(l,i)=v(l,i)**4
 101     continue
         call vspown1 (v(1,i),v(1,i),4.,nl)

         ip=i+1
         if(i.eq.nn)go to 40

         do 11 j=ip,nn
!           indx1=int(0.1 / (s(j-1)+epsil) ) + 1
            jind=j-2
            jind=max0(jind,1)

            do 111 l=1,nl

! transmission through cloud layer xnu
            xnu=1.-fnuage(l,j-1)
            if(fnuage(l,jind).lt.0.01)then
! entering a new cloud isolated from upper one
! keep in memory transmission down to top of new cloud tmem
! trmin is minimum transmission in cloud layer processed
! basic idea is random overlap (hence tmem X xnu = cloud transmittance)
! for cloud layers separated by clear ones; but maximum overlap
! (hence cloud transmittance is min (tmem,xnu))for adjacent cloud layers.
              tmem(l)= f(l,i,j-1)
              trmin(l)= xnu
            else
! inside a cloud use maximum overlap
! compute minimum transmission between adjacent layers
              trmin(l)=min(trmin(l),xnu)
            endif

            f(l,i,j)= tmem(l) * trmin(l)
!     elimine nuages si opnua false
               f(l,i,j)=max(f(l,i,j),trnuage)

               au(l)=au(l)+sudelq (l,j)
               az(l)=az(l)+sudeloz(l,j)
               azzt(l)=(az(l)*ps2pograv(l)+epsil)
               auzt(l)=(au(l)*ps2pograv(l)+epsil)
 111        continue
             call vslog(azzt,azzt,nl*2)
            do l=1,nl
               bz=azzt(l)*ln10inv
               bz=amax1(bz,xlo3)
               ku(l,i,j)=(mx*bz +rido3mx)
!               ku(l,i,j)=nint(mx*bz +itmp)
!     pour ozone indice 3 superieur a indice 2
               a=auzt(l)*ln10inv
               a=amax1(a,uu(1))
               a=amin1(a,borne)
               m=(mx*a +ridwatmx)
               ku(l,j,i)=m
!     pour transmission h2o-co2 indice 2 superieur a indice 3
               f(l,j,i)=tco2j(l,i,j) * th2o(m,icen)
             enddo

 11      continue

 1    continue

 40   continue
!  FLUX MONTANT *** FLUX UP

      call vspown1 (ts4,ts,4.,nl)

      DO 130 L=1,NL
!     FM(L,NN)=STEFAN*TS(L)**4
      FM(L,NN)=STEFAN*TS4(L)
      FD(L,1)=0.
!     IT=MAX0(1,NINT(AA*TS(L)+BB))
      I_TEMPO=AA*TS(L)+BB+0.5
      IT=MIN(MAX(1,I_TEMPO),NT)
      AZ(L)=FLOAT(IT)
 130  CONTINUE

      I1=1
      DO 2 I=I1,NK
!    CALCUL DE L'INTEGRALE PROCHE.

      DO 102 L=1,NL
      aco2wng=uco2j(l,i,1)
      aco2cen=uco2j(l,i,2)
      TX=T(L,I)
!     IT=MAX0(1,NINT(AA*TX+BB))
      I_TEMPO=AA*TX+BB+0.5
      IT=MIN(MAX(1,I_TEMPO),NT)
      WC=ALOG10(ACO2WNG+EPSIL)
      WC=AMAX1(WC,xlco2)
      WC=AMIN1(WC,.9999)
      WC2=ALOG10(ACO2CEN+EPSIL)
      WC2=AMAX1(WC2,xlco2)
      WC2=AMIN1(WC2,.9999)
      M=KU(L,I+1,I)
!     NWNG=NINT(MX*WC+ idco2*MX+1)
      NWNG=MX*WC+ idco2*MX+1.5
!     NCEN=NINT(MX*WC2+idco2*MX+1)
      NCEN=MX*WC2+idco2*MX+1.5
      A=DBCN(IT,ICEN) * YG3(M,NCEN,ICEN)
      Z=DBCN(IT,IWING)*YG3(M,NWNG,IWING)
      ZX=G3(M,IT)+A +DBO3(IT)*TO3(KU(L,I,I+1))*TRO3(M)+Z
      XIN(L,I)=(V(L,I)-V(L,I+1))*ZX*F(L,I,I+1)
      XPRO(L,I)=-XIN(L,I)
 102  CONTINUE


      IF(I.EQ.NK)GO TO 50

      IP=I+1
      DO 104 L=1,NL
      M=KU(L,IP,I)
      IT=ITEMP(L,IP)
      A=DBCN(IT,ICEN)
      Z=DBCN(IT,IWING)*TH2O(M,IWING) *TCO2J(L,IP,I)
      AU(L)=G2(M,IT)+A*F(L,IP,I)   +Z &
        +DBO3(IT)*TRO3(M)*TO3(KU(L,I,IP))
 104   CONTINUE

!    CALCUL DE L'INTEGRALE LOINTAINE.

      DO 3 J=IP,NK

      DO 103 L=1,NL
      R=AU(L)
      IT=ITEMP(L,J+1)
      M=KU(L,J+1,I)
      A=DBCN(IT,ICEN)
      Z=DBCN(IT,IWING)*TH2O(M,IWING)*TCO2J(L,J+1,I)
      AU(L)=G2(M,IT) +A *F(L,J+1,I) +Z &
        +DBO3(IT)*TRO3(M)*TO3(KU(L,I,J+1))
      XIN(L,I)=XIN(L,I)+.5*(R+AU(L))*(V(L,J)-V(L,J+1)) &
           *F(L,I,J+1)
 103  CONTINUE

  3   CONTINUE

  50  CONTINUE

      DO 140 L=1,NL

          if (TS_FLXIR) then

             IT=ITEMP(L,NN)
             M=KU(L,NN,I)
             IV=MAX(1,NINT(AA*TS(L)+BB))
             IV=MIN(IV,NT)

!  discontinuite a la surface entre TS**4 et V(NN)

             DDD=TS(L)**4 *( G1(M,IV) +BCN(IV,ICEN)*F(L,NN,I) &
             + BCN(IV,IWING) *TH2O(M,IWING) *TCO2J(L,I,NN) &
             + BO3(IV) *TRO3(M) *TO3(KU(L,I,NN)) ) &
             - V(L,NN)  *(G1(M,IT)  +BCN(IT,ICEN)*F(L,NN,I) &
             + BCN(IT,IWING) *TH2O(M,IWING) *TCO2J(L,I,NN) &
             + BO3(IT) *TRO3(M) *TO3(KU(L,I,NN)) )


             FM(L,I)=STEFAN* V(L,I) - XIN(L,I) + DDD

          else

             FM(L,I)=STEFAN* V(L,I) - XIN(L,I)

          endif


 140  CONTINUE

  2   CONTINUE


! FLUX DESCENDANT

      if (radfix) then
!   CORRECTION TEMPERATURE AU TOIT
      do L=1,NL
      Z1D(L)= T(L,1)+CORFAC*S(L,1)
      end do
      call vspown1(z1d4,z1d,4.,nl)

      DO 125 L=1,NL
      I_TEMPO=AA*Z1D(L)+BB+0.5
      IT=MIN(MAX(1,I_TEMPO),NT)
      ITEMP(L,1)=IT
      xpro(l,1)=xpro(l,1)*(z1d4(L)-v(l,2))/ &
                sign(max(abs(v(l,1)-v(l,2)),epsil),v(l,1)-v(l,2))
!     XPRO(L,1)=XPRO(L,1)*(Z**4-V(L,2))/(V(L,1)-V(L,2)+EPSIL)
      V(L,1)=Z1D4(L)
 125  CONTINUE
      endif


      I2=2
      DO 55 I=I2,NN

      DO 105 L=1,NL
      XIN(L,I)=-XPRO(L,I-1)
 105  CONTINUE
  55   CONTINUE

      DO 5 I=I2,NN
      IF(I.EQ.I2)GO TO 60

      DO 106 L=1,NL
      IT=ITEMP(L,1)
      M=KU(L,I,1)
      A=DBCN(IT,ICEN)
      Z=DBCN(IT,IWING)*TH2O(M,IWING)  *TCO2J(L,I,1)
      AU(L)=G2(M,IT)+A*F(L,I,1)  +Z &
        +DBO3(IT)* TRO3(M) *TO3(KU(L,1,I))
 106  CONTINUE

      DO 9 J=1,I-2

      DO 109 L=1,NL
      R=AU(L)
      IT=ITEMP(L,J+1)
      M=KU(L,I,J+1)
      A=DBCN(IT,ICEN)
      z=0.
      if(tco2j(l,i,j+1).gt.1.e-4)then
      Z=DBCN(IT,IWING)*TH2O(M,IWING)*TCO2J(L,I,J+1)
       endif
      AU(L)=G2(M,IT) +A*F(L,I,J+1) +Z &
       +DBO3(IT) *TRO3(M) * TO3(KU(L,J+1,I))
      XIN(L,I)=XIN(L,I)+.5*(R+AU(L))*(V(L,J)-V(L,J+1)) &
       *F(L,J+1,I)
 109  CONTINUE

   9  CONTINUE

  60  CONTINUE

!    CALCUL DU TERME DE REFROIDISSEMENT VERS L'ESPACE.

      DO 160 L=1,NL
      IT=ITEMP(L,1)
      M=KU(L,I,1)
      A=BCN(IT,ICEN)
      Z=BCN(IT,IWING)*TH2O(M,IWING)  *TCO2J(L,I,1)
      Y=G1(M,IT) +A*F(L,I,1) + Z &
       +BO3(IT)*TRO3(M)*TO3(KU(L,1,I))
      FD(L,I)=XIN(L,I)+ STEFAN*V(L,I) -V(L,1) *Y *F(L,1,I)
 160  CONTINUE

   5  CONTINUE


!cc   if(radfix) then
!     ajustement au toit fd(1) est extrapole et non zero
      do 175 l=1,nl
         fd(l,1)=amax1(0.,fd(l,2)-(fd(l,3)-fd(l,2))*del(l,1)/del(l,2))
         fd(l,1)=amin1(fd(l,1),fd(l,2))
 175  continue
!cc   else
!cc   do  l=1,nl
!     Utiliser la moyenne des deux approches - BD, 26 novembre 1993.
!cc      fd(l,1)=amax1( 0.,
!cc  +                  fd(l,2)-( fd(l,3)-fd(l,2) )*del(l,1)/del(l,2)
!cc  +                )
!cc  +          /                         2.0
!cc      fd(l,1)=amin1(fd(l,1),fd(l,2))
!     Ou bien, mettre le flux descendant a zero - BD, 4 juin 1997.
!cc      fd(l,1)=0.
!cc   enddo
!cc   endif

! recherche du jump dans transmission h2o (indice de tableau =1)
! elimination par interpolation du flux vers le bas a ce niveau
! ou les indices de tableau passent de 1 a > 1.
! ceci se produit entre 5 et 10 mb.
         do 671 i=3,nk
         do 672 l=1,nl
         if(ku(l,i,1) .gt.1. .and. ku(l,i-1,1).eq.1) &
          fd(l,i)=fd(l,i-1)+(fd(l,i+1)-fd(l,i-1))/(s(l,i+1)-s(l,i-1))* &
                  (s(l,i)-s(l,i-1))
 672     continue
 671     continue



!     calcul des indices IH et IB pour nuages 2-D
!     IH = niveau le plus pres de sigma=0.4
!     IB = niveau le plus pres de sigma=0.7
      do j=1,nk
         do l=1,nl
            if (sh(l,j).le.0.4) ih(l) = j
            if (sh(l,j).le.0.7) ib(l) = j
         end do
      end do

      do 208 l=1,nl
!     nuages totaux
         az(l)=1.-f(l,1,nn)

!        diagnostics de nuages 2-D
         nhaut(l) = 1. - f(l, 1   ,IH(l))
         nmoy (l) = 1. - f(l,IH(l),IB(l))
         nbas (l) = 1. - f(l,IB(l),nn   )

 208  continue


!     En mode reduction, interpoler fd et fm aux niveaux de flux complets

      if( reduc ) then

        call intchamps(flfm,fm,flss,s,nl,flnn,nn)
        call intchamps(flfd,fd,flss,s,nl,flnn,nn)

        do i=1,flnn
          do l=1,nl
            fm(l,i) = flfm(l,i)
            fd(l,i) = flfd(l,i)
          enddo
        enddo
      endif

      do 2005 i=1,flnn
         do 2006 l=1,nl
            xin(l,i)=fm(l,i)
            v(l,i)=fd(l,i)
 2006    continue
 2005 continue

!     smoothing 1/4-1/2-1/4 des flux
!************************************

      if ( radfltr ) then

           do  i=2,flnk
           do  l=1,nl
             fm(l,i)=0.25*(xin(l,i-1)+xin(l,i+1)) + 0.5*xin(l,i)
             fd(l,i)=0.25*(v(l,i-1)+v(l,i+1)) + 0.5*v(l,i)
           enddo
           enddo

      else

           do i=2,flnk
           do l=1,nl
              fm(l,i)=xin(l,i)
              fd(l,i)=v(l,i)
           enddo
           enddo

      endif
!**************************************

      secday=1./86400.

!     taux de refroidissement
      do 8 i=1,flnk

!        taumax = (( -877.193 * (flsh(i)-0.05)**2 ) - 0.5)*secday

         if(radfix) then
! >>>>>> Seulement si RADFIX est vrai <<<<<<<
         do 108 l=1,nl
!           taumax = (( -877.193 * (flsh(l,i)*ps(l)*1.E-5-0.05)**2) - 0.5)*secday
            taumax= (( -877.193 *  (flsh(l,i)*ps(l)*1.E-5-0.05)* &
                        (flsh(l,i)*ps(l)*1.E-5-0.05))  - 0.5)*secday

            taux(l,i)=-(fm(l,i)-fd(l,i)-fm(l,i+1)+fd(l,i+1)) / &
                 (ps(l)*cpd*fldel(l,i)) *grav
!     fit du taux au dessus de sigma=0.011 en fonction de temperature

!     legere correction des taux au-dessus de 100 mb de facon a
!     obtenir un equilibre quasi-parfait avec les taux visibles
!     de sun1.  correction inferieure a 0.15 k/j. le taux au toit
!     est fixe a -1.9 k/jour.
      r = flsh(l,i)*ps(l)*0.01
      if (r.lt.100.) taux(l,i)= taux(l,i)+0.15*secday

      if (i.eq.1) taux(l,i)= -1.9*secday
!     taux(l,i)=cvmgt((-1.9*secday),taux(l,i),(i.eq.1))
!     taux(l,i)=cvmgt((-1.9*secday),taux(l,i),flsh(l,i).lt.0.011)

!     limite sur les taux (en general trop froids) lorsque sigma <= 0.05
      if ((flsh(l,i)*ps(l)*1.E-5.le.0.05).and.(i.ne.1)) then
        taux(l,i) = max(taumax,taux(l,i))
! >>>>>> Seulement si RADFIX est vrai <<<<<<<
      endif
 108     continue

         else

!          calcul des taux radiatifs sans corrections
           do l=1,nl
              taux(l,i)=-(fm(l,i)-fd(l,i)-fm(l,i+1)+fd(l,i+1)) / &
                   (ps(l)*cpd*fldel(l,i)) *grav
           enddo


        endif
 8    continue

!     aux deux bouts taux calcules avec flux sans smoothing
      do 1088 l=1,nl
         y=-(xin(l,1)-v(l,1)-xin(l,2)+v(l,2)) / &
              (ps(l)*cpd*fldel(l,1)) *grav
!        si radfix est vrai, le taux est fixe par fit en temperature
         if ((flsh(l,1).ge.0.011) .and. .not.radfix ) taux(l,1) = y
         taux(l,flnk)=-(xin(l,flnk)-v(l,flnk)-xin(l,flnn)+v(l,flnn)) / &
              (ps(l)*cpd*fldel(l,flnk)) *grav

 1088 continue


      return
      end
