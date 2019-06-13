!-------------------------------------- LICENCE BEGIN ------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html

!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END --------------------------

subroutine vco2inf2(uco2,tco2,nl,nn,nk,nls,ni1,nmax, &
     sh,t,ps,s,sc,del,co2ppmv)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>

   real a1d,a1g,a2d,a2g,awing
   integer nl,nn,nk,nls,nmax,ni1
   real eco2,a1c,a2c,qco2,elsa,z,co2ppmv
   real uco2(nls,nmax,2),del(nls,nk),sc(nls,nmax),tco2(nls,nmax,nmax)
   real s(nls,nmax),sh(nls,nk),t(ni1,nk),ps(nls)

   !@author l.garand (1989)
   !@revision
   ! 001      g.pellerin(mar90)standard documentation
   ! 002      louis garand -add co2 wing bands
   ! 003      Y. Chartier (dec93) add first dimension nls to compute
   !          transmissivity in local sigma
   ! 004      l. garand (march 94) add temperature effects on transmission
   !                     these are important above 50 mb
   ! 005      l. garand (april 96) transition from Lorentz to Voigt line shape
   !                     following Giorgetta & Morcrette, MWR, 1995, p. 3381-3383
   ! 006      l. garand (november 97) - change CO2 concentration from 330 ppm
   !                                    to 360 ppm
   ! 007      b. dugas (sep 2002)     - CO2 ppmv concentration is passed
   !                                    as input parametre co2ppmv
   ! 008      b. bilodeau (april 2003) - IBM conversion
   !             - calls to vsexp  routine (from massvp4 library)
   !             - calls to vssqrt routine (from massvp4 library)
   !             - removal of loop 120 on k index; code now has a cost that
   !               is proportional to nk**2 instead of nk**3
   !             - removal of useless exponentiations
   ! 009      m. desgagne and m. valin (april 2005) - optimization for OpenMP
   !@object
   !          to precalculate the quantities of co2 and the
   !          transmissivity from level to level
   !@arguments
   !          - output -
   ! uco2     amount of co2 in each layer of thickness del (kg/m2)
   !          third index: 1: wings, 2: center
   ! tco2     precalculated transmissivity of co2 from level to level
   !          upper triangle of tco2 is used for the (strong) central band.
   !          the lower triangle of tco2 is used for the average of the
   !          right and left wings)
   !          - input -
   ! nl       number of profiles to process
   ! nn       number of levels (nk+1)
   ! nk       number of layers
   ! nls      1st dimension of uco2 and vco2
   ! ni1      1st dimension of t
   ! nmax     maximum number of flux layers
   ! sh       sigma levels at the centre of layers
   ! t        temperature (k)
   ! ps       surface pressure (newton/m2) for each profile
   ! s        sigma levels at the borders of the layers
   ! sc       work space
   ! del      sigma thickness from level to level
   ! ps       surface pressure (N/m2)
   ! co2ppmv  co2 concentration in ppmv

   real aprimec,aprimew
   integer k,kk,kkk,k2,l
   real voigt
   ! all these parameters in table 1 of internal publication
   ! by garand and mailhot (1990)
   !     Beware of eco2! See comment in loop 60
   parameter (eco2=1.00)
   parameter (a1c=198.0)
   parameter (a2c=0.9768)
   parameter (a1d=4.035)
   parameter (a2d=0.8224)
   parameter (a1g=5.439)
   parameter (a2g=0.9239)
   parameter (voigt = 60.)
   ! temperature coefficienta a' for center and wings in k-1
   ! especially important for wings; b' factor neglected
   parameter (aprimec= 3.1e-3)
   parameter (aprimew= 15.8e-3)

   real, dimension(NL     ) :: TRAPEZ1
   real, dimension(NL     ) :: TRAPEZ2
   real, dimension(NL     ) :: XT
   real, dimension(NL     ) :: XX
   real(REAL64), dimension(NL     ) :: XP
   real(REAL64), dimension(NL,NN  ) :: TRAPEZE1
   real(REAL64), dimension(NL,NN  ) :: TRAPEZE2
   real, dimension(NL,NN  ) :: XTK
   real, dimension(NL,NN  ) :: XXK

   !***********************************************************************


      qco2 = nint( 547. * ( co2ppmv / 360.d0 ) ) * 1.e-6

!  a1d et a2d sont les parametres de l'aile droite du co2
!  a1g et a2g """"""""""""""""""""""""""""" gauche """""""
!  a1c et a2c """"""""""""""""""" de la bande centrale (forte) du co2
!*


      awing= (a1g*a2g + a1d*a2d)/2.
!     parametre d'absortion moyen pour les deux ailes

      do 10 l=1,nl
         sc(L,nn)=qco2*(1.+ voigt/ps(L))
         s(L,nn)=1.
         s(L,1)=2.*sh(L,1)-((sh(L,1)+sh(L,2))*0.5)
!     cette definition du premier niveau de flux doit etre la meme
!     que dans le code de radiation
         s(L,1)=amax1(s(L,1),sh(L,1)/2.)
!        s(L,1)=amax1(s(L,1),0.0003)
         s(L,nn)=1.
 10   continue

      do 20 k=2,nk
         do 30 l=1,nl
            s(L,k)=(sh(L,k)+sh(L,k-1))*0.5
            del(L,k-1)=s(L,k)-s(L,k-1)
 30      continue
 20   continue

      do k=1,nk
         do l=1,nl
            xtk(L,k)=(t(L,k)-260.)*aprimec
            xxk(L,k)=(t(L,k)-260.)*aprimew
         enddo
      enddo


      call vsexp(xtk,xtk,nl*nk)
      call vsexp(xxk,xxk,nl*nk)

      do 40 l=1,nl
         del(L,nk)=1.-s(L,nk)
         uco2(L,nn,1)=0.
         uco2(L,nn,2)=0.
 40      continue


      do 50 k=1,nk
         do 60 l=1,nl
!           sc(L,k)=qco2*(s(L,k)+voigt/ps(L))**eco2
!           Beware! Attention! The following line of code
!           is valid only if eco2 = 1.
            sc(L,k)=qco2*(s(L,k)+voigt/ps(L))
            trapeze1(L,k) = 0.0d0
            trapeze2(L,k) = 0.0d0
 60      continue
 50   continue

      do l=1,nl
         trapeze2(l,1) = trapeze2(l,1) + &
                         xtk(l,1)*(sc(l,1)+sc(l,2))*del(l,1)
         trapeze1(l,1) = trapeze1(l,1) + &
                         xxk(l,1)*(sc(l,1)+sc(l,2))*del(l,1)
      end do

      do k=2,nk
         do l=1,nl
            trapeze2(l,k) = trapeze2(l,k-1) + &
                            xtk(l,k)*(sc(l,k)+sc(l,k+1))*del(l,k)
            trapeze1(l,k) = trapeze1(l,k-1) + &
                            xxk(l,k)*(sc(l,k)+sc(l,k+1))*del(l,k)
         enddo
      enddo

      elsa=1.66

      z=1./(101325.*9.80616)
      do l=1,nl
       xp(L)=(z*elsa*0.5)*ps(L)*ps(L)
      enddo
      do 70 k=1,nn
         do 80 l=1,nl
            tco2(L,k,k)=1.
 80      continue
 70   continue

      do 90 k=1,nn-1

         kk=k+1
         do 100 k2=kk,nn
            kkk=k2-k+1

            if (k.eq.1) then
               do l=1,nl
                  trapez1(L) = max(trapeze1(L,k2-1)*xp(L),0.0d0)
                  trapez2(L) = max(trapeze2(L,k2-1)*xp(L),0.0d0)
               end do
            else
               do l=1,nl
                  trapez1(L) = max((trapeze1(L,k2-1)-trapeze1(L,k-1))*xp(L),0.0d0)
                  trapez2(L) = max((trapeze2(L,k2-1)-trapeze2(L,k-1))*xp(L),0.0d0)
               end do
            endif

            if (k2.eq.kk) then
               do 150 l=1,nl
                  uco2(L,k,1)=trapez1(L)
                  uco2(L,k,2)=trapez2(L)
 150           continue
            endif

            do 160 l=1,nl
               xx(L)=a1c*a2c*trapez2(L)
               xt(L)=awing*trapez1(L)
 160        continue

            call vssqrt(xx,xx,nl)
            call vssqrt(xt,xt,nl)
            xx=-xx
            xt=-xt
            call vsexp(xx,xx,nl)
            call vsexp(xt,xt,nl)
            do l=1,nl
               tco2(L,k,k2)= xx(L)
               tco2(L,k2,k)= xt(L)
            end do

 100     continue
 90   continue

      return
      end
