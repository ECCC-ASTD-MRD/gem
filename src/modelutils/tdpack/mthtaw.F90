!-------------------------------------- LICENCE BEGIN ------------------------------------
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
!-------------------------------------- LICENCE END --------------------------------------
!**s/r mthtaw3  -  calcule tw ou thetaw
!
      Subroutine mthtaw3(tw,hu,tt,ps,swph,swth,ni,nk,n)
      use tdpack 
      implicit none
!!!#include <arch_specific.hf>
!
      Integer ni, nk, n
      Real tw(ni,nk), hu(ni,nk), tt(ni,nk)
      Real ps(ni,*)
      Logical swph, swth
!
!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      N. Brunet (May 1994) - change to l2ocprv in order to
!          obtain results closer to those of tephigram
! 003      B. Bilodeau (Jan 2001) - Automatic arrays
! 004      M. Lepine (March 2003) -  CVMG... Replacements
!
!Object
!          to calculate TW or THETAW (according to the value SWTT)
!          from specific humidity, temperature and pressure
!
!Arguments
!
!          - Output -
! tw       tw or thetaw in K
!
!          - input -
! hu       specific humidity in kg/kg
! tt       temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! swth     .true. to calculate thetaw
!          .false. to calculate tw
! ni       horizontal dimension
! nk       vertical dimension
! n        number of treated points
!*
!--------------------------------------------------------------------
      Real q1, dq1, q0, th, h, ft0, dft0, prn
      Real dlp, d, l2ocprv, qp, fac
      Integer iter, i, k, nn, j
!--------------------------------------------------------------------
!
      l2ocprv = 1.35e+07
!
      Do k=1,nk
         Do i=1,n
!
!     on trouve 'tw' en solutionnant par methode
!     de newton du 1er ordre avec 4 iterations.
!
            q0      = hu(i,k)
            tw(i,k) = tt(i,k)
!
            Do iter=1,4
               h = htvocp(tw(i,k))
               If(swph)Then
                  q1 = foqst(tw(i,k),ps(i,k))
                  dq1 = fodqs(q1,tw(i,k))
               Else
                  q1 = foqsa(tw(i,k),ps(i,k))
                  dq1 = fodqa(q1,tw(i,k))
               End If
               ft0 = -h*(q1 - q0)
               dft0 = -1.0 -(h*dq1)
               tw(i,k) = tw(i,k) - ft0/dft0
               q0 = q0 + ft0/(dft0*h)
            Enddo
!
            If(swth)Then
!
!     on trouve thetaw en utilisant la pente de l'adiabatique
!     mouillee et en faisant les calculs par tranche de
!     1000 pa ('dlp') en passant de 'ps' a 100000.
!
               If(ps(i,k).Ne.100000.)Then
                  dlp = 1000.
                  d = 100000. - ps(i,k)
                  nn = Int(Abs(d/dlp))
                  If(d.Lt.0.)dlp = -dlp
                  d = d - float(nn)*dlp
                  If(d.Ne.0.)Then
                     nn = nn + 1
                  Else
                     d = dlp
                  End If
!
                  prn = ps(i,k)
                  th = tw(i,k)
!
                  Do j=1,nn
                     qp = d/prn
                     q1 = foqst(th,prn)
                     If(.Not.swph)q1 = foqsa(th,prn)
                     h = htvocp(th)
                     fac = (cappa + (h*q1/th))/(1. + &
                     l2ocprv*q1/(th*th))
                     th = th*(1. + qp*fac)
!
                     prn = prn + d
                     If(j.Eq.1)d = dlp
                  Enddo
!
                  tw(i,k) = th
!
               End If
!
            End If
         Enddo
      Enddo
!
      Contains
!     --------------------------------------------------------
!     fonction calculant  h/cp
!     ttp= temperature en k
!     ai,bi,aw,bw,t1s,t2s,slp = voir common /ctesdyn/
!
!      htvocp(ttp) = cvmgt((ai-bi*ttp),cvmgt((aw-bw*ttp),
!     x             slp*((aw-bw*ttp)*(ttp-t2s)+(ai-bi*ttp)*(t1s-ttp))
!     y             ,ttp.ge.t1s),ttp.le.t2s)
!     -----------------------------------------------------------
!
!  internal function definition
!
      Real Function htvocp(ttp)
      Real ttp

      If      (ttp.Le.t2s) Then
         htvocp = ai-bi*ttp
      Else If (ttp.Ge.t1s) Then
            htvocp = aw-bw*ttp
      Else
            htvocp = slp*((aw-bw*ttp)*(ttp-t2s)+(ai-bi*ttp)*(t1s-ttp))
      Endif
      End Function htvocp

      End Subroutine mthtaw3


!**s/r mthtaw4  -  calcule tw ou thetaw
!
      Subroutine mthtaw4(tw,hu,tt,ps,swph,swth,ti,ni,nk,n)
      use tdpack
      implicit none
!!!#include <arch_specific.hf>
      Integer ni, nk, n
      Real tw(ni,nk), hu(ni,nk), tt(ni,nk)
      Real ps(ni,*), ti
      Real temp1
      Logical swph, swth
!
!Author
!          N. Brunet  (Jan91)
!
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      N. Brunet (May 1994) - change to l2ocprv in order to
!          obtain results closer to those of tephigram
! 003      N. Brunet (sept 2000) adaptation to new functions
!          sgamasp and schal.
!
!Object
!          to calculate tw or thetaw (according to the value swtt)
!          from specific humidity, temperature and pressure
!
!Arguments
!
!          - output -
! tw       tw or thetaw in K
!
!          - input -
! hu       specific humidity in kg/kg
! ti       temperature in K at which we start calculating
!            latent heat of sublimation
!            if swph=false, ti is n/a
!            ti must be .le. trpl
!
! tt       temperature in K
! ps       perssure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! swth     .true. to calculate thetaw
!          .false. to calculate tw
! ni       horizontal dimension
! nk       vertical dimension
! n        number of treated points
!*
!--------------------------------------------------------------------
      Real q1, dq1, q0, th, ft0, dft0, prn
      Real dlp, d, dt, dtpr, latheat, hscp
      Integer iter, i, k, nn, j
!--------------------------------------------------------------------
!
      Do k=1,nk
         Do i=1,n
!
!     on trouve 'tw' en solutionnant par methode
!     de newton du 1er ordre avec 4 iterations.
!
            q0 = hu(i,k)
            tw(i,k) = tt(i,k)
!
            Do iter=1,4
               latheat = schal(tw(i,k), ti, swph)
               hscp = latheat / cpd
               temp1 = ps(i,k)
               If(swph)Then
                  q1 = foqst(tw(i,k),temp1)
                  dq1 = fodqs(q1,tw(i,k))
               Else
                  q1 = foqsa(tw(i,k),temp1)
                  dq1 = fodqa(q1,tw(i,k))
               End If
               ft0 = -hscp*(q1 - q0)
               dft0 = -1.0 -(hscp*dq1)
               tw(i,k) = tw(i,k) - ft0/dft0
               q0 = q0 + ft0/(dft0*hscp)
            Enddo
!
            If(swth)Then
!
!     on trouve thetaw en utilisant la pente de l'adiabatique
!     mouillee et en faisant les calculs par tranche de
!     1000 Pa ('dlp') en passant de 'ps' a 100000 Pa.
!
               If(ps(i,k).Ne.100000.)Then
                  dlp = 1000.
                  d = 100000. - ps(i,k)
                  nn = Int(Abs(d/dlp))
                  If(d.Lt.0.)dlp = -dlp
                  d = d - float(nn)*dlp
                  If(d.Ne.0.)Then
                     nn = nn + 1
                  Else
                     d = dlp
                  End If
!
                  prn = ps(i,k)
                  th = tw(i,k)
!
                  Do j=1,nn
!                    --- calcule dt / dp
                     dt = sgamasp(th, prn, swph, ti)
!                    --- multiplie par delta p --> donne delta t
                     dtpr = dt * d
!                    --- mise a jour de la temp et de la pression
                     th = th + dtpr
                     prn = prn + d
                     If(j.Eq.1)d = dlp
                  Enddo
!
                  tw(i,k) = th
!
               End If
!
            End If
         Enddo
      Enddo
!
      End Subroutine mthtaw4
