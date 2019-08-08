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
!**s/r fonction sthtaw3  -  calcule tw ou thetaw
!
      Function sthtaw3(hu,tt,ps,swph,swth)
      use tdpack, only: foqst, foqsa, fodqs, fodqa, cappa, ai, aw, bi, bw, cpd, slp, t1s, t2s
      implicit none
!!!#include <arch_specific.hf>
      Real sthtaw3, hu, tt, ps
      Logical swph, swth
!
!author
!          n. brunet  (jan91)
!revision
! 001      b. bilodeau  (august 1991)- adaptation to unix
! 002      n. brunet  (may 1994) - change to l2ocprv in order
!          to obtain results closer to those of tephigram
! 003      m. lepine (march 2003) -  cvmg... replacements
!
!object
!          to return tw or thetaw by calculating from specific
!          humidity, temperature and pressure. result returned is in
!          kelvins
!
!arguments
!
!          - input -
! hu       specific humidity in kg/kg
! tx       temperature in K
! ps       pressure in Pa
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! swth     .true. to calculate theta
!          .false. to calculate tw
!*
!--------------------------------------------------------------------
      Real q1, dq1, q0, th, h, ft0, dft0
      Real d, dlp, qp, pr, fac, l2ocprv
      Integer iter, n, j
!--------------------------------------------------------------------
!     trouve d'abord  tw
!     solutionne par methode de newton du 1er ordre.
!
      th = tt
      q0 = hu
      Do iter=1,4
         h = htvocp(th)
            q1=foqst(th,ps)
            If(.Not.swph)q1=foqsa(th,ps)
            dq1=fodqs(q1,th)
            If(.Not.swph)dq1=fodqa(q1,th)
         ft0 = -h * (q1-q0)
         dft0 = -1.0 - (h*dq1)
         th = th - ft0/dft0
         q0 = q0 + ft0/(dft0*h)
      Enddo
!
      If(.Not.swth)Then
         sthtaw3 = th
         Return
      End If
!
!     trouve thetaw en utilisant la pente de l'adiabatique
!     mouillee et en faisant les calculs par tranche de
!     1000 pa (passant de pn a 100000 pa).
!
      If(ps.Eq.100000.)Then
         sthtaw3 = th
      Else
         dlp=1000.
         l2ocprv = 1.35e+07
         d = 100000. - ps
         n = Int(Abs(d/dlp))
         If(d.Lt.0.)dlp = -dlp
         d = d - float(n)*dlp
         If(d.Ne.0.)Then
            n = n + 1
         Else
            d = dlp
         End If
!
         pr = ps
!
         Do j=1,n
            qp = d/pr
            q1 = foqst(th,pr)
            If(.Not.swph)q1 = foqsa(th,pr)
            h = htvocp(th)
            fac = (cappa + (h*q1/th))/(1. + l2ocprv*q1/(th*th))
            th = th*(1. + qp*fac)
            pr = pr + d
            If(j.Eq.1)d = dlp
         Enddo
!
         sthtaw3 = th
      End If
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

      End Function sthtaw3


!**s/r fonction sthtaw4  -  calcule tw ou thetaw
!
      Function sthtaw4(hu,tt,ps,swph,swth,ti)
      use tdpack, only: schal, foqst, foqsa, fodqs, fodqa, sgamasp, cappa, ai, aw, bi, bw, cpd, slp, t1s, t2s
      implicit none
!!!#include <arch_specific.hf>
      Real sthtaw4, hu, tt, ps, ti
      Logical swph, swth
!
!Author
!          N. Brunet  (Jan91)
!Revision
! 001      B. Bilodeau  (August 1991)- Adaptation to UNIX
! 002      N. Brunet  (May 1994) - change to l2ocprv in order
!          to obtain results closer to those of tephigram
! 003      N. Brunet (sept 2000) adaptation to new functions
!          sgamasp and schal.
!
!Object
!          to compute tw or thetaw, function of specific
!          humidity, temperature and pressure. Result returned is in K
!
!Arguments
!
!          - Input -
! hu       specific humidity in kg/kg
! tt       temperature in K
! ps       pressure in K
! swph     .true. to consider water and ice phase
!          .false. to consider water phase only
! swth     .true. to calculate theta
!          .false. to calculate tw
! ti       temperature (k) at which we start considering only
!          latent heat of sublimation
!          if swph=false, ti is n/a
!          ti must be .le. trpl
!*
!--------------------------------------------------------------------
      Real q1, dq1, q0, th, ft0, dft0
      Real latheat, dt, dtpr
      Real hscp
      Real d, dlp, pr
      Integer iter, n, j
!--------------------------------------------------------------------
!
!
!     trouve d'abord  tw
!
      q0 = hu
      th = tt
!
!     solutionne par methode de newton du 1er ordre.
!
      Do iter=1,4
!        --- calcule la chaleur latente
         latheat = schal(th, ti, swph)
!
!        --- methode de newton
!
         hscp = latheat / cpd
!
         q1=foqst(th,ps)
         If(.Not.swph)q1=foqsa(th,ps)
         dq1=fodqs(q1,th)
         If(.Not.swph)dq1=fodqa(q1,th)
         ft0 = -hscp * (q1-q0)
         dft0 = -1.0 - (hscp*dq1)
         th = th - ft0/dft0
         q0 = q0 + ft0/(dft0*hscp)
      Enddo
!
      If(.Not.swth)Then
         sthtaw4 = th
         Return
      End If
!
!     trouve thetaw en utilisant la pente de l'adiabatique
!     mouillee et en faisant les calculs par tranche de
!     1000 Pa (passant de ps a 100000 Pa).
!
      If(ps.Eq.100000.)Then
         sthtaw4 = th
      Else
         dlp=1000.
         d = 100000. - ps
         n = Int(Abs(d/dlp))
         If(d.Lt.0.)dlp = -dlp
         d = d - float(n)*dlp
         If(d.Ne.0.)Then
            n = n + 1
         Else
            d = dlp
         End If
!
         pr = ps
!
         Do j=1,n
!
!           --- calcule dt/dp
            dt = sgamasp(th, pr, swph, ti)
!           --- multiplie par delta p --> donne delta t
            dtpr = dt * d
!           --- mise a jour de la temp et de la pression
            th = th + dtpr
            pr = pr + d
!
            If(j.Eq.1)d = dlp
!
         Enddo
!
         sthtaw4 = th
!
      End If
!
      End Function sthtaw4
