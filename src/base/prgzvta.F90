!---------------------------------- LICENCE BEGIN -------------------------------
! GEM - Library of kernel routines for the GEM numerical atmospheric model
! Copyright (C) 1990-2010 - Division de Recherche en Prevision Numerique
!                       Environnement Canada
! This library is free software; you can redistribute it and/or modify it
! under the terms of the GNU Lesser General Public License as published by
! the Free Software Foundation, version 2.1 of the License. This library is
! distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
! PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
! You should have received a copy of the GNU Lesser General Public License
! along with this library; if not, write to the Free Software Foundation, Inc.,
! 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
!---------------------------------- LICENCE END ---------------------------------

!**s/r prgzvta - interpolation of geopotential and virtual temperature on
!                given pressure levels

      subroutine prgzvta(F_gzout, F_vtout,  F_pres, Nkout, &
                         F_gzin,  F_vtin,   F_wlnph, F_la, &
                         F_vtund, F_gzund,  F_nundr, &
                         F_cubzt_L, F_linbot, &
                         Minx,Maxx,Miny,Maxy, F_Nk)
      use tdpack
      use glb_ld
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical F_cubzt_L
      integer F_nundr, F_linbot
      integer Minx,Maxx,Miny,Maxy,F_Nk,Nkout
      real    F_pres(Nkout)
      real    F_gzout(Minx:Maxx,Miny:Maxy,Nkout),F_vtout(Minx:Maxx,Miny:Maxy,Nkout)
      real    F_gzin (Minx:Maxx,Miny:Maxy,F_Nk), F_vtin (Minx:Maxx,Miny:Maxy,F_Nk)
      real    F_wlnph (Minx:Maxx,Miny:Maxy,F_Nk), F_la   (Minx:Maxx,Miny:Maxy)
      real    F_vtund(Minx:Maxx,Miny:Maxy,F_nundr), F_gzund(Minx:Maxx,Miny:Maxy,F_nundr)


      integer i, j, k, kk, pnk, pnkm, pnindex(l_ni), pnund,   pn1
      real    prlprso
      real    prd, pre, prr
      real    prfm0, prfm1, prfm2, prfm3, prfl2
      real    prl, prsmall
      real    prlptop, prvttop, prfitop
      real    prlpbot, prvtbot, prfibot
      real    logpres(Nkout)
      real(kind=REAL64)  invprd
!
!-------------------------------------------------------------------
!
      prsmall= .001
      logpres(1:Nkout) = log(F_pres(1:Nkout))

      do j  = 1, l_nj
         do kk = 1, Nkout
            pnindex = 0
            prlprso = logpres(kk)

            do k= 1, F_nk
               do i= 1, l_ni
                  if ( prlprso > F_wlnph(i,j,k) ) pnindex(i) = k
               end do
            end do

            do i= 1, l_ni
!******************************************************************************
!                                                                             *
! If:    output pressure   <   hydrostatic pressure on the                    *
!                              first level of the model                       *
!                                                                             *
! Then:  upward extrapolation                                                 *
!                                                                             *
!******************************************************************************

               if ( pnindex(i) == 0 ) then
                  prd = prlprso - F_wlnph(i,j,1)
                  F_vtout(i,j,kk) = F_vtin(i,j,1) + prd &
                                 * (F_vtin(i,j,1)-F_vtin(i,j,2)) &
                                / (F_wlnph(i,j,1)-F_wlnph(i,j,2))

                  F_gzout(i,j,kk) = F_gzin(i,j,1) - prd * rgasd_8 &
                                 * (F_vtin(i,j,1) + F_vtout(i,j,kk)) * 0.5
!
!******************************************************************************
!                                                                             *
! If:    output pressure   >   hydrostatic pressure on the                    *
!                              last level of the model                        *
!                                                                             *
! Then:  downward extrapolation                                               *
!                                                                             *
! The hypsometric equation is used:                                           *
!                                                                             *
!                         /    \                                              *
!                         | p  |                                              *
!                    _    |  t |                                              *
!  fi  - fi   = - R  T ln |----|                                          (1) *
!    t     b       d      | p  |                                              *
!                         |  b |                                              *
!                         \    /                                              *
!                                                                             *
!  Here the subscript t and b stand respectively for top and bottom of the    *
!  considered layer.                                                          *
!                                                dT                           *
!  We consider a constant temperature lapse rate --- = - L                    *
!                                            _   dfi                          *
!  (e.g. L = STLO) and use the definition of T:                               *
!                                                                             *
!          /                \                                                 *
!          |   fi  - fi     |                                                 *
!  _       |     t     b    |                                                 *
!  T = - L |----------------| ,                                           (2) *
!          |    / T   /   \ |                                                 *
!          | ln |  t / T  | |                                                 *
!          |    \   /   b / |                                                 *
!          \                /                                                 *
!                                                                             *
!  into expression (1) and get an expression for T :                          *
!                                                 b                           *
!               /                     \                                       *
!               |         / p   /   \ |                                       *
!  T  = T   exp | R  L ln |  b / p  | |                                   (3) *
!   b    t      |  d      \   /   t / |                                       *
!               \                     /                                       *
!                                                                             *
!  Then, we use the definition of L, to get an expression for fi :            *
!                                                               b             *
!              /         \                                                    *
!              | T  - T  |                                                    *
!              \  t    b /                                                    *
!  fi  = fi  + ----------- .                                              (4) *
!    b     t        L                                                         *
!                                                                             *
! In the case where L -> 0, we have to revert to expression (1) in which      *
!        _                                                                    *
! we use T = T .                                                              *
!             t                                                               *
!                                                                             *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                                                                             *
! At points where we want to use underground temperatures for extrapolation,  *
! we first determine the layer bottom pressure using (3) rearranged:          *
!                                                                             *
!             /    / T   /   \ \                                              *
!             | ln |  b / T  | |                                              *
!             |    \   /   t / |                                              *
! p  = p  exp | -------------- |                                              *
!  b    t     |      R  L      |                                              *
!             |       d        |                                              *
!             \                /                                              *
!                                                                             *
! In the case where L -> 0, we have to revert to the expression (1) in which  *
!        _                                                                    *
! we use T = T .                                                              *
!             t                                                               *
!                                                                             *
! Then, if the layer bottom pressure is larger than the destination pressure, *
! we proceed with calculation (3) and (4). Otherwise, we update the variables *
! at top and bottom for next layer calculation and iterate.                   *
!                                                                             *
!******************************************************************************
               else if ( pnindex(i) == F_nk ) then

                  do pnund=1,F_nundr
                     if ( F_gzin(i,j,F_nk) > F_gzund(i,j,pnund) ) exit
                  end do
                  prlptop = F_wlnph(i,j,F_nk)
                  prvttop = F_vtin (i,j,F_nk)
                  prfitop = F_gzin (i,j,F_nk)

                  do pn1=pnund,F_nundr

                     prvtbot = F_vtund (i,j,pn1)
                     prfibot = F_gzund (i,j,pn1)

                     if ( abs(prvtbot-prvttop) <= prsmall ) then
                        prlpbot = prlptop + (prfitop-prfibot)/(rgasd_8*prvttop)
                        if ( prlpbot >= prlprso ) then
                           F_vtout(i,j,kk) = prvttop
                           F_gzout(i,j,kk) = prfitop + rgasd_8*prvttop*(prlptop-prlprso)
                           go to 300
!                           cycle ! not good to go to next pn1...
                        end if
                     else
                        prl     = - ( prvttop - prvtbot ) / ( prfitop - prfibot )
                        prlpbot = prlptop + (log(prvtbot/prvttop)) / (rgasd_8*prl)
                        if ( prlpbot >= prlprso ) then
                           F_vtout(i,j,kk) = prvttop * &
                           exp ( rgasd_8 * prl * (prlprso-prlptop))
                           F_gzout(i,j,kk) = prfitop + (prvttop-F_vtout(i,j,kk)) / prl
                           go to 300
!                           cycle ! not good
                        end if
                     end if

                     prlptop = prlpbot
                     prvttop = prvtbot
                     prfitop = prfibot
                  end do

                  prl = stlo_8
                  if ( abs (F_la(i,j)*180./pi_8) >= 49.0 ) prl = .0005
                  F_vtout(i,j,kk) = prvttop * &
                     exp ( rgasd_8 * prl * (prlprso-prlptop))
                  F_gzout(i,j,kk) = prfitop + (prvttop-F_vtout(i,j,kk)) / prl

!******************************************************************************
!  Else, interpolate between appropriate levels                               *
!******************************************************************************
               else
!
!        **********************************************************************
!        *                                                                    *
!        * NOTE ABOUT "F_linbot"                                              *
!        *             --------                                               *
!        *                                                                    *
!        * this parameter is used to force a linear interpolation in a        *
!        * certain number of layers (equal to F_linbot) close to the bottom   *
!        * of the model even if F_cubzt_L is .true.                           *
!        *                                                                    *
!        * it has no effect if F_cubzt_L is .false.                           *
!        *                                                                    *
!        **********************************************************************
!
                  pnkm = pnindex(i)
                  pnk  = pnindex(i) + 1
                  prd  = F_wlnph(i,j,pnk) - F_wlnph(i,j,pnkm)
                  invprd = 1.0/prd
                  pre = prlprso - 0.5 * ( F_wlnph(i,j,pnk) + F_wlnph(i,j,pnkm) )

                  if ( F_cubzt_L .and. ( pnk < F_nk+1-F_linbot ) ) then
                     prr = 0.125 * prd * prd - 0.5 * pre * pre
                     prfm0 = 0.5 * ( F_gzin(i,j,pnk) + F_gzin(i,j,pnkm) )
                     prfm1 = ( F_gzin(i,j,pnk) - F_gzin(i,j,pnkm) ) * invprd
                     prfm2 = - rgasd_8 &
                             * ( F_vtin(i,j,pnk) - F_vtin(i,j,pnkm) ) * invprd
                     prfm3 = - rgasd_8 * ( F_vtin(i,j,pnk) + F_vtin(i,j,pnkm) )
                     prfm3 = ( prfm3 - prfm1 - prfm1 ) * invprd * invprd
                     prfl2 = prfm2 + 2.0 * pre * prfm3
                     F_gzout(i,j,kk) = prfm0 + pre * prfm1 - prr * prfl2
                     F_vtout(i,j,kk) = prfm1 + pre * prfl2 - 2.0 * prr * prfm3
                     F_vtout(i,j,kk) = - F_vtout(i,j,kk) / rgasd_8
                  else
                     prfm0 = 0.5 * ( F_gzin(i,j,pnk) + F_gzin(i,j,pnkm) )
                     prfm1 = ( F_gzin(i,j,pnk) - F_gzin(i,j,pnkm) ) * invprd
                     F_gzout(i,j,kk)= prfm0 + pre * prfm1
                     prfm0 = 0.5 * ( F_vtin(i,j,pnk) + F_vtin(i,j,pnkm) )
                     prfm1 = ( F_vtin(i,j,pnk) - F_vtin(i,j,pnkm) ) * invprd
                     F_vtout(i,j,kk)= prfm0 + pre * prfm1
                  end if
               end if

 300        end do
         end do
      end do
!
!-------------------------------------------------------------------
!
      return
      end
