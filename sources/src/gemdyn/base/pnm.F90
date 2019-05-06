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

!**s/r pnm - calculates MSL pressure
!
      subroutine pnm( F_pnm,   F_vts,   F_fis, F_lnps, F_la, &
                      F_vtund, F_fiund, F_und, &
                      Minx, Maxx, Miny, Maxy, Nk)

      use tdpack
      use glb_ld
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk, F_und
      real, dimension(Minx:Maxx, Miny:Maxy), intent(out) :: F_pnm
      real, dimension(Minx:Maxx, Miny:Maxy), intent(in) :: F_la
      real, dimension(Minx:Maxx, Miny:Maxy, Nk), intent(in) :: F_vts, F_fis, F_lnps
      real, dimension(Minx:Maxx,Miny:Maxy,F_und), intent(in) :: F_vtund, F_fiund
!
!object
!******************************************************************************
!                                                                             *
! The hypsometric equation is used:                                           *
!                                                                             *
!                        /    \                                               *
!                        | p  |                                               *
!                   _    |  t |                                               *
! fi  - fi   = - R  T ln |----|                                          (1)  *
!   t     b       d      | p  |                                               *
!                        |  b |                                               *
!                        \    /                                               *
!                                                                             *
! Here the subscript t and b stand respectively for top and bottom of the     *
! considered layer.                                                           *
!                                               dT                            *
! We consider a constant temperature lapse rate --- = - L                     *
!                                           _   dfi                           *
! (e.g. L = STLO) and use the definition of T:                                *
!                                                                             *
!         /                \                                                  *
!         |   fi  - fi     |                                                  *
! _       |     t     b    |                                                  *
! T = - L |----------------| ,                                           (2)  *
!         |    / T   /   \ |                                                  *
!         | ln |  t / T  | |                                                  *
!         |    \   /   b / |                                                  *
!         \                /                                                  *
!                                                                             *
! into expression (1) and get an expression for p :                           *
!                                                b                            *
!                                                                             *
!             /    / T   /   \ \                                              *
!             | ln |  b / T  | |                                              *
!             |    \   /   t / |                                              *
! p  = p  exp | -------------- |                                              *
!  b    t     |      R  L      |                                              *
!             |       d        |                                              *
!             \                /                                              *
!                                                                             *
! In the case where L -> 0, we have to revert to expression (1) in which      *
!        _                                                                    *
! we use T = T .                                                              *
!             t                                                               *
!                                                                             *
! At points where we want to use underground temperatures for calculation,    *
! we recursively compute the pressure at the bottom of each layer.            *
!                                                                             *
! * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!                                                                             *
! The temperature lapse rate in each virtual layer is either computed using   *
! the provided temperatures or assumed to be the Shumman-Newel lapse rate     *
! for the layer near ground except, when the temperatures exceeds given       *
! critical values:                                                            *
!                                                                             *
! T  = 301.75 - lat / 4                                                       *
!  c                                                                          *
!                                                                             *
! If T  is lower than T , the algorithm ensures that the bottom temperature   *
!     t                c                                                      *
! is not greater than T .                                                     *
!                      c                                                      *
! Else, if T  exceeds T , then the bottom temperature is set to:              *
!           t          c                                                      *
!                                                                             *
!                              2                                              *
! T  = T  - 0.005 * ( T  - T  ) .                                             *
!  b    c              t    c                                                 *
!                                                                             *
!******************************************************************************
!
!arguments
!  Name                        Description
!------------------------------------------------------------
! F_pnm         - MSL pressure
! F_vts         - surface virtual temperature
! F_fis         - surface geopotential height
! F_lnps        - surface log of hydrostatic pressure
! F_la          - geographical latitude (radian)
! F_vtund       - virtual temperatures for underground extrapolation
! F_fiund       - geopotential levels for which virtual temperature is
!                 given for underground extrapolation
! F_und         - number of virtual temperature levels for underground
!                 extrapolation
!               = 0 if no underground temperature is used and the
!                   the traditional scheme will be used
!
!notes
!   All fields in arguments are assumed to be workable on the same grid
!   (fni x fnj). This grid could be the staggered or the non staggered.
!
!   It is important that the data stored in F_vtund and F_fiund be ordered
!   in the proper manner:
!   F_vtund(i,j,1) and F_fiund(i,j,1) --> highest level
!   F_vtund(i,j,2) and F_fiund(i,j,2) --> second highest level
!   ......................................and so on
!
      integer :: i, j, pnund,   pn1
      real :: prl, prvtc, prsmall
      real :: prlptop, prvttop, prfitop
      real :: prlpbot, prvtbot, prfibot

      prsmall = .001

      do j= 1, l_nj
         do i= 1, l_ni
!
!           calculation of critical temperature
!                       --------------------
!
            prvtc = 301.75 - abs( (F_la(i,j) * 180.) / ( 4. * pi_8) )

            do pnund=1,F_und+1
               if ( pnund > F_und ) exit
               if ( F_fis(i,j,nk) > F_fiund(i,j,pnund) ) exit
            end do

            prlptop = F_lnps(i,j,nk)
            prvttop = F_vts(i,j,nk)
            prfitop = F_fis(i,j,nk)

            do pn1=pnund,F_und

               if ( prvttop <= prvtc ) then
                    prvtbot = min( F_vtund(i,j,pn1),  prvtc )
               else
                    prvtbot = prvtc - 0.005 * ( prvttop - prvtc ) **2
               end if

               prfibot  = F_fiund (i,j,pn1)

               if ( abs(prvtbot-prvttop) <= prsmall ) then
                  prlpbot = prlptop + (prfitop-prfibot)/(rgasd_8*prvttop)
               else
                  prl     = - ( prvttop - prvtbot ) / ( prfitop - prfibot )
                  prlpbot = prlptop + (log(prvtbot/prvttop)) / (rgasd_8*prl)
               end if

               prlptop = prlpbot
               prvttop = F_vtund(i,j,pn1)
               prfitop = prfibot
            end do

            if ( prvttop <= prvtc ) then
                 prvtbot = min( 1.0d0*prvttop + stlo_8 * 1.0d0*prfitop,  1.0d0*prvtc)
            else
                 prvtbot = prvtc - 0.005 * ( prvttop - prvtc ) **2
            end if
!
!           calculation of MSL pressure
!                       ------------
!
            if ((abs(prvtbot-prvttop) <= prsmall) .or. (prfitop <= 0.0)) then
                 F_pnm(i,j) = exp (prlptop+prfitop/(rgasd_8*prvttop))
            else
                 prl = - ( prvttop - prvtbot ) / ( prfitop )
                 F_pnm(i,j)=exp (prlptop+(log(prvtbot/prvttop))/(rgasd_8*prl))
            end if

         end do
      end do

      return
      end
