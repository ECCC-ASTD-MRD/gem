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

!**s/r adj_ss2topo- Obtain new surface pressure by comparing orographies and
!                   using given calculated 3d pressure field F_pres

      subroutine adj_ss2topo (F_ssq0, F_newtopo, F_pres, F_oldtopo, F_vt,&
                              Minx, Maxx, Miny, Maxy, nk, F_i0, F_in, F_j0, F_jn)
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
!object
!       computes hydrostatic surface pressure over model topography,
!
!  Assuming hydrostatic equilibrium and linear temperature lapse
!  rate in a layer, one can obtain analiticaly the following equation
!  by vertical integration:
!
!             /    / T   /   \ \                   1/
!             | ln |  b / T  | |                   / R  L
!             |    \   /   t / |        / T   /   \   d
! p  = p  exp | -------------- |  =  p  |  b / T  |          (1)
!  b    t     |      R  L      |      t \   /   t /
!             |       d        |
!             \                /
!
!
!  where the subscript t and b stand respectively for top and bottom of
!  the considered layer and L is the temperature lapse rate in the
!  layer defined as follow:
!
!              T  - T
!               t    b
!        L  =  -------                                       (2)
!              gz - gz
!                t    b
!
!  The use of equation (1) and (2) is not convenient when the lapse
!  rate is very small (nearly isothermal conditions) since the exponent
!  in (1) becomes infinite.
!
!  In this case, the hypsometric equation is used:
!
!              / gz  -  gz  \
!              |   t      b  |
!  p  = p  exp | ----------- |                                (3)
!   b    t     |   R  T      |
!              \    d       /
!
!       where T is the mean temperature in the layer.
!  Since equation (3) is used when T - T  --->  0   ,
!                                   t   b
!  the mean temperature is then taken from T  .
!                                           t
!
!  The algorithm is first looking for the closest analysis layer that
!  is found just above the destination terrain.  From that point, this
!  level is considered as the top of the layer.
!
!  At this point , the idea is to compute p  using equation (1) and (2)
!  except when L --> 0.                    b
!
!  When L --> 0,  equation (3) is used where T= T   .
!                                                t
!  Now if the found closest source layer is the lowest analysis level,
!  (this is where the destination model terrain is under the analysis
!  terrain) then there is no known bottom layers in (1) to (3).
!  In this case T and  L  are obtained assuming Schuman-Newel lapse
!  rate under analysis ground.
!
!arguments
!  Name                        Description
!------------------------------------------------------------
! F_newtopo           destination surface geopotential
! F_oldtopo           source surface geopotential
! F_pres              destination (3d pressure field)
! F_ssqn              source      (surface pressure)
! F_vt                source virtual temperature
! NK                  number of levels to look through
!------------------------------------------------------------------
      integer, intent(in) :: Minx,Maxx,Miny,Maxy,nk,F_i0,F_in,F_j0,F_jn
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in) :: F_oldtopo, F_newtopo
      real, dimension(Minx:Maxx,Miny:Maxy,nk), intent(in) :: F_vt, F_pres
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out) :: F_ssq0

      integer :: i,j,k,ik
      real(kind=REAL64) :: difgz, lapse, ttop, tbot, cons, q1, q2, q3, x0, xm, xp, &
                           aa, bb, cc, dd, zak, zbk, zck, invdet
      real, dimension(Minx:Maxx,Miny:Maxy,Nk) :: vma, vmb, vmc, gz_temp
!
!     ---------------------------------------------------------------
!

      do j=F_j0,F_jn
         do i=F_i0,F_in
            gz_temp(i,j,nk) = F_oldtopo(i,j)
         end do
      end do

      do k = 1,Nk
         do j = F_j0, F_jn
            do i = F_i0, F_in
               x0=F_pres(i,j,k)
               if (k == 1) then
                  xm = F_pres(i,j,1)
                  xp = F_pres(i,j,2)
                  aa = F_pres(i,j,3)-x0
                  bb = F_pres(i,j,2)-x0
               else if (k == nk) then
                  xm = F_pres(i,j,Nk-1)
                  xp = F_pres(i,j,Nk)
                  aa = F_pres(i,j,Nk-1)-x0
                  bb = F_pres(i,j,Nk-2)-x0
               else
                  xm = F_pres(i,j,k-1)
                  xp = F_pres(i,j,k+1)
                  aa = xm-x0
                  bb = xp-x0
               end if

               q1 = log(xp/xm)
               q2 = xp-xm
               q3 = (xp*xp - xm*xm)*0.5

               q3 = q3-x0*(2.0*q2-x0*q1)
               q2 = q2-x0*q1
               cc = aa*aa
               dd = bb*bb
               invdet = aa*dd-bb*cc
               invdet = 0.5/invdet
               vma(i,j,k) = (dd*q2-bb*q3)*invdet
               vmc(i,j,k) = (aa*q3-cc*q2)*invdet
               vmb(i,j,k) = q1*0.5-vma(i,j,k)-vmc(i,j,k)
            end do
         end do
      end do

      do j= F_j0, F_jn

         do i= F_i0, F_in
            zak = 2.0 * rgasd_8 * vma(i,j,nk)
            zbk = 2.0 * rgasd_8 * vmb(i,j,nk)
            zck = 2.0 * rgasd_8 * vmc(i,j,nk)
            gz_temp(i,j,nk-1) = zak * F_vt(i,j,nk-1) + zbk * F_vt(i,j,nk) + &
                                zck * F_vt(i,j,nk-2) + gz_temp(i,j,nk)
         end do

         do k = 1, nk-2
            ik  = nk-1-k
            do i= F_i0, F_in
               zak = 2.0 * rgasd_8 * vma(i,j,ik+1)
               zbk = 2.0 * rgasd_8 * vmb(i,j,ik+1)
               zck = 2.0 * rgasd_8 * vmc(i,j,ik+1)
               gz_temp(i,j,ik) = zak* F_vt(i,j,ik  ) + zbk* F_vt(i,j,ik+1) + &
                                 zck* F_vt(i,j,ik+2) + gz_temp(i,j,ik+2)
            end do
         end do

         do i= F_i0, F_in

            difgz = gz_temp(i,j,nk) - F_newtopo(i,j)

            if ( (abs(difgz) < 1.e-5) .or. (difgz > 0.) ) then

!          surface of target grid is below the surface of source grid
!          we assume SCHUMAN-NEWELL Lapse rate under ground to obtain
!          an estimates of the temperature at the target grid surface

               lapse = stlo_8
               k     = nk

            else

!          surface of target grid is above the surface of source grid
!          Then we are looking for the level in the source grid that
!          is just above the surface of the target grid...

               do k = nk, 2, -1
                  difgz = gz_temp(i,j,k) - F_newtopo(i,j)
                  if ( difgz > 0. ) exit
               end do

               lapse = - ( F_vt(i,j,k)-F_vt(i,j,k+1) ) / &
                         ( gz_temp(i,j,k)-gz_temp(i,j,k+1) )

            end if

            ttop = F_vt(i,j,k)
            tbot = ttop + lapse * difgz

            if ( abs(lapse) < 1E-10 ) then
               F_ssq0(i,j) = F_pres(i,j,k) * exp ( difgz/(rgasd_8*ttop) )
            else
               cons = 1. / ( rgasd_8 * lapse )
               F_ssq0(i,j) = F_pres(i,j,k) * ( tbot/ttop ) ** cons
            end if
         end do
      end do

!
!     ---------------------------------------------------------------
!
      return
      end
