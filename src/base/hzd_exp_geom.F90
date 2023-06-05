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

!**s/r hzd_exp_geom

      subroutine hzd_exp_geom()
      use gem_options
      use geomh
      use glb_ld
      use HORgrid_options
      use hzd_mod
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer :: j
      real(kind=REAL64) :: aaa,bbb,ccc,ddd,dx_8
!
!     ---------------------------------------------------------------
!
      dx_8 = Grd_dx * pi_8 / 180.0d0

      allocate ( Hzd_geom_q (l_miny:l_maxy,5),&
                 Hzd_geom_u (l_miny:l_maxy,5),&
                 Hzd_geom_v (l_miny:l_maxy,5) )

      do j= 2-G_haloy, l_nj+G_haloy-1

         aaa = (sin(geomh_yv_8(j))-sin(geomh_yv_8(j-1))) * geomh_invcy2_8(j)
         bbb = (geomh_sy_8(j  ) - geomh_sy_8(j-1)) / geomh_cyv2_8(j-1)
         ccc = (geomh_sy_8(j+1) - geomh_sy_8(j  )) / geomh_cyv2_8(j  )
         ddd = 1.0d0 / (dx_8 * (sin(geomh_yv_8 (j)) - sin(geomh_yv_8 (j-1))))

         Hzd_geom_q(j,2)= ddd*aaa/dx_8
         Hzd_geom_q(j,3)= Hzd_geom_q(j,2)
         Hzd_geom_q(j,4)= ddd*dx_8/bbb
         Hzd_geom_q(j,5)= ddd*dx_8/ccc
         Hzd_geom_q(j,1)=-(Hzd_geom_q(j,2)+Hzd_geom_q(j,3)+ &
                             Hzd_geom_q(j,4)+Hzd_geom_q(j,5))

         aaa= (geomh_sy_8(j+1) - geomh_sy_8(j)) * geomh_invcy2_8(j)
         bbb= geomh_cy2_8(j  ) / (sin(geomh_yv_8 (j)) - sin(geomh_yv_8 (j-1)))
         ccc= geomh_cy2_8(j+1) / (sin(geomh_yv_8 (j+1)) - sin(geomh_yv_8 (j)))
         ddd= 1.d0 / ( dx_8 * (geomh_sy_8(j+1) - geomh_sy_8(j)) )

         Hzd_geom_v(j,2)=ddd*aaa/dx_8
         Hzd_geom_v(j,3)=Hzd_geom_v(j,2)
         Hzd_geom_v(j,4)=ddd*dx_8*bbb
         Hzd_geom_v(j,5)=ddd*dx_8*ccc
         Hzd_geom_v(j,1)=-(Hzd_geom_v(j,2)+Hzd_geom_v(j,3)+ &
                             Hzd_geom_v(j,4)+Hzd_geom_v(j,5))

      end do

      Hzd_geom_u = Hzd_geom_q
!
!     ---------------------------------------------------------------
!
      return
      end
