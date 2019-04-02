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

!**s/p calomeg_w - compute vertical velocity in pressure coordinates
!                    from advection
!
      subroutine calomeg_w (F_ww,F_st1,F_sl,F_wt1,F_tt1,F_lnp,Minx,Maxx,Miny,Maxy,Nk)
      use dynkernel_options
      use tdpack
      use glb_ld
      use cstv
      use ver
      use type_mod
      implicit none
#include <arch_specific.hf>
!
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out) :: F_ww
      real, dimension(Minx:Maxx,Miny:Maxy),    intent(in)  :: F_st1
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in)  :: F_wt1, F_tt1, F_lnp
      real, dimension(Minx:Maxx,Miny:Maxy),    intent(in)  :: F_sl
!
!author
!     Claude Girard et Andre Plante avril 2008.
!
!object
!	compute vertical velocity in hydrostatic pressure coordinates
!
!        omega = dpi/dt ~ -g*rau*dz/dt = -g*rau*w
!
!        where : pi  is the hydrostatic pressure
!                rau the density
!                w   real vertical wind
!
!arguments
!  Name               Description
!---------------------------------------------------
! F_ww               dpi/dt (Pa/s)
! F_st1               s at time t1 = exp(PIt1/Zsruf)
! F_wt1               real vertical wind
! F_tt1               virtual temperature
!


      integer :: i,j,k
      real, dimension(l_ni,l_nj) :: t1,t2
!     __________________________________________________________________
!

!$omp parallel private(t1,t2,i,j,k)
!$omp do
      do k=1,l_nk
         if(trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then
            do j=1,l_nj
               do i=1,l_ni
                  t1(i,j)= F_lnp(i,j,k)
               end do
            end do
         else
            do j=1,l_nj
               do i=1,l_ni
                  t1(i,j)= Ver_a_8%t(k) + Ver_b_8%t(k)*F_st1(i,j) + Ver_c_8%t(k)*F_sl(i,j)
               end do
            end do
         end if
         call vsexp(t2,t1,l_ni*l_nj)
         do j=1,l_nj
            do i=1,l_ni
               F_ww (i,j,k) = -grav_8 *F_wt1(i,j,k)*t2(i,j)/ &
                              (rgasd_8*F_tt1(i,j,k))
            end do
         end do
      end do
!$omp enddo
!$omp end parallel
!     __________________________________________________________________
!
      return
      end
