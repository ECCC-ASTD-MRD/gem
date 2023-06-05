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

!**s/r HOR_bndry

      subroutine HOR_bndry_hlt ()
      use ctrl
      use glb_ld
      use HORgrid_options
      use mem_nest
      use mem_tracers
      use gmm_pw
      use tr3d
      implicit none
#include <arch_specific.hf>

      integer i,j,k,n
      integer jin, jj
!
!     ---------------------------------------------------------------
!
      if (Grd_yinyang_L) then
!$omp single
         do n=1, Tr3d_ntr
            call yyg_xchng (tracers_P(n)%pntr, l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                            G_nk, .true., 'CUBIC', .false.)
         end do
!$omp end single
      else
         call nest_HOR_gwa_hlt ()
         call pressure_hlt ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,pw_log_pt, &
                             pw_pm_plus_8,pw_p0_plus_8, &
                             l_minx,l_maxx,l_miny,l_maxy,l_nk,1 )
      endif

      if ( Ctrl_theoc_L ) then
      do n=1,Tr3d_ntr
         if (l_north) then
!$omp do
           do k=1,G_nk
               jin = l_nj-pil_n
               do j=1,pil_n
                  jj  = l_nj-pil_n+j
                  do i=1,l_ni
                     tracers_P(n)%pntr(i,jj,k)=tracers_P(n)%pntr(i,jin,k)
                  end do
               end do
            end do
!$omp end do
         end if
         if (l_south) then
!$omp do
            do k=1,G_nk
               jin = pil_s+1
               do j=1,pil_s
                  jj  = pil_s-j+1
                  do i=1,l_ni
                     tracers_P(n)%pntr(i,jj,k) = tracers_P(n)%pntr(i,jin,k)
                  end do
               end do
            end do
!$omp end do
         end if
      end do
      endif
!
!     ---------------------------------------------------------------
!
      return
      end subroutine HOR_bndry_hlt
