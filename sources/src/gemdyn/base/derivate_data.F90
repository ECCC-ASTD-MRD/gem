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

!** s/r derivate_data

      subroutine derivate_data ( F_zd, F_w, F_u, F_v, F_t, F_s, F_q,&
                                 Minx, Maxx, Miny, Maxy, Nk        ,&
                                 F_zd_L, F_w_L )
      use glb_ld
      use cstv
      use dynkernel_options
      implicit none

      integer, intent(in) ::  Minx, Maxx, Miny, Maxy, Nk
      logical, intent(in) ::  F_zd_L, F_w_L
      real, dimension(Minx:Maxx,Miny:Maxy   ),  intent(inout ):: F_s
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_zd, F_w,F_q
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_u, F_v, F_t
!      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(inout) :: F_q
!
!     ________________________________________________________________
!
      if( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' .or. &
          trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H' ) then
         call fislh_pres ( F_q, F_s, F_t, &
                           Minx, Maxx, Miny, Maxy, Nk)
         call fislh_diag_zd_w( F_zd, F_w, F_u, F_v, F_t, F_s, F_q,&
                               Minx, Maxx, Miny, Maxy, Nk,&
                               F_zd_L, F_w_L)
      end if

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then
         F_q= 0.
         F_s(1:l_ni,1:l_nj) = log(F_s(1:l_ni,1:l_nj)/Cstv_pref_8)
         call diag_zd_w ( F_zd, F_w, F_u, F_v, F_t, F_s,&
                          Minx, Maxx, Miny, Maxy, Nk   ,&
                          F_zd_L, F_w_L )
      end if
!
!     ________________________________________________________________
!
      return
      end subroutine derivate_data
