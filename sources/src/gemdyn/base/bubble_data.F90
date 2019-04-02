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

!**s/r bubble_data - generates initial condition for Robert's bubble
!                    experiment (Robert 1993 JAS)
!
      subroutine bubble_data ( F_u, F_v, F_t, F_s, F_q, F_topo, F_sls,&
                               Mminx,Mmaxx,Mminy,Mmaxy,nk )
      use dynkernel_options
      implicit none
#include <arch_specific.hf>

      integer Mminx,Mmaxx,Mminy,Mmaxy,nk
      real F_u      (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_v      (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_t      (Mminx:Mmaxx,Mminy:Mmaxy,nk), &
           F_s      (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_topo   (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_sls    (Mminx:Mmaxx,Mminy:Mmaxy   ), &
           F_q      (Mminx:Mmaxx,Mminy:Mmaxy,nk+1)

!
!     ---------------------------------------------------------------
!

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

         call bubble_fislP_data(F_u, F_v, F_t, F_s, F_q, F_topo, F_sls,&
                               Mminx,Mmaxx,Mminy,Mmaxy,nk )

      else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

         call bubble_fislH_data(F_u, F_v, F_t, F_s, F_q, F_topo, F_sls,&
                               Mminx,Mmaxx,Mminy,Mmaxy,nk )

      else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H') then

         call bubble_expoH_data(F_u, F_v, F_t, F_s, F_q, F_topo, F_sls,&
                               Mminx,Mmaxx,Mminy,Mmaxy,nk )

      end if
!
!     -----------------------------------------------------------------
!
      return
      end subroutine bubble_data

