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

      subroutine glbstat ( F_field, F_var_S, F_from_S   , &
                           Minx,Maxx,Miny,Maxy,Mink,Maxk, &
                           F_i0,F_in,F_j0,F_jn,F_k0,F_kn )
      use stat_mpi, only: statf_dm
      use gem_options
      use step_options
      use glb_ld
      use ptopo
      implicit none
#include <arch_specific.hf>

      character(len=*) F_var_S,F_from_S
      integer Minx,Maxx,Miny,Maxy,Mink,Maxk
      integer F_i0,F_in,F_j0,F_jn,F_k0,F_kn
      real F_field (Minx:Maxx,Miny:Maxy,Mink:Maxk)

!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_field       I         Field to be operated on
! F_var_S       I         User provided string to define F_field
! F_from_S      I         User provided string to comment the statfld
! F_i0,F_j0     I         Global lower-left indexes of the sub-domain
!                            on which to perform statistics
! F_in,F_jn     I         Global upper-right indexes of the sub-domain
!                            on which to perform statistics
! F_k0,F_kn     I         Range of levels on which to perform statistics
!----------------------------------------------------------------


      integer nk,rx
      real, dimension(:,:), allocatable :: wk1
!
!----------------------------------------------------------------------
!
      read (Lctl_rxstat_S(5:),*) rx

      if (Lctl_rxstat_S(1:3)=='GLB') then

         nk = Maxk-Mink+1
         if (Ptopo_myproc == 0) then
            allocate (wk1(G_ni*G_nj,Mink:Maxk))
         else
            allocate (wk1(1,1))
         end if
         call glbcolc (wk1,G_ni,G_nj,F_field,Minx,Maxx,Miny,Maxy,nk)

         if (Ptopo_myproc == 0)  then
            call statfld (wk1 ,F_var_S, Lctl_step, F_from_S, &
                           1,G_ni, 1,G_nj, Mink,Maxk, &
                           F_i0,F_j0,F_k0, F_in,F_jn,F_kn,rx)
         end if
         deallocate (wk1)

      else

         call statf_dm ( F_field,F_var_S,Lctl_step,F_from_S,&
                         Minx,Maxx,Miny,Maxy,Mink,Maxk     ,&
                         F_i0,F_j0,F_k0,F_in,F_jn,F_kn,rx)
      end if
!
!----------------------------------------------------------------------
!
      return
      end
