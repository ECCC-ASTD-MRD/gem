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

      subroutine set_num_threads ( F_nthreads, kount )
      use lun
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer F_nthreads, kount

!
!-------------------------------------------------------------------
!
      if (F_nthreads /= Ptopo_npeOpenMP) then

         call omp_set_num_threads (F_nthreads)

         Ptopo_npeOpenMP = F_nthreads

         if ((Lun_out > 0).and.(kount == 0)) then
            write (Lun_out,9000) Ptopo_npeOpenMP
         end if

      end if

 9000 format (/' Ptopo_npeOpenMP reset to: ',I7)
!
!-------------------------------------------------------------------
!
      return
      end
