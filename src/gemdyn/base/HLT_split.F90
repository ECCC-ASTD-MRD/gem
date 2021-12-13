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

      subroutine HLT_split ( F_i0, F_in, F_local_ni, F_start, F_end )
      use omp_lib
      implicit none

      integer, intent(IN)  :: F_i0, F_in
      integer, intent(OUT) :: F_local_ni, F_start, F_end

      integer nthreads, tid!, omp_get_num_threads, &
                            ! omp_get_thread_num
      integer mpx,irest,glb_np
!
!     ---------------------------------------------------------------
!
      tid = omp_get_thread_num()
      nthreads = omp_get_num_threads()
      glb_np   = F_in-F_i0+1

      mpx        = mod( tid, nthreads)
      F_local_ni = glb_np / nthreads
      irest  = glb_np  - f_local_ni * nthreads
      F_start = mpx * F_local_ni + F_i0

      if ( mpx < irest ) then
         F_local_ni   = F_local_ni + 1
         F_start = F_start + mpx
      else
         F_start = F_start + irest
      endif
      F_end= F_start + F_local_ni - 1
!      print*, tid,glb_np,NTHREADS,F_start,F_end,F_local_ni

!
!     ---------------------------------------------------------------
!
      return
      end subroutine HLT_split
