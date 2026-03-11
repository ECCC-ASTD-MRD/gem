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

!**s/r iau_set_mem - Allocate memory and set pointers
!		               for iau variables

      subroutine iau_set_mem
      use mem_iau
      use ctrl
      use init_options
      use glb_ld
      use lun
      use tr3d
      implicit none

      character(len=24) :: nomvar
      integer j,nv,ivar,dimTot,dimHor,dim3d
!
!     ---------------------------------------------------------------
!
      Ctrl_iau_L= (Iau_period > 0. .and. Iau_interval > 0.)
      if (.not. Ctrl_iau_L) return
      
      if (Lun_out > 0) write (Lun_out,1000)

      ivar= 1 ; nv= 0
      do while (len_trim(Iau_tracers_S(ivar)) > 0)
         nomvar= Iau_tracers_S(ivar)
         do j=1,Tr3d_ntr
            if (trim(Tr3d_name_S(j))=='HU') cycle
            if (trim(Tr3d_name_S(j))==nomvar) then
               nv= nv + 1
               IAU_trname(nv)= nomvar
               IAU_trindx(nv)= j
            endif
         end do
         ivar= ivar + 1
      end do
      IAU_ntr= nv

      dimHor = (l_maxx-l_minx+1) * (l_maxy-l_miny+1)
      dim3d  = dimHor * l_nk
      dimTot = (4+IAU_ntr)*dim3d + dimHor

      allocate (IAU_current(dimTot))

      iau_u (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> IAU_current(        1:)
      iau_v (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> IAU_current(dim3d  +1:)
      iau_t (l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> IAU_current(dim3d*2+1:)
      iau_hu(l_minx:l_maxx,l_miny:l_maxy,1:l_nk  )=> IAU_current(dim3d*3+1:)
      iau_p0(l_minx:l_maxx,l_miny:l_maxy         )=> IAU_current(dim3d*4+1:)
      iau_tr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk,1:IAU_ntr)=> IAU_current(dim3d*4+dimHor+1:)
      
 1000 format( &
      /,'INITIALIZATION OF IAU VARIABLE MEMORY', &
      /,'==============================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
