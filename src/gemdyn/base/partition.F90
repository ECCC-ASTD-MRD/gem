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
!**s/r partition
!
      integer function partition ( F_npts, F_alongxy, F_halo, &
                                   F_min, F_max, F_ni, F_i0 )
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
      character(len=*) :: F_alongxy
      integer F_npts, F_halo, F_min, F_max, F_ni, F_i0

      character(len=24) :: Part_type_S
      integer npe,val1,val2,irest,mpx,lstp
!
!-------------------------------------------------------------------
!
      Part_type_S = 'ALLEQUALLBUT1'
!      Part_type_S = 'MOSTUNIFORM'

      if (F_alongxy == 'NPEX') then
         npe = Ptopo_npex
         mpx = Ptopo_mycol
      else
         npe = Ptopo_npey
         mpx = Ptopo_myrow
      end if

      partition = -1
      F_i0 = F_npts + 1

      val1 = 0

      if (Part_type_S == 'ALLEQUALLBUT1') then
         val1 = (F_npts + npe - 1)/npe
         val2 =  F_npts -(npe - 1)*val1
         if (val2 > 0) F_i0 = mpx*val1 + 1
      end if

      if (Part_type_S == 'MOSTUNIFORM') then
         val1  = F_npts / npe
         irest = F_npts - val1 * npe
         F_i0  = mpx * val1 + 1
         if ( mpx < irest ) then
            val1 = val1 + 1
            F_i0 = F_i0 + mpx
         else
            F_i0 = F_i0 + irest
         end if
      end if

      if (F_i0 > F_npts)  then
         write (output_unit,9001) F_npts,npe,trim(Part_type_S)
         return
      end if

      partition = 0

      lstp  =  min(F_i0 + val1 - 1, F_npts)
      F_ni  =  lstp - F_i0 + 1

      F_min =  1    - F_halo
      F_max =  F_ni + F_halo

      F_max =  F_max + 1 - mod(val1,2)

 9001 format (/'  Invalid data decomposition (S/R PARTITION)'         /,&
               '  Cannot partition',i6,' points into',i6,' subdomains'/,&
               '  Under ',a,' method'/)
!
!-------------------------------------------------------------------
!
      return
      end

