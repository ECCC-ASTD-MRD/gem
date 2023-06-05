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

      logical function decomp ( F_npts, F_min, F_max, F_lni, F_npartiel, F_halo, F_start,  &
                                F_alongx_L, F_fill_L, F_npe, F_lowestsize, F_checkparti_L, &
                                F_relax)
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(in) :: F_alongx_L, F_fill_L, F_checkparti_L
      integer, intent(inout) :: F_min, F_max, F_Npartiel, F_start
      integer, intent(in) :: F_npts, F_halo, F_npe
      integer, intent(inout):: F_lni(F_npe)
      integer, intent(in) :: F_lowestsize, F_relax

      logical  :: check_parti
      integer  :: RPN_COMM_limit_2,rpn_comm_topo_2
      external :: RPN_COMM_limit_2,rpn_comm_topo_2,check_parti

      integer :: istat
      integer, dimension(F_npe) :: depl
!
!-------------------------------------------------------------------
!
      decomp = .false. ; F_lni= 0

      if (F_checkparti_L ) then

         istat= RPN_COMM_limit_2 (0, F_npe, 1, F_npts, F_min, F_max, &
                                  F_lni, depl, F_relax)
         if (istat >= 0) decomp = (  minval(F_lni) >= F_lowestsize )
                                  
      else

         istat= rpn_comm_topo_2( F_npts, F_min, F_max, F_lni(1), F_npartiel, &
                  F_halo, F_start, F_alongx_L, F_fill_L, F_relax, .false. )

         decomp = .not.(istat < 0)
         if (F_lowestsize > 0) decomp = decomp .and. (F_lni(1) >= F_lowestsize)
            
      end if
      
      if (.not.decomp) write (output_unit,1001) F_npts, F_npe
 1001 format(/' DECOMP: illegal partitionning ====> ',i7,' / ',i7)
!
!-------------------------------------------------------------------
!
      return
      end

