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

      subroutine adz_od_count ( F_export )
      use adz_mem
      use step_options
      implicit none
#include <arch_specific.hf>

      type(ADZ_SLOD), intent(IN) :: F_export

      integer lcl_sum,dim,err
      real rt
!
!     ---------------------------------------------------------------
!
      dim = l_ni*l_nj*l_nk
      lcl_sum= sum(F_export%stk(2,:))
      rt= real(lcl_sum)/real(dim)*100.
      call rpn_comm_REDUCE ( rt, nexports(Lctl_step), 1,&
                            "MPI_REAL","MPI_SUM",0,"grid",err )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_count
