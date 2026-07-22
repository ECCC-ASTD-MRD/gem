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

module cpl_update_mod

   private
   public :: cpl_update

contains

      subroutine cpl_update (F_f2u, F_name_S, F_pos, F_ni, rho, u, v, vmod, cmu, cplu)
      use cpl_mod
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_name_S
      integer, intent(in)          :: F_ni
      integer, dimension(2,F_ni),           intent(in)    :: F_pos
      real   , dimension  (F_ni),           intent(inout) :: F_f2u
      real   , dimension  (F_ni), optional, target,&
                                            intent(in)    :: rho, u, v, vmod
      real   , dimension  (F_ni), optional, target,&
                                            intent(inout) :: cmu
      logical, optional, intent(out) :: cplu

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version

!Purpose
! Update field F_f2u

      real   , dimension(F_ni), target  :: rho_w, u_w, v_w, vmod_w, cmu_w
      real   , dimension(:)   , pointer :: F_rho, F_u, F_v, F_vmod, F_cmu
      logical cplu_w

!     ---------------------------------------------------------------

      F_rho  => rho_w
      F_u    => u_w
      F_v    => v_w
      F_vmod => vmod_w
      F_cmu  => cmu_w
      cplu_w = .false.

      if (present(rho))  F_rho => rho
      if (present(u))    F_u   => u
      if (present(v))    F_v   => v
      if (present(vmod)) F_vmod => vmod
      if (present(cmu))  F_cmu => cmu

      if ( cpl_ocn_L ) &
        call cplocn_update(F_f2u, F_name_S, F_pos, F_ni, cplu_w, &
                           F_rho, F_u, F_v, F_vmod, F_cmu)

      if (present(cplu)) cplu=cplu_w
!
!     ---------------------------------------------------------------
!
      return
      end subroutine cpl_update

end module cpl_update_mod

