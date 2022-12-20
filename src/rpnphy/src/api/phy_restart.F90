!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

module phy_restart_mod
   use clib_itf_mod, only: clib_toupper
  use cpl_itf, only: cpl_restart
  implicit none
  private
  public :: phy_restart

  character(len=4), parameter, public :: PHY_RESTART_READ  = 'R'
  character(len=4), parameter, public :: PHY_RESTART_WRITE = 'W'

contains

  !/@*
  function phy_restart(F_WorR_S, F_spin_L) result(F_istat)
    implicit none
 !#TODO: F_spin_L never used
    character(len=*),intent(in) :: F_WorR_S
    logical,         intent(in) :: F_spin_L
    integer :: F_istat  !Return status (RMN_OK or RMN_ERR)

    !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2015
    !@revision
    !*@/

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

    integer :: istat
    character(len=32) :: WorR_S

    ! ------------------------------------------------------------------
    WorR_S = adjustl(F_WorR_S)
    istat = clib_toupper(WorR_S)

! physics has currently nothing to do for restart

    if (WorR_S == PHY_RESTART_READ) then
    endif

    if (WorR_S == PHY_RESTART_WRITE) then
    endif

    if (F_spin_L) then
       print *, 'F_spin_L ignored'
    endif

! coupling may have something to do for restart

    call cpl_restart(WorR_S)

    F_istat = RMN_OK

    ! ------------------------------------------------------------------
    return
  end function phy_restart

end module phy_restart_mod
