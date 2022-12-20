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

module phy_terminate_mod

  private
  public :: phy_terminate

contains

  !/@*
  function phy_terminate() result(F_istat)
      use cpl_itf, only: cpl_terminate
      use series_mod, only: series_terminate
      use phy_status, only: PHY_NONE, PHY_CTRL_INI_OK, phy_error_L, phy_init_ctrl
      implicit none

    integer :: F_istat  !Return status (RMN_OK or RMN_ERR)

    !@authors Desgagne, Chamberland, McTaggart-Cowan, Spacek -- Spring 2014
    !@revision

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

    integer :: istat
    ! ---------------------------------------------------------------------
    F_istat = RMN_ERR
    if (phy_init_ctrl == PHY_NONE) then
       F_istat = PHY_NONE
       return
    else if (phy_init_ctrl /= PHY_CTRL_INI_OK) then
       call msg(MSG_ERROR,'(phy_terminate) Physics not properly initialized.')
       return
    endif

    ! Set error return value
    F_istat = RMN_ERR

    ! Shut down time series
    istat = series_terminate()

    ! Shut down coupler
    call cpl_terminate (.true.)

    if (phy_error_L) return

    ! Successful completion
    F_istat = RMN_OK

    ! ---------------------------------------------------------------------
    return
  end function phy_terminate

end module phy_terminate_mod
