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
module dynkernel_options
   implicit none
   public
   save

   !# Main selector for dynamical kernel
   !# * "DYNAMICS_FISL_P" : Fully implicit SL in pressure
   !# * "DYNAMICS_FISL_H" : Fully implicit SL in height
   character(len=32) :: Dynamics_Kernel_S = 'DYNAMICS_FISL_P'
   namelist /dyn_kernel/ Dynamics_Kernel_S

   !# * True-> hydrostatic
   !# * False-> non-hydrostatic
   logical :: Dynamics_hydro_L = .false.
   namelist /dyn_kernel/ Dynamics_hydro_L

   !# True-> auto barotropic option
   logical :: Schm_autobar_L = .false.
   namelist /dyn_kernel/ Schm_autobar_L

   !# True-> Implicit metric terms in the LHS (needs 3D iterative solver)
   !# False-> Simplified approach (allows both direct and iterative solvers)
   logical :: FISLH_LHS_metric_L = .false.

   logical :: Dynamics_hauteur_L, Dynamics_FISL_L, Dynamics_autobar_L

contains

!**s/r dynkernel_nml - Read namelist dyn_kernel

      integer function dynkernel_nml (F_unf)
      use lun
      use clib_itf_mod
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
      integer err
!
!-------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         dynkernel_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=dyn_kernel)
         return
      end if

      dynkernel_nml= -1 ; nml_must= .false. ; nml_S= 'dyn_kernel'

      rewind(F_unf)
      read (F_unf, nml=dyn_kernel, end= 1001, err=1003)
      dynkernel_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         dynkernel_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (dynkernel_nml < 0 ) return
      if ((Lun_out>=0).and.(dynkernel_nml==0)) write (Lun_out, 6004) trim(nml_S)
      dynkernel_nml= 1

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      err = clib_toupper ( Dynamics_Kernel_S )
!
!-------------------------------------------------------------------
!
      return
      end function dynkernel_nml

   !/@*
   function dynkernel_options_init() result(F_istat)
      implicit none
      !@object Additional initialisation steps before reading the nml
      integer :: F_istat
      !*@/
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      !----------------------------------------------------------------
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      !----------------------------------------------------------------
      return
   end function dynkernel_options_init

end module dynkernel_options
