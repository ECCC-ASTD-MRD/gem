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
module dcmip_options
      use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   !# Dcmip case selector
   !# * Dcmip_case=0 (NONE)
   !# * Dcmip_case=11 (3D deformational flow)
   !# * Dcmip_case=12 (3D Hadley-like meridional circulation)
   !# * Dcmip_case=13 (2D solid-body rotation of thin cloud-like tracer in the presence of orography)
   !# * Dcmip_case=20 (Steady-state at rest in presence of oro.)
   !# * Dcmip_case=21 (Mountain waves over a Schaer-type mountain)
   !# * Dcmip_case=22 (As 21 but with wind shear)
   !# * Dcmip_case=31 (Gravity wave along the equator)
   !# * Dcmip_case=41X (Dry Baroclinic Instability Small Planet)
   !# * Dcmip_case=43 (Moist Baroclinic Instability Simple physics)
   !# * Dcmip_case=161 (Baroclinic wave with Toy Terminal Chemistry)
   !# * Dcmip_case=162 (Tropical cyclone)
   !# * Dcmip_case=163 (Supercell Small Planet)
   integer :: Dcmip_case = 0
   namelist /dcmip/ Dcmip_case

   !# Type of precipitation/microphysics
   !# * Dcmip_prec_type=-1 (none)
   !# * Dcmip_prec_type=0 (Kessler Microphysics)
   !# * Dcmip_prec_type=1 (Reed-Jablonowski Large-scale precipitation)
   integer :: Dcmip_prec_type = -1
   namelist /dcmip/ Dcmip_prec_type

   !# Type of planetary boundary layer
   !# * Dcmip_pbl_type=-1 (none)
   !# * Dcmip_pbl_type=0 (Reed-Jablonowski Boundary layer)
   !# * Dcmip_pbl_type=1 (Georges Bryan Boundary Layer)
   integer :: Dcmip_pbl_type = -1
   namelist /dcmip/ Dcmip_pbl_type

   !# Reset in Dcmip_set as: Dcmip_physics_L = Dcmip_pbl_type/=-1 or Dcmip_prec_type/=-1
   logical :: Dcmip_physics_L = .false.

   !# Account for moisture
   !# * Dcmip_moist=0 (dry)
   !# * Dcmip_moist=1 (moist)
   integer :: Dcmip_moist = 1
   namelist /dcmip/ Dcmip_moist

   !# Set lower value of Tracer in Terminator
   !# * Dcmip_lower_value=0 (free)
   !# * Dcmip_lower_value=1 (0)
   !# * Dcmip_lower_value=2 (1.0e-15)
   integer :: Dcmip_lower_value = 0
   namelist /dcmip/ Dcmip_lower_value

   !# Do Terminator chemistry if T
   logical :: Dcmip_Terminator_L = .false.
   namelist /dcmip/ Dcmip_Terminator_L

   !# Vertical Diffusion Winds (if <0,we remove REF)
   real :: Dcmip_nuZ_wd = 0.
   namelist /dcmip/ Dcmip_nuZ_wd

   !# Vertical Diffusion Theta (if <0,we remove REF)
   real :: Dcmip_nuZ_th = 0.
   namelist /dcmip/ Dcmip_nuZ_th

   !# Vertical Diffusion Tracers (if <0,we remove REF)
   real :: Dcmip_nuZ_tr = 0.
   namelist /dcmip/ Dcmip_nuZ_tr

   !# Earth's radius reduction factor
   real(kind=REAL64)  :: Dcmip_X = 1.d0
   namelist /dcmip/ Dcmip_X

   !# Reset in Dcmip_set as: Dcmip_vrd_L = Dcmip_nuZ_wd/=0.or.Dcmip_nuZ_tr/=0.or.Dcmip_nuZ_th/=0
   logical :: Dcmip_vrd_L = .false.

contains

!**s/r dcmip_nml - Read namelist dcmip

      integer function dcmip_nml (F_unf, F_dcmip_L)
      use lun
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      logical, intent(out) :: F_dcmip_L
      integer, intent(in) :: F_unf

      character(len=64) :: nml_S
      logical nml_must
!
!-------------------------------------------------------------------
!
! boiler plate - start

      F_dcmip_L = .false.

      if ( F_unf < 0 ) then
         dcmip_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=dcmip)
         return
      endif

      dcmip_nml= -1 ; nml_must= .false. ; nml_S= 'dcmip'
      rewind(F_unf)
      read (F_unf, nml=dcmip, end= 1001, err=1003)
      dcmip_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         dcmip_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      endif
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (dcmip_nml < 0 ) return
      if ((Lun_out>=0).and.(dcmip_nml==0)) write (Lun_out, 6004) trim(nml_S)
      dcmip_nml= 1
      F_dcmip_L = Dcmip_case > 0

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

!
!-------------------------------------------------------------------
!
      return
      end function dcmip_nml

!**s/r dcmip_set -  Setup for parameters DCMIP 2012/2016

      subroutine dcmip_set (F_adv_L,F_unout)

      use dcst
      use step_options
      use tdpack, only : rayt_8, omega_8
      use dyn_fisl_options
      use hvdif_options
      use ctrl
      use glb_ld
      use lun

      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>

      logical F_adv_L
      integer F_unout

      !object
      !=======================================
      !   Setup for parameters DCMIP 2012/2016
      !=======================================

      real(kind=REAL64) Rotation

      if (Dcmip_case==0) return

      !Identify 3D Advection runs
      !--------------------------
      F_adv_L = Dcmip_case>=11.and.Dcmip_case<=13
      if (F_unout>0.and.F_adv_L) write (F_unout, 7000)

      if (Ctrl_phyms_L.and.Dcmip_case/=162) call handle_error(-1,'SET_DCMIP','SET_DCMIP: Turn OFF GEM physics')

      if (Ctrl_phyms_L.and.Dcmip_case==162.and.(Dcmip_prec_type/=-1.or.Dcmip_pbl_type/=-1)) &
         call handle_error(-1,'SET_DCMIP','SET_DCMIP: T162: GEM physics + DCMIP2016 physics not allowed')

      !-----------------------------------
      !Set Earth's radius reduction factor
      !-----------------------------------

        if (Dcmip_case== 20.and.Dcmip_X/=1.   ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 20: Set Dcmip_X=1   ')
        if (Dcmip_case== 21.and.Dcmip_X/=500. ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 21: Set Dcmip_X=500 ')
        if (Dcmip_case== 22.and.Dcmip_X/=500. ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 22: Set Dcmip_X=500 ')
        if (Dcmip_case== 31.and.Dcmip_X/=125. ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 31: Set Dcmip_X=125 ')
        if (Dcmip_case== 43.and.Dcmip_X/=1.   ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 43: Set Dcmip_X=1   ')
        if (Dcmip_case==410.and.Dcmip_X/=1.   ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=410: Set Dcmip_X=1   ')
        if (Dcmip_case==411.and.Dcmip_X/=10.  ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=411: Set Dcmip_X=10  ')
        if (Dcmip_case==412.and.Dcmip_X/=100. ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=412: Set Dcmip_X=100 ')
        if (Dcmip_case==413.and.Dcmip_X/=1000.) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=413: Set Dcmip_X=1000')
        if (Dcmip_case==161.and.Dcmip_X/=1.   ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=161: Set Dcmip_X=1   ')
        if (Dcmip_case==162.and.Dcmip_X/=1.   ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=162: Set Dcmip_X=1   ')
        if (Dcmip_case==163.and.Dcmip_X/=120. ) call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case=163: Set Dcmip_X=120 ')

        !Reset Earth's radius
        !--------------------
        Dcst_rayt_8     = rayt_8/Dcmip_X ! rayt_8 = Reference Earth's Radius (m)
        Dcst_inv_rayt_8 = Dcmip_X/rayt_8

        !Reset time step
        !---------------
        Step_dt = Step_dt/Dcmip_X

        !Reset Vspng_coeftop
        !-------------------
        Vspng_coeftop = Vspng_coeftop/Dcmip_X

      !--------------------
      !Set Earth's rotation
      !--------------------

        Rotation = Dcmip_X

        !No Rotation
        !-----------
        if (Dcmip_case==163.or. &
            Dcmip_case== 20.or. &
            Dcmip_case== 21.or. &
            Dcmip_case== 22.or. &
            Dcmip_case== 31) Rotation = 0.

        !Reset Earth's angular velocity
        !------------------------------
        Dcst_omega_8 = Rotation * omega_8 ! omega_8 = Reference rotation rate of the Earth (s^-1)

      !------------------------
      !Toy Chemistry Terminator
      !------------------------
      if (Dcmip_case==161.and..NOT.Dcmip_Terminator_L) &
         call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 161: Set Dcmip_Terminator_L')

      !--------------
      !Riley Friction
      !--------------
      if ((Dcmip_case==21.or.Dcmip_case==22).and..NOT.Vspng_riley_L) &
         call handle_error(-1,'SET_DCMIP','SET_DCMIP Dcmip_case= 21/22: Set Vspng_riley_L')

      !------------------------------------------------------------------------------------
      !Vertical Diffusion: CAUTION: WE ASSUME COEFFICIENTS ARE ALREADY SET FOR SMALL PLANET
      !------------------------------------------------------------------------------------
      Dcmip_vrd_L = Dcmip_nuZ_wd/=0.or.Dcmip_nuZ_tr/=0.or.Dcmip_nuZ_th/=0

      !--------------------------
      !DCMIP 2016 Physics Package
      !--------------------------
      Dcmip_physics_L = Dcmip_prec_type/=-1.or.Dcmip_pbl_type/=-1

      if (Dcmip_pbl_type/=-1.and.Dcmip_case==163) &
         call handle_error(-1,'SET_DCMIP', 'SET_DCMIP: DONT activate Planetary Boundary Layer when Supercell')

      !--------------------------------------------------
      !Moist=1/Dry=0 Initial conditions (case=161/41X/43)
      !--------------------------------------------------
      if (Dcmip_moist==0.and.Dcmip_case==43) &
         call handle_error(-1,'SET_DCMIP', 'SET_DCMIP: Set Dcmip_moist==1 when Dcmip_case= 43' )
      if (Dcmip_moist==1.and.Dcmip_case>=410.and.Dcmip_case<=413) &
         call handle_error(-1,'SET_DCMIP', 'SET_DCMIP: Set Dcmip_moist==0 when Dcmip_case= 41X')

      if (Lun_out>0) write (Lun_out,1000)

      if (Lun_out>0) write (Lun_out,1001) Dcmip_X,Dcst_rayt_8,Dcst_omega_8,Dcmip_prec_type,Dcmip_pbl_type,Dcmip_Terminator_L,Dcmip_vrd_L

      return

 1000 format(                                                                    /, &
      '!----------------------------------------------------------------------|',/, &
      '!DESCRIPTION of DCMIP_2016_PHYSICS                                     |',/, &
      '!----------------------------------------------------------------------|',/, &
      '!  prec_type         | Type of precipitation/microphysics              |',/, &
      '!                    | ------------------------------------------------|',/, &
      '!                    |  0: Large-scale precipitation (Kessler)         |',/, &
      '!                    |  1: Large-scale precipitation (Reed-Jablonowski)|',/, &
      '!                    | -1: NONE                                        |',/, &
      '!----------------------------------------------------------------------|',/, &
      '!  pbl_type          | Type of planetary boundary layer                |',/, &
      '!                    | ------------------------------------------------|',/, &
      '!                    |  0: Reed-Jablonowski Boundary layer             |',/, &
      '!                    |  1: Georges Bryan Planetary Boundary Layer      |',/, &
      '!                    | -1: NONE                                        |',/, &
      '!----------------------------------------------------------------------|')

 1001 format( &
      /,'SETUP FOR PARAMETERS DCMIP_2016: (S/R DCMIP_2016_SET)',   &
      /,'=====================================================',/, &
        ' X Scaling Factor for Small planet  = ',F7.2          ,/, &
        ' Revised radius   for Small planet  = ',E14.5         ,/, &
        ' Revised angular velocity           = ',E14.5         ,/, &
        ' Precipitation/microphysics type    = ',I2            ,/, &
        ' Planetary boundary layer type      = ',I2            ,/, &
        ' Toy Chemistry                      = ',L2            ,/, &
        ' Vertical Diffusion                 = ',L2            ,   &
      /,'=====================================================',/)
 7000 format (//'  ====================='/&
                '  ACADEMIC 3D Advection'/&
                '  ====================='//)

      end subroutine dcmip_set

   function dcmip_options_init() result(F_istat)
      use, intrinsic :: iso_fortran_env
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.

      return
   end function dcmip_options_init

end module dcmip_options
