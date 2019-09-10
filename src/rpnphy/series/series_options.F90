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

module series_options
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   implicit none
   public
   save

   !#TODO: once validated, use a single string len
   integer, parameter :: SER_STRLEN = 32
   integer, parameter :: SER_STRLEN_STN = 128
   integer, parameter :: SER_STRLEN_VAR = 4

   type SER_STN_LALO_T
      sequence
      character(SER_STRLEN_STN) :: name
      real :: lat,lon
   end type SER_STN_LALO_T

   type SER_STN_T
      sequence
      character(SER_STRLEN_STN) :: name
      real    :: lat, lon
      integer :: ig, jg, gidx  !# Full grid i,j; stn Full index
      integer :: il, jl        !# Lcl  grid i,j
      integer :: iphy, jphy    !# Lcl  grid i,j phy (folded) coor
   end type SER_STN_T

   integer, parameter :: NSTATMAX = 5000
   integer, parameter :: NVARMAX  = 1000
   integer, parameter :: NVARGEO  = 12
   integer, parameter :: STN_MISSING = -9999
   character(len=*), parameter :: STN_MISSING_S = 'UNDEFINED'
   character(len=*), parameter :: PKGNAME_S = '(series) '
   character(len=32), parameter :: VERSION_S = 'rpnphy_series_6'
   !#Note: keep P_serg_srgeo_s in sync with series_geop_mod
   character(len=SER_STRLEN_VAR), parameter :: P_serg_srgeo_s(NVARGEO) = &
        (/'MA', 'LA', 'LO', 'ZP', 'MG', 'LH', &
          'AL', 'SD', 'TM', 'TP', 'GL', 'HS'/)

   type(SER_STN_LALO_T), parameter :: SER_STN_LALO_DEFAULT = &
        SER_STN_LALO_T(STN_MISSING_S, real(STN_MISSING), real(STN_MISSING))

   type(SER_STN_T), parameter :: SER_STN_DEFAULT =  &
        SER_STN_T( &
        STN_MISSING_S, real(STN_MISSING), real(STN_MISSING), &
        STN_MISSING, STN_MISSING, 0, &
        STN_MISSING, STN_MISSING, &
        STN_MISSING, STN_MISSING)

   !---- Namelist vars

   !# List of time series for surface variables
   character(len=SER_STRLEN_VAR) :: P_serg_srsrf_s(NVARMAX) = ' '
   namelist /series/ P_serg_srsrf_s

   !# List of time series for profile variables
   character(len=SER_STRLEN_VAR) :: P_serg_srprf_s(NVARMAX) = ' '
   namelist /series/ P_serg_srprf_s

   !# Number of timesteps between time-series writeout
   integer :: P_serg_srwri = 1
   namelist /series/ P_serg_srwri
   integer :: series_out_nsteps = 1
   !#TODO: rename to series_xtr_interval_S
   !#TODO: new to series_write_interval_S (or factor of series_xtr_interval_S)

   !# Times series package stops at this timestep
   integer :: P_serg_serstp = huge(1)
   namelist /series/ P_serg_serstp
   !#TODO: rename to series_laststep

   !# Stations chosen in lat,lon for time-series
   !# Format: "STN1_NAME",lat1,lon1, "STN2_NAME",lat2,lon2, ...
   type(SER_STN_LALO_T) :: xst_stn_latlon(NSTATMAX)
   namelist /series/ xst_stn_latlon


   !---- set in ser_nml
   logical :: series_on_L = .false.

   type(SER_STN_T) :: series_stng(NSTATMAX)
   integer :: series_nstng = 0

   integer :: series_ngeo  = NVARGEO
   integer :: series_nsurf = 0
   integer :: series_nprof = 0

   !---- set in ser_init
   logical :: series_initok_L = .false.
   integer :: series_master_pe = 0
   integer :: series_moyhr = 0
   real    :: series_delt = 0.
   integer :: series_interval_sec = 0
   integer(INT64) :: series_jdate = 0
   integer :: series_cmcdateo14(14) = 0
   logical :: series_satuco_L = .false.
   logical :: series_satues_L = .false.

   integer :: series_nstnb = 0
   integer :: series_nstnl = 0
   integer :: series_stnl_idxb(NSTATMAX) = STN_MISSING
   type(SER_STN_T) :: series_stnb(NSTATMAX)

   real, pointer :: series_gdata(:,:) => NULL()
   real, pointer :: series_sdata(:,:,:) => NULL()
   real, pointer :: series_pdata(:,:,:,:) => NULL()
   real, pointer :: series_gdata_out(:,:) => NULL()
   real, pointer :: series_sdata_out(:,:,:) => NULL()
   real, pointer :: series_pdata_out(:,:,:,:) => NULL()
   integer, pointer :: series_gdone(:) => NULL()
   integer, pointer :: series_sdone(:) => NULL()
   integer, pointer :: series_pdone(:) => NULL()
   integer, pointer :: series_sdone2(:,:) => NULL()
   integer, pointer :: series_pdone2(:,:) => NULL()
   integer, pointer :: series_sdone_out(:) => NULL()
   integer, pointer :: series_pdone_out(:) => NULL()

   !---- set in ser_step
   logical :: series_end_L = .false.
   integer :: series_kount = 0
   real(REAL64) :: series_heure_8 = 0.d0

   !----
   logical :: series_paused_L = .false.

   integer :: series_fileid = -1    !# Output file unit nb

contains

   function series_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      integer :: i
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      P_serg_srprf_s(:) = ' '
      P_serg_srsrf_s(:) = ' '
      series_stnl_idxb(:) = STN_MISSING
      do i = 1, NSTATMAX
         xst_stn_latlon(i) =  SER_STN_LALO_DEFAULT
         series_stng(i)    =  SER_STN_DEFAULT
         series_stnb(i)    =  SER_STN_DEFAULT
      enddo
      return
   end function series_options_init

end module series_options
