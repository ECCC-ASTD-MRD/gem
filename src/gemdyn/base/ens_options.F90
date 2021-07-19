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
module ens_options
   use wb_itf_mod
   use ens_spp, only: MAX_NSPP
   implicit none
   public
   save

   ! Internal configurations
   integer, parameter, public :: MAX2DC=6

   ! Internal variables
   integer, public :: nchains
   logical, public :: spp_L

   !# Switch to activate generation of Markov chains, use of SKEB
   !# and use of PTP
   logical :: Ens_conf = .false.
   namelist /ensembles/ Ens_conf

   !# Switch to activate recycling  of Markov chains parameters
   !# in PTP & SKEB
   logical :: Ens_recycle_mc = .false.
   namelist /ensembles/ Ens_recycle_mc

   !# Switch to activate SKEB
   !# (3D MARKOV CHAINES)
   logical :: Ens_skeb_conf = .false.
   namelist /ensembles/ Ens_skeb_conf

   !# switch to print global stat related to Markov chains, SKEB and PTP
   !# (3D MARKOV CHAINES)
   logical :: Ens_stat = .false.
   namelist /ensembles/ Ens_stat

   !# switch to do the calculation of the divergence due to SKEB forcing
   !# (3D MARKOV CHAINES)
   logical :: Ens_skeb_div = .false.
   namelist /ensembles/ Ens_skeb_div

   !# switch to do SKEB calculation based on diffusion
   !# (3D MARKOV CHAINES)
   logical :: Ens_skeb_dif = .false.
   namelist /ensembles/ Ens_skeb_dif

   !# switch to do SKEB calculation based on gravity wave drag
   !# (3D MARKOV CHAINES)
   logical :: Ens_skeb_gwd = .false.
   namelist /ensembles/ Ens_skeb_gwd

   !# Seed of the random number generator usually we put DAY and member number
   !# (3D MARKOV CHAINES)
   integer :: Ens_mc_seed = -1
   namelist /ensembles/ Ens_mc_seed

   !# number of longitudes of the gaussian grid used for the 3D Markov chains
   !# (in the SKEB calculation)
   !# (3D MARKOV CHAINES)
   integer :: Ens_skeb_nlon = 16
   namelist /ensembles/ Ens_skeb_nlon

   !# number of latitudes of the gaussian grid used  for the 3D Markov chains
   !# (used in the SKEB calculation)
   !# (3D MARKOV CHAINES)
   integer :: Ens_skeb_nlat = 8
   namelist /ensembles/ Ens_skeb_nlat

   !#
   !# (3D MARKOV CHAINES)
   integer :: Ens_skeb_ncha = 1
   namelist /ensembles/ Ens_skeb_ncha

   !# low wave number truncation limit used in 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   integer :: Ens_skeb_trnl = 2
   namelist /ensembles/ Ens_skeb_trnl

   !# high wave number truncation limit used in 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   integer :: Ens_skeb_trnh = 8
   namelist /ensembles/ Ens_skeb_trnh

   !# maximum value of the 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_max = 0.
   namelist /ensembles/ Ens_skeb_max

   !# minimum value of the 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_min = 0.
   namelist /ensembles/ Ens_skeb_min

   !# std. dev. value for the 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_std = 0.
   namelist /ensembles/ Ens_skeb_std

   !# decorrelation time (seconds) for 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_tau = 0.
   namelist /ensembles/ Ens_skeb_tau

   !# value of stretch for 3D Markov chain (used by SKEB)
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_str = 0.
   namelist /ensembles/ Ens_skeb_str

   !# coefficient Alpha for momentum in SKEB
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_alph = 0.
   namelist /ensembles/ Ens_skeb_alph

   !# coefficient Alpha for temperature in SKEB
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_alpt = 0.
   namelist /ensembles/ Ens_skeb_alpt

   !# coefficient for Gaussian filter used in SKEB
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_bfc = 1.0e-01
   namelist /ensembles/ Ens_skeb_bfc

   !# wavelength for Gaussian filter in SKEB
   !# (3D MARKOV CHAINES)
   real :: Ens_skeb_lam = 2.0e+05
   namelist /ensembles/ Ens_skeb_lam

   !# no. of longitudes for 2D Markov chains
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_nlon(MAX2DC) = 16
   namelist /ensembles/Ens_ptp_nlon

   !# no. of latitudes for 2D Markov chains
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_nlat(MAX2DC) = 8
   namelist /ensembles/ Ens_ptp_nlat

   !# number of 2d Markov chains
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_ncha = 0
   namelist /ensembles/ Ens_ptp_ncha

   !# low wave number horizontal truncation limit for 2D Markov chains
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_trnl(MAX2DC) = 1
   namelist /ensembles/ Ens_ptp_trnl

   !# high wave number horizontal truncation limit for 2D Markov chains
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_trnh(MAX2DC) = 8
   namelist /ensembles/ Ens_ptp_trnh

   !# (ignored) Ens_ptp_l = Ens_ptp_trnh-Ens_ptp_trnl+1
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_l(MAX2DC) = 0
   namelist /ensembles/ Ens_ptp_l

   !# (ignored) Ens_ptp_m = Ens_ptp_trnh+1
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_m(MAX2DC) = 0
   namelist /ensembles/ Ens_ptp_m

   !# (ignored) Ens_ptp_lmax = maxval(Ens_ptp_l)
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_lmax = 0
   namelist /ensembles/ Ens_ptp_lmax

   !# (ignored) Ens_ptp_mmax = maxval(Ens_ptp_m)
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_mmax = 0
   namelist /ensembles/ Ens_ptp_mmax

   !# minimum value of the 2D Markov chain
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_min(MAX2DC) = 0.0
   namelist /ensembles/ Ens_ptp_min

   !# maximum value of the 2D Markov chains
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_max(MAX2DC) = 0.0
   namelist /ensembles/ Ens_ptp_max

   !# (ignored) Ens_ptp_latmax = maxval(Ens_ptp_nlat)
   !# (2D MARKOV CHAINES)
   integer :: Ens_ptp_latmax = 0
   namelist /ensembles/ Ens_ptp_latmax

   !# standard deviation value for 2D Markov chains
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_std(MAX2DC) = 0.0
   namelist /ensembles/ Ens_ptp_std

   !# decorrelation time (seconds) for 2D Markov chains
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_tau(MAX2DC) = 0.0
   namelist /ensembles/ Ens_ptp_tau

   !# value of stretch for Markov chains
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_str(MAX2DC) = 0.0
   namelist /ensembles/ Ens_ptp_str

   !# switch to activate PTP (perturb tendencies of physics)
   !# (2D MARKOV CHAINES)
   logical :: Ens_ptp_conf = .false.
   namelist /ensembles/ Ens_ptp_conf

   !# upper value of transition zone of vertical envelope in sigma for PTP
   !# (above that full perturbation)
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_env_u = 1.0
   namelist /ensembles/ Ens_ptp_env_u

   !# bottom value of transition zone of vertical envelope in sigma for PTP
   !# (below that no perturbation)
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_env_b = 1.0
   namelist /ensembles/ Ens_ptp_env_b

   !# CAPE value in Kain-Fritsch scheme to stop perturbing the physical
   !# tendencies
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_cape = 0.0
   namelist /ensembles/ Ens_ptp_cape

   !# vertical velocity value (m/s) above which we stop perturbing the
   !# physical tendencies
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_crit_w = 100.0
   namelist /ensembles/ Ens_ptp_crit_w

   !# factor of reduction of the perturbation the physical tendencies (in PTP)
   !# when convection occurs
   !# (2D MARKOV CHAINES)
   real :: Ens_ptp_fac_reduc = 0.0
   namelist /ensembles/ Ens_ptp_fac_reduc

   !# Activate SPP for specified parameters
   character(len=32), dimension(MAX_NSPP) :: Ens_spplist = ''
   namelist /ensembles/ Ens_spplist

   !# Set cutoff latitude (degrees from equator) for adv_rhsint SPP perturbations
   real :: Ens_spp_rhsint_lat = 45.
   namelist /ensembles/ Ens_spp_rhsint_lat

contains

!**s/r ens_nml - read parametres for ensemble forecast

      integer function ens_nml (F_unf)
      use ens_param
      use HORgrid_options
      use lun
      implicit none
#include <arch_specific.hf>

      integer F_unf

      character(len=64) :: nml_S
      logical nml_must, ptp_L
      integer ier,i
!
!--------------------------------------------------------------------
!
! boiler plate - start
      if ( F_unf < 0 ) then
         ens_nml= 0
         if ( Lun_out >= 0) write (Lun_out,nml=ensembles)
         return
      end if

      ens_nml= -1 ; nml_must= .false. ; nml_S= 'ensembles'

      rewind(F_unf)
      read (F_unf, nml=ensembles, end= 1001, err=1003)
      ens_nml= 0 ; goto 1000
 1001 if (Lun_out >= 0) write (Lun_out, 6005) trim(nml_S)
      if (.not.nml_must) then
         ens_nml= 1
         if (Lun_out >= 0) write (Lun_out, 6002) trim(nml_S)
      end if
      goto 1000
 1003 if (Lun_out >= 0) write (Lun_out, 6007) trim(nml_S)

 1000 if (ens_nml < 0 ) return

 6002 format (' Skipping reading of namelist ',A)
 6004 format (' Reading of namelist ',A,' is successful')
 6005 format (' Namelist ',A,' NOT AVAILABLE')
 6007 format (/,' NAMELIST ',A,' IS INVALID'/)
! boiler plate - end

      ptp_L = .false.

      if (ens_nml == 0) then

         if (Lun_out>=0) write (Lun_out, 6004) trim(nml_S)

         if (Ens_ptp_conf) then
            if (Ens_ptp_ncha < 1) then
               if (Lun_out >= 0) write(Lun_out,"(a)") 'Ens_nml: Ens_ptp_ncha must be >= 1 when PTP is activated with Ens_ptp_conf'
               ens_nml = -1
            endif
         else
            Ens_ptp_ncha = 0
         endif

         if ((Grd_typ_S /= 'GU'.and.Grd_typ_S /= 'GY').and.Ens_skeb_conf) then
            if(Lun_out >= 0) write(Lun_out,"(a,a)" ) 'Ens_nml: ','Ens_skeb only available with grid GU/GY'
            ens_nml = -1
         end if

         if ((Ens_mc_seed < 0))then
            if(Lun_out >= 0)write(Lun_out,*)'You have to provide a positive integer as seed see Ens_mc_seed in NAMELIST'
            ens_nml = -1
          end if

            if (Ens_skeb_nlon /= 2*Ens_skeb_nlat)then
              if(Lun_out >= 0)write(Lun_out,*)' Nlon must equal 2*nlat'
              ens_nml = -1
            end if

         if (Ens_ptp_ncha > MAX2DC) then
            if(Lun_out >= 0)write(Lun_out,*)'Ens_ptp_ncha must be <=9'
            ens_nml = -1
         end if

         ens_ptp_latmax=-1
         do i=1,Ens_ptp_ncha
            if (Ens_ptp_nlon(i) /= 2*Ens_ptp_nlat(i))then
                if(Lun_out >= 0)write(Lun_out,*)'Nlon2 must equal 2*nlat2'
                ens_nml = -1
            end if
            Ens_ptp_latmax=max(Ens_ptp_nlat(i),ens_ptp_latmax)
         enddo

         if (ens_nml < 0) return

         Ens_skeb_conf  =  Ens_skeb_conf.and.Ens_conf
         Ens_skeb_l     =  Ens_skeb_trnh-Ens_skeb_trnl+1
         Ens_skeb_m     =  Ens_skeb_trnh+1
         Ens_skeb_div   =  Ens_skeb_div .and.Ens_conf
         Ens_stat       =  Ens_stat.and.Ens_conf
         Ens_ptp_l      =  Ens_ptp_trnh-Ens_ptp_trnl+1
         Ens_ptp_lmax   =  maxval(Ens_ptp_l)
         Ens_ptp_m      =  Ens_ptp_trnh+1
         Ens_ptp_mmax   =  maxval(Ens_ptp_m)
         ptp_L          =  Ens_conf .and. Ens_ptp_conf
         spp_L          =  Ens_conf .and. len_trim(Ens_spplist(1)) > 0

         if(Lun_out >= 0) then
            write(Lun_out,"(a,i8)" )'Ens_mc_seed   = ',Ens_mc_seed
            write(Lun_out,'(a,l5)' )'Ens_skeb_conf  = ',Ens_skeb_conf
            write(Lun_out,'(a,l5)' )'Ens_stat  = ',Ens_stat
            write(Lun_out,'(a,l5)' )'Ens_skeb_div   = ',Ens_skeb_div
            write(Lun_out,'(a,10i5)')'Ens_ptp_l     = ',Ens_ptp_l
            write(Lun_out,'(a,l5)' )'Ens_recycle_mc  = ',Ens_recycle_mc
            write(Lun_out,'(a,i5)' )'Ens_ptp_lmax = ',Ens_ptp_lmax
            write(Lun_out,'(a,10i5)')'Ens_ptp_m     = ',Ens_ptp_m
            write(Lun_out,'(a,i5)' )'Ens_ptp_mmax = ',Ens_ptp_mmax
            write(Lun_out,'(a,10i5)')'Ens_skeb_l     = ',Ens_skeb_l
            write(Lun_out,'(a,10i5)')'Ens_skeb_m     = ',Ens_skeb_m
            write(Lun_out,'(a,l5)' )'Ens_ptp_L = ',ptp_L
            write(Lun_out,'(a,l5)' )'Ens_spp_L      = ',spp_L
            write(Lun_out,'(a,i5)' )'Ens_ptp_nc     = ',Ens_ptp_ncha
            write(Lun_out,'(a,f8.5)' )'Ens_ens_ptp_env_u = ',Ens_ptp_env_u
            write(Lun_out,'(a,f8.5)' )'Ens_ens_ptp_env_b = ',Ens_ptp_env_b
            write(Lun_out,'(a,f8.5)' )'Ens_ens_ptp_cape = ',Ens_ptp_cape
!!$            write(Lun_out,'(a,f8.5)' )'Ens_ens_ptp_tlc = ',Ens_ptp_tlc
            write(Lun_out,'(a,f12.5)' )'Ens_ens_ptp_crit_w = ',Ens_ptp_crit_w
            write(Lun_out,'(a,f8.5)' )'Ens_ens_ptp_fac_reduc = ',Ens_ptp_fac_reduc
         end if

      end if

      ier= WB_OK
      if (ptp_L .or. spp_L) then
         ier= min(wb_put('ens/PTP_NC'     , Ens_ptp_ncha     , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/PTP'        , ptp_L            , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/SPP'        , spp_L            , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/PTPENVU'    , Ens_ptp_env_u    , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/PTPENVB'    , Ens_ptp_env_b    , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/PTPCAPE'    , Ens_ptp_cape     , WB_REWRITE_MANY),ier)
!!$         ier= min(wb_put('ens/PTPTLC'     , Ens_ptp_tlc      , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/PTPCRITW'   , Ens_ptp_crit_w   , WB_REWRITE_MANY),ier)
         ier= min(wb_put('ens/PTPFACREDUC', Ens_ptp_fac_reduc, WB_REWRITE_MANY),ier)
      end if

      if (.not.WB_IS_OK(ier)) ens_nml= -1

      if (ens_nml == 1) then
         ier= min(wb_put('ens/PTP'        , Ens_ptp_conf     , WB_REWRITE_MANY),ier)
      end if
!
!--------------------------------------------------------------------
!
      return
      end

   function ens_options_init() result(F_istat)
      implicit none
      integer :: F_istat
#include <rmnlib_basics.hf>
      logical, save :: init_L = .false.
      F_istat = RMN_OK
      if (init_L) return
      init_L = .true.
      return
   end function ens_options_init

end module ens_options
