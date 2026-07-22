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

      integer function cplocn_init (F_path_S, F_print_L, F_unout, F_dateo, &
                                    F_z0mtype, F_z0lat, F_z0ttype)
use iso_c_binding
use rmn_gmm
use cpl_mod
use cplocn_mod
use rpn_comm
use phygridmap, only: phy_yinyang_L, phy_yinyang_S
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_path_S,F_z0mtype,F_z0ttype
      logical, intent(in)          :: F_print_L
      integer, intent(in)          :: F_dateo, F_unout
      real   , intent(in)          :: F_z0lat(2)

!authors    Francois Roy -- spring 2015
! 
!revision
! v4_80 - Roy, F.  - initial version
! v5.0b1 - JM Belanger (dec 2017) - z0ttype

!Purpose
! Initialize coupling with ocean

#include <rmn/WhiteBoard.hf>

      integer ier,err
      integer :: unf, ios, nrec

      character(len=255) :: cplaofnml

      logical :: l_lam
!
!     ---------------------------------------------------------------
!
      cplocn_init=-1

      if (F_print_L) write(F_unout,2000)

      err= RPN_COMM_mype (cplocn_myproc, cplocn_mycol, cplocn_myrow)
      if ( err < 0 ) then
        write(F_unout,9900) 'PROBLEMS WHEN CALLING RPN_COMM_mype'
        call flush(F_unout)
        goto 995
      endif
      call RPN_COMM_size ('GRID', cplocn_nprocs, err)
      if (err < 0) then
        write(F_unout,9900) 'PROBLEMS WHEN CALLING RPN_COMM_size'
        call flush(F_unout)
        goto 995
      endif

      err = 0

      ! read cplao_settings.nml
      cplaofnml = trim(F_path_S)//'/MODEL_INPUT/cplao_settings.nml'
      call read_cplao_nml( cplaofnml, F_unout, F_print_L, F_z0mtype, F_z0lat, F_z0ttype )

      ! set start datestamp of coupling
      call datf2p (cplocn_runstrt_S,F_dateo)

      l_lam = .true.  ! LU or yinyang only now

      cplocn_it = 0
      if ( abs( mod( real(cplao_dt), cpl_drv_delt )) > 1e-5) then
         write(F_unout,*) 'cplao_dt is not dividable by model time step', cplao_dt, cpl_drv_delt
         err = -1
         goto 995
      endif
      cplocn_rap_dt = nint( real(cplao_dt) / cpl_drv_delt )
      if ( cplocn_rap_dt > 1 .and. cplao_xchg_mode == 0 ) then
         write(F_unout,*) 'cplao_dt =', cplao_dt, 'cpl_drv_delt =', cpl_drv_delt
         write(F_unout,*) 'cplao_dt has to be the same value as cpl_drv_delt for cplao_xchg_mode = 0'
         err = -1
         goto 995
      endif

      ! prepare some arrays
      call cplocn_init_busou
      call cplocn_init_busin( F_print_L, F_unout )

      ! prepare iris exchange
      call cpl_declare_grid_information()
      call cpl_declare_exchanged_fields()

      if (F_print_L) write(F_unout,2001)

      cplocn_init = 1
      return
 995  call handle_error(err,'cplocn_init','Problems with cplocn_init')

 2000 format( &
      /,'INITIALIZATION OF COUPLING INTERFACE S/R cplocn_init', &
      /,'=====================================================')
 2001 format( &
      /,'INITIALIZATION OF COUPLING INTERFACE ENDED S/R cplocn_init', &
      /,'===========================================================')
 9900 format (/,1x,a)
!
!     ---------------------------------------------------------------
!
      return
      end function cplocn_init


      subroutine cplocn_init_busou
      use cpl_mod
      use cplocn_mod
      implicit none
      ! locals
      integer :: ivar, nk, istat

!     ________________________________________________________________

      nk=cpl_drv_gnk-1

      allocate ( ocn_busou(cpl_drv_lni,cpl_drv_lnj,cplocn_n_fldou) )
      ocn_busou = 0.
      cplocn_it = 0

      !* FLUX COUPLING OCEAN BUS OUT
      !* Ice-Ocean independant
      !* 1- FB  - SW down                          p
      !* 2- FI  - LW down                          p
      !* 3- RT  - Precipitation                    p
      !* Ocean model flux calculation
      !* 4- TT  - Air temperature                  d
      !* 5- UU  - Wind x component                 d
      !* 6- VV  - Wind y component                 d
      !* 7- QA  - Specific humidity                d
      !* 8- PX  - First momentum level pressure    d
      !* 9- PX  - First thermo   level pressure    d
      !*10- P0  - Ground level pressure            d
      !* Atmospheric model flux calculated (DE-ACTIVATED)
      !*11- SHO - Sensible heat flux over water    ?
      !*12- SHI - Sensible heat flux over ice      ?
      !*13- LHO - Latent heat flux over water      ?
      !*14- LHI - Latent heat flux over ice        ?
      !*15- TXO - Wind stress x component (water)  ?
      !*16- TYO - Wind stress y component (water)  ?
      !*17- TXI - Wind stress x component (ice)    ?
      !*18- TYI - Wind stress y component (ice)    ?

      DO ivar = 1, cplocn_n_fldou

         SELECT CASE ( cplocn_cvou_S(ivar) )

         CASE ( 'FBA' )

           cplocn_cvou_N(ivar) = 'flusolis'
           cplocn_cvou_G(ivar) = 'P'
           cplocn_cvou_K(:,ivar) = (/ 1, 1 /)

         CASE ( 'FIA' )

           cplocn_cvou_N(ivar) = 'fdsi'
           cplocn_cvou_G(ivar) = 'P'
           cplocn_cvou_K(:,ivar) = (/ 1, 1 /)

         CASE ( 'RTA' )

           cplocn_cvou_N(ivar) = 'rt'
           cplocn_cvou_G(ivar) = 'P'
           cplocn_cvou_K(:,ivar) = (/ 1, 1 /)

         CASE ( 'TTA' )

           cplocn_cvou_N(ivar) = 'PW_TT:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ nk, nk /)

         CASE ( 'UUA' )

           cplocn_cvou_N(ivar) = 'PW_UU:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ nk, nk /)

         CASE ( 'VVA' )

           cplocn_cvou_N(ivar) = 'PW_VV:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ nk, nk /)

         CASE ( 'QQA' )

           cplocn_cvou_N(ivar) = 'TR/HU:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ nk, nk /)

         CASE ( 'PMA' )

           cplocn_cvou_N(ivar) = 'PW_PM:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ nk, nk /)

         CASE ( 'PTA' )

           cplocn_cvou_N(ivar) = 'PW_PT:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ nk, nk /)

         CASE ( 'P0A' )

           cplocn_cvou_N(ivar) = 'PW_P0:P'
           cplocn_cvou_G(ivar) = 'D'
           cplocn_cvou_K(:,ivar) = (/ 1, 1 /)

         CASE DEFAULT

           istat=-1
           call handle_error(istat,'cplocn_init_busou','WRONG NAME TAG')

         END SELECT
      ENDDO
!     ________________________________________________________________
!
      end subroutine cplocn_init_busou


      subroutine cplocn_init_busin( F_print_L, F_unout )
      use iso_c_binding
      use rmn_gmm
      use cpl_mod
      use cplocn_mod
      use, intrinsic :: iso_fortran_env
      implicit none

      logical, intent(in)          :: F_print_L
      integer, intent(in)          :: F_unout

      ! locals
      integer :: ivar, ier
      type(gmm_metadata) :: meta3d_ocn_busin, meta2d_surf0


      gmmk_ocn_busin_s = 'ocn_busin'
      gmmk_gli_0_s     = 'gli_0'
      gmmk_i8i_0_s     = 'i8i_0'
      gmmk_sdi_0_s     = 'sdi_0'
      gmmk_tmo_0_s     = 'tmo_0'

      nullify(ocn_busin,gli_0,i8i_0,sdi_0,tmo_0)
      call gmm_build_meta3D(meta3d_ocn_busin,                    &
                            1,cpl_drv_lni,0,0,cpl_drv_lni,       &
                            1,cpl_drv_lnj,0,0,cpl_drv_lnj,       &
                            1,cplocn_n_fldin,0,0,cplocn_n_fldin, &
                            0,GMM_NULL_FLAGS)
      call gmm_build_meta2D(meta2d_surf0,                        &
                            1,cpl_drv_lni,0,0,cpl_drv_lni,       &
                            1,cpl_drv_lnj,0,0,cpl_drv_lnj,       &
                            0,GMM_NULL_FLAGS)

      ier = gmm_create(gmmk_ocn_busin_s , ocn_busin , meta3d_ocn_busin , GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_gli_0_s ,     gli_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_i8i_0_s ,     i8i_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_sdi_0_s ,     sdi_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)
      ier = gmm_create(gmmk_tmo_0_s ,     tmo_0 ,     meta2d_surf0 ,     GMM_FLAG_RSTR)

      ier = gmm_get(gmmk_ocn_busin_s , ocn_busin)
      ier = gmm_get(gmmk_gli_0_s ,     gli_0)
      ier = gmm_get(gmmk_i8i_0_s ,     i8i_0)
      ier = gmm_get(gmmk_sdi_0_s ,     sdi_0)
      ier = gmm_get(gmmk_tmo_0_s ,     tmo_0)

      if ( .NOT. cpl_rstn_L ) then
         ocn_busin(:,:,:) = 0.
         gli_0(:,:) = cplocn_missval
         i8i_0(:,:) = cplocn_missval
         sdi_0(:,:) = cplocn_missval
         tmo_0(:,:) = cplocn_missval
       else
         if (F_print_L) write (F_unout,*) 'CPLOCN_INIT: RESTART MODE'
      endif

      do ivar=1,cplocn_n_fldin

         SELECT CASE ( cplocn_cvin_S(ivar) )

         CASE ( 'MCP' ) ; icvin_MCP=ivar
         CASE ( 'ALO' ) ; icvin_ALO=ivar
         CASE ( 'ALI' ) ; icvin_ALI=ivar
         CASE ( 'T4O' ) ; icvin_T4O=ivar
         CASE ( 'T4I' ) ; icvin_T4I=ivar
         CASE ( 'SHO' ) ; icvin_SHO=ivar
         CASE ( 'SHI' ) ; icvin_SHI=ivar
         CASE ( 'LHO' ) ; icvin_LHO=ivar
         CASE ( 'LHI' ) ; icvin_LHI=ivar
         CASE ( 'TXO' ) ; icvin_TXO=ivar
         CASE ( 'TYO' ) ; icvin_TYO=ivar
         CASE ( 'TXI' ) ; icvin_TXI=ivar
         CASE ( 'TYI' ) ; icvin_TYI=ivar
         CASE ( 'ZTO' ) ; icvin_ZTO=ivar
         CASE ( 'ZQO' ) ; icvin_ZQO=ivar
         CASE ( 'ZUO' ) ; icvin_ZUO=ivar
         CASE ( 'ZVO' ) ; icvin_ZVO=ivar
         CASE ( 'ZTI' ) ; icvin_ZTI=ivar
         CASE ( 'ZQI' ) ; icvin_ZQI=ivar
         CASE ( 'ZUI' ) ; icvin_ZUI=ivar
         CASE ( 'ZVI' ) ; icvin_ZVI=ivar
         CASE ( 'GLI' ) ; icvin_GLI=ivar
         CASE ( 'I8I' ) ; icvin_I8I=ivar
         CASE ( 'SDI' ) ; icvin_SDI=ivar
         CASE ( 'TMO' ) ; icvin_TMO=ivar
         CASE ( 'UUO' ) ; icvin_UUO=ivar
         CASE ( 'VVO' ) ; icvin_VVO=ivar
         CASE ( 'I7I' ) ; icvin_I7I=ivar
         CASE ( 'UUI' ) ; icvin_UUI=ivar
         CASE ( 'VVI' ) ; icvin_VVI=ivar
         CASE ( 'QSO' ) ; icvin_QSO=ivar
         CASE ( 'QSI' ) ; icvin_QSI=ivar
         CASE ( 'ILO' ) ; icvin_ILO=ivar
         CASE ( 'ILI' ) ; icvin_ILI=ivar
         CASE ( 'ZMO' ) ; icvin_ZMO=ivar
         CASE ( 'ZMI' ) ; icvin_ZMI=ivar
         CASE ( 'ZHO' ) ; icvin_ZHO=ivar
         CASE ( 'ZHI' ) ; icvin_ZHI=ivar

         CASE DEFAULT

           ier=-1
           call handle_error(ier,'cplocn_init_busin','WRONG NAME TAG')

         END SELECT

      enddo

      end subroutine cplocn_init_busin


      subroutine read_cplao_nml (F_nmlf_S, F_unout, F_comproc_L, &
                                 F_z0mtype, F_z0lat, F_z0ttype)
      use cpl_mod, only: cpl_drv_delt
      use cplocn_mod
      implicit none

      character(len=*), intent(in) :: F_nmlf_S
      integer, intent(in)          :: F_unout
      logical, intent(in)          :: F_comproc_L
      character(len=*), intent(in) :: F_z0mtype,F_z0ttype
      real   , intent(in)          :: F_z0lat(2)
      integer  fnom
      external fnom

      integer unf,nrec

      character(len=16) :: z0mtype,z0ttype
      real              :: z0tlat(2), z0tlat_w(2)


      namelist /cplao_cfgs/ z0mtype, z0ttype, z0tlat

      namelist /cplao_step/   cplao_dt, cplao_xchg_mode

!
!-------------------------------------------------------------------
!
! Defaults values are taken from arguments if not provided in the namelist
!
      cplao_dt    = cpl_drv_delt
      z0tlat_w(:) = F_z0lat(:)*180./3.14159265
      z0tlat(:)   = z0tlat_w(:)
      z0mtype     = TRIM(F_z0mtype)
      z0ttype     = TRIM(F_z0ttype)
      cplao_xchg_mode = 0
!
      if ((F_nmlf_S.eq.'print').or.(F_nmlf_S.eq.'PRINT')) then
         if (F_unout.ge.0.and.F_comproc_L) write (F_unout,nml=cplao_cfgs)
         return
      elseif (F_nmlf_S .ne. '') then
!
         unf = 0
         if (fnom (unf,F_nmlf_S, 'SEQ+OLD', nrec) .ne. 0) goto 9110

         rewind(unf)
         read (unf, nml=cplao_cfgs, end = 9120, err=9120)
         if(F_comproc_L) write( F_unout, cplao_cfgs )

         rewind(unf)
         read (unf, nml=cplao_step, end = 9121, err=9121)
         if(F_comproc_L) write( F_unout, cplao_step )

         call fclos (unf)

      endif
      !
      ! check consistency of GEM/CPL bulk formula
      !
      if ( ABS(z0tlat_w(1)-z0tlat(1)) > 0.01 .or. &
           ABS(z0tlat_w(2)-z0tlat(2)) > 0.01 ) then
         write (F_unout,9900) 'INCONSISTENT Z0TLAT VALUES GEM VS CPL'
         write (F_unout,*) 'GEM Z0TLAT(1)=',z0tlat_w(1)
         write (F_unout,*) 'GEM Z0TLAT(2)=',z0tlat_w(2)
         write (F_unout,*) 'CPL Z0TLAT(1)=',z0tlat(1)
         write (F_unout,*) 'CPL Z0TLAT(2)=',z0tlat(2)
         write (F_unout, 8000)
         goto 9998
      endif

      if (trim(z0mtype) /= trim(F_z0mtype)) then
         write (F_unout,9900) 'INCONSISTENT z0mtype VALUES GEM VS CPL'
         write (F_unout,*) 'GEM z0mtype  =',F_z0mtype
         write (F_unout,*) 'CPL z0mtype  =',z0mtype
         write (F_unout, 8000)
         goto 9998
      endif

      if (trim(z0ttype) /= trim(F_z0ttype)) then
         write (F_unout,9900) 'INCONSISTENT z0ttype VALUES GEM VS CPL'
         write (F_unout,*) 'GEM z0ttype  =',F_z0ttype
         write (F_unout,*) 'CPL z0ttype  =',z0ttype
         write (F_unout, 8000)
         goto 9998
      endif

      ! if all good go to return
      goto 9999
!
 9110 if (F_comproc_L.and.F_unout.ge.0) write (F_unout, 9050) trim( F_nmlf_S )
      if (F_comproc_L.and.F_unout.ge.0) write (F_unout, 8000)
      goto 9998
!
 9120 call fclos (unf)
      if (F_comproc_L.and.F_unout.ge.0) write (F_unout, 9150) 'cplao_cfgs',trim( F_nmlf_S )
      if (F_comproc_L.and.F_unout.ge.0) write (F_unout, 8000)
      goto 9998
!
 9121 call fclos (unf)
      if (F_comproc_L.and.F_unout.ge.0) write (F_unout, 9150) 'cplao_step',trim( F_nmlf_S )
      if (F_comproc_L.and.F_unout.ge.0) write (F_unout, 8000)
      goto 9998
!
 8000 format (/,'========= ABORT IN S/R read_cplao_nml.ft90 ============='/)
 9050 format (/,' FILE: ',A,' NOT AVAILABLE'/)
 9150 format (/,' NAMELIST ',A,' INVALID IN FILE: ',A/)
 9900 format (/,1x,a)
!
!-------------------------------------------------------------------
!
 9998  call handle_error(-1,'read_cplao_nml','Problems with read_cplao_nml')
 9999 return
      end subroutine read_cplao_nml
