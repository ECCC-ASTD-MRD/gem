!---------------------------------- LICENCE BEGIN -------------------------------

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

!**s/r set_world_view

      subroutine set_world_view()
      use iso_c_binding

      use adz_options
      use coriolis
      use dcmip_options
      use dyn_expo_options
      use gem_options
      use hvdif_options
      use init_options
      use inp_options
      use lam_options
      use out_options
      use spn_options
      use step_options
      use theo_options

      use geomh
      use inp_mod
      use out_collector, only: block_collect_set, Bloc_me

      use glb_ld
      use lun
      use tr3d
      use out3
      use out_meta
      use path
      use outgrid
      use wb_itf_mod
      use ptopo
      use tracers_attributes_mod, only: tracers_attributes
      implicit none
#include <arch_specific.hf>

#include <rmnlib_basics.hf>
      include "rpn_comm.inc"

      integer, external :: adv_config, gemdm_config, set_io_pes,&
                           domain_decomp,sol_transpose

      character(len=50) :: LADATE
      integer :: istat,options,wload,hzd,monot,massc,unf
      integer :: f1,f2,f3,f4
      integer, dimension(20) :: err
      real :: vmin,vmax
      character(len=12) :: intp_S
!
!-------------------------------------------------------------------
!
      err(:) = 0
      err(1) = wb_put( 'model/Hgrid/is_yinyang',Grd_yinyang_L,&
                       WB_REWRITE_NONE+WB_IS_LOCAL )
      if (Grd_yinyang_L) then
         Path_ind_S=trim(Path_input_S)//'/MODEL_INPUT/'&
                                      //trim(Grd_yinyang_S)
         err(2) = wb_put( 'model/Hgrid/yysubgrid',Grd_yinyang_S,&
                          WB_REWRITE_NONE+WB_IS_LOCAL )
      else
         Path_ind_S=trim(Path_input_S)//'/MODEL_INPUT'
      end if
      Path_phy_S=trim(Path_input_S)//'/'

! Read namelists from file Path_nml_S
      unf= 0
      if (fnom (unf,Path_nml_S, 'SEQ+OLD', 0) == 0) then
         if (Lun_out >= 0) write (Lun_out, 6000) trim( Path_nml_S )
         err( 3) = theocases_nml   (unf)
         err( 4) = HORgrid_nml     (unf)
         err( 5) = VERgrid_nml     (unf)
         err( 6) = step_nml        (unf)
         err( 7) = dynKernel_nml   (unf)
         err( 8) = gem_nml         (unf)
         err( 9) = init_nml        (unf)
         err(10) = inp_nml         (unf)
         err(11) = lam_nml         (unf)
         err(12) = out_nml         (unf)
         err(13) = spn_nml         (unf)
         if ( Dynamics_Kernel_S(1:13) == 'DYNAMICS_FISL' ) then
            err(14) = dyn_fisl_nml (unf)
            err(15) = adz_nml      (unf)
            err(16) = hvdif_nml    (unf)
         else if ( Dynamics_Kernel_S(1:13) == 'DYNAMICS_EXPO' ) then
            err(14) = dyn_expo_nml (unf)
         else
            if (lun_out > 0) then
               write (lun_out,1010) trim(Dynamics_Kernel_S)
            end if
            err(1)= -1
         end if
         istat= fclos(unf)
      else
         if (Lun_out >= 0) write (Lun_out, 6001) trim( Path_nml_S )
         err(1)= -1
      end if

      call gem_error ( minval(err(:)),'set_world_view',&
                       'Error reading nml or with wb_put' )

! Read physics namelist

      call itf_phy_nml()

!     Setup for parameters DCMIP
!     --------------------------
      call dcmip_set (Ctrl_testcases_adv_L, Lun_out)

! Establish final configuration

      err(:) = 0
      err(1) = HORgrid_config (adz_maxcfl_fact)
      call theo_cfg() !must absolutely be done here
      err(2) = VERgrid_config ()

! Establish domain decomposition (mapping subdomains and processors)
      err(3) = domain_decomp (Ptopo_npex, Ptopo_npey, .false.)

! Initialize GMM
      call set_gmm()

! Final configuration
      if (minval(err(:))>=0) err(4) = gemdm_config()

      call canonical_cases ("SET_ZETA")

      call gem_error ( minval(err(:)),'CONFIGURATION ERROR', &
                      'ABORT in set_world_view' )

      err(1) = dynKernel_nml   (-1)
      err(1) = step_nml        (-1)
      err(1) = gem_nml         (-1)
      err(1) = HORgrid_nml     (-1)
      err(1) = VERgrid_nml     (-1)
      err(1) = init_nml        (-1)
      err(1) = inp_nml         (-1)
      err(1) = lam_nml         (-1)
      err(1) = out_nml         (-1)
      err(1) = spn_nml         (-1)
      if ( Dynamics_Kernel_S(1:13) == 'DYNAMICS_FISL' ) then
         err(1) = dyn_fisl_nml (-1)
         err(1) = adz_nml      (-1)
         err(1) = hvdif_nml    (-1)
      else if ( Dynamics_Kernel_S(1:13) == 'DYNAMICS_EXPO' ) then
         err(1) = dyn_expo_nml (-1)
      end if

      if (lun_out > 0) then
         f1 = G_ni/Ptopo_npex + min(1,mod(G_ni,Ptopo_npex))
         f2 = G_ni-f1*(Ptopo_npex-1)
         f3 = G_nj/Ptopo_npey + min(1,mod(G_nj,Ptopo_npey))
         f4 = G_nj-f3*(Ptopo_npey-1)
         write (lun_out,1001) Grd_typ_S,G_ni,G_nj,G_nk,f1,f3,f2,f4
         LADATE='RUNSTART='//Step_runstrt_S(1:8)//Step_runstrt_S(10:11)
         call write_status_file3 (trim(LADATE))
         call write_status_file3 ( 'communications_established=YES' )
         if (Grd_yinyang_L) then
            call write_status_file3 ('GEM_YINYANG=1')
         end if
      end if

! Master output PE for all none distributed components

      options = WB_REWRITE_NONE+WB_IS_LOCAL
      f1= 0
      istat = wb_put('model/outout/pe_master', f1,options)
      istat = min(wb_put('model/l_minx',l_minx,options),istat)
      istat = min(wb_put('model/l_maxx',l_maxx,options),istat)
      istat = min(wb_put('model/l_miny',l_miny,options),istat)
      istat = min(wb_put('model/l_maxy',l_maxy,options),istat)
      call gem_error ( istat,'set_world_view', &
                       'Problem with min-max wb_put')

! Establish a grid id for RPN_COMM package and obtain Out3_iome,Inp_iome

      Out3_npes= max(1,min(Out3_npes,min(Ptopo_npex,Ptopo_npey)**2))
      Inp_npes = max(1,min(Inp_npes ,min(Ptopo_npex,Ptopo_npey)**2))

      err= 0
      if (lun_out > 0) write (lun_out,1002) 'Output',Out3_npes
      err(1)= set_io_pes (Out3_comm_id,Out3_comm_setno,Out3_iome,&
                          Out3_comm_io,Out3_iobcast,Out3_npes)
      if (lun_out > 0) write (lun_out,1002) 'Input',Inp_npes
      err(2)= set_io_pes (Inp_comm_id ,Inp_comm_setno ,Inp_iome ,&
                          Inp_comm_io ,Inp_iobcast ,Inp_npes )
      call gem_error ( min(err(1),err(2)),'set_world_view', &
                       'IO pes config is invalid' )

      Out3_ezcoll_L= .true.
      if ( (Out3_npex > 0) .and. (Out3_npey > 0) ) then
         Out3_npex= min(Out3_npex,Ptopo_npex)
         Out3_npey= min(Out3_npey,Ptopo_npey)
         call block_collect_set ( Out3_npex, Out3_npey )
         Out3_npes= Out3_npex * Out3_npey
         Out3_iome= -1
         if (Bloc_me == 0) Out3_iome= 0
         Out3_ezcoll_L= .false.
      end if
      out_stk_size= Out3_npes*2

      istat = tracers_attributes( 'DEFAULT,'//trim(Tr3d_default_s), &
                                  wload, hzd, monot, massc, vmin, vmax, intp_S )

      if (Lun_out > 0) then
         if (trim(Tr3d_default_s)=='') then
            write (Lun_out,'(/a)') &
            ' SYSTEM DEFAULTS FOR TRACERS ATTRIBUTES:'
         else
            write (Lun_out,'(/a)') &
            ' USER DEFAULTS FOR TRACERS ATTRIBUTES:'
         end if
         write (Lun_out,2001)
         write (Lun_out,2002) wload,hzd,monot,massc,vmin,vmax,intp_S
      end if

      call heap_paint()

! Initialize geometry of the model
      call set_geomh

      if (trim(Dynamics_Kernel_S(1:13)) == 'DYNAMICS_FISL') then
         err(1)= sol_transpose ( Ptopo_npex, Ptopo_npey, .false. )
         call gem_error (err(1), 'SET_WORLD_VIEW', 'sol_transpose -- ABORTING')
      end if

      call set_coriolis_shallow ( geomh_x_8, geomh_y_8, geomh_xu_8, geomh_yv_8, &
                                  Grd_rot_8, l_minx, l_maxx, l_miny, l_maxy )

      call set_opr(' ')

      call set_params()

      call set_sor()

 1001 format (' GRID CONFIG: GRTYP=',a,5x,'GLB=(',i5,',',i5,',',i5,')    maxLCL(',i4,',',i4,')    minLCL(',i4,',',i4,')')
 1002 format (/ ' Creating IO pe set for ',a,' with ',i4,' Pes')
 1010 format (/' Invalid choice of dynamic kernel: ',a,/)
 2001 format ( ' DEFAULT tracers attributes:'/3x,'Wload  Hzd   Mono  Mass    Min        Max     Intp')
 2002 format (4i6,3x,2(e10.3,1x),a12)
 6000 format (' READING configuration namelists from FILE: '/A)
 6001 format (/,' Namelist FILE: ',A,' NOT AVAILABLE'/)
!
!-------------------------------------------------------------------
!
      return
      end
