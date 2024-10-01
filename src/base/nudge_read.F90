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

subroutine nudge_read (F_step_kount, F_Lctl_step)
   use iso_c_binding

   use timestr_mod, only: timestr2step
   use clib_itf_mod, only: clib_toupper, clib_tolower
   use rmn_gmm
   use input_mod, only: input_new, input_add, input_nbvar, input_set_basedir, &
        input_set_filename, input_setgridid, input_isvarstep, input_meta, input_get
   use inputio_mod, only: inputio_new, inputio_add, inputio_nbvar, &
        inputio_set_filename, inputio_nbvar, inputio_isvarstep, inputio_meta, &
        inputio_get, INPUT_FILES_ANAL, INPUTIO_T
   use mu_jdate_mod, only: jdate_from_cmc
   use ptopo_utils, only: ptopo_io_set, ptopo_iotype, PTOPO_BLOC, PTOPO_IODIST
   use statfld_dm_mod, only: statfld_dm
   use tdpack
   use vGrid_Descriptors, only: vgrid_descriptor, vgd_free
   use vgrid_wb, only: vgrid_wb_get, vgrid_wb_put
   
   use cstv
   use gem_options
   use spn_options
   use inp_options
   use glb_ld
   use glb_pil
   use HORgrid_options, only: Grd_global_gid, &
        Grd_glbcore_gid, Grd_yinyang_L
   use inp_mod, only: Inp_comm_id
   use path
   use step_options
   use var_gmm
   use VERgrid_options, only: VGRID_M_S, VGRID_T_S

   use omp_timing
   implicit none
   !@params
   integer :: F_step_kount, F_Lctl_step
   !@author 
   !@object
   !  Read fields for nudging

#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   include "rpn_comm.inc"

   !# The following parameter should be promoted to nml options like iau_* ones
   !# To turn off nudge_read, just set period(1) > period(2)
   character(len=32), parameter  :: nudg_period_S(2) = (/ &
        "0h   ", &
        "9999h"/)
   
   integer, parameter :: STATS_PRECISION = 8
   character(len=2), parameter :: NUDG_PREFIX='N_'
   character(len=6), parameter :: NUDG_FILE = 'NUDGE '
   character(len=32), parameter  :: NUDG_VGRID_M_S = 'nudg-m'
   character(len=32), parameter  :: NUDG_VGRID_T_S = 'nudg-t'
   character(len=32), parameter  :: NUDG_RFLD_S = 'NUDGRFLD:P'
   character(len=32), parameter  :: NUDG_RFLD_LS_S = 'NUDGRFLDLS:P'
   character(len=32), parameter  :: NUDG_ALTFLD_M_S = 'NUDGALTFLDM:P'
   character(len=32), parameter  :: NUDG_ALTFLD_T_S = 'NUDGALTFLDT:P'
   logical, parameter :: ALLOWDIR = .true.
   logical, parameter :: DO_TT2VT = .true.
   
   logical, save :: is_converted = .false.
   integer, save :: nudg_period(2) = (/-1.,-1./) !# start step, end step

   logical, save :: is_init_L = .false.
   integer, save :: nbvar = 0
   type(INPUTIO_T), save :: inputobj
   real, pointer, save :: vt(:,:,:) => null()

   real, pointer, dimension(:,:) :: nudg_rfld, pw_rfld
   real, pointer, dimension(:,:,:) :: nudg_altfld, pw_altfld
   character(len=256) :: incfg0_S, incfg_S, vgrid_S, msg_S
   character(len=32)  :: rfld_S, rfldls_S, altfld_M_S, altfld_T_S
   character(len=16)  :: iname0_S, iname1_S, datev_S
   integer :: istat, dateo, datev, nudg_vtime, step_freq, ivar, ni1, nj1, &
        i, j, k, n, nw, add, lijk(3), uijk(3), step_0, ipass
   integer(IDOUBLE) :: jdateo
   real, pointer, dimension(:,:,:) :: data0, data1
   real, pointer, dimension(:,:,:) :: myptr0, myptr1
   real, pointer, dimension(:,:) :: myptr2d
   type(gmm_metadata) :: mymeta
   type(vgrid_descriptor) :: vgridm, vgridt
   integer, pointer :: ip1list_m(:), ip1list_t(:), ip1listref(:)
   logical :: vert_updated_L
   
   !--------------------------------------------------------------------------
!!$   * le step number a partir de 0z = lctl_step
!!$   * le step number a partir du debut (-3h) = step_kount
!!$   * Cstv_dt_8 = dble(Step_dt)
!!$   * la date de validite?  - voir au début de gem_run.F90     
!!$      dayfrac = dble(Step_kount) * Cstv_dt_8 / sec_in_day
!!$      call incdatsd (datev,Step_runstrt_S,dayfrac)

   if (.not. Grd_yinyang_L .or. spn_freq < 0 ) return

   if (.not.is_converted) then
      is_converted = .true.
      istat = timestr2step(nudg_period(1), nudg_period_S(1), Cstv_dt_8)
      istat = timestr2step(nudg_period(2), nudg_period_S(2), Cstv_dt_8)
   endif
   
   call datp2f(dateo, Step_runstrt_S)
   nudg_vtime = F_step_kount*Cstv_dt_8
   call incdatr(datev, dateo, dble(nudg_vtime)/3600.d0)
   call datf2p(datev_S, datev)

   if (F_lctl_step < nudg_period(1) .or. F_lctl_step > nudg_period(2)) then
!!$      write(msg_S, '(i6,i6)') F_step_kount, F_lctl_step
!!$      call msg(MSG_INFO, ' No NUDGE_READ : '//trim(Step_runstrt_S)//' + '//trim(msg_S)//" = "//trim(datev_S))
      return
   endif
   write(msg_S, '(i6,i6)') F_step_kount, F_lctl_step
   call msg(MSG_INFO, ' NUDGE_READ : '//trim(Step_runstrt_S)//' + '//trim(msg_S)//" = "//trim(datev_S))

   call gtmg_start(50, 'NUDG', 1)

   ptopo_iotype = PTOPO_IODIST

   call nudge_read_init()

   vert_updated_L = .false.
   DO_IVAR:  do ivar = 1, nbvar
      istat = inputio_meta(inputobj%cfg, ivar, iname0_S, iname1_S)
      
      nullify(data0, data1)
      istat = gmm_get(NUDG_PREFIX//trim(iname0_S), data0)
      if (iname1_S /= '') istat = gmm_get(NUDG_PREFIX//trim(iname1_S), data1)
      if (.not.associated(data0) .or. &
           (iname1_S /= '' .and..not.associated(data1))) then
         call msg(MSG_ERROR, '(nudge_read) Problem getting ptr for: '// &
              trim(iname0_S)//' '//trim(iname1_S)//' -- Skipping')
         cycle DO_IVAR
      end if

      istat = inputio_isvarstep(inputobj, ivar, F_step_kount)
      IF_READ: if (RMN_IS_OK(istat) .or. F_Lctl_step == nudg_period(1)) then

         call nudge_read_update_vert(vert_updated_L)

         write(msg_S, '(i9,a,i6,i6,a)') F_step_kount*nint(Cstv_dt_8), 's [Step=',F_step_kount,F_lctl_step,']'
         call msg(MSG_INFO, '(nudge_read) Getting: '//trim(iname0_S)// &
              ' at t0+'//trim(msg_S))
                  
         !# Get interpolted data from file
         vgrid_S = NUDG_VGRID_T_S
         if (iname0_S == 'uu') vgrid_S = NUDG_VGRID_M_S
         nullify(myptr0, myptr1, myptr2d)
         istat = inputio_get(inputobj, ivar, F_step_kount, myptr0, myptr1, &
              F_vgrid_S=vgrid_S)
         if (.not.RMN_IS_OK(istat) .or. .not.associated(myptr0) .or. &
           (iname1_S /= '' .and..not.associated(myptr1))) then
            call msg(MSG_ERROR, '(nudge_read) Problem getting data for: '// &
                 trim(iname0_S)//' '//trim(iname1_S))
            istat = RMN_ERR
         end if
         call handle_error(istat, 'nudge_read', 'Problem getting data for: '// &
              trim(iname0_S)//' '//trim(iname1_S))

         !# Print input fields statistics
         if (Spn_yy_nudge_data_stats_L) then
            if (associated(myptr0)) &
                 call statfld_dm(myptr0, iname0_S, F_step_kount, 'nudge_read', STATS_PRECISION)
            if (associated(myptr1)) &
                 call statfld_dm(myptr1, iname1_S, F_step_kount, 'nudge_read', STATS_PRECISION)
         end if

         !# Move data to grid with halos; save in GMM
         if (associated(myptr0)) then
            data0(1:l_ni,1:l_nj,:) = myptr0(1:l_ni,1:l_nj,:)
            deallocate(myptr0, stat=istat)
         end if
         if (associated(myptr1) .and. associated(data1)) then
            data1(1:l_ni,1:l_nj,:) = myptr1(1:l_ni,1:l_nj,:)
            deallocate(myptr1, stat=istat)
         end if

         !# Adapt units and horizontal positioning
         lijk = lbound(data0) ; uijk = ubound(data0)
         if (iname0_S == 'uu') then
            data0 = data0 * knams_8
            if (associated(data1)) data1 = data1 * knams_8
         elseif (iname0_S == 'tt') then
            data0 = data0 + TCDK
         end if
         if (iname0_S == 'p0') data0 = 100.*data0

!!$         if (Grd_yinyang_L .and. (iname0_S == 'hu' .or. &
!!$              .not.any(iname0_S == Nudg_tracers_S))) then
         if (Grd_yinyang_L) then ! .and. iname0_S == 'hu') then

            if (iname1_S /= ' ' .and. associated(data1)) then
               call yyg_xchng_vec_q2q ( data0, data1, l_minx, l_maxx, l_miny, l_maxy, G_nk)
               ! Exchange halos in the pilot zone
               if (Glb_pilotcirc_L) then
                   call rpn_comm_propagate_pilot_circular(data0, &
                        lijk(1),uijk(1),lijk(2),uijk(2), &
                        l_ni,l_nj,uijk(3),Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
                   call rpn_comm_propagate_pilot_circular(data1, &
                        lijk(1),uijk(1),lijk(2),uijk(2), &
                        l_ni,l_nj,uijk(3),Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
               else
                 call rpn_comm_xch_halo(data0,lijk(1),uijk(1),lijk(2),uijk(2), &
                 l_ni,l_nj,uijk(3),G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
                 call rpn_comm_xch_halo(data1,lijk(1),uijk(1),lijk(2),uijk(2),&
                 l_ni,l_nj,uijk(3),G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
               endif
            else
               call yyg_int_xch_scal ( data0, uijk(3), .false., 'CUBIC', .true.)
            end if

         end if

      end if IF_READ
      
   end do DO_IVAR

   !# Convert tt 2 vt
   nullify(data0, data1)
   istat = gmm_get(NUDG_PREFIX//'tt', data0)
   istat = gmm_get(NUDG_PREFIX//'hu', data1)
   call mfottv2(data0, data0, data1, &
                l_minx, l_maxx, l_miny, l_maxy, G_nk, &
                1-G_halox, l_ni+G_halox, 1-G_haloy, l_nj+G_haloy, DO_TT2VT)
!!$   ! or 
!!$   call mfottvh2(tt, tv, qq, qh, &
!!$                 l_minx, l_maxx, l_miny, l_maxy, G_nk, &
!!$                 1-G_halox, l_ni+G_halox, 1-G_haloy, l_nj+G_haloy, DO_TT2VT)
   
   call gtmg_stop(50)
   !--------------------------------------------------------------------------
   return

contains

   subroutine nudge_read_init()
      
      if (is_init_L) return
      is_init_L = .true.

      !# Set up input module
      jdateo = jdate_from_cmc(dateo)
      istat = ptopo_io_set(Inp_npes)
      istat = inputio_new(inputobj, jdateo, nint(Cstv_dt_8), " ", &
           Path_input_S, Grd_global_gid, Grd_glbcore_gid, Inp_comm_id, &
           F_li0=1, F_lin=l_ni, F_lj0=1, F_ljn=l_nj, F_iotype=ptopo_iotype)
      if (RMN_IS_OK(istat)) then
         istat = inputio_set_filename(inputobj, NUDG_FILE, NUDG_FILE, &
              ALLOWDIR, INPUT_FILES_ANAL)
      end if
      step_0 = nudg_period(1) - F_lctl_step + F_step_kount
      if (spn_yy_nudge_data_tint == 'near-cst-p') then
         step_freq = nint(max(spn_freq, spn_yy_nudge_data_freq)/Cstv_dt_8)
         step_0 = step_0 + step_freq/2  !#??? +1 ???
      else
         step_freq = nint(spn_freq/Cstv_dt_8)
      endif
      
      write(incfg0_S, '(a,i4,a,i4,a,a,a)') 'freq=', step_0, ',', step_freq, &
           '; search=', trim(NUDG_FILE), &
           '; tinterp='//spn_yy_nudge_data_tint(1:4)//'; hinterp=cubic; vinterp=c-cond; levels=-1'

      incfg_S = 'in=TT;'//trim(incfg0_S)
      call msg(MSG_INFO, '(nudge_read) add input: '//trim(incfg_S))
      istat = min(istat, inputio_add(inputobj%cfg, incfg_S))
      
      incfg_S = 'in=HU;'//trim(incfg0_S)
      call msg(MSG_INFO, '(nudge_read) add input: '//trim(incfg_S))
      istat = min(istat, inputio_add(inputobj%cfg, incfg_S))
      
      incfg_S = 'in=UU; IN2=VV;'//trim(incfg0_S)
      call msg(MSG_INFO, '(nudge_read) add input: '//trim(incfg_S))
      istat = min(istat, inputio_add(inputobj%cfg, incfg_S))
      
!!$      write(incfg0_S, '(a,i4,a,i4,a,a,a)') 'freq=', step_0, ',', step_freq, &
!!$           '; search=', trim(NUDG_FILE), &
!!$           '; typvar=P; tinterp=linear; hinterp=cubic'
!!$      
!!$      incfg_S = 'in=P0;'//trim(incfg0_S)
!!$      call msg(MSG_INFO, '(nudge_read) add input: '//trim(incfg_S))
!!$      istat = min(istat, inputio_add(inputobj%cfg, incfg_S))
      
      nbvar = inputio_nbvar(inputobj)
      call handle_error_l(RMN_IS_OK(istat).and.nbvar>0, 'nudge_read', &
           'Problem initializing the input module')

      !# Create data space to save inc values between read-incr
      DO_IVAR0: do ivar = 1, nbvar
         istat = inputio_meta(inputobj%cfg, ivar, iname0_S, iname1_S)
         nullify(data0, data1)
         if (iname0_S == 'p0') then
            mymeta = meta2d ; mymeta%l(3) = gmm_layout(1,1,0,0,1)
            istat = gmm_create(NUDG_PREFIX//trim(iname0_S), data0, &
                 mymeta, GMM_FLAG_IZER)
         else
            istat = gmm_create(NUDG_PREFIX//trim(iname0_S), data0, &
                 meta3d_nk, GMM_FLAG_IZER)
            if (iname1_S /= ' ') then
               istat = gmm_create(NUDG_PREFIX//trim(iname1_S), data1, &
                                  meta3d_nk, GMM_FLAG_IZER)
            end if
         end if
      end do DO_IVAR0
      
      !# define a vert coor with ref on l_ni/j
      nullify(ip1list_m, ip1list_t)
      istat = vgrid_wb_get(VGRID_M_S, vgridm, ip1list_m, F_sfcfld_S=rfld_S, &
           F_sfcfld2_S=rfldls_S, F_altfld_S=altfld_M_S)
      istat = vgrid_wb_get(VGRID_T_S, vgridt, ip1list_t, F_altfld_S=altfld_T_S)

      mymeta = GMM_NULL_METADATA
      mymeta%l(1) = gmm_layout(1,l_ni,0,0,l_ni)
      mymeta%l(2) = gmm_layout(1,l_nj,0,0,l_nj)
      
      if (rfld_S /= '') rfld_S = NUDG_RFLD_S
      if (rfldls_S /= '') rfldls_S = NUDG_RFLD_LS_S
      if (altfld_M_S /= '') altfld_M_S = NUDG_ALTFLD_M_S
      if (altfld_T_S /= '') altfld_T_S = NUDG_ALTFLD_T_S
      ip1listref => ip1list_m(1:G_nk)
      istat = vgrid_wb_put(NUDG_VGRID_M_S, vgridm, ip1listref, rfld_S, &
           rfldls_S, F_overwrite_L=.true., F_altfld_S=altfld_M_S)
      ip1listref => ip1list_t(1:G_nk)
      istat = vgrid_wb_put(NUDG_VGRID_T_S, vgridt, ip1listref, rfld_S, &
           rfldls_S, F_overwrite_L=.true., F_altfld_S=altfld_T_S)
     
      if (rfld_S /= '') then
         nullify(nudg_rfld)
         istat = gmm_create(NUDG_RFLD_S, nudg_rfld, mymeta)
         if (rfldls_S /= '') then
            nullify(nudg_rfld)
            istat = gmm_create(NUDG_RFLD_LS_S, nudg_rfld, mymeta)
         end if
      endif
      if (altfld_M_S /= '') then
         mymeta%l(3) = gmm_layout(1,G_nk,0,0,G_nk)
         nullify(nudg_altfld)
         istat = gmm_create(NUDG_ALTFLD_M_S, nudg_altfld, mymeta)
      endif
      if (altfld_T_S /= '') then
         mymeta%l(3) = gmm_layout(1,G_nk,0,0,G_nk)
         nullify(nudg_altfld)
         istat = gmm_create(NUDG_ALTFLD_T_S, nudg_altfld, mymeta)
      endif
      
      istat = vgd_free(vgridm)
      istat = vgd_free(vgridt)
      if (associated(ip1list_m)) deallocate(ip1list_m, stat=istat)
      if (associated(ip1list_t)) deallocate(ip1list_t, stat=istat)
      
      return
   end subroutine nudge_read_init

   
   subroutine nudge_read_update_vert(vert_updated_L)
      logical, intent(inout) :: vert_updated_L
      integer :: istat

      if (vert_updated_L) return
      vert_updated_L = .true.
      
      !# Update reference surface field for vgrid
      istat = vgrid_wb_get(VGRID_M_S, vgridm, F_sfcfld_S=rfld_S, &
           F_sfcfld2_S=rfldls_S, F_altfld_S=altfld_M_S)
      istat = vgd_free(vgridm)
      istat = vgrid_wb_get(VGRID_T_S, vgridt, F_altfld_S=altfld_T_S)
      istat = vgd_free(vgridt)

      if (rfld_S /= '') then 
         nullify(pw_rfld, nudg_rfld)
         istat = gmm_get(rfld_S, pw_rfld)
         istat = gmm_get(NUDG_RFLD_S, nudg_rfld)
         if (associated(nudg_rfld) .and. associated(pw_rfld)) then
            nudg_rfld(:,:) = pw_rfld(1:l_ni,1:l_nj)
         end if

         if (rfldls_S /= '') then
            nullify(pw_rfld, nudg_rfld)
            istat = gmm_get(rfldls_S, pw_rfld)
            istat = gmm_get(NUDG_RFLD_LS_S, nudg_rfld)
            if (associated(nudg_rfld) .and. associated(pw_rfld)) then
               nudg_rfld(:,:) = pw_rfld(1:l_ni,1:l_nj)
            end if
         end if
      end if

      if (altfld_M_S /= '') then
         nullify(pw_altfld, nudg_altfld)
         istat = gmm_get(altfld_M_S, pw_altfld)
         istat = gmm_get(NUDG_ALTFLD_M_S, nudg_altfld)
         if (associated(pw_altfld) .and. associated(nudg_altfld)) then
            nudg_altfld(:,:,:) = pw_altfld(1:l_ni,1:l_nj,1:G_nk)
         endif
      end if
      if (altfld_T_S /= '') then
         nullify(pw_altfld, nudg_altfld)
         istat = gmm_get(altfld_T_S, pw_altfld)
         istat = gmm_get(NUDG_ALTFLD_T_S, nudg_altfld)
         if (associated(pw_altfld) .and. associated(nudg_altfld)) then
            nudg_altfld(:,:,:) = pw_altfld(1:l_ni,1:l_nj,1:G_nk)
         endif
      end if
      return
   end subroutine nudge_read_update_vert
   
end subroutine nudge_read
