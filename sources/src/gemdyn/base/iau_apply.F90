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

subroutine iau_apply2(F_kount)
   use iso_c_binding

   use clib_itf_mod, only: clib_toupper, clib_tolower
   use gmm_itf_mod
   use input_mod, only: input_new, input_add, input_nbvar, input_set_basedir, &
        input_set_filename, input_setgridid, input_isvarstep, input_meta, input_get
   use inputio_mod, only: inputio_new, inputio_add, inputio_nbvar, &
        inputio_set_filename, inputio_nbvar, inputio_isvarstep, inputio_meta, &
        inputio_get, INPUT_FILES_ANAL, INPUTIO_T
   use mu_jdate_mod, only: jdate_from_cmc
   use ptopo_utils, only: ptopo_io_set, ptopo_iotype, PTOPO_BLOC, PTOPO_IO
   use statfld_dm_mod, only: statfld_dm
   use tdpack
   use vGrid_Descriptors, only: vgrid_descriptor, vgd_free
   use vgrid_wb, only: vgrid_wb_get, vgrid_wb_put

   use cstv
   use gem_options
   use init_options
   use inp_options
   use gmm_vt1
   use gmm_pw
   use glb_ld
   use glb_pil
   use HORgrid_options, only: Grd_local_gid, Grd_lclcore_gid, Grd_global_gid, &
        Grd_glbcore_gid, Grd_yinyang_L
   use inp_mod, only: Inp_comm_id
   use path
   use step_options
   use var_gmm

   use gem_timing
   implicit none
   !@params
   integer, intent(in) :: F_kount !step_kound
   !@author
   !       R. McTaggart-Cowan - Summer 2013
   !@revision
   !       S.Chamberland, 2014-07: use input_mod, allow vinterp
   !       S.Chamberland, 2017-09: use inputio_mod
   !@object
   !  Add an analysis increments to the model state (IAU).

#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>
   include "rpn_comm.inc"

   integer, parameter :: STATS_PRECISION = 8
   character(len=2), parameter :: IAU_PREFIX='I_'
   character(len=6), parameter :: IAU_FILE = 'IAUREP'
   character(len=32), parameter  :: VGRID_M_S = 'ref-m'
   character(len=32), parameter  :: VGRID_T_S = 'ref-t'
   character(len=32), parameter  :: IAU_VGRID_M_S = 'iau-m'
   character(len=32), parameter  :: IAU_VGRID_T_S = 'iau-t'
   character(len=32), parameter  :: IAU_REFP0_S = 'IAUREFP0:P'
   character(len=32), parameter  :: IAU_REFP0_LS_S = 'IAUREFP0LS:P'
   logical, parameter :: UVSCAL2WGRID = .false.
   logical, parameter :: ALLOWDIR = .true.
   logical, parameter :: DO_TT2TV = .true.
   logical, save :: is_init_L = .false.
   integer, save :: inputid = -1
   integer, save :: nbvar = 0
   integer, save :: step_freq2 = 0
   integer, save :: kount = 0
   integer, save :: rpncomm_gridid = -1
   type(INPUTIO_T), save :: inputobj
   real, pointer, contiguous, save :: weight(:)
   real, pointer, dimension(:,:), contiguous :: refp0, pw_p0
   character(len=256) :: incfg_S, vgrid_S, msg_S
   character(len=32)  :: refp0_S, refp0ls_S
   character(len=16)  :: iname0_S, iname1_S, datev_S
   integer :: istat, dateo, datev, iau_vtime, step_freq, ivar, ni1, nj1, &
        i, j, n, nw, add, lijk(3), uijk(3), step_0
   integer(IDOUBLE) :: jdateo
   real, pointer, dimension(:,:,:), contiguous :: data0, data1
   real, pointer, dimension(:,:,:) :: myptr0, myptr1
   real, pointer, dimension(:,:), contiguous :: myptr2d
   type(gmm_metadata) :: mymeta
   type(vgrid_descriptor) :: vgridm, vgridt
   integer, pointer, contiguous :: ip1list(:), ip1listref(:)
   logical :: input_use_old_l
   integer :: input_iotype
   !--------------------------------------------------------------------------
!!$   write(msg_S,'(l,i4,a,i7,a,i7)') (Cstv_dt_8*F_kount > Iau_period .or. Iau_interval<=0.),F_kount,'; t=',nint(Cstv_dt_8*F_kount),'; p=',nint(Iau_period)
!!$   call msg(MSG_INFO,'IAU YES/NO?: '//trim(msg_S))

   if (Cstv_dt_8*F_kount > Iau_period .or. Iau_interval<=0.) return
   call gemtime_start(50, 'IAU', 1)

   input_use_old_l = .true.
   input_iotype = PTOPO_BLOC
   if (Iau_input_type_S == 'DIST') then
      input_use_old_l = .false.
      input_iotype = PTOPO_IO
   else if (Iau_input_type_S == 'BLOC') then
      input_use_old_l = .false.
      input_iotype = PTOPO_BLOC
   end if
   ptopo_iotype = input_iotype

   if (input_use_old_l .or. input_iotype == PTOPO_BLOC) then
        istat = rpn_comm_bloc(Iau_ninblocx, Iau_ninblocy)
   end if

   call datp2f(dateo,Step_runstrt_S)
   iau_vtime = -Step_delay*Cstv_dt_8 + Iau_interval &
        * nint((Lctl_step)*Cstv_dt_8/Iau_interval-epsilon(1.))
   call incdatr(datev,dateo,dble(iau_vtime)/3600.d0)
   call datf2p(datev_S,datev)

   IF_INIT: if (.not.is_init_L) then
      is_init_L = .true.

      !# Set up input module
      OLD: if (input_use_old_l) then
         inputid = input_new(dateo, nint(Cstv_dt_8))
         istat = input_setgridid(inputid, Grd_lclcore_gid)
         istat = min(input_set_basedir(inputid, Path_input_S), inputid)
         istat = min(input_set_filename(inputid, IAU_FILE, IAU_FILE, &
              ALLOWDIR, INPUT_FILES_ANAL), istat)
      else !OLD
         jdateo = jdate_from_cmc(dateo)
         if (input_iotype == PTOPO_BLOC) then
            istat = inputio_new(inputobj, jdateo, nint(Cstv_dt_8), " ", &
                 Path_input_S, Grd_local_gid, Grd_lclcore_gid, rpncomm_gridid, &
                 F_li0=1, F_lin=l_ni, F_lj0=1, F_ljn=l_nj, F_iotype=input_iotype)
         else
            istat = ptopo_io_set(Inp_npes)
            istat = inputio_new(inputobj, jdateo, nint(Cstv_dt_8), " ", &
                 Path_input_S, Grd_global_gid, Grd_glbcore_gid, Inp_comm_id, &
                 F_li0=1, F_lin=l_ni, F_lj0=1, F_ljn=l_nj, F_iotype=input_iotype)
         end if
         if (RMN_IS_OK(istat)) then
            istat = inputio_set_filename(inputobj, IAU_FILE, IAU_FILE, &
                                         ALLOWDIR, INPUT_FILES_ANAL)
         end if
      end if OLD
      step_freq  = nint(Iau_interval/Cstv_dt_8)
      step_freq2 = nint((Iau_interval/2)/Cstv_dt_8) !TODO: check
      step_0 = mod((nint(Iau_period)/nint(Iau_interval)), 2)*step_freq2
      write(incfg_S, '(a,i4,a,i4,a,a,a)') 'freq=', step_0, ',', step_freq, &
           '; search=', trim(IAU_FILE), &
           '; typvar=R; hinterp=cubic; vinterp=c-cond'
      call msg(MSG_INFO, &
           '(iau_apply) add input: in=TT; levels=-1;'//trim(incfg_S))
      if (input_use_old_l) then
         istat = min(istat, &
              input_add(inputid, 'in=TT; levels=-1;'//trim(incfg_S)))
      else
         istat = min(istat, &
              inputio_add(inputobj%cfg, 'in=TT; levels=-1;'//trim(incfg_S)))
      end if
      call msg(MSG_INFO, &
           '(iau_apply) add input: in=HU; levels=-1;'//trim(incfg_S))
      if (input_use_old_l) then
         istat = min(istat, &
              input_add(inputid,'in=HU; levels=-1;'//trim(incfg_S)))
      else
         istat = min(istat, &
              inputio_add(inputobj%cfg, 'in=HU; levels=-1;'//trim(incfg_S)))
      end if
      call msg(MSG_INFO, &
           '(iau_apply) add input: in=UU; IN2=VV; levels=-1;'//trim(incfg_S))
      if (input_use_old_l) then
         istat = min(istat, &
              input_add(inputid, 'in=UU; IN2=VV; levels=-1;'//trim(incfg_S)))
      else
         istat = min(istat, &
              inputio_add(inputobj%cfg, 'in=UU; IN2=VV; levels=-1;'// &
              &           trim(incfg_S)))
      end if
      write(incfg_S, '(a,i4,a,i4,a,a,a)') 'freq=', step_0, ',', step_freq, &
           '; search=', trim(IAU_FILE), &
           '; typvar=R; hinterp=cubic'
      call msg(MSG_INFO, '(iau_apply) add input: in=P0;'//trim(incfg_S))
      if (input_use_old_l) then
         istat = min(input_add(inputid, 'in=P0; '//trim(incfg_S)), istat)
      else
         istat = min(inputio_add(inputobj%cfg, 'in=P0; '//trim(incfg_S)), istat)
      end if
      ivar = 1
      do while (len_trim(Iau_tracers_S(ivar)) > 0)
         istat = min(clib_tolower(Iau_tracers_S(ivar)), istat)
         if (Iau_tracers_S(ivar) /= 'hu') then
            call msg(MSG_INFO, &
                 '(iau_apply) add input: in='//trim(Iau_tracers_S(ivar))// &
                 '; levels=-1;'//trim(incfg_S))
            if (input_use_old_l) then
               istat = min(istat, &
                    input_add(inputid, 'in='//trim(Iau_tracers_S(ivar))// &
                    &         '; levels=-1;'//trim(incfg_S)))
            else
               istat = min(istat, &
                    inputio_add(inputobj%cfg, 'in='//trim(Iau_tracers_S(ivar))//&
                    &           '; levels=-1;'//trim(incfg_S)))
            end if
         end if
         ivar = ivar+1
         if (ivar > size(Iau_tracers_S)) exit
      end do
      if (input_use_old_l) then
         nbvar = input_nbvar(inputid)
      else
         nbvar = inputio_nbvar(inputobj)
      end if
      call handle_error_l(RMN_IS_OK(istat).and.nbvar>0, 'iau_apply', &
           'Problem initializing the input module')

      !# Create data space to save inc values between read-incr
      DO_IVAR0: do ivar = 1, nbvar
         if (input_use_old_l) then
            istat = input_meta(inputid, ivar, iname0_S, iname1_S)
         else
            istat = inputio_meta(inputobj%cfg, ivar, iname0_S, iname1_S)
         end if
         nullify(data0, data1)
         if (iname0_S == 'p0') then
            mymeta = meta2d ; mymeta%l(3) = gmm_layout(1,1,0,0,1)
            istat = gmm_create(IAU_PREFIX//trim(iname0_S), data0, &
                 mymeta, GMM_FLAG_IZER)
         else
            istat = gmm_create(IAU_PREFIX//trim(iname0_S), data0, &
                 meta3d_nk, GMM_FLAG_IZER)
            if (iname1_S /= ' ') then
               istat = gmm_create(IAU_PREFIX//trim(iname1_S), data1, &
                                  meta3d_nk, GMM_FLAG_IZER)
            end if
         end if
      end do DO_IVAR0

      !# Precompute filter coefficients on initialization
      nw = nint(Iau_period/Cstv_dt_8)
      add = 0; if (mod(nw, 2) == 0) add = 1
      nw = nw+add
      allocate(weight(nw))
      call msg(MSG_INFO, '(iau_apply) Precompute filter coefficients - '// &
           trim(Iau_weight_S))
      istat = clib_tolower(Iau_weight_S)
      select case (Iau_weight_S)
      case ('constant')
         weight = Cstv_dt_8 / Iau_period
      case ('spike')
         weight = 0.
         weight(nw/2) = 1.
      case ('sin')
         call handle_error_l(Iau_cutoff>0., 'iau_apply', &
              'Cutoff period must be greater than 0')
         j = 0 ; i = nw ; n = nw/2
         do while (j < nw)
            i = j-n; j = j+1
            if (i == 0) then
               weight(j) = 2.*Cstv_dt_8/(Iau_cutoff*3600.)
            else
               weight(j) = sin(i*pi_8/(n+1)) / (i*pi_8/(n+1)) * &
                    sin(i*(2.*pi_8*Cstv_dt_8/(Iau_cutoff*3600.))) / &
                    (i*pi_8)
            end if
         end do
         weight = weight/sum(weight(1:size(weight)-add))
      case default
         call handle_error(RMN_ERR, 'iau_apply', &
              'Unknown Iau_weight_S='//trim(Iau_weight_S))
      end select

      !# define a vert coor with ref on l_ni/j
      nullify(ip1list)
      istat = vgrid_wb_get(VGRID_M_S, vgridm, ip1list, F_sfcfld2_S=refp0ls_S)
      if (refp0ls_S /= '') refp0ls_S = IAU_REFP0_LS_S
      ip1listref => ip1list(1:G_nk)
      istat = vgrid_wb_put(IAU_VGRID_M_S, vgridm, ip1listref, IAU_REFP0_S, &
              refp0ls_S, F_overwrite_L=.true.)
      istat = vgd_free(vgridm)
      if (associated(ip1list)) deallocate(ip1list,stat=istat)
      nullify(ip1list)
      istat = vgrid_wb_get(VGRID_T_S, vgridt, ip1list, F_sfcfld2_S=refp0ls_S)
      if (refp0ls_S /= '') refp0ls_S = IAU_REFP0_LS_S
      ip1listref => ip1list(1:G_nk)
      istat = vgrid_wb_put(IAU_VGRID_T_S, vgridt, ip1listref, IAU_REFP0_S, &
                           refp0ls_S, F_overwrite_L=.true.)
      istat = vgd_free(vgridt)
      if (associated(ip1list)) deallocate(ip1list,stat=istat)
      mymeta = GMM_NULL_METADATA
      mymeta%l(1) = gmm_layout(1,l_ni,0,0,l_ni)
      mymeta%l(2) = gmm_layout(1,l_nj,0,0,l_nj)
      nullify(refp0)
      istat = gmm_create(IAU_REFP0_S, refp0, mymeta)
      if (refp0ls_S /= '') then
         nullify(refp0)
         istat = gmm_create(IAU_REFP0_LS_S, refp0, mymeta)
      end if

   end if IF_INIT

   !# Update reference surface field for vgrid
   istat = vgrid_wb_get(VGRID_M_S, vgridm, F_sfcfld_S=refp0_S, &
                        F_sfcfld2_S=refp0ls_S)
   istat = vgd_free(vgridm)

   nullify(pw_p0, refp0)
   istat = gmm_get(refp0_S, pw_p0)
   istat = gmm_get(IAU_REFP0_S, refp0)
   if (associated(refp0) .and. associated(pw_p0)) then
      refp0(:,:) = pw_p0(1:l_ni,1:l_nj)
   end if

   if (refp0ls_S /= '') then
      nullify(pw_p0, refp0)
      istat = gmm_get(refp0ls_S, pw_p0)
      istat = gmm_get(IAU_REFP0_LS_S, refp0)
      if (associated(refp0) .and. associated(pw_p0)) then
         refp0(:,:) = pw_p0(1:l_ni,1:l_nj)
      end if
   end if

   if (F_kount > 0) kount = F_kount+step_freq2-1
   DO_IVAR: do ivar = 1, nbvar
      if (input_use_old_l) then
         istat = input_meta(inputid, ivar, iname0_S, iname1_S)
      else
         istat = inputio_meta(inputobj%cfg, ivar, iname0_S, iname1_S)
      end if
      nullify(data0, data1)
      istat = gmm_get(IAU_PREFIX//trim(iname0_S), data0)
      if (iname1_S /= '') istat = gmm_get(IAU_PREFIX//trim(iname1_S), data1)
      if (.not.associated(data0) .or. &
           (iname1_S /= '' .and..not.associated(data1))) then
         call msg(MSG_ERROR, '(iau_apply) Problem getting ptr for: '// &
              trim(iname0_S)//' '//trim(iname1_S)//' -- Skipping')
         cycle DO_IVAR
      end if

      if (input_use_old_l) then
         istat = input_isvarstep(inputid, ivar, kount)
      else
         istat = inputio_isvarstep(inputobj, ivar, kount)
      end if
      IF_READ: if (RMN_IS_OK(istat)) then

         write(msg_S, '(i7)') kount*nint(Cstv_dt_8)
         call msg(MSG_INFO, '(iau_apply) Reading: '//trim(iname0_S)// &
              ' at t0+'//trim(msg_S)//'s')

         !# Get interpolted data from file
         vgrid_S = IAU_VGRID_T_S
         if (iname0_S == 'uu') vgrid_S = IAU_VGRID_M_S
         nullify(myptr0, myptr1, myptr2d)
         if (input_use_old_l) then
            istat = input_get(inputid, ivar, kount, Grd_local_gid, vgrid_S, &
                 myptr0, myptr1)
         else
            istat = inputio_get(inputobj, ivar, kount, myptr0, myptr1, &
                 F_vgrid_S=vgrid_S)
         end if
         if (.not.RMN_IS_OK(istat) .or. .not.associated(myptr0) .or. &
           (iname1_S /= '' .and..not.associated(myptr1))) then
            call msg(MSG_ERROR, '(iau_apply) Problem getting data for: '// &
                 trim(iname0_S)//' '//trim(iname1_S))
            istat = RMN_ERR
         end if
         call handle_error(istat, 'iau_apply', 'Problem getting data for: '// &
              trim(iname0_S)//' '//trim(iname1_S))

         !# Print input fields statistics
         if (iau_stats_L) then
            if (associated(myptr0)) &
                 call statfld_dm(myptr0, iname0_S, kount, 'iau_apply', STATS_PRECISION)
            if (associated(myptr1)) &
                 call statfld_dm(myptr1, iname1_S, kount, 'iau_apply', STATS_PRECISION)
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
         end if
         if (iname0_S == 'p0') data0 = 100.*data0

         if (Grd_yinyang_L .and. (iname0_S == 'hu' .or. &
              .not.any(iname0_S == Iau_tracers_S))) then

            if (iname1_S /= ' ' .and. associated(data1)) then
               call yyg_xchng_vec_q2q ( data0, data1, l_minx, l_maxx, l_miny, l_maxy, G_nk)
               ! Exchange halos in the pilot zone
               call rpn_comm_propagate_pilot_circular(data0,lijk(1),uijk(1),lijk(2),uijk(2), &
                                     l_ni,l_nj,uijk(3),Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
               call rpn_comm_propagate_pilot_circular(data1,lijk(1),uijk(1),lijk(2),uijk(2), &
                                     l_ni,l_nj,uijk(3),Glb_pil_e,Glb_pil_s,G_halox,G_haloy)
            else
               call yyg_xchng ( data0 ,lijk(1),uijk(1),lijk(2),uijk(2),l_ni,l_nj,uijk(3), &
                                .false., 'CUBIC', .true.)
            end if

         end if

      end if IF_READ

      IF_KOUNT0: if (F_kount > 0 .and. weight(max(1, F_kount)) > 0.) then

         !# Add increments to model and tracer states
         ni1 = l_ni ; nj1 = l_nj
         nullify(myptr0, myptr1, myptr2d)
         select case(iname0_S)
         case('tt')
            istat = gmm_get(gmmk_pw_tt_plus_s, myptr0)
!!$            print '(a,i4," : ",f10.5," : ",f10.5," < ",f10.5)','[iau_apply] ', F_kount, data0(l_ni/2, l_nj/2, l_nk/2), minval(data0), maxval(data0)
         case('uu')
            ni1 = l_ni ; nj1 = l_nj
            istat = gmm_get(gmmk_pw_uu_plus_s, myptr0)
            istat = gmm_get(gmmk_pw_vv_plus_s, myptr1)
         case('p0')
            istat = gmm_get(gmmk_st1_s, myptr2d)
         case default
            istat = clib_toupper(iname0_S)
            istat = gmm_get('TR/'//trim(iname0_S)//':P', myptr0)
         end select

         write(msg_S, '(a,i6,3(a,f12.6))') '; step=', F_kount, '; w=', &
              weight(F_kount), '; min=', minval(data0(1:ni1,1:l_nj,:)), &
              '; max=', maxval(data0(1:ni1,1:l_nj,:))
         call msg(MSG_INFO, '(iau_apply) Add increments (PE0): '// &
              trim(iname0_S)//trim(msg_S))
         if (associated(data1)) then
            write(msg_S, '(a,i6,3(a,f12.6))') '; step=', F_kount, '; w=', &
                 weight(F_kount), '; min=', minval(data1(1:ni1,1:l_nj,:)), &
                 '; max=', maxval(data1(1:ni1,1:l_nj,:))
            call msg(MSG_INFO, '(iau_apply) Add increments (PE0): '// &
                 trim(iname1_S)//trim(msg_S))
         end if

         if (associated(myptr0)) myptr0(1:ni1,1:l_nj,:) = &
              myptr0(1:ni1,1:l_nj,:) + weight(F_kount) * data0(1:ni1,1:l_nj,:)
         if (associated(myptr1) .and. associated(data1)) &
              myptr1(1:l_ni,1:nj1,:) = &
              myptr1(1:l_ni,1:nj1,:) + weight(F_kount) * data1(1:l_ni,1:nj1,:)
         if (associated(myptr2d)) then
            if (iname0_S == 'p0') then
               myptr2d(1:l_ni,1:l_nj) = myptr2d(1:l_ni,1:l_nj) + &
                    log( 1 + weight(F_kount)*data0(1:l_ni,1:l_nj,1) / &
                    (Cstv_pref_8*exp(myptr2d(1:l_ni,1:l_nj))) )
            else
               myptr2d(1:l_ni,1:l_nj) = myptr2d(1:l_ni,1:l_nj) + &
                    weight(F_kount) * data0(1:l_ni,1:l_nj,1)
            end if
         end if
      end if IF_KOUNT0

   end do DO_IVAR

   if (F_kount > 0) then
      call msg(MSG_INFO, ' IAU_APPLY - APPLIED ANALYSIS INCREMENTS VALID AT '//&
           trim(datev_S))
   end if
   call gemtime_stop(50)
   !--------------------------------------------------------------------------
   return
end subroutine iau_apply2
