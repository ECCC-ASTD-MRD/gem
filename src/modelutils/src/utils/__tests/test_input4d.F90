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

#include <rmn/msg.h>

!/@
subroutine test_input4d()
   use, intrinsic :: iso_fortran_env, only: INT64, REAL64
   use iso_c_binding
   use vGrid_Descriptors
   use clib_itf_mod
   use ezgrid_mod
   use inputio_mod
   use input_mod
   use mu_jdate_mod
   use ptopo_utils
   use statfld_dm_mod
   use testutils
   use vgrid_wb
   use rmn_gmm
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04, 2017-09
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"

   integer,parameter :: NVAR = 34
   integer,parameter :: NSKIP = 2

   integer,parameter :: NDIGITS = 9

   integer,parameter :: NSTEP = 1
   integer,parameter :: DELTAT = 9600  !#21600

   character(len=512) :: dir_S,filename_S,dfiles_S, name_S, &
        anal_S,anal2_S,clim_S,clim2_S,geop_S,geop2_S,inrep_S,inrep2_S
   integer :: istat, myproc, ndomains, idomain, ngrids, igrid
   ! ---------------------------------------------------------------------
   myproc = testutils_initmpi()
   ndomains = 1
   call ptopo_init_var(ndomains, idomain, ngrids, igrid)
   istat = ptopo_io_set(testutils_npeio)

   call testutils_set_name('test_input4d')
   istat = fstopc('MSGLVL', 'SYSTEM', RMN_OPT_SET)
   call msg_set_p0only(0)

   istat = clib_getenv('ATM_MODEL_DFILES',dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_ERROR,'ATM_MODEL_DFILES not defined')
      return
   endif

   write(dir_S,'(a,i4.4)') '__to-delete_test_input4d-',myproc
   istat = clib_mkdir(trim(dir_S))
   istat = clib_chdir(trim(dir_S))

!!$   anal_S = trim(dfiles_S)//'/bcmk/2009042700_000'
   anal_S = trim(dfiles_S)//'/bcmk_toctoc/2009042700_000'
!!$!   anal_S = '/users/dor/armn/sch/storage_model/zeta/ANALYSIS-mid'
!!$   anal_S = '/users/dor/armn/sch/storage_model/zeta/ANALYSIS-mid2'
   anal2_S = './ANALYSIS' !'./analysis'
   clim_S = trim(dfiles_S)//'/bcmk/climato'
   clim2_S = './CLIMATO'
   geop_S = trim(dfiles_S)//'/bcmk/geophy/Gem_geophy.fst'
   geop2_S = './GEOPHY'
!!$   inrep_S = trim(dfiles_S)//'/bcmk'
   inrep_S = trim(dfiles_S)//'/bcmk_toctoc'
   inrep2_S = './INREP'
   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))
   istat = clib_symlink(trim(anal_S),trim(anal2_S))
   istat = min(clib_symlink(trim(clim_S),trim(clim2_S)),istat)
   istat = min(clib_symlink(trim(geop_S),trim(geop2_S)),istat)
   istat = min(clib_symlink(trim(inrep_S),trim(inrep2_S)),istat)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_WARNING,'problem creating symlink')
   endif

   filename_S='../test_input4d_table'

   istat = clib_getenv('TEST_SUBNAME', name_S)
!!$   if (name_S == 'disthalo') then
!!$      call testutils_set_name('test_input4d_disthalo')
!!$      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_IO)
!!$   else if (name_S == 'blochalo') then
!!$      call testutils_set_name('test_input4d_blochalo')
!!$      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_BLOC)
!!$   else 
   if (name_S == 'dist') then
      call testutils_set_name('test_input4d_dist')
      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_IO)
   else if (name_S == 'bloc') then
      call testutils_set_name('test_input4d_bloc')
      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_BLOC)
    else
      call testutils_set_name('test_input4d_old')
      call test_input2_get_data_old(trim(filename_S), 5, -99)
   endif

   istat = clib_unlink(trim(anal2_S))
   istat = clib_unlink(trim(clim2_S))
   istat = clib_unlink(trim(geop2_S))
   istat = clib_unlink(trim(inrep2_S))

   istat = clib_chdir('./..')
   istat = clib_unlink(trim(dir_S))

   call testutils_stats_print()
   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return

contains

   !/@
   subroutine test_input2_get_data_even(F_filename_S, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer, intent(in) :: iotype
      !@/
      integer,parameter :: GNI = 8
      integer,parameter :: GNJ = 4
      integer, parameter :: HALO = 0
      integer,parameter :: MAXNK = 26
      integer,parameter :: NVALS = 10
 
      integer :: allvals2(NVALS, NSTEP+1, NVAR, MAXNK)
      ! ---------------------------------------------------------------------
      allvals2 = 0
      call test_input4d_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input2_get_data_even


   !/@
   subroutine test_input2_get_data_odd(F_filename_S, halo, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer, intent(in) :: halo, iotype
      !@/
      integer,parameter :: GNI = 81
      integer,parameter :: GNJ = 41
      integer,parameter :: MAXNK = 26
      integer,parameter :: NVALS = 10

      integer :: allvals2(NVALS, NSTEP+1, NVAR, MAXNK)
       ! ---------------------------------------------------------------------
      allvals2 = 0
      call test_input4d_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input2_get_data_odd

   !/@
   subroutine test_input2_get_data_old(F_filename_S, halo, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer, intent(in) :: halo, iotype
      !@/
      integer,parameter :: GNI = 81
      integer,parameter :: GNJ = 41
      integer,parameter :: MAXNK = 26
      integer,parameter :: NVALS = 10

      integer :: allvals2(NVALS, NSTEP+1, NVAR, MAXNK)
       ! ---------------------------------------------------------------------
      allvals2 = 0
      call test_input4d_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input2_get_data_old


   !/@
   subroutine test_input4d_run(F_filename_S, GNI, GNJ, HALO, allvals2, iotype)
      implicit none
      character(len=*),intent(in) :: F_filename_S
      integer,intent(in) :: GNI, GNJ, HALO
      integer :: allvals2(:,:,:,:), iotype
      !@/
      integer,parameter :: NK = 2
      logical, parameter :: ALONGX = .true.
      logical, parameter :: FILL = .false. !.true.
      real,parameter :: hyb(4) = (/0.5492443,0.7299818,0.8791828,0.9950425/)
      character(len=512) :: dateo_S,varname_S,varname2_S,skip_list_S(4),step_S
      ! character(len=512) :: dummy_S
      character(len=8) :: result_S(0:NSTEP,NVAR)
      integer :: istat,dateo,istep,dt,ii,gridid,nbvar,ivar,nhyb,k,ikind,ip1list(NK)
      integer :: mini, maxi, lni, lnimax, li0, minj, maxj, lnj, lnjmax, lj0, rpncomm_gridid,lijk0(3),lijk1(3),uijk0(3),uijk1(3), inputid
      logical :: ok_L(0:NSTEP,NVAR), ok1_L, input_use_old_l
      real :: hyblist(NK)
      real(8) :: ptop_8
      real, pointer :: data(:,:,:),data2(:,:,:),p0data(:,:)
      integer,pointer :: ip1listt(:), ip1listm(:)
      type(vgrid_descriptor) :: vgrid1
      type(gmm_metadata) :: p0meta

      type(INPUTIO_T) :: inputobj
      integer(INT64) :: jdateo

      integer :: ijkmin(3),ijkmax(3),ilvl,vals2(10), lclgridid
      real(REAL64) :: mean,var,rmin,rmax, fact
      ! ---------------------------------------------------------------------
      result_S = ' '
      skip_list_S = ' '

      ok_L = .false.

      gridid = ezqkdef(GNI,GNJ, 'G', 0,0,0,0,0)


      istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
           ALONGX, FILL)
      istat = rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
           .not.ALONGX, FILL)
      input_use_old_l = (.not.(iotype == PTOPO_IO .or. iotype == PTOPO_BLOC))
      if (iotype == PTOPO_IO) then
         rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)
         call testutils_assert_ok(RMN_IS_OK(rpncomm_gridid), 'rpncomm_gridid')
      else  !# if (iotype == PTOPO_BLOC) then
         rpncomm_gridid = -1
         lclgridid = ezgrid_sub(gridid, li0, lj0, li0+lni-1, lj0+lnj-1)
         if (.not.RMN_IS_OK(lclgridid)) then
!!$            write(dummy_S, '(i7,a,i7,4i6,a)') lclgridid,"= ezgrid_sub(",gridid,li0, lj0, li0+lni-1, lj0+lnj-1,")"
!!$            call msg(MSG_ERROR,'test_input7 - cannot creat lclgridid: '//trim(dummy_S))
!!$            print *,"ERROR: ",trim(dummy_S)
!!$            return
            lclgridid = gridid
         endif
         gridid = lclgridid
      endif

      lijk0(1:3) = (/1, 1, 1/)
      uijk0(1:3) = (/lni, lnj, NK/)

      call gmm_build_meta2D(p0meta, &
           1,lni,HALO,HALO,lni, &
           1,lnj,HALO,HALO,lnj, &
           0,GMM_NULL_FLAGS)
      nullify(p0data)
      istat = gmm_create('P0',p0data,p0meta)
      p0data = 99900.

      ptop_8 = 9575.0d0
      istat = vgd_new(vgrid1, &
           kind     = VGRID_HYBS_KIND, &
           version  = VGRID_HYBS_VER, &
           hyb      = hyb, &
           ptop_8   = ptop_8, &
           pref_8   = 100000.d0, &
           rcoef1   = 1., &
           rcoef2   = 1.)
      nullify(ip1listt, ip1listm)
      istat = vgd_get(vgrid1,key='VIPT',value=ip1listt)
      istat = vgrid_wb_put('gemthlvl',vgrid1,ip1listt,'P0')
      istat = vgd_get(vgrid1,key='VIPM',value=ip1listm)
      istat = vgrid_wb_put('gemmolvl',vgrid1,ip1listm,'P0')

      dateo_S = '20090427.000000'
      jdateo = jdate_from_print(dateo_S)
      dateo = jdate_to_cmc(jdateo)
      nhyb = 2
      hyblist(1:2) = (/0.985,0.995/)
      do k=1,nhyb
         ikind = RMN_CONV_HY
         call convip_plus(ip1list(k),hyblist(k),ikind,RMN_CONV_P2IPNEW,'',.false.)
      enddo
      dt = DELTAT
      if (input_use_old_l) then
         inputid = input_new(dateo, dt, F_filename_S)
         istat = input_setgridid(inputid, lclgridid)
         istat = min(input_set_basedir(inputid, '.'), inputid)
      else
         istat = inputio_new(inputobj, jdateo, dt, F_filename_S, '.', &
              gridid, gridid, rpncomm_gridid, &
              F_li0=1, F_lin=lni, F_lj0=1, F_ljn=lnj, F_iotype=iotype)
      endif
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'test_input3 - cannot init inputio_mod for: '//trim(F_filename_S))
         return
      endif
      if (input_use_old_l) then
         nbvar = input_nbvar(inputid)
      else
         nbvar = inputio_nbvar(inputobj)
      endif
!!$      call testutils_assert_eq(nbvar, NVAR+NSKIP, 'inputio_nbvar')
         write(step_S, '(i6)') nbvar
         call msg(MSG_INFO, 'nbvar='//trim(step_S))

      STEPLOOP: do istep = 0, NSTEP
         write(step_S, '(i2.2)') istep
         call msg(MSG_INFO, '============================================ '//trim(step_S))
         ii = 0
         VARLOOP: do ivar=1,nbvar
            if (input_use_old_l) then
               istat = input_isvarstep(inputid, ivar, istep)
            else
               istat = inputio_isvarstep(inputobj, ivar, istep)
            endif
            if (.not.RMN_IS_OK(istat)) cycle VARLOOP

            if (input_use_old_l) then
               istat = input_meta(inputid, ivar, varname_S, varname2_S)
            else
               istat = inputio_meta(inputobj%cfg, ivar, varname_S, varname2_S)
            endif
            if (any(varname_S == skip_list_S)) cycle VARLOOP
            call msg(MSG_INFO, '-------------------------------------------- '//trim(varname_S)//' : '//trim(step_S))

            nullify(data,data2)

            if (input_use_old_l) then
               if (varname_S == 'uu') then
                  istat = input_get(inputid, ivar, istep, lclgridid, 'gemmolvl', data, data2)
               else
                  istat = input_get(inputid, ivar, istep, lclgridid, 'gemthlvl', data, data2)
               endif
            else
               if (varname_S == 'uu') then
                  istat = inputio_get(inputobj, ivar, istep, data, data2, &
                       F_vgrid_S='gemmolvl')
               else
                  istat = inputio_get(inputobj, ivar, istep, data, data2, &
                       F_vgrid_S='gemthlvl')
               endif
            endif
            if (.not.(RMN_IS_OK(istat) .and. associated(data))) then
               call msg(MSG_WARNING,'test_input3 - var not found: '//trim(varname_S))
               cycle VARLOOP !exit STEPLOOP
            endif
!!$            print *,'====',inputobj%fid%files(1:inputobj%fid%nfiles)%unit

            ii = ii + 1
!!$            call testutils_assert_eq(varname_S, result_S(istep, ii), 'var '//trim(varname_S)//':'//trim(step_S))

            lijk1 = lbound(data)
            uijk1 = ubound(data)
            call testutils_assert_ok( &
                 all(lijk1(1:2) == lijk0(1:2)) .and. &
                 all(uijk1(1:2) == uijk0(1:2)), 'bounds '//trim(varname_S)//':'//trim(step_S))
            if (.not.(all(lijk1(1:2) == lijk0(1:2)) .and. &
                 all(uijk1(1:2) == uijk0(1:2)))) then
               print *,'bounds  '//trim(varname_S)//':'//trim(step_S)//':',maxi,maxj
               print *,'bounds0 '//trim(varname_S)//':'//trim(step_S)//':',uijk0(1:2)
               print *,'bounds1 '//trim(varname_S)//':'//trim(step_S)//':',uijk1(1:2)
               call flush(6)
            endif
            if (varname_S == 'tt') then
               call testutils_assert_ok( &
                    lijk1(3) == 1 .and. uijk1(3)==size(ip1listt), 'bounds k '//trim(varname_S)//':'//trim(step_S))
            elseif (varname_S == 'uu') then
               call testutils_assert_ok( &
                    lijk1(3) == 1 .and. uijk1(3)==size(ip1listm), 'bounds k '//trim(varname_S)//':'//trim(step_S))
            endif

            call testutils_assert_not_naninf(data, 'validity '//trim(varname_S)//':'//trim(step_S))
            if (varname2_S/=' ') &
                 call testutils_assert_not_naninf(data2, 'validity '//trim(varname2_S)//':'//trim(step_S))

            do ilvl=1,uijk1(3)
               call statfld_dm(data(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
               ok1_L = .true.
               if (ptopo_grid_ipe == RPN_COMM_MASTER) then
                  fact = 10.**nint(float(NDIGITS)-max(-20.,log10(maxval((/abs(mean), abs(var), abs(rmin), abs(rmax)/)) + tiny(rmax))))
                  vals2 = (/nint(fact*mean), nint(fact*var), nint(fact*rmin), nint(fact*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
                  if (all(vals2(:) == allvals2(:,istep+1,ii,ilvl))) then
                     ok1_L = .true.
                  else
                     ok1_L = .false.
                     print &
                          '(a,i1,a,i2,a,i2,a,i11,",",i11,",",i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)', &
                          "allvals2(:,", istep+1, ",", ii, ",", ilvl, ") = (/", vals2, "/) !"//varname_S(1:2)
!!$                     print &
!!$                          '(a,i1,a,i2,a,i2,a,i11,",",i11,",",i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)', &
!!$                          "allvals0(:,", istep+1, ",", ii, ",", ilvl, ") = (/", allvals2(:,istep+1,ii,ilvl), "/) !"//varname_S(1:2)
                  endif
               endif
!!$               write(dummy_S,*) ilvl
!!$               call testutils_assert_ok(ok1_L, 'stats '//varname_S(1:2)//':'//trim(step_S)//' '//trim(dummy_S))
            enddo
            if (varname2_S/=' ') then
               ii = ii + 1
               do ilvl=1,uijk1(3)
                  call statfld_dm(data2(:,:,ilvl), mean, var, rmin, rmax, ijkmin, ijkmax)
                  ok1_L = .true.
                  if (ptopo_grid_ipe == RPN_COMM_MASTER) then
                     fact = 10.**nint(5.-max(-20.,log10(maxval((/abs(mean), abs(var), abs(rmin), abs(rmax)/)) + tiny(rmax))))
                     vals2 = (/nint(fact*mean), nint(fact*var), nint(fact*rmin), nint(fact*rmax), ijkmin(1), ijkmin(2), ijkmin(3), ijkmax(1), ijkmax(2), ijkmax(3)/)
                     if (all(vals2(:) == allvals2(:,istep+1,ii,ilvl))) then
                        ok1_L = .true.
                     else
                        ok1_L = .false.
                        print &
                             '(a,i1,a,i2,a,i2,a,i11,",",i11,",",i11,",",i11,",",i2,",",i2,",",i2,",",i2,",",i2,",",i2,a)', &
                             "allvals2(:,", istep+1, ",", ii, ",", ilvl, ") = (/", vals2, "/) !"//varname2_S(1:2)
                     endif
                  endif
!!$                  write(dummy_S,*) ilvl
!!$                  call testutils_assert_ok(ok1_L, 'stats '//varname2_S(1:2)//':'//trim(step_S)//' '//trim(dummy_S))
               enddo
            endif
 
            if (associated(data)) deallocate(data,stat=istat)
            if (associated(data2)) deallocate(data2,stat=istat)
         enddo VARLOOP
!!$         if (ii < NVAR) then
!!$            ii = ii + 1
!!$            if (result_S(istep, ii) == '') then
!!$               call testutils_assert_ok(.true., 'nvar'//trim(step_S))
!!$            else
!!$               call testutils_assert_eq(' ', result_S(istep, ii), 'missing "'//trim(result_S(istep, ii))//'":'//trim(step_S))
!!$            endif
!!$         else
!!$            call testutils_assert_ok(.true., 'nvar'//trim(step_S))
!!$         endif
      enddo STEPLOOP
      if (input_use_old_l) then
         istat = input_close_files()
      else
         istat = inputio_close_files(inputobj)
      endif
      ! ---------------------------------------------------------------------
      return
   end subroutine test_input4d_run


end subroutine test_input4d
