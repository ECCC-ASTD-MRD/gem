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
subroutine test_inputio8()
   use, intrinsic :: iso_fortran_env, only: REAL64
   use iso_c_binding
   use clib_itf_mod
   use rmn_gmm
   use env_utils, only: env_get
   use testutils
   use vGrid_Descriptors
   use ezgrid_mod
   use input_mod, only: input_new, input_nbvar, input_set_basedir, &
        input_set_filename, input_setgridid, input_isvarstep, input_meta, input_get
   use inputio_mod, only: inputio_new, inputio_nbvar, inputio_set_filename, &
        inputio_nbvar, inputio_isvarstep, inputio_meta, inputio_get, &
        INPUT_FILES_GEOP, INPUT_FILES_CLIM, INPUTIO_T
   use mu_jdate_mod
   use ptopo_utils
!!$   use statfld_dm_mod
   use vgrid_wb
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2011-04, 2017-09, 2020-10
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"
   integer,parameter :: NDIGITS = 4
   integer,parameter :: NSTEP0 = 6
   integer,parameter :: NHR0 = 6
   integer,parameter :: NTEST = 1
   integer,parameter :: NVAR = 34
   integer,parameter :: NSKIP = 2
   integer,parameter :: GNI = 300 !# 198 !# 80
   integer,parameter :: GNJ = 240 !# 100 !# 40
   integer,parameter :: HALO = 0 !# 5
   integer,parameter :: NVALS = 10
   
!!$   integer,parameter :: NLVLS = 4
!!$   real,parameter :: HYBLVL_P(NLVLS) = (/0.5492443,0.7299818,0.8791828,0.9950425/)
   integer,parameter :: NLVLS = 35
   real,parameter :: HYBLVL_P(NLVLS) = (/ &
        0.196513996, &
        0.213311002, 0.2307     , 0.249127999, 0.2685     , 0.289252996, &
        0.3110     , 0.333811015, 0.3570     , 0.380953997, 0.4060     , &
        0.4320     , 0.4592     , 0.487251997, 0.515807986, 0.5453     , &
        0.575564981, 0.6070     , 0.639265001, 0.671959996, 0.704868972, &
        0.735844016, 0.765922010, 0.792918026, 0.818956017, 0.844021022, &
        0.865904987, 0.886768997, 0.906602025, 0.924284995, 0.940909982, &
        0.956465006, 0.970943987, 0.983220994, 0.994401991 &
        /)
!!$   real,parameter :: HYBLVL_P(80) = (/ &
!!$        0.000100000,  1.600285e-04, 2.561664e-04, 4.070487e-04, 6.320755e-04, &
!!$        9.528077e-04, 1.385376e-03, 1.964803e-03, 2.714585e-03, 3.643780e-03, &
!!$        4.794698e-03, 6.179091e-03, 7.825953e-03, 9.725254e-03, 1.192403e-02, &
!!$        0.01442,  &
!!$        0.017209800, 0.020292001, 0.023637200, 0.027247399, 0.031040501, &
!!$        0.035017598, 0.039265200, 0.043699902, 0.048151501, 0.052791599, &
!!$        0.057448599, 0.062121999, 0.066724300, 0.071254298, 0.075623199, &
!!$        0.079916999, 0.083958700, 0.0879     , 0.092072703, 0.096499197, &
!!$        0.10156    , 0.1073     , 0.1138     , 0.1213     , 0.1300     , &
!!$        0.14000    , 0.1516     , 0.1652     , 0.1802     , 0.196513996, &
!!$        0.213311002, 0.2307     , 0.249127999, 0.2685     , 0.289252996, &
!!$        0.3110     , 0.333811015, 0.3570     , 0.380953997, 0.4060     , &
!!$        0.4320     , 0.4592     , 0.487251997, 0.515807986, 0.5453     , &
!!$        0.575564981, 0.6070     , 0.639265001, 0.671959996, 0.704868972, &
!!$        0.735844016, 0.765922010, 0.792918026, 0.818956017, 0.844021022, &
!!$        0.865904987, 0.886768997, 0.906602025, 0.924284995, 0.940909982, &
!!$        0.956465006, 0.970943987, 0.983220994, 0.994401991/)

   real,parameter :: HYBLVL_H(NLVLS) = (/ &
        6800., 6600., 6400., 6200., 6000., &
        5800., 5600., 5400., 5200., 5000., &
        4800., 4600., 4400., 4200., 4000., &
        3800., 3600., 3400., 3200., 3000., &
        2800., 2600., 2400., 2200., 2000., &
        1800., 1600., 1400., 1200., 1000., &
        800., 600., 400., 200., 100. &
       /)
   
   character(len=32),parameter :: VGDTAGPT = 'gemthlvlp'
   character(len=32),parameter :: VGDTAGPM = 'gemmolvlp'
   character(len=32),parameter :: VGDTAGHT = 'gemthlvlh'
   character(len=32),parameter :: VGDTAGHM = 'gemmolvlh'
   character(len=32),parameter :: DATEO_S = '20090427.000000'

   character(len=512) :: dir0_S, dir1_S, dfiles_S, input_type_S, step_S, &
        incfg_S, incfg2_S, anal_S, anal2_S, clim_S, clim2_S, &
        geop_S, geop2_S, inrep_S, inrep2_S, time_S, mem_S
   integer :: istat, myproc, ndomains, idomain, ngrids, igrid, nstep, istep
   integer :: mini, maxi, lni, li0, minj, maxj, lnj, lj0, idt, nhr
   integer :: glbgridid, lclgridid

   integer, save :: inputid = -1
   type(INPUTIO_T), save :: inputobj
   integer, save :: nbvar = 0
   
   integer, external :: get_max_rss
   integer :: memuse
   real(kind=REAL64), external :: omp_get_wtime
   real(kind=REAL64) :: time0, time1
   ! ---------------------------------------------------------------------
   ndomains = 1
   ngrids = 1
   myproc = testutils_initmpi(F_ngrids=ngrids)
   call ptopo_init_var(ndomains, idomain, ngrids, igrid)

   call testutils_set_name('test_inputio8')
   istat = fstopc('MSGLVL', 'SYSTEM', RMN_OPT_SET)
   ! call msg_set_p0only(0)

   nstep = NSTEP0
   istat = env_get('TEST_NSTEP', nstep, NSTEP0, 1, 21600)
   istat = env_get('TEST_NHR', nhr, NHR0)
   idt = nhr*3600 / nstep

   istat = env_get('TEST_SUBNAME', input_type_S, F_default='dist', &
        F_normalize_L=.true., F_okvalues=(/'gem48', 'bloc ', 'dist '/))
    if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_ERROR,'TEST_SUBNAME should be one of: gem48, bloc, dist')
      return
   endif
  
   istat = env_get('ATM_MODEL_DFILES', dfiles_S)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_ERROR,'ATM_MODEL_DFILES not defined')
      return
   endif

!!$   dir0_S = '__to-delete_test_inputio8'
   dir0_S = './'
   write(dir1_S,'(a,i4.4)') trim(dir0_S)//'/', myproc
   istat = clib_mkdir(trim(dir1_S))
   istat = clib_chdir(trim(dir1_S))

   anal_S = trim(dfiles_S)//'/2009042700_000'
   anal2_S = './ANALYSIS' !'./analysis'
   clim_S = trim(dfiles_S)//'/climato'
   clim2_S = './CLIMATO'
   geop_S = trim(dfiles_S)//'/geophy/Gem_geophy.fst'
   geop2_S = './GEOPHY'
   inrep_S = trim(dfiles_S)
   inrep2_S = './INREP'
   incfg_S = trim(dfiles_S)//'/test_inputio8_table'
   incfg2_S = './physics_input_table'

   istat = clib_symlink(trim(anal_S), trim(anal2_S))
   istat = min(clib_symlink(trim(clim_S), trim(clim2_S)),istat)
   istat = min(clib_symlink(trim(geop_S), trim(geop2_S)),istat)
   istat = min(clib_symlink(trim(inrep_S), trim(inrep2_S)),istat)
   istat = min(clib_symlink(trim(incfg_S), trim(incfg2_S)),istat)
   if (.not.RMN_IS_OK(istat)) then
      call msg(MSG_WARNING, 'problem creating symlink')
   endif

   istat = priv_init_grid(mini, maxi, lni, li0, minj, maxj, lnj, lj0, &
        glbgridid, lclgridid)
   
   if (RMN_IS_OK(istat)) then
      istat = priv_init_io(inputid, inputobj, nbvar, idt, &
           incfg2_S, '.', geop2_S, &
           mini, maxi, lni, li0, minj, maxj, lnj, lj0, &
           glbgridid, lclgridid)
   endif
   
   if (RMN_IS_OK(istat)) then
      do istep = 0,nstep
         write(step_S, '(i5.5)') istep
         call msg(MSG_INFO, '============================================ '//trim(step_S))
         time0 = omp_get_wtime()
!!$         istat = priv_input(istep, input_type_S, inputid, inputobj, lclgridid, step_S, VGDTAGPT, VGDTAGPM)
         istat = priv_input(istep, input_type_S, inputid, inputobj, lclgridid, step_S, VGDTAGHT, VGDTAGHM)
         time1 = omp_get_wtime()
         memuse = get_max_rss()
         write(time_S,'(1pe13.6)') time1 - time0
         write(mem_S,'(i7)') memuse
         write(RMN_STDOUT,'(a)') 'Step '//trim(step_S)//' stats: '// &
              trim(time_S)//' sec ; '//trim(mem_S)//' Kbytes/PE]'
      enddo
   endif

!!$   istat = clib_unlink(trim(anal2_S))
!!$   istat = clib_unlink(trim(clim2_S))
!!$   istat = clib_unlink(trim(geop2_S))
!!$   istat = clib_unlink(trim(inrep2_S))
!!$
!!$   istat = clib_chdir('./../..')
!!$   istat = clib_unlink(trim(dir1_S))
!!$   istat = clib_unlink(trim(dir0_S))

   call testutils_stats_print()
   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
   call rpn_comm_finalize(istat)

!!$   ----------------------
!!$
!!$   istat = ptopo_io_set(testutils_npeio)
!!$
!!$
!!$   filename_S='../test_inputio8_table'
!!$
!!$   istat = clib_getenv('TEST_SUBNAME', input_type_S)
!!$   if (input_type_S == 'gem48') then
!!$      call testutils_set_name('test_inputio8_gem48')
!!$      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_GEM48)
!!$   else if (input_type_S == 'bloc') then
!!$      call testutils_set_name('test_inputio8_bloc')
!!$      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_BLOC)
!!$   else
!!$      call testutils_set_name('test_inputio8_dist')
!!$      call test_input2_get_data_odd(trim(filename_S), 5, PTOPO_IO)
!!$   endif

   ! ---------------------------------------------------------------------
   return

contains

   !/@*
   function priv_init_grid(mini, maxi, lni, li0, minj, maxj, lnj, lj0, &
        glbgridid, lclgridid) result(F_istat)
      implicit none
      integer, intent(OUT) :: mini, maxi, lni, li0, minj, maxj, lnj, lj0
      integer, intent(OUT) :: glbgridid, lclgridid
      integer :: F_istat
      !*@/
      logical, parameter :: ALONGX = .true.
      logical, parameter :: FILL = .false. !.true.
      real, parameter :: XSPAN = 300.
      real, parameter :: YSPAN = 120.
      real, parameter :: Hyb_rcoef(4) = [ 1., 1., -1., -1. ]

      character(len=2) :: grtyp, grref
      integer :: lnimax, lnjmax
      integer :: istat, ig1, ig2, ig3, ig4, i,j
      integer,pointer :: ip1listt(:), ip1listm(:)
      real :: xlat12, xlon1, xlon2, dx, dy, ax(GNI), ay(GNJ), ax0, ay0
      real, pointer :: data2d(:,:)
      real(REAL64) :: ptop_8
      type(gmm_metadata) :: meta2d
      type(vgrid_descriptor) :: vgrid1, vgrid2
      ! ---------------------------------------------------------------------
      !# horizontal
      F_istat = rpn_comm_topo(GNI, mini, maxi, lni, lnimax, HALO, li0, &
           ALONGX, FILL)
      F_istat = min(F_istat, &
           rpn_comm_topo(GNJ, minj, maxj, lnj, lnjmax, HALO, lj0, &
           &             .not.ALONGX, FILL))
      call testutils_assert_ok(RMN_IS_OK(F_istat), '(init_grid) rpn_comm_topo')

!!$      glbgridid = ezqkdef(GNI,GNJ, 'G', 0,0,0,0,0)
      grtyp='Z'
      
      grref='E'
      xlat12=0.
      xlon1=180.
      xlon2=270.
      dx = XSPAN / GNI
      dy = YSPAN / GNJ
      ax0 = (360. - XSPAN) / 2.
      ay0 = -1.* (YSPAN/ 2.)
      do i=1,GNI
         ax(i) = ax0 + float(i-1)*dx
      enddo
      do j=1,GNJ
         ay(j) = ay0 + float(j-1)*dy
      enddo
      call cxgaig(grref, ig1, ig2, ig3, ig4, xlat12, xlon1, xlat12, xlon2)
      glbgridid = ezgdef_fmem(GNI, GNJ, grtyp, grref, ig1, ig2, ig3, ig4, ax, ay)
      call testutils_assert_ok(RMN_IS_OK(glbgridid), '(init_grid) glbgridid')
      F_istat = min(F_istat, glbgridid)

      lclgridid = RMN_ERR
      if (input_type_S /= 'dist') then
         lclgridid = ezgrid_sub(glbgridid, li0, lj0, li0+lni-1, lj0+lnj-1)
         call testutils_assert_ok(RMN_IS_OK(lclgridid), '(init_grid) lclgridid')
         F_istat = min(F_istat, lclgridid)
      endif

      call gmm_build_meta2D(meta2d, &
           1, lni, HALO, HALO, lni, &
           1, lnj, HALO, HALO, lnj, &
           0, GMM_NULL_FLAGS)

      !# vertical Press
      nullify(data2d)
      istat = gmm_create('P0', data2d, meta2d)
      data2d = 99900.
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) gmm_create p0')
      F_istat = min(F_istat, istat)
      
      ptop_8 = 9575.0d0
      istat = vgd_new(vgrid1, &
           kind     = VGRID_HYBS_KIND, &
           version  = VGRID_HYBS_VER, &
           hyb      = HYBLVL_P, &
           ptop_8   = ptop_8, &
           pref_8   = 100000.d0, &
           rcoef1   = 1., &
           rcoef2   = 1.)
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) vgd_new HYBS ')
      F_istat = min(F_istat, istat)
      
      nullify(ip1listt, ip1listm)
      istat = vgd_get(vgrid1, key='VIPT', value=ip1listt)
      istat = vgrid_wb_put(VGDTAGPT, vgrid1, ip1listt, 'P0')
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) vgrid_wb_put P T')
      F_istat = min(F_istat, istat)
      
      istat = vgd_get(vgrid1, key='VIPM', value=ip1listm)
      istat = vgrid_wb_put(VGDTAGPM, vgrid1, ip1listm, 'P0')
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) vgrid_wb_put P M')
      F_istat = min(F_istat, istat)
      
      !# vertical Height
      nullify(data2d)
      istat = gmm_create('ME', data2d, meta2d)
      data2d = 10.
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) gmm_create ME')
      F_istat = min(F_istat, istat)
      
      istat = vgd_new(vgrid2, &
           kind    = VGRID_GC_KIND, &
           version = VGRID_GC_VER, &
           hyb     = HYBLVL_H, &
           rcoef1  = Hyb_rcoef(1), &
           rcoef2  = Hyb_rcoef(2), &
           rcoef3  = Hyb_rcoef(3), &
           rcoef4  = Hyb_rcoef(4), &
           dhm     = 0., &
           dht     = 0., &
           hyb_flat = HYBLVL_H(1))
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) vgd_new GC')
      F_istat = min(F_istat, istat)
      
      nullify(ip1listt, ip1listm)
      istat = vgd_get(vgrid2, key='VIPT', value=ip1listt)
      istat = vgrid_wb_put(VGDTAGHT, vgrid2, ip1listt, 'ME')
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) vgrid_wb_put H T')
      F_istat = min(F_istat, istat)
      
      istat = vgd_get(vgrid2, key='VIPM', value=ip1listm)
      istat = vgrid_wb_put(VGDTAGHM, vgrid2, ip1listm, 'ME')
      call testutils_assert_ok(RMN_IS_OK(istat), '(init_grid) vgrid_wb_put H M')
      F_istat = min(F_istat, istat)
      ! ---------------------------------------------------------------------
      return
   end function priv_init_grid
   
   !/@*
   function priv_init_io(F_inputid, F_inputobj, F_nbvar, F_idt, &
        F_incfg_S, F_basedir_S, F_geoname_S, &
        mini, maxi, lni, li0, minj, maxj, lnj, lj0, &
        glbgridid, lclgridid) result(F_istat)
      implicit none
      integer, intent(inout) :: F_inputid
      type(INPUTIO_T), intent(inout) :: F_inputobj
      integer, intent(inout) :: F_nbvar
      integer, intent(in) :: F_idt
      integer, intent(in) :: mini, maxi, lni, li0, minj, maxj, lnj, lj0
      integer, intent(in) :: glbgridid, lclgridid
      character(len=*), intent(in) :: F_incfg_S, F_basedir_S, F_geoname_S
      integer :: F_istat
      !*@/
      logical, parameter :: IS_DIR = .true.

      logical, save :: is_init_L = .false.
      integer, save :: istatus = RMN_ERR

      integer :: istat, iotype, rpncomm_gridid
      integer(INT64) :: jdateo
      ! ---------------------------------------------------------------------
      F_istat = istatus
      if (input_type_S == 'gem48') ptopo_iotype = PTOPO_BLOC
      if (is_init_L) return
      is_init_L = .true.

      jdateo = jdate_from_print(DATEO_S)

      rpncomm_gridid = -1
      F_istat = RMN_OK
      IF_GEM48: if (input_type_S == 'gem48') then
         F_inputid = input_new(jdateo, F_idt, F_incfg_S)
         istat = F_inputid
         if (RMN_IS_OK(istat)) then
            F_nbvar = input_nbvar(F_inputid)
            istat = input_set_basedir(F_inputid, F_basedir_S)
            istat = min(input_setgridid(F_inputid, lclgridid), istat)
         endif
         if (RMN_IS_OK(istat)) then
            istat = input_set_filename(F_inputid, 'geop', F_geoname_S, &
                 IS_DIR, INPUT_FILES_GEOP)
         endif
      else
         if (input_type_S == 'bloc') then
            iotype = PTOPO_BLOC
            istat = inputio_new(F_inputobj, jdateo, F_idt, F_incfg_S, &
                 F_basedir_S, lclgridid, lclgridid, rpncomm_gridid, &
                 F_li0=1, F_lin=lni, F_lj0=1, F_ljn=lnj, &
                 F_iotype=iotype)
         else
            iotype = PTOPO_IODIST
            istat = ptopo_io_set(testutils_npeio)
            rpncomm_gridid = rpn_comm_create_2dgrid(GNI, GNJ, mini, maxi, minj, maxj)
            istat = inputio_new(F_inputobj, jdateo, F_idt, F_incfg_S, &
                 F_basedir_S, glbgridid, glbgridid, rpncomm_gridid, &
!!$                 F_li0=li0, F_lin=li0+lni-1, F_lj0=lj0, F_ljn=lj0+lnj-1, &
                 F_li0=1, F_lin=lni, F_lj0=1, F_ljn=lnj, &
                 F_iotype=iotype)
         endif
         if (RMN_IS_OK(istat)) then
            F_nbvar = inputio_nbvar(F_inputobj)
            istat = inputio_set_filename(F_inputobj, 'geop', F_geoname_S, &
                 IS_DIR, INPUT_FILES_GEOP)
         endif
      endif IF_GEM48
      if (.not.RMN_IS_OK(istat)) &
           call msg(MSG_ERROR, '(phy_input) input module initialization problem')
      F_istat = min(istat, F_istat)

      call collect_error(F_istat)
      call testutils_assert_ok(RMN_IS_OK(F_istat), '(init_io) end')
      istatus = F_istat
      ! ---------------------------------------------------------------------
      return
   end function priv_init_io

   !/@
   function priv_input(F_step, input_type_S, inputid, inputobj, lclgridid, step_S, vgdtag_t_S, vgdtag_m_S) result(F_istat)
      implicit none
      integer, intent(IN) :: F_step, inputid, lclgridid
      character(len=*), intent(IN) :: input_type_S, step_S, vgdtag_t_S, vgdtag_m_S
      type(INPUTIO_T), intent(INOUT) :: inputobj
      integer :: F_istat
      !@/
      integer :: ivar, istat
      character(len=32) :: inname_S, inname2_S, vgrid_S
      real, pointer, dimension(:,:,:) :: data, data2
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      VARLOOP: do ivar=1,nbvar
         if (input_type_S == 'gem48') then
            istat = input_isvarstep(inputid, ivar, F_step)
         else
            istat = inputio_isvarstep(inputobj, ivar, F_step)
         endif
         if (.not.RMN_IS_OK(istat)) then
            cycle VARLOOP !var not requested at this step
         endif
         if (input_type_S == 'gem48') then
            istat = input_meta(inputid, ivar, inname_S, inname2_S)
         else
            istat = inputio_meta(inputobj%cfg, ivar, inname_S, inname2_S)
         endif
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(priv_input) problem getting input varname')
            cycle VARLOOP
         endif
         call msg(MSG_INFO, '-------------------------------------------- '//trim(inname_S))
         vgrid_S = vgdtag_t_S
         if (inname_S == 'uu') vgrid_S = vgdtag_m_S
         nullify(data, data2)
         if (input_type_S == 'gem48') then
            istat = input_get(inputid, ivar, F_step, lclgridid, vgrid_S, data, data2, F_ovname1_S=inname_S, F_ovname2_S=inname2_S)
         else
            istat = inputio_get(inputobj, ivar, F_step, data, data2, &
                 F_vgrid_S=vgrid_S, F_ovname1_S=inname_S, F_ovname2_S=inname2_S)
         endif
         call testutils_assert_ok(RMN_IS_OK(istat).and.associated(data), '(input) '//trim(step_S)//' get: '//trim(inname_S)//' '//(inname2_S))
         if (.not.RMN_IS_OK(istat)) then
            call msg(MSG_ERROR,'(priv_input) problem getting data')
            print *,'ERROR: RMN_IS_OK(istat),associated(data): ',RMN_IS_OK(istat),associated(data)
            if (associated(data)) deallocate(data, stat=istat)
            if (associated(data2)) deallocate(data2, stat=istat)
            cycle VARLOOP
         endif

         !#TODO: do something
         
         if (associated(data)) deallocate(data, stat=istat)
         if (associated(data2)) deallocate(data2, stat=istat)
         
      enddo VARLOOP
      F_istat = RMN_OK
      ! ---------------------------------------------------------------------
      return
   end function priv_input

end subroutine test_inputio8
