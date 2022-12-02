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
subroutine test_vinterp()
   use, intrinsic :: iso_fortran_env, only: REAL64
   use iso_c_binding
   use testutils
   use rmn_gmm
   use vGrid_Descriptors
   use vgrid_wb
   use vinterp_mod
   implicit none
   !@objective 
   !@author Stephane Chamberland, 2021-04
!@/
#include <rmnlib_basics.hf>
   include "rpn_comm.inc"

   integer, parameter :: GNI = 1
   integer, parameter :: GNJ = 1
   integer, parameter :: GNK_I = 4
   integer, parameter :: GNK_O = 5
   integer, parameter :: HALO = 0

   logical, parameter :: IS_PRESS = .true.
   logical, parameter :: IS_METER= .not.IS_PRESS

   real, parameter :: MB2PA = 100.

   real, parameter :: HY_REFVAL_I = 99900.
   real, parameter :: HY_REFVAL_O = 99905.
   real, parameter :: HY_HYB_I(GNK_I) = [ 0.5,  0.7, 0.8, 0.99 ]
   real, parameter :: HY_HYB_O(GNK_O) = [ 0.51, 0.6, 0.75, 0.9, 0.98 ]
   character(len=*), parameter :: HY_VNAME_I = 'LVL_HY_I'
   character(len=*), parameter :: HY_VNAME_O = 'LVL_HY_O'
   character(len=*), parameter :: HY_REFNAME_I = 'P0'
   character(len=*), parameter :: HY_REFNAME_O = 'P1'
   character(len=*), parameter :: HY_ALTNAME_I = ' '
   character(len=*), parameter :: HY_ALTNAME_O = ' '

   real, parameter :: HY2_REFVAL_I = HY_REFVAL_I
   real, parameter :: HY2_REFVAL_O = HY_REFVAL_O
   real, parameter :: HY2_HYB_I(GNK_I) = HY_HYB_I
   real, parameter :: HY2_HYB_O(GNK_O) = HY_HYB_O
   character(len=*), parameter :: HY2_VNAME_I = 'LVL_HY2_I'
   character(len=*), parameter :: HY2_VNAME_O = 'LVL_HY2_O'
   character(len=*), parameter :: HY2_REFNAME_I = 'P20'
   character(len=*), parameter :: HY2_REFNAME_O = 'P21'
   character(len=*), parameter :: HY2_ALTNAME_I = 'GZ_I'
   character(len=*), parameter :: HY2_ALTNAME_O = 'GZ_O'
   
   real, parameter :: PR_REFVAL_I = 99900.
   real, parameter :: PR_REFVAL_O = 99905.
   real, parameter :: PR_HYB_I(GNK_I) = [ 500., 700., 800., 999. ]
   real, parameter :: PR_HYB_O(GNK_O) = [ 510., 600., 750., 900., 980. ]
   character(len=*), parameter :: PR_VNAME_I = 'LVL_PR_I'
   character(len=*), parameter :: PR_VNAME_O = 'LVL_PR_O'
   character(len=*), parameter :: PR_REFNAME_I = ' '
   character(len=*), parameter :: PR_REFNAME_O = ' '
   character(len=*), parameter :: PR_ALTNAME_I = ' '
   character(len=*), parameter :: PR_ALTNAME_O = ' '

   real, parameter :: GC_REFVAL_I = 10.
   real, parameter :: GC_REFVAL_O = 11.
   real, parameter :: GC_HYB_I(GNK_I) = [ 6010., 4000., 2000., 10.]
   real, parameter :: GC_HYB_O(GNK_O) = [ 6000., 5000., 3000., 1000., 15. ]
   character(len=*), parameter :: GC_VNAME_I = 'LVL_GC_I'
   character(len=*), parameter :: GC_VNAME_O = 'LVL_GC_O'
   character(len=*), parameter :: GC_REFNAME_I = 'M0'
   character(len=*), parameter :: GC_REFNAME_O = 'M1'
   character(len=*), parameter :: GC_ALTNAME_I = 'PR_I'
   character(len=*), parameter :: GC_ALTNAME_O = 'PR_O'

   real, parameter :: GC2_REFVAL_I = GC_REFVAL_I
   real, parameter :: GC2_REFVAL_O = GC_REFVAL_O
   real, parameter :: GC2_HYB_I(GNK_I) = GC_HYB_I
   real, parameter :: GC2_HYB_O(GNK_O) = GC_HYB_O
   character(len=*), parameter :: GC2_VNAME_I = 'LVL_GC2_I'
   character(len=*), parameter :: GC2_VNAME_O = 'LVL_GC2_O'
   character(len=*), parameter :: GC2_REFNAME_I = 'M20'
   character(len=*), parameter :: GC2_REFNAME_O = 'M21'
   character(len=*), parameter :: GC2_ALTNAME_I = ' '
   character(len=*), parameter :: GC2_ALTNAME_O = ' '
   
   integer :: istat, hy_nk_i, hy_nk_o, gc_nk_i, gc_nk_o, pr_nk_i, pr_nk_o
!!$   integer :: myproc, ndomains, idomain, ngrids, igrid
   ! ---------------------------------------------------------------------
!!$   myproc = testutils_initmpi()
!!$   ndomains = 1
!!$   call ptopo_init_var(ndomains, idomain, ngrids, igrid)
!!$   istat = ptopo_io_set(testutils_npeio)

   !# export TEST_VERBOSITY_PROC=0
   !# export TEST_VERBOSITY=d
   call testutils_verbosity()

   call testutils_set_name('test_vinterp')
   
   hy_nk_i = def_vcoor(HY_HYB_I, HY_VNAME_I, HY_REFNAME_I, HY_REFVAL_I, IS_PRESS, HY_ALTNAME_I)
   call testutils_assert_ok(RMN_IS_OK(hy_nk_i), 'def_vcoor HY I')
   hy_nk_o = def_vcoor(HY_HYB_O, HY_VNAME_O, HY_REFNAME_O, HY_REFVAL_O, IS_PRESS, HY_ALTNAME_O)
   call testutils_assert_ok(RMN_IS_OK(hy_nk_o), 'def_vcoor HY O')

   pr_nk_i = def_vcoor(PR_HYB_I, PR_VNAME_I, PR_REFNAME_I, PR_REFVAL_I, IS_PRESS, PR_ALTNAME_I)
   call testutils_assert_ok(RMN_IS_OK(pr_nk_i), 'def_vcoor PR I')
   pr_nk_o = def_vcoor(PR_HYB_O, PR_VNAME_O, PR_REFNAME_O, PR_REFVAL_O, IS_PRESS, PR_ALTNAME_O)
   call testutils_assert_ok(RMN_IS_OK(pr_nk_o), 'def_vcoor PR O')

   gc_nk_i = def_vcoor(GC_HYB_I, GC_VNAME_I, GC_REFNAME_I, GC_REFVAL_I, IS_METER, GC_ALTNAME_I)
   call testutils_assert_ok(RMN_IS_OK(gc_nk_i), 'def_vcoor GC I')
   gc_nk_o = def_vcoor(GC_HYB_O, GC_VNAME_O, GC_REFNAME_O, GC_REFVAL_O, IS_METER, GC_ALTNAME_O)
   call testutils_assert_ok(RMN_IS_OK(gc_nk_o), 'def_vcoor GC O')

   hy_nk_i = def_vcoor(HY2_HYB_I, HY2_VNAME_I, HY2_REFNAME_I, HY2_REFVAL_I, IS_PRESS, HY2_ALTNAME_I)
   call testutils_assert_ok(RMN_IS_OK(hy_nk_i), 'def_vcoor HY2 I')
   hy_nk_o = def_vcoor(HY2_HYB_O, HY2_VNAME_O, HY2_REFNAME_O, HY2_REFVAL_O, IS_PRESS, HY2_ALTNAME_O)
   call testutils_assert_ok(RMN_IS_OK(hy_nk_o), 'def_vcoor HY2 O')
   
   gc_nk_i = def_vcoor(GC2_HYB_I, GC2_VNAME_I, GC2_REFNAME_I, GC2_REFVAL_I, IS_METER, GC2_ALTNAME_I)
   call testutils_assert_ok(RMN_IS_OK(gc_nk_i), 'def_vcoor GC2 I')
   gc_nk_o = def_vcoor(GC2_HYB_O, GC2_VNAME_O, GC2_REFNAME_O, GC2_REFVAL_O, IS_METER, GC2_ALTNAME_O)
   call testutils_assert_ok(RMN_IS_OK(gc_nk_o), 'def_vcoor GC2 O')

   
   istat = set_altcoor(GC_ALTNAME_I, PR_HYB_I, gc_nk_i, -100.)
   call testutils_assert_ok(RMN_IS_OK(istat), 'set_altcoor gc in')
   
   istat = set_altcoor(GC_ALTNAME_O, PR_HYB_O, gc_nk_o, -100.)
   call testutils_assert_ok(RMN_IS_OK(istat), 'set_altcoor gc out')

   istat = set_altcoor(HY2_ALTNAME_I, GC_HYB_I, hy_nk_i, 1000.)
   call testutils_assert_ok(RMN_IS_OK(istat), 'set_altcoor hy in')
   
   istat = set_altcoor(HY2_ALTNAME_O, GC_HYB_O, hy_nk_o, 1000.)
   call testutils_assert_ok(RMN_IS_OK(istat), 'set_altcoor hy out')

   
   call testutils_set_name('test_vinterp_hy-hy')
   call test_vinterp_run(HY_VNAME_I,HY_VNAME_O,hy_nk_i,hy_nk_o,'hy-hy')
   
   call testutils_set_name('test_vinterp_pr-pr')
   call test_vinterp_run(PR_VNAME_I,PR_VNAME_O,pr_nk_i,pr_nk_o,'pr-pr')
   
   call testutils_set_name('test_vinterp_gc-gc')
   call test_vinterp_run(GC_VNAME_I,GC_VNAME_O,gc_nk_i,gc_nk_o,'gc-gc')

   
   call testutils_set_name('test_vinterp_hy-pr')
   call test_vinterp_run(HY_VNAME_I,PR_VNAME_O,hy_nk_i,pr_nk_o,'hy-pr')
   
   call testutils_set_name('test_vinterp_pr-hy')
   call test_vinterp_run(PR_VNAME_I,HY_VNAME_O,pr_nk_i,hy_nk_o,'pr-hy')

   
   call testutils_set_name('test_vinterp_hy-gc')
   call test_vinterp_run(HY_VNAME_I,GC_VNAME_O,hy_nk_i,gc_nk_o,'hy-gc')
   
   call testutils_set_name('test_vinterhy_gc-hy')
   call test_vinterp_run(GC_VNAME_I,HY_VNAME_O,gc_nk_i,hy_nk_o,'gc-hy')

   
   call testutils_set_name('test_vinterp_pr-gc')
   call test_vinterp_run(PR_VNAME_I,GC_VNAME_O,pr_nk_i,gc_nk_o,'pr-gc')
   
   call testutils_set_name('test_vinterp_gc-pr')
   call test_vinterp_run(GC_VNAME_I,PR_VNAME_O,gc_nk_i,pr_nk_o,'gc-pr')

   call testutils_set_name('test_vinterp_hy2-gc2')
   call test_vinterp_run(HY2_VNAME_I,GC2_VNAME_O,hy_nk_i,gc_nk_o,'hy2-gc2')
   
   call testutils_set_name('test_vinterhy_gc2-hy2')
   call test_vinterp_run(GC2_VNAME_I,HY2_VNAME_O,gc_nk_i,hy_nk_o,'gc2-hy2')

   call testutils_stats_print()
!!$   call rpn_comm_barrier(RPN_COMM_WORLD, istat)
!!$   call rpn_comm_finalize(istat)
   ! ---------------------------------------------------------------------
   return

contains
   
   !/@
   function def_vcoor(F_hyb, F_name_S, F_ref_S, F_ref_val, F_press_L, F_alt_S) &
        result(F_nk)
      implicit none
      real, intent(in) :: F_hyb(:)
      character(len=*), intent(in) :: F_name_S, F_ref_S, F_alt_S
      real, intent(in) :: F_ref_val
      logical, intent(in) :: F_press_L
      integer :: F_nk
      !@/
      real(REAL64), parameter :: PTOP_8 = 9575.0d0
      real(REAL64), parameter :: PREF_8 = 100000.d0
      real, parameter :: Hyb_rcoef(4) = [ 1., 1., -1., -1. ]

      type(vgrid_descriptor) :: vgrid
      integer, pointer :: ip1listt(:)
      type(gmm_metadata) :: p0meta
      real, pointer :: p0data(:,:)
      integer :: istat
      ! ---------------------------------------------------------------------
      F_nk = RMN_ERR

      if (F_press_L) then
         if (F_ref_S == ' ') then
            istat = vgd_new(vgrid, &
                 kind    = VGRID_PRES_KIND , &
                 version = VGRID_PRES_VER, &
                 hyb     = F_hyb)
         else
            istat = vgd_new(vgrid, &
                 kind    = VGRID_HYBS_KIND, &
                 version = VGRID_HYBS_VER, &
                 hyb     = F_hyb, &
                 ptop_8  = PTOP_8, &
                 pref_8  = PREF_8, &
                 rcoef1  = Hyb_rcoef(1), &
                 rcoef2  = Hyb_rcoef(2))
         endif
      else
         istat = vgd_new(vgrid, &
              kind    = VGRID_GC_KIND, &
              version = VGRID_GC_VER, &
              hyb     = F_hyb, &
              rcoef1  = Hyb_rcoef(1), &
              rcoef2  = Hyb_rcoef(2), &
              rcoef3  = Hyb_rcoef(3), &
              rcoef4  = Hyb_rcoef(4), &
              dhm     = 0., &
              dht     = 0., &
              hyb_flat = F_hyb(1))
      endif
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(def_vcoor) Problem in vgd_new for: '//trim(F_name_S))
         return
      endif

      nullify(ip1listt)
      istat = vgd_get(vgrid,key='VIPT', value=ip1listt)
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(def_vcoor) Problem in vgd_get for: '//trim(F_name_S))
         return
      endif

      if (F_alt_S /= ' ') then
         istat = vgrid_wb_put(trim(F_name_S), vgrid, ip1listt, F_ref_S, F_altfld_S=F_alt_S)
      else
         istat = vgrid_wb_put(trim(F_name_S), vgrid, ip1listt, F_ref_S)
      endif
      if (.not.RMN_IS_OK(istat)) then
         call msg(MSG_ERROR,'(def_vcoor) Problem in vgrid_wb_put for: '//trim(F_name_S))
         return
      endif

      nullify(p0data)
      if (F_ref_S /= ' ') then
         call gmm_build_meta2D(p0meta, &
              1,GNI,HALO,HALO,GNJ, &
              1,GNI,HALO,HALO,GNJ, &
              0,GMM_NULL_FLAGS)
         istat = gmm_create(F_ref_S, p0data, p0meta)
         if (.not.(RMN_IS_OK(istat) .and. associated(p0data))) then
            call msg(MSG_ERROR,'(def_vcoor) Problem in gmm_create for: '//trim(F_name_S))
            return
         endif
         p0data = F_ref_val
      endif

      deallocate(ip1listt, stat=istat)
      F_nk = size(ip1listt)
      ! ---------------------------------------------------------------------
      return
   end function def_vcoor

   
   !/@
   function set_altcoor(F_alt_S, F_hyb_pp, F_nk, F_offset) &
        result(F_istat)
      implicit none
      character(len=*), intent(in) :: F_alt_S
      real, intent(in) :: F_hyb_pp(:)
      integer, intent(in) :: F_nk
      real, intent(in) :: F_offset
      integer :: F_istat
      !@/
      type(gmm_metadata) :: meta3d
      real, pointer :: data3d(:,:,:)
      integer :: istat, k
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR

      if (F_nk /= size(F_hyb_pp) + 2) then
         call msg(MSG_ERROR,'(set_altcoor) wrong size for: '//trim(F_alt_S))
         return
      endif

      call gmm_build_meta3D(meta3d, &
           1,GNI,HALO,HALO,GNI, &
           1,GNJ,HALO,HALO,GNJ, &
           1,F_nk,HALO,HALO,F_nk, &
           0,GMM_NULL_FLAGS)
      nullify(data3d)
      istat = gmm_create(F_alt_S, data3d, meta3d)
      if (.not.(RMN_IS_OK(istat) .and. associated(data3d))) then
         call msg(MSG_ERROR,'(def_vcoor) Problem in gmm_create for: '//trim(F_alt_S))
         return
      endif

      data3d(:,:,1) = F_hyb_pp(1) + (F_offset * 2.)
      data3d(:,:,2) = F_hyb_pp(2) + F_offset
      do k=1,size(F_hyb_pp)
         data3d(:,:,2+k) = F_hyb_pp(k)
      enddo
      data3d = data3d * MB2PA

      F_istat = RMN_OK
      ! ---------------------------------------------------------------------
      return
   end function set_altcoor

   
   !/@
   subroutine test_vinterp_run(F_vname_i, F_vname_o, F_nk_i, F_nk_o, F_msg_S)
      implicit none
      character(len=*), intent(in) :: F_vname_i, F_vname_o, F_msg_S
      integer, intent(in) :: F_nk_i, F_nk_o
      !@/
      integer :: istat, k
      real, pointer :: data_i(:,:,:), data_o(:,:,:)
      ! ---------------------------------------------------------------------
      allocate(data_i(GNI,GNJ,F_nk_i), data_o(GNI,GNJ,F_nk_o))
      
      do k=1,F_nk_i
         data_i(:,:,k) = float(k)
         ! print *,'i',k,data_i(1,1,k)
      enddo
      
      data_o = -1.
      istat = vinterp(data_o, F_vname_o, data_i, F_vname_i, F_msg_S=F_msg_S)
      call testutils_assert_ok(RMN_IS_OK(istat), 'vinterp')
      call testutils_assert_ok(all(data_o >= 0.), 'vinterp vall')

      do k=1,F_nk_o
         print *,'o',k,data_o(1,1,k)
      enddo

      deallocate(data_i, data_o, stat=istat)
      ! ---------------------------------------------------------------------
      return
   end subroutine test_vinterp_run
   

end subroutine test_vinterp
