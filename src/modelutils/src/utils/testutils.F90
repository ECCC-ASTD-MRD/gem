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
!NOTE: testutils crash on AIX-powerpc7 with MAXLEN >= 1024
#define DEF_MSG_MAXLEN 512
#define DEF_MSG_OK 'OK   '
#define DEF_MSG_FAIL 'FAIL '

module testutils_stats
   use, intrinsic :: iso_fortran_env, only: REAL64
   public
   integer, save :: ntot = 0
   integer, save :: nok = 0

contains

   subroutine testutils_stats_print()
      implicit none
#include <rmnlib_basics.hf>
      logical :: canWrite_L
      integer :: msgLevelMin, msgUnit
      character(len=64) :: msgFormat_S
      ! ---------------------------------------------------------------------
      call msg_getInfo(canWrite_L, msgLevelMin, msgUnit, msgFormat_S)
      if (.not.canWrite_L) return
      write(RMN_STDERR, '(a)') '============================================='
      write(RMN_STDERR, '(a,i6,a,i6,a,i6,a)') 'Successfull assertions:',nok, ' /', ntot, ' (failed=', ntot-nok,')'
      write(RMN_STDERR, '(a)') '============================================='
      call flush(RMN_STDERR)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_stats_print

end module testutils_stats


module testutils
   use, intrinsic :: iso_fortran_env, only: INT64
   use iso_c_binding
   use ieee_arithmetic
   use rpn_comm_itf_mod
   use clib_itf_mod, only: clib_getenv, clib_isdigit, clib_tolower
   use wb_itf_mod
   use str_mod, only: str_tab2space
   use testutils_stats
   implicit none
   private
   public :: testutils_initmpi, testutils_getenv_int,testutils_verbosity, &
        testutils_set_name,testutils_set_tolerence,testutils_assert_eq, &
        testutils_assert_ok, testutils_stats_print, testutils_assert_neq, &
        testutils_assert_not_naninf

   integer, public :: testutils_myproc = 0
   integer, public :: testutils_npeio  = 1

#include <rmnlib_basics.hf>

   interface testutils_assert_ok
      module procedure testutils_assert_ok1
      module procedure testutils_assert_ok2
   end interface

   interface testutils_assert_not_naninf
      module procedure testutils_assert_not_naninf_r4
      module procedure testutils_assert_not_naninf_r4_1d
      module procedure testutils_assert_not_naninf_r4_2d
      module procedure testutils_assert_not_naninf_r4_3d
   end interface

   interface testutils_assert_neq
      module procedure testutils_assert_neq_i4
      module procedure testutils_assert_neq_r4
   end interface

   interface testutils_assert_eq
      module procedure testutils_assert_eq_s
      module procedure testutils_assert_eq_s_1d
      module procedure testutils_assert_eq_i4
      module procedure testutils_assert_eq_i4_1d
      module procedure testutils_assert_eq_r4
      module procedure testutils_assert_eq_r4_1d
      module procedure testutils_assert_eq_r8
      module procedure testutils_assert_eq_r8_1d
      module procedure testutils_assert_eq_L
      module procedure testutils_assert_eq_L_1d
   end interface

   character(len=5),parameter :: TESTOK_S   = DEF_MSG_OK
   character(len=5),parameter :: TESTFAIL_S = DEF_MSG_FAIL

   character(len=64),save :: m_name_S
   integer,save :: m_myproc = 0
   integer,save :: m_bloc_myproc = 0
   real,save :: m_r_tolerence_8 = 1.e-5

contains

   !/@*
   function testutils_initmpi(F_ngrids,F_npex,F_npey,F_nblocx,F_nblocy,F_npeio) result(F_myproc)
      implicit none
      !@objective init mpi for tests
      !@arguments
      integer,intent(in),optional :: F_ngrids,F_npex,F_npey,F_nblocx,F_nblocy,F_npeio
      !@return
      integer :: F_myproc
      !@author  Stephane Chamberland, 2010-05
      !*@/
      integer :: err,ngrids,igrid,myproc,numproc,npex,npey,mycol,myrow,nblocx,nblocy,mydomain,npeio
      external :: testutils_get_doms,testutils_get_npxy
      !---------------------------------------------------------------------
      call rpn_comm_mydomain(testutils_get_doms,mydomain)

      ngrids = 1
      if (present(F_ngrids)) ngrids = max(1,F_ngrids)
      npex = 0
      if (present(F_npex)) npex = F_npex
      npey = 0
      if (present(F_npey)) npey = F_npey
      igrid = rpn_comm_init_multigrid(testutils_get_npxy,myproc,numproc,npex,npey,ngrids)

      if (present(F_nblocx)) then
         nblocx = F_nblocx
      else
         err = testutils_getenv_int('MPI_NBLOCX',nblocx)
         if (.not.RMN_IS_OK(err)) nblocx = npex
      endif

      if (present(F_nblocy)) then
         nblocy = F_nblocy
      else
         err = testutils_getenv_int('MPI_NBLOCY',nblocy)
         if (.not.RMN_IS_OK(err)) nblocy = npey
      endif
      nblocx = min(nblocx,npex)
      nblocy = min(nblocy,npey)
!!$      err = rpn_comm_bloc(nblocx,nblocy)
      err = rpn_comm_mype(myproc,mycol,myrow)
!!$      call rpn_comm_rank(RPN_COMM_BLOC_COMM,m_bloc_myproc,err)
      testutils_myproc = myproc

      if (present(F_npeio)) then
         npeio = F_npeio
      else
         err = testutils_getenv_int('MPI_NPEIO',npeio)
!!$         print *,'MPI_NPEIO=', err, npeio
         if (.not.RMN_IS_OK(err)) npeio = 1
      endif
      testutils_npeio = npeio

      err = wb_put('ptopo/npx',npex)
      err = wb_put('ptopo/npy',npey)
      err = wb_put('ptopo/numproc',numproc)
      err = wb_put('ptopo/nblocx',nblocx)
      err = wb_put('ptopo/nblocy',nblocy)
      err = wb_put('ptopo/npeio',npeio)
      err = wb_put('ptopo/ngrids',ngrids)

      err = wb_put('ptopo/igrid',igrid)
      err = wb_put('ptopo/myproc',myproc)
      err = wb_put('ptopo/mycol',mycol)
      err = wb_put('ptopo/myrow',myrow)
!!$      err = wb_put('ptopo/mybloc',mybloc)
!!$      err = wb_put('ptopo/myblocx',myblocx)
!!$      err = wb_put('ptopo/myblocy',myblocy)

      call testutils_verbosity()

      m_myproc = myproc
      F_myproc = myproc
      !---------------------------------------------------------------------
      return
   end function testutils_initmpi


   !/@*
   function testutils_getenv_int(F_name_S,F_ival) result(F_istat)
      implicit none
      !@objective
      !@arguments
      character(len=*),intent(in) :: F_name_S
      integer,intent(out) :: F_ival
      !@returns
      integer :: F_istat
      !@author  Stephane Chamberland, 2010-05
      !*@/
      integer :: err
      character(len=256) :: tmp_S
      !---------------------------------------------------------------------
      F_istat = clib_getenv(trim(F_name_S),tmp_S)
      tmp_S = adjustl(tmp_S)
      F_istat = min(clib_isdigit(tmp_S(1:1)),F_istat)
      if (RMN_IS_OK(F_istat)) then
         read(tmp_S,fmt=*,iostat=err) F_ival
         if (err /= 0) F_istat = RMN_ERR
      endif
      !---------------------------------------------------------------------
      return
   end function testutils_getenv_int


   !/@*
   subroutine testutils_verbosity()
      implicit none
      !@objective Set verbosity
      !@author  Stephane Chamberland, 2011-09
      !*@/
      integer :: istat,myproc
      character(len=256) :: tmp_S
      !---------------------------------------------------------------------
      myproc = m_myproc
      istat = clib_getenv('TEST_VERBOSITY_PROC',tmp_S)
      if (RMN_IS_OK(istat)) then
         call str_tab2space(tmp_S)
         tmp_S = adjustl(tmp_S)
         istat = clib_tolower(tmp_S)
         select case(tmp_S(1:1))
         case('0') !- print from pe0 only
            myproc = m_myproc
         case('a') !- print from all proc
            myproc = 0
         case('b') !- print from bloc master only
            myproc = m_bloc_myproc
         end select
      endif
      call msg_set_p0only(myproc)

      call msg_set_minMessageLevel(MSG_CRITICAL)
      istat = wb_verbosity(WB_MSG_FATAL)
      istat = clib_getenv('TEST_VERBOSITY',tmp_S)
      if (RMN_IS_OK(istat)) then
         call str_tab2space(tmp_S)
         tmp_S = adjustl(tmp_S)
         istat = clib_tolower(tmp_S)
         select case(tmp_S(1:1))
         case('d')
            call msg_set_minMessageLevel(MSG_DEBUG)
            istat = wb_verbosity(WB_MSG_WARN) !WB_MSG_DEBUG)
         case('p')
            call msg_set_minMessageLevel(MSG_INFOPLUS)
            istat = wb_verbosity(WB_MSG_WARN) !WB_MSG_INFO)
         case('i')
            call msg_set_minMessageLevel(MSG_INFO)
            istat = wb_verbosity(WB_MSG_WARN) !WB_MSG_INFO)
         case('w')
            call msg_set_minMessageLevel(MSG_WARNING)
            istat = wb_verbosity(WB_MSG_WARN)
            istat= fstopc('MSGLVL','SYSTEM', RMN_OPT_SET)
        case('e')
            call msg_set_minMessageLevel(MSG_ERROR)
            istat = wb_verbosity(WB_MSG_ERROR)
            istat= fstopc('MSGLVL','SYSTEM', RMN_OPT_SET)
         end select
      else
         istat= fstopc('MSGLVL','SYSTEM', RMN_OPT_SET)
      endif
      !---------------------------------------------------------------------
      return
   end subroutine testutils_verbosity


   !/@
   subroutine testutils_set_name(F_name_S)
      implicit none
      character(len=*), intent(in) :: F_name_S
      !@/
      ! ---------------------------------------------------------------------
      m_name_S = F_name_S
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_set_name


   !/@
   subroutine testutils_set_tolerence(F_r_tolerence_8)
      implicit none
      real(REAL64), intent(in) :: F_r_tolerence_8
      !@/
      ! ---------------------------------------------------------------------
      m_r_tolerence_8 = F_r_tolerence_8
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_set_tolerence


   !/@
   subroutine testutils_assert_not_naninf_r4(F_userval,F_msg_S)
      implicit none
      real, intent(in) :: F_userval
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
!!$      real :: infinity_4
      ! ---------------------------------------------------------------------
      !#uses fortran 2003 ieee_is_nan() and ieee_is_finite()
!!$      infinity_4 = huge(infinity_4)
!!$      if (F_userval > infinity_4 .or. F_userval < -infinity_4) then
      if (.not.ieee_is_finite(F_userval)) then
         ok_L = .false.
         write(msg_S, '(a)') trim(F_msg_S)//' - Value is Inf'
      else
         ok_L = .true.
         msg_S = F_msg_S
      endif
      if (ok_L) then
!!$         if (F_userval /= F_userval) then
      if (ieee_is_nan(F_userval)) then
            ok_L = .false.
            write(msg_S, '(a)') trim(F_msg_S)//' - Value is NaN'
         else
            ok_L = .true.
            write(msg_S, '(a)') trim(F_msg_S)//' - Value is a valid number'
         endif
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_not_naninf_r4


   !/@
   subroutine testutils_assert_not_naninf_r4_1d(F_userval,F_msg_S)
      implicit none
      real, intent(in) :: F_userval(:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: i
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      ok_L = .true.
      msg_S = F_msg_S
      do i = lbound(F_userval,1), ubound(F_userval,1)
         if (.not.ieee_is_finite(F_userval(i))) then
            ok_L = .false.
            write(msg_S, '(a,1x,i0)') trim(F_msg_S)//' - Inf value found at: ',i
         endif
         if (ieee_is_nan(F_userval(i))) then
            ok_L = .false.
            write(msg_S, '(a,1x,i0)') trim(F_msg_S)//' - NaN value found at: ',i
         endif
         if (.not.ok_L) exit
      enddo
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_not_naninf_r4_1d


   !/@
   subroutine testutils_assert_not_naninf_r4_2d(F_userval,F_msg_S)
      implicit none
      real, intent(in) :: F_userval(:,:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: i, j
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      ok_L = .true.
      msg_S = F_msg_S
      do j = lbound(F_userval,2), ubound(F_userval,2)
         do i = lbound(F_userval,1), ubound(F_userval,1)
            if (.not.ieee_is_finite(F_userval(i,j))) then
               ok_L = .false.
               write(msg_S, '(a,1x,i0,1x,i0)') trim(F_msg_S)//' - Inf value found at: ',i,j
            endif
            if (ieee_is_nan(F_userval(i,j))) then
               ok_L = .false.
               write(msg_S, '(a,1x,i0,1x,i0)') trim(F_msg_S)//' - NaN value found at: ',i,j
            endif
            if (.not.ok_L) exit
         enddo
         if (.not.ok_L) exit
      enddo
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_not_naninf_r4_2d


   !/@
   subroutine testutils_assert_not_naninf_r4_3d(F_userval,F_msg_S)
      implicit none
      real, intent(in) :: F_userval(:,:,:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: i, j, k
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      ok_L = .true.
      msg_S = F_msg_S
      do k = lbound(F_userval,3), ubound(F_userval,3)
         do j = lbound(F_userval,2), ubound(F_userval,2)
            do i = lbound(F_userval,1), ubound(F_userval,1)
               if (.not.ieee_is_finite(F_userval(i,j,k))) then
                  ok_L = .false.
                  write(msg_S, '(a,1x,i0,1x,i0,1x,i0)') trim(F_msg_S)//' - Inf value found at: ',i,j,k
               endif
               if (ieee_is_nan(F_userval(i,j,k))) then
                  ok_L = .false.
                  write(msg_S, '(a,1x,i0,1x,i0,1x,i0)') trim(F_msg_S)//' - NaN value found at: ',i,j,k
               endif
               if (.not.ok_L) exit
            enddo
            if (.not.ok_L) exit
         enddo
      enddo
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_not_naninf_r4_3d


   !/@
   subroutine testutils_assert_eq_s(F_userval,F_expected,F_msg_S)
      implicit none
      character(len=*), intent(in) :: F_userval,F_expected
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S,v0,vsep,vu,v1
      logical :: ok_L
      integer :: mymaxlen
      ! ---------------------------------------------------------------------
      if (F_userval == F_expected) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         vsep = ' != '
         v0 = trim(F_msg_S)//' - got,exp:'
         mymaxlen = len(msg_S)/2 - len_trim(v0) - len_trim(vsep)
         vu = F_userval(1:min(len_trim(F_userval),mymaxlen))
         v1 = F_expected(1:min(len_trim(F_expected),mymaxlen))
         write(msg_S, '(4a)') trim(v0),trim(vu),trim(vsep),trim(v1)
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_s


   !/@
   subroutine testutils_assert_eq_s_1d(F_userval,F_expected,F_msg_S)
      implicit none
      character(len=*), intent(in) :: F_userval(:),F_expected(:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: n
      character(len=DEF_MSG_MAXLEN) :: msg_S,v0,vsep
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (all(F_userval == F_expected)) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         n = size(F_userval)
         if (n /= size(F_expected)) then
            write(msg_S, '(a,1x,i0,1x,a,1x,i0)') trim(F_msg_S)//' - nb values differ, got,exp:',size(F_userval),' != ',size(F_expected)
         else
            v0 = trim(F_msg_S)//' - got,exp:'
            vsep = ' != '
            !TODO-later: make sure what is written in msg_S is not too long
            if (n == 1) then
               write(msg_S, '(4a)') trim(v0),trim(F_userval(1)),trim(vsep),trim(F_expected(1))

            else if (n == 1) then
               write(msg_S,*) trim(v0),trim(F_userval(1)),', ', &
                    trim(F_userval(n)),trim(vsep),trim(F_expected(1)), &
                    ', ',trim(F_expected(n))
            else
               write(msg_S,*) trim(v0),trim(F_userval(1)),'...', &
                    trim(F_userval(n)),trim(vsep),trim(F_expected(1)), &
                    ', ...,',trim(F_expected(n))
            endif
         endif
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_s_1d


   !/@
   subroutine testutils_assert_neq_i4(F_userval,F_expected,F_msg_S)
      implicit none
      integer, intent(in) :: F_userval,F_expected
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (F_userval /= F_expected) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval,' == ',F_expected
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_neq_i4


   !/@
   subroutine testutils_assert_eq_i4(F_userval,F_expected,F_msg_S)
      implicit none
      integer, intent(in) :: F_userval,F_expected
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (F_userval == F_expected) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval,' != ',F_expected
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_i4


   !/@
   subroutine testutils_assert_eq_i4_1d(F_userval,F_expected,F_msg_S)
      implicit none
      integer, intent(in) :: F_userval(:),F_expected(:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: n
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (all(F_userval == F_expected)) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         n = size(F_userval)
         if (n /= size(F_expected)) then
            write(msg_S,*) trim(F_msg_S)//' - nb values differ, got,exp:',size(F_userval),' != ',size(F_expected)
         else
            if (n <= 3) then
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1:n),' != ',F_expected(1:n)
            else
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1),'...',F_userval(n),' != ',F_expected(1),'...',F_expected(n)
            endif
         endif
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_i4_1d


   !/@
   subroutine testutils_assert_neq_r4(F_userval,F_expected,F_msg_S)
      implicit none
      real, intent(in) :: F_userval,F_expected
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      !TODO-later: shoud we check relative error?
      if (abs(dble(F_userval) - dble(F_expected)) > m_r_tolerence_8) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval,' == ',F_expected
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_neq_r4

   !/@
   subroutine testutils_assert_eq_r4(F_userval,F_expected,F_msg_S)
      implicit none
      real, intent(in) :: F_userval,F_expected
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      !TODO-later: shoud we check relative error?
      if (abs(dble(F_userval) - dble(F_expected)) <= m_r_tolerence_8) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval,' != ',F_expected
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_r4


   !/@
   subroutine testutils_assert_eq_r4_1d(F_userval,F_expected,F_msg_S)
      implicit none
      real, intent(in) :: F_userval(:),F_expected(:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: n
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (all(abs(dble(F_userval) - dble(F_expected)) <= m_r_tolerence_8)) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         n = size(F_userval)
         if (n /= size(F_expected)) then
            write(msg_S,*) trim(F_msg_S)//' - nb values differ, got,exp:',size(F_userval),' != ',size(F_expected)
         else
            if (n <= 3) then
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1:n),' != ',F_expected(1:n)
            else
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1),'...',F_userval(n),' != ',F_expected(1),'...',F_expected(n)
            endif
         endif
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_r4_1d


   !/@
   subroutine testutils_assert_eq_r8(F_userval,F_expected,F_msg_S)
      implicit none
      real(REAL64), intent(in) :: F_userval,F_expected
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (abs(F_userval - F_expected) <= m_r_tolerence_8) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval,' != ',F_expected
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_r8


   !/@
   subroutine testutils_assert_eq_r8_1d(F_userval,F_expected,F_msg_S)
      implicit none
      real(REAL64), intent(in) :: F_userval(:),F_expected(:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: n
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (all(abs(F_userval - F_expected) <= m_r_tolerence_8)) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         n = size(F_userval)
         if (n /= size(F_expected)) then
            write(msg_S,*) trim(F_msg_S)//' - nb values differ, got,exp:',size(F_userval),' != ',size(F_expected)
         else
            if (n <= 3) then
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1:n),' != ',F_expected(1:n)
            else
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1),'...',F_userval(n),' != ',F_expected(1),'...',F_expected(n)
            endif
         endif
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_r8_1d


   !/@
   subroutine testutils_assert_eq_L(F_userval_L,F_expected_L,F_msg_S)
      implicit none
      logical, intent(in) :: F_userval_L,F_expected_L
      character(len=*), intent(in) :: F_msg_S
      !@/
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (F_userval_L .eqv. F_expected_L) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval_L,F_expected_L
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_l


   !/@
   subroutine testutils_assert_eq_l_1d(F_userval,F_expected,F_msg_S)
      implicit none
      logical, intent(in) :: F_userval(:),F_expected(:)
      character(len=*), intent(in) :: F_msg_S
      !@/
      integer :: n
      character(len=DEF_MSG_MAXLEN) :: msg_S
      logical :: ok_L
      ! ---------------------------------------------------------------------
      if (all(F_userval .eqv. F_expected)) then
         ok_L = .true.
         msg_S = F_msg_S
      else
         ok_L = .false.
         n = size(F_userval)
         if (n /= size(F_expected)) then
            write(msg_S,*) trim(F_msg_S)//' - nb values differ, got,exp:',size(F_userval),' != ',size(F_expected)
         else
            if (n <= 3) then
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1:n),' != ',F_expected(1:n)
            else
               write(msg_S,*) trim(F_msg_S)//' - got,exp:',F_userval(1),'...',F_userval(n),' != ',F_expected(1),'...',F_expected(n)
            endif
         endif
      endif
      call testutils_assert_ok0(ok_L,m_name_S,msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_eq_l_1d

   !/@
   subroutine testutils_assert_ok1(F_ok_L,F_msg_S)
      implicit none
      logical, intent(in) :: F_ok_L
      character(len=*), intent(in) :: F_msg_S
      !@/
      ! ---------------------------------------------------------------------
      call testutils_assert_ok0(F_ok_L,m_name_S,F_msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_ok1

   !/@
   subroutine testutils_assert_ok2(F_ok_L,F_name_S,F_msg_S)
      implicit none
      logical, intent(in) :: F_ok_L
      character(len=*), intent(in) :: F_name_S,F_msg_S
      !@/
      ! ---------------------------------------------------------------------
      call testutils_assert_ok0(F_ok_L,F_name_S,F_msg_S)
      ! ---------------------------------------------------------------------
      return
   end subroutine testutils_assert_ok2

   !TODO-later: assert gt, lt, ne for all type
end module testutils


!/@
subroutine testutils_assert_ok0(F_ok_L,F_name_S,F_msg_S)
   use testutils_stats
   implicit none
   logical, intent(in) :: F_ok_L
   character(len=*), intent(in) :: F_name_S,F_msg_S
!@/
#include <rmnlib_basics.hf>
   character(len=5),parameter :: TESTOK_S   = DEF_MSG_OK
   character(len=5),parameter :: TESTFAIL_S = DEF_MSG_FAIL
   integer :: istat,msgLevelMin,msgUnit
   logical :: canWrite_L
   character(len=64) :: msgFormat_S
   character(len=DEF_MSG_MAXLEN) :: msg2_S
   ! ---------------------------------------------------------------------
   istat = RMN_ERR
   if (F_ok_L) istat = RMN_OK
   call msg_getInfo(canWrite_L,msgLevelMin,msgUnit,msgFormat_S)
   call collect_error(istat)
   ntot = ntot + 1
   if (RMN_IS_OK(istat)) nok = nok + 1
   if (.not.canWrite_L) return
   msg2_S = adjustl(F_msg_S)
   if (RMN_IS_OK(istat)) then
      write(RMN_STDERR,'(a)') TESTOK_S//trim(F_name_S)//' - '//trim(msg2_S)
   else
      write(RMN_STDERR,'(a)') TESTFAIL_S//trim(F_name_S)//' - '//trim(msg2_S)
   endif
   call flush(RMN_STDERR)
   ! ---------------------------------------------------------------------
   return
end subroutine testutils_assert_ok0

