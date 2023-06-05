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
module fst_read_mod
   use, intrinsic :: iso_fortran_env, only: REAL64, INT64
   use vGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use vgrid_wb
   use vgrid_from_file_mod
   use cmcdate_mod
   use ezgrid_mod
   use hinterp4yy_mod
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2011-04
   !@description
   ! Public functions
   public :: fst_find, fst_read, fst_rdhint, fst_get_hgridid, fst_get_gridid, &
        fst_get_vgrid, fst_getmeta, fst_checkalloc, fst_checkalloc2
   public :: fst_find_vect, fst_find_0, fst_find_3d_vect, fst_find_3d_0
   ! Public constants
   integer, parameter, public :: FST_FIND_LT    = -2 !find nearest datev .lt. F_datev
   integer, parameter, public :: FST_FIND_LE    = -1 !find nearest datev .le. F_datev
   integer, parameter, public :: FST_FIND_NEAR  = 0  !find nearest datev
   integer, parameter, public :: FST_FIND_GE    = 1  !find nearest datev .ge. F_datev
   integer, parameter, public :: FST_FIND_GT    = 2  !find nearest datev .gt. F_datev
   integer, parameter, public :: FST_FIND_DIAG_T = -999  !find thermo diag lvl
   integer, parameter, public :: FST_FIND_DIAG_M = -998  !find mom diag lvl
!@/

#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

   interface fst_get_gridid !# Name kept for backward compatibility
      module procedure fst_get_hgridid
   end interface

   interface fst_find
      module procedure fst_find_vect
      module procedure fst_find_0
      module procedure fst_find_3d_vect
      module procedure fst_find_3d_0
   end interface fst_find

   interface fst_read
      module procedure fst_read_3d_r4_vect
      module procedure fst_read_3d_r4
   end interface

   interface fst_rdhint
      module procedure fst_rdhint_3d_r4_vect
      module procedure fst_rdhint_3d_r4
   end interface fst_rdhint

   integer,parameter :: MYMSG_QUIET = 99
   integer, parameter :: FST_MIN_TIME_RES = 40  !fstinf min time/datev (sec) resolution in search
   integer, parameter :: FST_MIN_TIME_RES_PRE80 = 1840  !fstinf min time/datev (sec) resolution in search
   real(REAL64), parameter :: SEC_PER_HR = 3600.d0

contains


   !/@
   function fst_find_0(F_fileid,F_nomvar,F_datev,F_ip1,F_ip2,F_ip3, &
        F_datevfuzz,F_fuzzopr,F_typvar_S) result(F_key)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      character(len=*),intent(in) :: F_nomvar
      integer,intent(inout) :: F_datev
      integer,intent(in) :: F_fileid,F_ip1,F_ip2,F_ip3
      integer,intent(in),optional :: F_datevfuzz,F_fuzzopr
      character(len=*),intent(inout),optional :: F_typvar_S
      !@return
      integer :: F_key
      !@author
      !@revision
      !@/
      integer,parameter :: NMAX = 9999
!!$      character(len=12) :: dummy_S
      character(len=1) :: grtyp_S
      character(len=2) :: typvar_S,typvar1_S
      character(len=4) :: nomvar_S
      character(len=12):: etiket_S
      character(len=256) :: msg_S
      integer :: i,istat,keylist(NMAX),nkeys,keys2(NMAX),nkeys2,&
           datevfuzz,fuzzopr,min_time_res
      integer(INT64) :: dist64, dist064
      integer :: ni1,nj1,nk1, datev, datev0, searchdate, &
           dateo,deet,npas,nbits, datyp, ip1, ip2, ip3, &
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
!!$      real :: zp1
      real(REAL64) :: nhours_8
      !---------------------------------------------------------------------
      nomvar_S = F_nomvar
      F_key = RMN_ERR
      if (F_fileid <= 0) then
         call msg(MSG_WARNING, '(fst) Find - Invalid file unit')
         return
      endif
      datevfuzz = 0
      fuzzopr = FST_FIND_NEAR
      typvar_S = RMN_ANY_TYP
      if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
      if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      write(msg_S,'(a,3i10,a,i10,i3,a)') &
           '(fst) Find, looking for: '//nomvar_S(1:4)// &
           ' [datev='//cmcdate_toprint(F_datev)// &
           '] [ip123=', F_ip1,F_ip2,F_ip3, &
           '] [typvar='//typvar_S(1:2)// &
           '] [datevfuzz=', datevfuzz,fuzzopr,']'
      call msg(MSG_DEBUG,msg_S)

      searchdate = F_datev
      min_time_res = FST_MIN_TIME_RES
      if (searchdate >= 0) then
         if (cmcdate_year(searchdate) <= 1980) &
              min_time_res = FST_MIN_TIME_RES_PRE80
         select case(fuzzopr)
         case(FST_FIND_GT)
            nhours_8 = dble(min_time_res)/SEC_PER_HR
            call incdatr(searchdate,F_datev,nhours_8)
            fuzzopr = FST_FIND_GE
         case(FST_FIND_LT)
            nhours_8 = dble(-min_time_res)/SEC_PER_HR
            call incdatr(searchdate,F_datev,nhours_8)
            fuzzopr = FST_FIND_LE
         end select
      endif

      ip1 = F_ip1
      if (ip1 == 1200) ip1 = 0
      if (ip1 <= 0) then
         F_key = fstinf(F_fileid,ni1,nj1,nk1,searchdate,RMN_ANY_ETK,  &
              ip1,F_ip2,F_ip3,typvar_S,F_nomvar)
         if (RMN_IS_OK(F_key)) then
            F_datev = searchdate
!!$            if (typvar_S /= RMN_ANY_TYP .and. present(F_typvar_S)) then
!!$               F_typvar_S =  !#TODO:
!!$            endif
            write(msg_S,'(a,3i10,a,i12,a)') '(fst) Found: '//nomvar_S(1:4)// &
                 ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=',F_ip1,F_ip2,F_ip3, &
                 '] [key=',F_key,']'
            call msg(MSG_DEBUG,msg_S)
            return
         endif
      endif

      ip1 = F_ip1
      if (ip1 == 0) ip1 = 1200
      if (F_ip1 >= 0) then
!!$         call convip_plus(ip1, zp1, kind, RMN_CONV_IP2P, dummy_S, .not.RMN_CONV_USEFORMAT_L)
!!$         ip1 = ip1_all(zp1,kind)
         ip1 = priv_ip1_all(ip1)
         F_key = fstinf(F_fileid,ni1,nj1,nk1,searchdate,RMN_ANY_ETK,  &
              ip1,F_ip2,F_ip3,typvar_S,F_nomvar)
!!$            F_key = fstinf(F_fileid,ni1,nj1,nk1,searchdate,RMN_ANY_ETK,  &
!!$                 ip1_all(zp1,kind),F_ip2,F_ip3,typvar_S,F_nomvar)
      endif

      if (RMN_IS_OK(F_key) .or. datevfuzz <= 0) then
         F_datev = searchdate
!!$            if (typvar_S /= RMN_ANY_TYP .and. present(F_typvar_S)) then
!!$               F_typvar_S =  !#TODO:
!!$            endif
         if (RMN_IS_OK(F_key)) then
            write(msg_S,'(a,3i10,a,i12,a)') '(fst) Found: '//nomvar_S(1:4)// &
                 ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=',F_ip1,F_ip2,F_ip3, &
                 '] [key=',F_key,']'
            call msg(MSG_DEBUG,msg_S)
         else
            write(msg_S,'(a,3i10,a)') '(fst) Not Found: '//nomvar_S(1:4)// &
                 ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=',F_ip1,F_ip2,F_ip3,']'
            call msg(MSG_DEBUG, msg_S)
         endif
         return
      endif

      nkeys = 0
      if (F_ip1 >= 0) then
!!$         ip1 = ip1_all(zp1,kind)
         ip1 = priv_ip1_all(ip1)
         istat = fstinl(F_fileid,ni1,nj1,nk1,RMN_ANY_DATE,RMN_ANY_ETK,  &
              ip1,F_ip2,F_ip3,typvar_S,F_nomvar, &
              keylist,nkeys,NMAX)
         if (F_ip1 == 0) then
            ip1 = F_ip1
            istat = fstinl(F_fileid,ni1,nj1,nk1,RMN_ANY_DATE,RMN_ANY_ETK,  &
                 ip1,F_ip2,F_ip3,typvar_S,F_nomvar, &
                 keys2,nkeys2,NMAX)
            do i=1,nkeys2
               if (nkeys >= NMAX) exit
               nkeys = nkeys + 1
               keylist(nkeys) = keys2(i)
            enddo
         endif
      else
         ip1 = F_ip1
         istat = fstinl(F_fileid,ni1,nj1,nk1,RMN_ANY_DATE,RMN_ANY_ETK,  &
              ip1,F_ip2,F_ip3,typvar_S,F_nomvar, &
              keylist,nkeys,NMAX)
      endif
      if (.not.RMN_IS_OK(istat) .or. nkeys == 0 .or. nkeys == NMAX) then
         write(msg_S,'(a,3i10,a,i12,a)') '(fst) Not Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=',F_ip1,F_ip2,F_ip3,'] [key=',F_key,']'
         call msg(MSG_DEBUG,msg_S)
         return
      endif

      dist64 = int(datevfuzz, kind=INT64)
      if (datevfuzz /= huge(datevfuzz)) &
           dist64 = int(datevfuzz+min_time_res, kind=INT64)
      do i=1,nkeys
         istat = fstprm(keylist(i), dateo,deet,npas, ni1,nj1,nk1, &
              nbits, datyp, ip1, ip2, ip3, &
              typvar1_S, nomvar_S, etiket_S, &
              grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
              ubc, extra1, extra2, extra3)
         if (.not.RMN_IS_OK(istat)) then
            F_key = RMN_ERR
            call msg(MSG_WARNING,'(fst_find) problem getting record meta')
            return
         endif
         nhours_8 = (DBLE(deet)*DBLE(npas))/SEC_PER_HR
         call incdatr(datev0,dateo,nhours_8)
         call difdatr(searchdate,datev0,nhours_8) !#TODO: what if searchdate=-1
         dist064 = nint(real(nhours_8*SEC_PER_HR), kind=INT64)
         select case(fuzzopr)
         case(FST_FIND_GE)
            dist064 = -dist064
         case(FST_FIND_NEAR)
            dist064 = abs(dist064)
         !case(FST_FIND_LE)
         end select
         if (dist064 >= 0 .and. dist064 < dist64) then
            F_key  = keylist(i)
            dist64 = dist064
            datev  = datev0
         endif
      end do
      if (typvar_S /= RMN_ANY_TYP .and. present(F_typvar_S)) then
         F_typvar_S =  typvar1_S
      endif

      if (RMN_IS_OK(F_key)) then
         F_datev = datev
         write(msg_S,'(a,3i10,a,i12,a)') '(fst) Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=',F_ip1,F_ip2,F_ip3, &
              '] [typvar='//typvar1_S(1:2)//'] [key=',F_key,']'
      else
         write(msg_S,'(a,3i10,a,i12,a)') '(fst) Not Found: '//nomvar_S(1:4)// &
              ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=',F_ip1,F_ip2,F_ip3, &
              '] [typvar='//typvar_S(1:2)//'] [key=',F_key,']'
      endif
      call msg(MSG_DEBUG,msg_S)
      !---------------------------------------------------------------------
      return
   end function fst_find_0


   !/@
   function fst_find_vect(F_key1, F_key2, F_fileid, F_nomvar1, F_nomvar2,  &
        F_datev, F_ip1, F_ip2, F_ip3, F_datevfuzz, F_fuzzopr, F_typvar_S) &
        result(F_istat)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, intent(out) :: F_key1, F_key2
      character(len=*), intent(in) :: F_nomvar1, F_nomvar2
      integer, intent(inout) :: F_datev
      integer, intent(in) :: F_fileid, F_ip1, F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S  !#TODO:inout
      !@return
      integer :: F_istat
      !@author
      !@revision
      !@/
      character(len=2) :: typvar_S
      character(len=4) :: nomvar1_S, nomvar2_S
      character(len=128) :: msg_S
      integer :: istat, datevfuzz, fuzzopr, ip1, ip2, ip3
      !---------------------------------------------------------------------
      nomvar1_S = F_nomvar1
      nomvar2_S = F_nomvar2
      write(msg_S, '(a,3i10,a)') &
           '(fst) Find vect, looking for: '//nomvar1_S(1:4)//' + '//nomvar2_S(1:4)// &
           ' [datev='//cmcdate_toprint(F_datev)//'] [ip123=', F_ip1, F_ip2, F_ip3, ']'
      call msg(MSG_DEBUG, msg_S)
      F_istat = RMN_ERR
      F_key1 = RMN_ERR
      F_key2 = RMN_ERR
      if (F_fileid <= 0) return
      datevfuzz = 0
      fuzzopr = FST_FIND_NEAR
      typvar_S = RMN_ANY_TYP
      if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
      if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
      if (present(F_typvar_S)) typvar_S = F_typvar_S

      F_key1 = fst_find(F_fileid, F_nomvar1, F_datev, F_ip1, F_ip2, F_ip3,  &
           datevfuzz, fuzzopr, typvar_S)
      if (RMN_IS_OK(F_key1)) then
         datevfuzz = 0
         fuzzopr = FST_FIND_NEAR
         istat = fst_getmeta(F_key1, F_ip1=ip1, F_ip2=ip2, F_ip3=ip3, &
              F_typvar_S=typvar_S)
         F_key2 = fst_find(F_fileid, F_nomvar2, F_datev, F_ip1, F_ip2, F_ip3,  &
              datevfuzz, fuzzopr, typvar_S)
         if (RMN_IS_OK(F_key2)) F_istat = RMN_OK
      endif
      !---------------------------------------------------------------------
      return
   end function fst_find_vect


   !/@
   function fst_find_3d_0(F_keys1, F_fileid, F_nomvar1, &
        F_datev, F_ip1s, F_ip2, F_ip3, F_datevfuzz, F_fuzzopr, F_typvar_S, &
        F_vgrid, F_lvltyp_S) result(F_nkeys)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, pointer :: F_keys1(:)  !# out
      integer, intent(in) :: F_fileid
      character(len=*), intent(in) :: F_nomvar1
      integer, intent(inout) :: F_datev
      integer, pointer :: F_ip1s(:)  !# inout
      integer, intent(in) :: F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S  !#TODO:inout
      type(vgrid_descriptor), intent(out), optional :: F_vgrid
      character(len=*), intent(out), optional :: F_lvltyp_S
      !@return
      integer :: F_nkeys
      !@/
      integer, pointer :: ip1s(:)
      integer :: key, istat, i, datevfuzz, fuzzopr, datev, vb0, ikind
      real :: zp1
      type(vgrid_descriptor) :: vgrid
      character(len=32) :: lvltyp_S, dummy_S
      character(len=2) :: typvar_s
      logical :: needvgd_L, needdiag_L
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) find_3d_0 [BEGIN]')
      F_nkeys = RMN_ERR
      istat = RMN_ERR
      datevfuzz = 0
      fuzzopr = FST_FIND_NEAR
      typvar_S = RMN_ANY_TYP
      if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
      if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_lvltyp_S)) F_lvltyp_S = ' '
      if (present(F_vgrid)) call vgrid_nullify(F_vgrid) 
      datev = F_datev
      if (datevfuzz > 0) datev = -1
      needvgd_L = .not.associated(F_ip1s)
      if (.not.needvgd_L) needvgd_L = (F_ip1s(1) < 0 .and. size(F_ip1s) > 1)
      if (.not.needvgd_L) needdiag_L = &
           (any(F_ip1s(1) == (/FST_FIND_DIAG_M, FST_FIND_DIAG_T/)))
      IF_VGRID: if (needvgd_L) then
         key = fst_find_0(F_fileid, F_nomvar1, datev, RMN_ANY_I, F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S) !# Only needed it datev = RMN_ANY_DATE
         istat = RMN_ERR
         if (key >= 0) istat = vgrid_from_file(&
              F_fileid, F_nomvar1, datev, vgrid, F_ip1s, lvltyp_S)
         if (RMN_IS_OK(istat)) then
            if (present(F_lvltyp_S)) F_lvltyp_S = lvltyp_S
            if (present(F_vgrid)) F_vgrid = vgrid
         endif
      elseif (needdiag_L .or. present(F_vgrid) .or. present(F_lvltyp_S)) then
         key = fst_find_0(F_fileid, F_nomvar1, datev, RMN_ANY_I, F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S)
         IF_KEY_OK: if (key >= 0) then
            nullify(ip1s)
            istat = vgrid_from_file(F_fileid, F_nomvar1, datev, vgrid, &
                 ip1s, lvltyp_S)
            IF_VGRID_OK: if (RMN_IS_OK(istat)) then
               if (present(F_lvltyp_S)) F_lvltyp_S = lvltyp_S
               if (present(F_vgrid)) F_vgrid = vgrid
               IF_DIAG: if (needdiag_L) then
                  call msg_verbosity_get(vb0)
                  call msg_verbosity(MYMSG_QUIET)
                  if (F_ip1s(1) == FST_FIND_DIAG_M) then
                     istat = vgd_get(vgrid, 'DIPM', F_ip1s(1))
                  else
                     istat = vgd_get(vgrid, 'DIPT', F_ip1s(1))
                  endif
                  call msg_verbosity(vb0)
               endif IF_DIAG
            endif IF_VGRID_OK
            IF_DIAG2: if (needdiag_L .and. &
                 .not.(RMN_IS_OK(istat) .and.  F_ip1s(1) >= 0)) then
               call msg(MSG_INFO,'(fst) find: Problem getting vgrid '// &
                    'diag for '//trim(F_nomvar1)//' - trying hyb 1')
               zp1 = 1.
               ikind = RMN_CONV_HY
               dummy_S = ' '
               call convip_plus(F_ip1s(1), zp1, ikind, RMN_CONV_P2IPNEW, &
                    dummy_S, .not.RMN_CONV_USEFORMAT_L)
               istat = RMN_ERR
               if (F_ip1s(1) >= 0) istat = RMN_OK
            endif IF_DIAG2
         endif IF_KEY_OK
      else
         istat = RMN_OK
      endif IF_VGRID
      if (.not.RMN_IS_OK(istat)) return
      if (associated(F_keys1)) then
         if (size(F_keys1) < size(F_ip1s)) then
            call msg(MSG_WARNING, '(fst) find - provided key list too small')
            return
         endif
      else
         allocate(F_keys1(size(F_ip1s)))
      endif

      F_nkeys = 0
      F_keys1 = RMN_ERR
      do i = 1, size(F_ip1s)
         if (i > 1 .and. F_ip1s(i) < 0) exit
         F_keys1(i) = fst_find_0(F_fileid, &
              F_nomvar1, F_datev, F_ip1s(i), F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S)
         F_nkeys = F_nkeys + 1
         if (RMN_IS_OK(F_keys1(i))) then
            datevfuzz = 0
            fuzzopr = FST_FIND_NEAR
            !#TODO: should we make sure all levels have the same typvar_S (if == "")
!!$            F_keys1(F_nkeys) = F_keys1(i)
!!$            F_ip1s(F_nkeys) = F_ip1s(i)
!!$            if (i /= F_nkeys) then
!!$               F_keys1(i) = RMN_ERR
!!$               F_ip1s(i) = RMN_ERR
!!$            endif
         ! else
         !    F_keys1(i) = RMN_ERR
         endif
      enddo
      if (maxval(F_keys1) < 0) F_nkeys = 0
!!$      if (present(F_typvar_S)) F_typvar_S = typvar_S !#TODO
      call msg(MSG_DEBUG, '(fst) find_3d_0 [END]')
      !---------------------------------------------------------------------
      return
   end function fst_find_3d_0


   !/@
   function fst_find_3d_vect(F_keys1, F_keys2, F_fileid, &
        F_nomvar1, F_nomvar2, F_datev, F_ip1s, F_ip2, F_ip3,  &
        F_datevfuzz, F_fuzzopr, F_typvar_S, F_vgrid, F_lvltyp_S) result(F_nkeys)
      implicit none
      !@objective Try to find rec key for given field params
      !@arguments
      integer, intent(in) :: F_fileid
      integer, pointer :: F_keys1(:), F_keys2(:)  !# out
      character(len=*), intent(in) :: F_nomvar1, F_nomvar2
      integer, intent(inout) :: F_datev
      integer, pointer :: F_ip1s(:)  !# inout
      integer, intent(in) :: F_ip2, F_ip3
      integer, intent(in), optional :: F_datevfuzz, F_fuzzopr
      character(len=*), intent(in), optional :: F_typvar_S  !#TODO:inout
      type(vgrid_descriptor), intent(out), optional :: F_vgrid  !# out
      character(len=*), intent(out), optional :: F_lvltyp_S
      !@return
      integer :: F_nkeys
      !@author
      !@revision
      !@/
      integer, pointer :: ip1s(:)
      integer :: key, istat, i, datevfuzz, fuzzopr, datev, vb0, ikind
      real :: zp1
      type(vgrid_descriptor) :: vgrid
      character(len=32) :: lvltyp_S, dummy_S
      character(len=2) :: typvar_s
      logical :: needvgd_L, needdiag_L
      !---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) find_3d_vect [BEGIN]')
      F_nkeys = RMN_ERR
      istat = RMN_ERR
      datevfuzz = 0
      fuzzopr = FST_FIND_NEAR
      typvar_S = RMN_ANY_TYP
      if (present(F_datevfuzz)) datevfuzz = F_datevfuzz
      if (present(F_fuzzopr)) fuzzopr = F_fuzzopr
      if (present(F_typvar_S)) typvar_S = F_typvar_S
      if (present(F_lvltyp_S)) F_lvltyp_S = ' '
      if (present(F_vgrid)) call vgrid_nullify(F_vgrid) 
      datev = F_datev
      if (datevfuzz > 0) datev = -1
      needvgd_L = .not.associated(F_ip1s)
      if (.not.needvgd_L) needvgd_L = (F_ip1s(1) < 0 .and. size(F_ip1s) > 1)
      if (.not.needvgd_L) needdiag_L = &
           (any(F_ip1s(1) == (/FST_FIND_DIAG_M, FST_FIND_DIAG_T/)))
      IF_VGRID: if (needvgd_L) then
         key = fst_find_0(F_fileid, F_nomvar1, datev, RMN_ANY_I, F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S) !# Only needed it datev = RMN_ANY_DATE
         istat = RMN_ERR
         if (key >= 0) istat = vgrid_from_file( &
              F_fileid, F_nomvar1, datev, vgrid, F_ip1s, lvltyp_S)
         if (RMN_IS_OK(istat)) then
            if (present(F_lvltyp_S)) F_lvltyp_S = lvltyp_S
            if (present(F_vgrid)) F_vgrid = vgrid
         endif
      elseif (needdiag_L .or. present(F_vgrid) .or. present(F_lvltyp_S)) then
         key = fst_find_0(F_fileid, F_nomvar1, datev, RMN_ANY_I, F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S)
         IF_KEY_OK: if (key >= 0) then
            nullify(ip1s)
            istat = vgrid_from_file(F_fileid, F_nomvar1, datev, vgrid, &
                 ip1s, lvltyp_S)
            IF_VGRID_OK: if (RMN_IS_OK(istat)) then
               if (present(F_lvltyp_S)) F_lvltyp_S = lvltyp_S
               if (present(F_vgrid)) F_vgrid = vgrid
               IF_DIAG: if (needdiag_L) then
                  call msg_verbosity_get(vb0)
                  call msg_verbosity(MYMSG_QUIET)
                  if (F_ip1s(1) == FST_FIND_DIAG_M) then
                     istat = vgd_get(vgrid, 'DIPM', F_ip1s(1))
                  else
                     istat = vgd_get(vgrid, 'DIPT', F_ip1s(1))
                  endif
                  call msg_verbosity(vb0)
               endif IF_DIAG
            endif IF_VGRID_OK
            IF_DIAG2: if (needdiag_L .and. &
                 .not.(RMN_IS_OK(istat) .and.  F_ip1s(1) >= 0)) then
               call msg(MSG_INFO,'(fst) find: Problem getting vgrid '// &
                    'diag for '//trim(F_nomvar1)//' - trying hyb 1')
               zp1 = 1.
               ikind = RMN_CONV_HY
               dummy_S = ' '
               call convip_plus(F_ip1s(1), zp1, ikind, RMN_CONV_P2IPNEW, &
                    dummy_S, .not.RMN_CONV_USEFORMAT_L)
               istat = RMN_ERR
               if (F_ip1s(1) >= 0) istat = RMN_OK
            endif IF_DIAG2
         endif IF_KEY_OK
      else
         istat = RMN_OK
      endif IF_VGRID
      if (.not.RMN_IS_OK(istat)) return
      if (associated(F_keys1)) then
         if (size(F_keys1) < size(F_ip1s)) then
            call msg(MSG_WARNING, '(fst) find - provided key list too small')
            return
         endif
      else
         allocate(F_keys1(size(F_ip1s)), stat=istat)
      endif
      if (associated(F_keys2)) then
         if (size(F_keys2) < size(F_ip1s)) then
            call msg(MSG_WARNING, '(fst) find - provided key list too small')
            return
         endif
      else
         allocate(F_keys2(size(F_ip1s)), stat=istat)
      endif

      F_nkeys = 0
      F_keys1 = RMN_ERR
      F_keys2 = RMN_ERR
      do i = 1, size(F_ip1s)
         if (i > 1 .and. F_ip1s(i) < 0) exit
         istat = fst_find_vect(F_keys1(i), F_keys2(i), F_fileid, &
              F_nomvar1, F_nomvar2, F_datev, F_ip1s(i), F_ip2, F_ip3, &
              datevfuzz, fuzzopr, typvar_S)
         F_nkeys = F_nkeys + 1
         if (RMN_IS_OK(F_keys1(i)) .and. &
              (F_nomvar2 == ' ' .or. RMN_IS_OK(F_keys2(i)))) then
            datevfuzz = 0
            fuzzopr = FST_FIND_NEAR
            !#TODO: should we make sure all levels have the same typvar_S (if == "")
         else
            F_keys1(i) = RMN_ERR
            F_keys2(i) = RMN_ERR
         endif
      enddo
      if (maxval(F_keys1) < 0) F_nkeys = 0
      call msg(MSG_DEBUG, '(fst) find_3d_vect [END]')
      !---------------------------------------------------------------------
      return
   end function fst_find_3d_vect


   !/@
   function fst_read_3d_r4(F_key,F_data,F_fileid,F_gridid,&
        F_nomvar_S,F_etiket_S,F_dateo,F_deet,F_npas,F_ip1,F_ip2,F_ip3,&
        F_typvar_S, F_realloc_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_key
      real,pointer :: F_data(:,:,:)
      integer,intent(in),optional :: F_fileid
      integer,intent(out),optional :: F_gridid
      character(len=*),intent(out),optional :: F_nomvar_S, F_etiket_S,F_typvar_S
      integer,intent(out),optional :: F_dateo,F_deet,F_npas,F_ip1, F_ip2, F_ip3
      logical,intent(in),optional :: F_realloc_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: istat
      character(len=2) :: typvar_S
      character(len=12) :: etiket_S
      character(len=4) :: nomvar_S
      character(len=2) :: grtyp_S
      character(len=128) :: msg_S
      integer :: ni1,nj1,nk1, &
           dateo,deet,npas,nbits,datyp,ip1,ip2,ip3,&
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
      logical :: realloc_L
      real:: nhours
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) read_3d_r4 [BEGIN]')
      F_istat = RMN_ERR
      if (present(F_gridid)) F_gridid = RMN_ERR
      if (F_key < 0) return
      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      istat = fstprm(F_key, dateo,deet,npas, ni1,nj1,nk1, &
           nbits, datyp, ip1, ip2, ip3, &
           typvar_S, nomvar_S, etiket_S, &
           grtyp_S(1:1), ig1, ig2, ig3, ig4, swa, lng, dltf, &
           ubc, extra1, extra2, extra3)
      if (.not.RMN_IS_OK(istat) .or. ni1<1 .or. nj1<1 .or. nk1<1) then
         call msg(MSG_WARNING,'(fst) read - Cannot get field dims')
         return
      endif

!!$      print *,'(fst_read) found:',trim(nomvar_S),ni1,nj1,':',ip1, ip2, ip3,':',grtyp_S(1:1),ig1, ig2, ig3, ig4

      istat = fst_checkalloc(F_data, ni1, nj1, nk1, realloc_L)
      if (.not.RMN_IS_OK(istat)) return

      F_istat = fstluk(F_data,F_key,ni1,nj1,nk1)

      nhours = float(deet)*float(npas)/3600.
      if (RMN_IS_OK(F_istat)) then
         write(msg_S, '(a,f7.1,a,3i10,a)') '(fst) Read: '//nomvar_S(1:4)// &
              ' [dateo='//cmcdate_toprint(dateo)//' +',nhours,'h] [ip123=', &
              ip1,ip2,ip3, ', tv='//trim(typvar_S)//']'
         call msg(MSG_INFO, msg_S)
      else
         write(msg_S, '(a,f7.1,a,3i10,a)') &
              '(fst) Problem Reading: '//nomvar_S(1:4)// &
              ' [dateo='//cmcdate_toprint(dateo)//' +',nhours,'h] [ip123=', &
              ip1,ip2,ip3, ', tv='//trim(typvar_S)//']'
         call msg(MSG_WARNING, msg_S)
      endif

      if (present(F_gridid).and.present(F_fileid)) then
!!$         print *,'(fst_read) looking for gridL:',grtyp_S(1:1),ig1,ig2,ig3,ig4 ; call flush(6)
         if (any(grtyp_S(1:1) == (/'u','U'/))) then
            ni1 = -1 ; nj1 = -1
         endif
         F_gridid = ezqkdef(ni1,nj1, grtyp_S(1:1), ig1, ig2, ig3, ig4, F_fileid)
      endif
      if (present(F_nomvar_S)) F_nomvar_S = nomvar_S
      if (present(F_etiket_S)) F_etiket_S = etiket_S
      if (present(F_typvar_S)) F_typvar_S = typvar_S
      if (present(F_dateo)) F_dateo = dateo
      if (present(F_deet)) F_deet = deet
      if (present(F_npas)) F_npas = npas
      if (present(F_ip1)) F_ip1 = ip1
      if (present(F_ip2)) F_ip2 = ip2
      if (present(F_ip3)) F_ip3 = ip3

!!$      write(msg_S,'(a,3i10,a)') '(fst) Read: '//nomvar_S(1:4)//' [dateo='//cmcdate_toprint(dateo)//'] [ip123=',ip1,ip2,ip3,']'
!!$      call msg(MSG_INFO, msg_S)
!!$      call msg(MSG_DEBUG, '(fst) Read [END]')
      call msg(MSG_DEBUG, '(fst) read_3d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_read_3d_r4


   !/@
   function fst_read_3d_r4_vect(F_key1, F_key2, F_data1, F_data2, F_fileid, &
        F_gridid, F_nomvar1_S, F_nomvar2_S, F_etiket_S, F_dateo, F_deet,  &
        F_npas, F_ip1, F_ip2, F_ip3, F_typvar_S, F_realloc_L) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer, intent(in) :: F_key1, F_key2
      real, pointer :: F_data1(:,:,:), F_data2(:,:,:)
      integer, intent(in), optional :: F_fileid
      integer, intent(out), optional :: F_gridid
      character(len=*), intent(out), optional :: F_nomvar1_S, F_nomvar2_S, F_etiket_S,F_typvar_S
      integer, intent(out), optional :: F_dateo, F_deet, F_npas, F_ip1, F_ip2, F_ip3
      logical, intent(in), optional :: F_realloc_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: istat
      character(len=2) :: typvar_S
      character(len=12) :: etiket_S
      character(len=4) :: nomvar1_S, nomvar2_S
      integer :: dateo, deet, npas, ip1, ip2, ip3
      logical :: realloc_L
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) read_3d_r4_vect [BEGIN]')
      F_istat = RMN_ERR
      if (present(F_gridid)) F_gridid = RMN_ERR
      if (F_key1 < 0 .or. F_key2 < 0) return
      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      if (present(F_gridid).and.present(F_fileid)) then
         istat = fst_read_3d_r4(F_key1, F_data1, F_fileid, F_gridid, &
              nomvar1_S, etiket_S, dateo, deet, npas, &
              ip1, ip2, ip3, typvar_S, F_realloc_L=realloc_L)
      else
         istat = fst_read_3d_r4(F_key1, F_data1, &
              F_nomvar_S=nomvar1_S, F_etiket_S=etiket_S, F_dateo=dateo, &
              F_deet=deet, F_npas=npas, F_ip1=ip1, F_ip2=ip2, F_ip3=ip3, &
              F_typvar_S=typvar_S, F_realloc_L=realloc_L)
      endif
      if (RMN_IS_OK(istat)) then
         F_istat = fst_read_3d_r4(F_key2, F_data2, &
              F_nomvar_S=nomvar2_S, F_etiket_S=etiket_S, F_dateo=dateo, &
              F_deet=deet, F_npas=npas, F_ip1=ip1, F_ip2=ip2, F_ip3=ip3, &
              F_typvar_S=typvar_S, F_realloc_L=realloc_L)
      endif

      if (present(F_nomvar1_S)) F_nomvar1_S = nomvar1_S
      if (present(F_nomvar2_S)) F_nomvar2_S = nomvar2_S
      if (present(F_etiket_S)) F_etiket_S = etiket_S
      if (present(F_typvar_S)) F_typvar_S = typvar_S
      if (present(F_dateo)) F_dateo = dateo
      if (present(F_deet)) F_deet = deet
      if (present(F_npas)) F_npas = npas
      if (present(F_ip1)) F_ip1 = ip1
      if (present(F_ip2)) F_ip2 = ip2
      if (present(F_ip3)) F_ip3 = ip3
      call msg(MSG_DEBUG, '(fst) read_3d_r4_vect [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_read_3d_r4_vect


   !/@
   function fst_rdhint_3d_r4(F_data, F_status, F_keylist, F_hintlist_S, &
        F_fileids, F_outgridid, F_coregridid, F_realloc_L) &
        result(F_istat)
      implicit none
      !@objective
      !@arguments
      real, pointer :: F_data(:,:,:)
      integer, intent(out) :: F_status(:)
      integer, intent(in) :: F_keylist(:)
      character(len=*), intent(in) :: F_hintlist_S(:)
      integer, intent(in) :: F_fileids(:)
      integer, intent(in) :: F_outgridid
      integer, intent(in), optional :: F_coregridid
      logical,intent(in), optional :: F_realloc_L
      !#TODO: alternate itf with hgridis_S so we have min, max as well
      !@author
      !@return
      integer :: F_istat
      !@/
      real, pointer :: indata(:,:,:)
      integer :: istat, coregridid, nij(2), nkeys, nhint, nfids, ikey, ingridid
      character(len=12) :: nomvar_S
      logical :: realloc_L
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) rdhint_3d_r4 [BEGIN]')
      F_istat = RMN_OK

      coregridid = F_outgridid
      if (present(F_coregridid)) coregridid = F_coregridid
      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      nkeys = size(F_keylist)
      nhint = size(F_hintlist_S)
      nfids = size(F_fileids)

      if (size(F_status) < nkeys) then
         F_istat = RMN_ERR
         call msg(MSG_WARNING, '(fst) rdhint: status array too small')
         return
      endif

      F_istat = ezgrid_params(F_outgridid, nij)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING, '(fst) rdhint: Problem getting grid params')
         return
      endif

      F_istat = fst_checkalloc(F_data, nij(1), nij(2), nkeys, realloc_L)
      if (.not.RMN_IS_OK(F_istat)) return

      F_status = RMN_ERR
      nullify(indata)
      DO_NKEYS: do ikey = 1, nkeys
 
         if (F_keylist(ikey) < 0) cycle

         ingridid = -1
         istat =  fst_read_3d_r4(F_keylist(ikey), indata, &
              F_fileids(min(ikey, nfids)), ingridid, nomvar_S, &
              F_realloc_L=realloc_L)
         if (.not.RMN_IS_OK(istat)) then
            F_status(ikey) = RMN_ERR
            call msg(MSG_WARNING, '(fst) rdhint: Problem reading: '//trim(nomvar_S))
            cycle
         endif
         !#TODO: allow for optional stats of read field before interp

         F_status(ikey) = hinterp4yy2d(F_data, indata, ikey, ingridid, &
              F_outgridid, coregridid, F_hintlist_S(min(ikey,nhint)), &
              nomvar_S)
         if (.not.RMN_IS_OK(F_status(ikey))) then
            call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in hinterp4yy2d for: '//trim(nomvar_S))
            cycle
         endif

      enddo DO_NKEYS
      if (associated(indata)) deallocate(indata, stat=istat)
      F_istat = maxval(F_status(1:nkeys))
      call msg(MSG_DEBUG, '(fst) rdhint_3d_r4 [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_rdhint_3d_r4


   !/@
   function fst_rdhint_3d_r4_vect(F_data1, F_data2, F_status, F_keys1, F_keys2,&
        F_hintlist_S, F_fileids, F_outgridid, F_coregridid, &
        F_realloc_L) result(F_istat)
      implicit none
      !@objective
      !@arguments
      real, pointer :: F_data1(:,:,:), F_data2(:,:,:)
      integer, intent(out) :: F_status(:)
      integer, intent(in) :: F_keys1(:), F_keys2(:)
      character(len=*), intent(in) :: F_hintlist_S(:)
      integer, intent(in) :: F_fileids(:)
      integer, intent(in) :: F_outgridid
      integer, intent(in), optional :: F_coregridid
      logical,intent(in), optional :: F_realloc_L
      !#TODO: alternate itf with hgridis_S so we have min, max as well
      !@author
      !@return
      integer :: F_istat
      !@/
      real, pointer :: indata1(:,:,:), indata2(:,:,:)
      integer :: istat, coregridid, nij(2), nkeys, nhint, nfids, ikey, ingridid
      character(len=12) :: nomvar1_S, nomvar2_S
      logical :: realloc_L
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) rdhint_3d_r4_vect [BEGIN]')
      F_istat = RMN_OK

      coregridid = F_outgridid
      if (present(F_coregridid)) coregridid = F_coregridid
      realloc_L = .false.
      if (present(F_realloc_L)) realloc_L = F_realloc_L

      nkeys = size(F_keys1)
      nhint = size(F_hintlist_S)
      nfids = size(F_fileids)

      if (nkeys /= size(F_keys2)) then
         call msg(MSG_WARNING, '(fst) rdhint: incompatible list size')
         F_istat = RMN_ERR
         return
      endif

      if (size(F_status) < nkeys) then
         F_istat = RMN_ERR
         call msg(MSG_WARNING, '(fst) rdhint: status array too small')
         return
      endif

      F_istat = ezgrid_params(F_outgridid, nij)
      if (.not.RMN_IS_OK(F_istat)) then
         call msg(MSG_WARNING, '(fst) rdhint: Problem getting grid params')
         return
      endif

      F_istat = fst_checkalloc(F_data1, nij(1), nij(2), nkeys, realloc_L)
      if (.not.RMN_IS_OK(F_istat)) return
      F_istat = fst_checkalloc(F_data2, nij(1), nij(2), nkeys, realloc_L)
      if (.not.RMN_IS_OK(F_istat)) return

      F_status = RMN_ERR
      nullify(indata1, indata2)
      DO_IKEY: do ikey = 1, nkeys

         ingridid = -1
         istat =  fst_read_3d_r4(F_keys1(ikey), indata1, &
              F_fileids(min(ikey, nfids)), ingridid, nomvar1_S, &
              F_realloc_L=realloc_L)
         if (.not.RMN_IS_OK(istat)) then
            F_status(ikey) = RMN_ERR
            call msg(MSG_WARNING, '(fst) rdhint: Problem reading: '//trim(nomvar1_S))
            cycle
         endif
         istat =  fst_read_3d_r4(F_keys2(ikey), indata2, &
              F_fileids(min(ikey, nfids)), ingridid, nomvar2_S, &
              F_realloc_L=realloc_L)
         if (.not.RMN_IS_OK(istat)) then
            F_status(ikey) = RMN_ERR
            call msg(MSG_WARNING, '(fst) rdhint: Problem reading: '//trim(nomvar2_S))
            cycle
         endif
         !#TODO: allow for optional stats of read field before interp

         F_status(ikey) = hinterp4yy2d(F_data1, F_data2, indata1, indata2, ikey, &
              ingridid, F_outgridid, coregridid, &
              F_hintlist_S(min(ikey, nhint)), nomvar1_S, nomvar2_S)
         if (.not.RMN_IS_OK(F_status(ikey))) then
            call msg(MSG_WARNING, '(fstmpio) rdhint: Problem in hinterp4yy2d for: '//trim(nomvar1_S)//' + '//trim(nomvar2_S))
            cycle
         endif

      enddo DO_IKEY
      if (associated(indata1)) deallocate(indata1, stat=istat)
      if (associated(indata2)) deallocate(indata2, stat=istat)
      F_istat = maxval(F_status(1:nkeys))
      call msg(MSG_DEBUG, '(fst) rdhint_3d_r4_vect [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_rdhint_3d_r4_vect


   !/@
   function fst_getmeta(F_key,F_nomvar_S,F_dateo,F_deet,F_npas,&
        F_ip1, F_ip2, F_ip3, F_etiket_S,F_typvar_S) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_key
      character(len=*),intent(out),optional :: F_nomvar_S, F_etiket_S,F_typvar_S
      integer,intent(out),optional :: F_dateo,F_deet,F_npas,F_ip1, F_ip2, F_ip3
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: istat
      character(len=1) :: grtyp_S
      character(len=2) :: typvar_S
      character(len=4) :: nomvar_S
      character(len=12):: etiket_S
      integer :: ni1,nj1,nk1, &
           dateo,deet,npas,nbits,datyp,ip1,ip2,ip3,&
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) getmet [BEGIN]')
      F_istat = RMN_ERR
      if (F_key < 0) return

      istat = fstprm(F_key, dateo,deet,npas, ni1,nj1,nk1, &
           nbits, datyp, ip1, ip2, ip3, &
           typvar_S, nomvar_S, etiket_S, &
           grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
           ubc, extra1, extra2, extra3)
      if (.not.RMN_IS_OK(istat) .or. ni1<1 .or. nj1<1 .or. nk1<1) return

      if (present(F_nomvar_S)) F_nomvar_S = nomvar_S
      if (present(F_etiket_S)) F_etiket_S = etiket_S
      if (present(F_typvar_S)) F_typvar_S = typvar_S
      if (present(F_dateo)) F_dateo = dateo
      if (present(F_deet)) F_deet = deet
      if (present(F_npas)) F_npas = npas
      if (present(F_ip1)) F_ip1 = ip1
      if (present(F_ip2)) F_ip2 = ip2
      if (present(F_ip3)) F_ip3 = ip3
      call msg(MSG_DEBUG, '(fst) getmet [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_getmeta


   !/@
   function fst_get_hgridid(F_fileid, F_key) result(F_gridid)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid,F_key
      !@author
      !@return
      integer :: F_gridid
      !@/
      integer :: istat
      character(len=1) :: grtyp_S
      character(len=2) :: typvar_S
      character(len=4) :: nomvar_S
      character(len=12):: etiket_S
      integer :: ni1,nj1,nk1, &
           dateo,deet,npas,nbits,datyp,ip1,ip2,ip3,&
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) get_hgridid [BEGIN]')
      F_gridid = RMN_ERR
      if (F_fileid <= 0 .or. F_key < 0) return

      istat = fstprm(F_key, dateo,deet,npas, ni1,nj1,nk1, &
           nbits, datyp, ip1, ip2, ip3, &
           typvar_S, nomvar_S, etiket_S, &
           grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
           ubc, extra1, extra2, extra3)

      if (any(grtyp_S(1:1) == (/'u','U'/))) then
         ni1 = -1 ; nj1 = -1
      endif
      F_gridid = ezqkdef(ni1,nj1, grtyp_S, ig1, ig2, ig3, ig4, F_fileid)
      call msg(MSG_DEBUG, '(fst) get_hgridid [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_get_hgridid


   !/@
   function fst_get_vgrid(F_fileid,F_key,F_vgrid,F_ip1s,F_lvltyp_S) result(F_istat)
      implicit none
      !@objective 
      !@arguments
      integer,intent(in) :: F_fileid,F_key
      type(vgrid_descriptor),intent(out) :: F_vgrid
      integer,pointer :: F_ip1s(:)
      character(len=*),intent(out) :: F_lvltyp_S
      !@author Ron McTaggartCowan, Aug 2012
      !@return
      integer :: F_istat
      !@/
      integer,parameter :: NMAX = 9999
      integer :: istat,keylist(NMAX),nkeys,k
      character(len=1) :: grtyp_S
      character(len=2) :: typvar_S
      character(len=4) :: nomvar_S
      character(len=12):: etiket_S
      integer :: ni1,nj1,nk1, &
           dateo,datev,deet,npas,nbits,datyp,ip1,ip2,ip3,&
           ig1, ig2, ig3, ig4, swa, lng, dltf, ubc, extra1, extra2, extra3
      ! ---------------------------------------------------------------------
      call msg(MSG_DEBUG, '(fst) get_vgrid [BEGIN]')
      F_istat = RMN_ERR
      call vgrid_nullify(F_vgrid) 
      F_lvltyp_S = ' '
      if (F_fileid <= 0 .or. F_key < 0) return

      istat = fstprm(F_key, dateo,deet,npas, ni1,nj1,nk1, &
           nbits, datyp, ip1, ip2, ip3, &
           typvar_S, nomvar_S, etiket_S, &
           grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
           ubc, extra1, extra2, extra3)
      if (.not.RMN_IS_OK(istat)) return

      istat = vgd_new(F_vgrid,unit=F_fileid,ip1=ig1,ip2=ig2)
      ! if (istat /= VGD_OK) &
      !      istat = vgd_new(F_vgrid,unit=F_fileid,ip1=ig4,ip2=-1)
      if (istat /= VGD_OK) &
           istat = vgd_new(F_vgrid,unit=F_fileid,ip1=-1,ip2=-1)
      if (istat /= VGD_OK) return

      F_lvltyp_S = 'M'
      istat = vgd_get(F_vgrid,'VIP'//trim(F_lvltyp_S),F_ip1s)
      if (.not.(istat == VGD_OK .and. associated(F_ip1s))) return
      if (any(ip1 == F_ip1s)) then
         F_istat = size(F_ip1s)
         return
      endif

      F_lvltyp_S = 'T'
      istat = vgd_get(F_vgrid,'VIP'//trim(F_lvltyp_S),F_ip1s)
      if (.not.(istat == VGD_OK .and. associated(F_ip1s))) return
      if (any(ip1 == F_ip1s)) then
         F_istat = size(F_ip1s)
         return
      endif

      F_lvltyp_S = 'SFC'
      call incdatr(datev,dateo,(dble(deet)*dble(npas))/SEC_PER_HR)
      istat = fstinl(F_fileid,ni1,nj1,nk1,datev,'',ip1,-1,-1,'',nomvar_S,keylist,nkeys,size(keylist))
      if (.not.RMN_IS_OK(istat) .or. nkeys < 1) return
      if (associated(F_ip1s)) deallocate(F_ip1s, stat=istat)
      allocate(F_ip1s(nkeys),stat=istat)
      if (.not.associated(F_ip1s)) return
      istat = RMN_OK
      do k=1,nkeys
         istat = min(fstprm(F_key, dateo,deet,npas, ni1,nj1,nk1, &
              nbits, datyp, F_ip1s(k), ip2, ip3, &
              typvar_S, nomvar_S, etiket_S, &
              grtyp_S, ig1, ig2, ig3, ig4, swa, lng, dltf, &
              ubc, extra1, extra2, extra3),istat)
      enddo
      if (.not.RMN_IS_OK(istat)) then
         if (associated(F_ip1s)) deallocate(F_ip1s, stat=istat)
         return
      endif
      !TODO: check that it's sfc (kind 3) levels
      !TODO: sort_unique
      F_istat = size(F_ip1s)
      call msg(MSG_DEBUG, '(fst) get_vgrid [END]')
      ! ---------------------------------------------------------------------
      return
   end function fst_get_vgrid


   !/@
   function fst_checkalloc(F_data, F_ni, F_nj, F_nk, F_realloc_L) &
        result(F_istat)
      implicit none
      !@objective Check F_data dims and optionally alloc/re-alloc
      !@arguments
      real, pointer :: F_data(:,:,:)
      integer, intent(in) :: F_ni, F_nj, F_nk
      logical, intent(in) :: F_realloc_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: istat
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (associated(F_data)) then
         if (size(F_data,1) /= F_ni .or. &
              size(F_data,2) /= F_nj .or. &
              size(F_data,3) /= F_nk) then
            if (F_realloc_L) then
               deallocate(F_data, stat=istat)
               nullify(F_data)
            else
               call msg(MSG_WARNING,'(fst) Wrong dims of Provided Field')
               return
            endif
         endif
      endif
      if (.not.associated(F_data)) then
         allocate(F_data(F_ni, F_nj, F_nk), stat=istat)
         if (istat /= 0) return
      endif
      F_istat = RMN_OK
      ! ---------------------------------------------------------------------
      return
   end function fst_checkalloc


   !/@
   function fst_checkalloc2(F_data, F_mini, F_maxi, F_minj, F_maxj, &
        F_nk, F_realloc_L) result(F_istat)
      implicit none
      !@objective Check F_data dims and optionally alloc/re-alloc
      !@arguments
      real, pointer :: F_data(:,:,:)
      integer, intent(in) :: F_mini, F_maxi, F_minj, F_maxj, F_nk
      logical, intent(in) :: F_realloc_L
      !@author
      !@return
      integer :: F_istat
      !@/
      integer :: istat
      ! ---------------------------------------------------------------------
      F_istat = RMN_ERR
      if (associated(F_data)) then
         if (any(lbound(F_data) /= (/F_mini, F_minj, 1/)) .or. &
              any(ubound(F_data) /= (/F_maxi, F_maxj, F_nk/))) then
            if (F_realloc_L) then
               deallocate(F_data, stat=istat)
               nullify(F_data)
            else
               call msg(MSG_WARNING,'(fst) Wrong shape of Provided Field')
               return
            endif
         endif
      endif
      if (.not.associated(F_data)) then
         allocate(F_data(F_mini:F_maxi, F_minj:F_maxj, F_nk), stat=istat)
         if (istat /= 0) return
      endif
      F_istat = RMN_OK
      ! ---------------------------------------------------------------------
      return
   end function fst_checkalloc2

   !==== Private Functions =================================================


   function priv_ip1_all(F_ip1) result(F_ip1out)
      implicit none
      integer,intent(in) :: F_ip1
      integer :: F_ip1out
!!$      integer,parameter :: NMAXIP1ALL = 1024
!!$      integer,save :: m_ip1all(NMAXIP1ALL,2) = -1
!!$      integer,save :: m_nip1all = 0
!!$      integer :: nn
      character(len=12) :: dummy_S = ' '
      integer :: ip1, ikind
      real :: zp1
      ! -----------------------------------------------------------------
!!$      do nn=1,m_nip1all
!!$         if (m_ip1all(nn,1) == F_ip1) then
!!$            F_ip1out = m_ip1all(nn,2)
!!$            return
!!$         endif
!!$      enddo
      ip1 = F_ip1
      call convip_plus(ip1, zp1, ikind, RMN_CONV_IP2P, dummy_S, .not.RMN_CONV_USEFORMAT_L)
      F_ip1out = ip1_all(zp1,ikind)
!!$      if (m_nip1all < NMAXIP1ALL)  then
!!$         m_nip1all = m_nip1all + 1
!!$         m_ip1all(m_nip1all,1) = F_ip1
!!$         m_ip1all(m_nip1all,2) = F_ip1out
!!$         if (all(m_ip1all(1:m_nip1all,1) /= F_ip1out) &
!!$              .and. m_nip1all < NMAXIP1ALL)  then
!!$            m_nip1all = m_nip1all + 1
!!$            m_ip1all(m_nip1all,1) = F_ip1out
!!$            m_ip1all(m_nip1all,2) = F_ip1out
!!$         endif
!!$      endif
      ! -----------------------------------------------------------------
      return
   end function priv_ip1_all

end module fst_read_mod
