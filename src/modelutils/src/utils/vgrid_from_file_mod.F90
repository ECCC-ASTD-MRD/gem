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

!/@*
module vgrid_from_file_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use iso_c_binding
   use rpn_comm_itf_mod
   use VGrid_Descriptors
   use vgrid_ov, only: vgrid_nullify
   use vgrid_wb
   use sort_mod
   use mu_jdate_mod
   use ptopo_utils
   implicit none
   private
   !@objective 
   !@author Stephane Chamberland,2013-03
   !@description
   ! Public functions
   public :: vgrid_from_file, vgrid_from_file_mpi, vgrid_from_file_rfld_key, &
        vgrid_from_file_rfld_key_mpi
   ! Public parameters
   public :: VGRID_FROM_FILE_NORFLD
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   integer, parameter :: VGRID_FROM_FILE_NORFLD = RMN_ERR - 1
   integer, parameter :: MAXLEV=1024, NLEVTYP=2, IKIND_SURF=3
   integer, parameter :: MYMSG_QUIET = 99

   interface vgrid_from_file
      module procedure vgrid_from_file_4
      module procedure vgrid_from_file_8
   end interface

   interface vgrid_from_file_mpi
      module procedure vgrid_from_file_mpi_4
      module procedure vgrid_from_file_mpi_8
   end interface

   interface vgrid_from_file_rfld_key_mpi
      module procedure vgrid_from_file_rfld_key_mpi_4
      module procedure vgrid_from_file_rfld_key_mpi_8
   end interface

   interface vgrid_from_file_rfld_key
      module procedure vgrid_from_file_rfld_key_4
      module procedure vgrid_from_file_rfld_key_8
   end interface


contains 


   !/@*
   function vgrid_from_file_8(F_unit, F_varname_S, F_jdatev, F_vgrid, &
        F_iplist, F_levtype_S, F_fstkeys, F_sfcRefKey, F_sfcRefKey2) &
        result(F_istat)
      implicit none
      !@objective Get vgrid and ip1 list for specific field in file
      !@author Ron McTaggart-Cowan, 2012-08
      !@revision
      ! 2014-11, S.Chamberland: sorting ip1,keys as in vgrid from files
      !@arguments
      integer,intent(in) :: F_unit                !Input file unit (opened)
      character(len=*),intent(in) :: F_varname_S  !Field name in FST file
      integer(INT64),intent(in) :: F_jdatev     !Valid date (jsec)
      type(vgrid_descriptor),intent(out) :: F_vgrid  !Vgrid of record
      integer,pointer              :: F_iplist(:)  !List of IP1s for the field
      character(len=*),intent(out) :: F_levtype_S  !Short level name
      integer,pointer,optional     :: F_fstkeys(:) !List of field's fst keys
      integer,intent(out),optional :: F_sfcRefKey   !Key for Sfc Ref Field
      integer,intent(out),optional :: F_sfcRefKey2  !Key for Sfc LS Ref Field
      !@return
      integer :: F_istat
      !*@/
      integer :: datev, sfcRefKey2
      !------------------------------------------------------------------
      call vgrid_nullify(F_vgrid) 
      datev = RMN_ANY_DATE
      if (F_jdatev /= MU_JDATE_ANY) datev = jdate_to_cmc(F_jdatev)
      sfcRefKey2 = RMN_ERR
      if (present(F_fstkeys) .and. present(F_sfcRefKey)) then
         F_istat = vgrid_from_file_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, F_fstkeys, F_sfcRefKey, &
              F_sfcRefKey2=sfcRefKey2)
      else if (present(F_fstkeys)) then
         F_istat = vgrid_from_file_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, F_fstkeys)
      else if (present(F_sfcRefKey)) then
         F_istat = vgrid_from_file_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, F_sfcRefKey=F_sfcRefKey, &
              F_sfcRefKey2=sfcRefKey2)
      else
         F_istat = vgrid_from_file_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S)
      endif
      if (present(F_sfcRefKey2)) F_sfcRefKey2 = sfcRefKey2
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_8


   !/@*
   function vgrid_from_file_4(F_unit, F_varname_S, F_datev, F_vgrid, &
        F_iplist, F_levtype_S, F_fstkeys, F_sfcRefKey, F_sfcRefKey2) &
        result(F_istat)
      implicit none
      !@objective Get vgrid and ip1 list for specific field in file
      !@author Ron McTaggart-Cowan, 2012-08
      !@revision
      ! 2014-11, S.Chamberland: sorting ip1,keys as in vgrid from files
      !@arguments
      integer,intent(in) :: F_unit                !Input file unit (opened)
      character(len=*),intent(in) :: F_varname_S  !Field name in FST file
      integer,intent(in) :: F_datev               !Validity CMC datestamp format
      type(vgrid_descriptor),intent(out) :: F_vgrid  !Vgrid of record
      integer,pointer              :: F_iplist(:)  !List of IP1s for the field
      character(len=*),intent(out) :: F_levtype_S  !Short level name
      integer,pointer,optional     :: F_fstkeys(:) !List of field's fst keys
      integer,intent(out),optional :: F_sfcRefKey  !Key for Sfc Ref Field
      integer,intent(out),optional :: F_sfcRefKey2 !Key for Sfc LS Ref Field
      !@return
      integer :: F_istat
      !*@/
      integer :: ier,i,j,ni,nj,nk,nkeys,dateo,deet,npas,nbits,datyp,ip2,ip3, &
           ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,ex1,ex2,ex3,ikind,nip1,nip1b, &
           msglevel0,i2
      integer, dimension(MAXLEV) :: keylist,fst_iplist,subip1list

      integer, pointer :: vg_iplist(:)
      real :: ip1r
      character(len=1) :: typvar,grtyp
      character(len=4) :: nomvar
      character(len=12) :: etiket
      character(len=1), dimension(NLEVTYP) :: levtype_S=(/'M','T'/)
      logical :: sfc_L, ok_L
      !------------------------------------------------------------------
      F_istat = RMN_ERR
      F_levtype_S = ' '
      call vgrid_nullify(F_vgrid) 
      if (present(F_fstkeys)) then
         !TODO: reuse provided space if any and big enough
         if (associated(F_fstkeys)) deallocate(F_fstkeys,stat=ier)
         nullify(F_fstkeys)
      endif
      if (present(F_sfcRefKey))  F_sfcrefkey = RMN_ERR
      if (present(F_sfcRefKey2)) F_sfcrefkey2 = RMN_ERR

      ier = fstinl(F_unit,ni,nj,nk,F_datev,' ',RMN_ANY_I,RMN_ANY_I,RMN_ANY_I,' ',F_varname_S,keylist,nkeys,size(keylist))
      if (ier < 0 .or. nkeys < 1) then
         call msg(MSG_WARNING,'(vgrid_from_file) Cannot find any record for: '//trim(F_varname_S))
         return
      endif

      do i=1,nkeys
         ier = fstprm(keylist(i),dateo,deet,npas,ni,nj,nk,nbits,datyp,fst_iplist(i),ip2, &
              ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,ex1,ex2,ex3)
         if (ier /= 0) then
            call msg(MSG_WARNING,'(vgrid_from_file) Problem in fstprm')
            return
         endif
      enddo
      nip1 = sort(fst_iplist(1:nkeys),keylist(1:nkeys),unique_L=.true.)

      call msg_verbosity_get(msglevel0)
      call msg_verbosity(MYMSG_QUIET)
      ier = vgd_new(F_vgrid,unit=F_unit,ip1=ig1,ip2=ig2)
      ! if (ier /= VGD_OK) &
      !     ier = vgd_new(F_vgrid,unit=F_unit,ip1=ig4,ip2=RMN_ANY_I)
      if (ier /= VGD_OK) &
           ier = vgd_new(F_vgrid,unit=F_unit,ip1=RMN_ANY_I,ip2=RMN_ANY_I)
      call msg_verbosity(msglevel0)
      if (ier /= VGD_OK) then
         call msg(MSG_WARNING,'(vgrid_from_file) Problem in vgd_new')
         return
      endif

      ! Cross-check ip1 list to determine level type
      nullify(vg_iplist)
      j = 0
      DO_TYPE: do while (F_levtype_S == ' ' .and. j < NLEVTYP)
         j = j+1
         if (associated(vg_iplist)) deallocate(vg_iplist,stat=ier)
         nullify(vg_iplist)
         ier = vgd_get(F_vgrid,'VIP'//trim(levtype_S(j)),vg_iplist)
         if (ier /= VGD_OK) then
            call msg(MSG_WARNING,'(vgrid_from_file) Problem with vgd_get for VIP'//trim(levtype_S(j)))
            return
         endif
         if (size(vg_iplist) < nip1) cycle
         do i=1,nip1
            if (all(vg_iplist /= fst_iplist(i))) cycle DO_TYPE
         enddo

         nip1b = 0
         DO_I2: do i2=1,size(vg_iplist)
            do i=1,nip1
               if (vg_iplist(i2) == fst_iplist(i)) then
                  nip1b = nip1b + 1
                  subip1list(nip1b) = i
                  cycle DO_I2
               endif
            enddo
         enddo DO_I2

         F_levtype_S = trim(levtype_S(j))
      enddo DO_TYPE
      if (associated(vg_iplist)) deallocate(vg_iplist, stat=ier)

      ! Check for surface coordinates
      SFC_SEARCH: if (F_levtype_S == ' ') then
         sfc_L = .false.
         do i=1,nip1
            call convip_plus(fst_iplist(i),ip1r,ikind,-1,'',.false.)
            subip1list(i) = i
            if (ikind == IKIND_SURF) sfc_L = .true.
         enddo
         if (sfc_L) F_levtype_S = 'SFC'
      endif SFC_SEARCH

      ! Check for identified type and return IP1 list of found variable
      if (F_levtype_S == ' ') then
         call msg(MSG_WARNING,'(vgrid_from_file) unable to find full input dataset for '//trim(F_varname_S))
!!$         return
      endif

      ! Fill output fields
      ok_L = associated(F_iplist)
      if (ok_L) ok_L = (size(F_iplist) >= nip1)
      !#TODO: Memory leak if F_iplist already allocated but size too small
      if (.not.ok_L) allocate(F_iplist(nip1))
      if (present(F_fstkeys)) then
         ok_L = associated(F_fstkeys)
         if (ok_L) ok_L = (size(F_fstkeys) >= nip1)
         !#TODO: Memory leak if F_fstkeys already allocated but size too small
         if (.not.ok_L) allocate(F_fstkeys(nip1))
      endif
      do i=1,nip1
         F_iplist(i) = fst_iplist(subip1list(i))
         if (present(F_fstkeys)) F_fstkeys(i) = keylist(subip1list(i))
      enddo

      if (F_levtype_S /= 'SFC' .and. present(F_sfcRefKey)) then
         F_sfcRefKey = vgrid_from_file_rfld_key(F_unit, F_datev, F_vgrid)
         if (present(F_sfcRefKey2)) then
            F_sfcRefKey2 = vgrid_from_file_rfld_key(F_unit, F_datev, F_vgrid, &
                 F_rfls_L=.true.)
         endif
      endif

      F_istat = RMN_OK
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_4


   !/@*
   function vgrid_from_file_mpi_8(F_unit, F_varname_S, F_jdatev, F_vgrid, &
        F_iplist, F_levtype_S, F_fstkeys, F_sfcRefKey, F_sfcRefKey2, &
        F_comm_S, F_ipe_master, F_ipe) result(F_istat)
      implicit none
      !@objective Get vgrid and ip1 list for specific field in file
      !@author Ron McTaggart-Cowan, 2012-08
      !@revision
      ! 2014-11, S.Chamberland: sorting ip1,keys as in vgrid from files
      !@arguments
      integer,intent(in) :: F_unit                !Input file unit (opened)
      character(len=*),intent(in) :: F_varname_S  !Field name in FST file
      integer(INT64),intent(in) :: F_jdatev     !Valid date (jsec)
      type(vgrid_descriptor),intent(out) :: F_vgrid  !Vgrid of record
      integer,pointer              :: F_iplist(:)  !List of IP1s for the field
      character(len=*),intent(out) :: F_levtype_S  !Short level name
      integer,pointer,optional     :: F_fstkeys(:) !List of field's fst keys
      integer,intent(out),optional :: F_sfcRefKey  !Key for Sfc Ref Field
      integer,intent(out),optional :: F_sfcRefKey2 !Key for LSSfc Ref Field
       character(len=*), intent(in), optional :: F_comm_S  !- RPN_COMM communicator name
      integer, intent(in), optional :: F_ipe_master !- Sending PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      integer, intent(in), optional :: F_ipe      !- PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
     !@return
      integer :: F_istat
      !*@/
      integer :: datev, sfcrefkey2, me, ipe_master, istat
      character(len=64) :: comm_S
      !------------------------------------------------------------------
      comm_S = RPN_COMM_BLOC_COMM
      if (present(F_comm_S)) comm_S = F_comm_S
      if (present(F_ipe)) then
         me = F_ipe
      else
         call rpn_comm_rank(comm_S, me, istat)
      endif
      ipe_master = RPN_COMM_MASTER
      if (present(F_ipe_master)) ipe_master = F_ipe_master
      call vgrid_nullify(F_vgrid) 

      datev = RMN_ANY_DATE
      if (F_jdatev /= MU_JDATE_ANY) datev = jdate_to_cmc(F_jdatev)
      sfcRefKey2 = RMN_ERR
      if (present(F_fstkeys) .and. present(F_sfcRefKey)) then
         F_istat = vgrid_from_file_mpi_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, F_fstkeys, F_sfcRefKey, sfcRefKey2, &
              F_comm_S=comm_S,  F_ipe_master=ipe_master, F_ipe=me)
      else if (present(F_fstkeys)) then
         F_istat = vgrid_from_file_mpi_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, F_fstkeys, &
              F_comm_S=comm_S,  F_ipe_master=ipe_master, F_ipe=me)
      else if (present(F_sfcRefKey)) then
         F_istat = vgrid_from_file_mpi_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, F_sfcRefKey=F_sfcRefKey, &
              F_sfcRefKey2=sfcRefKey2, &
              F_comm_S=comm_S,  F_ipe_master=ipe_master, F_ipe=me)
      else
         F_istat = vgrid_from_file_mpi_4(F_unit, F_varname_S, datev, F_vgrid, &
              F_iplist, F_levtype_S, &
              F_comm_S=comm_S,  F_ipe_master=ipe_master, F_ipe=me)
      endif
      if (present(F_sfcRefKey2)) F_sfcRefKey2 = sfcRefKey2
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_mpi_8


   !/@*
   function vgrid_from_file_mpi_4(F_unit, F_varname_S, F_datev, F_vgrid, &
        F_iplist, F_levtype_S, F_fstkeys, F_sfcRefKey, F_sfcRefKey2, &
        F_comm_S, F_ipe_master, F_ipe) result(F_istat)
      implicit none
      !@objective Get vgrid and ip1 list for specific field in file
      !@author Stephane Chamberland, 2012-08
      !@arguments
      integer,intent(in) :: F_unit                !Input file unit (opened)
      character(len=*),intent(in) :: F_varname_S  !Field name in FST file
      integer,intent(in) :: F_datev               !Validity CMC datestamp format
      type(vgrid_descriptor),intent(out) :: F_vgrid  !Vgrid of record
      integer,pointer              :: F_iplist(:)  !List of IP1s for the field
      character(len=*),intent(out) :: F_levtype_S  !Short level name
      integer,pointer,optional     :: F_fstkeys(:) !List of field's fst keys
      integer,intent(out),optional :: F_sfcRefKey  !Key for Sfc Ref Field
      integer,intent(out),optional :: F_sfcRefKey2 !Key for Sfc LS Ref Field
      character(len=*), intent(in), optional :: F_comm_S  !- RPN_COMM communicator name
      integer, intent(in), optional :: F_ipe_master !- Sending PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      integer, intent(in), optional :: F_ipe      !- PE number in RPN_COMM communicator, default RPN_COMM_MASTER=0
      !@return
      integer :: F_istat
      !*@/
      integer :: levtype, istat, me, ipe_master
      character(len=16) :: sfcfld_S
      character(len=64) :: comm_S
      logical :: ismaster_L
      !------------------------------------------------------------------
      F_istat = RMN_OK
      if (present(F_fstkeys)) then
         !TODO: reuse provided space if any and big enough
         if (associated(F_fstkeys)) deallocate(F_fstkeys,stat=istat)
         nullify(F_fstkeys)
      endif
      if (present(F_sfcRefKey))  F_sfcrefkey  = RMN_ERR
      if (present(F_sfcRefKey2)) F_sfcrefkey2 = RMN_ERR
      call vgrid_nullify(F_vgrid) 

      comm_S = RPN_COMM_BLOC_COMM
      if (present(F_comm_S)) comm_S = F_comm_S
      if (present(F_ipe)) then
         me = F_ipe
      else
         call rpn_comm_rank(comm_S, me, istat)
      endif
      ipe_master = RPN_COMM_MASTER
      if (present(F_ipe_master)) ipe_master = F_ipe_master
      ismaster_L = (me == ipe_master)

      if (ismaster_L) then
         if (present(F_fstkeys)) then
            F_istat = vgrid_from_file(F_unit, F_varname_S, F_datev, F_vgrid, &
                 F_iplist, F_levtype_S, F_fstkeys)
         else
            F_istat = vgrid_from_file(F_unit, F_varname_S, F_datev, F_vgrid, &
                 F_iplist, F_levtype_S)
         endif
         levtype = VGRID_SURF_TYPE
         if (F_levtype_S(1:1) == 'M') levtype = VGRID_UPAIR_M_TYPE
         if (F_levtype_S(1:1) == 'T') levtype = VGRID_UPAIR_T_TYPE
         sfcfld_S = ' '
         if (F_levtype_S /= 'SFC' .and. present(F_sfcRefKey)) then
            F_sfcRefKey = vgrid_from_file_rfld_key(F_unit, F_datev, F_vgrid)
            if (present(F_sfcRefKey2)) then
               F_sfcRefKey2 = vgrid_from_file_rfld_key(F_unit, F_datev, &
                    F_vgrid, F_rfls_L=.true.)
            endif
         endif
      endif
      call collect_error(F_istat)
      if (.not.RMN_IS_OK(F_istat)) return
      F_istat = vgrid_wb_bcast(F_vgrid, F_iplist, levtype, sfcfld_S, &
           comm_S, ipe_master, me)
      !TODO: ?need to distribute F_sfcRefKey?
      F_levtype_S = 'SFC'
      if (levtype == VGRID_UPAIR_M_TYPE) F_levtype_S = 'M'
      if (levtype == VGRID_UPAIR_T_TYPE) F_levtype_S = 'T'
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_mpi_4


   !/@*
   function vgrid_from_file_rfld_key_mpi_8(F_unit, F_jdatev, F_vgrid, &
        F_rfls_L) result(F_sfcRefKey)
      implicit none
      !@objective Get Key for Sfc Ref Field
      !@author Stephane Chamberland, 2013-11
      !@arguments
      integer,intent(in) :: F_unit                   !Input file unit (already opened)
      integer(INT64),intent(in) :: F_jdatev        !Valid date (jsec)
      type(vgrid_descriptor),intent(in) :: F_vgrid   !Vertical grid descriptor of record
      logical,intent(in),optional :: F_rfls_L        !it .T. get Key for Sfc LS Ref Field
      !@return
      integer :: F_sfcRefKey
      !*@/
      integer :: datev
      !------------------------------------------------------------------
      datev = RMN_ANY_DATE
      if (F_jdatev /= MU_JDATE_ANY) datev = jdate_to_cmc(F_jdatev)
      if (present(F_rfls_L)) then
         F_sfcRefKey = vgrid_from_file_rfld_key_mpi_4(F_unit, datev, F_vgrid, F_rfls_L)
      else
         F_sfcRefKey = vgrid_from_file_rfld_key_mpi_4(F_unit, datev, F_vgrid)
      endif
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_rfld_key_mpi_8


   !/@*
   function vgrid_from_file_rfld_key_mpi_4(F_unit, F_datev, F_vgrid, &
        F_rfls_L) result(F_sfcRefKey)
      implicit none
      !@objective Get Key for Sfc Ref Field
      !@author Stephane Chamberland, 2013-11
      !@arguments
      integer,intent(in) :: F_unit                   !Input file unit (already opened)
      integer,intent(in) :: F_datev                  !Validity date (CMC datestamp format)
      type(vgrid_descriptor),intent(in) :: F_vgrid   !Vertical grid descriptor of record
      logical,intent(in),optional :: F_rfls_L        !it .T. get Key for Sfc LS Ref Field
      !@return
      integer :: F_sfcRefKey
      !*@/
      integer :: me,istat
      !------------------------------------------------------------------
      call rpn_comm_rank(RPN_COMM_BLOC_COMM,me,istat)
      F_sfcrefkey = RMN_OK
      if (me == RPN_COMM_MASTER) then
         if (present(F_rfls_L)) then
            F_sfcRefKey = vgrid_from_file_rfld_key(F_unit, F_datev, F_vgrid, F_rfls_L)
         else
            F_sfcRefKey = vgrid_from_file_rfld_key(F_unit, F_datev, F_vgrid)
         endif
      endif
      istat = F_sfcRefKey
      call collect_error(istat)
      if (.not.RMN_IS_OK(istat)) F_sfcRefKey = istat
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_rfld_key_mpi_4


   !/@*
   function vgrid_from_file_rfld_key_8(F_unit, F_jdatev, F_vgrid, F_rfls_L) &
        result(F_sfcRefKey)
      implicit none
      !@objective Get Key for Sfc Ref Field
      !@author Stephane Chamberland, 2013-11
      !@arguments
      integer,intent(in) :: F_unit                   !Input file unit (already opened)
      integer(INT64),intent(in) :: F_jdatev        !Valid date (jsec)
      type(vgrid_descriptor),intent(in) :: F_vgrid   !Vertical grid descriptor of record
      logical,intent(in),optional :: F_rfls_L        !it .T. get Key for Sfc LS Ref Field
      !@return
      integer :: F_sfcRefKey
      !*@/
      integer :: datev
      !------------------------------------------------------------------
      datev = RMN_ANY_DATE
      if (F_jdatev /= MU_JDATE_ANY) datev = jdate_to_cmc(F_jdatev)
      if (present(F_rfls_L)) then
         F_sfcRefKey = vgrid_from_file_rfld_key_4(F_unit, datev, F_vgrid, F_rfls_L)
      else
         F_sfcRefKey = vgrid_from_file_rfld_key_4(F_unit, datev, F_vgrid)
      endif
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_rfld_key_8


   !/@*
   function vgrid_from_file_rfld_key_4(F_unit, F_datev, F_vgrid, F_rfls_L) &
        result(F_sfcRefKey)
      implicit none
      !@objective Get Key for Sfc Ref Field
      !@author Stephane Chamberland, 2013-11
      !@arguments
      integer,intent(in) :: F_unit                   !Input file unit (already opened)
      integer,intent(in) :: F_datev                  !Validity date (CMC datestamp format)
      type(vgrid_descriptor),intent(in) :: F_vgrid   !Vertical grid descriptor of record
      logical,intent(in),optional :: F_rfls_L        !it .T. get Key for Sfc LS Ref Field
      !@return
      integer :: F_sfcRefKey
      !*@/
      integer :: ier,i,ni,nj,nk,nkeys,dateo,deet,npas,nbits,datyp,ip1,ip2,ip3, &
           ig1,ig2,ig3,ig4,swa,lng,dltf,ubc,ex1,ex2,ex3,istat,ikind,ivers
      integer, dimension(MAXLEV) :: keylist
      character(len=1) :: typvar,grtyp
      character(len=4) :: nomvar,sfcfld_S,rfld_s
      character(len=12) :: etiket
      !------------------------------------------------------------------
      F_sfcRefKey = RMN_ERR
      rfld_s = 'RFLD'
      if (present(F_rfls_L)) then
         if (F_rfls_L) rfld_s = 'RFLS'
      endif
      istat = vgd_get(F_vgrid, key=rfld_s, value=sfcfld_S)
      if (istat /= VGD_OK) then
         istat = vgd_get(F_vgrid, key='KIND', value=ikind)
         istat = vgd_get(F_vgrid, key='VERS', value=ivers)
         if ((rfld_s == 'RFLD' .and. ikind == 2 .and. ivers == 1) .or. &
              (rfld_s == 'RFLS' .and. ikind /= 5 .and. ivers /= 100)) then
            F_sfcRefKey = VGRID_FROM_FILE_NORFLD
         else
            call msg(MSG_WARNING,'(vgrid_from_file) Problem with vgd_get for '//trim(rfld_s))
         endif
         return
      endif
      ier = fstinl(F_unit, ni, nj, nk, F_datev, '', RMN_ANY_I, RMN_ANY_I, &
           RMN_ANY_I, '', sfcfld_S, keylist, nkeys, size(keylist))
      if (rfld_s == 'RFLS' .and. (ier < 0 .or. nkeys < 1)) then
         ier = fstinl(F_unit, ni, nj, nk, RMN_ANY_I, '', RMN_ANY_I ,RMN_ANY_I, &
              RMN_ANY_I, '', sfcfld_S, keylist, nkeys, size(keylist))
      endif
      if (ier < 0 .or. nkeys < 1) then
         call msg(MSG_WARNING,'(vgrid_from_file) Cannot find any record for '//trim(rfld_s)//': '//trim(sfcfld_S))
         return
      endif
      do i=1,nkeys
         ier = fstprm(keylist(i),dateo,deet,npas,ni,nj,nk, &
              nbits,datyp,ip1,ip2, &
              ip3,typvar,nomvar,etiket,grtyp,ig1,ig2,ig3,ig4, &
              swa,lng,dltf,ubc,ex1,ex2,ex3)
         if (ier /= 0) then
            call msg(MSG_WARNING,'(vgrid_from_file) Problem in fstprm')
            cycle
         endif
         if (all(typvar /= (/'r','R'/))) then
            F_sfcRefKey = keylist(i)
            exit
         endif
      enddo
      !------------------------------------------------------------------
      return
   end function vgrid_from_file_rfld_key_4

end module vgrid_from_file_mod
