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
module statfld_mod
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   !@author Stephane Chamberland, 2012-05
   !@Objective calcule la moyenne, la variance, le minimum et 
   !           le maximum d un champs et imprime le resultat
   private
   public :: statfld,statfld_print
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   interface statfld
      module procedure statfld_r4_2d
      module procedure statfld_r4_2d_e
      module procedure statfld_r4_3d
      module procedure statfld_r4_3d_e
      module procedure statfld_r8_2d
      module procedure statfld_r8_2d_e
      module procedure statfld_r8_3d
      module procedure statfld_r8_3d_e
   end interface

   character(len=256),parameter :: STATFORMAT_R4 = &
        "(i4,1x,a16,' Mean:',1pe14.6,' Var:',1pe14.6, "//&
        "' Min:[(',i4,',',i4,',',i3,')',1pe14.6,']',"//&
        "' Max:[(',i4,',',i4,',',i3,')',1pe14.6,'] ',a)"
   character(len=256),parameter :: STATFORMAT_R8 = &
        "(i4,1x,a16,' Mean:',1pe22.12,' Var:',1pe22.12, "//&
        "' Min:[(',i4,',',i4,',',i3,')',1pe22.12,']',"//&
        "' Max:[(',i4,',',i4,',',i3,')',1pe22.12,'] ',a)"

   !TODO: integer types

contains


   !/@*
   subroutine statfld_print(F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax,F_nv_S,F_no,F_from_S,F_rx) 
      implicit none
      !@arguments
      real(REAL64),intent(in) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(in) :: F_ijkmin(3),F_ijkmax(3),F_no,F_rx
      character(len=*),intent(in) :: F_nv_S, F_from_S
      !*@/
      integer,external :: msg_getUnit
      character(len=256) :: msg_S
      integer :: msgUnit,n
      !---------------------------------------------------------------
      msgUnit = msg_getUnit(MSG_INFO)
      if (RMN_IS_OK(msgUnit)) then
         n=min(max(1,len_trim(F_from_S)),32)
         if (F_rx < 8) then
            write(msg_S,STATFORMAT_R4) F_no,F_nv_S,F_mean,F_var,&
                 F_ijkmin(1),F_ijkmin(2),F_ijkmin(3),F_rmin,&
                 F_ijkmax(1),F_ijkmax(2),F_ijkmax(3),F_rmax,F_from_S(1:n)
         else
            write(msg_S,STATFORMAT_R8) F_no,F_nv_S,F_mean,F_var,&
                 F_ijkmin(1),F_ijkmin(2),F_ijkmin(3),F_rmin,&
                 F_ijkmax(1),F_ijkmax(2),F_ijkmax(3),F_rmax,F_from_S(1:n)
         endif
         call msg(MSG_INFO,msg_S)
      endif
      !---------------------------------------------------------------
      return
   end subroutine statfld_print


   !/@*
   subroutine statfld_r4_2d_e(F_fld,F_nv_S,F_no,F_from_S,F_rx,F_ijk0) 
      implicit none
      !@arguments
      ! F_fld         I  Field to be operated on (or array section of it)
      ! F_nv_S        I  User provided string to define F_field
      ! F_no          I  Usually the timestep # or level (if k0==kn)
      ! F_from_S      I  Usually the name of the calling subroutine
      ! F_rx          I  Stat precision 4/8
      ! F_ijk0        I  For partial F_fld, start index of the field crop
      real,intent(in) :: F_fld(:,:)
      character(len=*),intent(in) :: F_nv_S, F_from_S
      integer,intent(in) :: F_no,F_rx
      integer,intent(in),optional :: F_ijk0(3)
      !*@/
      integer :: ijk0(3),ijkmin(3),ijkmax(3)
      real(REAL64) :: mean,var,rmin,rmax
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      call statfld(F_fld,mean,var,rmin,rmax,ijkmin,ijkmax,ijk0) 
      call statfld_print(mean,var,rmin,rmax,ijkmin,ijkmax,F_nv_S,F_no,F_from_S,F_rx)
      !---------------------------------------------------------------
      return
   end subroutine statfld_r4_2d_e


   !/@*
   subroutine statfld_r4_2d(F_fld,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax,F_ijk0) 
      implicit none
      !@arguments
      real,intent(in) :: F_fld(:,:)
      real(REAL64),intent(out) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(out) :: F_ijkmin(3),F_ijkmax(3)
      integer,intent(in),optional :: F_ijk0(3) !for partial F_fld, start index of the field crop
      !*@/
      integer  :: i,j,ijk0(3),lijk(2),uijk(2)
      real(REAL64) :: mysum,rnijk
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      ijk0 = ijk0 - 1
      rnijk = dble(size(F_fld))
      lijk = lbound(F_fld)
      uijk = ubound(F_fld)

      mysum = 0.
      do j=lijk(2),uijk(2)
         do i=lijk(1),uijk(1)
            mysum = mysum + dble(F_fld(i,j))
         enddo
      enddo

      F_mean = mysum / rnijk

      mysum = 0.
      do j=lijk(2),uijk(2)
         do i=lijk(1),uijk(1)
            mysum = mysum + (dble(F_fld(i,j)) - F_mean)**2.
         enddo
      enddo

      F_var = sqrt(mysum / rnijk)

      F_ijkmin(1:2) = minloc(F_fld) ; F_ijkmin(3) = 1
      F_ijkmax(1:2) = maxloc(F_fld) ; F_ijkmax(3) = 1
      F_rmin = dble(F_fld(F_ijkmin(1),F_ijkmin(2)))
      F_rmax = dble(F_fld(F_ijkmax(1),F_ijkmax(2)))
      F_ijkmin = F_ijkmin + ijk0
      F_ijkmax = F_ijkmax + ijk0
      !---------------------------------------------------------------
      return
   end subroutine statfld_r4_2d


   !/@*
   subroutine statfld_r4_3d_e(F_fld,F_nv_S,F_no,F_from_S,F_rx,F_ijk0) 
      implicit none
      !@arguments
      ! F_fld         I  Field to be operated on (or array section of it)
      ! F_nv_S        I  User provided string to define F_field
      ! F_no          I  Usually the timestep # or level (if k0==kn)
      ! F_from_S      I  Usually the name of the calling subroutine
      ! F_rx          I  Stat precision 4/8
      ! F_ijk0        I  For partial F_fld, start index of the field crop
      real,intent(in) :: F_fld(:,:,:)
      character(len=*),intent(in) :: F_nv_S, F_from_S
      integer,intent(in) :: F_no,F_rx
      integer,intent(in),optional :: F_ijk0(3)
      !*@/
      integer :: ijk0(3),ijkmin(3),ijkmax(3)
      real(REAL64) :: mean,var,rmin,rmax
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      call statfld(F_fld,mean,var,rmin,rmax,ijkmin,ijkmax,ijk0) 
      call statfld_print(mean,var,rmin,rmax,ijkmin,ijkmax,F_nv_S,F_no,F_from_S,F_rx)
      !---------------------------------------------------------------
      return
   end subroutine statfld_r4_3d_e


   !/@*
   subroutine statfld_r4_3d(F_fld,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax,F_ijk0) 
      implicit none
      !@arguments
      real,intent(in) :: F_fld(:,:,:)
      real(REAL64),intent(out) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(out) :: F_ijkmin(3),F_ijkmax(3)
      integer,intent(in),optional :: F_ijk0(3) !for partial F_fld, start index of the field crop
      !*@/
      integer  :: i,j,k,ijk0(3),lijk(3),uijk(3)
      real(REAL64) :: mysum,rnijk
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      ijk0 = ijk0 - 1
      rnijk = dble(size(F_fld))
      lijk = lbound(F_fld)
      uijk = ubound(F_fld)

      mysum = 0.
      do k=lijk(3),uijk(3)
         do j=lijk(2),uijk(2)
            do i=lijk(1),uijk(1)
               mysum = mysum + dble(F_fld(i,j,k))
            enddo
         enddo
      enddo

      F_mean = mysum / rnijk

      mysum = 0.
      do k=lijk(3),uijk(3)
         do j=lijk(2),uijk(2)
            do i=lijk(1),uijk(1)
               mysum = mysum + (dble(F_fld(i,j,k)) - F_mean)**2.
            enddo
         enddo
      enddo

      F_var = sqrt(mysum / rnijk)

      F_ijkmin = minloc(F_fld)
      F_ijkmax = maxloc(F_fld)
      F_rmin = dble(F_fld(F_ijkmin(1),F_ijkmin(2),F_ijkmin(3)))
      F_rmax = dble(F_fld(F_ijkmax(1),F_ijkmax(2),F_ijkmax(3)))
      F_ijkmin = F_ijkmin + ijk0
      F_ijkmax = F_ijkmax + ijk0
      !---------------------------------------------------------------
      return
   end subroutine statfld_r4_3d


   !/@*
   subroutine statfld_r8_2d_e(F_fld,F_nv_S,F_no,F_from_S,F_rx,F_ijk0) 
      implicit none
      !@arguments
      ! F_fld         I  Field to be operated on (or array section of it)
      ! F_nv_S        I  User provided string to define F_field
      ! F_no          I  Usually the timestep # or level (if k0==kn)
      ! F_from_S      I  Usually the name of the calling subroutine
      ! F_rx          I  Stat precision 4/8
      ! F_ijk0        I  For partial F_fld, start index of the field crop
      real(REAL64),intent(in) :: F_fld(:,:)
      character(len=*),intent(in) :: F_nv_S, F_from_S
      integer,intent(in) :: F_no,F_rx
      integer,intent(in),optional :: F_ijk0(3)
      !*@/
      integer :: ijk0(3),ijkmin(3),ijkmax(3)
      real(REAL64) :: mean,var,rmin,rmax
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      call statfld(F_fld,mean,var,rmin,rmax,ijkmin,ijkmax,ijk0) 
      call statfld_print(mean,var,rmin,rmax,ijkmin,ijkmax,F_nv_S,F_no,F_from_S,F_rx)
      !---------------------------------------------------------------
      return
   end subroutine statfld_r8_2d_e


   !/@*
   subroutine statfld_r8_2d(F_fld,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax,F_ijk0) 
      implicit none
      !@arguments
      real(REAL64),intent(in) :: F_fld(:,:)
      real(REAL64),intent(out) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(out) :: F_ijkmin(3),F_ijkmax(3)
      integer,intent(in),optional :: F_ijk0(3) !for partial F_fld, start index of the field crop
      !*@/
      integer  :: i,j,ijk0(3),lijk(2),uijk(2)
      real(REAL64) :: mysum,rnijk
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      ijk0 = ijk0 - 1
      rnijk = dble(size(F_fld))
      lijk = lbound(F_fld)
      uijk = ubound(F_fld)

      F_mean = sum(F_fld) / rnijk
      
      mysum = 0.
      do j=lijk(2),uijk(2)
         do i=lijk(1),uijk(1)
            mysum = mysum + (F_fld(i,j) - F_mean)**2.
         enddo
      enddo

      F_var = sqrt(mysum / rnijk)

      F_ijkmin(1:2) = minloc(F_fld) ; F_ijkmin(3) = 1
      F_ijkmax(1:2) = maxloc(F_fld) ; F_ijkmax(3) = 1
      F_rmin = F_fld(F_ijkmin(1),F_ijkmin(2))
      F_rmax = F_fld(F_ijkmax(1),F_ijkmax(2))
      F_ijkmin = F_ijkmin + ijk0
      F_ijkmax = F_ijkmax + ijk0
      !---------------------------------------------------------------
      return
   end subroutine statfld_r8_2d


   !/@*
   subroutine statfld_r8_3d_e(F_fld,F_nv_S,F_no,F_from_S,F_rx,F_ijk0) 
      implicit none
      !@arguments
      ! F_fld         I  Field to be operated on (or array section of it)
      ! F_nv_S        I  User provided string to define F_field
      ! F_no          I  Usually the timestep # or level (if k0==kn)
      ! F_from_S      I  Usually the name of the calling subroutine
      ! F_rx          I  Stat precision 4/8
      ! F_ijk0        I  For partial F_fld, start index of the field crop
      real(REAL64),intent(in) :: F_fld(:,:,:)
      character(len=*),intent(in) :: F_nv_S, F_from_S
      integer,intent(in) :: F_no,F_rx
      integer,intent(in),optional :: F_ijk0(3)
      !*@/
      integer :: ijk0(3),ijkmin(3),ijkmax(3)
      real(REAL64) :: mean,var,rmin,rmax
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      call statfld(F_fld,mean,var,rmin,rmax,ijkmin,ijkmax,ijk0) 
      call statfld_print(mean,var,rmin,rmax,ijkmin,ijkmax,F_nv_S,F_no,F_from_S,F_rx)
      !---------------------------------------------------------------
      return
   end subroutine statfld_r8_3d_e


   !/@*
   subroutine statfld_r8_3d(F_fld,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax,F_ijk0) 
      implicit none
      !@arguments
      real(REAL64),intent(in) :: F_fld(:,:,:)
      real(REAL64),intent(out) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(out) :: F_ijkmin(3),F_ijkmax(3)
      integer,intent(in),optional :: F_ijk0(3) !for partial F_fld, start index of the field crop
      !*@/
      integer  :: i,j,k,ijk0(3),lijk(3),uijk(3)
      real(REAL64) :: mysum,rnijk
      !---------------------------------------------------------------
      ijk0 = (/1,1,1/)
      if (present(F_ijk0)) ijk0 = F_ijk0
      ijk0 = ijk0 - 1
      rnijk = dble(size(F_fld))
      lijk = lbound(F_fld)
      uijk = ubound(F_fld)

      F_mean = sum(F_fld) / rnijk

      mysum = 0.
      do k=lijk(3),uijk(3)
         do j=lijk(2),uijk(2)
            do i=lijk(1),uijk(1)
               mysum = mysum + (F_fld(i,j,k) - F_mean)**2.
            enddo
         enddo
      enddo

      F_var = sqrt(mysum / rnijk)

      F_ijkmin = minloc(F_fld)
      F_ijkmax = maxloc(F_fld)
      F_rmin = F_fld(F_ijkmin(1),F_ijkmin(2),F_ijkmin(3))
      F_rmax = F_fld(F_ijkmax(1),F_ijkmax(2),F_ijkmax(3))
      F_ijkmin = F_ijkmin + ijk0
      F_ijkmax = F_ijkmax + ijk0
      !---------------------------------------------------------------
      return
   end subroutine statfld_r8_3d


end module statfld_mod
