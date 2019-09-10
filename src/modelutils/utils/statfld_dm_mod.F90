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
module statfld_dm_mod
   use, intrinsic :: iso_fortran_env, only: REAL64
   use iso_c_binding
   use rpn_comm_itf_mod
   use ptopo_utils
   use statfld_mod
   implicit none
   !@author Stephane Chamberland, 2012-05
   !@Objective calcule la moyenne, la variance, le minimum et 
   !           le maximum d un champs et imprime le resultat avec GLBCOLC
   !           de sorte que le calcule ne depend pas de la topologie
   private
   public :: statfld_dm,statfld_dm_print
   !*@/
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <msg.h>

   interface statfld_dm_print
      module procedure statfld_print
   end interface

   interface statfld_dm
      module procedure statfld_dm_r4_2d
      module procedure statfld_dm_r4_2d_e
      module procedure statfld_dm_r4_3d
      module procedure statfld_dm_r4_3d_e
   end interface

contains

   !/@*
   subroutine statfld_dm_r4_2d_e(F_fld,F_nv_S,F_no,F_from_S,F_rx) 
      implicit none
      !@arguments
      ! F_fld         I  Field to be operated on (w/o halos)
      ! F_nv_S        I  User provided string to define F_field
      ! F_no          I  Usually the timestep # or level (if k0==kn)
      ! F_from_S      I  Usually the name of the calling subroutine
      ! F_rx          I  Stat precision 4/8
      real,intent(in) :: F_fld(:,:)
      character(len=*),intent(in) :: F_nv_S, F_from_S
      integer,intent(in) :: F_no,F_rx
      !*@/
      integer :: ijkmin(3),ijkmax(3)
      real(REAL64) :: mean,var,rmin,rmax
      !---------------------------------------------------------------
      call statfld_dm(F_fld,mean,var,rmin,rmax,ijkmin,ijkmax) 
      call statfld_print(mean,var,rmin,rmax,ijkmin,ijkmax,F_nv_S,F_no,F_from_S,F_rx)
      !---------------------------------------------------------------
      return
   end subroutine statfld_dm_r4_2d_e


   !/@*
   subroutine statfld_dm_r4_2d(F_fld,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax) 
      implicit none
      !@arguments
      real,intent(in),target :: F_fld(:,:)
      real(REAL64),intent(out) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(out) :: F_ijkmin(3),F_ijkmax(3)
      !*@/
      real,pointer :: wk(:,:),pfld(:,:)
      integer :: istat,lni,lnj,gni,gnj,gi0,gj0
      !---------------------------------------------------------------
      call ptopo_init_var()
      F_mean = 0.
      F_var = 0.
      F_rmin = 0.
      F_rmax = 0.
      F_ijkmin = 0
      F_ijkmax = 0
 
      lni = size(F_fld,1)
      lnj = size(F_fld,2)
      istat = ptopo_collect_dims(RPN_COMM_GRID,lni,lnj,gni,gnj,gi0,gj0)
      nullify(wk)
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         allocate(wk(gni,gnj),stat=istat)
      endif
      pfld => F_fld
      istat = ptopo_collect(wk,pfld,RPN_COMM_GRID,gi0,gj0,lni,lnj)
 
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         call statfld(wk,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax) 
!!$      else !#TODO: bcast
      endif
      if (associated(wk)) deallocate(wk,stat=istat)
      !---------------------------------------------------------------
      return
   end subroutine statfld_dm_r4_2d


   !/@*
   subroutine statfld_dm_r4_3d_e(F_fld,F_nv_S,F_no,F_from_S,F_rx) 
      implicit none
      !@arguments
      ! F_fld         I  Field to be operated on (w/o halos)
      ! F_nv_S        I  User provided string to define F_field
      ! F_no          I  Usually the timestep # or level (if k0==kn)
      ! F_from_S      I  Usually the name of the calling subroutine
      ! F_rx          I  Stat precision 4/8
      real,intent(in) :: F_fld(:,:,:)
      character(len=*),intent(in) :: F_nv_S, F_from_S
      integer,intent(in) :: F_no,F_rx
      !*@/
      integer :: ijkmin(3),ijkmax(3)
      real(REAL64) :: mean,var,rmin,rmax
      !---------------------------------------------------------------
      call statfld_dm(F_fld,mean,var,rmin,rmax,ijkmin,ijkmax) 
      call statfld_print(mean,var,rmin,rmax,ijkmin,ijkmax,F_nv_S,F_no,F_from_S,F_rx)
      !---------------------------------------------------------------
      return
   end subroutine statfld_dm_r4_3d_e


   !/@*
   subroutine statfld_dm_r4_3d(F_fld,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax) 
      implicit none
      !@arguments
      real,intent(in),target :: F_fld(:,:,:)
      real(REAL64),intent(out) :: F_mean,F_var,F_rmin,F_rmax
      integer,intent(out) :: F_ijkmin(3),F_ijkmax(3)
      !*@/
      real,pointer :: wk(:,:,:),pfld(:,:,:)
      integer :: istat,lni,lnj,lnk,gni,gnj,gi0,gj0
      !---------------------------------------------------------------
      call ptopo_init_var()
      F_mean = 0.
      F_var = 0.
      F_rmin = 0.
      F_rmax = 0.
      F_ijkmin = 0
      F_ijkmax = 0
 
      lni = size(F_fld,1)
      lnj = size(F_fld,2)
      lnk = size(F_fld,3)
      istat = ptopo_collect_dims(RPN_COMM_GRID,lni,lnj,gni,gnj,gi0,gj0)
      nullify(wk)
      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         allocate(wk(gni,gnj,lnk),stat=istat)
      endif
      pfld => F_fld
      istat = ptopo_collect(wk,pfld,RPN_COMM_GRID,gi0,gj0,lni,lnj)

      if (ptopo_grid_ipe == RPN_COMM_MASTER) then
         call statfld(wk,F_mean,F_var,F_rmin,F_rmax,F_ijkmin,F_ijkmax) 
!!$      else !#TODO: bcast
      endif
      if (associated(wk)) deallocate(wk,stat=istat)
      !---------------------------------------------------------------
      return
   end subroutine statfld_dm_r4_3d


end module statfld_dm_mod
