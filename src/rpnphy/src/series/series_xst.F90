!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

module series_xst_mod
   use clib_itf_mod, only: clib_toupper
   use series_options
   use phygridmap, only: phydim_ni, phydim_nj, phydim_nk
   implicit none
   private

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   public :: series_xst_geo, series_xst, series_isstep, series_isvar

   interface series_xst
      module procedure series_xst2d
      module procedure series_xst3d
   end interface series_xst
 
   !@Author R. Benoit (RPN 1984)
   !@Revision
   ! 001      J. Cote RPN(January 1985)
   !                - Recoding compatible SEF/RFE version
   !                - Documentation
   ! 002      M. Lepine  -  RFE model code revision project (Feb 87)
   ! 002      M. Lepine  -  Ensuring that the code is re-entrant (Oct 87)
   ! 003      R. Benoit  -  Extraction by levels for the PROFILS
   ! 004      B. Reid  (June 89) - Zonal diagnostics
   ! 005      B. Bilodeau (Mar 91) - Eliminate the entrance point
   !                VSER XST and the call to ZONXST.
   ! 006      B. Bilodeau  (July 1991)- Adaptation to UNIX
   ! 007      N. Ek (Mar 1995) - output only every SERINT time-steps
   ! 008      B. Bilodeau (Nov 1995) - KAM
   ! 009      B. Bilodeau and M. Desgagne (March 2001) - Build lists
   !            surface(m,2) and profils(m,2) even if nstat=0 because in
   !            MPI mode, processor 0 needs full lists in call to serwrit2
   ! 010      B. Bilodeau (Jan 2006) - Variable NOM converted in upper case
   ! 2017-01, S.Chamberland - Major revision

contains

   !/@*
   function series_isstep(F_varname_S) result(F_istat_L)
      implicit none
      !@objective 
      !@arguments
      character(len=*), intent(in), optional :: F_varname_S
      !@return
      logical :: F_istat_L
      !*@/
      !---------------------------------------------------------------
      F_istat_L = .false.
      if (series_paused_L .or. .not.(series_initok_L .and. series_on_L)) return
!!$      if (series_kount /= 1 .and. (mod(series_kount, p_serg_srwri) /= 0)) return
      if (mod(series_kount, p_serg_srwri) /= 0) return
      F_istat_L = .true.
      if (present(F_varname_S)) F_istat_L = series_isvar(F_varname_S)
      !---------------------------------------------------------------
      return
   end function series_isstep


   !/@*
   function series_isvar(F_varname_S) result(F_istat_L)
      implicit none
      !@objective 
      !@arguments
      character(len=*), intent(in) :: F_varname_S
      !@return
      logical :: F_istat_L
      !*@/
      character(len=SER_STRLEN) :: varname_S
      integer :: istat
      !---------------------------------------------------------------
      varname_S = F_varname_S
      istat = clib_toupper(varname_S)
      F_istat_L = .false.
      if (any(p_serg_srsrf_s(1:series_nsurf) == varname_S) .or. &
           any(P_serg_srprf_s(1:series_nprof) == varname_S)) &
           F_istat_L = .true.
      !---------------------------------------------------------------
      return
   end function series_isvar


   !/@*
   subroutine series_xst_geo(F_data, F_varname_S, F_j)
      implicit none
      !@objective 
      !@arguments
      real,    intent(in) :: F_data(:)  !# Slice J data(1:phydim_ni)
      character(len=*), intent(in) :: F_varname_S
      integer, intent(in) :: F_j        !# Physics slice number
      !*@/
      integer :: nn, nnl, i, idxvar
      !---------------------------------------------------------------
      if (series_nstnl == 0) return

      idxvar = priv_varidx(F_varname_S, p_serg_srgeo_s(1:series_ngeo))
      if (RMN_IS_OK(idxvar)) then
         do nnl=1,series_nstnl
            nn = series_stnl_idxb(nnl)
            if (nn > 0 .and. nn <= series_nstnb .and. &
                 series_stnb(nn)%jphy == F_j) then
               !#TODO: review index order
               i = series_stnb(nn)%iphy
               series_gdata(idxvar, nn) = F_data(i)
               series_gdone(idxvar) = 1
            endif
         enddo
      endif
      !---------------------------------------------------------------
      return
   end subroutine series_xst_geo


   !/@*
   subroutine series_xst2d(F_data, F_varname_S, F_j, F_facs, F_facf, &
        F_overwrite_L, F_warnprof_L)
      implicit none
      !@objective 
      !@arguments
      real, intent(in) :: F_data(:)  !# Slice J data(1:phydim_ni)
      character(len=*), intent(in) :: F_varname_S
      integer, intent(in) :: F_j  !# Physics slice number
      real, intent(in), optional :: F_facs  !# mul factor on the time-series
      real, intent(in), optional :: F_facf  !# mul factor before extraction
      logical, intent(in), optional :: F_overwrite_L  !# Allow overrite
      logical, intent(in), optional :: F_warnprof_L   !# Warn if profile
      !*@/
      real :: facs, facf
      logical :: overwrite_L, warnprof_L
      integer :: nn, nnl, i, idxvar, kount2
      !---------------------------------------------------------------
      if (series_nstnl == 0) return
      if (.not.series_isstep(F_varname_S)) return

      facs = 0.
      if (present(F_facs)) facs = F_facs
      facf = 1.
      if (present(F_facf)) facf = F_facf
      overwrite_L = .true.
      if (present(F_overwrite_L)) overwrite_L = F_overwrite_L
      warnprof_L = .true.
      if (present(F_warnprof_L)) warnprof_L = F_warnprof_L

!!$      kount2 = 1 + mod(series_kount, series_out_nsteps)  !#TODO:
      kount2 = 1

      idxvar = priv_varidx(F_varname_S, p_serg_srsrf_s(1:series_nsurf))
      IF_SURF: if (RMN_IS_OK(idxvar)) then
         do nnl=1,series_nstnl
            nn = series_stnl_idxb(nnl)
            if (nn > 0 .and. nn <= series_nstnb .and. &
                 series_stnb(nn)%jphy == F_j) then
               if (series_sdone2(idxvar,nnl) == series_kount .and. &
                    .not.overwrite_L) then
                  if (series_kount == P_serg_srwri) &
                       call msg(MSG_INFOPLUS, PKGNAME_S//'skipped, called twice in same step for var: '//trim(F_varname_S))
                  return
               endif
               !#TODO: review index order
               i = series_stnb(nn)%iphy
               series_sdata(idxvar, nn, kount2) = &
                    facs * series_sdata(idxvar, nn, kount2) + &
                    facf * F_data(i)
               series_sdone(idxvar) = 1
               series_sdone2(idxvar,nnl) = series_kount
            endif
         enddo
      endif IF_SURF

      !# Profiles
      if (warnprof_L) then
         idxvar = priv_varidx(F_varname_S, P_serg_srprf_s(1:series_nprof))
         if (RMN_IS_OK(idxvar) .and. series_kount <= P_serg_srwri .and. F_j == 1) then
            call msg(MSG_WARNING, PKGNAME_S//'series_xst, Profile requested but only 2d data is provided for '//trim(F_varname_S)//' - Ignoring')
         endif
      endif
      !---------------------------------------------------------------
      return
   end subroutine series_xst2d


   !/@*
   subroutine series_xst3d(F_data, F_varname_S, F_j, F_facs, F_facf, &
        F_overwrite_L)
      implicit none
      !@objective 
      !@arguments
      real, intent(in) :: F_data(:,:)  !# Slice j data(1:phydim_ni, 1:phydim_nk)
      character(len=*), intent(in) :: F_varname_S
      integer, intent(in) :: F_j  !# Physics slice number
      real, intent(in), optional :: F_facs  !# mul factor on the time-series
      real, intent(in), optional :: F_facf  !# mul factor before extraction
      logical, intent(in), optional :: F_overwrite_L  !# Allow overrite
      !*@/
      real :: facs, facf
      logical :: overwrite_L
      integer :: nn, nnl, idxvar, kount2, i, k, nk
      !---------------------------------------------------------------
      if (series_nstnl == 0) return
      if (.not.series_isstep(F_varname_S)) return
      facs = 0.
      if (present(F_facs)) facs = F_facs
      facf = 1.
      if (present(F_facf)) facf = F_facf
      overwrite_L = .true.
      if (present(F_overwrite_L)) overwrite_L = F_overwrite_L

!!$      kount2 = 1 + mod(series_kount, series_out_nsteps)  !#TODO:
      kount2 = 1

      !# Surface
      idxvar = priv_varidx(F_varname_S, p_serg_srsrf_s(1:series_nsurf))
      IF_SURF: if (RMN_IS_OK(idxvar)) then
         k = min(phydim_nk, ubound(F_data,2))
         call series_xst2d(F_data(:,k), F_varname_S, F_j, facs, facf, &
              overwrite_L, F_warnprof_L=.false.)
      endif IF_SURF

      !# Profiles
      idxvar = priv_varidx(F_varname_S, P_serg_srprf_s(1:series_nprof))
      IF_PROF: if (RMN_IS_OK(idxvar)) then
         do nnl=1,series_nstnl
            nn = series_stnl_idxb(nnl)
            if (nn > 0 .and. nn <= series_nstnb .and. &
                 series_stnb(nn)%jphy == F_j) then
               if (series_pdone2(idxvar,nnl) == 1 .and. &
                    .not.overwrite_L) then
                  if (series_kount == P_serg_srwri) &
                       call msg(MSG_INFOPLUS, PKGNAME_S//'skipped, called twice in same step for var: '//trim(F_varname_S))
                  return
               endif
               !#TODO: review index order
               i = series_stnb(nn)%iphy
               nk = min(phydim_nk, ubound(F_data,2))
               series_pdata(idxvar, nn, 1:nk, kount2) = &
                    facs * series_pdata(idxvar, nn, 1:nk, kount2) + &
                    facf * F_data(i,1:nk)
               series_pdone(idxvar) = 1
               series_pdone2(idxvar,nnl) = 1
            endif
         enddo
      endif IF_PROF
      !---------------------------------------------------------------
      return
   end subroutine series_xst3d

   !#TODO: ??? series_xst4d

   !#---- Private functions ------------------------------------------

   function priv_varidx(F_varname_S, F_varlist_S) result(F_idx)
      character(len=*), intent(in) :: F_varname_S
      character(len=*), intent(in) :: F_varlist_S(:)
      integer :: F_idx, istat
      character(len=SER_STRLEN) :: varname_S
      varname_S = F_varname_S
      istat = clib_toupper(varname_S)
      do F_idx=1,size(F_varlist_S)
         if (varname_S == F_varlist_S(F_idx)) then
            return
         endif
      enddo
      F_idx = RMN_ERR
      return
   end function priv_varidx

end module series_xst_mod
