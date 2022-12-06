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
module mu_fglob_mod
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_isdir, clib_isfile, clib_glob
   implicit none
   private
   !@objective 
   !@author  Stephane Chamberland, 2011-04
   !@description
   ! Public functions
   public :: mu_fglob
   ! Public constants
   integer,parameter,public :: MU_FGLOB_FILE_ONLY = 1
   integer,parameter,public :: MU_FGLOB_DIR_ONLY  = 2
!@/

#include <rmnlib_basics.hf>
!!!#include <arch_specific.hf>

contains

   !/@
   recursive function mu_fglob(F_pathlist_S, F_path_S, F_recurse, F_filter) result(F_npaths)
      implicit none
      !@objective Open rpn std file
      !@arguments
      character(len=*), intent(out) :: F_pathlist_S(:)
      character(len=*), intent(in) :: F_path_S
      integer, intent(in),optional :: F_recurse
      integer, intent(in),optional :: F_filter
      !@author
      !@return
      integer :: F_npaths
      !@/
      integer,parameter :: NMAXFILES = 2048
      logical :: ok_L
      integer :: recurse, filter, n, npaths, npaths2, ifile, nmax, istat
      character(len=RMN_PATH_LEN) :: path_S, pathlist_S(NMAXFILES)
      ! ---------------------------------------------------------------------
      F_npaths = -1
      recurse = 0
      filter = 0
      if (present(F_recurse)) recurse = F_recurse
      if (present(F_filter)) filter = F_filter

      path_S = F_path_S
      if (RMN_IS_OK(clib_isdir(trim(F_path_S)))) then
         n = len_trim(path_S)
         path_S = trim(path_S)//'/*'
         if (path_S(n:n) == '/') path_S = trim(path_S)//'*'
      endif
      nmax = min(size(F_pathlist_S), NMAXFILES)
      npaths = 0
      istat = clib_glob(pathlist_S, npaths, trim(path_S), nmax)
      if (.not.(RMN_IS_OK(istat) .and. npaths >=0)) return 
      ifile = 0
      do n = 1,npaths
         ok_L = .true.
         if (filter == MU_FGLOB_FILE_ONLY) &
              ok_L = RMN_IS_OK(clib_isfile(trim(pathlist_S(n))))
         if (filter == MU_FGLOB_DIR_ONLY) &
              ok_L = RMN_IS_OK(clib_isdir(trim(pathlist_S(n))))
         if (ok_L .and. ifile < nmax) then
            ifile = ifile + 1
            F_pathlist_S(ifile) = pathlist_S(n)
         endif
      end do
      do n = 1,npaths
         if (RMN_IS_OK(clib_isdir(trim(pathlist_S(n)))) &
              .and. recurse > 0) then
            npaths2 = mu_fglob(F_pathlist_S(ifile+1:), pathlist_S(n), recurse-1, filter)
            if (RMN_IS_OK(npaths2)) &
                 ifile = ifile + npaths2
            if (ifile >= nmax) then
               call msg(MSG_WARNING,'(mu_fglob) Too many paths for provided list')
               F_npaths = -nmax
               return
            endif
         endif
      end do
      F_npaths = ifile
      ! ---------------------------------------------------------------------
      return
   end function mu_fglob

end module mu_fglob_mod
