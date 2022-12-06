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
function file_open_existing(F_fileName_S,F_fileType_S) result(F_fileUnit)
   use, intrinsic :: iso_fortran_env, only: INT64
   use clib_itf_mod, only: clib_isfile, CLIB_OK
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>
   !@objective Open a file, handle error and return unit
   !@arguments
   character(len=*) :: F_fileName_S,F_fileType_S
   !@returns
   integer :: F_fileUnit
   !@author Stephane Chamberland, Nov 2009
   !@revisions
   !  YYYY-MM, MY_AUTHOR: original code
!*@/
   integer :: iErr,nrec=0
   !---------------------------------------------------------------------
   iErr = clib_isfile(F_fileName_S)
   if (iErr == CLIB_OK) then
      F_fileUnit = 0
      iErr = fnom(F_fileUnit, F_fileName_S, F_fileType_S, nrec)
      if (iErr<0) then
         F_fileUnit = RMN_ERR
         call msg(MSG_ERROR,'Unable to open file: '//trim(F_fileName_S))
      endif
   else
      F_fileUnit = RMN_ERR
      call msg(MSG_ERROR,'File does not exists when it should: '//trim(F_fileName_S))
   endif
   !TODO: register a call back to close file on error
   !TODO: if filetype == STD call fstouv() as well
   !---------------------------------------------------------------------
   return
end function file_open_existing

!TODO: create file_close: close file (and FSTFRM) and remove call back
