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
function serdim(nstn,n,nk) result(F_dim)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
   !@object Determine vector dimensions for time-series
   !@params
   ! nstn     number of stations
   ! n        number of surface or profile variables requested
   ! nk       vertical dimension
   integer,intent(in) :: nstn,n,nk
   !@return
   ! F_dim    dimension to be allocated by the dynamics
   integer :: F_dim
   !@author m. desgagne (mar 99)
   !*@/
   include "series.cdk"
   !---------------------------------------------------------------
   F_dim = 0
   if (nstn <= 0) return
   mxstt = nstn
   if (nk == 1) then
      mxsrf = n
      F_dim = mxstt * mxsrf
   else
      mxprf = n
      mxnvo = nk
      F_dim = mxstt * mxprf * mxnvo
   endif
   !---------------------------------------------------------------
   return
end function serdim
