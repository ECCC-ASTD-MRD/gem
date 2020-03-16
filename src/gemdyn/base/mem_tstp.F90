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
module mem_tstp
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      real, allocatable,  dimension (:,:  ) :: rhsb
      real, allocatable,  dimension (:,:,:) :: rhsu,rhsv,rhst,rhsc,rhsw,rhsf,rhsp
!      real, allocatable,  dimension (:,:,:) :: ruw1 
!      real, allocatable,  dimension (:,:,:) :: rvw1 
!      real, allocatable,  dimension (:,:,:) :: ruw2 
!      real, allocatable,  dimension (:,:,:) :: rvw2 
!      real, allocatable,  dimension (:,:,:) :: rhsx 
!      real, allocatable,  dimension (:,:,:) :: rhsq 

      real, allocatable,  dimension (:,:,:) :: orhsu,orhsv,orhst,orhsc,orhsw,orhsf

      real, allocatable,  dimension (:,:  ) :: nl_b
      real, allocatable,  dimension (:,:,:) :: nl_u,nl_v,nl_t,nl_c,nl_w,nl_f

      real(kind=REAL64), allocatable, dimension (:,:,:) :: rhs_sol, lhs_sol

end module mem_tstp
