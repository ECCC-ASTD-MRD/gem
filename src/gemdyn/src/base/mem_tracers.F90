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
module mem_tracers
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

      integer :: tracers_block_size = 4
      integer :: tracers_nblocks
      real, dimension    (:), pointer :: trt2, trt1, trt0, trdf, trtb
      real(kind=REAL64), dimension(:,:,:), pointer :: sumq_8

      type :: memTR_pntrs
         real, dimension(:,:,:), pointer :: pntr
      end type memTR_pntrs

      type(memTR_pntrs), allocatable :: tracers_P (:)
      type(memTR_pntrs), allocatable :: tracers_M (:)
      type(memTR_pntrs), allocatable :: tracers_t2 (:)
      type(memTR_pntrs), allocatable :: tracers_B (:)

contains
      integer function tr_get(F_name_S, F_pntr)
      use glb_ld
      use tr3d
      implicit none

      character(len=*), intent(IN) :: F_name_S
      real, pointer, dimension (:,:,:), intent(OUT) :: F_pntr

      character(len=2 ) tf
      character(len=64) name
      integer n, indx,dim

      dim = (l_maxx-l_minx+1) * (l_maxy-l_miny+1) * l_nk
      tr_get = -99999 ; nullify (F_pntr)
      indx = index(F_name_s,":")
      name= F_name_s ; tf= 'P'
      if (indx > 0) then
         name = F_name_S(1:indx-1)
         tf   = F_name_S(indx+1: )
      endif
      do n=1,Tr3d_ntr
         if (trim(Tr3d_name_S(n)) == trim(name) )then
            tr_get=n
            select case(tf)
            case('M' , 'm')
               F_pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt0((n-1)*dim+1:)
            case('P' , 'p')
               F_pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt1((n-1)*dim+1:)
            case('T2' , 't2')
               F_pntr(l_minx:l_maxx,l_miny:l_maxy,1:l_nk) => trt2((n-1)*dim+1:)
            case default
               tr_get= -99999
               print*, 'Timeframe ',tf,' for tracer ',trim(name),' NOT FOUND'
            end select
         endif
      end do
      if (tr_get<0) then
         print*, 'Tracer ',trim(name),' NOT FOUND'
      endif

      return
      end function tr_get

end module mem_tracers
