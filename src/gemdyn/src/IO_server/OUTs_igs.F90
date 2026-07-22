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

!**s/p IO_igs - initialize the ig1 and ig2 parameters for a grid

subroutine OUTs_igs ( F_ig1, F_ig2     ,&
     F_longs, F_latgs, F_ni, F_nj      ,&
     F_ig1ro,F_ig2ro, F_ig3ro, F_ig4ro ,&
     F_i0, F_in, F_is, F_j0, F_jn, F_js )
   implicit none

   integer, intent(IN ) :: F_ni, F_nj, F_i0, F_in, F_is, &
                           F_j0, F_jn, F_js            , &
                           F_ig1ro,F_ig2ro,F_ig3ro,F_ig4ro
   integer, intent(OUT) :: F_ig1, F_ig2
   real,    intent(IN ) ::  F_longs(F_ni), F_latgs(F_nj)
   !object
   !     use GRID definition in the interface to define a set of ig1/ig2 values
   !     using a cyclic redundancy check that uniquely define all output grids.
   !
   ! Internal variables
   integer :: i,cnt,crc,dgid,err
   integer, parameter :: ELEM=4
   real :: xlat1,xlon1,xlat2,xlon2
   real, dimension(F_ni*2 + F_nj*2 + ELEM) :: identity_vec
   real, dimension(F_ni,2) :: latlon_we
   real, dimension(F_nj,2) :: latlon_sn

#include <rmnlib_basics.hf>

   integer, external :: f_crc32
!
!----------------------------------------------------------------------
!
   ! Get lat/lon values for the horizontal grid descriptors
   dgid = ezgdef_fmem ( F_ni, 1  , 'Z', 'E',                &
        F_ig1ro, F_ig2ro, F_ig3ro, F_ig4ro, &
        F_longs, F_latgs(F_nj/2) )
   err  = gdll ( dgid, latlon_we, latlon_we(1,2) )

   dgid = ezgdef_fmem ( 1   , F_nj, 'Z', 'E',               &
        F_ig1ro, F_ig2ro, F_ig3ro, F_ig4ro, &
        F_longs(F_ni/2), F_latgs )
   err  = gdll ( dgid, latlon_sn, latlon_sn(1,2) )

   ! Set unique values of F_ig1 and F_ig2 for the descriptors

   F_ig1 = -1; F_ig2 = -1

   cnt = 1
   do i=F_j0,F_jn,F_js
      identity_vec(cnt) = latlon_sn(i,1) ; cnt=cnt+1
      identity_vec(cnt) = latlon_sn(i,2) ; cnt=cnt+1
   enddo
   do i=F_i0,F_in,F_is
      identity_vec(cnt) = latlon_we(i,1) ; cnt=cnt+1
      identity_vec(cnt) = latlon_we(i,2) ; cnt=cnt+1
   enddo

   call cigaxg ( 'E',xlat1,xlon1,xlat2,xlon2, &
                  F_ig1ro,F_ig2ro,F_ig3ro,F_ig4ro )

   identity_vec(cnt:cnt+ELEM-1) = (/xlat1,xlon1,xlat2,xlon2/)
   cnt= cnt+ELEM-1

   crc = f_crc32 (0., identity_vec(1:cnt), cnt*4)
   ! before rmn_011 convip was bugged for 3200 < ip1 < 32768
   ! we therefore add 32768 for now

   F_ig1 = ibits(crc,0,16)  + 32768
   F_ig2 = ibits(crc,16,16) + 32768
!
!----------------------------------------------------------------------
!
   return
   end subroutine OUTs_igs
