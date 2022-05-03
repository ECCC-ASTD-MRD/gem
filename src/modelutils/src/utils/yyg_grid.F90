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


subroutine grid_info(F_ni, F_nj, F_div, F_arakawa, F_gid, &
     F_searchkey, F_iun)
   implicit none
   !@object Decoding of the g4 parameter that contains the nomvar
   !        of the positional parameter and the Arakawa C position.
   !@arguments
   integer :: F_ni, F_nj, F_div, F_arakawa, F_gid, F_searchkey, F_iun
   !@author Michel Desgagne  -  Spring 2012

   integer, external :: fstprm, fstinl
   integer, parameter :: nlis = 1024
   integer, parameter :: ngrids = 100

   character(len=2)  :: typ_S, grd_S
   character(len=4)  :: var_S
   character(len=12) :: lab_S

   integer :: dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
        dty, swa, lng, dlf, ubc, ex1, ex2, ex3
   integer :: ni1, nj1, nk1, err,lislon
   integer :: liste(nlis)

   character(len=12) :: grid_type(ngrids)
   integer :: grid_key(ngrids)
   common /grid_c/ grid_type
   common /grid_i/ grid_key
   !---------------------------------------------------------------
   F_gid = -1
   err= fstprm (F_searchkey, dte, det, ipas, ni1, nj1, nk1, bit, &
        dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1, &
        g2, g3, g4, swa, lng, dlf, ubc, ex1, ex2, ex3)

   ! Normally var_S and F_arakawa would be encoded into g4
   ! (encoding/decoding function are required)

   if (g4 == 1) then
      F_gid             = 1
      grid_type(F_gid)  = 'YINYANG'
      var_S             = '<Y0>'
      F_arakawa         = -1
      F_div             = 2
      F_ni              = ni1
      F_nj              = nj1/2
      err               = fstinl (F_iun,ni1,nj1,nk1,-1,' ',g1,g2,g3,' ', &
           var_S,liste,lislon,nlis)
      grid_key(F_gid)   = liste(1)
   endif
   !---------------------------------------------------------------
   return
end subroutine grid_info


subroutine grid_getll(F_posx, F_posy, F_corners, F_yin_rot, F_yan_rot, F_gid)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@object Obtain positional parameters
   !@arguments
   integer :: F_yin_rot(4), F_yan_rot(4), F_gid
   real(REAL64) :: F_posx(*), F_posy(*)
   real :: F_corners(4)
   !@author Michel Desgagne  -  Spring 2012

!!$   integer, external :: fstprm, fstluk
   integer, parameter :: ngrids = 100

   character(len=2)  :: typ_S, grd_S
   character(len=4)  :: var_S
   character(len=12) :: lab_S

   integer :: dte, det, ipas, p1, p2, p3, g1, g2, g3, g4, bit, &
        dty, swa, lng, dlf, ubc, ex1, ex2, ex3
   integer ni1, nj1, nk1, err
   real, dimension(:), allocatable :: yy

   character(len=12) :: grid_type(ngrids)
   integer :: grid_key(ngrids)
   common /grid_c/ grid_type
   common /grid_i/ grid_key
   !---------------------------------------------------------------
   err = fstprm(grid_key(F_gid), dte, det, ipas, ni1, nj1, nk1, bit, &
        dty, p1, p2, p3, typ_S, var_S, lab_S, grd_S, g1, &
        g2, g3, g4, swa,  lng, dlf, ubc, ex1, ex2, ex3)

   allocate(yy(ni1))
   err = fstluk(yy,grid_key(F_gid),ni1,nj1,nk1)

   if (grid_type(F_gid) == 'YINYANG') then

      ni1 = nint (yy(1))
      nj1 = nint (yy(2))

      F_corners(1) = yy(3)
      F_corners(2) = yy(4)
      F_corners(3) = yy(5)
      F_corners(4) = yy(6)
      F_yin_rot(1) = nint (yy(7 ))
      F_yin_rot(2) = nint (yy(8 ))
      F_yin_rot(3) = nint (yy(9 ))
      F_yin_rot(4) = nint (yy(10))
      F_yan_rot(1) = nint (yy(11))
      F_yan_rot(2) = nint (yy(12))
      F_yan_rot(3) = nint (yy(13))
      F_yan_rot(4) = nint (yy(14))

      print*, 'YY decoded value:'
      print*, 'GNI=', nint (yy(1) )
      print*, 'GNJ=', nint (yy(2) )
      print*, 'BOTTOM/LEF=', yy(3),yy(4)
      print*, 'UPPER/RIGHT', yy(5),yy(6)
      print*, 'YIN_ROT=', F_yin_rot
      print*, 'YAN_ROT=', F_yan_rot
      F_posx(1:ni1) = yy(15:14+ni1)
      F_posy(1:nj1) = yy(15+ni1:)

      ! call UOG_parpos(F_posx, F_posy,                                     &
      !                 dble(yy(3)), dble(yy(4)), dble(yy(5)), dble(yy(6)), &
      !                 ni1,nj1)

   endif

   deallocate(yy, stat=err)
   !---------------------------------------------------------------
   return
end subroutine grid_getll


subroutine UOG_parpos(F_x_8, F_y_8, &
     F_xbeg, F_ybeg, F_xend, F_yend, &
     NX, NY)
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@object Compute Yin-Yang grid positional parameters
   !@arguments
   ! F_x_8        O     - positions alog-x
   ! F_y_8        O     - positions alog-y
   ! F_xbeg       I     - starting longitude degree of the x-axis
   ! F_ybeg       I     - starting latitude  degree of the y-axis
   ! F_xend       I     - ending longitude degree of the x-axis
   ! F_yend       I     - ending latitude  degree of the y-axis
   ! NX           I     - total number of points on the x-axis
   ! NY           I     - total number of points on the y-axis
   integer :: NX, NY
   real(REAL64) ::  F_x_8(NX), F_y_8(NY),F_xbeg, F_ybeg, F_xend, F_yend
   !@author  Michel Desgagne  - Spring 2012

   integer i
   real(REAL64) :: delta_8
   !-----------------------------------------------------------
   delta_8 = (F_xend-F_xbeg)/dble(NX-1)
   F_x_8(1)  = F_xbeg
   F_x_8(NX) = F_xend

   do i=2,NX-1
      F_x_8(i)= F_xbeg + (i-1)*delta_8
   end do

   delta_8 = (F_yend-F_ybeg)/dble(NY-1)
   F_y_8(1 ) = F_ybeg
   F_y_8(NY) = F_yend

   do i=2,NY-1
      F_y_8(i)= F_ybeg + (i-1)*delta_8
   end do
   !-----------------------------------------------------------
   return
end subroutine UOG_parpos


subroutine yyg_yangrot ( F_yinlat1, F_yinlon1, F_yinlat2, F_yinlon2, &
     F_yanlat1, F_yanlon1, F_yanlat2, F_yanlon2  )
   use, intrinsic :: iso_fortran_env, only: REAL64
   implicit none
   include "rmnlib_basics.inc"
   !@object return the rotation for the Yang grid using the rotation from Yin
   !@arguments
   real(REAL64) :: F_yinlat1, F_yinlon1, F_yinlat2, F_yinlon2, &
        F_yanlat1, F_yanlon1, F_yanlat2, F_yanlon2
   !@author V. Lee/A. Qaddouri - April 2011
   !@revision
   ! v4_40 - Qaddouri/Lee      - To obtain Yang grid

   real(REAL64) :: xlat1,xlon1,xlat2,xlon2
   real(REAL64) :: a_8,b_8,c_8,d_8,xyz1_8(3),xyz2_8(3)
   real(REAL64) :: rot_8(3,3),invrot_8(3,3),xyz3_8(3),xyz4_8(3)
   integer :: i,j
   !-------------------------------------------------------------------

   !Get the rotation for the input grid (F_yinlat1,F_yinlon1,F_yinlat2,F_yinlon2)
   call llacar_8(xyz1_8, F_yinlon1, F_yinlat1, 1, 1)
   call llacar_8(xyz2_8, F_yinlon2, F_yinlat2, 1, 1)
   a_8 = (xyz1_8(1)*xyz2_8(1)) + (xyz1_8(2)*xyz2_8(2))  &
        + (xyz1_8(3)*xyz2_8(3))
   b_8 = sqrt (((xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3)))**2 &
        +  ((xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3)))**2 &
        +  ((xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2)))**2)
   c_8 = sqrt ( xyz1_8(1)**2 + xyz1_8(2)**2 + xyz1_8(3)**2 )
   d_8 = sqrt ( ( ( (a_8*xyz1_8(1)) - xyz2_8(1) ) / b_8 )**2 + &
        ( ( (a_8*xyz1_8(2)) - xyz2_8(2) ) / b_8 )**2 + &
        ( ( (a_8*xyz1_8(3)) - xyz2_8(3) ) / b_8 )**2  )
   rot_8(1,1)=  -xyz1_8(1)/c_8
   rot_8(1,2)=  -xyz1_8(2)/c_8
   rot_8(1,3)=  -xyz1_8(3)/c_8
   rot_8(2,1)=  ( ((a_8*xyz1_8(1)) - xyz2_8(1)) / b_8)/d_8
   rot_8(2,2)=  ( ((a_8*xyz1_8(2)) - xyz2_8(2)) / b_8)/d_8
   rot_8(2,3)=  ( ((a_8*xyz1_8(3)) - xyz2_8(3)) / b_8)/d_8
   rot_8(3,1)=   &
        ( (xyz1_8(2)*xyz2_8(3)) - (xyz2_8(2)*xyz1_8(3)))/b_8
   rot_8(3,2)=   &
        ( (xyz2_8(1)*xyz1_8(3)) - (xyz1_8(1)*xyz2_8(3)))/b_8
   rot_8(3,3)=   &
        ( (xyz1_8(1)*xyz2_8(2)) - (xyz2_8(1)*xyz1_8(2)))/b_8

   !Get transpose of rotation
   do i=1,3
      do j=1,3
         invrot_8(i,j)=rot_8(j,i)
      enddo
   enddo

   ! Find the centre of Yang grid through Yin by setting
   xlon1 = 0.0d0
   xlat1 = 0.0d0
   ! And set the rotation for Yang grid with respect to Yin
   xlon2 = 0.0d0
   xlat2 = 90.0d0

   ! Obtain the cartesian coordinates
   call llacar_8(xyz1_8, xlon1, xlat1, 1, 1)
   call llacar_8(xyz2_8, xlon2, xlat2, 1, 1)
   !     call mxma8(invrot_8,1,3,xyz1_8,1,3,xyz3_8,1,3,3,3,1)
   !     call mxma8(invrot_8,1,3,xyz2_8,1,3,xyz4_8,1,3,3,3,1)
   !     call dgemm('N','N',1,3,3,1._REAL64,xyz1_8,1,invrot_8,1, 0._REAL64,xyz3_8,1)
   !     call dgemm('N','N',1,3,3,1._REAL64,xyz2_8,1,invrot_8,1, 0._REAL64,xyz4_8,1)
   ! Only multiply a 2d 3x3 matrice to a single vector:
   do i=1,3
      xyz3_8(i)=0.0d0
      xyz4_8(i)=0.0d0
      do j=1,3
         xyz3_8(i) = xyz3_8(i) + invrot_8(i,j)*xyz1_8(j)
         xyz4_8(i) = xyz4_8(i) + invrot_8(i,j)*xyz2_8(j)
      enddo
   enddo

   !Obtain the real geographic coordinates
   call cartall_8(xlon1,xlat1, xyz3_8,1)
   call cartall_8(xlon2,xlat2, xyz4_8,1)
   if (xlon1 >= 360.) xlon1=xlon1-360.0
   if (xlon2 >= 360.) xlon2=xlon2-360.0

   F_yanlat1= xlat1
   F_yanlon1= xlon1
   F_yanlat2= xlat2
   F_yanlon2= xlon2
   !-------------------------------------------------------------------
   return
end subroutine yyg_yangrot
