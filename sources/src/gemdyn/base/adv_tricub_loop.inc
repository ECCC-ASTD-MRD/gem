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

      subroutine adv_tricub_loop()

      logical :: zcubic_L

      integer :: id, n, n0, nx, ny, nz, m1, nit, nij, kkmax, &
                 ni_of, nj_of, o1, o2
      real*8  :: ra, rb, rc, rd, dz
      real*8, dimension(4) :: px, py, pz

      real*8 :: triprd, zb, zc, zd
      real :: za
      real :: linz_cubyx, tricub_zyx
      triprd(za,zb,zc,zd)=(za-zb)*(za-zc)*(za-zd)

!--------------------------------------------------------

      kkmax = F_nk - 1

      nit = adv_nit
      nij = adv_nijag
      ni_of = adv_int_i_off
      nj_of = adv_int_j_off

!$omp parallel private(zcubic_L,n0,nx,ny,nz,m1,o1,o2,&
!$omp                  n, id ,px, py, pz, dz, ra,rb,rc,rd) &
!$omp          shared ( kkmax , ii,nit, nij        , &
!$omp                 p_zabcd_8,p_zbacd_8,p_zcabd_8, &
!$omp                 p_zdabc_8,p_bsz_8,p_zbc_8,     &
!$omp                 ni_of, nj_of)
!$omp do
      do id=1,F_nind

         m1=id*4
         n0=m1-3
         nx=m1-2
         ny=m1-1
         nz=m1
         n=ii(n0)

         zcubic_L = (ii(nz) > 0) .and. (ii(nz) < kkmax-1)

         ra = adv_bsx_8(ii(nx)-1)
         rb = adv_bsx_8(ii(nx)  )
         rc = adv_bsx_8(ii(nx)+1)
         rd = adv_bsx_8(ii(nx)+2)
         px(1) = triprd(F_x(n),rb,rc,rd)*adv_xabcd_8
         px(2) = triprd(F_x(n),ra,rc,rd)*adv_xbacd_8
         px(3) = triprd(F_x(n),ra,rb,rd)*adv_xcabd_8
         px(4) = triprd(F_x(n),ra,rb,rc)*adv_xdabc_8

         ra = adv_bsy_8(ii(ny)-1)
         rb = adv_bsy_8(ii(ny)  )
         rc = adv_bsy_8(ii(ny)+1)
         rd = adv_bsy_8(ii(ny)+2)
         py(1) = triprd(F_y(n),rb,rc,rd)*adv_yabcd_8
         py(2) = triprd(F_y(n),ra,rc,rd)*adv_ybacd_8
         py(3) = triprd(F_y(n),ra,rb,rd)*adv_ycabd_8
         py(4) = triprd(F_y(n),ra,rb,rc)*adv_ydabc_8

         o1 = (ii(nz)-1)*nij + (ii(ny)-nj_of-1)*nit + (ii(nx)-ni_of-1) - nit

         if (zcubic_L) then
            ra = p_bsz_8(ii(nz)-1)
            rb = p_bsz_8(ii(nz)  )
            rc = p_bsz_8(ii(nz)+1)
            rd = p_bsz_8(ii(nz)+2)
            pz(1) = triprd(F_z(n),rb,rc,rd)*p_zabcd_8(ii(nz)+1)
            pz(2) = triprd(F_z(n),ra,rc,rd)*p_zbacd_8(ii(nz)+1)
            pz(3) = triprd(F_z(n),ra,rb,rd)*p_zcabd_8(ii(nz)+1)
            pz(4) = triprd(F_z(n),ra,rb,rc)*p_zdabc_8(ii(nz)+1)
            F_cub(n)=  tricub_zyx(F_in(o1),px(1),py(1),pz(1),nit,nij)
         else
            dz = (F_z(n)-p_bsz_8(ii(nz)))*p_zbc_8(ii(nz)+1)
            o2 = o1 + nij
            F_cub(n) = linz_cubyx (F_in(o2),px(1),py(1),dz,nit,nij)
         endif

   end do
!$omp enddo
!$omp end parallel

      end subroutine adv_tricub_loop
