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

      subroutine adv_BC_LAM_flux_loop_o ()

!$omp parallel private(id,n,n0,nx,ny,nz,m1,zcubic_L,&
!$omp            a1, a2, a3, a4, b1, b2, b3, b4,&
!$omp            c1, c2, c3, c4, d1, d2, d3, d4,&
!$omp            p1, p2, p3, p4, ra, rb, rc, rd,&
!$omp            o1, o2, o3, o4) &
!$omp          shared (p_bsz_8,p_zbc_8,p_zabcd_8, &
!$omp            p_zbacd_8,p_zcabd_8,p_zdabc_8,kkmax, &
!$omp            nind_o,ii_o)

!$omp do

      do id=1,nind_o

            m1=id*4
            n0=m1-3
            nx=m1-2
            ny=m1-1
            nz=m1

            n=ii_o(n0)

            zcubic_L = (ii_o(nz) > 0) .and. (ii_o(nz) < kkmax-1)

            !- x interpolation
            ra = adv_bsx_8(ii_o(nx)-1)
            rb = adv_bsx_8(ii_o(nx)  )
            rc = adv_bsx_8(ii_o(nx)+1)
            rd = adv_bsx_8(ii_o(nx)+2)
            p1 = triprd(F_x(n),rb,rc,rd)*adv_xabcd_8
            p2 = triprd(F_x(n),ra,rc,rd)*adv_xbacd_8
            p3 = triprd(F_x(n),ra,rb,rd)*adv_xcabd_8
            p4 = triprd(F_x(n),ra,rb,rc)*adv_xdabc_8

            o2 = (ii_o(nz)-1)*adv_nijag + (ii_o(ny)-adv_int_j_off-1)*adv_nit + (ii_o(nx)-adv_int_i_off)
            o1 = o2-adv_nit
            o3 = o2+adv_nit
            o4 = o3+adv_nit

            if (zcubic_L) then
               a1 = p1 * F_in_o(o1-1) + p2 * F_in_o(o1) + p3 * F_in_o(o1+1) + p4 * F_in_o(o1+2)
               a2 = p1 * F_in_o(o2-1) + p2 * F_in_o(o2) + p3 * F_in_o(o2+1) + p4 * F_in_o(o2+2)
               a3 = p1 * F_in_o(o3-1) + p2 * F_in_o(o3) + p3 * F_in_o(o3+1) + p4 * F_in_o(o3+2)
               a4 = p1 * F_in_o(o4-1) + p2 * F_in_o(o4) + p3 * F_in_o(o4+1) + p4 * F_in_o(o4+2)
            endif

            o1 = o1 + adv_nijag
            o2 = o2 + adv_nijag
            o3 = o3 + adv_nijag
            o4 = o4 + adv_nijag

            b1 = p1 * F_in_o(o1-1) + p2 * F_in_o(o1) + p3 * F_in_o(o1+1) + p4 * F_in_o(o1+2)
            b2 = p1 * F_in_o(o2-1) + p2 * F_in_o(o2) + p3 * F_in_o(o2+1) + p4 * F_in_o(o2+2)
            b3 = p1 * F_in_o(o3-1) + p2 * F_in_o(o3) + p3 * F_in_o(o3+1) + p4 * F_in_o(o3+2)
            b4 = p1 * F_in_o(o4-1) + p2 * F_in_o(o4) + p3 * F_in_o(o4+1) + p4 * F_in_o(o4+2)

            o1 = o1 + adv_nijag
            o2 = o2 + adv_nijag
            o3 = o3 + adv_nijag
            o4 = o4 + adv_nijag

            c1 = p1 * F_in_o(o1-1) + p2 * F_in_o(o1) + p3 * F_in_o(o1+1) + p4 * F_in_o(o1+2)
            c2 = p1 * F_in_o(o2-1) + p2 * F_in_o(o2) + p3 * F_in_o(o2+1) + p4 * F_in_o(o2+2)
            c3 = p1 * F_in_o(o3-1) + p2 * F_in_o(o3) + p3 * F_in_o(o3+1) + p4 * F_in_o(o3+2)
            c4 = p1 * F_in_o(o4-1) + p2 * F_in_o(o4) + p3 * F_in_o(o4+1) + p4 * F_in_o(o4+2)

            if (zcubic_L) then
               o1 = o1 + adv_nijag
               o2 = o2 + adv_nijag
               o3 = o3 + adv_nijag
               o4 = o4 + adv_nijag

               d1 = p1 * F_in_o(o1-1) + p2 * F_in_o(o1) + p3 * F_in_o(o1+1) + p4 * F_in_o(o1+2)
               d2 = p1 * F_in_o(o2-1) + p2 * F_in_o(o2) + p3 * F_in_o(o2+1) + p4 * F_in_o(o2+2)
               d3 = p1 * F_in_o(o3-1) + p2 * F_in_o(o3) + p3 * F_in_o(o3+1) + p4 * F_in_o(o3+2)
               d4 = p1 * F_in_o(o4-1) + p2 * F_in_o(o4) + p3 * F_in_o(o4+1) + p4 * F_in_o(o4+2)
            endif

            !- y interpolation
            ra = adv_bsy_8(ii_o(ny)-1)
            rb = adv_bsy_8(ii_o(ny)  )
            rc = adv_bsy_8(ii_o(ny)+1)
            rd = adv_bsy_8(ii_o(ny)+2)
            p1 = triprd(F_y(n),rb,rc,rd)*adv_yabcd_8
            p2 = triprd(F_y(n),ra,rc,rd)*adv_ybacd_8
            p3 = triprd(F_y(n),ra,rb,rd)*adv_ycabd_8
            p4 = triprd(F_y(n),ra,rb,rc)*adv_ydabc_8

            b1 = p1 * b1 + p2 * b2 + p3 * b3 + p4 * b4
            c1 = p1 * c1 + p2 * c2 + p3 * c3 + p4 * c4

            !- z interpolation
            if (zcubic_L) then
               a1 = p1 * a1 + p2 * a2 + p3 * a3 + p4 * a4
               d1 = p1 * d1 + p2 * d2 + p3 * d3 + p4 * d4
               ra = p_bsz_8(ii_o(nz)-1)
               rb = p_bsz_8(ii_o(nz)  )
               rc = p_bsz_8(ii_o(nz)+1)
               rd = p_bsz_8(ii_o(nz)+2)
               p1 = triprd(F_z(n),rb,rc,rd)*p_zabcd_8(ii_o(nz)+1)
               p2 = triprd(F_z(n),ra,rc,rd)*p_zbacd_8(ii_o(nz)+1)
               p3 = triprd(F_z(n),ra,rb,rd)*p_zcabd_8(ii_o(nz)+1)
               p4 = triprd(F_z(n),ra,rb,rc)*p_zdabc_8(ii_o(nz)+1)

               F_cub_o(n) = p1 * a1 + p2 * b1 + p3 * c1 + p4 * d1
            else
               p3 = (F_z(n)-p_bsz_8(ii_o(nz)))*p_zbc_8(ii_o(nz)+1)
               p2 = 1. - p3

               F_cub_o(n) = p2 * b1 + p3 * c1
            endif
      end do
!$omp enddo
!$omp end parallel

      end subroutine adv_BC_LAM_flux_loop_o
