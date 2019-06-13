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

      subroutine adv_BC_LAM_flux_loop_iH ()

      integer :: m, midxk, idxk, midxjk, idxjk, &
                 o1, o2, o3, o4 , ii, jj, kk, ii1, jj1, kk1, sig

      real(kind=REAL64)  :: a1, a2, a3, a4, b1, b2, b3, b4, &
                 c1, c2, c3, c4, d1, d2, d3, d4, &
                 p1, p2, p3, p4, &
                 ra, rb, rc, rd, rri, rrj, rrk, p_z00_8

      logical :: zcubic_L

      real(kind=REAL64) :: triprd, zb, zc, zd
      real :: za
      triprd(za,zb,zc,zd)=(za-zb)*(za-zc)*(za-zd)

      !------------------------------------------------------

      !Vertical variable type:  Height--> sig <0 , Pressure --> sig >0
      sig=int((Ver_z_8%m(l_nk)-Ver_z_8%m(1))/(abs(  Ver_z_8%m(l_nk)-Ver_z_8%m(1) )))

      p_z00_8 = Ver_z_8%m(0)

!$omp parallel private(m,midxk,idxk,midxjk,idxjk,zcubic_L,&
!$omp            a1, a2, a3, a4, b1, b2, b3, b4,&
!$omp            c1, c2, c3, c4, d1, d2, d3, d4,&
!$omp            p1, p2, p3, p4, ra, rb, rc, rd,&
!$omp            o1, o2, o3, o4, ii, jj, kk,&
!$omp            ii1,jj1,kk1,rri,rrj,rrk)&
!$omp          shared (p_bsz_8,p_zbc_8,p_zabcd_8,i_h,p_lcz,&
!$omp            p_zbacd_8,p_zcabd_8,p_zdabc_8,kkmax,p_z00_8,sig,&
!$omp            i0_i,in_i,j0_i,jn_i,minx_e,maxx_e,miny_e,maxy_e,&
!$omp            px,py,pz)

!$omp do

      do k=F_k0,F_nk

         midxk = (k-1)*(maxx_e-minx_e+1)*(maxy_e-miny_e+1)
          idxk = (k-1)*l_ni*l_nj

         do j=j0_i,jn_i

            midxjk = midxk + ((j-miny_e)*(maxx_e-minx_e+1))
             idxjk =  idxk + ((j-1)*l_ni)

            do i=i0_i,in_i

               m = midxjk + i-minx_e + 1 !Location in MIN/MAX grid of a point in I0/IN/J0/JN grid

               rri = px(m)
               ii1 = 1 + (rri - adv_x00_8) * adv_ovdx_8
               ii  = max(2,min(ii1,adv_iimax))

               rrj = py(m)
               jj1 = 1 + (rrj - adv_y00_8) * adv_ovdy_8
               jj  = max(G_haloy,min(jj1,adv_jjmax))

               rrk = pz(m)
               kk1 = (rrk - p_z00_8) * adv_ovdz_8*sig
               kk1 = p_lcz (kk1+1)
               if ( real(sig) * (rrk - p_bsz_8(kk1)) < 0. ) kk1 = kk1 - 1
               kk  = min(kkmax-1,max(0,kk1))

               zcubic_L = (kk > 0) .and. (kk < kkmax-1)

               !- x interpolation
               ra = adv_bsx_8(ii-1)
               rb = adv_bsx_8(ii  )
               rc = adv_bsx_8(ii+1)
               rd = adv_bsx_8(ii+2)
               p1 = triprd(px(m),rb,rc,rd)*adv_xabcd_8
               p2 = triprd(px(m),ra,rc,rd)*adv_xbacd_8
               p3 = triprd(px(m),ra,rb,rd)*adv_xcabd_8
               p4 = triprd(px(m),ra,rb,rc)*adv_xdabc_8

               o2 = (kk-1)*adv_nijag + (jj-adv_int_j_off-1)*adv_nit + (ii-adv_int_i_off)
               o1 = o2-adv_nit
               o3 = o2+adv_nit
               o4 = o3+adv_nit

               if (zcubic_L) then
                  a1 = p1 * F_in_i(o1-1) + p2 * F_in_i(o1) + p3 * F_in_i(o1+1) + p4 * F_in_i(o1+2)
                  a2 = p1 * F_in_i(o2-1) + p2 * F_in_i(o2) + p3 * F_in_i(o2+1) + p4 * F_in_i(o2+2)
                  a3 = p1 * F_in_i(o3-1) + p2 * F_in_i(o3) + p3 * F_in_i(o3+1) + p4 * F_in_i(o3+2)
                  a4 = p1 * F_in_i(o4-1) + p2 * F_in_i(o4) + p3 * F_in_i(o4+1) + p4 * F_in_i(o4+2)
               end if

               o1 = o1 + adv_nijag
               o2 = o2 + adv_nijag
               o3 = o3 + adv_nijag
               o4 = o4 + adv_nijag

               b1 = p1 * F_in_i(o1-1) + p2 * F_in_i(o1) + p3 * F_in_i(o1+1) + p4 * F_in_i(o1+2)
               b2 = p1 * F_in_i(o2-1) + p2 * F_in_i(o2) + p3 * F_in_i(o2+1) + p4 * F_in_i(o2+2)
               b3 = p1 * F_in_i(o3-1) + p2 * F_in_i(o3) + p3 * F_in_i(o3+1) + p4 * F_in_i(o3+2)
               b4 = p1 * F_in_i(o4-1) + p2 * F_in_i(o4) + p3 * F_in_i(o4+1) + p4 * F_in_i(o4+2)

               o1 = o1 + adv_nijag
               o2 = o2 + adv_nijag
               o3 = o3 + adv_nijag
               o4 = o4 + adv_nijag

               c1 = p1 * F_in_i(o1-1) + p2 * F_in_i(o1) + p3 * F_in_i(o1+1) + p4 * F_in_i(o1+2)
               c2 = p1 * F_in_i(o2-1) + p2 * F_in_i(o2) + p3 * F_in_i(o2+1) + p4 * F_in_i(o2+2)
               c3 = p1 * F_in_i(o3-1) + p2 * F_in_i(o3) + p3 * F_in_i(o3+1) + p4 * F_in_i(o3+2)
               c4 = p1 * F_in_i(o4-1) + p2 * F_in_i(o4) + p3 * F_in_i(o4+1) + p4 * F_in_i(o4+2)

               if (zcubic_L) then
                  o1 = o1 + adv_nijag
                  o2 = o2 + adv_nijag
                  o3 = o3 + adv_nijag
                  o4 = o4 + adv_nijag

                  d1 = p1 * F_in_i(o1-1) + p2 * F_in_i(o1) + p3 * F_in_i(o1+1) + p4 * F_in_i(o1+2)
                  d2 = p1 * F_in_i(o2-1) + p2 * F_in_i(o2) + p3 * F_in_i(o2+1) + p4 * F_in_i(o2+2)
                  d3 = p1 * F_in_i(o3-1) + p2 * F_in_i(o3) + p3 * F_in_i(o3+1) + p4 * F_in_i(o3+2)
                  d4 = p1 * F_in_i(o4-1) + p2 * F_in_i(o4) + p3 * F_in_i(o4+1) + p4 * F_in_i(o4+2)
               end if

               !- y interpolation
               ra = adv_bsy_8(jj-1)
               rb = adv_bsy_8(jj  )
               rc = adv_bsy_8(jj+1)
               rd = adv_bsy_8(jj+2)
               p1 = triprd(py(m),rb,rc,rd)*adv_yabcd_8
               p2 = triprd(py(m),ra,rc,rd)*adv_ybacd_8
               p3 = triprd(py(m),ra,rb,rd)*adv_ycabd_8
               p4 = triprd(py(m),ra,rb,rc)*adv_ydabc_8

               b1 = p1 * b1 + p2 * b2 + p3 * b3 + p4 * b4
               c1 = p1 * c1 + p2 * c2 + p3 * c3 + p4 * c4

               !- z interpolation
               if (zcubic_L) then
                  a1 = p1 * a1 + p2 * a2 + p3 * a3 + p4 * a4
                  d1 = p1 * d1 + p2 * d2 + p3 * d3 + p4 * d4
                  ra = p_bsz_8(kk-1)
                  rb = p_bsz_8(kk  )
                  rc = p_bsz_8(kk+1)
                  rd = p_bsz_8(kk+2)
                  p1 = triprd(pz(m),rb,rc,rd)*p_zabcd_8(kk+1)
                  p2 = triprd(pz(m),ra,rc,rd)*p_zbacd_8(kk+1)
                  p3 = triprd(pz(m),ra,rb,rd)*p_zcabd_8(kk+1)
                  p4 = triprd(pz(m),ra,rb,rc)*p_zdabc_8(kk+1)

                  i_h(i,j,k) = p1 * a1 + p2 * b1 + p3 * c1 + p4 * d1
               else
                  p3 = (pz(m)-p_bsz_8(kk))*p_zbc_8(kk+1)
                  p2 = 1. - p3

                  i_h(i,j,k) = p2 * b1 + p3 * c1
               end if

            end do
         end do
      end do
!$omp end do
!$omp end parallel

      end subroutine adv_BC_LAM_flux_loop_iH
