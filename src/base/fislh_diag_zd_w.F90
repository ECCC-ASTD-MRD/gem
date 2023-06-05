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

!** s/r - Computes model vertical velocities zd and w diagnostically.
!         Height-type vertical coordinate

      subroutine fislh_diag_zd_w (F_zd, F_w, F_u, F_v, F_t, F_q, &
              F_metric,Minx, Maxx, Miny, Maxy, Nk, F_zd_L, F_w_L )
      use dyn_fisl_options
      use gem_options
      use geomh
      use glb_ld
      use lun
      use metric
      use tdpack
      use ver
      implicit none

      integer, intent(in) ::  Minx, Maxx, Miny, Maxy, Nk
      logical, intent(in) ::  F_zd_L, F_w_L
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(out)   :: F_zd, F_w
      real, dimension(Minx:Maxx,Miny:Maxy,Nk),  intent(inout) :: F_u, F_v, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1),intent(inout) :: F_q
      type(Vmetric) , intent(in ) :: F_metric

      integer :: i, j, k, kp, km, i0, in, j0, jn
      real, dimension(Minx:Maxx,Miny:Maxy)    :: rJzX, rJzY, rJzZ
!     ________________________________________________________________
!
      if ( (F_zd_L.or.F_w_L) ) then
      call rpn_comm_xch_halo (F_u,  l_minx,l_maxx,l_miny,l_maxy, l_niu,l_nj, Nk,   &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo (F_v,  l_minx,l_maxx,l_miny,l_maxy, l_ni,l_njv, Nk,   &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo( F_t,  l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj ,Nk,   &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo( F_q,  l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj ,Nk+1, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      endif
                           
!     local grid setup for final results
      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj
      if (l_west)  i0 = 2    - G_halox
      if (l_east)  in = in-1 + G_halox
      if (l_south) j0 = 2    - G_haloy
      if (l_north) jn = jn-1 + G_haloy

      if (F_zd_L) then
         if (Lun_debug_L.and.F_zd_L) write (Lun_out,1000)

! Compute Zdot

         F_zd=0.0
         rJzZ=0.0
         do k=Nk,2,-1
            km=max(k-1,1)
            do j=j0,jn
               do i=i0-1,in
                  rJzX(i,j) = 0.5d0*((F_metric%ztht_8(i+1,j,k)-F_metric%ztht_8(i+1,j,km)) &
                            + (F_metric%ztht_8(i  ,j,k)-F_metric%ztht_8(i  ,j,km)))*Ver_idz_8%m(k)
               end do
            end do
            do j=j0-1,jn
               do i=i0,in
                  rJzY(i,j) = 0.5d0*((F_metric%ztht_8(i,j+1,k)-F_metric%ztht_8(i,j+1,km)) &
                            + (F_metric%ztht_8(i,j  ,k)-F_metric%ztht_8(i,j  ,km)))*Ver_idz_8%m(k)
               end do
            end do
            do j=j0,jn
               do i=i0,in
                  F_zd(i,j,k-1) = rJzZ(i,j)*F_zd(i,j,k)*Ver_zeronk(k) + Ver_dz_8%m(k)*( &
                                (rJzX(i  ,j)*F_u(i  ,j,k)                   &
                               -rJzX(i-1,j)*F_u(i-1,j,k))*geomh_invDX_8(j) &
                              +(rJzY(i,j  )*F_v(i,j  ,k)*geomh_cyV_8(j  )  &
                               -rJzY(i,j-1)*F_v(i,j-1,k)*geomh_cyV_8(j-1))*geomh_invDYM_8(j) )
                  rJzZ(i,j) = (F_metric%zmom_8(i,j,k)-F_metric%zmom_8(i,j,k-1))*Ver_idz_8%t(k-1)
                  F_zd(i,j,k-1) = F_zd(i,j,k-1)/rJzZ(i,j)
               end do
            end do
         end do

      end if

      if(F_w_L) then
         if (Lun_debug_L.and.F_w_L ) write (Lun_out,1001)

! Compute W (which depends on previous computation of Zdot)

         i0 = 1
         in = l_ni
         j0 = 1
         jn = l_nj
         if (l_west)  i0 = 3    - G_halox
         if (l_east)  in = in-1 + G_halox
         if (l_south) j0 = 3    - G_haloy
         if (l_north) jn = jn-1 + G_haloy
         F_w=0.
         do k=1,Nk
            km=max(k-1,1)
            kp=min(k+1,Nk)
            do j=j0,jn
               do i=i0,in
                  F_w(i,j,k) = 0.25d0* ( &
                              (F_u(i,j,kp)*F_metric%mc_Jx_8(i,j,kp)+F_u(i-1,j,kp)*F_metric%mc_Jx_8(i-1,j,kp))   &
                             +(F_u(i,j,k )*F_metric%mc_Jx_8(i,j,k )+F_u(i-1,j,k )*F_metric%mc_Jx_8(i-1,j,k ))   &
                             +(F_v(i,j,kp)*F_metric%mc_Jy_8(i,j,kp)+F_v(i,j-1,kp)*F_metric%mc_Jy_8(i,j-1,kp))   &
                             +(F_v(i,j,k )*F_metric%mc_Jy_8(i,j,k )+F_v(i,j-1,k )*F_metric%mc_Jy_8(i,j-1,k )) ) &
                             +(Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km))  &
                             *(F_metric%zmom_8(i,j,k+1)-F_metric%zmom_8(i,j,k))*Ver_idz_8%t(k)
               end do
            end do
         end do
      end if

 1000 format(3X,'COMPUTE DIAGNOSTIC ZDT1: (S/R DIAG_ZD_W_H)')
 1001 format(3X,'COMPUTE DIAGNOSTIC WT1:  (S/R DIAG_ZD_W_H)')
!
!     ________________________________________________________________
!
      return
      end subroutine fislh_diag_zd_w
