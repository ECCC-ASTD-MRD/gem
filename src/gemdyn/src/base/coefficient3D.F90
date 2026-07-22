!**s/r coefficient3D  -
!
      subroutine coefficient3D (F_s, F_u1, F_v1, F_w1,  lminx, lmaxx, lminy, lmaxy, nk)

!      use cstv
!      use HORgrid_options
!      use dynkernel_options
!      use hvdif_options
!      use gmm_itf_mod
!      use gmm_smag
!      use hzd_mod
!      use tdpack

      use gem_options
      use geomh
      use glb_ld
      use cstv
      use ver
      use metric
      use hzd_mod
      use dcst
!
      use ptopo
      use stat_mpi, only: statf_dm
!
      use, intrinsic :: iso_fortran_env
      implicit none


      integer, intent(in) :: lminx,lmaxx,lminy,lmaxy, nk
      real, dimension(lminx:lmaxx,lminy:lmaxy,nk), intent(in) :: F_u1, F_v1, F_w1
      real, dimension(lminx:lmaxx,lminy:lmaxy,nk),    intent(out) :: F_s
      real(kind=REAL64), parameter :: half=0.5d0

      integer i0,in,j0,jn,i,j,k
      real(kind=REAL64)  C1, C2, C3, C4
      real(kind=REAL64) UT1_t, UT1_tmi, UT1_phi, UT1p1_phi
      real(kind=REAL64) VT1_t, VT1_tmj, VT1_phi, VT1p1_phi,VT1_t_mj, VT1_t_phi
      real(kind=REAL64) WT1_u, WT1_um1, WT1_m,  WT1p1_m 
      real(kind=REAL64) WT1_v,WT1_vm1
      real(kind=REAL64) VT1_tr, VT1_tr_p1, VT1_lr, VT1_lr_p1, VT1_tl, VT1_tl_p1, VT1_ll, VT1_ll_p1
      real(kind=REAL64) UT1_tr, UT1_tr_p1, UT1_lr, UT1_lr_p1, UT1_tl, UT1_tl_p1, UT1_ll, UT1_ll_p1
!      real, dimension(lminx:lmaxx,lminy:lmaxy,nk) :: k1_h, k2_h, k3_h, k4_h, k5_h 
      real(kind=REAL64) k1_h, k2_h, k3_h, k4_h, k5_h
      real, dimension(lminx:lmaxx,lminy:lmaxy,0:nk+1) :: F_u, F_v ,F_w
      real(kind=REAL64) Cjx(lminx:lmaxx,lminy:lmaxy),Cjy(lminx:lmaxx,lminy:lmaxy)
!
      i0  = 1    + pil_w
      in  = l_ni - pil_e

      j0  = 1    + pil_s
      jn  = l_nj - pil_n


      F_u = 0.0
      F_v = 0.0
      F_w = 0.0

      F_u(:,:,1:nk)= F_u1(:,:,1:nk)
      F_v(:,:,1:nk)= F_v1(:,:,1:nk)
      F_w(:,:,1:nk)=F_w1(:,:,1:nk)
!  other boundary condition 
       F_u(:,:,0)=F_u(:,:,1)
       F_v(:,:,0)=F_v(:,:,1)
!
       F_u(:,:,nk+1)= F_u1(:,:,nk)
       F_v(:,:,nk+1)= F_v1(:,:,nk)
       do j=j0-1, jn+1
          do i=i0-1, in+1
       Cjx(i,j)=  (GVM%zmom_8(i+1,j,nk+1)-GVM%zmom_8(i,j,nk+1))*geomh_invDX_8(j)
       Cjy(i,j)=  (GVM%zmom_8(i,j+1,nk+1)-GVM%zmom_8(i,j,nk+1))*geomh_invDY_8
       F_w(i,j,nk+1)=Cjx(i,j) *F_u(i,j,nk+1)+Cjy(i,j)*F_v(i,j,nk+1)
          enddo
       enddo


! UT1(i,j,k): U-wind component on u-point and m-levels
! VT1(i,j,k): V-wind component on v-point and m-levels
! WT1(i,j,k): dz/dt (not d\zeta / dt) on phi-points and t-levels

! Computation on phi-points and t-levels of SQRT(S^2)

      F_s=0.0

      do k=1,nk
         do j=j0-1, jn+1
            do i=i0-1, in+1
! u on k t-levels and (i,j) u-points
      UT1_t  = half*(F_u(i,j,k)+F_u(i,j,k+1))
! u on k t-levels and (i-1,j) u-points
      UT1_tmi = half*(F_u(i-1,j,k)+F_u(i-1,j,k+1))
! du/dx on k t-levels and (i,j) phi-points
      C1  = ( UT1_t - UT1_tmi ) * geomh_invDX_8(j)
! u on k m-levels and (i,j) phi-points
      UT1_phi = half * (F_u(i,j,k)+F_u(i-1,j,k))
! u on k+1 m-levels and (i,j) phi-points
      UT1p1_phi = half * (F_u(i,j,k+1)+F_u(i-1,j,k+1))
! (dz/dx)*(du/dz) on k t-levels and (i,j) phi-points
      C2 = half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k))*GVM%mc_iJz_8(i,j,k)*( UT1p1_phi - UT1_phi )
! v on k t-levels and (i,j) v-points
      VT1_t = half*(F_v(i,j,k)+F_v(i,j,k+1))
! v on k t-levels and (i,j-1) v-points
      VT1_t_mj = half*(F_v(i,j-1,k)+F_v(i,j-1,k+1))
! v on k t-levels and (i,j) phi-points
      VT1_t_phi = half*(VT1_t+VT1_t_mj)
! (dz/dx)*(du/dz) + v\tan\phi/a on k t-levels and (i,j) phi-points
      C2 = C2 + VT1_t_phi * geomh_tyoa_8(j) /  Dcst_rayt_8 
! 2 ( du/dx - (dz/dx)*(du/dz) - v\tan\phi/a )^2 on k t-levels and (i,j) phi-points
      k1_h = 2.D0*(C1-C2)*(C1-C2)

! dv/dy on k t-levels and (i,j) phi-points
      C1  = ( VT1_t - VT1_t_mj ) * geomh_invDY_8
! 
! v on k m-levels and (i,j) phi-points
      VT1_phi = half * (F_v(i,j,k)+F_V(i,j-1,k))
! v on k+1 m-levels and (i,j) phi-points
      VT1p1_phi = half * (F_v(i,j,k+1)+F_v(i,j-1,k+1))

! (dz/dy)*(dv/dz) on k t-levels and (i,j) phi-points
      C2 = half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k))*GVM%mc_iJz_8(i,j,k)*( VT1p1_phi - VT1_phi ) 
! 2 ( dv/dy - (dz/dy)*(dv/dz) )^2 on k t-levels and (i,j) phi-points
      k2_h = 2.D0*(C1-C2)*(C1-C2)
!
! w on k t-levels and (i,j) u-points
      WT1_u = half*(F_w(i,j,k)+F_w(i+1,j,k))
! w on k t-levels and (i-1,j) u-points
      WT1_um1 = half*(F_w(i-1,j,k)+F_w(i,j,k))
! dw/dx on k t-levels and (i,j) phi-points
      C1  = ( WT1_u - WT1_um1 ) * geomh_invDX_8(j)
! w on k m-levels and (i,j) phi-points
       WT1_m =   Ver_wp_8%m(k)*F_w(i  ,j,k)+Ver_wm_8%m(k)*F_w(i  ,j,k-1)
! w on k+1 m-levels and (i,j) phi-points
       if (k.eq.NK)  then 
          WT1p1_m= F_w(i  ,j,k+1)
       else
          WT1p1_m =  Ver_wp_8%m(k+1)*F_w(i  ,j,k+1)+Ver_wm_8%m(k+1)*F_w(i  ,j,k)
       endif    
! (dz/dx)*(dw/dz) on k t-levels and (i,j) phi-points
      C2 = half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k))*GVM%mc_iJz_8(i,j,k)*( WT1p1_m - WT1_m ) 
      k3_h = (C1-C2)*(C1-C2)

!
! w on k t-levels and (i,j) v-points
      WT1_v = half*(F_w(i,j,k)+F_w(i,j+1,k))
! w on k t-levels and (i,j-1) v-points
      WT1_vm1 = half*(F_w(i,j-1,k)+F_w(i,j,k))
! dw/dy on k t-levels and (i,j) phi-points
      C1  = ( WT1_v - WT1_vm1 ) * geomh_invDY_8
! (dz/dy)*(dw/dz) on k t-levels and (i,j) phi-points
      C2 = half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k))*GVM%mc_iJz_8(i,j,k)*( WT1p1_m - WT1_m ) 
      k4_h = (C1-C2)*(C1-C2)
!
! v on k m-level and top-right corner
      VT1_tr = half*(F_v(i,j,k)+F_v(i+1,j,k))
! v on k+1 m-level and top-right corner
      VT1_tr_p1 = half*(F_v(i,j,k+1)+F_v(i+1,j,k+1))
! v on k m-level and lower-right corner
      VT1_lr = half*(F_v(i,j-1,k)+F_v(i+1,j-1,k))
! v on k+1 m-level and lower-right corner
      VT1_lr_p1 = half*(F_v(i,j-1,k+1)+F_v(i+1,j-1,k+1))
! v on k m-level and top-left corner
      VT1_tl = half*(F_v(i,j,k)+F_v(i-1,j,k))
! v on k+1 m-level and top-left corner
      VT1_tl_p1 = half*(F_v(i,j,k+1)+F_v(i-1,j,k+1))
! v on k m-level and lower-left corner
      VT1_ll = half*(F_v(i,j-1,k)+F_v(i-1,j-1,k))
! v on k+1 m-level and lower-left corner
      VT1_ll_p1 = half*(F_v(i,j-1,k+1)+F_v(i-1,j-1,k+1))
! dv/dx on k m-level and (i,j) phi-points
      C1 = ( half*(VT1_tr+VT1_lr) - half*(VT1_tl+VT1_ll) ) * geomh_invDX_8(j)
! dv/dx on k+1 m-level and (i,j) phi-points
      C2 = ( half*(VT1_tr_p1+VT1_lr_p1) - half*(VT1_tl_p1+VT1_ll_p1) ) * geomh_invDX_8(j)
! dv/dx on k t-level and (i,j) phi-points
      C1 = half*(C1+C2)

! u on k m-level and top-right corner
      UT1_tr = half*(F_u(i,j,k)+F_u(i,j+1,k))
! u on k+1 m-level and top-right corner
      UT1_tr_p1 = half*(F_u(i,j,k+1)+F_u(i,j+1,k+1))
! u on k m-level and lower-right corner
      UT1_lr = half*(F_u(i,j-1,k)+F_u(i,j,k))
! u on k+1 m-level and lower-right corner
      UT1_lr_p1 = half*(F_u(i,j-1,k+1)+F_u(i,j,k+1))
! u on k m-level and top-left corner
      UT1_tl = half*(F_u(i-1,j,k)+F_u(i-1,j+1,k))
! u on k+1 m-level and top-left corner
      UT1_tl_p1 = half*(F_u(i-1,j,k+1)+F_u(i-1,j+1,k+1))
! u on k m-level and lower-left corner
      UT1_ll = half*(F_u(i-1,j-1,k)+F_u(i-1,j,k))
! u on k+1 m-level and lower-left corner
      UT1_ll_p1 = half*(F_u(i-1,j-1-1,k+1)+F_u(i-1,j,k+1))
! du/dy on k m-level and (i,j) phi-points
      C2 = ( half*(UT1_tr+UT1_tl) - half*(UT1_lr+UT1_ll) ) * geomh_invDY_8
! du/dy on k+1 m-level and (i,j) phi-points
      C3 = ( half*(UT1_tr_p1+UT1_tl_p1) - half*(UT1_lr_p1+UT1_ll_p1) ) * geomh_invDY_8
! du/dy on k t-level and (i,j) phi-points
      C2 = half*(C2+C3)
! (dz/dx)*(dv/dz) on k t-levels and (i,j) phi-points
      C3 = half*(GVM%mc_Jxt_8(i,j,k)+GVM%mc_Jxt_8(i-1,j,k))*GVM%mc_iJz_8(i,j,k)*( VT1p1_phi - VT1_phi )
! (dz/dy)*(du/dz) on k t-levels and (i,j) phi-points
      C4 = half*(GVM%mc_Jyt_8(i,j,k)+GVM%mc_Jyt_8(i,j-1,k))*GVM%mc_iJz_8(i,j,k)*( UT1p1_phi - UT1_phi ) 
! ( dv/dx - (dz/dx)*(dv/dz) + du/dy - (dz/dy)*(du/dz) )^2 on k t-levels and (i,j) phi-points
      k5_h = (C1+C2-C3-C4)*(C1+C2-C3-C4)     
!
      F_s(i,j,k) = ( k1_h+k2_h+k3_h+k4_h+k5_h)**half
!       F_s(i,j,k)= sqrt(k1_h+k2_h+k3_h+k4_h+k5_h)
    enddo
  enddo
enddo
!
!      F_s(:,:,nk)= F_s(:,:,nk-1)
!
      return
      end

