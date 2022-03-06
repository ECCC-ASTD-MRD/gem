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
!
!*s/r fislh_rhs - compute the right-hand sides

      subroutine fislh_rhs ( F_dt_8 )
      use HORgrid_options
      use gem_options
      use dyn_fisl_options
      use dynkernel_options
      use coriolis
      use geomh
      use gmm_geof
      use tdpack
      use gmm_phy
      use glb_ld
      use gmm_contiguous
      use gmm_vt0
      use gmm_vt1
      use adz_mem
      use mem_tstp
      use mem_tracers
      use cstv
      use ver
      use metric
      use, intrinsic :: iso_fortran_env
      implicit none
      
      real(kind=REAL64), intent(IN) :: F_dt_8
      
!Author: Claude Girard, July 2017
!        Syed Husain, Jun 2019 (new dyn-phy coupling and
!                               new elimination of variables)

      integer :: HLT_j0, HLT_jn, HLT_nj, HLT_nk, &
                 HLT_np, HLT_start, HLT_end
      integer :: i, j, k, km, kp, n
      real, dimension(:,:,:), pointer :: logT, logP
      real(kind=REAL64) :: div, barz, barzp, u_interp, v_interp,&
               t_interp, w2, w3, w4, invT_8, invT_nh_8,invT_m_8
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
!
!     ---------------------------------------------------------------
!
!$omp do
      do n= 1, ubound(dynt0,1)
         dynt0(n) = dynt1(n)
      end do
!$omp end do nowait
!$omp do
      do n= 1, ubound(trt0,1)
         trt0(n) = trt1(n)
      end do
!$omp end do nowait
      
      logT(1:l_ni,1:l_nj,1:l_nk) => WS1(1:)
      logP(1:l_ni,1:l_nj,1:l_nk) => WS1(l_ni*l_nj*l_nk+1:)
!$omp do collapse(2)
      do k=1, l_nk
         do j=1, l_nj
         do i= 1, l_ni
            logT(i,j,k)= log (tt1(i,j,k)/Cstv_Tstr_8)
            logP(i,j,k)= log ((qt1(i,j,k) - fis0(i,j))*Cstv_invFI_8 + 1.d0)
        end do
      end do
      end do
!$omp enddo

      HLT_j0 = 1
      HLT_jn = l_nj
      HLT_nk = l_nk
      HLT_nj = HLT_jn - HLT_j0 + 1
      call HLT_split (1, HLT_nj*HLT_nk, HLT_np, HLT_start, HLT_end)

!**********************************************
! Compute Ru, Rv : RHS of U, V equations      *
! Compute Rw & Rt: RHS of w & T equations     *
! Compute Rc     : RHS of continuity equation *
!**********************************************

      w3=Cstv_bar1_8*half/(cpd_8*Cstv_Tstr_8)
      w4=Cstv_bar1_8*epsi_8/grav_8
      invT_8   = one/F_dt_8/Cstv_bA_8
      invT_m_8 = one/F_dt_8/Cstv_bA_m_8
      invT_nh_8= one/F_dt_8/Cstv_bA_nh_8    

      do n= HLT_start, HLT_end
         k= (n-1)/HLT_nj
         j= n - k*HLT_nj + HLT_j0 - 1
         k= k+1
         
         km=max(k-1,1)
         kp=min(k+1,l_nk)

         do i= 1, l_ni

            !assuming T at momentum level 1 equal to T at thermo level 3/2
            barz  = Ver_wp_8%m(k)*tt1(i  ,j,k)+Ver_wm_8%m(k)*tt1(i  ,j,km)
            barzp = Ver_wp_8%m(k)*tt1(i+1,j,k)+Ver_wm_8%m(k)*tt1(i+1,j,km)
            t_interp = (barz + barzp)*half/Cstv_Tstr_8
            v_interp = 0.25d0*(vt1(i,j,k)+vt1(i,j-1,k)+vt1(i+1,j,k)+vt1(i+1,j-1,k))
            orhsu_ext(i,j,k) = invT_m_8  * ut1(i,j,k) - Cstv_Beta_m_8 * ( &
                           t_interp * ( qt1(i+1,j,k) - qt1(i,j,k) ) * geomh_invDX_8(j)        &
                         - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * ut1(i,j,k) ) * v_interp )  &
                         + Cstv_Beta_m_8 * t_interp * mc_Jx_8(i,j,k) * ( &
                           Ver_wp_8%m(k)*half*( (qt1(i+1,j,k+1)-qt1(i+1,j,k ))*mc_iJz_8(i+1,j,k )   &
                                              + (qt1(i  ,j,k+1)-qt1(i  ,j,k ))*mc_iJz_8(i  ,j,k ) ) &
                         + Ver_wm_8%m(k)*half*( (qt1(i+1,j,k  )-qt1(i+1,j,km))*mc_iJz_8(i+1,j,km)   &
                                              + (qt1(i  ,j,k  )-qt1(i  ,j,km))*mc_iJz_8(i  ,j,km) ) )&
                         + ((1d0-phy_cplm(i,j))/Cstv_bA_m_8) * rhs_uu_tend(i,j,k)

            barz  = Ver_wp_8%m(k)*tt1(i,j  ,k)+Ver_wm_8%m(k)*tt1(i,j  ,km)
            barzp = Ver_wp_8%m(k)*tt1(i,j+1,k)+Ver_wm_8%m(k)*tt1(i,j+1,km)
            t_interp = ( barz + barzp)*half/Cstv_Tstr_8
            u_interp = 0.25d0*(ut1(i,j,k)+ut1(i-1,j,k)+ut1(i,j+1,k)+ut1(i-1,j+1,k))
            orhsv_ext(i,j,k) = invT_m_8  * vt1(i,j,k) - Cstv_Beta_m_8 * ( &
                           t_interp * ( qt1(i,j+1,k) - qt1(i,j,k) ) * geomh_invDY_8           &
                         + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp )   &
                         + Cstv_Beta_m_8 * t_interp * mc_Jy_8(i,j,k) * ( &
                           Ver_wp_8%m(k)*half*( (qt1(i,j+1,k+1)-qt1(i,j+1,k ))*mc_iJz_8(i,j+1,k )   &
                                              + (qt1(i,j  ,k+1)-qt1(i,j  ,k ))*mc_iJz_8(i,j  ,k ) ) &
                         + Ver_wm_8%m(k)*half*( (qt1(i,j+1,k  )-qt1(i,j+1,km))*mc_iJz_8(i,j+1,km)   &
                                              + (qt1(i,j  ,k  )-qt1(i,j  ,km))*mc_iJz_8(i,j  ,km) ) )&
                         + ((1d0-phy_cplm(i,j))/Cstv_bA_m_8) * rhs_vv_tend(i,j,k)

            orhst_ext(i,j,k) = invT_8 * ( logT(i,j,k) - w3*(qt1(i,j,k+1)+qt1(i,j,k)) ) &
                             - Cstv_Beta_8 * mu_8 * wt1(i,j,k) &
               + rhs_phytv*((1d0-phy_cplt(i,j))/Cstv_bA_8) * 1./tt1(i,j,k) * phy_tv_tend(i,j,k)

            orhsf_ext(i,j,k) = invT_nh_8 * (ztht_8(i,j,k)-Ver_z_8%t(k)) * Cstv_bar1_8 &
                         - Cstv_Beta_nh_8 * ( Ver_wpstar_8(k)*zdt1(i,j,k)+Ver_wmstar_8(k)*zdt1(i,j,km) - wt1(i,j,k) )
            div = (ut1 (i,j,k)-ut1 (i-1,j,k))*geomh_invDXM_8(j)     &
                   + (vt1 (i,j,k)*geomh_cyM_8(j)-vt1 (i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j) &
                   + (zdt1(i,j,k)*Ver_wpstar_8(k)+(Ver_wmstar_8(k)-Ver_onezero(k))*zdt1(i,j,km))*Ver_idz_8%m(k) &
                   + half * ( mc_Ix_8(i,j,k)*(ut1(i,j,k)+ut1(i-1,j,k))        &
                   + mc_Iy_8(i,j,k)*(vt1(i,j,k)+vt1(i,j-1,k)) )      &
                   + mc_Iz_8(i,j,k)*(Ver_wp_8%m(k)*zdt1(i,j,k)+Ver_wm_8%m(k)*Ver_onezero(k)*zdt1(i,j,km))

            orhsc_ext (i,j,k) = invT_8 *  w4 * qt1(i,j,k) +   invT_8 * mc_logJz_8(i,j,k)  &
                          - Cstv_Beta_8 * ( div-epsi_8*(Ver_wp_8%m(k)*wt1(i,j,k)+Ver_onezero(k)*Ver_wm_8%m(k)*wt1(i,j,km)) ) &
                          + (1d0/Cstv_bA_8) * &
                          (              Ver_wp_8%m(k)*phy_tv_tend(i,j,k )/tt1(i,j,k ) + &
                          Ver_onezero(k)*Ver_wm_8%m(k)*phy_tv_tend(i,j,km)/tt1(i,j,km) ) &
                          + (1.0d0-Cstv_bar1_8) * invT_8 * logP(i,j,k)
 
        end do
      end do
      
      if (.not.Dynamics_hydro_L) then
         do n= HLT_start, HLT_end
            k= (n-1)/HLT_nj
            j= n - k*HLT_nj + HLT_j0 - 1
            k= k+1
            do i= 1, l_ni
               w2 = tt1(i,j,k)/Cstv_Tstr_8
               orhsw_ext(i,j,k) = invT_nh_8 * wt1(i,j,k) - Cstv_Beta_nh_8 *w2* &
                              ( (qt1(i,j,k+1)-qt1(i,j,k))*mc_iJz_8(i,j,k)  -grav_8*(one-one/w2) )
            end do
         end do
      endif

!$OMP BARRIER
!     
!     ---------------------------------------------------------------
!
      return
      end subroutine fislh_rhs
