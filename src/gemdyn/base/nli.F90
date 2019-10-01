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
!**s/r nlip   - compute non-linear terms:  Nu, Nv, Nt, Nc, Nw, Nf,
!             - compute full right-hand side of Helmholtz eqn: Rp=Rc-Nc
!
!**********************************************************************
!
      subroutine nli ( F_nu , F_nv , F_nt   , F_nc , F_nw , F_nf  , &
                       F_u  , F_v  , F_t    , F_s  , F_zd , F_q   , &
                       F_rhs, F_rc , F_sl   , F_fis, F_nb , F_hu  , &
                       Minx, Maxx, Miny, Maxy, Nk , ni, nj,         &
                       i0, j0, in, jn, k0, icln )
      use coriolis
      use cstv
      use dcst
      use gem_options
      use dynkernel_options
      use dyn_fisl_options
      use geomh
      use glb_ld
      use HORgrid_options
      use lun
      use tdpack
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk, ni, nj, i0, j0, in, jn, k0, icln

      real, dimension(Minx:Maxx,Miny:Maxy), intent(inout) :: F_s
      real, dimension(Minx:Maxx,Miny:Maxy), intent(in)    :: F_fis, F_sl
      real, dimension(Minx:Maxx,Miny:Maxy), intent(out)   :: F_nb
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout)   :: F_u, F_v, F_t
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(in)      :: F_zd, F_hu, F_rc
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(out)     :: F_nu, F_nv, F_nt, F_nc, F_nw, F_nf
      real, dimension(Minx:Maxx,Miny:Maxy,Nk+1), intent(inout) :: F_q
      real(kind=REAL64), dimension(ni,nj,Nk), intent(out) :: F_rhs

!

      integer :: i, j, k, km, i0u, inu, j0v, jnv, nij, k0t, onept
      real    :: w_nt
      real(kind=REAL64)  :: c1,qbar,ndiv,w1,w2,w3,w4,barz,barzp,MUlin, &
                 t_interp, mu_interp, u_interp, v_interp
      real(kind=REAL64) , dimension(i0:in,j0:jn) :: xtmp_8, ytmp_8
      real(kind=REAL64), parameter :: one=1.d0, half=0.5d0
      real, dimension(Minx:Maxx,Miny:Maxy,l_nk) :: MU
      real, dimension(Minx:Maxx,Miny:Maxy,l_nk+1) :: BsPq, BsPrq, FI
!     __________________________________________________________________
!
      if (Lun_debug_L) then
         write(Lun_out,1000)
      end if

      if (icln > 1) then
         call rpn_comm_xch_halo( F_u   ,l_minx,l_maxx,l_miny,l_maxy,l_niu,l_nj ,G_nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         call rpn_comm_xch_halo( F_v   ,l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_njv,G_nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         call rpn_comm_xch_halo( F_t   ,l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,G_nk, &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         call rpn_comm_xch_halo( F_s   ,l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,1   , &
                                 G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         if (.not.Dynamics_hydro_L) then
            call rpn_comm_xch_halo( F_q,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk+1, &
                                    G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         end if
      end if

      c1 = Dcst_rayt_8**2

      k0t=k0
      if (Schm_opentop_L) k0t=k0-1
      nij = (in - i0 +1)*(jn - j0 +1)

      onept= 0
      if (Grd_yinyang_L) onept=1
!
!***********************************************************
! The nonlinear deviation of horizontal momentum equations *
!***********************************************************
!
!     Indices

      i0u = i0-1
      j0v = j0-1
      inu = l_niu-pil_e
      jnv = l_njv-pil_n

      if (l_west ) i0u=i0 -onept
      if (l_south) j0v=j0 -onept
      if (l_east ) inu=inu+onept
      if (l_north) jnv=jnv+onept

      call diag_fip( FI, F_s, F_sl, F_t, F_q, F_fis, l_minx,l_maxx,l_miny,l_maxy,&
                                               l_nk, i0u, inu+1, j0v, jnv+1 )

      if (Dynamics_hydro_L) then
         MU = 0.
      else
         call diag_mu ( MU, F_q, F_s, F_sl, l_minx,l_maxx,l_miny,l_maxy, l_nk,&
                                                 i0u, inu+1, j0v, jnv+1 )
      end if


      do k=1,l_nk+1
          do j=j0v,jnv+1
             do i=i0u,inu+1
                BsPq(i,j,k)  = Ver_b_8%m(k) *F_s(i,j)  &
                             + Ver_c_8%m(k) *F_sl(i,j) + F_q(i,j,k)
                BsPrq(i,j,k) = Ver_b_8%m(k) *F_s(i,j)  &
                             + Ver_c_8%m(k) *F_sl(i,j) + F_q(i,j,k)

             end do
          end do
      end do

      do k=k0,l_nk
      km=max(k-1,1)

!     Compute Nu
!     ~~~~~~~~~~

!     V barY stored in wk2
!     ~~~~~~~~~~~~~~~~~~~~

      do j = j0, jn
         do i = i0u, inu

   !        mu barXZ
   !        ~~~~~~~~
            barz  = Ver_wpM_8(k)*MU(i  ,j,k)+Ver_wmM_8(k)*MU(i  ,j,km)
            barzp = Ver_wpM_8(k)*MU(i+1,j,k)+Ver_wmM_8(k)*MU(i+1,j,km)
            mu_interp = (barz+barzp)*half

   !        Pressure gradient and mu terms: RT barXZ * dBsPq/dX + mu barXZ * dfi'/dX
   !                                        - RTstr barZ * dBsPrq/dX
   !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            barz  = Ver_wpM_8(k)*(F_t(i  ,j,k)-Cstv_Tstr_8)+Ver_wmM_8(k)*(F_t(i  ,j,km)-Cstv_Tstr_8)
            barzp = Ver_wpM_8(k)*(F_t(i+1,j,k)-Cstv_Tstr_8)+Ver_wmM_8(k)*(F_t(i+1,j,km)-Cstv_Tstr_8)
            t_interp = (barz+barzp)*half

            barz  = Ver_wpM_8(k)*Cstv_Tstr_8+Ver_wmM_8(k)*Cstv_Tstr_8

            w1 = ( BsPq(i+1,j,k)  - BsPq(i,j,k)  ) * geomh_invDXMu_8(j)
            w2 = (   FI(i+1,j,k)  -   FI(i,j,k)  ) * geomh_invDXMu_8(j)

            F_nu(i,j,k) = rgasd_8 * t_interp * w1 + mu_interp * w2

   !        Coriolis term & metric terms: - (f + tan(phi)/a * U ) * V barXY
   !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            v_interp = 0.25d0*(F_v(i,j,k)+F_v(i,j-1,k)+F_v(i+1,j,k)+F_v(i+1,j-1,k))

            F_nu(i,j,k) = F_nu(i,j,k) - ( Cori_fcoru_8(i,j) + geomh_tyoa_8(j) * F_u(i,j,k) ) * v_interp

         end do
      end do

!     Compute Nv
!     ~~~~~~~~~~
      do j = j0v, jnv
         do i = i0, in

   !        mu barYZ
   !        ~~~~~~~~
            barz  = Ver_wpM_8(k)*MU(i,j  ,k)+Ver_wmM_8(k)*MU(i,j  ,km)
            barzp = Ver_wpM_8(k)*MU(i,j+1,k)+Ver_wmM_8(k)*MU(i,j+1,km)
            mu_interp = (barz+barzp)*half

   !        Pressure gradient and Mu term: RT' barYZ * dBsPq/dY + mu barYZ * dfi'/dY
   !                                     - RTstr barZ * dBsPrq/dY
   !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            barz  = Ver_wpM_8(k)*(F_t(i,j  ,k)-Cstv_Tstr_8)+Ver_wmM_8(k)*(F_t(i,j  ,km)-Cstv_Tstr_8)
            barzp = Ver_wpM_8(k)*(F_t(i,j+1,k)-Cstv_Tstr_8)+Ver_wmM_8(k)*(F_t(i,j+1,km)-Cstv_Tstr_8)
            t_interp = (barz+barzp)*half

            barz  = Ver_wpM_8(k)*Cstv_Tstr_8+Ver_wmM_8(k)*Cstv_Tstr_8

            w1 = (  BsPq(i,j+1,k) -  BsPq(i,j,k) ) * geomh_invDYMv_8(j)
            w2 = (    FI(i,j+1,k) -    FI(i,j,k) ) * geomh_invDYMv_8(j)

            F_nv(i,j,k) = rgasd_8 * t_interp * w1 + mu_interp * w2

   !        Coriolis term & metric terms: + f * U barXY + tan(phi)/a * (U barXY)^2
   !        ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            u_interp = 0.25d0*(F_u(i,j,k)+F_u(i-1,j,k)+F_u(i,j+1,k)+F_u(i-1,j+1,k))

            F_nv(i,j,k) = F_nv(i,j,k) + ( Cori_fcorv_8(i,j) + geomh_tyoav_8(j) * u_interp ) * u_interp

         end do
      end do

      end do


!     Set  Nu=0  on the east and west boundaries of the LAM grid
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      if (.not.Grd_yinyang_L) then

         if (l_west) then

            do k=1,l_nk
               do j=j0,jn
                  F_nu(pil_w,j,k) = 0.
               end do
            end do

         end if
         if (l_east) then

            do k=1,l_nk
               do j=j0,jn
                  F_nu(l_ni-pil_e,j,k) = 0.
               end do
            end do

         end if

!     Set  Nv=0  on the north and south boundaries  of the LAM grid
!     and        at the north and south poles       of the GLOBAL grid
!     ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

         if (l_south) then

            do k=1,l_nk
               do i=i0,in
                  F_nv(i,pil_s,k) = 0.
               end do
            end do

         end if
         if (l_north) then

            do k=1,l_nk
               do i=i0,in
                  F_nv(i,l_nj-pil_n,k) = 0.
               end do
            end do

         end if

      end if


      do k=k0t,l_nk
         km=max(k-1,1)
!**********************************
!   The nonlinear deviation of    *
! the thermodynamic equation: Nt' *
!**********************************

!        Compute Nw and Nt' (Nf=0)
!        ~~~~~~~~~~~~~~~~~~~~~~~~~
         w1 = one / Cstv_Tstr_8
         do j= j0, jn
            do i= i0, in
               xtmp_8(i,j) = F_t(i,j,k) * w1
            end do
         end do
         call vlog ( ytmp_8, xtmp_8, nij )
         if(Schm_opentop_L.and.k == k0t) then
            do j= j0, jn
               do i= i0, in
                  F_nb(i,j) = Cstv_invT_8*(ytmp_8(i,j)-xtmp_8(i,j)+one)
               end do
            end do
         end if
         w1 = Ver_idz_8%t(k) / Rgasd_8 / Cstv_Tstr_8
         w2 = one / Cstv_Tstr_8 * Ver_idz_8%t(k) / Rgasd_8
         do j= j0, jn
            do i= i0, in
               w4=Ver_wpstar_8(k)*F_zd(i,j,k)+Ver_wmstar_8(k)*F_zd(i,j,km)
               qbar=Ver_wpstar_8(k)*F_q(i,j,k+1)+Ver_wmstar_8(k)*half*(F_q(i,j,k)+F_q(i,j,km))
               qbar=Ver_wp_8%t(k)*qbar+Ver_wm_8%t(k)*F_q(i,j,k)
               MUlin=Ver_idz_8%t(k)*(F_q(i,j,k+1)-F_q(i,j,k)) + qbar
               F_nw(i,j,k) = - grav_8 * ( MU(i,j,k) - MUlin )
               F_nt(i,j,k) = Cstv_invT_8*(ytmp_8(i,j)  &
                                    + w2*( FI(i,j,k+1)+Rgasd_8*Cstv_Tstr_8*BsPrq(i,j,k+1) &
                                     - FI(i,j,k  )-Rgasd_8*Cstv_Tstr_8*BsPrq(i,j,k  ) ) &
                                     - MUlin )
               F_nf(i,j,k) = 0.0
            end do
         end do

!        Compute Nt" and Nf"
!        ~~~~~~~~~~~~~~~~~~~
         w1 = cappa_8/ ( Rgasd_8 * Cstv_Tstr_8 )
         w2 = Cstv_invT_m_8 / ( cappa_8 + epsi_8 )
         do j= j0, jn
            do i= i0, in
               w_nt = F_nt(i,j,k) + Ver_igt_8 * F_nw(i,j,k)
               F_nt(i,j,k) = w2 * ( w_nt + Ver_igt2_8 * F_nf(i,j,k) )
               F_nf(i,j,k) = w2 * ( w_nt - w1 * F_nf(i,j,k) )
            end do
         end do

      end do

!***************************************
!     The nonlinear deviation of       *
!   the continuity equation: Nc and    *
! the horizontal Divergence of (Nu,Nv) *
!   combined with Nc (stored in Nc)    *
!***************************************


      do k=k0,l_nk

!        Compute Nc
!        ~~~~~~~~~~
         km=max(k-1,1)
         do j = j0, jn
            do i = i0, in
               xtmp_8(i,j) = one + Ver_dbdz_8%m(k)*F_s(i,j)  &
                                 + Ver_dcdz_8%m(k)*F_sl(i,j)
            end do
         end do
         call vlog(ytmp_8, xtmp_8, nij)
         do j = j0, jn
            do i = i0, in
               F_nc(i,j,k) = Cstv_invT_8 * ( ytmp_8(i,j) +  &
                             ( Cstv_bar1_8*(Ver_b_8%m(k)-Ver_bzz_8(k)) - Ver_dbdz_8%m(k) )*F_s(i,j)  + &
                             ( Cstv_bar1_8*(Ver_c_8%m(k)-Ver_czz_8(k)) - Ver_dcdz_8%m(k) )*F_sl(i,j) ) &
                                      + (Ver_wpC_8(k)-Ver_wp_8%m(k)) * F_zd(i,j,k) &
                                      + (Ver_wmC_8(k)-Ver_wm_8%m(k)) * Ver_onezero(k) * F_zd(i,j,km)
            end do
         end do

      end do

      do k=k0,l_nk

!        Compute Nc"
!        ~~~~~~~~~~~
         km=max(k-1,1)
         w1=Ver_igt_8*Ver_wpA_8(k)
         w2=Ver_igt_8*Ver_wmA_8(k)*Ver_onezero(k)
         do j = j0, jn
            do i = i0, in
               ndiv = (F_nu(i,j,k)-F_nu(i-1,j,k)) * geomh_invDXM_8(j) &
                  + (F_nv(i,j,k)*geomh_cyM_8(j)-F_nv(i,j-1,k)*geomh_cyM_8(j-1))*geomh_invDYM_8(j)
               F_nc(i,j,k) = ndiv  - Cstv_invT_m_8 * ( F_nc(i,j,k) - w1*F_nw(i,j,k) - w2*F_nw(i,j,km) )
            end do
         end do

      end do


!**********************************************************
! The full contributions to the RHS of Helmholtz equation *
!**********************************************************

!     Finish computations of NP (combining Nc", Nt", Nf")
!     Substract NP from RP(Rc") and store result(RP-NP) in RP


      do j = j0, jn
         do i = i0, in
            F_nt(i,j,l_nk) = (F_nt(i,j,l_nk) - Ver_wmstar_8(l_nk)*F_nt(i,j,l_nk-1)) &
                             /Ver_wpstar_8(l_nk)
         end do
      end do

      do k=k0,l_nk
         km=max(k-1,1)
         w1=(Ver_idz_8%m(k) + Ver_wp_8%m(k))
         w2=(Ver_idz_8%m(k) - Ver_wm_8%m(k))*Ver_onezero(k)
         w3=Ver_wpA_8(k)*epsi_8
         w4=Ver_wmA_8(k)*epsi_8*Ver_onezero(k)
         do j = j0, jn
            do i = i0, in
               F_rhs(i,j,k) =  c1 * ( F_rc(i,j,k) - F_nc(i,j,k) &
                                     + w1 * F_nt(i,j,k) - w2 * F_nt(i,j,km)  &
                                     + w3 * F_nf(i,j,k) + w4 * F_nf(i,j,km)  )
            end do
         end do
      end do


!     Apply boundary conditions
!     ~~~~~~~~~~~~~~~~~~~~~~~~~

      if(Schm_opentop_L) then
         F_rhs(:,:,1:k0t) = 0.0

         do j= j0, jn
            do i= i0, in
               F_nb(i,j)    = F_nt(i,j,k0t)-Ver_ikt_8*F_nb(i,j)
               F_rhs(i,j,k0)= F_rhs(i,j,k0) + c1 * Ver_cstp_8 * F_nb(i,j)
            end do
         end do

      end if

      do j= j0, jn
         do i= i0, in
             F_nt(i,j,l_nk) =  Ver_wpstar_8(l_nk) * F_nt(i,j,l_nk)
             F_rhs(i,j,l_nk) = F_rhs(i,j,l_nk) - c1 * Ver_cssp_8 * F_nt(i,j,l_nk)
         end do
      end do


1000 format(/,5X,'COMPUTE NON-LINEAR RHS: (S/R NLI)')
!     __________________________________________________________________
!
      return
      end

