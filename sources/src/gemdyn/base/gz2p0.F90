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

!**s/r gz2p0 - Compute P0 from GZ from pressure coordinate
!
      subroutine gz2p0 ( F_ps, F_gz, F_topo, F_rna ,&
                         Mminx,Mmaxx,Mminy,Mmaxy,Nk,&
                         F_i0,F_in,F_j0,F_jn )
      use dyn_fisl_options
      use dynkernel_options
      use geomh
      use tdpack
      use glb_ld
      use cstv
      use ver
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Nk,Mminx,Mmaxx,Mminy,Mmaxy,F_i0,F_in,F_j0,F_jn
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk), intent(in) :: F_gz
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy), intent(out) :: F_ps
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy), intent(inout) :: F_topo
      real, dimension(Nk), intent(in) :: F_rna

!arguments
!  Name                       Description
!-----------------------------------------------------------
! F_ps         - ln(pi_s/z_s)
! F_gz         - geopotential height
! F_topo       - topography
! F_rna        - pressure levels from pressure analyse
! Nk           - number of levels from the pressure analyse


      integer i,j,k,m,NN
      real, allocatable, dimension(:  ) :: guess,a,topo
      real, allocatable, dimension(:,:) :: zcol,tcol
      real lna(Nk),sdd(Nk),conv,acc
!
!     ---------------------------------------------------------------
!
      if (Schm_autobar_L) then
         F_topo = 0.
         do j=1,l_nj
            do i=1,l_ni
               F_ps(i,j) =  (grav_8*F_gz(i,j,1)-F_topo(i,j)) &
                           /(Rgasd_8*Cstv_Tstr_8) &
                           +Ver_z_8%m(1)-Cstv_Zsrf_8
               F_ps(i,j) = exp(F_ps(i,j)) * Cstv_pref_8
            end do
         end do
         return
      end if

      if (Nk < 2) then
         call gem_error(-1,'gz2p0','Impossible to proceed: NK < 2')
      end if

      acc = .1 * grav_8
      conv = log(100.)
!     Convert millibar to log of pascal unit - Pressure Analysis
      do k=1,Nk
         lna(k) = log(F_rna(k))
      end do
      do k=1,Nk-1
         sdd(k) = 1./(lna(k+1)-lna(k))
      end do

      NN=(F_jn-F_j0+1)*(F_in-F_i0+1)
      allocate(guess(NN),a(NN),topo(NN))
      allocate(zcol(NN,Nk),tcol(NN,Nk))

      m=0
      do j=F_j0,F_jn
         do i=F_i0,F_in
            m=m+1
            topo(m)=F_topo(i,j)
         end do
      end do
      do k=1,Nk
         m=0
         do j=F_j0,F_jn
            do i=F_i0,F_in
               m=m+1
               zcol(m,k) = grav_8*F_gz(i,j,k)
            end do
         end do
      end do
!
!     Compute derivative of geopotential (from vdfds)
!
      do k=1,Nk-1
         do i=1,NN
            tcol(i,k+1) = sdd(k)*(zcol(i,k+1)-zcol(i,k))
         end do
      end do

      do i=1,NN
        a(i) = tcol(i,2)
      end do

      do k=2,Nk-1
         do i=1,NN
            tcol(i,k) = (sdd(k)*tcol(i,k+1)+sdd(k-1)*tcol(i,k)) &
                        /(sdd(k)+sdd(k-1))
         end do
      end do

!     BOUNDARIES
      do i=1,NN
         tcol(i,1)  = a(i)
         tcol(i,Nk) = tcol(i,Nk)
      end do

!     Compute pressure at the surface (PS)
      do i=1,NN
         guess(i) = lna(Nk)-topo(i)/(rgasd_8*250.)
      end do

      call vterp (guess,topo,zcol,tcol,lna,acc,NN,Nk)
      m=0
      do j=F_j0,F_jn
         do i=F_i0,F_in
            m=m+1
            F_ps(i,j) = exp (guess(m) + conv)
         end do
      end do

      deallocate(guess,a,topo,zcol,tcol)
!
!     ---------------------------------------------------------------
!
      return
      end
