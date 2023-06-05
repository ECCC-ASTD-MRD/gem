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

!**s/r nest_HOR_gwa

      subroutine nest_HOR_gwa()
      use dynkernel_options
      use lam_options
      use gmm_vt1
      use mem_nest
      use mem_tracers
      use glb_ld
      use tr3d
      implicit none
#include <arch_specific.hf>

      integer i,j,k
      logical :: using_qt1
      integer :: n,deb
!
!----------------------------------------------------------------------
!
      if ( (Lam_blend_Hx <= 0).and.(Lam_blend_Hy <= 0) ) goto 999

      using_qt1 = ( .not. Dynamics_hydro_L ) .or. Dynamics_hauteur_L

      do k=1,G_nk
      do j=1,l_nj
      do i=1,l_niu
         ut1(i,j,k) = ut1(i,j,k)*(1.-nest_weightu(i,j,k)) + nest_u(i,j,k)*nest_weightu(i,j,k)
      enddo
      enddo

      do j=1,l_njv
      do i=1,l_ni
         vt1(i,j,k) = vt1(i,j,k)*(1.-nest_weightv(i,j,k)) + nest_v(i,j,k)*nest_weightv(i,j,k)
      enddo
      enddo
      do j=1,l_nj
      do i=1,l_ni
         tt1(i,j,k) =  tt1(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_t (i,j,k)*nest_weightm(i,j,k)
         wt1(i,j,k) =  wt1(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_w (i,j,k)*nest_weightm(i,j,k)
         zdt1(i,j,k)= zdt1(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_zd(i,j,k)*nest_weightm(i,j,k)
      enddo
      enddo
      enddo

      if ( using_qt1 ) then
      do k=1,G_nk+1
      do j=1,l_nj
      do i=1,l_ni
         qt1(i,j,k) = qt1(i,j,k)*(1.-nest_weightm(i,j,k)) + nest_q(i,j,k)*nest_weightm(i,j,k)
      enddo
      enddo
      enddo
      endif

      do j=1,l_nj
      do i=1,l_ni
        st1(i,j)= st1(i,j)*(1.-nest_weightm(i,j,G_nk+1)) + nest_s(i,j)*nest_weightm(i,j,G_nk+1)
      enddo
      enddo

      do n=1,Tr3d_ntr
         deb= (n-1)*G_nk + 1
         tracers_P(n)%pntr(1:l_ni,1:l_nj,1:G_nk) = &
         tracers_P(n)%pntr(1:l_ni,1:l_nj,1:G_nk)*(1.-nest_weightm(1:l_ni,1:l_nj,1:G_nk)) + &
                   nest_tr(1:l_ni,1:l_nj,deb:deb+G_nk-1)*nest_weightm(1:l_ni,1:l_nj,1:G_nk)
      enddo
      
 999  call spn_main()
!
!----------------------------------------------------------------------
!
      return
      end
