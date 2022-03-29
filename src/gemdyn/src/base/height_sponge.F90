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

!**s/r  pressure_sponge -  Performs vertical blending
!
      subroutine height_sponge()
      use gmm_vt0
      use gmm_geof
      use gem_options
      use theo_options
      use glb_ld
      use dcst
      implicit none
#include <arch_specific.hf>

      real betav_m(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
           betav_t(l_minx:l_maxx,l_miny:l_maxy,l_nk),ubar

!----------------------------------------------------------------------

      if (Dynamics_hauteur_L) then
         call height_spongeH ()
         return
      endif
      
!$omp single
      call set_betav (betav_m, betav_t, st0, sls, l_minx, l_maxx,&
                      l_miny, l_maxy, l_nk)

      ubar= mtn_flo

      if(Theo_case_S /= 'MTN_SCHAR' ) then
         call apply (ut0, ubar, betav_m, l_minx,l_maxx,l_miny,l_maxy, l_nk)
      end if

      call apply (wt0, 0., betav_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)
!$omp end single

!----------------------------------------------------------------------
      return
      end

!=======================================================================


      subroutine apply(ff,valu,betav, Minx,Maxx,Miny,Maxy, Nk)
      use glb_ld
      implicit none
#include <arch_specific.hf>

      integer :: Minx,Maxx,Miny,Maxy, Nk

      real ff(Minx:Maxx,Miny:Maxy,Nk),valu,betav(Minx:Maxx,Miny:Maxy,Nk)

      integer i,j,k,i0,in,j0,jn

      i0 = 1
      in = l_ni
      j0 = 1
      jn = l_nj
      if (l_west ) i0 = 1+pil_w
      if (l_east ) in = l_ni-pil_e
      if (l_south) j0 = 1+pil_s
      if (l_north) jn = l_nj-pil_n

      do k=1,Nk
         do j=j0,jn
            do i=i0,in
               ff(i,j,k)=(1.-betav(i,j,k))*ff(i,j,k)+betav(i,j,k)*valu
            end do
         end do
      end do

      return

      end
!=======================================================================
