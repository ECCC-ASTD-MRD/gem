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
      use gmm_vt1
      use gmm_geof
      use gem_options
      use theo_options
      use gmm_itf_mod
      use glb_ld
      use dcst
      implicit none
#include <arch_specific.hf>
!author
!     Plante A.           - May 2004
!
!revision
!

      type(gmm_metadata) :: mymeta
      integer :: istat
      real betav_m(l_minx:l_maxx,l_miny:l_maxy,l_nk),betav_t(l_minx:l_maxx,l_miny:l_maxy,l_nk),ubar

!----------------------------------------------------------------------

      istat = gmm_get(gmmk_ut1_s,ut1,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(ut1)'
      istat = gmm_get(gmmk_vt1_s,vt1,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(vt1)'
      istat = gmm_get(gmmk_wt1_s,wt1,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(wt1)'
      istat = gmm_get(gmmk_tt1_s,tt1,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(tt1)'
      istat = gmm_get(gmmk_st1_s,st1,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(st1)'
      istat = gmm_get(gmmk_qt1_s,qt1,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(qt1)'
      istat = gmm_get(gmmk_sls_s,sls,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(sls)'
      istat = gmm_get(gmmk_fis0_s,fis0,mymeta)
      if (GMM_IS_ERROR(istat)) print *,'height_sponge ERROR at gmm_get(fis0)'

      call set_betav(betav_m,betav_t,st1,sls,fis0, &
                     l_minx,l_maxx,l_miny,l_maxy,l_nk)

      ubar= mtn_flo

      if(Theo_case_S /= 'MTN_SCHAR' ) then
         call apply (ut1, ubar, betav_m, l_minx,l_maxx,l_miny,l_maxy, l_nk)
      end if

      call apply (wt1, 0., betav_t, l_minx,l_maxx,l_miny,l_maxy, l_nk)

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
