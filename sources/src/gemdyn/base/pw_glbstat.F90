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

!**s/r pw_glbstat - Global statistics on physical (PW) quantities
!
      subroutine pw_glbstat (name_S)
      use step_options
      use gmm_pw
      use gem_options
      use glb_ld
      use tr3d
      use gmm_itf_mod
      implicit none

      character(len=*) name_S
#include <arch_specific.hf>

      type(gmm_metadata) :: meta
      integer istat,fld,i,err
      integer, parameter :: nvar=17
      character(len=GMM_MAXNAMELENGTH) laliste(nvar)
      real, dimension(:,:  ), pointer :: wk2d
      real, dimension(:,:,:), pointer :: wk3d
!     ________________________________________________________________
!
      laliste( 1) = gmmk_pw_uu_moins_s
      laliste( 2) = gmmk_pw_uu_plus_s
      laliste( 3) = gmmk_pw_vv_moins_s
      laliste( 4) = gmmk_pw_vv_plus_s
      laliste( 5) = gmmk_pw_wz_plus_s
      laliste( 6) = gmmk_pw_tt_moins_s
      laliste( 7) = gmmk_pw_tt_plus_s
      laliste( 8) = gmmk_pw_pm_moins_s
      laliste( 9) = gmmk_pw_pm_plus_s
      laliste(10) = gmmk_pw_pt_moins_s
      laliste(11) = gmmk_pw_pt_plus_s
      laliste(12) = gmmk_pw_gz_moins_s
      laliste(13) = gmmk_pw_gz_plus_s
      laliste(14) = gmmk_pw_me_moins_s
      laliste(15) = gmmk_pw_me_plus_s
      laliste(16) = gmmk_pw_p0_moins_s
      laliste(17) = gmmk_pw_p0_plus_s

      do fld=1,nvar
         err = gmm_getmeta (trim(laliste(fld)), meta)
         if (meta%l(3)%high <= 1) then
            nullify(wk2d)
            istat= gmm_get(laliste(fld),wk2d)
            call glbstat2 (wk2d, laliste(fld)(4:), name_S,&
                       l_minx,l_maxx, l_miny,l_maxy, 1,1,&
                       1,G_ni,1,G_nj,1, 1)
         else
            nullify(wk3d)
            istat= gmm_get(laliste(fld),wk3d)
            call glbstat2 (wk3d, laliste(fld)(4:), name_S,&
                       l_minx,l_maxx, l_miny,l_maxy, 1,meta%l(3)%high,&
                       1,G_ni,1,G_nj,1, meta%l(3)%high)
         end if
      end do

      do i=1,Tr3d_ntr
         nullify(wk3d)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(i))//':M',wk3d)
         call glbstat2 (wk3d,trim(Tr3d_name_S(i))//':M', name_S,&
                        l_minx,l_maxx, l_miny,l_maxy, 1,G_nk,&
                        1,G_ni,1,G_nj,1, G_nk)
         nullify(wk3d)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(i))//':P',wk3d)
         call glbstat2 (wk3d,trim(Tr3d_name_S(i))//':P', name_S,&
                        l_minx,l_maxx, l_miny,l_maxy, 1,G_nk,&
                        1,G_ni,1,G_nj,1, G_nk)
      end do
!     ________________________________________________________________
!
      return
      end
