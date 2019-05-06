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

!**s/r pw_update - Update physical quantities WZ, GZ, PM and PT
!
      subroutine fislh_pw_update_GPW(F_q)
      use dynkernel_options
      use metric
      use gem_timing
      use gmm_vt1
      use gmm_vt0
      use gmm_pw
      use gmm_geof
      use glb_ld
      use ver
      use gmm_itf_mod
      use tdpack
      implicit none
#include <arch_specific.hf>

      real, dimension(l_minx:l_maxx,l_miny:l_maxy), intent(in) :: F_q

      integer :: k, istat
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: delps, pw_pm
!     ________________________________________________________________
!
      call gemtime_start ( 5, 'PW_UPDATE_GPW', 0)

      istat = gmm_get(gmmk_pw_wz_plus_s , pw_wz_plus )
      istat = gmm_get(gmmk_pw_gz_plus_s , pw_gz_plus )
      istat = gmm_get(gmmk_pw_pm_plus_s , pw_pm_plus )
      istat = gmm_get(gmmk_pw_pt_plus_s , pw_pt_plus )
      istat = gmm_get(gmmk_pw_me_plus_s , pw_me_plus )
      istat = gmm_get(gmmk_pw_p0_plus_s , pw_p0_plus )
      istat = gmm_get(gmmk_pw_log_pm_s  , pw_log_pm  )
      istat = gmm_get(gmmk_pw_log_pt_s  , pw_log_pt  )

      istat = gmm_get(gmmk_tt1_s, tt1)
      istat = gmm_get(gmmk_wt1_s, wt1)
      istat = gmm_get(gmmk_st1_s, st1)
      istat = gmm_get(gmmk_fis0_s, fis0)
      istat = gmm_get(gmmk_qt1_s, qt1)
      istat = gmm_get(gmmk_tt0_s, tt0)

      delps(1:l_ni,1:l_nj)=exp(lg_pstar(1:l_ni,1:l_nj,l_nk+1))*&
               (exp(qt1(1:l_ni,1:l_nj,l_nk+1)/(rgasd_8*Ver_Tstar_8%m(l_nk+1)))-&
                exp(F_q(1:l_ni,1:l_nj)/(rgasd_8*Ver_Tstar_8%m(l_nk+1))))
      do k=l_nk,1,-1
         pw_pm(1:l_ni,1:l_nj)=exp(qt1(1:l_ni,1:l_nj,k)/(rgasd_8*Ver_Tstar_8%m(k))+lg_pstar(1:l_ni,1:l_nj,k))

         qt1(1:l_ni,1:l_nj,k)=rgasd_8*Ver_Tstar_8%m(k)*&
                     (log(pw_pm(1:l_ni,1:l_nj)+delps(1:l_ni,1:l_nj)*exp(lg_pstar(1:l_ni,1:l_nj,k))/&
                     exp(lg_pstar(1:l_ni,1:l_nj,l_nk+1)))-lg_pstar(1:l_ni,1:l_nj,k))
      end do

      do k=1,l_nk
         pw_wz_plus(:,:,k) = wt1(:,:,k)
         pw_gz_plus(1:l_ni,1:l_nj,k)= grav_8*zmom(1:l_ni,1:l_nj,k)
         if(k == 1) then
            pw_me_plus(1:l_ni,1:l_nj)= fis0(1:l_ni,1:l_nj)
            pw_log_pm(1:l_ni,1:l_nj,k)=(qt1(1:l_ni,1:l_nj,k)/(rgasd_8*Ver_Tstar_8%m(k))+lg_pstar(1:l_ni,1:l_nj,k))
         end if
         pw_pm_plus(1:l_ni,1:l_nj,k)=exp(pw_log_pm(1:l_ni,1:l_nj,k))
         pw_log_pm(1:l_ni,1:l_nj,k+1)=(qt1(1:l_ni,1:l_nj,k+1)/(rgasd_8*Ver_Tstar_8%m(k+1))+lg_pstar(1:l_ni,1:l_nj,k+1))
         if(k == l_nk) then
            pw_p0_plus(1:l_ni,1:l_nj)=exp(pw_log_pm(1:l_ni,1:l_nj,l_nk+1))
         end if
      end do

      do k=1,l_nk
         pw_log_pt(1:l_ni,1:l_nj,k)=0.5*(pw_log_pm(1:l_ni,1:l_nj,k+1)+pw_log_pm(1:l_ni,1:l_nj,k))
      end do

      pw_log_pt(1:l_ni,1:l_nj,l_nk+1)=pw_log_pm(1:l_ni,1:l_nj,l_nk+1)

      do k=1,l_nk
         pw_pt_plus(1:l_ni,1:l_nj,k)=exp(pw_log_pt(1:l_ni,1:l_nj,k))
      end do

      call timing_stop (5)
!     ________________________________________________________________
!
      return
      end
