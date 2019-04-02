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

!**   s/r metric_H - calculate metric coefficients for GEM in height-base coordinates
   subroutine metric_H
      use mtn_options
      use HORgrid_options
      use gmm_geof
      use geomh
      use tdpack
      use gmm_itf_mod
      use glb_ld
      use lun
      use metric
      use cstv
      use ver
      implicit none
#include <arch_specific.hf>

      integer :: istat,i,j,k
      real, parameter :: zero=0.d0, one=1.d0, half=.5d0
      logical,save :: done = .false.
!
!     ---------------------------------------------------------------

      if(.not.done) then
         allocate ( zmom(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    ztht(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                lg_pstar(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )

         allocate ( mc_Jx (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Jy (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_iJz(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                  mc_logJz(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Ix (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Iy (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    mc_Iz (l_minx:l_maxx,l_miny:l_maxy,G_nk) )

         allocate ( mc_css_H_8   (l_minx:l_maxx,l_miny:l_maxy), &
                    mc_alfas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                    mc_betas_H_8 (l_minx:l_maxx,l_miny:l_maxy))
       done=.true.
      end if

      istat = gmm_get (gmmk_fis0_s,fis0)
      istat = gmm_get (gmmk_sls_s,sls)

      ztht(:,:,0)=ver_z_8%m(0)
      zmom(:,:,0)=ver_z_8%m(0)

      do k=1,G_nk
         zmom(:,:,k)=ver_z_8%m(k)+(Ver_b_8%m(k)*fis0(:,:)+Ver_c_8%m(k)*sls(:,:))/grav_8
         ztht(:,:,k)=ver_z_8%t(k)+(Ver_b_8%t(k)*fis0(:,:)+Ver_c_8%t(k)*sls(:,:))/grav_8
      end do
      zmom(:,:,G_nk+1)=fis0(:,:)/grav_8
      ztht(:,:,G_nk+1)=fis0(:,:)/grav_8

      lg_pstar(:,:,G_nk+1)=log(1.d5)-grav_8*zmom(:,:,G_nk+1)/(rgasd_8*Cstv_Tstr_8)

      do k=G_nk,1,-1
         lg_pstar(:,:,k)=lg_pstar(:,:,k+1)+grav_8*(zmom(:,:,k+1)-zmom(:,:,k))/(rgasd_8*Cstv_Tstr_8)
      enddo

      do k=1,G_nk
         do j=l_miny+1,l_maxy-1
         do i=l_minx+1,l_maxx-1
            mc_Jx (i,j,k)=(zmom(i+1,j,k)-zmom(i,j,k))*geomh_invDX_8(j)
            mc_Jy (i,j,k)=(zmom(i,j+1,k)-zmom(i,j,k))*geomh_invDY_8
            mc_iJz(i,j,k)=one/(zmom(i,j,k+1)-zmom(i,j,k))
         enddo
         enddo
      enddo
      ztht(:,:,G_nk)=fis0(:,:)/grav_8
      do k=1,G_nk
         do j=l_miny+1,l_maxy-1
         do i=l_minx+1,l_maxx-1
      !     mc_Ix(i,j,k)=log( (ztht(i+1,j,k)-ztht(i+1,j,k-1))/(ztht(i-1,j,k)-ztht(i-1,j,k-1)) )*0.5d0*geomh_invDX_8(j)
      !     mc_Iy(i,j,k)=log( (ztht(i,j+1,k)-ztht(i,j+1,k-1))/(ztht(i,j-1,k)-ztht(i,j-1,k-1)) )*0.5d0*geomh_invDY_8
      !     mc_Iz(i,j,k)=log( (zmom(i,j,k+1)-zmom(i,j,k))/(Ver_z_8%m(k+1)-Ver_z_8%m(k)) &
      !                      /(zmom(i,j,k)-zmom(i,j,k-1))*(Ver_z_8%m(k)-Ver_z_8%m(k-1)) )*Ver_idz_8%m(k)
           mc_Ix(i,j,k)=0.0
           mc_Iy(i,j,k)=0.0
           mc_Iz(i,j,k)=0.0
           mc_logJz(i,j,k)=log( (ztht(i,j,k)-ztht(i,j,k-1))/(Ver_z_8%x(k)-Ver_z_8%x(k-1)) )
       !    mc_logJz(i,j,k)=0.0
         enddo
         enddo
      enddo
      ztht(:,:,G_nk)=ver_z_8%t(G_nk)+(Ver_b_8%t(G_nk)*fis0(:,:)+Ver_c_8%t(G_nk)*sls(:,:))/grav_8

      if (trim(sol_type_S) == 'ITERATIVE_3D') then

         do j=l_miny+1,l_maxy-1
            do i=l_minx+1,l_maxx-1

               mc_css_H_8(i,j)   = one/(gama_8*(mc_iJz(i,j,G_nk)-half*mu_8))

               mc_alfas_H_8(i,j) = mc_css_H_8(i,j)*(gama_8*(mc_iJz(i,j,G_nk  )+half*mu_8) &
                          +Ver_wmstar_8(G_nk)*gama_8*(mc_iJz(i,j,G_nk-1)-half*mu_8) )

               mc_betas_H_8(i,j) = mc_css_H_8(i,j)* &
                           Ver_wmstar_8(G_nk)*gama_8*(mc_iJz(i,j,G_nk-1)+half*mu_8)
            enddo
         enddo
      else

         do j=l_miny+1,l_maxy-1
            do i=l_minx+1,l_maxx-1
                mc_css_H_8(i,j)   = one/(gama_8*(Ver_idz_8%t(G_nk)-half*mu_8))

                mc_alfas_H_8(i,j) = mc_css_H_8(i,j)*(gama_8*(Ver_idz_8%t(G_nk  )+half*mu_8)   &
                          +Ver_wmstar_8(G_nk)*gama_8*(Ver_idz_8%t(G_nk-1)-half*mu_8) )

                mc_betas_H_8(i,j) = mc_css_H_8(i,j)* &
                           Ver_wmstar_8(G_nk)*gama_8*(Ver_idz_8%t(G_nk-1)+half*mu_8)
            enddo
         enddo
      endif


!     ---------------------------------------------------------------
!
      return
   end subroutine metric_H
