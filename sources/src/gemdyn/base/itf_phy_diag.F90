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
!/@*
subroutine itf_phy_diag ()
   use phy_itf, only: phy_get
   use gmm_pw
   use HORgrid_options
   use ctrl
   use gmm_itf_mod
   use glb_ld
   use outp
      use wb_itf_mod
   implicit none
#include <arch_specific.hf>

   !@objective
   ! To compute diagnostic level values using the surface layer module.
   !@arguments
   !@author  Ron McTaggart-Cowan - Winter 2015
   !*@/

#include <rmnlib_basics.hf>
#include <msg.h>

   logical, save :: init_L = .false., dodiag_L = .true.

   ! Local variable declarations
   integer :: istat
   real :: zu,zt
   real, dimension(:,:), pointer :: tdiag,qdiag,udiag,vdiag,ptr2d
   real, dimension(:,:,:), pointer :: tt,uu,vv,hu
!
!----------------------------------------------------------------------
!
   if (.not. Ctrl_phyms_L) return
   if (.not. init_L) then
      istat = wb_get('phy/zu', zu)
      istat = wb_get('phy/zt', zt)
      dodiag_L = (zu >= 0. .and. zt >= 0.)
      init_L = .true.
   end if
   if (.not. dodiag_L) return

   ! Retrieve diagnostic level
   nullify(tdiag,qdiag,udiag,vdiag)
   istat = gmm_get(gmmk_diag_tt_s,tdiag)
   istat = gmm_get(gmmk_diag_uu_s,udiag)
   istat = gmm_get(gmmk_diag_vv_s,vdiag)
   istat = gmm_get(gmmk_diag_hu_s,qdiag)

   ! Retrieve physical world and humidity state information
   nullify(tt,uu,vv,hu)
   istat = gmm_get(gmmk_pw_tt_plus_s,tt)
   istat = gmm_get(gmmk_pw_uu_plus_s,uu)
   istat = gmm_get(gmmk_pw_vv_plus_s,vv)
   istat = gmm_get('TR/HU:P'        ,hu)

   ! Pre-fill diagnostic fields with copy-down for outside phy domain
   tdiag = tt(:,:,G_nk)
   udiag = uu(:,:,G_nk)
   vdiag = vv(:,:,G_nk)
   qdiag = hu(:,:,G_nk)

   ! Retrieve physics information only if physics is active
   if (.not.Ctrl_phyms_L) return

   ! Retrieve physics internal diagnostics

   ptr2d => tdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
   istat = phy_get(ptr2d,'tdiag',F_npath='V')
   ptr2d => udiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
   istat = phy_get(ptr2d,'udiag',F_npath='V')
   ptr2d => vdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
   istat = phy_get(ptr2d,'vdiag',F_npath='V')
   ptr2d => qdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn)
   istat = phy_get(ptr2d,'qdiag',F_npath='V')

   Outp_diag_S= ''
   return

!!$   ! Compute diagnostic level values from external surface layer module.  This
!!$   ! option is currently disabled while a better design is developed.
!!$
!!$   RECOMPUTE_DIAG: if (len_trim(Outp_diag_S) > 0) then
!!$   ! Retrieve physical world and humidity state information
!!$   nullify(gz,p0,me)
!!$   istat = gmm_get(gmmk_pw_gz_plus_s ,gz)
!!$   istat = gmm_get(gmmk_pw_p0_plus_s ,p0)
!!$   istat = gmm_get(gmmk_pw_me_moins_s,me)
!!$
!!$   hghtm = ( gz(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,G_nk) - &
!!$        me(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn) ) / grav_8
!!$   hghtt = hghtm * .5
!!$
!!$
!!$      call msg(MSG_ERROR,'itf_phy_diag','Recalculation of diagnostic quantities disabled')
!!$
!!$      ! Retrieve near-surface information from the physics
!!$      nullify(tsurf,qsurf,z0m3d,z0t3d,fcor,dlat)
!!$      istat = RMN_OK
!!$      istat = min(istat,phy_get(tsurf,'tsurf',F_npath='V'))
!!$      istat = min(istat,phy_get(qsurf,'qsurf',F_npath='V'))
!!$      istat = min(istat,phy_get(z0m3d,'z0',F_npath='V'))
!!$      z0m = z0m3d(:,:,size(z0m3d,dim=3))
!!$      istat = min(istat,phy_get(z0t3d,'z0t',F_npath='V'))
!!$      z0t = z0t3d(:,:,size(z0t3d,dim=3))
!!$      istat = phy_get(fcor,'fcor',F_npath='V')
!!$      istat = phy_get(dlat,'dlat',F_npath='V')
!!$      istat = min(istat,wb_get('phy/zu',zu))
!!$      istat = min(istat,wb_get('phy/zt',zt))
!!$      if (.not.RMN_IS_OK(istat)) call msg(MSG_ERROR,'itf_phy_diag ERROR retrieving physics information')
!!$
!!$      ! Compute diagnostic level quantities
!!$      istat = SL_OK
!!$      do j=Grd_lphy_j0,Grd_lphy_jn
!!$         jphy = j-Grd_lphy_j0+1
!!$         istat = min(istat,sl_prelim( &
!!$              t_air=tt(Grd_lphy_i0:Grd_lphy_in,j,G_nk), &
!!$              q_air=hu(Grd_lphy_i0:Grd_lphy_in,j,G_nk), &
!!$              u_air=uu(Grd_lphy_i0:Grd_lphy_in,j,G_nk), &
!!$              v_air=vv(Grd_lphy_i0:Grd_lphy_in,j,G_nk), &
!!$              p_sfc=p0(Grd_lphy_i0:Grd_lphy_in,j), &
!!$              hghtm_air=hghtm(:,jphy), &
!!$              spd_air=vmod,dir_air=vdir,tv_air=tv,rho_air=rho))
!!$         istat = min(istat,sl_sfclayer( &
!!$              t_air=tt(Grd_lphy_i0:Grd_lphy_in,j,G_nk), &
!!$              q_air=hu(Grd_lphy_i0:Grd_lphy_in,j,G_nk), &
!!$              spd_air=vmod, &
!!$              dir_air=vdir, &
!!$              hghtm_air=hghtm(:,jphy), &
!!$              hghtt_air=hghtt(:,jphy), &
!!$              t_sfc=tsurf(:,jphy,1), &
!!$              q_sfc=qsurf(:,jphy,1), &
!!$              z0m=z0m(:,jphy), &
!!$              z0t=z0t(:,jphy), &
!!$              lat=dlat(:,jphy,1), &
!!$              fcor=fcor(:,jphy,1), &
!!$              hghtm_diag=zu, &
!!$              hghtt_diag=zt, &
!!$              t_diag=my_tdiag(:,jphy), &
!!$              q_diag=my_qdiag(:,jphy), &
!!$              u_diag=my_udiag(:,jphy), &
!!$              v_diag=my_vdiag(:,jphy)))
!!$      end do
!!$      if (.not.(istat==SL_OK)) call msg(MSG_ERROR,'itf_phy_diag ERROR computing surface layer diagnostics')
!!$
!!$      ! Free space required for surface layer calculation
!!$      deallocate(ptr2d,tsurf,qsurf,z0m3d,z0t3d,fcor,dlat)
!!$
!!$   ! Update diagnostic values at user request
!!$   if (index(Outp_diag_S,'TT') > 0) tdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn) = my_tdiag
!!$   if (index(Outp_diag_S,'HU') > 0) qdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn) = my_qdiag
!!$   if (index(Outp_diag_S,'UU') > 0) udiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn) = my_udiag
!!$   if (index(Outp_diag_S,'VV') > 0) vdiag(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn) = my_vdiag
!!$   end if RECOMPUTE_DIAG
!
!----------------------------------------------------------------------
!
   return
end subroutine itf_phy_diag
