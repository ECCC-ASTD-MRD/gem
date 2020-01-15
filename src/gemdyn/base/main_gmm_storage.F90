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

!**s/r main_gmm_storage - Allocate model gmm storage

      subroutine main_gmm_storage()
      use dynkernel_options
      use glb_ld
      use gmm_geof
      use gmm_itf_mod
      use HORgrid_options
      use init_options
      use lun
      use var_gmm
      use mem_tstp
      use ldnh
      implicit none
#include <arch_specific.hf>

      integer :: istat
!
!-------------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,2000)

!     Initialize the time-dependent variables modules
!     -------------------------------------------------
      call heap_paint()

      call set_vt()

      if (Grd_yinyang_L) then
         call yyg_init()
         call yyg_initstencils()
         call yyg_rhs_initscalbc()
      else
         call nest_set_gmmvar()
      end if

!     Initialize digital filter variables modules
!     --------------------------------------------
      if ( Init_mode_L ) call set_vta()

      gmmk_fis0_s      = 'FIS0'
      gmmk_sls_s       = 'SLS'
      gmmk_topo_low_s  = 'TOPOLOW'
      gmmk_topo_high_s = 'TOPOHIGH'

      nullify (fis0, sls, topo_low, topo_high)
      istat = gmm_create(gmmk_fis0_s,fis0,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_sls_s,sls,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_get (gmmk_fis0_s, fis0)
      istat = gmm_get (gmmk_sls_s , sls )

      istat = gmm_create(gmmk_topo_low_s ,topo_low ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_topo_high_s,topo_high,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_get (gmmk_topo_low_s , topo_low )
      istat = gmm_get (gmmk_topo_high_s, topo_high)

      if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' .or. &
           trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' ) then

      allocate (rhsu(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhsv(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhst(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhsc(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhsw(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhsf(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhsp(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                rhsb(l_minx:l_maxx,l_miny:l_maxy))
      rhsu=0.;rhsv=0.;rhst=0.;rhsc=0.;rhsw=0.;rhsf=0.;rhsp=0.;rhsb=0.

      allocate (orhsu(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsv(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhst(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsc(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsw(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsf(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      orhsu=0.;orhsv=0.;orhst=0.;orhsc=0.;orhsw=0.;orhsf=0.

      allocate (nl_u(l_maxx-l_minx+1,l_maxy-l_miny+1,l_nk),&
                nl_v(l_maxx-l_minx+1,l_maxy-l_miny+1,l_nk),&
                nl_t(l_maxx-l_minx+1,l_maxy-l_miny+1,l_nk),&
                nl_c(l_maxx-l_minx+1,l_maxy-l_miny+1,l_nk),&
                nl_w(l_maxx-l_minx+1,l_maxy-l_miny+1,l_nk),&
                nl_f(l_maxx-l_minx+1,l_maxy-l_miny+1,l_nk),&
                nl_b(l_maxx-l_minx+1,l_maxy-l_miny+1))

      allocate (rhs_sol(ldnh_maxx-ldnh_minx+1,ldnh_maxy-ldnh_miny+1,l_nk),&
                lhs_sol(ldnh_maxx-ldnh_minx+1,ldnh_maxy-ldnh_miny+1,l_nk))

      endif

 2000 format( /,'INITIALIZATION OF MAIN GMM VARIABLES S/R MAIN_GMM_STORAGE', &
              /,'====================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
