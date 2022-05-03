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
      use adz_mem
      use dynkernel_options
      use dyn_fisl_options  
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_itf_mod
      use HORgrid_options
      use init_options
      use lun
      use metric
      use var_gmm
      use gmm_phy
      use mem_tstp
      use ldnh
      use tr3d
      implicit none
#include <arch_specific.hf>

#include "gmm_gem_flags.hf"
#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      type(gmm_metadata) :: meta, mymeta
      integer(kind=INT64) :: flag_m_f
      integer :: istat,dim
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
         call nest_set_mem
      end if

!     Initialize digital filter variables modules
!     --------------------------------------------
      if ( Init_mode_L ) call set_vta()

      nullify (fis0, orols, sls, topo_low, topo_high, me_full, me_large)
      istat = gmm_create(gmmk_fis0_s ,fis0 ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_orols_s,orols,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_sls_s  ,sls  ,meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)

      call gmm_build_meta3D(meta, &
                            l_minx,l_maxx,G_halox,G_halox,l_ni, &
                            l_miny,l_maxy,G_haloy,G_haloy,l_nj, &
                            1,2,0,0,2, 0,GMM_NULL_FLAGS)
      flag_m_f = FLAG_LVL_M
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)

      istat = gmm_create(gmmk_topo_low_s ,topo_low ,mymeta,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_topo_high_s,topo_high,mymeta,GMM_FLAG_RSTR+GMM_FLAG_IZER)

      istat = gmm_create(gmmk_me_full_s, me_full  , meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      istat = gmm_create(gmmk_me_large_s, me_large, meta2d,GMM_FLAG_RSTR+GMM_FLAG_IZER)
      
      allocate (rhs_zero(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      allocate (rhsb(l_minx:l_maxx,l_miny:l_maxy))
      rhsb = 0.
      if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' .or. &
           trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' ) then

      if (Schm_phycpl_S == 'RHS' .or. Schm_phycpl_S == 'AVG') then
         rhs_phytv= 1.d0
         rhs_uu_tend => phy_uu_tend
         rhs_vv_tend => phy_vv_tend
      else
         rhs_phytv= 0.d0
         rhs_uu_tend => rhs_zero
         rhs_vv_tend => rhs_zero
      endif
      dim= (l_maxx-l_minx+1)*(l_maxy-l_miny+1)*G_nk
      allocate (rhs(6*dim))
      rhsu (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(      1:)
      rhsv (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(  dim+1:)
      rhst (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(2*dim+1:)
      rhsc (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(3*dim+1:)
      rhsw (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(4*dim+1:)
      rhsf (l_minx:l_maxx,l_miny:l_maxy,1:l_nk)=>rhs(5*dim+1:)

      rhsu=0.;rhsv=0.;rhst=0.;rhsc=0.;rhsw=0.;rhsf=0.;rhs_zero=0.

      allocate (orhsu(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsv(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhst(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsc(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsw(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                orhsf(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      orhsu=0.;orhsv=0.;orhst=0.;orhsc=0.;orhsw=0.;orhsf=0.

      allocate (nl_u(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_v(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_t(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_c(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_w(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_f(l_minx:l_maxx,l_miny:l_maxy,l_nk),&
                nl_b(l_minx:l_maxx,l_miny:l_maxy))
      nl_b = 0.
      allocate (rhs_sol(ldnh_maxx-ldnh_minx+1,ldnh_maxy-ldnh_miny+1,l_nk),&
                lhs_sol(ldnh_maxx-ldnh_minx+1,ldnh_maxy-ldnh_miny+1,l_nk))
      rhs_sol= 0.

      endif
             
      if (Dynamics_hauteur_L) then
         allocate ( GVM%zmom_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    GVM%ztht_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    GVM%zmom_u(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    GVM%ztht_u(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    GVM%zmom_v(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    GVM%ztht_v(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1), &
                    GVM%lg_pstar_8(l_minx:l_maxx,l_miny:l_maxy,0:G_nk+1) )

         allocate ( GVM%mc_Jx_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    GVM%mc_Jy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    GVM%mc_iJz_8(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                  GVM%mc_logJz_8(l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    GVM%mc_Ix_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    GVM%mc_Iy_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk), &
                    GVM%mc_Iz_8 (l_minx:l_maxx,l_miny:l_maxy,G_nk) )

         allocate ( GVM%mc_css_H_8   (l_minx:l_maxx,l_miny:l_maxy), &
                    GVM%mc_alfas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                    GVM%mc_betas_H_8 (l_minx:l_maxx,l_miny:l_maxy), &
                    GVM%mc_cssp_H_8  (l_minx:l_maxx,l_miny:l_maxy) )

         allocate ( GVM%mc_cst_8   (l_minx:l_maxx,l_miny:l_maxy), &
                    GVM%mc_alfat_8 (l_minx:l_maxx,l_miny:l_maxy), &
                    GVM%mc_cstp_8  (l_minx:l_maxx,l_miny:l_maxy) )
         GVM%mc_cst_8= 0. ; GVM%mc_alfat_8= 0. ; GVM%mc_cstp_8= 0.
     endif

 2000 format( /,'INITIALIZATION OF MAIN GMM VARIABLES S/R MAIN_GMM_STORAGE', &
              /,'====================================================')
!
!     ---------------------------------------------------------------
!
      return
      end
