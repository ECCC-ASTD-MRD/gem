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

      subroutine adz_set()
      use adz_mem
      use adz_options
      use ctrl
      use gem_options
      use glb_pil
      use HORgrid_options
      use geomh
      use lam_options
      use gmm_itf_mod
      use gmm_pw
      use ptopo
      use rstr
      use tr3d
      use ver
      use, intrinsic :: iso_fortran_env
      implicit none

#include "gmm_gem_flags.hf"
#define SET_GMMUSR_FLAG(MYMETA,MYFLAG) gmm_metadata(MYMETA%l,gmm_attributes(MYMETA%a%key,ior(MYMETA%a%uuid1,MYFLAG),MYMETA%a%uuid2,MYMETA%a%initmode,MYMETA%a%flags))

      include "tricublin_f90.inc"

      integer  i, j, k, n, k0, BCS_BASE, pnz, ext, halom
      real(kind=REAL64), parameter :: EPS_8= 1.D-5
      real(kind=REAL64), dimension(:), allocatable :: lat
      real :: verysmall
      real(kind=REAL128) :: smallest_dz,posz

      type(gmm_metadata) :: meta, mymeta
      integer(kind=INT64) :: flag_m_f
      integer :: dim, kp1, flag_r_n, istat,istatu,istatv
!
!---------------------------------------------------------------------
!
      Adz_maxcfl= max(1,Grd_maxcfl)
      Adz_halox = Adz_maxcfl + 1
      Adz_haloy = Adz_halox

      Adz_lminx = 1    - Adz_halox
      Adz_lmaxx = l_ni + Adz_halox
      Adz_lminy = 1    - Adz_haloy
      Adz_lmaxy = l_nj + Adz_haloy
! Above lines should be replaced by the following when SLOD
! is completed. Adz_maxcfl should also no longer be needed.
!!$      Adz_halox = G_halox
!!$      Adz_haloy = G_haloy
!!$      Adz_lminx = l_minx
!!$      Adz_lmaxx = l_maxx
!!$      Adz_lminy = l_miny
!!$      Adz_lmaxy = l_maxy

      Adz_nit = Adz_lmaxx - Adz_lminx + 1
      Adz_njt = Adz_lmaxy - Adz_lminy + 1
      Adz_nij = Adz_nit * Adz_njt
      Adz_2dnh = l_ni * l_nj
      Adz_3dnh = Adz_2dnh * l_nk

      Adz_ioff = l_i0 - 2
      Adz_joff = l_j0 - 2

      ext=1
      if (Grd_yinyang_L) ext=2

      Adz_i0=   1  + pil_w - (ext+1)*west
      Adz_in= l_ni - pil_e +  ext   *east
      Adz_j0=   1  + pil_s - (ext+1)*south
      Adz_jn= l_nj - pil_n +  ext   *north

      Adz_i0u = Adz_i0 ;  Adz_inu = Adz_in
      Adz_j0v = Adz_j0 ;  Adz_jnv = Adz_jn
      if (.not. Grd_yinyang_L) then
         Adz_i0u= 1 + pil_w     ; Adz_j0v= 1 + pil_s
         Adz_inu= l_niu - pil_e ; Adz_jnv= l_njv - pil_n
      end if

      Adz_k0 = Lam_gbpil_t+1
      Adz_k0t= Adz_k0
      if(Lam_gbpil_t > 0) Adz_k0t= Adz_k0 - 1
      Adz_k0m= max(Adz_k0t-2,1)

      Adz_num_u= (Adz_inu-Adz_i0u+1)*(Adz_jn -Adz_j0 +1)*(l_nk-Adz_k0 +1)
      Adz_num_v= (Adz_in -Adz_i0 +1)*(Adz_jnv-Adz_j0v+1)*(l_nk-Adz_k0 +1)
      Adz_num_q= (Adz_in -Adz_i0 +1)*(Adz_jn -Adz_j0 +1)*(l_nk-Adz_k0 +1)
      Adz_num_t= (Adz_in -Adz_i0 +1)*(Adz_jn -Adz_j0 +1)*(l_nk-Adz_k0t+1)

      Adz_kkmax= l_nk-2

      Adz_cpntr_q= vsearch_setup_plus(Ver_z_8%m(1:G_nk), G_nk,&
                                 Adz_nit,Adz_njt,1-l_i0,1-l_j0)
      Adz_cpntr_t= vsearch_setup_plus(Ver_z_8%t(1:G_nk), G_nk,&
                                 Adz_nit,Adz_njt,1-l_i0,1-l_j0)

      BCS_BASE= 4
      if (Grd_yinyang_L) BCS_BASE = 3

      allocate (Adz_Xlim(2,0:Ptopo_npex), Adz_Ylim(2,0:Ptopo_npey))
      do i=0, Ptopo_npex-1
         Adz_Xlim(1,i) = Ptopo_gindx(1,Ptopo_colrow(Ptopo_couleur,i,0)+1)
         Adz_Xlim(2,i) = Ptopo_gindx(2,Ptopo_colrow(Ptopo_couleur,i,0)+1)
      end do
      do i=0, Ptopo_npey-1
         Adz_Ylim(1,i) = Ptopo_gindx(3,Ptopo_colrow(Ptopo_couleur,0,i)+1)
         Adz_Ylim(2,i) = Ptopo_gindx(4,Ptopo_colrow(Ptopo_couleur,0,i)+1)
      end do
      Adz_Xlim(1,0) = 1 +BCS_BASE
      Adz_Xlim(2,Ptopo_npex-1) = Ptopo_gindx&
                            (2,Ptopo_colrow(Ptopo_couleur,Ptopo_npex-1,0)+1) -BCS_BASE
      Adz_Ylim(1,0) = 1 +BCS_BASE
      Adz_Ylim(2,Ptopo_npey-1) = Ptopo_gindx&
                            (4,Ptopo_colrow(Ptopo_couleur,0,Ptopo_npey-1)+1) -BCS_BASE
      Adz_Xlim(1,Ptopo_npex) = Adz_Xlim(2,Ptopo_npex-1)
      Adz_Ylim(1,Ptopo_npey) = Adz_Ylim(2,Ptopo_npey-1)

! FOR SLOD
!!$      Adz_iminposx = l_i0+   1-Adz_halox   + EPS_8
!!$      Adz_imaxposx = l_i0+l_ni+Adz_halox-2 - EPS_8
!!$      Adz_iminposy = l_j0+   1-Adz_haloy   + EPS_8
!!$      Adz_imaxposy = l_j0+l_nj+Adz_haloy-2 - EPS_8

      Adz_iminposx = l_i0+adz_lminx   + EPS_8
      Adz_imaxposx = l_i0+adz_lmaxx-2 - EPS_8
      Adz_iminposy = l_j0+adz_lminy   + EPS_8
      Adz_imaxposy = l_j0+adz_lmaxy-2 - EPS_8
      if (l_west ) Adz_iminposx = l_i0+     BCS_BASE   + EPS_8
      if (l_east ) Adz_imaxposx = l_i0+l_ni-BCS_BASE-1 - EPS_8
      if (l_south) Adz_iminposy = l_j0+     BCS_BASE   + EPS_8
      if (l_north) Adz_imaxposy = l_j0+l_nj-BCS_BASE-1 - EPS_8

      Adz_yyminposx = 2    + Glb_pil_w     + EPS_8
      Adz_yymaxposx = G_ni - Glb_pil_e - 1 - EPS_8
      Adz_yyminposy = 2    + Glb_pil_s     + EPS_8
      Adz_yymaxposy = G_nj - Glb_pil_n - 1 - EPS_8

      ext = Grd_maxcfl + 1

      Adz_i0b =    1 + BCS_BASE*west
      Adz_inb = l_ni - BCS_BASE*east
      Adz_j0b =    1 + BCS_BASE*south
      Adz_jnb = l_nj - BCS_BASE*north

      Adz_num_b= (Adz_inb -Adz_i0b +1)*(Adz_jnb -Adz_j0b +1)*(l_nk-Adz_k0+1)

      allocate (  Adz_delz_m(0:l_nk),  Adz_delz_t(0:l_nk), &
                 Adz_odelz_m(0:l_nk), Adz_odelz_t(0:l_nk) )

! Vert coord: sig= 1.d0 for P coord, sig= -1.d0 for H coord
      sig=       (Ver_z_8%m(l_nk)-Ver_z_8%m(1)) / &
            ( abs(Ver_z_8%m(l_nk)-Ver_z_8%m(1)) )

      Adz_delz_m= 9.e33 ; Adz_delz_t= 9.e33
      verysmall= tiny(verysmall)*1.e5
      do k = 0,l_nk
         kp1 = k+1
         if (abs(Ver_z_8%m(kp1)-Ver_z_8%m(k))>verysmall) &
         Adz_delz_m(k) = Ver_z_8%m(kp1)-Ver_z_8%m(k)
         if (abs(Ver_z_8%t(kp1)-Ver_z_8%t(k))>verysmall) &
         Adz_delz_t(k) = Ver_z_8%t(kp1)-Ver_z_8%t(k)
         Adz_odelz_m(k) = 1.0d0 / Adz_delz_m(k)
         Adz_odelz_t(k) = 1.0d0 / Adz_delz_t(k)
      end do

      smallest_dz= minval(sig*Adz_delz_m(1:l_nk))
      adz_ovdzm_8 = 1.0d0/smallest_dz
      pnz = nint(1.0+ sig*(Ver_z_8%m(l_nk+1)-Ver_z_8%m(0))*adz_ovdzm_8)
      adz_ovdzm_8 = sig/smallest_dz
      allocate ( Adz_search_m(pnz) )
      k0 = 0
      do k= 1, pnz
         posz = sig*Ver_z_8%m(0) + dble(k-1) * smallest_dz
         if (posz > sig*Ver_z_8%m(k0+1)) k0 = min(l_nk+1, k0+1)
         Adz_search_m(k) = k0
      end do

      smallest_dz= minval(sig*Adz_delz_t(1:l_nk))
      adz_ovdzt_8 = 1.0d0/smallest_dz
      pnz = nint(1.0+ sig*(Ver_z_8%t(l_nk+1)-Ver_z_8%t(0))*adz_ovdzt_8)
      adz_ovdzt_8 = sig/smallest_dz
      allocate ( Adz_search_t(pnz) )
      k0 = 0
      do k= 1, pnz
         posz = sig*Ver_z_8%t(0) + dble(k-1) * smallest_dz
         if (posz > sig*Ver_z_8%t(k0+1)) k0 = min(l_nk+1, k0+1)
         Adz_search_t(k) = k0
      end do

      allocate ( Adz_cy_8(Adz_lminy:Adz_lmaxy) )
      halom= max(Adz_haloy,G_haloy)
      allocate ( lat(1-halom:G_nj+halom+1))
      lat(1-G_haloy:G_nj+G_haloy+1)= G_yg_8(1-G_haloy:G_nj+G_haloy+1)
      do j= Adz_lminy, -G_haloy
         lat(j)= G_yg_8(1-G_haloy) + (j-1+G_haloy) * geomh_hy_8
      end do
      do j = G_nj+G_haloy+2, Adz_lmaxy
         lat(j)= G_yg_8(G_nj+G_haloy) + (j-G_nj-G_haloy) * geomh_hy_8
      end do

      do j = Adz_lminy, Adz_lmaxy
         Adz_cy_8(j) = 1.d0 / cos(lat(l_j0+j-1))
      end do
      deallocate(lat)

      allocate (&
         Adz_uu_ext( Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_vv_ext( Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_ww_ext( Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_uvw_d(3,Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk), &
         Adz_uu_arr(l_ni,l_nj,l_nk), Adz_vv_arr(l_ni,l_nj,l_nk), &
         Adz_ww_arr(l_ni,l_nj,l_nk), Adz_uvw_dep(3,l_ni,l_nj,l_nk) )

      Adz_uu_ext=0. ; Adz_vv_ext=0. ; Adz_ww_ext=0.

      allocate(Adz_pm (3,Adz_i0 :Adz_in , Adz_j0 :Adz_jn ,Adz_k0 :l_nk),&
               Adz_pmu(3,Adz_i0u:Adz_inu, Adz_j0 :Adz_jn ,Adz_k0 :l_nk),&
               Adz_pmv(3,Adz_i0 :Adz_in , Adz_j0v:Adz_jnv,Adz_k0 :l_nk),&
               Adz_pt (3,Adz_i0 :Adz_in , Adz_j0 :Adz_jn ,Adz_k0t:l_nk),&
               Adz_pb (3,Adz_i0b:Adz_inb, Adz_j0b:Adz_jnb,Adz_k0t:l_nk) )

      allocate ( Adz_pxyzt (3,l_ni,l_nj,l_nk),&
                 Adz_wpxyz(-1:l_ni+2,-1:l_nj+2,l_nk,3) )
      Adz_wpxyz(-1:0,:,:,:)=0. ; Adz_wpxyz(l_ni+1:l_ni+2,:,:,:)=0.
      Adz_wpxyz(:,-1:0,:,:)=0. ; Adz_wpxyz(:,l_nj+1:l_nj+2,:,:)=0.

      flag_m_f = FLAG_LVL_M
      flag_r_n = GMM_FLAG_RSTR+GMM_FLAG_IZER

      call gmm_build_meta4D (meta,  1,3   ,0,0,   3, &
                                    1,l_ni,0,0,l_ni, &
                                    1,l_nj,0,0,l_nj, &
                                    1,l_nk,0,0,l_nk, &
                                    0,GMM_NULL_FLAGS)
      mymeta= SET_GMMUSR_FLAG(meta, flag_m_f)

      istat= gmm_create('ADZ_PXYZM',Adz_pxyzm, mymeta, flag_r_n)
      istat= gmm_get('ADZ_PXYZM',Adz_pxyzm)

!     Wind information at the lowest thermodynamic level
!     from the physics surface layer scheme
      nullify(Adz_uslt,Adz_vslt)
      if ( .not. Ctrl_phyms_L ) Adz_slt_winds= .false.
      if (Adz_slt_winds) then
         istatu= gmm_get(gmmk_pw_uslt_s, Adz_uslt)
         istatv= gmm_get(gmmk_pw_vslt_s, Adz_vslt)
         if ((istatu/=0).or.(istatv/=0)) Adz_slt_winds= .false.
      end if

      allocate (Adz_stack(max(10,Tr3d_ntr)))

     !if (.not.Grd_yinyang_L) allocate(Adz_stack(max(1,Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP))%pil &
     !                                 (l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      !This works but not the previous when Tr3d_ntrTRICUB_WP>2 (Dont know why...)
      !---------------------------------------------------------------------------
      if (.not.Grd_yinyang_L) then

         do n = 1,max(1,Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)

            allocate(Adz_stack(n)%pil(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

         end do

      end if

      !Allocation Conservation Postprocessing
      !--------------------------------------
      allocate (adz_flux     (max(1,Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)))
      allocate (adz_flux_3CWP(max(1,Tr3d_ntrTRICUB_WP)))
      allocate (adz_flux_BQWP(max(1,Tr3d_ntrBICHQV_WP)))

      do n = 1,max(1,Tr3d_ntrTRICUB_WP)

         allocate(adz_flux_3CWP(n)%fo(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_flux_3CWP(n)%fi(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      do n = 1,max(1,Tr3d_ntrBICHQV_WP)

         allocate(adz_flux_BQWP(n)%fo(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_flux_BQWP(n)%fi(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      allocate (adz_post     (max(1,Tr3d_ntrTRICUB_WP,Tr3d_ntrBICHQV_WP)))
      allocate (adz_post_3CWP(max(1,Tr3d_ntrTRICUB_WP)))
      allocate (adz_post_BQWP(max(1,Tr3d_ntrBICHQV_WP)))

      do n = 1,max(1,Tr3d_ntrTRICUB_WP)

         allocate(adz_post_3CWP(n)%lin(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_3CWP(n)%min(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_3CWP(n)%max(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      do n = 1,max(1,Tr3d_ntrBICHQV_WP)

         allocate(adz_post_BQWP(n)%lin(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_BQWP(n)%min(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))
         allocate(adz_post_BQWP(n)%max(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      allocate (adz_bc(1+Tr3d_ntrTRICUB_WP+Tr3d_ntrBICHQV_WP))

      do n = 1,1+Tr3d_ntrTRICUB_WP+Tr3d_ntrBICHQV_WP

         allocate(adz_bc(n)%w(l_minx:l_maxx,l_miny:l_maxy,1:l_nk))

      end do

      Adz_niter = Adz_itraj
      if ( .not. Rstri_rstn_L ) call adz_inittraj

      dim = l_ni*l_nj*l_nk
      allocate ( Adz_expq%ijk (  dim  ), Adz_expq%gpos(3*dim,8),&
                 Adz_expq%resu(  dim  ), Adz_expq%req ( dim   ),&
                 Adz_expq%intp(  dim  )                        ,&
                 Adz_expu%ijk (  dim  ), Adz_expu%gpos(3*dim,8),&
                 Adz_expu%resu(  dim  ), Adz_expu%req ( dim   ),&
                 Adz_expu%intp(  dim  )                        ,&
                 Adz_expv%ijk (  dim  ), Adz_expv%gpos(3*dim,8),&
                 Adz_expv%resu(  dim  ), Adz_expv%req ( dim   ),&
                 Adz_expv%intp(  dim  )                        ,&
                 Adz_expt%ijk (  dim  ), Adz_expt%gpos(3*dim,8),&
                 Adz_expt%resu(  dim  ), Adz_expt%req ( dim   ),&
                 Adz_expt%intp(  dim  )                        ,&
                 Adz_expq%COMM_handle(Ptopo_numproc)           ,&
                 Adz_expu%COMM_handle(Ptopo_numproc)           ,&
                 Adz_expv%COMM_handle(Ptopo_numproc)           ,&
                 Adz_expt%COMM_handle(Ptopo_numproc)            )

      Adz_exppe = -1
! later in the one-sided MPI in GY
!!$      if (Ptopo_myrow>0) then
!!$         if (Ptopo_mycol<Ptopo_npex-1) Adz_exppe(1)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow-1)
!!$         if (Ptopo_mycol>0           ) Adz_exppe(3)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow-1)
!!$                                       Adz_exppe(2)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol  ,Ptopo_myrow-1)
!!$      endif
!!$      if (Ptopo_myrow<Ptopo_npey-1) then
!!$         if (Ptopo_mycol<Ptopo_npex-1) Adz_exppe(7)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow+1)
!!$         if (Ptopo_mycol>0           ) Adz_exppe(5)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow+1)
!!$                                       Adz_exppe(6)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol  ,Ptopo_myrow+1)
!!$      endif
!!$      if (Ptopo_mycol>0              ) Adz_exppe(4)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow  )
!!$      if (Ptopo_mycol<Ptopo_npex-1   ) Adz_exppe(8)=Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow  )

      if (Ptopo_myrow>0) then
         if (Ptopo_mycol<Ptopo_npex-1) Adz_exppe(1)=Ptopo_colrow(0,Ptopo_mycol+1,Ptopo_myrow-1)
         if (Ptopo_mycol>0           ) Adz_exppe(3)=Ptopo_colrow(0,Ptopo_mycol-1,Ptopo_myrow-1)
                                       Adz_exppe(2)=Ptopo_colrow(0,Ptopo_mycol  ,Ptopo_myrow-1)
      endif
      if (Ptopo_myrow<Ptopo_npey-1) then
         if (Ptopo_mycol<Ptopo_npex-1) Adz_exppe(7)=Ptopo_colrow(0,Ptopo_mycol+1,Ptopo_myrow+1)
         if (Ptopo_mycol>0           ) Adz_exppe(5)=Ptopo_colrow(0,Ptopo_mycol-1,Ptopo_myrow+1)
                                       Adz_exppe(6)=Ptopo_colrow(0,Ptopo_mycol  ,Ptopo_myrow+1)
      endif
      if (Ptopo_mycol>0              ) Adz_exppe(4)=Ptopo_colrow(0,Ptopo_mycol-1,Ptopo_myrow  )
      if (Ptopo_mycol<Ptopo_npex-1   ) Adz_exppe(8)=Ptopo_colrow(0,Ptopo_mycol+1,Ptopo_myrow  )
!
!---------------------------------------------------------------------
!
      return

 1000 format(/,'=================================================',/,A31,I2,/,&
               '=================================================' )

      end subroutine adz_set
