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

!**s/r out_dyn_casc - model output for cascade

      subroutine out_dyn_casc()
      use dyn_fisl_options
      use dynkernel_options
      use glb_ld
      use HORgrid_options
      use gmm_geof
      use gmm_itf_mod
      use gmm_pw
      use gmm_vt1
      use grdc_options
      use levels
      use metric
      use out3
      use out_mod
      use out_options
      use out_vref
      use outp
      use tdpack
      use phy_itf, only: phy_get
      use vertical_interpolation
      use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_free,VGD_OK,VGD_ERROR
      use vgrid_wb, only: vgrid_wb_get
      implicit none
#include <arch_specific.hf>


      character(len=512) varname
      integer :: k, istat, indo(G_nk+2)
      integer, dimension(:), pointer  :: ip1m
      real, dimension(:    ), pointer :: hybm,hybt,hybt_w
      real, dimension(:,:  ), pointer :: tdiag,udiag,vdiag,qdiag
      real, dimension(:,:,:), pointer :: ptr3d,tr2,wlnph_ta,wlnph_m

      real, dimension(l_minx:l_maxx,l_miny:l_maxy,1:G_nk+1),target :: tr1
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk+1) :: gzm, gzt
      type(vgrid_descriptor) :: vcoord
      real, dimension(1) :: hybm_gnk2, hybt_gnk2, hyb0
      integer, dimension(1) :: ind0
      real :: conv
!
!------------------------------------------------------------------
!
      nullify (pw_tt_plus,pw_uu_plus,pw_vv_plus,tdiag,udiag,vdiag)
      istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)
      istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get (gmmk_pw_p0_plus_s, pw_p0_plus)
      istat = gmm_get (gmmk_diag_tt_s   , tdiag     )
      istat = gmm_get (gmmk_diag_uu_s   , udiag     )
      istat = gmm_get (gmmk_diag_vv_s   , vdiag     )

      istat = gmm_get(gmmk_zdt1_s,zdt1)
      istat = gmm_get(gmmk_fis0_s,fis0)
      istat = gmm_get(gmmk_wt1_s ,wt1 )

      do k=1,G_nk+2
         indo(k) = k
      end do

      gzm = 0.
      gzt = 0.

      nullify (ip1m,hybm,hybt,hybt_w)
      istat = vgrid_wb_get('ref-m',vcoord,ip1m)
      deallocate(ip1m); nullify(ip1m)
      istat = vgd_get (vcoord,'VCDM - vertical coordinate (m)',hybm)
      istat = vgd_get (vcoord,'VCDT - vertical coordinate (t)',hybt)
      istat = vgd_free(vcoord)

      hyb0(1)=0.0
      hybm_gnk2(1)=hybm(G_nk+2)
      hybt_gnk2(1)=hybt(G_nk+2)
      ind0(1)=1

      call itf_phy_diag()

      Out_reduc_l = .true.

      call out_href ( 'Mass_point', Grdc_gid, Grdc_gif, 1,&
                                     Grdc_gjd, Grdc_gjf, 1 )

      call out_vref_itf ( etiket=Out3_etik_S )

      conv = -tcdk_8

      call out_fstecr ( pw_tt_plus,l_minx,l_maxx,l_miny,l_maxy,hybt,&
                         'TT  ',1., conv,Level_kind_ip1,-1,G_nk,indo,&
                         G_nk,Grdc_nbits,.false. )
      if (Out3_sfcdiag_L) then
         call out_fstecr ( tdiag,l_minx,l_maxx,l_miny,l_maxy ,&
                            hybt_gnk2, 'TT  ', 1., conv,4,-1,1,&
                            ind0,1,Grdc_nbits, .false. )
      end if

      if (.not. Schm_autobar_L) then
         call out_fstecr ( pw_p0_plus,l_minx,l_maxx,l_miny,l_maxy,hyb0,&
                   'P0  ',.01, 0., 2,-1,1, ind0, 1, Grdc_nbits, .false. )
      end if

      call out_fstecr ( wt1 ,l_minx,l_maxx,l_miny,l_maxy, hybt,&
                         'WT1 ',1., 0.,Level_kind_ip1,-1,G_nk  ,&
                          indo,G_nk,Grdc_nbits,.false. )
      call out_fstecr ( zdt1,l_minx,l_maxx,l_miny,l_maxy, hybt,&
                         'ZDT1',1., 0.,Level_kind_ip1,-1,G_nk  ,&
                         indo,G_nk,Grdc_nbits,.false. )

      if ( Dynamics_hauteur_L ) then
         istat = gmm_get(gmmk_qt1_s ,qt1 )
         call out_fstecr ( qt1,l_minx,l_maxx,l_miny,l_maxy, hybm,&
                         'QT1',1., 0.,Level_kind_ip1,-1,G_nk+1 ,&
                         indo,G_nk+1,Grdc_nbits,.false. )
         if (.not. Schm_autobar_L) then
            conv = 0.1d0
            call out_fstecr ( zmom(l_minx,l_miny,1), l_minx,l_maxx,&
                    l_miny,l_maxy,hybm,'GZ',conv, 0.,Level_kind_ip1,&
                    -1,G_nk+1,indo,G_nk,Grdc_nbits,.false. )
            call out_fstecr ( ztht(l_minx,l_miny,1), l_minx,l_maxx,&
                    l_miny,l_maxy,hybt,'GZ',conv, 0.,Level_kind_ip1,&
                    -1,G_nk+1,indo,G_nk+1,Grdc_nbits,.false. )
         else
            conv = 0.1d0 / grav_8
            istat = gmm_get(gmmk_qt1_s,qt1)
            gzm(:,:,1:G_nk+1) = qt1(:,:,1:G_nk+1) + 1.0d0/Cstv_invFI_8
            call out_fstecr ( gzm(l_minx,l_miny,1), l_minx,l_maxx ,&
                    l_miny,l_maxy,hybm,'GZ',conv, 0.,Level_kind_ip1,&
                    -1,G_nk+1,indo,G_nk,Grdc_nbits,.false. )
         end if
      else
         istat = gmm_get(gmmk_qt1_s      , qt1     )
         istat = gmm_get(gmmk_pw_log_pt_s, wlnph_ta)
         istat = gmm_get(gmmk_pw_log_pm_s, wlnph_m )
         call diag_fi (gzm, st1, tt1, qt1, l_minx,l_maxx,l_miny,l_maxy,&
                       G_nk, 1, l_ni, 1, l_nj)
         gzt(:,:,l_nk+1)= gzm(:,:,l_nk+1)
         call vertint2 ( gzt, wlnph_ta,G_nk, gzm, wlnph_m,G_nk+1  ,&
                         l_minx,l_maxx,l_miny,l_maxy,1,l_ni,1,l_nj )
         conv = 0.1d0 / grav_8
         call out_fstecr ( gzm(l_minx,l_miny,1), l_minx,l_maxx ,&
                 l_miny,l_maxy,hybm,'GZ',conv, 0.,Level_kind_ip1,&
                 -1,G_nk+1,indo,G_nk,Grdc_nbits,.false. )
         if (.not.Schm_autobar_L) then
         call out_fstecr ( gzt(l_minx,l_miny,1), l_minx,l_maxx ,&
                 l_miny,l_maxy,hybt,'GZ',conv, 0.,Level_kind_ip1,&
                 -1,G_nk+1,indo,G_nk,Grdc_nbits,.false. )
         call out_fstecr ( gzt(l_minx,l_miny,1), l_minx,l_maxx ,&
                 l_miny,l_maxy,hybt,'GZ',conv, 0.,Level_kind_ip1,&
                 -1,G_nk+1,indo(G_nk+1),1,Grdc_nbits,.false. )
         end if
      endif

      do k=1,Grdc_ntr
         nullify (tr2)
         varname = 'TR/'//trim(Grdc_trnm_S(k))//':P'
         istat= gmm_get (varname,tr2)
         call out_fstecr ( tr2 ,l_minx,l_maxx,l_miny,l_maxy,hybt, &
                            Grdc_trnm_S(k),1.,0.,Level_kind_ip1,-1,&
                            G_nk, indo, G_nk, Grdc_nbits,.false. )
         if ( Out3_sfcdiag_L ) then
            if (trim(varname)=='TR/HU:P') then
               istat = gmm_get(gmmk_diag_hu_s,qdiag)
               if (istat == 0) &
               call out_fstecr ( qdiag ,l_minx,l_maxx,l_miny,l_maxy, &
                                  hybt_gnk2,Grdc_trnm_S(k),1.,0.,4, &
                                  -1,1,ind0,1,Grdc_nbits,.false. )
            else
               tr1(:,:,G_nk+1) = tr2(:,:,G_nk)
               ptr3d => tr1(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:G_nk+1)
               istat = phy_get(ptr3d, varname, F_npath='VO', F_bpath='D')
               if (istat == 0) &
               call out_fstecr ( tr1(l_minx,l_miny,G_nk+1), &
                                l_minx,l_maxx,l_miny,l_maxy, &
                                hybt_gnk2,Grdc_trnm_S(k),1.,0.,4, &
                                -1,1,ind0,1,Grdc_nbits,.false. )
            endif
         end if
      end do

      conv = 1.d0 / knams_8
      call out_fstecr ( pw_uu_plus, l_minx,l_maxx,l_miny,l_maxy, hybm,&
                         'UU  ',conv, 0., Level_kind_ip1,-1,G_nk,indo ,&
                         G_nk,Grdc_nbits,.false. )
      call out_fstecr ( pw_vv_plus, l_minx,l_maxx,l_miny,l_maxy, hybm,&
                         'VV  ',conv, 0., Level_kind_ip1,-1,G_nk,indo ,&
                         G_nk,Grdc_nbits,.false. )
      if (Out3_sfcdiag_L) then
         call out_fstecr ( udiag,l_minx,l_maxx,l_miny,l_maxy,hybm_gnk2,&
                    'UU  ' , conv, 0., 4,-1,1,ind0,1,Grdc_nbits,.false. )
         call out_fstecr (vdiag, l_minx,l_maxx,l_miny,l_maxy,hybm_gnk2,&
                    'VV  ' , conv, 0., 4,-1,1,ind0,1,Grdc_nbits,.false. )
      end if

      deallocate (hybm,hybt)
!
!------------------------------------------------------------------
!
      return
      end

