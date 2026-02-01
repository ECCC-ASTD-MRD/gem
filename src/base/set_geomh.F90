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
!**s/r set_geomh - initialize model geometry

      subroutine set_geomh()
      use dcst
      use gem_options
      use ctrl
      use geomh
      use glb_ld
      use glb_pil
      use HORgrid_options
      use hgc
      use hgrid_wb, only: hgrid_wb_put
      use lun
      use ptopo
      use tdpack
      use wb_itf_mod
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, external :: ezgdef_fmem,gdll
      integer offi,offj,indx,err,dgid,hgc_array(4), nicore, njcore
      integer gphy_i0, gphy_in, gphy_j0, gphy_jn, gphy_ni, gphy_nj, gphy_nicore, gphy_njcore
      integer i,j,istat,ni,nj,offset, glbphy_gid, glbphycore_gid, dimx,dimy
      real xfi(0:l_ni+1),yfi(0:l_nj+1)
      real gxfi(G_ni),gyfi(G_nj)
      real(kind=REAL64) posx_8(1-G_halox:G_ni+G_halox+1), posy_8(1-G_haloy:G_nj+G_haloy+1)
      real(kind=REAL64) rad2deg_8,deg2rad_8,x0,xl,y0,yl
      real(kind=REAL64) scale_factor, scale_factor_v
      real(kind=REAL64) Del_xg , Del_yg
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,1000)

      rad2deg_8 = 180.0d0/pi_8 !Should not need this..
      deg2rad_8 = pi_8/180.0d0

      geomh_minx= west
      geomh_miny= south
      geomh_maxx= l_ni + 1 - east
      geomh_maxy= l_nj + 1 - north
      
      allocate (G_xg_8(1-G_halox:G_ni+G_halox+1) , G_yg_8(1-G_haloy:G_nj+G_haloy+1) )
      allocate (geomh_latrx(l_ni,l_nj), geomh_lonrx(l_ni,l_nj))
      allocate (geomh_lonQ(1-G_halox:G_ni+G_halox), geomh_latQ(1-G_haloy:G_nj+G_haloy),&
                geomh_lonF(1-G_halox:G_ni+G_halox), geomh_latF(1-G_haloy:G_nj+G_haloy) )
      allocate (geomh_latij(geomh_minx:geomh_maxx,geomh_miny:geomh_maxy), &
                geomh_lonij(geomh_minx:geomh_maxx,geomh_miny:geomh_maxy))

      dimx= G_ni+2*G_halox ; dimy= G_nj+2*G_haloy
      allocate (geomh_4output(8+2*dimx+2*dimy))
      geomh_4output(1) = Grd_rot_xg(1,1)
      geomh_4output(2) = Grd_rot_xg(2,1)
      geomh_4output(3) = Grd_rot_xg(3,1)
      geomh_4output(4) = Grd_rot_xg(4,1)
      geomh_4output(5) = Grd_rot_xg(1,2)
      geomh_4output(6) = Grd_rot_xg(2,2)
      geomh_4output(7) = Grd_rot_xg(3,2)
      geomh_4output(8) = Grd_rot_xg(4,2)

      geomh_latgs(1-G_haloy:G_nj+G_haloy) => geomh_4output(            9:)
      geomh_longs(1-G_halox:G_ni+G_halox) => geomh_4output(dimy       +9:)
      geomh_latgv(1-G_haloy:G_nj+G_haloy) => geomh_4output(dimx+  dimy+9:)
      geomh_longu(1-G_halox:G_ni+G_halox) => geomh_4output(dimx+2*dimy+9:)
      
      ni= Grd_ni + 2*G_halox + 1
      nj= Grd_nj + 2*G_haloy + 1
      x0= Grd_x0_8 - dble(G_halox)  *dble(Grd_dx)
      y0= Grd_y0_8 - dble(G_haloy)  *dble(Grd_dy)
      xl= Grd_xl_8 + dble(G_halox+1)*dble(Grd_dx)
      yl= Grd_yl_8 + dble(G_haloy+1)*dble(Grd_dy)

      call set_gemHgrid4 ( posx_8, posy_8, ni, nj, &
                           Grd_dx, Grd_dy, x0,xl,y0,yl, Grd_yinyang_L )

      G_xg_8(1-G_halox:G_ni+G_halox+1) = posx_8(1-G_halox:G_ni+G_halox+1)*deg2rad_8
      G_yg_8(1-G_haloy:G_nj+G_haloy+1) = posy_8(1-G_haloy:G_nj+G_haloy+1)*deg2rad_8
      if (posy_8(1-G_haloy) <= -90.0.or.posy_8(G_nj+G_haloy+1) >= 90.0) then
          if (Lun_out > 0) write (Lun_out,1001) posy_8(1-G_haloy),posy_8(G_nj+G_haloy+1)
          call gem_error ( -1,'set_geomh','' )
      end if

      geomh_lonQ(1-G_halox:G_ni+G_halox) = posx_8(1-G_halox:G_ni+G_halox)
      geomh_latQ(1-G_haloy:G_nj+G_haloy) = posy_8(1-G_haloy:G_nj+G_haloy)

      do i= 1-G_halox, G_ni+G_halox
         geomh_lonF(i) = (posx_8(i+1) + posx_8(i)) * 0.5d0
      end do
      do j= 1-G_haloy, G_nj+G_haloy
         geomh_latF(j) = (posy_8(j+1) + posy_8(j)) * 0.5d0
      end do

      do i= 1-G_halox, G_ni+G_halox
         geomh_longs(i) =  posx_8(i)
         geomh_longu(i) = geomh_lonF(i)
      end do

      do i= 1-G_haloy, G_nj+G_haloy
         geomh_latgs(i) =  posy_8(i)
         geomh_latgv(i) = geomh_latF(i)
      end do

      allocate (geomh_x_8 (l_minx:l_maxx),geomh_xu_8 (l_minx:l_maxx),&
                geomh_sx_8(l_minx:l_maxx),geomh_sy_8 (l_miny:l_maxy),&
                geomh_cx_8(l_minx:l_maxx),geomh_cy_8 (l_miny:l_maxy),&

                geomh_y_8      (l_miny:l_maxy),geomh_yv_8     (l_miny:l_maxy),&
                geomh_cy2_8    (l_miny:l_maxy),geomh_cyv_8    (l_miny:l_maxy),&
                geomh_tyoa_8   (l_miny:l_maxy),geomh_tyoav_8  (l_miny:l_maxy),&
                geomh_cyv2_8   (l_miny:l_maxy),geomh_cyM_8    (l_miny:l_maxy),&
                geomh_invDYM_8 (l_miny:l_maxy),geomh_invDYMv_8(l_miny:l_maxy),&
                geomh_invcy_8  (l_miny:l_maxy),geomh_invcyv_8 (l_miny:l_maxy),&
                geomh_invcy2_8 (l_miny:l_maxy),geomh_invcyv2_8(l_miny:l_maxy),&
                geomh_invDX_8  (l_miny:l_maxy),geomh_invDXM_8 (l_miny:l_maxy),&
                geomh_invDXMu_8(l_miny:l_maxy),geomh_invDXv_8 (l_miny:l_maxy),&
                geomh_area_8(l_ni,l_nj),geomh_mask_8(l_ni,l_nj),geomh_area_mask_8(l_ni,l_nj))

      offi = Ptopo_gindx(1,Ptopo_myproc+1)-1
      offj = Ptopo_gindx(3,Ptopo_myproc+1)-1

      Del_xg= G_xg_8(2) - G_xg_8(1)
      Del_yg= G_yg_8(2) - G_yg_8(1)

      geomh_hx_8 = Del_xg
      geomh_hy_8 = Del_yg
      geomh_inv_hx_8 = 1.d0 / geomh_hx_8
      geomh_inv_hy_8 = 1.d0 / geomh_hy_8

      geomh_invDY_8  = 1.0d0 / ( Dcst_rayt_8 * Del_yg )

      do i=1-G_halox, l_ni+G_halox
         indx = offi + i
         geomh_x_8 (i) =  G_xg_8(indx)
         geomh_xu_8(i) = (G_xg_8(indx+1) + G_xg_8(indx)) * 0.5d0
         geomh_sx_8(i) = sin( geomh_x_8(i) )
         geomh_cx_8(i) = cos( geomh_x_8(i) )
      end do

      do j=1-G_haloy, l_nj+G_haloy
         indx = offj + j
         geomh_y_8  (j) =  G_yg_8(indx)
         geomh_yv_8 (j) = (G_yg_8(indx+1) + G_yg_8(indx)) * 0.5d0
         geomh_sy_8  (j)= sin( geomh_y_8 (j) )
         geomh_cy_8  (j)= cos( geomh_y_8 (j) )
         geomh_cy2_8 (j)= cos( geomh_y_8 (j) )**2
         geomh_cyv_8 (j)= cos( geomh_yv_8(j) )
         geomh_cyv2_8(j)= cos( geomh_yv_8(j) )**2
         geomh_cyM_8 (j)= geomh_cyv_8(j)
         geomh_invcy2_8 (j) = 1.d0/geomh_cy2_8  (j)
         geomh_invcyv2_8(j) = 1.d0/geomh_cyv2_8 (j)
         geomh_invcy_8  (j) = 1.d0/geomh_cy_8   (j)
         geomh_invcyv_8 (j) = 1.d0/geomh_cyv_8  (j)
      end do

      do j=1-G_haloy, l_nj+G_haloy
         indx = offj + j
         geomh_invDYMv_8(j)= geomh_invDY_8
         geomh_invDYM_8 (j)= geomh_invDY_8*geomh_invcy_8(j)
         geomh_tyoa_8   (j)= 0.d0
         geomh_tyoav_8  (j)= 0.d0
      end do

      do j=1-G_haloy, l_nj+G_haloy
         geomh_tyoa_8(j) = tan(geomh_y_8 (j)) * Dcst_inv_rayt_8
         geomh_tyoav_8(j)= tan(geomh_yv_8(j)) * Dcst_inv_rayt_8
      end do

      do j=1-G_haloy, l_nj+G_haloy
         scale_factor   = Dcst_rayt_8 * geomh_cy_8(j)
         scale_factor_v = Dcst_rayt_8 * geomh_cyv_8(j)
         geomh_invDX_8  (j) = 1.0d0/(scale_factor   * Del_xg )
         geomh_invDXv_8(j)  = 1.0d0/(scale_factor_v * Del_xg )
         geomh_invDXM_8 (j) = geomh_invDX_8(j)
         geomh_invDXMu_8(j) = geomh_invDX_8(j)
      end do

      do i=1,l_ni
         indx = offi + i
         xfi(i) = posx_8(indx)
      end do
      do i=1,l_nj
         indx = offj + i
         yfi(i) = posy_8(indx)
      end do
      gxfi(1:G_ni) = posx_8(1:G_ni)
      gyfi(1:G_nj) = posy_8(1:G_nj)

      Grd_global_gid = ezgdef_fmem (G_ni , G_nj , 'Z', 'E', Hgc_ig1ro, &
                       Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, gxfi(1), gyfi(1))
      Grd_local_gid  = ezgdef_fmem (l_ni , l_nj , 'Z', 'E', Hgc_ig1ro, &
                       Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,  xfi(1) ,  yfi(1))
      Grd_lclcore_gid= ezgdef_fmem (l_ni-pil_w-pil_e, l_nj-pil_s-pil_n,&
                       'Z', 'E', Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
                       xfi(1+pil_w), yfi(1+pil_s) )

      if ( .not. Ctrl_theoc_L) then
          nicore = G_ni - Glb_pil_w - Glb_pil_e
          njcore = G_nj - Glb_pil_s - Glb_pil_n

          Grd_glbcore_gid = ezgdef_fmem (nicore , njcore , 'Z', 'E', Hgc_ig1ro, &
           Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, gxfi(1+Glb_pil_w), gyfi(1+Glb_pil_s))
      end if

      offset= 2
      Grd_lphy_i0 = 1 + offset*west  ; Grd_lphy_in = l_ni-offset*east
      Grd_lphy_j0 = 1 + offset*south ; Grd_lphy_jn = l_nj-offset*north
      Grd_lphy_ni = Grd_lphy_in - Grd_lphy_i0 + 1
      Grd_lphy_nj = Grd_lphy_jn - Grd_lphy_j0 + 1

      Grd_lphy_gid  = ezgdef_fmem ( Grd_lphy_ni, Grd_lphy_nj, 'Z', 'E', &
                      Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro      , &
                      xfi(Grd_lphy_i0), yfi(Grd_lphy_j0) )

      if ( .not. Ctrl_theoc_L) then
        gphy_i0 = 1 + offset
        gphy_in = G_ni - offset
        gphy_j0 = 1 + offset
        gphy_jn = G_nj - offset
        gphy_ni = gphy_in - gphy_i0 + 1
        gphy_nj = gphy_jn - gphy_j0 + 1
        glbphy_gid = ezgdef_fmem (gphy_ni, gphy_nj,'Z', 'E', &
           Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, &
           gxfi(gphy_i0), gyfi(gphy_j0))

        gphy_i0 = 1 + max(offset, Glb_pil_w)
        gphy_j0 = 1 + max(offset, Glb_pil_s)
        gphy_nicore = gphy_ni - max(offset, Glb_pil_w) - max(offset, Glb_pil_e)
        gphy_njcore = gphy_nj - max(offset, Glb_pil_s) - max(offset, Glb_pil_n)
        glbphycore_gid = ezgdef_fmem(gphy_nicore, gphy_njcore, &
           'Z', 'E', Hgc_ig1ro, Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro,&
           gxfi(gphy_i0), gyfi(gphy_j0))

        istat= hgrid_wb_put ('model/Hgrid/glbcore', Grd_glbcore_gid, &
                            F_lni=nicore, F_lnj=njcore, F_rewrite_L=.true.)
        istat= hgrid_wb_put ('model/Hgrid/glbphy', glbphy_gid, &
                            F_lni=gphy_ni, F_lnj=gphy_nj, F_rewrite_L=.true.)
        istat= hgrid_wb_put ('model/Hgrid/glbphycore', glbphycore_gid, &
                            F_lni=gphy_nicore, F_lnj=gphy_njcore, F_rewrite_L=.true.)
      end if


      istat= hgrid_wb_put ('model/Hgrid/global', Grd_global_gid  , &
                            F_lni=G_ni, F_lnj=G_nj, F_rewrite_L=.true.)
      istat= hgrid_wb_put ('model/Hgrid/local', Grd_local_gid   , &
                            F_lni=l_ni, F_lnj=l_nj, F_rewrite_L=.true.)
      istat= hgrid_wb_put ('model/Hgrid/lclcore', Grd_lclcore_gid , &
                  F_lni=l_ni-pil_w-pil_e, F_lnj=l_nj-pil_s-pil_n , &
                  F_i0=1+pil_w, F_j0=1+pil_s, F_rewrite_L=.true.)
      istat= hgrid_wb_put ('model/Hgrid/lclphy' ,Grd_lphy_gid    , &
                  F_lni=Grd_lphy_ni, F_lnj=Grd_lphy_nj           , &
                  F_i0=Grd_lphy_i0, F_j0=Grd_lphy_j0, F_rewrite_L=.true.)


      hgc_array(1)= Hgc_ig1ro
      hgc_array(2)= Hgc_ig2ro
      hgc_array(3)= Hgc_ig3ro
      hgc_array(4)= Hgc_ig4ro

      istat= wb_put('model/Hgrid/hgcrot',hgc_array,WB_REWRITE_NONE+WB_IS_LOCAL)

      err = gdll (Grd_local_gid, geomh_latrx, geomh_lonrx)
      do j=1,l_nj
      do i=1,l_ni
         if (geomh_lonrx(i,j) >= 180.0) geomh_lonrx(i,j)=geomh_lonrx(i,j)-360.0
      end do
      end do

      do i=geomh_minx, geomh_maxx
         indx = offi + i
         xfi(i) = posx_8(indx)
      end do
      do i=geomh_miny, geomh_maxy
         indx = offj + i
         yfi(i) = posy_8(indx)
      end do
      ni = geomh_maxx-geomh_minx+1
      nj = geomh_maxy-geomh_miny+1
      dgid = ezgdef_fmem (ni , nj , 'Z', 'E', Hgc_ig1ro, &
                Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, xfi(geomh_minx) , yfi(geomh_miny) )
      err = gdll (dgid,geomh_latij,geomh_lonij)

 1000 format(/,'INITIALIZATION OF MODEL HORIZONTAL GEOMETRY (S/R set_geomh)', &
             /'===============================================')
 1001 format(/,'ERROR in (S/R set_geomh) Y AXIS in grid touch the poles', &
             /,'from:',F14.5,' to: ',F14.5, &
             /'======== Delta Y    is too large  =============')
!
!     ---------------------------------------------------------------
!
      return
      end subroutine set_geomh
