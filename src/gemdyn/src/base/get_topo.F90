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

!**s/r get_topo - Obtain topography from geophysical file

      subroutine get_topo ()
      use dynkernel_options
      use dyn_fisl_options
      use cstv
      use gem_options
      use geomh
      use glb_ld
      use lun
      use path
      use gmm_geof
      use gmm_pw
      use ptopo
      use tdpack
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      character(len=8) :: inttyp
      character(len=1024) :: fn
      integer :: stats, istat, zlist
      real, dimension(1-G_halox:G_ni+G_halox) :: xfi
      real, dimension(1-G_haloy:G_nj+G_haloy) :: yfi
      real, dimension(1-G_halox:G_ni+G_halox,1-G_haloy:G_nj+G_haloy) ::&
                                                  topo_destination, h1
      real(kind=REAL64) :: oneoRT
!
!-----------------------------------------------------------------------
!
      stats= 0 ; topo_high= 0. ; topo_low= 0. ; sls= 0.
      if (.not.Schm_topo_L) return

      if (.not.Lun_debug_L) istat= fstopc ('MSGLVL','SYSTEM',RMN_OPT_SET)
      zlist= -1
      inttyp      = 'LINEAR'
      
      If (Ptopo_myproc==0) Then
         xfi(1-G_halox:G_ni+G_halox)= &
                    real(geomh_longs(1-G_halox:G_ni+G_halox))
         yfi(1-G_haloy:G_nj+G_haloy)= &
                    real(geomh_latgs(1-G_haloy:G_nj+G_haloy))
         fn = trim(Path_input_S)//'/GEOPHY/Gem_geophy.fst'
         call get_field (topo_destination,G_ni+2*G_halox,G_nj+2*G_haloy,&
                         'ME', trim(fn), inttyp, xfi, yfi, stats)
         zlist=1
      end if
      call gem_error (stats,'GET_TOPO','Topography "ME" NOT specified')
      
      call glbdist_os (topo_destination, topo_high(l_minx,l_miny,1),&
                       l_minx,l_maxx,l_miny,l_maxy,1,&
                       G_ni+G_halox,G_nj+G_haloy,zlist,1,grav_8,0.d0)

      if (Schm_sleve_L) then
         if ( Schm_orols_fromgeophy_L ) then
            if (Ptopo_myproc==0) Then
               call get_field (h1, G_ni+2*G_halox, G_nj+2*G_haloy,&
                               'MELS', trim(fn), inttyp, xfi, yfi, stats)
            end if
            call gem_error (stats,'GET_TOPO',&
                           'Topography large scale "MESL" NOT specified')
            
            call glbdist_os (h1, topo_high(l_minx,l_miny,2),&
                             l_minx,l_maxx,l_miny,l_maxy,1 ,&
                G_ni+G_halox,G_nj+G_haloy,zlist,1,grav_8,0.d0)
         else
            call mc2_topols (topo_high(l_minx,l_miny,2),&
                             topo_high(l_minx,l_miny,1),&
                             l_minx,l_maxx,l_miny,l_maxy,Schm_orols_np)
         endif
         
         pw_p0_ls= 0.
         oneoRT=1.d0 / (rgasd_8 * Tcdk_8)
         if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' ) then
            pw_p0_ls(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy) = Cstv_pref_8 * exp(-topo_high(1-G_halox:l_ni+G_halox,1-G_haloy:l_nj+G_haloy,2) * oneoRT)
         endif
      end if
      
      topo_low = topo_high

      istat = fstopc ('MSGLVL','INFORM',RMN_OPT_SET)
!
!-----------------------------------------------------------------------
!
      return
      end subroutine get_topo
