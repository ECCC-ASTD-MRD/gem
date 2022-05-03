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

!**s/r nest_indata - Read and process nesting data during LAM
!                    integration for LBC.

      subroutine nest_indata ( F_u , F_v, F_w, F_t , F_q, F_zd, F_s, &
                               F_tr, F_topo, F_stag_L, F_datev_S    ,&
                               Mminx, Mmaxx, Mminy, Mmaxy, Nk, Ntr )
      use gmm_geof
      use dynkernel_options
      use dyn_fisl_options
      use gem_options
      use inp_mod
      use mem_nest
      use glb_ld
      use tr3d
      implicit none
      
      character(len=*), intent(in):: F_datev_S
      logical, intent(in) :: F_stag_L
      integer, intent(in) :: Mminx, Mmaxx, Mminy, Mmaxy, Nk, Ntr
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk),   intent(out) :: F_u, F_v, F_w, F_t, F_zd
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk+1), intent(out) :: F_q
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy),      intent(out) :: F_s
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,2),    intent(out) :: F_topo
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk*Ntr),intent(out) :: F_tr
      
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: uu, vv, tt, sumpqj
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: nest_sls
!     
!     ---------------------------------------------------------------
!
      call inp_data (F_u , F_v, F_w, F_t , F_q, F_zd, F_s, F_tr      ,&
                     F_topo(l_minx,l_miny,1),F_topo(l_minx,l_miny,2),&
              F_stag_L,F_datev_S,l_minx,l_maxx,l_miny,l_maxy,G_nk,Ntr)

      if (Schm_sleve_L) then
         call update_sls (F_topo(l_minx,l_miny,2),nest_sls,&
                          l_minx,l_maxx,l_miny,l_maxy)
      else
         nest_sls= 0.
      endif

      if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' )  &
      call vertical_metric (nest_metric, F_topo, nest_sls, &
                                l_minx,l_maxx,l_miny,l_maxy)
                   
      call canonical_indata()
      
      if (.not.F_stag_L) then
         uu= F_u ; vv= F_v ; tt= F_t ; sumpqj= 0.
         call mfottvh2 (tt, F_t, F_tr(l_minx,l_minx,(Tr3d_hu-1)*l_nk+1),&
                        sumpqj,l_minx, l_maxx, l_miny, l_maxy, l_nk    ,&
                   1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy,.true.)
         call hwnd_stag2 ( F_u, F_v,uu,vv,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk   ,&
                         1-G_halox*west ,l_niu+G_halox*east ,&
                         1-G_haloy*south,l_njv+G_haloy*north, .true. )
      endif
                       
      call derivate_data ( F_zd, F_w, F_u, F_v, F_t , F_s, F_q         ,&
                           F_topo(l_minx,l_miny,1),nest_sls,nest_metric,&
                           l_minx,l_maxx,l_miny,l_maxy, G_nk           ,&
                           .not.Inp_zd_L, .not.Inp_w_L, .not. Inp_qt_L )
!
!     ---------------------------------------------------------------
!
      return
      end subroutine nest_indata

