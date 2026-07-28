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

!**s/r hzd_exp_set

      subroutine hzd_exp_set
      use hzd_mod
      use hvdif_options
      use tdpack
      use glb_ld
      use cstv
      use lun
      use ver
      use ptopo
      implicit none
#include <arch_specific.hf>

!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write(Lun_out,1002)

      hzd_hyb_top = hzd_hyb_lev(1)
      hzd_hyb_bot = hzd_hyb_lev(2)

      ! Constant z diffusion
      if ( Hzd_pwr_z < 0) Hzd_pwr_z = 2
      if ( Hzd_pwr_theta_z < 0) Hzd_pwr_theta_z = Hzd_pwr_z 
      ! Hybrid diffusion
      if (hzd_hyb_bot >0) then
         if ( Hzd_pwr_z < 0) Hzd_pwr_z = 2
         if ( Hzd_pwr_theta_z < 0) Hzd_pwr_theta_z = Hzd_pwr_z 

         if ( Hzd_pwr < 0) Hzd_pwr = Hzd_pwr_z
         if ( Hzd_pwr_theta < 0) Hzd_pwr_theta =  Hzd_pwr_theta_z

         if ( Hzd_lnr_theta < 0.) Hzd_lnr_theta = Hzd_lnr_theta_z
         if ( Hzd_lnr < 0.) Hzd_lnr = Hzd_lnr_z
      endif
 
      if(Hzd_lnr_theta_z >0. .OR. Hzd_lnr_z > 0.)  Hzd_alh_L=.true.

      if((Hzd_lnr_z > 0.).and.(Hzd_lnr_theta_z > 0.))then

        if(hzd_hyb_bot >0)then
           if (Lun_out > 0) then
              write(Lun_out,1005) Hzd_lnr_z,Hzd_lnr_theta_z
           end if
        else
           if (Lun_out > 0) then
              write(Lun_out,1006) Hzd_lnr_z,Hzd_lnr_theta_z
           end if
        end if
      end if

      Hzd_lnr = min(max(0.,Hzd_lnr),0.9999999)
      Hzd_lnr_z = min(max(0.,Hzd_lnr_z),0.9999999)
      Hzd_pwr = Hzd_pwr / 2
      Hzd_pwr = min(max(2,Hzd_pwr*2),8)

      Hzd_lnr_theta= min(max(0.,Hzd_lnr_theta),0.9999999)
      Hzd_lnr_theta_z= min(max(0.,Hzd_lnr_theta_z),0.9999999)
      Hzd_pwr_theta= Hzd_pwr_theta / 2
      Hzd_pwr_theta= min(max(2,Hzd_pwr_theta*2),8)

      if (Hzd_lnr_tr < 0.) Hzd_lnr_tr = Hzd_lnr
      if (Hzd_pwr_tr < 0 ) Hzd_pwr_tr = Hzd_pwr
      Hzd_lnr_tr = min(max(0.,Hzd_lnr_tr),0.9999999)
      Hzd_pwr_tr = Hzd_pwr_tr / 2
      Hzd_pwr_tr = min(max(2,Hzd_pwr_tr*2),8)

      if ((Hzd_lnr <= 0.).and.(Hzd_lnr_theta <= 0.)  &
                         .and.(Hzd_lnr_tr <= 0.)) then
         if((Hzd_smago_param <= 0.).and.(Hzd_smago_lnr(2) == 0.) &
                         .and.(Hzd_lnr_z <= 0).and.(Hzd_lnr_theta_z <= 0)) then
            if (Lun_out > 0) write(Lun_out,1003)
         elseif((Hzd_smago_param > 0.).and.(Hzd_smago_lnr(2) > 0.)) then
            if (Lun_out > 0) then
               write(Lun_out,1004) Hzd_smago_param,100*Hzd_smago_lnr(2)
            end if
         endif
      endif

      allocate(pres_pt (l_minx:l_maxx,l_miny:l_maxy,1:l_nk), &
               theta(l_minx:l_maxx,l_miny:l_maxy,1:l_nk),    & 
               theta0(l_minx:l_maxx,l_miny:l_maxy,1:l_nk),   &
               wk1(l_minx:l_maxx,l_miny:l_maxy,1:l_nk),      &
               wk2(l_minx:l_maxx,l_miny:l_maxy,1:l_nk),  &
               sfd(l_minx:l_maxx,l_miny:l_maxy,1:l_nk),      &
               sfd1(l_minx:l_maxx,l_miny:l_maxy,1:l_nk),  &
               wrkt1(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot) , &
               wrkd1(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot), &
               wrkt2(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top), & 
               wrkd2(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top), & 
               u_wrk(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot),   &
               v_wrk(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot),   &
               w_wrk(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot),   &
               zdt_wrk(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_bot), &
               u_wrk1(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top),  &
               v_wrk1(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top),  &
               w_wrk1(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top),  &
               zdt_wrk1(l_minx:l_maxx,l_miny:l_maxy,1:hzd_hyb_top) ,& 
               fdg2_4(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               Afdg1(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               Bfdg1(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               Afdg2(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               Bfdg2(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               bdd_v8(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               add_v8(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               cdd_v8(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               cdd_v82(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               cdd_v81(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               cflux(l_minx:l_maxx, l_miny:l_maxy,l_nk+1),&
               stencilV(l_minx:l_maxx, l_miny:l_maxy,3,l_nk),&
       	       a_th(l_minx:l_maxx, l_miny:l_maxy,l_nk) ,&
       	       b_th(l_minx:l_maxx, l_miny:l_maxy,l_nk) ,&
       	       d_th(l_minx:l_maxx, l_miny:l_maxy,l_nk))

      pres_pt=0. ; theta=0.; theta0=0. ; wk1=0. ; wk2=0. 
      wrkt1=0. ; wrkd1=0. ; wrkt2=0.; wrkd2=0.
      u_wrk=0.  ; v_wrk  =0. ; w_wrk  =0. ; zdt_wrk  =0.
      u_wrk1=0. ; v_wrk1 =0. ; w_wrk1 =0. ; zdt_wrk1 =0.
      sfd=0. ; sfd1=0.;fdg2_4=0.;
      Afdg1=0.d0; Bfdg1=0.d0; Afdg2=0.d0; Bfdg2=0.d0
      bdd_v8=0.d0; add_v8=0.d0; cdd_v8=0.d0; cdd_v82=0.d0; cdd_v81=0.d0
      a_th=0.d0; b_th=0.d0; d_th=0.d0;
      stencilV=0.d0; cflux=0.d0

      call hzd_exp_geom ()

      call hzd_exp5p_set ()

 1002 format(/,'INITIALIZATING HIGH ORDER HORIZONTAL DIFFUSION ',  &
               '(S/R HZD_SET)',/,60('='))
 1003 format(/,'NO HORIZONTAL DIFFUSION REQUESTED',/,33('='))
 1004 format(/,'  HORIZONTAL DIFFUSION A LA SMAGORINSKY',/,2x,37('=')// &
              ,'  PARAMETER =',f5.2,'  BACKGROUND =',f4.1,' %/TIMESTEP')
 1005 format(/,'  HORIZONTAL HYBRID  DIFFUSION ',/,2x,37('=')// &
              ,'  HZD_LNR_Z =',f6.3,'  HZD_LNR_THETA_Z =',f6.3)
 1006 format(/,'  HORIZONTAL DIFFUSION ALONG CONSTANT Z',/,2x,37('=')// &
             ,'  HZD_LNR_Z =',f5.2,'  HZD_LNR_THETA_Z =',f5.1)
!
!     ---------------------------------------------------------------
!
      return
      end
