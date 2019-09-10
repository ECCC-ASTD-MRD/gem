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
!**s/r nest_init -- Initializes nesting data for LAM configuration

      subroutine nest_init ()
      use lam_options
      use gmm_vt1
      use gmm_nest
      use gmm_geof
      use lun
      use gmm_itf_mod
      use step_options
      use tr3d
      implicit none
#include <arch_specific.hf>

      character(len=GMM_MAXNAMELENGTH) :: tr_name
      integer n,yy,mo,dd,hh,mm,ss,dum,istat
      real, pointer, dimension(:,:,:) :: tr1,trf
!
!     ---------------------------------------------------------------
!
      if (Lun_debug_L) write(Lun_out,1000)

      istat = gmm_get(gmmk_ut1_s ,ut1 )
      istat = gmm_get(gmmk_vt1_s ,vt1 )
      istat = gmm_get(gmmk_tt1_s ,tt1 )
      istat = gmm_get(gmmk_st1_s ,st1 )
      istat = gmm_get(gmmk_wt1_s ,wt1 )
      istat = gmm_get(gmmk_qt1_s ,qt1 )
      istat = gmm_get(gmmk_zdt1_s,zdt1)
      istat = gmm_get(gmmk_fis0_s,fis0)

!     copying values from UT1 to nest_u variables

      if (Lam_ctebcs_L) then
!     LAM with same (constant) pilot conditions
         istat = gmm_get(gmmk_nest_u_s ,nest_u )
         istat = gmm_get(gmmk_nest_v_s ,nest_v )
         istat = gmm_get(gmmk_nest_t_s ,nest_t )
         istat = gmm_get(gmmk_nest_s_s ,nest_s )
         istat = gmm_get(gmmk_nest_w_s ,nest_w )
         istat = gmm_get(gmmk_nest_q_s ,nest_q )
         istat = gmm_get(gmmk_nest_zd_s,nest_zd)
         istat = gmm_get(gmmk_nest_fullme_s,nest_fullme)
         nest_u  = ut1
         nest_v  = vt1
         nest_t  = tt1
         nest_s  = st1
         nest_w  = wt1
         nest_q  = qt1
         nest_zd = zdt1
         nest_fullme = fis0

         do n=1,Tr3d_ntr
            tr_name = 'TR/'//trim(Tr3d_name_S(n))//':P'
            istat = gmm_get(tr_name,tr1)
            tr_name = 'NEST/'//trim(Tr3d_name_S(n))//':C'
            istat = gmm_get(tr_name,trf)
            trf = tr1
         end do

      else
!     ordinary LAM with future pilot conditions
         istat = gmm_get(gmmk_nest_u_fin_s ,nest_u_fin )
         istat = gmm_get(gmmk_nest_v_fin_s ,nest_v_fin )
         istat = gmm_get(gmmk_nest_t_fin_s ,nest_t_fin )
         istat = gmm_get(gmmk_nest_s_fin_s ,nest_s_fin )
         istat = gmm_get(gmmk_nest_w_fin_s ,nest_w_fin )
         istat = gmm_get(gmmk_nest_q_fin_s ,nest_q_fin )
         istat = gmm_get(gmmk_nest_zd_fin_s,nest_zd_fin)
         istat = gmm_get(gmmk_nest_fullme_fin_s,nest_fullme_fin)
         nest_u_fin  = ut1
         nest_v_fin  = vt1
         nest_t_fin  = tt1
         nest_s_fin  = st1
         nest_w_fin  = wt1
         nest_q_fin  = qt1
         nest_zd_fin = zdt1
         nest_fullme_fin = fis0
         do n=1,Tr3d_ntr
            tr_name = 'TR/'//trim(Tr3d_name_S(n))//':P'
            istat = gmm_get(tr_name,tr1)
            tr_name = 'NEST/'//trim(Tr3d_name_S(n))//':F'
            istat = gmm_get(tr_name,trf)
            trf = tr1
         end do

      end if

      Lam_previous_S= ''
      Lam_current_S = Step_runstrt_S
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,Lam_current_S)
      call pdfjdate2 (Lam_tfin, yy,mo,dd,hh,mm,ss)
      Lam_tdeb      = Lam_tfin
!
!     ---------------------------------------------------------------
!
 1000 format(3X,'NESTING INITIALIZATION (NEST_INIT)')
      return
      end
