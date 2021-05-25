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

!**s/r init_bar - prepare data for autobarotropic runs (Williamson cases)

      subroutine init_bar ()

      use cstv
      use ctrl
      use dyn_fisl_options
      use dynkernel_options
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use inp_mod
      use mem_tracers
      use step_options
      use tr3d
      use wil_options

      implicit none

#include <arch_specific.hf>

      !object
      !============================================================
      !Prepare data for autobarotropic runs (Williamson cases)
      !------------------------------------------------------------
      !NOTE: U,V output on Staggered grids
      !============================================================

      integer :: k,i0,in,j0,jn
      real, dimension (:,:,:), pointer :: hu
      real, dimension (l_minx:l_maxx,l_miny:l_maxy,G_nk) :: gz_t
!
!-------------------------------------------------------------------
!
      if (Vtopo_L)      call gem_error (-1,'INIT_BAR','Vtopo_L not available YET')

      if (Schm_sleve_L) call gem_error (-1,'INIT_BAR','  SLEVE not available YET')

      orols = 0. ; sls = 0. ; topo_low = 0. ; topo_high = 0.

      i0= 1-G_halox ; in= l_ni+G_halox
      j0= 1-G_haloy ; jn= l_nj+G_haloy

      !Setup Williamson Case 7: The 21 December 1978 Initial conditions are read
      !-------------------------------------------------------------------------
      if (Williamson_case==7) then

         call inp_data ( ut1, vt1, wt1, tt1, qt1, zdt1,st1,trt1,&
                         fis0, orols, .true., Step_runstrt_S   ,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr)

         if (Schm_sleve_L) then
            call update_sls (orols,sls,l_minx,l_maxx,l_miny,l_maxy)
         endif

         !Initialize log(surface pressure)
         !--------------------------------
         if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then

            st1(i0:in,j0:jn) = log(st1(i0:in,j0:jn)/Cstv_pref_8)

         !Initialize PHI perturbation in q
         !--------------------------------
         else if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') then

            do k=1,G_nk+1
               qt1(i0:in,j0:jn,k) = st1(i0:in,j0:jn) - 1.0d0/Cstv_invFI_8
            end do

         end if

      end if

      !Initialize T/ZD/W/Q
      !-------------------
      tt1 = Cstv_Tstr_8 ; zdt1 = 0. ; wt1 = 0.

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') qt1 = 0.

      !Initialize HU
      !-------------
      hu => tracers_P(Tr3d_hu)%pntr

      hu = 0.

      !Initialize d(Zeta)dot and dz/dt
      !-------------------------------
      Inp_zd_L = .true.
      Inp_w_L  = .true.

      !Prepare initial conditions (staggered u-v,gz,s,topo) for Williamson cases
      !-------------------------------------------------------------------------
      call wil_init (ut1,vt1,gz_t,st1,fis0,qt1,l_minx,l_maxx,l_miny,l_maxy,G_nk,.true.)

      !Required for CASE5 LAM version (NESTING)
      !----------------------------------------
      topo_high(i0:in,j0:jn,1) =      fis0(i0:in,j0:jn)
      topo_low (i0:in,j0:jn,1) = topo_high(i0:in,j0:jn,1)
      topo_high(i0:in,j0:jn,2) =     orols(i0:in,j0:jn)
      topo_low (i0:in,j0:jn,2) = topo_high(i0:in,j0:jn,2)

      if (Schm_sleve_L) then
         call update_sls (orols,sls,l_minx,l_maxx,l_miny,l_maxy)
      endif

      !Estimate U-V and T on scalar grids
      !----------------------------------
      call hwnd_stag2( pw_uu_plus,pw_vv_plus,ut1,vt1      ,&
                       l_minx,l_maxx,l_miny,l_maxy,G_nk   ,&
                       1-G_halox*west ,l_niu+G_halox*east ,&
                       1-G_haloy*south,l_njv+G_haloy*north, .false. )

      pw_tt_plus = tt1

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H' .and. .not.Ctrl_testcases_adv_L) then
         do k=1,G_nk
            qt1(i0:in,j0:jn, k) = max(gz_t(i0:in,j0:jn,1), 0.)
         end do
      end if
!
!-------------------------------------------------------------------
!
      return
      end subroutine init_bar
