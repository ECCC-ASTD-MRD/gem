!--------------------------------- LICENCE BEGIN -------------------------------
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

!**s/r dcmip_vrd_set - Setting for vertical diffusion (Based on eqspng_set)

      subroutine dcmip_vrd_set()

      use step_options
      use dcmip_options
      use dcmip_vrd_coef
      use glb_ld
      use cstv
      use lun
      use ver
      use dynkernel_options

      implicit none

#include <arch_specific.hf>

      !object
      !=========================================================
      !     Setting for vertical diffusion (Based on eqspng_set)
      !=========================================================

      !--------------------------------------------------------------------------------

      real    :: Href2,Sponge,lnR_wd,lnR_tr,lnR_th,nuZeta_wd,nuZeta_th,nuZeta_tr,dZeta
      integer :: k,km,kp,istat

      !--------------------------------------------------------------------------------

      if (Dcmip_case==0) return

      Dcmip_wd_L = .false.
      Dcmip_th_L = .false.
      Dcmip_tr_L = .false.

      if (Dcmip_nuZ_wd/=0.) Dcmip_wd_L = .true.
      if (Dcmip_nuZ_th/=0.) Dcmip_th_L = .true.
      if (Dcmip_nuZ_tr/=0.) Dcmip_tr_L = .true.

      if (.NOT.Dcmip_wd_L.and..NOT.Dcmip_th_L.and..NOT.Dcmip_tr_L) return

      if (Dcmip_case/=163) call handle_error(-1,'DCMIP_VRD_SET','UREF etc. need to be prescribed')

      if (Lun_out > 0) write(Lun_out,1000)

      Dcmip_ref_wd = 0.
      Dcmip_ref_th = 0.
      Dcmip_ref_tr = 0.

      if (Dcmip_nuZ_wd<0.) Dcmip_ref_wd = 1.
      if (Dcmip_nuZ_th<0.) Dcmip_ref_th = 1.
      if (Dcmip_nuZ_tr<0.) Dcmip_ref_tr = 1.

      if (Lun_out>0) then
      if (Dcmip_nuZ_wd==0.) write(Lun_out,*) ' NO VERTICAL DIFFUSION WD'
      if (Dcmip_nuZ_th==0.) write(Lun_out,*) ' NO VERTICAL DIFFUSION TH'
      if (Dcmip_nuZ_tr==0.) write(Lun_out,*) ' NO VERTICAL DIFFUSION TR'
      if (Dcmip_nuZ_wd> 0.) write(Lun_out,*) '    VERTICAL DIFFUSION WD WITHOUT REF'
      if (Dcmip_nuZ_th> 0.) write(Lun_out,*) '    VERTICAL DIFFUSION TH WITHOUT REF'
      if (Dcmip_nuZ_tr> 0.) write(Lun_out,*) '    VERTICAL DIFFUSION TR WITHOUT REF'
      if (Dcmip_nuZ_wd< 0.) write(Lun_out,*) '    VERTICAL DIFFUSION WD WITH    REF'
      if (Dcmip_nuZ_th< 0.) write(Lun_out,*) '    VERTICAL DIFFUSION TH WITH    REF'
      if (Dcmip_nuZ_tr< 0.) write(Lun_out,*) '    VERTICAL DIFFUSION TR WITH    REF'
      end if

   !!!Href2   = (16000./log(10.))**2 ! (~6950 m)**2
      Href2   = (8780.2)**2          !DCMIP_height (as in SET_GEOM)
      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H') Href2 = 1.0d0

      nuZeta_wd = abs(Dcmip_nuZ_wd)/Href2
      nuZeta_th = abs(Dcmip_nuZ_th)/Href2
      nuZeta_tr = abs(Dcmip_nuZ_tr)/Href2

      dZeta = Ver_dz_8%m(G_nk/2)

      if (Dcmip_nuZ_wd/=0.) lnR_wd = nuZeta_wd*(2./Dzeta)**2*Cstv_dt_8
      if (Dcmip_nuZ_th/=0.) lnR_th = nuZeta_th*(2./Dzeta)**2*Cstv_dt_8
      if (Dcmip_nuZ_tr/=0.) lnR_tr = nuZeta_tr*(2./Dzeta)**2*Cstv_dt_8

      if (Lun_out>0.and.Dcmip_nuZ_wd/=0.) write(Lun_out,*) ' DCMIP lnR_wd = ',lnR_wd*100.,'% (Damping of the smallest wave)'
      if (Lun_out>0.and.Dcmip_nuZ_th/=0.) write(Lun_out,*) ' DCMIP lnR_th = ',lnR_th*100.,'% (Damping of the smallest wave)'
      if (Lun_out>0.and.Dcmip_nuZ_tr/=0.) write(Lun_out,*) ' DCMIP lnR_tr = ',lnR_tr*100.,'% (Damping of the smallest wave)'

      Sponge = max(abs(Dcmip_nuZ_wd),abs(Dcmip_nuZ_th),abs(Dcmip_nuZ_tr))

      !Check stability criteria
      !------------------------
      istat = 0
      do k= 1,G_nk
         if (Step_dt*Sponge/(Href2*Ver_dz_8%m(k)*Ver_dz_8%m(k))>.25) then
            istat = -1
            exit
         end if
      end do

      call handle_error(istat, 'DCMIP_VERTICAL_DIFFUSION', &
           'Selected diffusion coefficients making the scheme unstable, aborting')

      allocate ( dcmip_coef_m(1:G_nk+1),dcmip_coef_t(1:G_nk+1))

      !------------
      !Momentum U,V
      !------------
      allocate ( dcmip_cm_wd_m(G_nk), dcmip_cp_wd_m(G_nk) )

      dcmip_coef_m(     1)= 0.
      dcmip_coef_m(G_nk+1)= 0.
      dcmip_coef_m(2:G_nk)= abs(Dcmip_nuZ_wd)

      do k=1,G_nk
         km=max(1,k-1)
         dcmip_cp_wd_m(k) = Step_dt*dcmip_coef_m(k+1)/(Href2*Ver_dz_8%m(k)*Ver_dz_8%t(k ))
         dcmip_cm_wd_m(k) = Step_dt*dcmip_coef_m(k  )/(Href2*Ver_dz_8%m(k)*Ver_dz_8%t(km))
      end do

      !-----------
      !Thermo W,ZD
      !-----------
      allocate ( dcmip_cm_wd_t(G_nk), dcmip_cp_wd_t(G_nk) )

      dcmip_coef_t(     1)= 0.
      dcmip_coef_t(G_nk+1)= 0.
      dcmip_coef_t(2:G_nk)= abs(Dcmip_nuZ_wd)

      do k=1,G_nk
         kp=min(G_nk,k+1)
         dcmip_cp_wd_t(k) = Step_dt*dcmip_coef_t(k+1)/(Href2*Ver_dz_8%t(k)*Ver_dz_8%m(kp))
         dcmip_cm_wd_t(k) = Step_dt*dcmip_coef_t(k  )/(Href2*Ver_dz_8%t(k)*Ver_dz_8%m(k ))
      end do

      !------------
      !Thermo Theta
      !------------
      allocate ( dcmip_cm_th_t(G_nk), dcmip_cp_th_t(G_nk) )

      dcmip_coef_t(     1)= 0.
      dcmip_coef_t(G_nk+1)= 0.
      dcmip_coef_t(2:G_nk)= abs(Dcmip_nuZ_th)

      do k=1,G_nk
         kp=min(G_nk,k+1)
         dcmip_cp_th_t(k) = Step_dt*dcmip_coef_t(k+1)/(Href2*Ver_dz_8%t(k)*Ver_dz_8%m(kp))
         dcmip_cm_th_t(k) = Step_dt*dcmip_coef_t(k  )/(Href2*Ver_dz_8%t(k)*Ver_dz_8%m(k ))
      end do

      !--------------
      !Thermo Tracers
      !--------------
      allocate ( dcmip_cm_tr_t(G_nk), dcmip_cp_tr_t(G_nk) )

      dcmip_coef_t(     1)= 0.
      dcmip_coef_t(G_nk+1)= 0.
      dcmip_coef_t(2:G_nk)= abs(Dcmip_nuZ_tr)

      do k=1,G_nk
         kp=min(G_nk,k+1)
         dcmip_cp_tr_t(k) = Step_dt*dcmip_coef_t(k+1)/(Href2*Ver_dz_8%t(k)*Ver_dz_8%m(kp))
         dcmip_cm_tr_t(k) = Step_dt*dcmip_coef_t(k  )/(Href2*Ver_dz_8%t(k)*Ver_dz_8%m(k ))
      end do

      return

 1000 format( &
      /,'INITIALIZATING VERTICAL DIFFUSION: (S/R DCMIP_VRD_SET)', &
      /,'======================================================')

      end subroutine dcmip_vrd_set
