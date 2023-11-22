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

module dcmip_vrd_coef
   implicit none
   public
   save

   !-------------------------------------------------------------------------------------------------------
   !Variables used in dcmip_vrd (DCMIP vertical diffusion)
   !-------------------------------------------------------------------------------------------------------
   logical dcmip_wd_L,dcmip_th_L,dcmip_tr_L

   real dcmip_ref_wd,dcmip_ref_th,dcmip_ref_tr

   real, dimension(:  ), pointer :: dcmip_coef_m,dcmip_cm_wd_m,dcmip_cp_wd_m,dcmip_cm_tr_t,dcmip_cp_tr_t, &
                                    dcmip_coef_t,dcmip_cm_wd_t,dcmip_cp_wd_t,dcmip_cm_th_t,dcmip_cp_th_t

end module dcmip_vrd_coef
