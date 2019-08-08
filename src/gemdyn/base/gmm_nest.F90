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

module gmm_nest
   implicit none
   public
   save
!
! v4_10 - Tanguay M.        - Adjust digital filter when LAM
!
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH NESTING for current timestep              |
!  For Nest_uf,Nest_vf - used for future timesteps with one extra level|
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Nest_u             | x component of velocity                         |
! Nest_v             | y component of velocity                         |
! Nest_t             | T (temperature)                                 |
! Nest_zd            |                                                 |
! Nest_s             | ln (dpi/dpi*)                                   |
!--------------------|-------------------------------------------------|
! Nest_w             | z component of velocity                         |
! Nest_q             |                                                 |
!--------------------|-------------------------------------------------|
! Nest_tr            | tracer 3d variables                             |
!----------------------------------------------------------------------|
!

      real, pointer, contiguous, dimension (:,:,:) :: nest_u_deb  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_v_deb  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_t_deb  => null()
      real, pointer, contiguous, dimension (:,:  ) :: nest_s_deb  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_w_deb  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_q_deb  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_zd_deb => null()
      real, pointer, contiguous, dimension (:,:  ) :: nest_fullme_deb => null()

      real, pointer, contiguous, dimension (:,:,:) :: nest_u      => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_v      => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_t      => null()
      real, pointer, contiguous, dimension (:,:  ) :: nest_s      => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_w      => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_q      => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_zd     => null()
      real, pointer, contiguous, dimension (:,:  ) :: nest_fullme     => null()

      real, pointer, contiguous, dimension (:,:,:) :: nest_u_fin  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_v_fin  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_t_fin  => null()
      real, pointer, contiguous, dimension (:,:  ) :: nest_s_fin  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_w_fin  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_q_fin  => null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_zd_fin => null()
      real, pointer, contiguous, dimension (:,:  ) :: nest_fullme_fin => null()

      real, pointer, contiguous, dimension (:,:,:) :: nest_weightm=> null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_weightq=> null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_weightu=> null()
      real, pointer, contiguous, dimension (:,:,:) :: nest_weightv=> null()

      integer, parameter :: MAXNAMELENGTH = 32

      character(len=MAXNAMELENGTH) :: &
               gmmk_nest_u_deb_s , gmmk_nest_v_deb_s , gmmk_nest_t_deb_s ,&
               gmmk_nest_s_deb_s , gmmk_nest_w_deb_s , gmmk_nest_q_deb_s ,&
               gmmk_nest_zd_deb_s, gmmk_nest_fullme_deb_s                ,&
               gmmk_nest_u_s , gmmk_nest_v_s , gmmk_nest_t_s ,&
               gmmk_nest_s_s , gmmk_nest_w_s , gmmk_nest_q_s ,&
               gmmk_nest_zd_s, gmmk_nest_fullme_s            ,&
               gmmk_nest_u_fin_s , gmmk_nest_v_fin_s , gmmk_nest_t_fin_s ,&
               gmmk_nest_s_fin_s , gmmk_nest_w_fin_s , gmmk_nest_q_fin_s ,&
               gmmk_nest_zd_fin_s, gmmk_nest_fullme_fin_s                ,&
               gmmk_nest_weightm_s, gmmk_nest_weightu_s, &
               gmmk_nest_weightv_s, gmmk_nest_weightq_s

end module gmm_nest
