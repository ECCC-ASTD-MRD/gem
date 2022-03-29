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

module mem_nest
   use metric
   implicit none
   public
   save

      real, pointer, dimension (:,:,:) :: nest_u_deb  => null()
      real, pointer, dimension (:,:,:) :: nest_v_deb  => null()
      real, pointer, dimension (:,:,:) :: nest_t_deb  => null()
      real, pointer, dimension (:,:  ) :: nest_s_deb  => null()
      real, pointer, dimension (:,:,:) :: nest_w_deb  => null()
      real, pointer, dimension (:,:,:) :: nest_q_deb  => null()
      real, pointer, dimension (:,:,:) :: nest_zd_deb => null()
      real, pointer, dimension (:,:,:) :: nest_tr_deb => null()
      real, pointer, dimension (:,:,:) :: nest_fullme_deb => null()

      real, pointer, dimension (:,:,:) :: nest_u      => null()
      real, pointer, dimension (:,:,:) :: nest_v      => null()
      real, pointer, dimension (:,:,:) :: nest_t      => null()
      real, pointer, dimension (:,:  ) :: nest_s      => null()
      real, pointer, dimension (:,:,:) :: nest_w      => null()
      real, pointer, dimension (:,:,:) :: nest_q      => null()
      real, pointer, dimension (:,:,:) :: nest_zd     => null()
      real, pointer, dimension (:,:,:) :: nest_tr     => null()
      real, pointer, dimension (:,:,:) :: nest_fullme => null()

      real, pointer, dimension (:,:,:) :: nest_u_fin  => null()
      real, pointer, dimension (:,:,:) :: nest_v_fin  => null()
      real, pointer, dimension (:,:,:) :: nest_t_fin  => null()
      real, pointer, dimension (:,:  ) :: nest_s_fin  => null()
      real, pointer, dimension (:,:,:) :: nest_w_fin  => null()
      real, pointer, dimension (:,:,:) :: nest_q_fin  => null()
      real, pointer, dimension (:,:,:) :: nest_zd_fin => null()
      real, pointer, dimension (:,:,:) :: nest_tr_fin => null()
      real, pointer, dimension (:,:,:) :: nest_fullme_fin => null()

      real, pointer, dimension (:,:,:) :: nest_weightm=> null()
      real, pointer, dimension (:,:,:) :: nest_weightq=> null()
      real, pointer, dimension (:,:,:) :: nest_weightu=> null()
      real, pointer, dimension (:,:,:) :: nest_weightv=> null()

      real, pointer, dimension (:) :: nest_deb, nest_now,&
                                      nest_fin

      type(Vmetric) :: nest_metric
end module mem_nest
