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
!---------------------------------- LICENCE END ----------------------

!*s/r spn_init - initialize spectral nudging profile, filter

      subroutine spn_init()
      use spn_work_mod
      use HORgrid_options
      use spn_options
      use glb_ld
      use cstv
      use lun
      use ver
      use ldnh
      use glb_pil
      implicit none
#include <arch_specific.hf>

!author
!     Minwei Qian (CCRD) & Bernard Dugas, Syed Husain  (MRB)  - summer 2015
!
!revision
! v4_80 - Qian, Dugas, Hussain            - initial version
! v4_80 - Baek - correction for NK


      integer :: k,err1
      real :: t_turn, b_turn, pi2, nudging_tau
!
!----------------------------------------------------------------------
!
      if ( Grd_yinyang_L .or. (Spn_nudging_S == ' ') ) return

      err1 = 0

      pi2 = atan( 1.d0 )*2.d0

      call up2low( Spn_nudging_S, Spn_nudging_S )
      call low2up( Spn_trans_shape_S, Spn_trans_shape_S )

      if (Lun_out > 0) write(Lun_out,1000) Spn_nudging_S

      ! Allocate once and for all

      allocate( prof(G_nk), &
                fxy(G_ni+2,G_nj+2), &
                Ldiff3D(ldnh_minx:ldnh_maxx,ldnh_miny:ldnh_maxy,G_nk), stat=err1 )

      if (err1 > 0) call stop_mpi(-1,'spn_init','Error in spn_init')

      call spn_calfiltre (G_ni-Lam_pil_w-Lam_pil_e, &
                          G_nj-Lam_pil_s-Lam_pil_n)

      prof=0.

      ! nudging_tau is nudging time scale in hours
      ! t_turn and b_turn are top and bottom turnning points

      t_turn= max( Spn_up_const_lev,Ver_hyb%m(  1 ) )
      b_turn= min( Spn_start_lev   ,Ver_hyb%m(G_nk) )
      nudging_tau = Spn_relax_hours

      if (Spn_trans_shape_S == 'COS2' ) then

         do k=1,G_nk
            if (Ver_hyb%m(k) <= b_turn .and. Ver_hyb%m(k) >= t_turn) then
               prof(k) = cos(pi2-pi2*(b_turn-Ver_hyb%m(k))/(b_turn-t_turn))
            elseif (Ver_hyb%m(k) < t_turn) then
               prof(k)=1.
            else
               prof(k)=0.
            end if
            prof(k) = prof(k)*prof(k)
         end do

      elseif (Spn_trans_shape_S == 'LINEAR' ) then

         do k=1,G_nk
            if (Ver_hyb%m(k) <= b_turn .and. Ver_hyb%m(k) >= t_turn) then
               prof(k) =  (b_turn-Ver_hyb%m(k))/(b_turn-t_turn)
            elseif (Ver_hyb%m(k) < t_turn) then
               prof(k)=1.
            else
               prof(k)=0.
            end if
         end do

      else

         if (Lun_out > 0) write(Lun_out,1001) Spn_trans_shape_S
         call stop_mpi(-1,'spn_init','Error in spn_init')

      end if

      do k=1,G_nk
         prof(k) = prof(k) * Cstv_dt_8/3600./nudging_tau
      end do

 1000 format(/' In SPN_INIT, Spn_nudging_S = ',A8/)
 1001 format(/' In SPN_INIT, unknown Spn_trans_shape_S ',A8/)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_init
