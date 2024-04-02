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
      use HORgrid_options
      use spn_options
      use glb_ld
      use glb_pil
      use cstv
      use lun
      use ver
      use ldnh
      use step_options
      use tdpack
      use ptopo
      implicit none
#include <arch_specific.hf>

      logical, external :: decomp
      integer, parameter :: lowest = 2
      integer :: minx, maxx, npartiel, n0
      integer :: k, err1, err2, tmdt
      integer, dimension(Ptopo_npex) :: lnis      
      real    :: t_turn, b_turn
!
!----------------------------------------------------------------------
!
      if ( Grd_yinyang_L .or. (Spn_freq<=0) ) return

      err1= 0 ; err2= 0

      call low2up( Spn_trans_shape_S, Spn_trans_shape_S )

      if (Lun_out > 0) write(Lun_out,1000)
      
      if (Lun_out > 0) write(Lun_out,1002) ' Transpose 1===>2 for SPN:', &
                 ' G_nk distributed on Ptopo_npex PEs', G_nk,Ptopo_npex
      
      if (.not. decomp (G_nk, minx, maxx, lnis, npartiel, 0, n0, &
                .true., .true., Ptopo_npex, -1, .false., 3 )) err1 = -1
      
      Spn_12smin = minx
      Spn_12smax = maxx
      Spn_12sn   = lnis(1)
      Spn_12sn0  = n0
      
      if (Lun_out > 0) write(Lun_out,1002) ' Transpose 2===>2 for SPN:', &
                 ' G_ni distributed on Ptopo_npey PEs', G_ni,Ptopo_npey

      if (.not. decomp (G_ni, minx, maxx, lnis, npartiel, 0, n0, &
                .false., .true., Ptopo_npex, lowest, .false., 0 )) err1 = -1

      Spn_22min = minx
      Spn_22max = maxx
      Spn_22n   = lnis(1)
      Spn_22n0  = n0

      Spn_22pil_w= 0 ;  Spn_22pil_e= 0
      if (Spn_22n0==1)              Spn_22pil_w= Grd_extension
      if (Spn_22n0+Spn_22n-1==G_ni) Spn_22pil_e= Grd_extension

      if (err1<0) goto 999
      
      allocate ( prof(G_nk), &
         Spn_fft(ldnh_maxy ,Spn_12smax,G_ni+2+Ptopo_npex),&
         Spn_fdg(Spn_12smax,Spn_22max ,G_nj  +Ptopo_npey),&
         Spn_wrk(ldnh_maxx,ldnh_maxy,l_nk) )
      prof=0. ; Spn_fft= 0. ; Spn_fdg= 0. ; Spn_wrk= 0.
     
      Spn_njnh  = ldnh_maxy-ldnh_miny+1
      Spn_nk12  = Spn_12smax-Spn_12smin+1
      Spn_ni22  = Spn_22max-Spn_22min+1
      tmdt      = int(Cstv_dt_8)
      Spn_interval = Spn_freq/tmdt
      Spn_interval = max(1,Spn_interval)
      Spn_ws = Step_nesdt/tmdt
      Spn_weight= 1.0
      
      call spn_calfiltre ()

      t_turn= max( Spn_up_const_lev,Ver_hyb%m(  1 ) )
      b_turn= min( Spn_start_lev   ,Ver_hyb%m(G_nk) )

      if (Spn_trans_shape_S == 'COS2' ) then

         do k=1,G_nk
            if (Ver_hyb%m(k) <= b_turn .and. Ver_hyb%m(k) >= t_turn) then
               prof(k) = cos(pi_8-pi_8*(b_turn-Ver_hyb%m(k))/(b_turn-t_turn))
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
         err2 = -1
      end if
      
 999  call gem_error ( min(err1,err2),'spn_init',&
          'Wrong choice for Spn_trans_shape_S or transpose problems')

      do k=1,G_nk
         prof(k) = prof(k) * Cstv_dt_8/(Spn_relax_hours*3600.)
      end do

      Spn_ON_L= .true.

 1000 format (/' SPN_INIT: Initialization of spectral nudging'/)
 1002 format (a/a45,i6,' /',i5)
!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_init
