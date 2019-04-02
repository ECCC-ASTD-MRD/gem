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

!**s/r eqspng_set

       subroutine eqspng_set

       use step_options
       use hvdif_options
       use geomh
       use glb_ld
       use ver
       use dynkernel_options
       implicit none
#include <arch_specific.hf>

       real    :: Href2
       real*8  :: pdb, pdtmp, pda,  pdc
       integer :: i,j,k,km,istat
!
!     ---------------------------------------------------------------
!
       eq_nlev = 0
       Href2   = (16000./log(10.))**2 ! (~6950 m)**2
       if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' ) Href2 = 1.0d0

! Check stability criteria and count level to diffuse.

       istat = 0
       do k= 1, size(Eq_sponge)
          if ( Eq_sponge(k) > epsilon(Eq_sponge(k)) )then
             if (Step_dt*Eq_sponge(k)/(Href2*Ver_dz_8%m(k)*Ver_dz_8%m(k))>.25)then
                istat = -1
                exit
             end if
             eq_nlev= eq_nlev+1
          else
             exit
          end if
       end do

       call handle_error(istat, 'equatorial_sponge', &
            'Selected diffusion coeficients making the scheme unsatable, aborting')

       if (eq_nlev <= 0) return

       eq_nlev= eq_nlev+1

       allocate ( coef(1:eq_nlev+1), cm(eq_nlev), cp(eq_nlev) )

       coef(        1)= 0.
       coef(eq_nlev+1)= 0.
       coef(2:eq_nlev)= Eq_sponge(1:eq_nlev-1)

       do k=1,eq_nlev
          km=max(1,k-1)
          cp(k)=Step_dt*coef(k+1)/(Href2*Ver_dz_8%m(k)*Ver_dz_8%t(k ))
          cm(k)=Step_dt*coef(k  )/(Href2*Ver_dz_8%m(k)*Ver_dz_8%t(km))
       end do

       allocate ( eponmod(l_ni,l_nj) )

       if (Eq_ramp_L) then
          pdb = P_lmvd_high_lat - P_lmvd_low_lat
          do j= 1, l_nj
            do i= 1, l_ni
               pdtmp = abs(geomh_latrx(i,j))
               pda   = min(1.0d0,max(0.0d0,(pdtmp-P_lmvd_low_lat)/pdb))
               pdc   = (3.-2.*pda)*pda*pda
               eponmod(i,j) = pdc  * P_lmvd_weigh_high_lat + &
                              (1. - pdc) * P_lmvd_weigh_low_lat
               eponmod(i,j) = max(0.0,min(1.0,eponmod(i,j)))
            end do
          end do
       else
          eponmod = 1.
       end if
!
!     ---------------------------------------------------------------
!
       return
       end subroutine eqspng_set
