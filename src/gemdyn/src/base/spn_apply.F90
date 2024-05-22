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
!---------------------------------- LICENCE END --------------------------------

!*s/r spn_fld - 2D forward FFT + filter + backward FFT 
!               and apply nudging tendency

      subroutine spn_apply ( F_ft1, F_nest, Minx, Maxx, Miny, Maxy, Nk )
      use spn_options
      use gmm_vt1
      use glb_ld
      use glb_pil
      use HORgrid_options
      use ldnh
      use, intrinsic :: iso_fortran_env
      implicit none

      integer, intent(in) :: Minx, Maxx, Miny, Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(INOUT) :: F_ft1
      real, dimension(Minx:Maxx,Miny:Maxy,  Nk), intent(IN   ) :: F_nest

      integer i,j,k
      real(kind=REAL64) pri
!
!----------------------------------------------------------------------
!      
      Spn_wrk(1:l_ni,1:l_nj,1:l_nk) = F_nest(1:l_ni,1:l_nj,1:l_nk) - F_ft1(1:l_ni,1:l_nj,1:l_nk)
      
      call rpn_comm_transpose ( Spn_wrk, ldnh_minx, ldnh_maxx, G_ni, Spn_njnh,&
                                Spn_12smin, Spn_12smax, G_nk, Spn_fft, 1, 2 )

      call itf_fft_set( G_ni-2*Grd_extension,'QCOS',pri )

      do k= 1, Spn_12sn                ! forward fft in X direction
         call itf_fft_drv( Spn_fft(1+pil_s,k,1+Lam_pil_w), &
              Spn_njnh*Spn_nk12,1,(Spn_njnh-pil_s-pil_n), -1 )
      end do

      call rpn_comm_transpose ( Spn_fft, ldnh_miny, ldnh_maxy, G_nj, Spn_nk12,&
                                Spn_22min, Spn_22max, G_ni, Spn_fdg,  2, 2 )
      
      call itf_fft_set( G_nj-2*Grd_extension,'QCOS',pri )

      do k=1, Spn_22n             ! forward fft in Y direction
         call itf_fft_drv( Spn_fdg(1,k,1+Lam_pil_s), &
                     Spn_nk12*Spn_ni22,1, Spn_nk12, -1 )
      end do

      do j= 1, G_nj            ! filter in X-Y direction
         do i= 1, Spn_22n
            do k= 1, Spn_12sn
               Spn_fdg(k,i,j)=Spn_fdg(k,i,j)*Spn_flt(i,j)
            end do
         end do
      end do

      do k=1, Spn_22n           ! backward fft in Y direction
         call itf_fft_drv( Spn_fdg(1,k,1+Lam_pil_s), &
                    Spn_nk12*Spn_ni22,1, Spn_nk12, +1 )
      end do

      call rpn_comm_transpose ( Spn_fft, ldnh_miny, ldnh_maxy, G_nj, Spn_nk12,&
                                Spn_22min, Spn_22max, G_ni, Spn_fdg, -2, 2 )
      
      call itf_fft_set( G_ni-2*Grd_extension,'QCOS',pri )

      do k= 1, Spn_12sn    ! backward fft in X direction
         call itf_fft_drv( Spn_fft(1+pil_s,k,1+Lam_pil_w), &
            Spn_njnh*Spn_nk12,1,(Spn_njnh-pil_s-pil_n), +1 )
      end do

      call rpn_comm_transpose ( Spn_wrk, ldnh_minx, ldnh_maxx, G_ni, Spn_njnh,&
                               Spn_12smin,Spn_12smax, G_nk, Spn_fft, -1, 2 )

      do k= 2, G_nk
         F_ft1(1:l_ni,1:l_nj,k)= F_ft1(1:l_ni,1:l_nj,k) + &
                prof(k)*Spn_wrk(1:l_ni,1:l_nj,k)*Spn_weight
      end do

!
!----------------------------------------------------------------------
!
      return
      end subroutine spn_apply

