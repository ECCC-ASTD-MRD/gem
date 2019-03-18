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

!**s/r ens_filter_gauss - Gaussian filter
!
      subroutine ens_filter_gauss(bfact, lambda, dsp_local)
      use dcst
      use geomh
      use glb_ld
      use lun
      use ptopo
      use tdpack
      use trp
      implicit none
#include <arch_specific.hf>
!
!author
!     Lubos Spacek - rpn - oct 2005
!
!revision
!
!object
!
!arguments
!     bfact - atenuation factor
!     lambda - wavelength
!     dsp_local - field to filter
!

      real*8, intent(in) :: bfact, lambda
      real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk), intent(inout) :: dsp_local

!
      integer kx,kz
      integer i, j, k
      real    err, min_local, min_global
      real*8  sigma, trx, try, pri8
!     Arrays
      integer,dimension(:),     allocatable  :: i1,j1
      real*8, dimension(:),     allocatable  :: freqy
      real*8, dimension(:,:),   allocatable  :: freqx
      real*8, dimension(:,:,:), allocatable  :: psi_local,work_local,ffx_local
      logical, save :: init_done
!
      allocate( freqx(l_nj,G_ni+2),freqy(G_nj+2))
      allocate( psi_local(l_minx:l_maxx,l_miny:l_maxy,l_nk))
      allocate( work_local(l_miny:l_maxy,Trp_12dmin:Trp_12dmax, &
                                                   G_ni+2+Ptopo_npex))
      allocate( ffx_local(Trp_12dmin:Trp_12dmax,Trp_22min:Trp_22max, &
                                                   G_nj+2+Ptopo_npey))
      work_local= 0.0d0
      ffx_local = 0.0d0

!
!     Determine sigma
!
!      bfact=1.0d-01 ; lambda=3.130793428843078e+05

      sigma=lambda*sqrt(-2.d0*log(bfact))/(2*pi_8);
!
!
!     Transfer to a real*8 array to apply 2D Fourier filter
!     ===================================================================
!
      psi_local=dble(dsp_local(:,:,1:l_nk))

      call rpn_comm_transpose( psi_local, l_minx, l_maxx, G_ni, &
                               (l_maxy-l_miny+1), &
                      Trp_12dmin, Trp_12dmax, G_nk, work_local, 1, 2 )
!
!     x-direction
!     ============
      call itf_fft_set(G_ni,'PERIODIC',pri8)
      do kz=1,Trp_12dn
        call ffft8(work_local(l_miny,kz,1),(l_maxy-l_miny+1) &
             *(Trp_12dmax-Trp_12dmin+1),1, (l_maxy-l_miny+1), -1 )
      end do
!
!     apply the Gaussian bell
!
      do j=1,l_nj
         trx=1.d0/(2*pi_8*Dcst_rayt_8*geomh_cy_8(j))**2
         do i=0,G_ni+1,2
           freqx(j,i+1:i+2)=-2.0d0*trx*(pi_8*sigma*dble(i/2))**2
         end do
      end do
      call vexp(freqx,freqx,l_nj*(G_ni+2))
!
         do i=1,G_ni+2
           do k=1,Trp_12dn
             do j=1,l_nj
               work_local(j,k,i)=work_local(j,k,i)*freqx(j,i)
             end do
           end do
         end do
!
      call rpn_comm_transpose( work_local, l_miny, l_maxy, G_nj, &
                               (Trp_12dmax-Trp_12dmin+1), &
                          Trp_22min , Trp_22max, G_ni, ffx_local, 2, 2 )
!
!     y-direction
!     ============
      call itf_fft_set(G_nj,'PERIODIC',pri8)
      do kx=1,Trp_22n
      call ffft8(ffx_local(Trp_12dmin,kx,1),(Trp_12dmax-Trp_12dmin+1) &
           *(Trp_22max-Trp_22min+1),1, (Trp_12dmax-Trp_12dmin+1), -1 )
      end do
!     Apply the Gaussian bell
      try=1.d0/(pi_8*Dcst_rayt_8)**2
        do j=0,G_nj+1,2
          freqy(j+1:j+2)=-2.d0*try*(pi_8*sigma*dble(j/2))**2
        end do
!
      call vexp(freqy,freqy,G_nj+2)
!
      do k=1,G_nj+2
        do j=1,Trp_22n
          do i=Trp_12dmin,Trp_12dmax
             ffx_local(i,j,k)=ffx_local(i,j,k)*freqy(k)
          end do
        end do
      end do
!
!     the way back to the original space
!     ==================================
      do kx=1,Trp_22n
      call ffft8(ffx_local(Trp_12dmin,kx,1),(Trp_12dmax-Trp_12dmin+1) &
           *(Trp_22max-Trp_22min+1),1, (Trp_12dmax-Trp_12dmin+1), +1 )
      end do
!
      call rpn_comm_transpose( work_local, l_miny, l_maxy, G_nj, &
                               (Trp_12dmax-Trp_12dmin+1), &
                          Trp_22min , Trp_22max, G_ni, ffx_local, -2, 2 )
      call setfft8(G_ni)
      do kz=1,Trp_12dn
        call ffft8(work_local(l_miny,kz,1),(l_maxy-l_miny+1) &
             *(Trp_12dmax-Trp_12dmin+1),1, (l_maxy-l_miny+1), +1 )
      end do
!
      call rpn_comm_transpose( psi_local, l_minx, l_maxx, G_ni, &
                               (l_maxy-l_miny+1), &
                      Trp_12dmin, Trp_12dmax, G_nk, work_local, -1, 2 )
!     Transfer back
!
      dsp_local(:,:,1:l_nk)=real(psi_local)

      min_local=minval(dsp_local)
      call rpn_comm_reduce (min_local,min_global, 1, &
                       "MPI_REAL","MPI_MIN",0,"grid",err )
      call rpn_comm_bcast(min_global, 1, "MPI_REAL", 0, "grid", err)

      if(min_global/=0.0)then
         where (dsp_local < -min_global)
            dsp_local=0.0
         end where
      end if
!
      INITIALIZE: if (.not.init_done) then
         if (Lun_out > 0) then
            write(*,*)'min_global',min_global
            write( Lun_out,1000)
            write( Lun_out,1009)lambda,bfact,sigma
            allocate(i1(5),j1(5))
            i1(1)=3 ;i1(2)=G_ni/4 ;i1(3)=G_ni/2 ;i1(4)=3*G_ni/4 ;i1(5)=G_ni
            j1(1)=3 ;j1(2)=G_nj/4 ;j1(3)=G_nj/2 ;j1(4)=3*G_nj/4 ;j1(5)=G_nj
            write(Lun_out,1010)(i1(j)/2,j=1,5)
            do i=1,l_nj
              write(Lun_out,1011)180.d0*geomh_y_8(i)/pi_8,(freqx(i,i1(j)),j=1,5)
            end do
            write(Lun_out,1012)(j1(j)/2,j=1,5),geomh_x_8(1), &
                                              (freqy(j1(j)),j=1,5)
            deallocate(i1,j1)
         end if
         init_done=.true.
      end if INITIALIZE
      deallocate(psi_local,work_local,ffx_local)
      deallocate(freqx,freqy)
!
 1000 format( &
         /,'INITIALIZE SCHEMES CONTROL PARAMETERS (S/R ENS_FILTER_GAUSS)', &
         /,'============================================================')
 1009 format( &
         /,'GAUSSIAN FILTER     ','     LAMBDA  ','         B    ', &
           '       SIGMA', &
         /,20X,3E14.6)
 1010 format( &
         /,'LATITUDE  FILTER VALUES FOR 1ST ROW WAVENUMBERS - DIRECTION G_NI', &
         /,64('='),/,7x,5i14,/)
 1011 format(f7.2,5e14.7)
 1012 format( &
         /,'LONGITUDE FILTER VALUES FOR WAVENUMBERS - DIRECTION G_NJ', &
         /,56('='),/,7x,5i14,/,f7.2,5e14.7)
      end subroutine ens_filter_gauss
