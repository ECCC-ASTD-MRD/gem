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

!**s/r ens_filter_ggauss - Gaussian grid space filter
!
      subroutine ens_filter_ggauss(bfact,lambda,dsp_local)
      use dcst
      use ens_options
      use gem_options
      use step_options
      use gmm_vt1
      use geomh
      use tdpack
      use glb_ld
      use lun
      use ldnh
      use glb_pil
      use trp
      use ptopo
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
  real(kind=REAL64)                            :: bfact, lambda
  real, dimension(l_minx:l_maxx,l_miny:l_maxy,l_nk) :: dsp_local
!
!author
!     Lubos Spacek - rpn - oct 2012
!
!revision
!
!object
!     Gaussian filter in gridpoint space
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! bfact         I         attenuation factor units [1]
! lambda        I         wavelength attenuated by the factor bfact
! dsp_local    I/O        filtered field
!


  integer i, j, k
  real(kind=REAL64)  sigma
!     Arrays
  real, dimension(7)                   :: fg
  real, dimension(:,:),   allocatable  :: f2, wk2
  real, dimension(:,:,:), allocatable  :: psi_local
  logical, save :: init_done=.false.
  real(kind=REAL64) half
      parameter( half  = 0.5 )

!
  allocate( psi_local(l_minx:l_maxx,l_miny:l_maxy,l_nk))
!
!     Determine sigma
!
  sigma=lambda*sqrt(-2.d0*log(bfact))/(2*pi_8);



!
  call rpn_comm_xch_halo (dsp_local,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk, &
               G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
  psi_local=dsp_local



! x-axis
  do k=1,l_nk
    do j=1,l_nj
!      fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hxu_8(1)*geomh_cy_8(j);fg(3)=-fg(5)
      fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hx_8*geomh_cy_8(j);fg(3)=-fg(5)
      fg(2)=2*fg(3);fg(1)=3*fg(3);fg(6)=2*fg(5);fg(7)=3*fg(5)
      fg=1/(sigma*sqrt(2*pi_8))*exp(-fg*fg/(2*sigma*sigma))
      fg=fg/sum(fg)
      do i=1,l_ni
        dsp_local(i,j,k)= sum(fg*psi_local(i-3:i+3,j,k))
      end do
    end do
  end do
  call rpn_comm_xch_halo (dsp_local,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,G_nk, &
               G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
  psi_local=dsp_local
!y-axis
!  fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hyv_8(1);fg(3)=-fg(5)
   fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hy_8;fg(3)=-fg(5)
  fg(2)=2*fg(3);fg(1)=3*fg(3);fg(6)=2*fg(5);fg(7)=3*fg(5)
  fg=1/(sigma*sqrt(2*pi_8))*exp(-fg*fg/(2*sigma*sigma))
  fg=fg/sum(fg)
  do k=1,l_nk
    do j=1,l_nj
      do i=1,l_ni
        dsp_local(i,j,k)= sum(fg*psi_local(i,j-3:j+3,k))
      end do
    end do
  end do
!
  INITIALIZE: if (.not.init_done) then
     allocate(f2(l_ni,l_nj),wk2(G_ni,G_nj))
     f2=1.0;wk2=2.0
     do j=1,l_nj
!        f2(4,j)=0.;f2(5,j)=Dcst_rayt_8*geomh_hxu_8(1)*geomh_cy_8(j)
        f2(4,j)=0.;f2(5,j)=Dcst_rayt_8*geomh_hx_8*geomh_cy_8(j)
        f2(3,j)=-f2(5,j)
        f2(2,j)=2*f2(3,j);f2(1,j)=3*f2(3,j)
        f2(6,j)=2*f2(5,j);f2(7,j)=3*f2(5,j)
        f2(1:7,j)=1/(sigma*sqrt(2*pi_8))* &
             exp(-f2(1:7,j)*f2(1:7,j)/(2*sigma*sigma))
        f2(1:7,j)=f2(1:7,j)/sum(f2(1:7,j))
        f2(8,j)=180.d0*geomh_y_8(j)/pi_8
     end do
     call glbcolc(wk2,G_ni,G_nj,f2,1,l_ni,1,l_nj,1)
!     fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hyv_8(1);fg(3)=-fg(5)
     fg(4)=0.;fg(5)=Dcst_rayt_8*geomh_hy_8;fg(3)=-fg(5)
     fg(2)=2*fg(3);fg(1)=3*fg(3);fg(6)=2*fg(5);fg(7)=3*fg(5)
     fg=1/(sigma*sqrt(2*pi_8))*exp(-fg*fg/(2*sigma*sigma))
     fg=fg/sum(fg)
     if (Lun_out > 0) then
        write(Lun_out,1000)
        write(Lun_out,1009)lambda,bfact,sigma
        do j=1,G_nj
           write(Lun_out,1011)wk2(8,j),(wk2(i,j),i=1,7)
        end do
        write(Lun_out,1011)0.d0,(fg(i),i=1,7)
        write(Lun_out,"(56('='),/)")
     end if
     deallocate(f2,wk2)
     init_done=.true.
  end if INITIALIZE

 deallocate(psi_local)
!
 1000 format( &
         /,'INITIALIZE SCHEMES CONTROL PARAMETERS (S/R ENS_FILTER_GAUSS)', &
         /,'============================================================')
 1009 format( &
         /,'GAUSSIAN FILTER     ','     LAMBDA  ','         B    ', &
           '       SIGMA', &
         /,64('='),/,20X,3E14.6)
 1011 format(f7.2,7e14.7)
      end subroutine ens_filter_ggauss
