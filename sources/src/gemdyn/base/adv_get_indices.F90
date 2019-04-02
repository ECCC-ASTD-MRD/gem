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

!@objective:  pre-compute indices to be used in tricub interp  adx_tricub_lag3d_loop.cdk

   subroutine adv_get_indices (ii,  F_x, F_y, F_z, F_num, &
                                 nind, F_i0, F_in, F_j0, F_jn, F_k0, F_nk, F_lev_S)
   use gem_options
   use glb_ld
   use ver
   use adv_grid
   use adv_interp
   use outgrid
   implicit none

#include <arch_specific.hf>

!@arguments
   character(len=*), intent(in) :: F_lev_S                          ! m/t : Momemtum/thermo level
   integer, intent(in) :: F_nk                                      ! number of vertical levels
   integer, intent(in) :: F_i0, F_in, F_j0, F_jn, F_k0              ! scope of operator
   integer, intent(in) :: F_num
   real,dimension(F_num), intent(in) :: F_x, F_y, F_z ! interpolation target x,y,z coordinates

   !
   ! author  Rabah Aider    August 2015
   !

   integer :: kkmax , idxk, idxjk, ii1, jj1,kk1
   integer :: i,j,k,n,n0,n1,n2,n3,m,nind
   integer :: midxk,midxjk,mni,mnj,mnk,nn
   real*8  :: p_z00_8
   real*8  :: rri,rrj,rrk
   integer, dimension(4*nind), intent(out) :: ii    ! index, to be used in tricubic lagrangian interpolations
   integer, dimension(:),pointer, contiguous :: p_lcz
   real*8,  dimension(:),pointer, contiguous :: p_bsz_8

!---------------------------------------------------------------------

! Vertical variable type:  Height--> sig <0 , Pressure --> sig >0
    sig=int((Ver_z_8%m(l_nk)-Ver_z_8%m(1))/(abs(  Ver_z_8%m(l_nk)-Ver_z_8%m(1) )))

    p_z00_8 = Ver_z_8%m(0)
    kkmax   = F_nk-1
    if (F_lev_S == 'm') then
      p_lcz     => adv_lcz%m
      p_bsz_8   => adv_bsz_8%m
    elseif  (F_lev_S == 't') then
      p_lcz     => adv_lcz%t
      p_bsz_8   => adv_bsz_8%t
    elseif  (F_lev_S == 'x') then
      p_lcz     => adv_lcz%x
      p_bsz_8   => adv_bsz_8%x
    end if

! pre-compute indices ii

    mni=F_in-F_i0+1
    mnj=F_jn-F_j0+1
    mnk=F_nk-F_k0+1

!$omp parallel private(i,j,k, ii1, jj1, kk1,&
!$omp                  n,n0,n1,n2,n3,nn,m,&
!$omp                  idxk,idxjk,midxk,midxjk,&
!$omp                  rri,rrj,rrk) &
!$omp          shared(ii,kkmax,sig,p_lcz,p_z00_8,p_bsz_8,&
!$omp                 mni,mnj,mnk)

!$omp do

   do k=F_k0,F_nk
         idxk  = (k-1)*l_ni*l_nj
         midxk = (k-F_k0)*mni*mnj
         do j=F_j0,F_jn
            idxjk = idxk + ((j-1)*l_ni)
            midxjk=midxk + ((j-F_j0)*mni)
            do i=F_i0,F_in
               n = idxjk + i
               nn= midxjk+ (i-F_i0+1)
               m=nn*4

               n0=m-3
               n1=m-2
               n2=m-1
               n3=m

               ii(n0)=n

               rri = F_x(n)
               ii1 = 1 + (rri - adv_x00_8) * adv_ovdx_8
               ii(n1) = max(2,min(ii1,adv_iimax))

               rrj = F_y(n)
               jj1 = 1 + (rrj - adv_y00_8) * adv_ovdy_8
               ii(n2) = max(G_haloy,min(jj1,adv_jjmax))

               rrk = F_z(n)
               kk1 = (rrk - p_z00_8) * adv_ovdz_8*sig
               kk1 = p_lcz (kk1+1)
               if ( real(sig) * (rrk - p_bsz_8(kk1)) < 0. ) kk1 = kk1 - 1
               ii(n3) = min(kkmax-1,max(0,kk1))
            end do
         end do
   end do

!$omp enddo
!$omp end parallel

end subroutine adv_get_indices
