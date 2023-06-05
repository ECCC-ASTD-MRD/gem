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

module adz_interp_hlt_mod
  use ISO_C_BINDING
  use adz_mem
  use adz_options
  use HORgrid_options
  use tr3d
  use mem_tstp
  use ens_gmm_var, only: mcrhsint
  implicit none
  public

contains
      subroutine adz_tricub_hlt ( F_stk,F_nptr,F_xyz,F_geom,F_num,&
                   F_i0,F_in,F_j0,F_jn,F_k0,F_ext_L,F_post,F_QV_L)
      implicit none

      logical, intent(IN), optional :: F_QV_L,F_ext_L
      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0
      real, dimension(*), intent(in ) :: F_xyz
      type(meta_tracers), dimension(F_nptr), intent(in), optional :: F_post
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(C_PTR), intent(in) :: F_geom

      include "tricublin_f90.inc"

      type(C_PTR),dimension(F_nptr) :: stkpntr
      logical linmima_l, sto_L, qv_L, ext_L
      integer ij,i,j,j1,j2,k,k1,k2,kk,ni,nj,nij,n,n1,n2,np,slc,dimext
      real :: pert
      integer :: slice
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
      real, dimension(:,:,:,:), pointer :: extended
!
!---------------------------------------------------------------------
!
      ni = F_in-F_i0+1
      nj = F_jn-F_j0+1
      nij= ni*nj
      if (nij<=0) return
      dimext= Adz_nij*l_nk
      
      ext_L = .true.
      if (present(F_ext_L)) ext_L= F_ext_L
      qv_L = .false.
      if (present(F_qv_L)) qv_L= F_qv_L

      if (ext_L) then !source is already extended
         if (qv_L) then
            print*, 'Source must be NOT expended when using quintic - abort'
            stop (1)
         endif
         do n=1,F_nptr
            stkpntr(n)= c_loc(F_stk(n)%src(1,1,1))
         end do
      else
         ! Only necessary for tracers since the source is NOT extended
         extended(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,1:l_nk,1:F_nptr) => WS1(1:)
         do n=1,F_nptr
!$omp do collapse(2)
            do k=1, l_nk
            do j=1, l_nj
            do i= 1, l_ni
               extended(i,j,k,n) = F_stk(n)%src(i,j,k)
            end do
            end do
            end do
!$omp end do
            stkpntr(n)= c_loc(extended(1,1,1,n))
         end do

!$omp single
         call rpn_comm_xch_halo (extended, Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,&
               l_ni,l_nj, F_nptr*l_nk, Adz_halox,Adz_haloy, .false.,.false., l_ni,0)
!$omp end single

         if (present(F_post)) then
            if (.not.Grd_yinyang_L.and.Adz_BC_LAM_flux==1) then
                call adz_BC_LAM_Aranami_hlt (extended,Adz_pb,Adz_num_b,1,Adz_lminx,Adz_lmaxx, &
                                             Adz_lminy,Adz_lmaxy,F_post,F_nptr)
            endif
         end if
      endif
      
      sto_L= associated(mcrhsint)
      linmima_l = present(F_post) .or. sto_L

      slice = ni*nj

!$omp do
      do slc= 1, F_num/slice + min(1,mod(F_num,slice))
         n1 = (slc-1)*slice + 1
         n2= min(n1+slice-1,F_num)
         np= n2-n1+1
         k1= n1/nij + F_k0
         k2= n2/nij + min(1,mod(n2,nij)) + F_k0 - 1
         ij= 0

         if (qv_L) call bicubHQV_lin_min_max_ntr ( wrkc,lin,mi,ma,extended ,&
                        F_xyz((n1-1)*3+1),np, Adz_lminx,Adz_lmaxx,Adz_lminy,&
                        Adz_lmaxy,dimext,F_nptr )

         if (linmima_l) then
            if (.not.qv_L) call tricublin_mono_zyx_m_n ( wrkc,lin,mi,ma,stkpntr,&
                                           F_xyz((n1-1)*3+1),F_geom ,np, F_nptr )
            do n=1,F_nptr
               do k=k1,k2
                  kk= (k-F_k0)*nj
                  j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
                  j2=min(nj,n2/ni  -kk) + F_j0 - 1
                  if (sto_L) then
                     do j=j1,j2
!DIR$ SIMD
                        do i=F_i0,F_in
                           ij= ij+1
                           pert = mcrhsint(i,j) * abs(wrkc(ij)-lin(ij))
                           F_stk(n)%dst(i,j,k)= wrkc(ij) + pert
                        end do
                     end do
                  else
                     do j=j1,j2
!DIR$ SIMD
                        do i=F_i0,F_in
                           ij= ij+1
                           F_stk(n)%dst(i,j,k)= wrkc(ij)
                        end do
                     end do
                  endif 
               end do
            end do
            ij= 0
            if (present(F_post)) then
               do n=1,F_nptr
                  do k=k1,k2
                     kk= (k-F_k0)*nj
                     j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
                     j2=min(nj,n2/ni  -kk) + F_j0 - 1
                     do j=j1,j2
!DIR$ SIMD
                        do i=F_i0,F_in
                           ij= ij+1
                           Adz_post(n)%lin(i,j,k) = max(mi(ij),min(ma(ij),lin(ij)))
                           Adz_post(n)%min(i,j,k) = mi  (ij)
                           Adz_post(n)%max(i,j,k) = ma  (ij)
                        end do
                     end do
                  end do
               enddo
            endif
         else
            if (.not.qv_L) call tricublin_zyx1_m_n ( wrkc, stkpntr, F_xyz((n1-1)*3+1),&
                                                     F_geom ,np, F_nptr )
            do n=1,F_nptr
            do k=k1,k2
               kk= (k-F_k0)*nj
               j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
               j2=min(nj,n2/ni  -kk) + F_j0 - 1
               do j=j1,j2
!DIR$ SIMD
               do i=F_i0,F_in
                  ij= ij+1
                  F_stk(n)%dst(i,j,k)= wrkc(ij)
               end do
               end do
            end do
            end do
         endif
      end do
!$omp enddo nowait
!---------------------------------------------------------------------
!
      return
      end subroutine adz_tricub_hlt
      
end module adz_interp_hlt_mod
