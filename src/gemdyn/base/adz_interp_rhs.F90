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

module adz_interp_rhs_mod
  use ISO_C_BINDING
  use adz_mem
  use adz_options
  use gem_timing
  use HORgrid_options
  use tr3d
  use ptopo
  use ens_gmm_var, only: mcrhsint
  implicit none
  public

contains
      subroutine adz_tricub_rhs ( F_stk,F_nptr,F_xyz,F_geom,F_num,&
                                  F_i0,F_in,F_j0,F_jn,F_k0,F_post )
      implicit none

      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0
      real, dimension(*), intent(in ) :: F_xyz
      type(meta_tracers), dimension(F_nptr), intent(in), optional :: F_post
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(C_PTR), intent(in) :: F_geom

      include "tricublin_f90.inc"

      type(C_PTR),dimension(F_nptr) :: stkpntr
      logical linmima_l, sto_L
      integer ij,i,j,j1,j2,k,k1,k2,kk,ni,nj,nij,n,n1,n2,np
      real :: pert
      real, dimension(Adz_lminx:Adz_lmaxx, Adz_lminy:Adz_lmaxy,l_nk,F_nptr), target :: extended
      integer :: slice
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
!
!---------------------------------------------------------------------
!
      Adz_cnt_int= Adz_cnt_int+1
! This next loop should use a normal rpn_comm_xch_halo
! once SLOD is completed
      do n=1,F_nptr
         call rpn_comm_xch_halox( F_stk(n)%src, l_minx,l_maxx    ,&
              l_miny,l_maxy, l_ni,l_nj,l_nk, Adz_halox,Adz_haloy ,&
              .false., .false., extended(Adz_lminx,Adz_lminy,1,n),&
              Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, 0)
         stkpntr(n)= c_loc(extended(1,1,1,n))
      end do
      call time_trace_barr(gem_time_trace, 10500+Adz_cnt_int,&
                           Gem_trace_barr, Ptopo_intracomm, MPI_BARRIER)

      if (present(F_post)) then
         call gemtime_start (83, 'C_BCFLUX_TR', 37)
         if (.not.Grd_yinyang_L.and.Adz_BC_LAM_flux==1) then
            call adz_BC_LAM_Aranami (extended,Adz_pb(1,Adz_i0b,Adz_j0b,Adz_k0),Adz_lminx,Adz_lmaxx, &
                                  Adz_lminy,Adz_lmaxy,F_post,F_nptr)
         endif
         call gemtime_stop  (83)
      end if
      
      sto_L= associated(mcrhsint)
      linmima_l = present(F_post) .or. sto_L

      if (linmima_l) call gemtime_start (84, 'TRICUBLIN', 37)

      ni = F_in-F_i0+1
      nj = F_jn-F_j0+1
      nij= ni*nj
      n= nj
      slice = n*ni

      n1=1
      do while (n1 <= F_num)
         n2= min(n1+slice-1,F_num)
         np= n2-n1+1
         k1= n1/nij + F_k0
         k2= n2/nij + min(1,mod(n2,nij)) + F_k0 - 1
         ij= 0

         if (linmima_l) then
            call tricublin_mono_zyx_m_n (wrkc,lin,mi,ma, stkpntr, &
                             F_xyz((n1-1)*3+1),F_geom ,np, F_nptr )
            do n=1,F_nptr
               do k=k1,k2
                  kk= (k-F_k0)*nj
                  j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
                  j2=min(nj,n2/ni  -kk) + F_j0 - 1
                  if (sto_L) then
                     do j=j1,j2
                        do i=F_i0,F_in
                           ij= ij+1
                           pert = mcrhsint(i,j) * abs(wrkc(ij)-lin(ij))
                           F_stk(n)%dst(i,j,k)= wrkc(ij) + pert
                        end do
                     end do
                  else
                     do j=j1,j2
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
            call tricublin_zyx1_m_n ( wrkc, stkpntr, F_xyz((n1-1)*3+1),&
                                      F_geom ,np, F_nptr )
            do n=1,F_nptr
            do k=k1,k2
               kk= (k-F_k0)*nj
               j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
               j2=min(nj,n2/ni  -kk) + F_j0 - 1
               do j=j1,j2
               do i=F_i0,F_in
                  ij= ij+1
                  F_stk(n)%dst(i,j,k)= wrkc(ij)
               end do
               end do
            end do
            end do
         endif
         n1=n2+1
      end do
      if (linmima_l) call gemtime_stop (84)
      call time_trace_barr(gem_time_trace, 10600+adz_cnt_int,&
                           Gem_trace_barr, Ptopo_intracomm, MPI_BARRIER)
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_tricub_rhs

      subroutine adz_bicubHQV_rhs ( F_stk,F_nptr,F_pxt,F_pyt,F_pzt,F_num,&
                                    F_i0,F_in,F_j0,F_jn,F_k0,F_post )
      implicit none

      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0
      real, dimension(*), intent(in) :: F_pxt,F_pyt,F_pzt
      type(meta_tracers), dimension(F_nptr), intent(in), optional :: F_post
      type(Adz_pntr_stack), dimension(F_nptr),  intent(inout) :: F_stk

      logical linmima_l, sto_L
      integer ijk,i,j,k,n
      real :: pert
      real, dimension(Adz_lminx:Adz_lmaxx, Adz_lminy:Adz_lmaxy,l_nk,F_nptr), target :: extended
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
      integer, pointer, contiguous, dimension(:) :: ii_w
!
!---------------------------------------------------------------------
!
      do n=1,F_nptr
         call rpn_comm_xch_halox( F_stk(n)%src, l_minx,l_maxx    ,&
              l_miny,l_maxy, l_ni,l_nj,l_nk, Adz_halox,Adz_haloy ,&
              .false., .false., extended(Adz_lminx,Adz_lminy,1,n),&
              Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, 0)
      end do
      
      if (present(F_post)) then
         if (.not.Grd_yinyang_L.and.Adz_BC_LAM_flux==1) then
            call adz_BC_LAM_Aranami (extended,Adz_pb(1,Adz_i0b,Adz_j0b,Adz_k0),Adz_lminx,Adz_lmaxx, &
                                  Adz_lminy,Adz_lmaxy,F_post,F_nptr)
         endif
      end if
      
      sto_L= associated(mcrhsint)
      linmima_l = present(F_post) .or. sto_L
      
      allocate (ii_w(4*F_num))

      call adv_get_indices (ii_w,F_pxt,F_pyt,F_pzt,l_ni*l_nj*l_nk,F_num, &
                            F_i0,F_in,F_j0,F_jn,F_k0,l_nk,'t')

      if (linmima_l) then

         do n=1,F_nptr

            call bicubHQV_lin_min_max (wrkc,lin,mi,ma,extended(Adz_lminx,Adz_lminy,1,n), &
                                       F_pxt,F_pyt,F_pzt,l_ni*l_nj*l_nk,F_num,ii_w,l_nk,'t')
            if (sto_L) then
               do k=F_k0,l_nk
                  do j=F_j0,F_jn
                     do i=F_i0,F_in
                        ijk= (k-1)*l_ni*l_nj + (j-1)*l_ni + i
                        pert = mcrhsint(i,j) * abs(wrkc(ijk)-lin(ijk))
                        F_stk(n)%dst(i,j,k)= wrkc(ijk) + pert
                     end do
                  end do
               end do
            else
               do k=F_k0,l_nk
                  do j=F_j0,F_jn
                     do i=F_i0,F_in
                        ijk= (k-1)*l_ni*l_nj + (j-1)*l_ni + i
                        F_stk(n)%dst(i,j,k)= wrkc(ijk)
                     end do
                  end do
               end do
            endif               
            if (present(F_post)) then
               do k=F_k0,l_nk
                  do j=F_j0,F_jn
                     do i=F_i0,F_in
                        ijk= (k-1)*l_ni*l_nj + (j-1)*l_ni + i
                        Adz_post(n)%lin(i,j,k) = max(mi(ijk),min(ma(ijk),lin(ijk)))
                        Adz_post(n)%min(i,j,k) = mi  (ijk)
                        Adz_post(n)%max(i,j,k) = ma  (ijk)
                     end do
                  end do
               end do
            endif
         end do

      else

         do n=1,F_nptr

            call bicubHQV_lin_min_max (wrkc,lin,mi,ma,extended(Adz_lminx,Adz_lminy,1,n), &
                                       F_pxt,F_pyt,F_pzt,l_ni*l_nj*l_nk,F_num,ii_w,l_nk,'t')

            do k=F_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     ijk= (k-1)*l_ni*l_nj + (j-1)*l_ni + i
                     F_stk(n)%dst(i,j,k)= wrkc(ijk)
                  end do
               end do
            end do

         end do

      end if

      deallocate (ii_w)
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_bicubHQV_rhs

end module adz_interp_rhs_mod
