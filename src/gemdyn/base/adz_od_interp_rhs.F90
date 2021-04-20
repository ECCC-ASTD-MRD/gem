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

module adz_od_interp_rhs
  use ISO_C_BINDING
  use gem_options
  use adz_mem
  use adz_options
  use ens_gmm_var, only: mcrhsint
  use HORgrid_options
  implicit none
  public
      
contains
      subroutine adz_od_tricub ( F_stk,F_nptr,F_xyz,F_geom,F_num ,&
                        F_export,F_i0,F_in,F_j0,F_jn,F_k0,F_spp_L )
      implicit none

      logical, intent(IN) :: F_spp_L
      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0
      real, dimension(*), intent(in ) :: F_xyz
      type(ADZ_SLOD), intent(INOUT) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(C_PTR), intent(in) :: F_geom

      include "tricublin_f90.inc"

      type(C_PTR),dimension(F_nptr) :: stkpntr
      integer ij,i,j,j1,j2,k,k1,k2,kk,ni,nj,nij,n,n1,n2,np
      integer :: slice
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc
!
!---------------------------------------------------------------------
!
      if (F_spp_L) then
         call adz_od_tricub_mono ( F_stk,F_nptr,F_xyz,F_geom,F_num,&
                         F_export,F_i0,F_in,F_j0,F_jn,F_k0,F_spp_L )
         return
      endif
      
      do n=1,F_nptr
         call rpn_comm_xch_halo( F_stk(n)%src,l_minx,l_maxx,l_miny,l_maxy,&
                                 l_ni,l_nj,l_nk,G_halox,G_haloy,&
         G_periodx,G_periody,G_ni,0)
         stkpntr(n)= c_loc(F_stk(n)%src(1,1,1))
      end do
      
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
         n1=n2+1
      end do

      call adz_od_corrhs (F_stk,stkpntr,F_nptr,F_geom,F_export)
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_od_tricub

      subroutine  adz_od_tricub_mono ( F_stk,F_nptr,F_xyz,F_geom,F_num, &
                       F_export,F_i0,F_in,F_j0,F_jn,F_k0,F_spp_L,F_post )
      use tr3d
      implicit none

      logical, intent(IN) :: F_spp_L
      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0
      type(post), dimension(F_nptr), intent(in), optional :: F_post
      real, dimension(*), intent(in ) :: F_xyz
      type(ADZ_SLOD), intent(INOUT) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(C_PTR), intent(in) :: F_geom

      include "tricublin_f90.inc"

      type(C_PTR),dimension(F_nptr) :: stkpntr
      integer ij,i,j,j1,j2,k,k1,k2,kk,ni,nj,nij,n,n1,n2,np
      integer :: slice
      real :: pert
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
!
!---------------------------------------------------------------------
!
      do n=1,F_nptr
         call rpn_comm_xch_halo( F_stk(n)%src,l_minx,l_maxx,l_miny,l_maxy,&
                                 l_ni,l_nj,l_nk,G_halox,G_haloy,&
                                 G_periodx,G_periody,G_ni,0)
         stkpntr(n)= c_loc(F_stk(n)%src(1,1,1))
      end do

      if (present(F_post).and..not.Grd_yinyang_L.and.Adz_BC_LAM_flux==1) then

         call adz_od_BC_LAM_Aranami (F_stk,Adz_pb,Adz_expt,Adz_num_b,1, &
                                     l_minx,l_maxx,l_miny,l_maxy,Tr_3CWP,F_nptr)

      end if
      
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

            call tricublin_mono_zyx_m_n (wrkc,lin,mi,ma, stkpntr, &
                             F_xyz((n1-1)*3+1),F_geom ,np, F_nptr )

            ij= 0
            do n=1,F_nptr
               do k=k1,k2
                  kk= (k-F_k0)*nj
                  j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
                  j2=min(nj,n2/ni  -kk) + F_j0 - 1
                  if (F_spp_L) then
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
                           F_post(n)%lin(i,j,k) = max(mi(ij),min(ma(ij),lin(ij)))
                           F_post(n)%min(i,j,k) = mi  (ij)
                           F_post(n)%max(i,j,k) = ma  (ij)
                        end do
                     end do
                  end do
               enddo
            endif
            
         n1=n2+1
      end do
         
      if (present(F_post)) then
         call adz_od_corrhs_mono (F_stk,stkpntr,F_nptr,F_geom,F_export,F_spp_L,F_post=F_post)
      else
         call adz_od_corrhs_mono (F_stk,stkpntr,F_nptr,F_geom,F_export,F_spp_L)
      endif
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_od_tricub_mono

      subroutine adz_od_corrhs_mono (F_stk,F_cptr,F_nptr,F_geom,&
                                     F_export,F_spp_L,F_post)
      use ptopo
      implicit none

      logical, intent(IN) :: F_spp_L
      integer, intent(in) :: F_nptr
      type(C_PTR), intent(in) :: F_geom
      type(C_PTR), dimension (F_nptr), intent(in) :: F_cptr
      type(ADZ_SLOD), intent(INOUT) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(post), dimension(F_nptr), intent(in), optional :: F_post

      include "tricublin_f90.inc"
      integer :: i,j,k,n, cnt,dim,nr,jj,km,err,ijk, req, err_dest
      integer, dimension(Ptopo_world_numproc) :: rank, nreq, dsp
      real, dimension(Adz_MAX_MPI_OS_SIZE*F_nptr) :: send,lin,mi,ma
      real, dimension(Adz_MAX_MPI_OS_SIZE*3     ) :: request_from
      real :: pert, cubic, linea
      real, dimension(:,:,:,:), allocatable :: cor
!
!     ---------------------------------------------------------------
!
! get requested positions from others and post return adresses
      dim= Adz_MAX_MPI_OS_SIZE*3
      call adz_traj_from (request_from,nr,F_export,F_nptr*4,dim)

! compute the requested interpolated values and put the results in
! Adz_cor which is exposed to all PEs through Adz_wincor window

      call tricublin_mono_zyx_m_n (send,lin,mi,ma,  F_cptr, &
                            request_from, F_geom ,nr, F_nptr)
      err_dest=0
      if (4*F_nptr*nr>size(Adz_cor)) then
         err_dest= -1
      else
         call xfer_2_win (Adz_cor, send,lin,mi,ma, nr,F_nptr)
      endif
      call gem_error(err_dest,'adz_od_corrhs_mono',&
           'NOT ENOUGH MEMORY - Adz_MAX_MPI_OS_SIZE for Adz_cor')

!     get requested positions from others

      ORIGIN_DATATYPE = REAL_DATATYPE
      TARGET_DATATYPE = REAL_DATATYPE
      dim= l_ni*l_nj

      call MPI_Win_fence(0, Adz_offs_win, err)
      call rpn_comm_barrier ("MULTIGRID", err)
      call MPI_Win_fence(0, Adz_wincor, err)

      cnt= 0
      do req=1,2*Ptopo_world_numproc,2
         if (Adz_offs(req) > 0) then
            cnt=cnt+1
            rank(cnt)= req/2
            nreq(cnt)= Adz_offs(req)
            dsp (cnt)= Adz_offs(req+1)
         endif
         nr= cnt
      end do
      
! fetch and plug interpolated requests from neighbors
      allocate (cor(4,F_nptr,Adz_MAX_MPI_OS_SIZE,nr))
      do req=1,nr
         ORIGIN_COUNT= nreq(req) * 4*F_nptr
         TARGET_RANK = rank(req)
         TARGET_DISP = dsp (req)
         TARGET_COUNT= ORIGIN_COUNT
         call MPI_Get( cor(1,1,1,req), ORIGIN_COUNT, ORIGIN_DATATYPE,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Adz_wincor, err )
      end do
      call MPI_Win_fence(0, Adz_wincor, err)
      do req=1,nr
         do km=1, nreq(req)
            ijk= F_export%dest(km,rank(req))
            k = ijk/dim - 1 + min(mod(ijk,dim),1)
            jj= ijk-k*dim
            j = jj/l_ni - 1 + min(mod(jj,l_ni),1)
            i = jj - j*l_ni ; j=j+1 ; k=k+1
            do n=1, F_nptr
               cubic= cor(1,n,km,req)
               linea= cor(2,n,km,req)
               if (F_spp_L) then
                  pert = mcrhsint(i,j) * abs(cubic-linea)
                  F_stk(n)%dst(i,j,k) = cubic + pert
               else
                  F_stk(n)%dst(i,j,k) = cubic
               endif
               if (present(F_post)) then
                  F_post(n)%lin(i,j,k) = max(cor(3,n,km,req),min(cor(4,n,km,req),linea))
                  F_post(n)%min(i,j,k) = cor(3,n,km,req)
                  F_post(n)%max(i,j,k) = cor(4,n,km,req)
               endif
            end do
         end do
      enddo
      deallocate (cor)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_corrhs_mono
      
      subroutine xfer_2_win (F_cor,send,lin,mi,ma,nr, F_nptr)
      implicit none
      integer nr, F_nptr
      real, dimension(*) :: send,lin,mi,ma
      real, dimension(4,F_nptr,nr) :: F_cor

      integer i,j
      do i=1,nr
         do j=1,F_nptr
            F_cor(1,j,i)= send(i+(j-1)*nr)
            F_cor(2,j,i)= lin (i+(j-1)*nr)
            F_cor(3,j,i)= mi  (i+(j-1)*nr)
            F_cor(4,j,i)= ma  (i+(j-1)*nr)
         end do
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine xfer_2_win

      subroutine adz_od_bicubHQV_rhs ( F_stk,F_nptr,F_xyz,F_num,F_export, &
                                       F_i0,F_in,F_j0,F_jn,F_k0,F_spp_L,F_post )

      use tr3d
      implicit none

      logical, intent(in) :: F_spp_L
      integer, intent(in) :: F_nptr,F_num,F_i0,F_in,F_j0,F_jn,F_k0
      type(post), dimension(F_nptr), intent(in), optional :: F_post
      real, dimension(*), intent(in ) :: F_xyz
      type(ADZ_SLOD), intent(inout) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk

      integer :: ijk,i,j,k,n
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
      real :: pert
!
!---------------------------------------------------------------------
!
      do n=1,F_nptr
         call rpn_comm_xch_halo( F_stk(n)%src,l_minx,l_maxx,l_miny,l_maxy,&
                                 l_ni,l_nj,l_nk,G_halox,G_haloy,&
                                 G_periodx,G_periody,G_ni,0)
      end do

      if (present(F_post).and..not.Grd_yinyang_L.and.Adz_BC_LAM_flux==1) then

         call adz_od_BC_LAM_Aranami (F_stk,Adz_pb,Adz_expt,Adz_num_b,1, &
                                     l_minx,l_maxx,l_miny,l_maxy,Tr_BQWP,F_nptr)

      end if

      do n=1,F_nptr

         call bicubHQV_lin_min_max (wrkc,lin,mi,ma,F_stk(n)%src,F_xyz,F_num, &
                                    l_minx,l_maxx,l_miny,l_maxy)

         if (F_spp_L) then
            ijk = 0
            do k=F_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     ijk = ijk + 1 
                     pert = mcrhsint(i,j) * abs(wrkc(ijk)-lin(ijk))
                     F_stk(n)%dst(i,j,k)= wrkc(ijk) + pert
                  end do
               end do
            end do
         else
            ijk = 0
            do k=F_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     ijk = ijk + 1 
                     F_stk(n)%dst(i,j,k)= wrkc(ijk)
                  end do
               end do
            end do
         end if

         if (present(F_post)) then
            ijk = 0
            do k=F_k0,l_nk
               do j=F_j0,F_jn
                  do i=F_i0,F_in
                     ijk = ijk + 1 
                     F_post(n)%lin(i,j,k) = max(mi(ijk),min(ma(ijk),lin(ijk)))
                     F_post(n)%min(i,j,k) = mi  (ijk)
                     F_post(n)%max(i,j,k) = ma  (ijk)
                  end do
               end do
            end do
         end if

      end do

      if (present(F_post)) then
         call adz_od_corrhs_mono_BQ (F_stk,F_nptr,F_export,F_spp_L,F_post=F_post)
      else
         call adz_od_corrhs_mono_BQ (F_stk,F_nptr,F_export,F_spp_L)
      end if
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_od_bicubHQV_rhs

      subroutine adz_od_corrhs_mono_BQ (F_stk,F_nptr,F_export,F_spp_L,F_post)

      use ptopo
      implicit none

      logical, intent(in) :: F_spp_L
      integer, intent(in) :: F_nptr
      type(ADZ_SLOD), intent(inout) :: F_export
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(post), dimension(F_nptr), intent(in), optional :: F_post

      integer :: i,j,k,n, cnt,dim,nr,jj,km,err,ijk, req, err_dest, pos
      integer, dimension(Ptopo_world_numproc) :: rank, nreq, dsp
      real, dimension(Adz_MAX_MPI_OS_SIZE*F_nptr) :: send,lin,mi,ma
      real, dimension(Adz_MAX_MPI_OS_SIZE*3     ) :: request_from
      real :: pert, cubic, linea
      real, dimension(:,:,:,:), allocatable :: cor
!
!     ---------------------------------------------------------------
!
! get requested positions from others and post return adresses
      dim= Adz_MAX_MPI_OS_SIZE*3
      call adz_traj_from (request_from,nr,F_export,F_nptr*4,dim)

! compute the requested interpolated values and put the results in
! Adz_cor which is exposed to all PEs through Adz_wincor window

      do n=1,F_nptr

         pos = 1+(n-1)*nr

         call bicubHQV_lin_min_max (send(pos),lin(pos),mi(pos),ma(pos), &
                                    F_stk(n)%src,request_from,nr,       &
                                    l_minx,l_maxx,l_miny,l_maxy)

      end do

      err_dest=0
      if (4*F_nptr*nr>size(Adz_cor)) then
         err_dest= -1
      else
         call xfer_2_win (Adz_cor, send,lin,mi,ma, nr,F_nptr)
      endif
      call gem_error(err_dest,'adz_od_corrhs_mono_BQ',&
           'NOT ENOUGH MEMORY - Adz_MAX_MPI_OS_SIZE for Adz_cor')

!     get requested positions from others

      ORIGIN_DATATYPE = REAL_DATATYPE
      TARGET_DATATYPE = REAL_DATATYPE
      dim= l_ni*l_nj

      call MPI_Win_fence(0, Adz_offs_win, err)
      call rpn_comm_barrier ("MULTIGRID", err)

      cnt= 0
      do req=1,2*Ptopo_world_numproc,2
         if (Adz_offs(req) > 0) then
            cnt=cnt+1
            rank(cnt)= req/2
            nreq(cnt)= Adz_offs(req)
            dsp (cnt)= Adz_offs(req+1)
         endif
         nr= cnt
      end do
      
! fetch and plug interpolated requests from neighbors
      allocate (cor(4,F_nptr,Adz_MAX_MPI_OS_SIZE,nr))
      do req=1,nr
         ORIGIN_COUNT= nreq(req) * 4*F_nptr
         TARGET_RANK = rank(req)
         TARGET_DISP = dsp (req)
         TARGET_COUNT= ORIGIN_COUNT
         call MPI_Get( cor(1,1,1,req), ORIGIN_COUNT, ORIGIN_DATATYPE,&
                       TARGET_RANK,TARGET_DISP, TARGET_COUNT,&
                       TARGET_DATATYPE, Adz_wincor, err )
      end do
      call MPI_Win_fence(0, Adz_wincor, err)
      do req=1,nr
         do km=1, nreq(req)
            ijk= F_export%dest(km,rank(req))
            k = ijk/dim - 1 + min(mod(ijk,dim),1)
            jj= ijk-k*dim
            j = jj/l_ni - 1 + min(mod(jj,l_ni),1)
            i = jj - j*l_ni ; j=j+1 ; k=k+1
            do n=1, F_nptr
               cubic= cor(1,n,km,req)
               linea= cor(2,n,km,req)
               if (F_spp_L) then
                  pert = mcrhsint(i,j) * abs(cubic-linea)
                  F_stk(n)%dst(i,j,k) = cubic + pert
               else
                  F_stk(n)%dst(i,j,k) = cubic
               endif
               if (present(F_post)) then
                  F_post(n)%lin(i,j,k) = max(cor(3,n,km,req),min(cor(4,n,km,req),linea))
                  F_post(n)%min(i,j,k) = cor(3,n,km,req)
                  F_post(n)%max(i,j,k) = cor(4,n,km,req)
               endif
            end do
         end do
      enddo
      deallocate (cor)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_corrhs_mono_BQ

end module adz_od_interp_rhs
