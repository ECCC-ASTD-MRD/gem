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

module adz_tricub
  use ISO_C_BINDING
  use adz_mem
  implicit none
  public

contains
      subroutine adz_cubiclag ( F_stk,F_nptr,F_xyz,F_geom,F_num,&
                                F_i0,F_in,F_j0,F_jn,F_k0,F_mono_L )
      implicit none

      logical, intent(in) :: F_mono_L
      integer, intent(in) :: F_nptr, F_num,F_i0,F_in,F_j0,F_jn,F_k0
      real, dimension(*), intent(in ) :: F_xyz
      type(Adz_pntr_stack), dimension(F_nptr), intent(inout) :: F_stk
      type(C_PTR), intent(in) :: F_geom

      include "tricublin_f90.inc"

      type(C_PTR),dimension(F_nptr) :: stkpntr
      integer ij,i,j,j1,j2,k,k1,k2,kk,ni,nj,nij,n,n1,n2,np
      real, dimension(Adz_lminx:Adz_lmaxx,Adz_lminy:Adz_lmaxy,l_nk,F_nptr), target :: extended
      integer :: slice
      real, dimension(l_ni*l_nj*l_nk*F_nptr) :: wrkc,lin,mi,ma
!
!---------------------------------------------------------------------
!
      do n=1,F_nptr
         call rpn_comm_xch_halox( F_stk(n)%src, l_minx,l_maxx    ,&
              l_miny,l_maxy, l_ni,l_nj,l_nk, Adz_halox,Adz_haloy ,&
              .false., .false., extended(Adz_lminx,Adz_lminy,1,n),&
              Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy, l_ni, 0)
         stkpntr(n)= c_loc(extended(1,1,1,n))
      end do

      ni = F_in-F_i0+1
      nj = F_jn-F_j0+1
      nij= ni*nj
      n= nj
      slice = n*ni

      n1=1
      do while (n1 <= F_num)
         n2=min(n1+slice-1,F_num)
         np=n2-n1+1
         k1= n1/nij + F_k0
         k2= n2/nij + min(1,mod(n2,nij)) + F_k0 - 1
         ij=0

         if (F_mono_L) then
            call tricublin_mono_zyx_m_n (wrkc,lin,mi,ma, stkpntr, &
                             F_xyz((n1-1)*3+1),F_geom ,np, F_nptr )
            do n=1,F_nptr
            do k=k1,k2
               kk= (k-F_k0)*nj
               j1=max(1 ,n1/ni+1-kk) + F_j0 - 1
               j2=min(nj,n2/ni  -kk) + F_j0 - 1
               do j=j1,j2
               do i=F_i0,F_in
                  ij= ij+1
                  F_stk(n)%dst(i,j,k)= max(mi(ij), &
                                       min(ma(ij), wrkc(ij)))
               end do
               end do
            end do
            end do
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
!
!---------------------------------------------------------------------
!
      return
      end subroutine adz_cubiclag

end module adz_tricub
