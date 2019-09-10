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

!**s/r adz_BC_LAM_Aranami - Estimate FLUX_out/FLUX_in based on Aranami et al. (2015)

      subroutine adz_BC_LAM_Aranami (F_extended,F_aminx,F_amaxx,F_aminy,F_amaxy,F_post,F_nptr)

      use adz_mem
      use adz_options
      use HORgrid_options
      use tr3d

      implicit none

#include <arch_specific.hf>

      !arguments
      !---------
      integer, intent(IN) :: F_aminx,F_amaxx,F_aminy,F_amaxy,F_nptr
      integer, dimension(F_nptr), intent(IN) :: F_post
      real, dimension(F_aminx:F_amaxx,F_aminy:F_amaxy,l_nk,F_nptr), intent(IN) :: F_extended

      include "tricublin_f90.inc"

      !object
      !=============================================================
      !     Estimate FLUX_out/FLUX_in based on Aranami et al. (2015)
      !=============================================================

      integer :: n,itr,k,nf,nptr_F,tr_num(F_nptr),&
                 ij,i,j,j1,j2,k1,k2,kk,ni,nj,nij,n1,n2,np,slice,nn,i0,in,j0,jn
      real, dimension(F_aminx:F_amaxx,F_aminy:F_amaxy,l_nk,F_nptr), target :: adv_o,adv_i
      real, dimension(l_ni*l_nj*l_nk*F_nptr*2) :: wrkc
      type(C_PTR),dimension(F_nptr*2) :: stkpntr
!
!---------------------------------------------------------------------
!
      if (Grd_yinyang_L.or.Adz_BC_LAM_flux==2) return

      i0 = F_aminx + Adz_halox*west
      in = F_amaxx - Adz_halox*east
      j0 = F_aminy + Adz_haloy*south
      jn = F_amaxy - Adz_haloy*north

      !Bermejo-Conde LAM: Apply mask_o/mask_i to Tracer
      !------------------------------------------------
      if (F_nptr<=MAXTR3D) then !Treat Tracer in Tr3d_list

         nf = 0

         do n=1,F_nptr

            itr = F_post(n)

            if (Tr3d_mass(itr)/=1) cycle

            nf = nf + 1

            tr_num(nf) = n

            do k=1,l_nk
            do j=j0,jn
            do i=i0,in
               adv_o(i,j,k,nf) = Adz_BC_LAM_mask_o(i,j,k) * F_extended(i,j,k,n)
               adv_i(i,j,k,nf) = Adz_BC_LAM_mask_i(i,j,k) * F_extended(i,j,k,n)
            end do
            end do
            end do

            stkpntr(2*(nf-1) + 1) = c_loc(adv_o(1,1,1,nf))
            stkpntr(2*(nf-1) + 2) = c_loc(adv_i(1,1,1,nf))

         end do

      !PSADJ LAM: Use mask_o/mask_i based on Tracer=1
      !----------------------------------------------
      else !Treat Tracer=1 used in PSADJ_LAM

         nf = 1

         tr_num(nf) = 1

         adv_o(i0:in,j0:jn,1:l_nk,1) = Adz_BC_LAM_mask_o(i0:in,j0:jn,1:l_nk)

         adv_i(i0:in,j0:jn,1:l_nk,1) = Adz_BC_LAM_mask_i(i0:in,j0:jn,1:l_nk)

         stkpntr(2*(nf-1) + 1) = c_loc(adv_o(1,1,1,nf))
         stkpntr(2*(nf-1) + 2) = c_loc(adv_i(1,1,1,nf))

      end if

      if (nf==0) return

      !Estimate FLUX_out/FLUX_in using TRICUB semi-Lagrangian advection
      !----------------------------------------------------------------
      nptr_F = nf

      ni = Adz_inb-Adz_i0b+1
      nj = Adz_jnb-Adz_j0b+1
      nij= ni*nj
      n= nj
      slice = n*ni

      n1=1
      do while (n1 <= Adz_num_b)
         n2= min(n1+slice-1,Adz_num_b)
         np= n2-n1+1
         k1= n1/nij + Adz_k0
         k2= n2/nij + min(1,mod(n2,nij)) + Adz_k0 - 1
         ij= 0
         call tricublin_zyx1_m_n ( wrkc, stkpntr, Adz_pb((n1-1)*3+1,Adz_i0b,Adz_j0b,Adz_k0),&
                                   Adz_cpntr_t, np, nptr_F*2 )
         do nf=1,nptr_F
         do nn=1,2
         do k=k1,k2
            kk= (k-Adz_k0)*nj
            j1=max(1 ,n1/ni+1-kk) + Adz_j0b - 1
            j2=min(nj,n2/ni  -kk) + Adz_j0b - 1
            do j=j1,j2
            do i=Adz_i0b,Adz_inb
               ij= ij+1
               Adz_flux(i,j,k,(tr_num(nf)-1)*2+nn)= wrkc(ij)
            end do
            end do
         end do
         end do
         end do
         n1=n2+1
      end do

      !Do the appropriate ZEROING to FLUX_out/FLUX_in
      !----------------------------------------------
      do nf=1,nptr_F

         do k=1,l_nk
            do j=1+pil_s,l_nj-pil_n
               do i=1+pil_w,l_ni-pil_e
                  Adz_flux(i,j,k,(tr_num(nf)-1)*2+1) = 0.
               end do
            end do
         end do

         if (l_north) then
            Adz_flux(1:l_ni,l_nj-pil_n+1:l_nj,1:l_nk,(tr_num(nf)-1)*2+2) = 0.
         end if

         if (l_east) then
            Adz_flux(l_ni-pil_e+1:l_ni,1:l_nj,1:l_nk,(tr_num(nf)-1)*2+2) = 0.
         end if

         if (l_south) then
            Adz_flux(1:l_ni,1:pil_s,1:l_nk,(tr_num(nf)-1)*2+2) = 0.
         end if

         if (l_west) then
            Adz_flux(1:pil_w,1:l_nj,1:l_nk,(tr_num(nf)-1)*2+2) = 0.
         end if

      end do
!
!---------------------------------------------------------------------
!
      return
      end
