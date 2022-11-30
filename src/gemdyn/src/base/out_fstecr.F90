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

!**s/r out_fstecr

      subroutine out_fstecr ( fa,lminx,lmaxx,lminy,lmaxy,rf,nomvar,&
                              mul,add,kind,lstep,nkfa,ind_o,nk_o  ,&
                              nbit,F_empty_stk_L )
      use ISO_C_BINDING
      use step_options
      use gem_options
      use glb_ld
      use out_mod
      use out_meta
      implicit none
#include <arch_specific.hf>

      character(len=*) nomvar
      logical F_empty_stk_L
      integer lminx,lmaxx,lminy,lmaxy,nkfa,nbit,nk_o,kind,lstep
      integer ind_o(nk_o)
      real fa (lminx:lmaxx,lminy:lmaxy,nkfa), rf(nkfa), mul,add

      include 'rmn/convert_ip123.inc'

Interface
subroutine out_stkecr ( fa,lminx,lmaxx,lminy,lmaxy,meta,nplans, &
                         g_id,g_if,g_jd,g_jf )
      use out_meta
      integer lminx,lmaxx,lminy,lmaxy,nplans
      integer g_id,g_if,g_jd,g_jf
      real fa(lminx:lmaxx,lminy:lmaxy,nplans)
      type (meta_fstecr), dimension(:), pointer :: meta
End Subroutine out_stkecr
End Interface

      type(FLOAT_IP) :: RP1,RP2,RP3
      real,   parameter :: eps=1.e-12
      character(len=8) dumc
      logical, save :: done= .false.
      integer modeip1,i,j,k,err
      integer, save :: istk = 0
      real, save, dimension (:,:,:), pointer :: f2c => null()
      type (meta_fstecr), save, dimension(:), pointer :: meta => null()
!
!----------------------------------------------------------------------
!
      if (.not.associated(meta)) then
         allocate (meta(out_stk_size))
      end if

      if (.not.associated(f2c )) then
         allocate (f2c(l_minx:l_maxx,l_miny:l_maxy,out_stk_size))
      end if

      if (.not.done) then
         out_stk_full = 0 ; out_stk_part = 0 ; done = .true.
      end if

      if (F_empty_stk_L) then
         if ( istk > 0) then
            call out_stkecr ( f2c,l_minx,l_maxx,l_miny,l_maxy ,&
                               meta,istk, Out_gridi0,Out_gridin,&
                                          Out_gridj0,Out_gridjn )
            if (istk < out_stk_size) out_stk_part= out_stk_part + 1
         end if
         istk= 0
         deallocate (meta, f2c) ; nullify (meta, f2c)
         return
      end if

      modeip1= 1
      if (kind == 2) modeip1= 3 !old ip1 style for pressure lvls output

      if ( lstep > 0 ) then
         RP2%lo  = dble(lstep           ) * dble(Step_dt) / 3600.d0
         RP2%hi  = dble(max(0,lctl_step)) * dble(Step_dt) / 3600.d0
         RP2%kind= KIND_HOURS
         RP3%lo= 0. ; RP3%hi= 0. ; RP3%kind= 0
      end if

      do k= 1, nk_o
         istk = istk + 1
         do j= 1, l_nj
         do i= 1, l_ni
            f2c(i,j,istk)= fa(i,j,ind_o(k))*mul + add
         end do
         end do
         if ( lstep > 0 ) then
            RP1%lo  = rf(ind_o(k))
            RP1%hi  = RP1%lo
            RP1%kind= kind
            err= encode_ip ( meta(istk)%ip1,meta(istk)%ip2,&
                             meta(istk)%ip3,RP1,RP2,RP3 )
         else
            meta(istk)%ip2= Out_ip2
            if (Out_ip2>32767) call convip( meta(istk)%ip2, real(Out_ip2), 10, +2, dumc,.false.)
            meta(istk)%ip3= Out_ip3
            call convip ( meta(istk)%ip1, rf(ind_o(k)), kind,&
                          modeip1,dumc,.false. )
         end if
         meta(istk)%nv   = nomvar
         meta(istk)%ig1  = Out_ig1
         meta(istk)%ig2  = Out_ig2
         meta(istk)%ig3  = Out_ig3
         meta(istk)%ni   = Out_gridin - Out_gridi0 + 1
         meta(istk)%nj   = Out_gridjn - Out_gridj0 + 1
         meta(istk)%nbits= nbit
         if (nbit <= 16) then
            meta(istk)%dtyp = 134
         else
            meta(istk)%dtyp = 133
         endif
         if (istk == out_stk_size) then
            call out_stkecr ( f2c,l_minx,l_maxx,l_miny,l_maxy ,&
                               meta,istk, Out_gridi0,Out_gridin,&
                                          Out_gridj0,Out_gridjn )
            istk=0
            out_stk_full= out_stk_full + 1
         end if
      end do
!
!--------------------------------------------------------------------
!
      return
      end
