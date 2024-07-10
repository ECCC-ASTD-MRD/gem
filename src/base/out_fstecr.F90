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
      implicit none
#include <arch_specific.hf>

      character(len=*) nomvar
      logical F_empty_stk_L
      integer lminx,lmaxx,lminy,lmaxy,nkfa,nbit,nk_o,kind,lstep
      integer ind_o(nk_o)
      real fa (lminx:lmaxx,lminy:lmaxy,nkfa), rf(nkfa), mul,add

      include 'rmn/convert_ip123.inc'

Interface
subroutine out_stkecr ( fa,lminx,lmaxx,lminy,lmaxy,metaf,nplans, &
                         g_id,g_if,g_jd,g_jf )
      use out_mod
      use rmn_fst24
      integer lminx,lmaxx,lminy,lmaxy,nplans
      integer g_id,g_if,g_jd,g_jf
      real fa(lminx:lmaxx,lminy:lmaxy,nplans)
      type (fst_record), dimension(:), pointer :: metaf
End Subroutine out_stkecr
End Interface

      type(FLOAT_IP) :: RP1,RP2,RP3
      real,   parameter :: eps=1.e-12
      character(len=8) dumc
      logical, save :: done= .false.
      integer modeip1,i,j,k,err
      integer, save :: istk = 0
      real, save, dimension (:,:,:), pointer :: f2c => null()
      type (fst_record), save, dimension(:), pointer :: metaf => null()
!
!----------------------------------------------------------------------
!
      if (.not.associated(metaf)) then
         allocate (metaf(out_stk_size))
      end if

      if (.not.associated(f2c )) then
         allocate (f2c(l_minx:l_maxx,l_miny:l_maxy,out_stk_size))
      end if

      if (.not.done) then
         Out_stk_full = 0 ; Out_stk_part = 0 ; done = .true.
      end if

      if (F_empty_stk_L) then
         if ( istk > 0) then
            call out_stkecr ( f2c,l_minx,l_maxx,l_miny,l_maxy ,&
                               metaf,istk, Out_gridi0,Out_gridin,&
                                          Out_gridj0,Out_gridjn )
            if (istk < out_stk_size) out_stk_part= out_stk_part + 1
         end if
         istk= 0
         deallocate (metaf, f2c) ; nullify (metaf, f2c)
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
            err= encode_ip ( metaf(istk)%ip1,metaf(istk)%ip2,&
                             metaf(istk)%ip3,RP1,RP2,RP3 )
         else
            metaf(istk)%ip2= Out_ip2
            if (Out_ip2>32767) call convip( metaf(istk)%ip2, real(Out_ip2), 10, +2, dumc,.false.)
            metaf(istk)%ip3= Out_ip3
            call convip ( metaf(istk)%ip1, rf(ind_o(k)), kind,&
                          modeip1,dumc,.false. )
         end if
         metaf(istk)%nomvar= nomvar
         metaf(istk)%ig1  = Out_ig1
         metaf(istk)%ig2  = Out_ig2
         metaf(istk)%ig3  = Out_ig3
         metaf(istk)%ni   = Out_gridin - Out_gridi0 + 1
         metaf(istk)%nj   = Out_gridjn - Out_gridj0 + 1
         metaf(istk)%data_bits= 32
         metaf(istk)%pack_bits= nbit
         if (nbit <= 16) then
            metaf(istk)%data_type = 134
         else
            metaf(istk)%data_type = 133
         endif
         if (istk == out_stk_size) then
            call out_stkecr ( f2c,l_minx,l_maxx,l_miny,l_maxy ,&
                               metaf,istk, Out_gridi0,Out_gridin,&
                                          Out_gridj0,Out_gridjn )
            istk=0
            Out_stk_full= Out_stk_full + 1
         end if
      end do
!
!--------------------------------------------------------------------
!
      return
      end
