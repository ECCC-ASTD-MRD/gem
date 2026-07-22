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
      subroutine OUTs_fstecr ( F_levels, F_indo, F_nk, F_k0, &
                             F_kn, F_offset, F_stag_S)
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use mpi_f08
      use IOs
      use OUTs
      use vGrid_Descriptors
      use rmn_fst24
      implicit none

      character(len=4), intent(IN) :: F_stag_S
      integer, intent(IN) :: F_nk, F_k0, F_kn, F_offset
      integer, intent(IN) :: F_indo(F_nk)
      real   , intent(IN) :: F_levels(*)

#include <rmnlib_basics.hf>

      character(len=1) :: gtyp
      character(len=12) :: etiket_S
      integer :: i,j,k,nis,njs,wk_njs,dtyp,ip1,ip2,ip3,err
      real, save, dimension (:,:,:), pointer :: wk2,wk2_uu
      real, dimension (:,:,:), pointer :: wk, wk_ptr
      real, dimension (:,:), pointer :: x,y,lat,lon
      real, save, dimension (:), pointer :: wk_target
      logical :: vector_int_L,wait_for_v_L
      logical, save :: done_L=.false.
      logical :: success
      type(wind_rec), save :: windu

!     
!--------------------------------------------------------------------
!
      if(F_k0 > F_kn)return
      if(.not.done_L)then
         done_L=.true.
         call OUTS_init_wind_rec(windu)
         ! Note: wk2, wk2_uu are reused if possible
         ! see subroutine OUTs_hor_int
         nullify(wk2,wk2_uu,wk_target)
      endif

      nullify(wk,wk_ptr)
      nis = Out_in - Out_i0 + 1
      njs = Out_jn - Out_j0 + 1
      wk_njs = njs ;  gtyp='Z'

      if (Grd_yinyang_L) then
         wk_njs = njs*2
         gtyp='U'
      endif
      if ( (nis < 1) .or. (njs < 1) ) return
      
      call OUTs_check_allocation_3D(wk,wk_target,nis,wk_njs,F_k0,F_kn)
      if (Grd_yinyang_L) then
         do k=F_k0, F_kn
            wk(1:nis,1:njs,k) = &
            IOs_glbdata(Out_i0:Out_in,Out_j0:Out_jn,k,1)
            wk(1:nis,njs+1:wk_njs,k) = &
            IOs_glbdata(Out_i0:Out_in,Out_j0:Out_jn,k,2)
         end do
      else
         wk(1:nis,1:njs,F_k0:F_kn) = &
         IOs_glbdata(Out_i0:Out_in,Out_j0:Out_jn,F_k0:F_kn,1)
      endif
      if (Out_nbit <= 16) then
         dtyp = 134
      else
         dtyp = 133
      endif
      ip1 = 0
      etiket_S=Out_etik_S
      wk_ptr=>wk
      vector_int_L=.false.
      if(OUTs_hor_int_L(F_stag_S))then
         if(.not.OUTs_hor_int(wk2,wk2_uu,wk,Out_nomvar,F_stag_S,windu,nis,wk_njs,F_k0,F_kn,&
              F_levels(F_k0),gtyp,Out_ig1,Out_ig2,Out_ig3,Out_ig4,etiket_S,vector_int_L,wait_for_v_L))return
         if(wait_for_v_L)then            
            ! Waiting for V record
            return
         endif
         wk_ptr=>wk2
      endif
      Out_rec%typvar=Out_typvar_S
      Out_rec%grtyp=gtyp
      Out_rec%data_type=dtyp
      Out_rec%etiket=etiket_S
      Out_rec%dateo=Out_dateo
      Out_rec%deet=Out_deet
      Out_rec%npas=Out_npas
      Out_rec%data_bits=Out_nbit
      Out_rec%pack_bits=Out_nbit
      Out_rec%nomvar=Out_nomvar
      Out_rec%ni=nis
      Out_rec%nj=wk_njs
      Out_rec%nk=1
      Out_rec%ig1=Out_ig1
      Out_rec%ig2=Out_ig2
      Out_rec%ig3=Out_ig3
      Out_rec%ig4=Out_ig4
      
      if ( F_stag_S(2:2) == 'S' ) then
         if (F_nk > 1) then
            do k=F_k0, F_kn
               call encode_ips (ip1,ip2,ip3,F_levels(k))
              ! err = fstecr ( wk_ptr(1,1,k),wk,-Out_nbit,Out_unf,&
              ! Out_dateo,Out_deet,Out_npas,nis,wk_njs,1  ,&
              ! ip1,ip2,ip3, Out_typvar_S,Out_nomvar,etiket_S,gtyp,&
              ! Out_ig1,Out_ig2,Out_ig3,Out_ig4,dtyp,Out_rewrit_L )

               Out_rec%ip1=ip1
               Out_rec%ip2=ip2
               Out_rec%ip3=ip3
               Out_rec%data=c_loc(wk_ptr(1,1,k))
               
               success = Out_file%write(Out_rec,rewrite=FST_SKIP)
            end do
         else
            if ((F_k0==1).and.(F_kn==1)) then
               call encode_ips (ip1,ip2,ip3,0.)
             !  err = fstecr ( wk_ptr,wk,-Out_nbit,Out_unf,&
             !  Out_dateo,Out_deet,Out_npas,nis,wk_njs,1  ,&
             !  ip1,ip2,ip3, Out_typvar_S,Out_nomvar,etiket_S,gtyp,&
             !  Out_ig1,Out_ig2,Out_ig3,Out_ig4,dtyp,Out_rewrit_L )

               Out_rec%data=c_loc(wk_ptr)
               Out_rec%ip1=ip1
               Out_rec%ip2=ip2
               Out_rec%ip3=ip3
               success = Out_file%write(Out_rec,rewrite=FST_SKIP)
            endif
         endif
      else
         do k=F_k0, F_kn
            call encode_ips (ip1,ip2,ip3,F_levels(F_indo(k)))
            !err = fstecr ( wk_ptr(1,1,k),wk,-Out_nbit,Out_unf,&
            !      Out_dateo,Out_deet,Out_npas,nis,wk_njs,1  ,&
            !      ip1,ip2,ip3, Out_typvar_S,Out_nomvar,etiket_S,gtyp,&
            !      Out_ig1,Out_ig2,Out_ig3,Out_ig4,dtyp,Out_rewrit_L )

            Out_rec%data=c_loc(wk_ptr(1,1,k))
            Out_rec%ip1=ip1
            Out_rec%ip2=ip2
            Out_rec%ip3=ip3
            success = Out_file%write(Out_rec,rewrite=FST_SKIP)
            if(vector_int_L)then
              ! err = fstecr ( wk2_uu(1,1,k),wk,-Out_nbit,Out_unf,&
              !      Out_dateo,Out_deet,Out_npas,nis,wk_njs,1  ,&
              !      ip1,ip2,ip3, Out_typvar_S,'UU',etiket_S,gtyp,&
              !      Out_ig1,Out_ig2,Out_ig3,Out_ig4,dtyp,Out_rewrit_L )

               Out_rec%data=c_loc(wk2_uu(1,1,k))
               Out_rec%nomvar='UU'
               success = Out_file%write(Out_rec,rewrite=FST_SKIP)
               ! Reset nomvar for next level, in this context Out_nomvar='VV'
               Out_rec%nomvar=Out_nomvar
            endif
         end do
      endif
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_fstecr
      
      subroutine OUTs_fstecr_sfc ( F_vsfc, F_i0, F_in, F_stag_S )
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use IOs
      use OUTs
      implicit none

      integer, intent(IN) :: F_i0, F_in
      type (sfc_var), intent(IN) :: F_vsfc(*)
      character(len=4), intent(IN) :: F_stag_S
#include <rmnlib_basics.hf>
      
      character(len=1) :: gtyp
      character(len=12) :: etiket_S
      integer :: j,k,nis,njs,wk_njs,dtyp,err,ip1,ip2,ip3,h_stag
      real, save, dimension (:,:), pointer :: wk2,wk_ptr
      real, dimension (:,:,:), pointer :: wk
      real, save, dimension (:), pointer :: wk_target
      logical, save :: done_L=.false.
      logical :: wind_cycke_L
      logical :: success
!     
!--------------------------------------------------------------------
!
      nis = Out_in - Out_i0 + 1
      njs = Out_jn - Out_j0 + 1
      wk_njs = njs ;  gtyp='Z'      
      
      if (Grd_yinyang_L) then
         wk_njs = njs*2
         gtyp='U'
      endif
      
      if ( (nis < 1) .or. (njs < 1) ) return
      
      if(.not.done_L)then
         done_L=.true.
         ! Note: wk2 is reused if possible
         ! see subroutine OUTs_hor_int_sfc
         nullify(wk2,wk_target)
      endif

      call OUTs_check_allocation_3D(wk,wk_target,nis,wk_njs,F_i0,F_in)
      
      if (Grd_yinyang_L) then
         do k=F_i0, F_in
            wk(1:nis,1:njs,k) = &
            IOs_glbdata(Out_i0:Out_in,Out_j0:Out_jn,k,1)
            wk(1:nis,njs+1:wk_njs,k) = &
            IOs_glbdata(Out_i0:Out_in,Out_j0:Out_jn,k,2)
         end do
      else
         wk(1:nis,1:njs,F_i0:F_in) = &
         IOs_glbdata(Out_i0:Out_in,Out_j0:Out_jn,F_i0:F_in,1)
      endif

      etiket_S=Out_etik_S

      if (Out_nbit <= 16) then
         dtyp = 134
      else
         dtyp = 133
      endif      

      Out_rec%typvar=Out_typvar_S
      Out_rec%data_type=dtyp
      Out_rec%dateo=Out_dateo
      Out_rec%deet=Out_deet
      Out_rec%npas=Out_npas
      Out_rec%data_bits=Out_nbit

      do k= F_i0, F_in
         wk_ptr=>wk(:,:,k)
         Out_kind=F_vsfc(k)%knd
         if(OUTs_hor_int_L(F_vsfc(k)%stag))then
            ! Find matching horizontal
            read(F_vsfc(k)%stag(4:4),*)h_stag
            do j=1,size(hgd_usr)
               if(hgd_usr(j)%usr_grid_index == h_stag)then               
                  call OUTs_wrtref_usr (F_stag_S,hgd_usr(j))
                  exit
               end if
            end do
            if(.not.OUTs_hor_int_sfc(wk2,wk_ptr,F_vsfc(k),nis,wk_njs,&
                 gtyp,Out_ig1,Out_ig2,Out_ig3,Out_ig4,etiket_S,wind_cycke_L))return
            if(wind_cycke_L)cycle
            wk_ptr=>wk2
         else
            call OUTs_wrtref ( F_vsfc(k)%stag )
         endif
         Out_rec%ig1=Out_ig1
         Out_rec%ig2=Out_ig2
         Out_rec%ig3=Out_ig3
         Out_rec%ig4=Out_ig4
         call encode_ips (ip1,ip2,ip3,F_vsfc(k)%lvl)
         !err = fstecr ( wk_ptr(1,1),wk,-F_vsfc(k)%nbits,Out_unf,&
         !          Out_dateo,Out_deet,Out_npas,nis,wk_njs,1  ,&
         !          ip1,ip2,ip3 ,&
         !          Out_typvar_S,F_vsfc(k)%nv,etiket_S,gtyp ,&
         !          Out_ig1,Out_ig2,Out_ig3,Out_ig4,dtyp,Out_rewrit_L )
         Out_rec%data=c_loc(wk_ptr(1,1))
         Out_rec%nomvar=F_vsfc(k)%nv
         Out_rec%etiket=etiket_S               
         Out_rec%pack_bits=F_vsfc(k)%nbits
         Out_rec%ni=nis
         Out_rec%nj=wk_njs
         Out_rec%nk=1
         Out_rec%ip1=ip1
         Out_rec%ip2=ip2
         Out_rec%ip3=ip3
         Out_rec%grtyp=gtyp
         success = Out_file%write(Out_rec,rewrite=FST_SKIP)
        
      end do
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_fstecr_sfc

      subroutine encode_ips (ip1,ip2,ip3,rf)
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use mpi_f08
      use OUTs
      use vGrid_Descriptors
      implicit none

      real   , intent(IN) :: rf
      integer, intent(OUT) :: ip1,ip2,ip3
      
      include 'rmn/convert_ip123.inc'

      character(len=8) dumc
      integer :: lstep,modeip1,err,kind
      type(FLOAT_IP) :: RP1,RP2,RP3
      real lvl
!     
!--------------------------------------------------------------------
!
      modeip1= 1
      if (Out_kind == 1 .or. Out_kind == 2) modeip1= 3 !old ip1 style for pressure and eta lvls output
      lstep= -1
      if ( lstep > 0 ) then
         RP2%lo  = dble(lstep          ) * dble(Out_deet) / 3600.d0
         RP2%hi  = dble(max(0,Out_npas)) * dble(Out_deet) / 3600.d0
         RP2%kind= KIND_HOURS
         RP3%lo= 0. ; RP3%hi= 0. ; RP3%kind= 0
      end if
      if ( lstep > 0 ) then
         RP1%lo  = rf
         RP1%hi  = RP1%lo
         RP1%kind= Out_kind
         err= encode_ip ( ip1,ip2,ip3,RP1,RP2,RP3 )
      else
         ip2 = Out_ip2
         if (Out_ip2>32767) call convip( ip2, real(Out_ip2), 10, +2, dumc,.false.)
         ip3 = Out_ip3
         lvl = rf ;kind=Out_kind
         call convip ( ip1, lvl, kind, modeip1,dumc,.false. )
      end if
!
!--------------------------------------------------------------------
!
    end subroutine encode_ips
