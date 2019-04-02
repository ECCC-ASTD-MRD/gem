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

module fislh_inp_base
   use iso_c_binding
   use vertical_interpolation
   use vGrid_Descriptors
   use inp_options
   use dyn_fisl_options
   use geomh
   use glb_ld
   use inp_base, only : inp_read, inp_read_uv, inp_3dpres
   use inp_mod
   use tdpack
   implicit none
#include <arch_specific.hf>

   public

contains

!**s/r inp_get_H - Read variable F_var_S valid at F_datev, perform hor.
!                interpolation to F_hgrid_S hor. grid and perform
!                vertical interpolation to F_vgrid_S vertical grid

      integer function inp_get_H ( F_var_S, F_hgrid_S, F_ip1_dst,      &
                                   F_ip1_gz_src, F_gz_src, F_vgd_dst , &
                                   F_sfc_dst, F_sfcLS_dst,             &
                                   F_dest , Minx, Maxx, Miny, Maxy    ,&
                                   F_nk, F_inttype_S )

      character(len=*)          , intent(in) :: F_var_S,F_hgrid_S
      character(len=*), optional, intent(in) :: F_inttype_S
      integer                   , intent(in) :: Minx,Maxx,Miny,Maxy, F_nk
      integer, dimension(:)  , pointer, intent(in) :: F_ip1_dst, F_ip1_gz_src
      real, dimension (:,:,:), pointer, intent(in) :: F_gz_src
      type(vgrid_descriptor) , intent(IN) :: F_vgd_dst
      real, dimension (:,:), pointer, intent(in) :: F_sfc_dst, F_sfcLS_dst
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out):: F_dest

!     local variables
      character(len=12) :: inttype
      integer :: nka
      integer, dimension (:    ), pointer :: ip1_list
      real   , dimension (:,:  ), pointer :: dummy
      real   , dimension (:,:,:), pointer :: wrkr,srclev,dstlev
      integer, parameter :: INP_OK = 0, INP_ERROR = -1
!
!---------------------------------------------------------------------
!
      inp_get_h= -1
      if ( any (Inp_blacklist_S(1:MAX_blacklist) == trim(F_var_S)) ) return

      nullify (ip1_list, wrkr, dummy)

      inp_get_h= inp_read ( F_var_S, F_hgrid_S, wrkr, ip1_list, nka )

      if (inp_get_h .lt. 0) then
         if (associated(ip1_list)) deallocate (ip1_list)
         if (associated(wrkr    )) deallocate (wrkr    )
         return
      end if

      allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,nka), dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk) )
      if ( inp_match_heights(srclev, F_gz_src, ip1_list, F_ip1_gz_src, nka, size(F_gz_src,dim=3)) &
           == INP_ERROR) call gem_error (-1,'inp_get_h','See error above.')
      call inp_3dpres ( F_vgd_dst, F_ip1_dst, F_sfc_dst, F_sfcLS_dst, dstlev, 1, G_nk)

      inttype= 'cubic'
      if (present(F_inttype_S)) inttype= F_inttype_S
      call vertint2 ( F_dest,dstlev,G_nk, wrkr, srclev,nka       , l_minx,l_maxx,l_miny,l_maxy, 1&
           &,l_ni,1,l_nj, varname=F_var_S, inttype= inttype)

      !open(unit=62,file='source.txt',status='unknown')
      !open(unit=63,file='desti.txt',status='unknown')
      !do k=1,nka
      !   write(62,*)srclev(1,1,k),wrkr(1,1,k)
      !end do
      !do k=1,G_nk
      !   write(63,*)dstlev(1,1,k),F_dest(1,1,k)
      !end do
      !close(62)
      !close(63)

      deallocate (ip1_list,wrkr,srclev,dstlev)
!
!---------------------------------------------------------------------
!
      return
      end function inp_get_H
!
!---------------------------------------------------------------------
!
      subroutine inp_3dhgts ( F_vgd, F_ip1, F_sfc, F_sfcL, F_dest, k0, kn)

      type(vgrid_descriptor) , intent(IN) :: F_vgd
      integer                , intent(IN) :: k0,kn
      integer, dimension(:)    , pointer, intent(in ) :: F_ip1
      real   , dimension(:,:  ), pointer, intent(in ) :: F_sfc,F_sfcL
      real   , dimension(:,:,:), pointer, intent(out) :: F_dest

!     local variables
      integer :: istat
      real, dimension (:,:  ), pointer :: hgt,hgtls
      real, dimension (:,:,:), pointer :: ptr3d
!
!---------------------------------------------------------------------
!
      hgt  => F_sfc (1:l_ni,1:l_nj      )
      ptr3d => F_dest(1:l_ni,1:l_nj,k0:kn)

      hgtls=> F_sfcL (1:l_ni,1:l_nj)
      istat= vgd_levels (F_vgd, F_ip1(k0:kn), ptr3d, sfc_field=hgt, sfc_field_ls=hgtls)
!
!---------------------------------------------------------------------
!
      return
    end subroutine inp_3dhgts

!**s/r inp_hwnd_H - Read horizontal winds UU,VV F valid at F_datev
!                 and perform vectorial horizontal interpolation
!                 and  vertical interpolation to momentum levels in heights

      subroutine inp_hwnd_H(F_u,F_v, F_gz, F_gz_ip1_list, F_vgd_dst, F_stag_L, &
                            F_sfc_dst, F_sfcLS_dst,                            &
                            Minx, Maxx, Miny, Maxy, F_nk )

        use ver

        implicit none
        integer                , intent(IN)  :: Minx,Maxx,Miny,Maxy, F_nk
        integer, dimension(:), pointer, intent(in) :: F_gz_ip1_list
        type(vgrid_descriptor) , intent(IN)  :: F_vgd_dst
        real, dimension (:,:), pointer, intent(in) :: F_sfc_dst, F_sfcLS_dst
        real, dimension (:,:,:), pointer :: F_gz
        real, dimension(Minx:Maxx,Miny:Maxy,F_nk), intent(out) :: F_u, F_v
        logical :: F_stag_L
        integer, parameter :: INP_OK = 0, INP_ERROR = -1

        ! Local variables
        integer :: nka
        integer, dimension (:    ), pointer :: ip1_list, ip1_target
        real   , dimension (:,:,:), pointer :: srclev,dstlev
        real   , dimension (:,:,:), pointer :: ur,vr

        nullify (ip1_list,ip1_target,srclev,dstlev,ur,vr)

        if (F_stag_L) then
           call inp_read_uv ( ur, vr, 'UV' , ip1_list, nka )
        else
           call inp_read_uv ( ur, vr, 'Q ' , ip1_list, nka )
        end if
        allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,nka) ,&
             dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk) )
        allocate (ip1_target(1:G_nk))
        ip1_target(1:G_nk)= Ver_ip1%m(1:G_nk)

        if ( inp_match_heights(srclev, F_gz, ip1_list, F_gz_ip1_list, nka, size(F_gz,dim=3)) &
             == INP_ERROR) call gem_error (-1,'inp_hwnd_H','See error above.')

        call inp_3dhgts ( F_vgd_dst, ip1_target, F_sfc_dst, F_sfcLS_dst,&
             dstlev, 1, G_nk)

        call vertint2 ( F_u,dstlev,G_nk, ur,srclev,nka,&
             l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )
        call vertint2 ( F_v,dstlev,G_nk, vr,srclev,nka,&
             l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )

      !open(unit=62,file='source.txt',status='unknown')
      !open(unit=63,file='desti.txt',status='unknown')
      !do k=1,nka
      !   write(62,*)srclev(1,1,k),sqrt(ur(1,1,k)*ur(1,1,k)+vr(1,1,k)*vr(1,1,k))
      !end do
      !do k=1,G_nk
      !   write(63,*)dstlev(1,1,k),sqrt(F_u(1,1,k)*F_u(1,1,k)+F_v(1,1,k)*F_v(1,1,k))
      !end do
      !close(62)
      !close(63)

      deallocate (ip1_list,ip1_target,ur,vr,srclev,dstlev)
      nullify    (ip1_list,ip1_target,ur,vr,srclev,dstlev)

      end subroutine inp_hwnd_H

!**s/r inp_match_heights - To match heights
!

      integer function inp_match_heights(F_ho, F_hi, F_ip1_list_o, &
                                         F_ip1_list_i, nko, nki) result(status)

      real, dimension (:,:,:), pointer, intent(out) :: F_ho
      real, dimension (:,:,:), pointer, intent(in)  :: F_hi
      integer, dimension (:)    , pointer, intent(out) :: F_ip1_list_o
      integer, dimension (:)    , pointer, intent(in)  :: F_ip1_list_i
      integer, intent(inout) :: nko
      integer, intent(in)    :: nki
      integer, parameter :: INP_OK = 0, INP_ERROR = -1

!     Local variables
      integer :: index_diag_AGL, ko, ki, kind, i, j, k
      logical :: found_L
      real :: lev
      character(len=128) :: message_S

      status = INP_ERROR

      index_diag_AGL = -1

      do ko=1, nko
         found_L = .false.
         call convip (F_ip1_list_o(ko), lev, kind, -1, message_S, .false.)
         if( kind .eq. 4 )then
            found_L = .true.
            ! Pour le moment sort_ip1 ne garde que le dernier kind = 4 de la liste.
            ! Donc cette erreur ne se produira pas.
            if( index_diag_AGL /= -1 )&
                 call gem_error (-1,'inp_match_heights','more than one diagnostic level, review code')
            index_diag_AGL = ko
            F_ho(:,:,ko) = F_hi(:,:,nki) + lev
         else
            do ki=1, nki
               if (F_ip1_list_o(ko) == F_ip1_list_i(ki)) then
                   found_L = .true.
                   F_ho(:,:,ko) = F_hi(:,:,ki)
                   exit
               end if
            end do
         end if
         if(.not. found_L)then
            write(message_S,*)'Missing field: GZ for ip1 = ',&
                 F_ip1_list_o(ko)
            call gem_error (-1,'inp_match_heights',message_S)
         end if
      end do

      if(index_diag_AGL /= -1)then
         ! Check monotonicity
         found_L = .false.
         do j = 1, l_nj
            do i = 1, l_ni
               if( F_ho(i,j, index_diag_AGL-1) <= F_ho(i,j, index_diag_AGL) )then
                  found_L = .true.
               end if
            end do
         end do
         if(index_diag_AGL < nko)then
            ! Check level below
            do j = 1, l_nj
               do i = 1, l_ni
                  if( F_ho(i,j, index_diag_AGL+1) >= F_ho(i,j, index_diag_AGL) )then
                     found_L = .true.
                  end if
               end do
            end do
         end if
         if(found_L)then
            ! Remove diag level
            do k = index_diag_AGL, nko-1
               F_ho(:,:,k) = F_ho(:,:,k+1)
            end do
            nko = nko - 1
         end if
      end if

      status = INP_OK
      return

    end function inp_match_heights

    integer function inp_comp_pres_from_vt_p0(F_q, F_ps, F_vt, F_zs, F_zsls, F_vgd, &
         Minx,Maxx,Miny,Maxy,F_nk) result(istat)
      integer                                    , intent(in)  :: Minx,Maxx,Miny,Maxy, F_nk
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk+1), intent(out) :: F_q
      real, dimension(Minx:Maxx,Miny:Maxy,F_nk)  , intent(in)  :: F_vt
      real   , dimension(:,:), pointer           , intent(in ) :: F_ps, F_zs,F_zsls
      type (vgrid_descriptor) :: F_vgd
      ! Local varibales
      integer :: i, j, k
      integer, dimension(:), pointer :: ip1_list
      real, dimension(:,:,:), pointer :: zz
      real*8 :: aaa
      integer, parameter :: INP_OK = 0, INP_ERROR = -1

      nullify(zz,ip1_list)

      istat = INP_ERROR

      ! TODO, trap error
      if( vgd_get(F_vgd,'VIPM - level ip1 list (m)',ip1_list) == VGD_ERROR)then
         return
      end if

      allocate(zz(l_minx:l_maxx,l_miny:l_maxy,size(ip1_list)))
      call inp_3dhgts ( F_vgd, ip1_list, F_zs, F_zsls, zz, 1, size(ip1_list))
      ! Calculate q=RTstar*log(surface pressure/1.e5) for level G_nk + 1
      aaa=rgasd_8*Cstv_tstr_8
      do j=1,l_nj
      do i=1,l_ni
         F_q(i,j,G_nk+1)=aaa*log(F_ps(i,j)/1.e5)
      end do
      end do
      ! Integrate hydrostatic equation, dq/dz=-gTstar/Tv, to obtain q at momentum levels
      aaa = grav_8*Cstv_tstr_8
      do k=G_nk,1,-1
         do j=1,l_nj
         do i=1,l_ni
            F_q(i,j,k) = F_q(i,j,k+1) + aaa*(zz(i,j,k+1) - zz(i,j,k))/F_vt(i,j,k)
         end do
         end do
      end do
      ! Add gz to obtain qprime, stored in F_q
      F_q(1:l_ni,1:l_nj,1:G_nk+1)=F_q(1:l_ni,1:l_nj,1:G_nk+1)+grav_8*zz(1:l_ni,1:l_nj,1:G_nk+1)
      deallocate(zz,ip1_list)
      istat = INP_OK
      return
    end function inp_comp_pres_from_vt_p0

end module fislh_inp_base
