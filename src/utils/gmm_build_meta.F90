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
subroutine gmm_build_meta1D(meta,lowx,highx,halox,halomaxx,nx,initmode,flags)
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use rmn_gmm
  implicit none

  type(gmm_metadata), intent(out) :: meta
  integer :: lowx,highx,halox,halomaxx,nx
  integer :: initmode,flags

  meta = GMM_NULL_METADATA
  meta%l(1) = gmm_layout(lowx,highx,halox,halomaxx,nx)
  meta%a%initmode = initmode
  meta%a%flags = flags

  return
end subroutine gmm_build_meta1D

subroutine gmm_build_meta2D(meta,lowx,highx,halox,halomaxx,nx,lowy,highy,haloy,&
     halomaxy,ny,initmode,flags)
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use rmn_gmm
  implicit none

  type(gmm_metadata), intent(out) :: meta
  integer :: lowx,highx,halox,halomaxx,nx
  integer :: lowy,highy,haloy,halomaxy,ny
  integer :: initmode,flags

  meta = GMM_NULL_METADATA
  meta%l(1) = gmm_layout(lowx,highx,halox,halomaxx,nx)
  meta%l(2) = gmm_layout(lowy,highy,haloy,halomaxy,ny)
  meta%a%initmode = initmode
  meta%a%flags = flags

  return
end subroutine gmm_build_meta2D

subroutine gmm_build_meta3D(meta,lowx,highx,halox,halomaxx,nx,lowy,highy,haloy,halomaxy,ny,&
                            lowz,highz,haloz,halomaxz,nz,initmode,flags)
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use rmn_gmm
  implicit none

  type(gmm_metadata), intent(out) :: meta
  integer :: lowx,highx,halox,halomaxx,nx
  integer :: lowy,highy,haloy,halomaxy,ny
  integer :: lowz,highz,haloz,halomaxz,nz
  integer :: initmode,flags

  meta = GMM_NULL_METADATA
  meta%l(1) = gmm_layout(lowx,highx,halox,halomaxx,nx)
  meta%l(2) = gmm_layout(lowy,highy,haloy,halomaxy,ny)
  meta%l(3) = gmm_layout(lowz,highz,haloz,halomaxz,nz)
  meta%a%initmode = initmode
  meta%a%flags = flags

  return
end subroutine gmm_build_meta3D

subroutine gmm_build_meta4D(meta,lowx,highx,halox,halomaxx,nx,lowy,highy,haloy,halomaxy,ny,&
                            lowz,highz,haloz,halomaxz,nz,low4,high4,halo4,halomax4,n4,initmode,flags)
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use rmn_gmm
  implicit none

  type(gmm_metadata), intent(out) :: meta
  integer :: lowx,highx,halox,halomaxx,nx
  integer :: lowy,highy,haloy,halomaxy,ny
  integer :: lowz,highz,haloz,halomaxz,nz
  integer :: low4,high4,halo4,halomax4,n4
  integer :: initmode,flags

  meta = GMM_NULL_METADATA
  meta%l(1) = gmm_layout(lowx,highx,halox,halomaxx,nx)
  meta%l(2) = gmm_layout(lowy,highy,haloy,halomaxy,ny)
  meta%l(3) = gmm_layout(lowz,highz,haloz,halomaxz,nz)
  meta%l(4) = gmm_layout(low4,high4,halo4,halomax4,n4)
  meta%a%initmode = initmode
  meta%a%flags = flags

  return
end subroutine gmm_build_meta4D

subroutine print_meta(meta)
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use rmn_gmm
  implicit none

  type(gmm_metadata), intent(in) :: meta

  write(6,'(a,5(2x,i6))') 'meta layout 1=',meta%l(1)
  write(6,'(a,5(2x,i6))') 'meta layout 2=',meta%l(2)
  write(6,'(a,5(2x,i6))') 'meta layout 3=',meta%l(3)
  write(6,'(a,5(2x,i6))') 'meta layout 4=',meta%l(4)
  write(6,'(a,i10)') 'meta key =',meta%a%key
  write(6,'(a,2(2x,i8))') 'meta uuid1,uuid2=',meta%a%uuid1,meta%a%uuid2
  write(6,'(a,2x,i6)') 'meta initmode =',meta%a%initmode
  write(6,'(a,2x,i8)') 'meta flags=',meta%a%flags
  return
end subroutine print_meta

subroutine gmm_print_keylist()
  use, intrinsic :: iso_fortran_env, only: REAL64, INT64
  use rmn_gmm
  implicit none

  character(len=GMM_MAXNAMELENGTH), dimension(:), allocatable :: keylist
  integer :: nkeys, i

  nkeys = gmm_nkeys()
  allocate(keylist(nkeys))
  nkeys = gmm_keys(keylist)
  write(6,'(a,2x,i3)') '(gmm_print_keylist) number of gmm keys is:',nkeys
  write(6,'(a)') '(gmm_print_keylist) key names are:'
  do i=1,nkeys
    write(6,'(a)') keylist(i)
  enddo
  deallocate(keylist, stat=i)
  return
end subroutine gmm_print_keylist
