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

module ens_spp
  implicit none
  private
  save

  ! External API
  public :: spp_options_init, spp_nml
  integer, parameter, public :: MAX_NSPP=50     !Maximum number of SPP instances (chains)
  integer, public :: spp_ncha, spp_lmax=-1, spp_mmax=-1, spp_latmax=-1
  character(len=64), dimension(:), pointer, public :: spp_list=>null()

  ! Internal parameters
  integer, parameter :: MAX_SPP_ELEMENTS=25     !Maximum dict/nml entries per SPP instance
  integer, parameter :: MAX_WB_ENTRIES=MAX_NSPP*MAX_SPP_ELEMENTS

contains

  function spp_options_init(F_spp_L, F_spplist) result(F_istat)
    logical :: F_spp_L                                          !Switch to use SPP
    character(len=*), dimension(:), intent(in) :: F_spplist     !List of SPP chain names
    integer :: F_istat
#include <rmnlib_basics.hf>
    integer :: i
    logical, save :: init_L = .false.
    F_istat = RMN_OK
    if (init_L) return
    init_L = .true.
    ! Determine the number of SPP chains required
    i = 1
    do while (F_spplist(i) /= '')
       i = i+1
       if (i > size(F_spplist)) exit
    enddo
    spp_ncha = i-1
    if (.not.F_spp_L) spp_ncha = 0
    if (.not.associated(spp_list)) then
       allocate(spp_list(spp_ncha))
       spp_list(:) = F_spplist(1:spp_ncha)
    endif
    return
  end function spp_options_init

  function spp_nml(F_fdout, F_inpath, F_cfgfile) result(F_istat)
    use wb_itf_mod
    use  clib_itf_mod, only: clib_toupper
    implicit none
    integer, intent(in) :: F_fdout
    character(len=*), intent(in) :: F_inpath, F_cfgfile
    integer :: F_istat
#include <rmnlib_basics.hf>
    integer :: i, j, k, stat, nlon, nlat, nkeys, dtype, tlen, dlen, opt
    integer, dimension(2) :: spp_trn
    integer, dimension(MAX_NSPP) :: ival
    real, dimension(MAX_NSPP) :: rval
    character(len=WB_MAXNAMELENGTH+1) :: prefix, section
    character(len=2048) :: dict_fname
    character(len=WB_MAXNAMELENGTH), dimension(MAX_WB_ENTRIES) :: keylist
    character(len=WB_MAXNAMELENGTH), dimension(MAX_NSPP) :: cval
    logical, dimension(MAX_NSPP) :: lval
    F_istat = RMN_ERR

    ! Only process configurations if SPP is active
    if (spp_ncha < 1) then
       F_istat = RMN_OK
       return
    endif

    ! Retrieve user configurations for all SPP entries
    dict_fname = trim(F_inpath)//'/spp.dict'
    do i=1,spp_ncha
       prefix = 'spp/'//trim(spp_list(i))//'/'
       HAVE_SPP_IN_WB: if (wb_check(prefix, WB_INITIALIZED) > 0) then
          if (WB_IS_ERROR(wb_get(trim(prefix)//'spp_nlat', nlat))) then
             write(F_fdout, *) 'Error retrieving stored lat definitions for '// &
                  trim(section)
             return
          endif
          if (WB_IS_ERROR(wb_get(trim(prefix)//'spp_trn', spp_trn, dlen))) then
             write(F_fdout, *) 'Error retrieving stored wavenumber definitions for '// &
                  trim(section)
             return
          endif
       else
          stat = wb_read(prefix, dict_fname, 'spp', WB_STRICT_DICTIONARY)
          if (WB_IS_ERROR(stat)) then
             write(F_fdout, *) 'Error reading spp entry from dictionary '// &
                  trim(dict_fname)//' into '//trim(spp_list(i))//': ',stat
             return
          endif
          section = 'spp_'//trim(spp_list(i))
          stat = wb_read(prefix, F_cfgfile, section, WB_FORBID_DEFINE)
          if (WB_IS_ERROR(stat)) then
             write(F_fdout, *) 'Error reading '//trim(section)// &
                  ' configuration from '//trim(F_cfgfile)
             return
          endif
          stat = WB_OK
          stat = min(wb_get(trim(prefix)//'spp_nlat', nlat), stat)
          nlon = 2*nlat
          stat = min(wb_put(trim(prefix)//'spp_nlon', nlon), stat)
          if (.not.WB_IS_OK(stat)) then
             write(F_fdout, *) 'Error processing SPP lat/lon definitions for '// &
                  trim(section)//' in '//trim(F_cfgfile)
             return
          endif
          if (WB_IS_ERROR(wb_get(trim(prefix)//'spp_trn', spp_trn, dlen))) then
             write(F_fdout, *) 'Error retrieving SPP wavenumber definitions for '// &
                  trim(section)//' in '//trim(F_cfgfile)
             return
          endif
       endif HAVE_SPP_IN_WB
       spp_lmax = max( maxval(spp_trn) - minval(spp_trn) + 1, spp_lmax )
       spp_mmax = max( maxval(spp_trn) + 1, spp_mmax )
       spp_latmax=max(nlat,spp_latmax)
    enddo

    ! Post required SPP information to whiteboard
    if (WB_IS_ERROR(wb_put('spp/list', spp_list(:), WB_REWRITE_MANY))) then
       write(F_fdout, *) 'Error posting SPP list to whiteboard'
       return
    endif

    ! Post-treatment of values and generation of listing on request
    if (WB_IS_ERROR(wb_keys(keylist, nkeys, 'spp'))) then
       write(F_fdout, *) 'Error retrieving spp entries in whiteboard'
       return
    endif
    stat = WB_OK
    do i=1,nkeys
       if (WB_IS_ERROR(wb_get_meta(keylist(i), dtype, tlen, dlen, opt))) then
          write(F_fdout, *) 'Error retrieving metadat for '//trim(keylist(i))
          return
       endif
       if (dlen > MAX_NSPP) then
          write(F_fdout, *) 'Too many entries for '//trim(keylist(i))// &
               '.  Increase MAX_NSPP in spp::spp_nml to at least ',dlen
          return
       endif
       stat = WB_OK
       select case (dtype)
       case (WB_FORTRAN_REAL)
          if (dlen == 0) then
             stat = min(wb_get(keylist(i), rval(1)), stat)
          else
             stat = min(wb_get(keylist(i), rval, dlen), stat)
          endif
          if (F_fdout >= 0) write(F_fdout, *) keylist(i)//' = ',(rval(j), j=1,max(dlen,1))
       case (WB_FORTRAN_INT)
          if (dlen == 0) then
             stat = min(wb_get(keylist(i), ival(1)), stat)
          else
             stat = min(wb_get(keylist(i), ival, dlen), stat)
          endif
          if (F_fdout >= 0) write(F_fdout, *) keylist(i)//' = ',(ival(j), j=1,max(dlen,1))
       case (WB_FORTRAN_CHAR)
          if (dlen == 0) then
             stat = min(wb_get(keylist(i), cval(1)), stat)
             j = clib_toupper(cval(1))
             stat = min(wb_put(keylist(i), cval(1)), stat)
          else
             stat = min(wb_get(keylist(i), cval, dlen), stat)
             do k=1,dlen
                j = clib_toupper(cval(k))
             enddo
             stat = min(wb_put(keylist(i), cval(1:dlen)), stat)
          endif
          if (F_fdout >= 0) write(F_fdout, *) keylist(i)//' = ',(trim(cval(j))//',', j=1,max(dlen,1))
       case (WB_FORTRAN_BOOL)
          if (dlen == 0) then
             stat = min(wb_get(keylist(i), lval(1)), stat)
          else
             stat = min(wb_get(keylist(i), lval, dlen), stat)
          endif
          if (F_fdout >= 0) write(F_fdout, *) keylist(i)//' = ',(lval(j), j=1,max(dlen,1))
       case DEFAULT
          continue
       end select
    enddo

    F_istat = RMN_OK
  end function spp_nml

end module ens_spp
