!---------------------------------- LICENCE BEGIN ------------------------------
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
!---------------------------------- LICENCE END --------------------------------

module series_nml_mod
   use clib_itf_mod, only: clib_isreadok, clib_toupper
   use wb_itf_mod
   use series_options
   private
   public :: series_nml

!!!#include <arch_specific.hf>
#include <rmn/msg.h>
#include <rmnlib_basics.hf>

   integer, parameter, public :: SER_NML_ERR  = RMN_ERR
   integer, parameter, public :: SER_NML_NONE = RMN_OK
   integer, parameter, public :: SER_NML_OK   = RMN_OK + 1

   !@author R. McTaggart-Cowan    - Apr 2009
   !@revisions
   ! v4_50 - Desgagne M.    - rename the routine & other minors mods
   ! 2017-01, S.Chamberland - Major revision

contains

   !/@*
   function series_nml(F_namelistf_S) result(F_istat)
      implicit none
      !@object: Default configuration and reading series namelist
      !@params
      character(len=*), intent(in) :: F_namelistf_S
      !@return
      integer :: F_istat
      !*@/

      integer :: istat
      !----------------------------------------------------------------
      F_istat = SER_NML_ERR

      istat = series_options_init()
      if (.not.RMN_IS_OK(istat)) return
      istat = ser_nml_read(F_namelistf_S)
      if (.not.RMN_IS_OK(istat)) return
      istat = ser_nml_check()
      if (.not.RMN_IS_OK(istat)) return
      istat = ser_nml_post_init()
      if (.not.RMN_IS_OK(istat)) return
      call ser_nml_print()

      F_istat = SER_NML_OK
      !----------------------------------------------------------------
      return
   end function series_nml

   !#---- Private functions -------------------------------------------

   function ser_nml_read(m_namelist) result(m_istat)
      implicit none
      character(len=*), intent(in) :: m_namelist
      integer :: m_istat
      character(len=1024) :: msg_S, namelist_S
      integer :: istat, unf
      !----------------------------------------------------------------
      m_istat = RMN_ERR
      namelist_S = '&series'

      unf = 0
      istat = clib_isreadok(m_namelist)
      if (RMN_IS_OK(istat)) istat = fnom(unf, m_namelist, 'SEQ+OLD', 0)
      if (.not.RMN_IS_OK(istat)) then
         write(msg_S,'(a,a,a,a)') PKGNAME_S, &
              'Using default config, ', &
              'Namelist file Not found/readable: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
         return
      endif

      read(unf, nml=series, iostat=istat)
      if (istat == 0) then     !# Read ok
         m_istat = RMN_OK + 1
      else if (istat < 0) then !# EOF, nml not found
         write(msg_S,'(a,a,a,a,a)') PKGNAME_S, &
              'No Series, Namelist ',&
              trim(namelist_S), ' not available in file: ', trim(m_namelist)
         call msg(MSG_INFO, msg_S)
         m_istat = RMN_OK
      else !# Error
         write(msg_S,'(a,a,a,a,a)') PKGNAME_S, &
              'Namelist', trim(namelist_S), &
              ' invalid in file: ', trim(m_namelist)
         call msg(MSG_ERROR, msg_S)
      endif
      istat = fclos(unf)
      !----------------------------------------------------------------
      return
   end function ser_nml_read


   subroutine ser_nml_print()
      implicit none
      integer, external :: msg_getUnit
      integer :: fileid, i
      !----------------------------------------------------------------
      fileid = msg_getUnit(MSG_INFO)
      if (fileid > 0 .and. series_on_L) then
!!$         write(fileid, nml=series)
         write(fileid,'(a)') "&series"
         write(fileid,'(a,i6)') "P_serg_serstp = ",P_serg_serstp
         write(fileid,'(a,i6)') "P_serg_srwri  = ",P_serg_srwri
         do i=1,size(P_serg_srprf_s)
            if (trim(P_serg_srprf_s(i)) /= ' ') then
               write (fileid,'(a,i3,a,a)') "P_serg_srprf_s(",i,") = ", &
                    trim(P_serg_srprf_s(i))
            else 
               exit
            endif
         enddo
         do i=1,size(P_serg_srsrf_s)
            if (trim(P_serg_srsrf_s(i)) /= ' ') then
               write(fileid,'(a,i3,a,a)') "P_serg_srsrf_s(",i,") = ", &
                    trim(P_serg_srsrf_s(i))
            else 
               exit
            endif
         enddo
         write(fileid,'(a)') "/"
      endif
      !----------------------------------------------------------------
      return
   end subroutine ser_nml_print


   function ser_nml_check() result(m_istat)
      implicit none
      integer :: m_istat
      !----------------------------------------------------------------
      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function ser_nml_check


   function ser_nml_post_init() result(m_istat)
      implicit none
      integer :: m_istat

      integer, external :: msg_getUnit
      integer :: istat, options
      character(len=1024) :: msg_S
      !----------------------------------------------------------------
      m_istat = RMN_ERR

      !# Determine if the user has specified any stations
      series_nstng = 0
      do while(series_nstng < NSTATMAX .and. &
           &   int(xst_stn_latlon(series_nstng+1)%lat) /= STN_MISSING .and. &
           &   int(xst_stn_latlon(series_nstng+1)%lon) /= STN_MISSING)
         series_nstng = series_nstng + 1 
         series_stng(series_nstng)%name = xst_stn_latlon(series_nstng)%name
         series_stng(series_nstng)%lat  = xst_stn_latlon(series_nstng)%lat
         series_stng(series_nstng)%lon  = xst_stn_latlon(series_nstng)%lon
         istat = clib_toupper(series_stng(series_nstng)%name)
      enddo

      series_on_L = .false.
      IF_NSTAT: if (series_nstng > 0) then
         write(msg_S,'(i5,a)') series_nstng, ' Stations specified'
         call msg(MSG_INFO, PKGNAME_S//msg_S)

         series_nsurf = 0
         do while(series_nsurf < NVARMAX .and. &
              P_serg_srsrf_s(series_nsurf+1) /= ' ') 
            istat = clib_toupper(P_serg_srsrf_s(series_nsurf+1))
            series_nsurf = series_nsurf+1
         enddo

         series_nprof = 0
         do while(series_nprof < NVARMAX .and. &
              P_serg_srprf_s(series_nprof+1) /= ' ')
            istat = clib_toupper(P_serg_srprf_s(series_nprof+1))
            series_nprof = series_nprof+1
         enddo

         series_on_L = (series_nsurf + series_nprof > 0)
      endif IF_NSTAT
      if (.not.series_on_L) then
         series_nstng = 0
         call msg(MSG_INFO, PKGNAME_S//'NO TIME SERIES REQUESTED')
      endif

      options = WB_REWRITE_NONE + WB_IS_LOCAL
      istat = wb_put('model/series/P_serg_srsus_L', series_on_L, options)

      m_istat = RMN_OK
      !----------------------------------------------------------------
      return
   end function ser_nml_post_init

end module series_nml_mod
