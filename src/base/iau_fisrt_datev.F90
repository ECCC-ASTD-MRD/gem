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
      subroutine iau_fisrt_datev (F_datev, F_ubstp, F_stepno)
      use iso_c_binding
      use rmn_fst24
      use inp_mod
      use step_options
      use init_options
      use cstv
      use path
      use ptopo
      implicit none
#include <rmnlib_basics.hf>

      character(len=*), intent(OUT) :: F_datev
      integer         , intent(OUT) :: F_ubstp
      integer         , intent(IN ) :: F_stepno

      include 'mpif.h'
      character(len=2048) root,list_files,cmd,filename,fn
      character(len=16) :: pdfdate
      logical :: success,found_L
      integer :: err,err_code,unf,n,valid_lbound, valid_ubound
      real(kind=REAL64) :: dayfrac
      real(kind=REAL64), parameter :: one=1.0d0, &
                                      sid=86400.0d0, rsid=one/sid
      type(fst_file) :: fst
      type(fst_record) :: my_record
      type(fst_query)  :: my_query
!     
!     ---------------------------------------------------------------
!
! By giving Iau_interval, the user is responsible for providing
! increments files valid for a certain number of timesteps around
! the validity time date0 of each files without overlap between
! files. The users will also provide enough of those files to cover
! the entire Iau_period. Hence the computing of a range of validdity
! timesteps below (valid_lbound,valid_ubound):
      
      found_L= .false. ; err_code= 0 ; F_datev= 'NUL'
      root= trim(Path_input_S)//'/IAUREP'
      list_files = 'IAUREP_files'

      if (Inp_iome >= 0) then
         cmd='ls -1 '//trim(root)//' > '//trim(list_files)
         call system(cmd)

         unf= 0 ; err_code= -1
         if (fnom( unf,trim(list_files),'SEQ+FMT+OLD',0 ) /= 0) unf= 0
         if (unf == 0) goto 33
 55      read (unf,*,end=33) filename
         fn= trim(root)//'/'//trim(filename)

         if (fst%open(trim(fn),'RND+OLD+R/O')) then
            my_query = fst%new_query(nomvar = 'TT')
            success  = my_query%find_next(my_record)
            if (success) then
               ! Determine range of validity timestep
               call datf2p ( pdfdate,my_record%datev )
               call difdatsd(dayfrac,Step_runstrt_S,pdfdate)
               n= Iau_interval/2.d0/Cstv_dt_8
               valid_ubound= nint(dayfrac*sid/Cstv_dt_8) + n
               valid_lbound= nint(dayfrac*sid/Cstv_dt_8) - n + 1
               if ((F_stepno>=valid_lbound).and.&
                   (F_stepno<=valid_ubound)) then
                   F_datev = pdfdate
                   F_ubstp = valid_ubound
                   found_L = .true.
                endif
             endif
             call my_query%free()
             success= fst%close()
          endif
          goto 55
          
 33       if (found_L) err_code= 0
       endif

       call gem_error ( err_code, 'iau_fisrt_datev', &
                       'Unable to find valid data in increments files' )
       call MPI_bcast (F_ubstp,1,MPI_INTEGER, 0,&
                       COMM_multigrid, err)
       call MPI_bcast (F_datev,len(F_datev),MPI_CHARACTER, 0,&
                       COMM_multigrid, err)
!     ---------------------------------------------------------------
!
      return
      end subroutine iau_fisrt_datev

