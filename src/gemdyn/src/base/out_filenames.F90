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

!**s/r out_filenames - liste of output files to open

      subroutine out_filenames ( F_sorties, F_levels, F_grid, F_prefix_S, F_usrdir, &
                                 d1,d2, F_nfiles, F_filenames_S)
      use timestr_mod
      use dimout
      use levels
      use step_options
      use gem_options
      use out_options
      use cstv
      use out_mod
      use out3
      use outgrid, only: OutGrid_hgrid_usr
      use outusrdir, only: OutUsrdir_name_S
      use, intrinsic :: iso_fortran_env
      implicit none

      character(len=*)       , intent (IN)  :: F_prefix_S
      character(len=*)       , intent (OUT) :: F_filenames_S(*)
      integer                , intent (IN)  :: d1,d2
      integer                , intent (INOUT)  :: F_nfiles
      integer, dimension(0:MAXSET,d1:d2), intent (IN) :: F_sorties
      integer, dimension(MAXSET)        , intent (IN) :: F_levels, F_grid, F_usrdir

      character(len=2 ) :: prefix_S, dumc
      character(len=4 ) :: unit_ext
      character(len=16) :: datev,fdate
      character(len=32) :: filen,myformat_S,my_hour
      logical keep
      integer prognum,err,i,jj,kk,indx,len0,len1,levset,gridset,usrdir
      real(kind=REAL64), parameter :: OV_day = 1.0d0/86400.0d0
      real(kind=REAL64)  dayfrac
!
!------------------------------------------------------------------
!
      do jj=1, F_sorties(0,Lctl_step)

            if ( trim(F_prefix_S) == 'casc' ) then

               dayfrac = dble(lctl_step-Step_delay) * Cstv_dt_8 * OV_day
               call incdatsd (datev,Step_runstrt_S,dayfrac)
               filen = trim(F_prefix_S)//'_'//trim(datev)

            else

               kk       = F_sorties(jj,Lctl_step)
               gridset  = f_grid(kk)
               levset   = F_levels(kk)
               usrdir   = F_usrdir(kk)
               if(usrdir /= -1)then
                  prefix_S(1:2)='u'//OutUsrdir_name_S(usrdir)(4:4)
               else
                  prefix_S(1:2) = F_prefix_S//Level_typ_S(levset)
               endif               
               call up2low (prefix_S ,dumc)
               prefix_S= dumc
               err= timestr_prognum (prognum, Out3_unit_S, Out3_close_interval,&
                          Out_dateo,sngl(Cstv_dt_8),lctl_step,Out_endstepno)
               unit_ext = ' '
               if (Out3_unit_S(1:3) == 'SEC') unit_ext = 's'
               if (Out3_unit_S(1:3) == 'MIN') unit_ext = 'm'
               if (Out3_unit_S(1:3) == 'DAY') unit_ext = 'd'
               if (Out3_unit_S(1:3) == 'STE') unit_ext = 'p'
               if (Out3_unit_S(1:3) == 'MON') unit_ext = 'n'

               len1 = max(3,Out3_ndigits)
               if (any(Out3_unit_s(1:3) == (/'SEC','STE'/))) len1 = max(6,len1)
               write(my_hour,'(i20)')  abs(prognum)
               indx=1
 77            i=index(my_hour(indx:20)," ")
               if ((i > 0).and.(indx <= 20)) then
                  indx=indx+1
                  goto 77
               end if
               len1 = max(len1,max(0,21-indx))
               len0 = len1
               if (prognum < 0) len0 = len0+1
               write(myformat_S,'(a,i1.1,a,i1.1,a)') '(a,i',len0,'.',len1,')'
               my_hour = ' '
               write(my_hour,trim(myformat_S)) '_', prognum

               call datf2p (fdate, Out3_date)
               filen= trim(prefix_S)//fdate(1:8)//fdate(10:11)// &
                           trim(my_hour)//trim(unit_ext)

            end if

            keep=.true.
            do i=1,F_nfiles
               if (filen == F_filenames_S(i)) keep=.false.
            end do
            if (keep) then
               F_nfiles= F_nfiles+1
               F_filenames_S(F_nfiles)= trim(filen)
            endif
         enddo
!
!------------------------------------------------------------------
!
      return
      end subroutine out_filenames

