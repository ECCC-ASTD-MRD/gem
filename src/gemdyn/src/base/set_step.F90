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

!**s/r set_step - initialization of common block TIMESTEP
!
      integer function set_step(F_argc, F_argv_S, F_cmdtyp, F_v1, F_v2)
      use ctrl
      use cstv
      use lun
      use step_options
      use timestep
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

        integer, intent(in) :: F_argc,F_v1,F_v2
        character(len=*), intent(in) :: F_argv_S(0:F_argc),F_cmdtyp
!
!object
!       initialization of the common block TIMESTEP. This function is
!       called when the keyword "steps" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
! ie:   steps=1,hour,0.,3.,6.,12.,24.,48.;
!
!       The "rpn_fortran_callback" routine will process the above
!       statement and return 5 arguments to this function. For more
!       information to how this is processed, see "SREQUET".
!
!arguments
!  Name                        Description
!------------------------------------------------------------
! F_argc        - number of elements in F_argv_S
! F_argv_S      - array of elements received
!                 if F_argv_S(ii) contains "[", the value in this
!                  argument indicates number of elements following it.
! F_cmdtyp      - character command type - not used
! F_v1          - integer parameter 1 - not used
! F_v2          - integer parameter 2 - not used
!
!Notes:
!
! steps=stepset#,[init],step/hour,{list};
!
! ie:   steps=2,step,-1;
!       steps=3,hour,<0.,48.,3.>;
!       steps=4,init,step,<0.,6.,1.>;
!       steps=5,step,[1,2,3,9,11];
!       steps=6,day,<0.,31.,1.>;
!       steps=7,month,<0.,1.,1.>;
!
!      Should label the stepset# sequentially: 1,2,3,....
!      'init' - means this command only applies to output during the
!               initialization period
!      'month'- output in months
!      'day'  - output in days
!      'hour' - output in hours
!      'step' - output in timesteps in the model
!      '-1' with "step" will give every timestep of the model
!      [a,b,c] means a,b and c are requested
!      <a,b,c> means a to b, incrementing every c are requested
!
! Concerning the 'month' option: all months end at 24h00 on the last
! day of a month so that, depending on the starting date, the first
! month could be shorter than expected.
!
      character(len=5) :: stuff
      integer i, j, k, ii, num, istep
      logical month_flag, day_flag, hour_flag, step_flag, found_L
      real frarg, month, day, hour, rstep
      integer transtep, transtep2, stepset, argc_out
      transtep(frarg) = nint(3600.0 * frarg / Cstv_dt_8)
      transtep2(frarg) = nint(86400.0 * frarg / Cstv_dt_8)
!
      integer,external :: month_since_start
      integer,external :: dcmip_div_X
!
      argc_out=min(F_argc,6)
      if (Lun_out > 0) then
          write(Lun_out,*)
          if (argc_out < F_argc) then
            write(Lun_out,*) F_argv_S(0),'=',F_argv_S(1),',',F_argv_S(2),',',(F_argv_S(i),i=3,argc_out),'...'
          else
            write(Lun_out,*) F_argv_S(0),'=',F_argv_S(1),',',F_argv_S(2),',',(F_argv_S(i),i=3,argc_out)
          end if
      end if
      set_step=0
      read(F_argv_S(1),*)stepset
      Timestep_sets = Timestep_sets + 1

      if (Timestep_sets > MAXSET) then
          if (Lun_out > 0) then
            write(Lun_out,*)'SET_STEP WARNING: too many TIMESTEP sets'
          end if
          Timestep_sets = Timestep_sets - 1
          set_step=1
          return
      end if

      j=Timestep_sets
      i=0
      month_flag= .false.
      day_flag  = .false.
      hour_flag = .false.
      step_flag = .false.
      Timestep_id(j)=stepset
      Timestep_init_L(j)=.false.
      do ii=2,F_argc
         if (index(F_argv_S(ii),'[') > 0) then
             stuff=F_argv_S(ii)
             read(stuff(2:4),*) num
         else if (F_argv_S(ii) == 'month' ) then
             month_flag= .true.
             day_flag  = .false.
             hour_flag = .false.
             step_flag = .false.
         else if (F_argv_S(ii) == 'day' ) then
             month_flag= .false.
             day_flag  = .true.
             hour_flag = .false.
             step_flag = .false.
         else if (F_argv_S(ii) == 'hour') then
             hour_flag = .true.
             step_flag = .false.
             day_flag  = .false.
             month_flag= .false.
         else if (F_argv_S(ii) == 'step') then
             step_flag = .true.
             hour_flag = .false.
             day_flag  = .false.
             month_flag= .false.
         else if (F_argv_S(ii) == 'init') then
             Timestep_init_L(j)=.true.
         else if (step_flag) then
             i = i+1
             read(F_argv_S(ii),*)rstep
             if (rstep == -1) then
             i = i-1
             do istep=lctl_step,Step_total
               i = i+1
               Timestep_tbl(i,j)=istep
             end do
             else
               Timestep_tbl(i,j)=int(rstep)
             end if
         else if (hour_flag) then
             i = i+1
             read(F_argv_S(ii),*)hour
             Timestep_tbl(i,j)=transtep(hour)
             if ( Ctrl_canonical_dcmip_L ) then
               Timestep_tbl(i,j)= dcmip_div_X(transtep(hour))
             end if
         else if (day_flag) then
             i = i+1
             read(F_argv_S(ii),*)day
             Timestep_tbl(i,j)=transtep2(day)
         else if (month_flag) then
             i = i+1
             read(F_argv_S(ii),*)month
             Timestep_tbl(i,j)=month_since_start(month)
         else
             if (Lun_out > 0) then
               write(Lun_out,*)'SET_STEP WARNING: Timestep type not recognizable'
             end if
             Timestep_sets = Timestep_sets - 1
             set_step=1
             return
         end if
      end do

      if (i > MAXSTEP) then
          if (Lun_out > 0) then
            write(Lun_out,*)'SET_STEP WARNING: Requested timesteps > MAXSTEP'
          end if
          Timestep_sets = Timestep_sets - 1
          set_step=1
          return
      end if
!
!     Eliminate repeated timesteps in one Timestep set
      istep = 1
      do ii = 2, i
         found_L = .false.
         do k = 1, ii-1
            if ( Timestep_tbl(ii,j) == Timestep_tbl(k,j) ) found_L = .true.
         end do
         if (.not. found_L) then
             istep = istep + 1
             Timestep_tbl(istep,j) = Timestep_tbl(ii,j)
         end if
      end do

      Timestep_max(Timestep_sets)=istep

      if (Lun_out > 0) then
         write(Lun_out,*) ' Timestep_set(',j,') : Timestep_id=',Timestep_id(j)
         write(Lun_out,*) ' Timestep_init_L=',Timestep_init_L(j)
         if (Timestep_max(j) > 30) then
            write(Lun_out,*) ' Timestep=', &
               (Timestep_tbl(i,j),i=1,30),',... up to ,',Timestep_tbl(Timestep_max(j),j)
         else
            write(Lun_out,*) ' Timestep=',(Timestep_tbl(i,j),i=1,Timestep_max(j))
         end if
      end if

      return
      end

      integer function month_since_start (month)
      use glb_ld
      use cstv
      use lun
      use dimout
      use out_mod
      use out3
      use, intrinsic :: iso_fortran_env
      implicit none

#include <arch_specific.hf>
      real     month

!arguments
! month        I    - months since startup (timestep 0)
!
! Note that all months end at 24h00 on the last day of a month
! so that, depending on the starting date, the first month could
! be shorter than expected.
!
      real(kind=REAL64)   hours
      logical, save :: First_call_L = .true.
      integer, save :: Dts_month1(14),Steps_month1
      integer  Start_date(14),Current_month(14)
      integer  ier, MOD2,MOD3, MOD

      integer  newdate
      external newdate,datmgp2

!----------------------------------------------------------------
      month_since_start = -1

      if(first_call_L) then

         First_call_L   = .false.
         Start_date(14) = Out3_date

!        decode Out3_date
         call datmgp2( Start_date )

!        define the end of the first month
         Dts_month1(2) = Start_date(2)+1  ! month (1-12)
         Dts_month1(3) = 1                ! day   (1-31)
         Dts_month1(4) = Start_date(4)    ! year
         Dts_month1(5) = 0                ! zulu  (0-23)
         Dts_month1(6) = 0                ! 100*seconds since last hour
         if(Dts_month1(2) == 13) then
            Dts_month1(2) = 1
            Dts_month1(4) = Dts_month1(4)+1
         end if

!        encode this information
         MOD2 = (Dts_month1(4)*100+Dts_month1(2))*100+1
         MOD3 = 0 ; MOD = +3
         ier = newdate( Dts_month1(14), MOD2,MOD3, MOD )
         if (ier /= 0) then
            if (Lun_out > 0) write(Lun_out,*)'month_since_start: error un newdate(+3), MOD2,MOD3=',MOD2,MOD3
            stop ' in month_since_start, called by set_step'
         end if

!        number of hours between that and the start
         call difdatr( Dts_month1(14),Start_date(14),hours )
         Steps_month1 = nint(3600.0 * hours / Cstv_dt_8)

      end if

      if(month == 0) then

         month_since_start = 0

      else if(month == 1) then

         month_since_start = Steps_month1

      else if(month > 1) then

!        define the current month
         Current_month(2) = Dts_month1(2)+nint( amod( month-1, 12. ) )
         Current_month(3) = 1
         Current_month(4) = Dts_month1(4)+          ( month-1)/12
         Current_month(5) = 0
         Current_month(6) = 0
         if(Current_month(2) > 12) then
            Current_month(2) = Current_month(2)-12
            Current_month(4) = Current_month(4)+1
         end if

!        encode the current month information
         MOD2 = (Current_month(4)*100+Current_month(2))*100+1
         MOD3 = 0 ; MOD = +3
         ier = newdate( Current_month(14), MOD2,MOD3, MOD )
         if (ier /= 0) then
            if (Lun_out > 0) write(Lun_out,*)'month_since_start: error un newdate(+3), MOD2,MOD3=',MOD2,MOD3
            stop ' in month_since_start, called by set_step'
         end if

!        number of hours between that and the start
         call difdatr( Current_month(14),Dts_month1(14),hours )
         month_since_start = Steps_month1 + nint(3600.0 * hours / Cstv_dt_8)

      end if

      return
      end


