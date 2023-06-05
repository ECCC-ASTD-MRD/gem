!/****************************************************************************/
!/****************       Diurnal FFMC Adjustment Program      ****************/
!/*                                                                          */
!/* Date: March 2, 1993                                                      */
!/*                                                                          */
!/* Authors: Brad Armitage                       Bruce Lawson                */
!/*          Fire Research Technician            Head, Fire Research Program */
!/*          Pacific Forestry Centre             Pacific Forestry Centre     */
!/*          506 W. Burnside Rd.                 506 W. Burnside Rd.         */
!/*          Victoria, B.C.                      Victoria, B.C.              */
!/*          Phone (604) 363-0693                Phone (604) 363-0710        */
!/*          E-mail barmitage@pfc.forestry.ca                                */
!/*              or   blawson@pfc.forestry.ca                                */ 
!/*                                                                          */
!/* Slightly modified by Kerry Anderson, programmer extraorinaire:           */
!/*                                                                          */
!/* 1. Calc_New_FFMC accepts a negative rh as missing (90 is substituted in) */
!/* 2. routines main, InputError, Input3_Output and Input4_Output            */
!/*    have been removed                                                     */
!/****************************************************************************/
 module mach_cffeps_diurnal_mod

   private
   public :: bdl_ffmc, eq_ffmc, hourly_ffmc

   contains

!/*----------------------------------------------------------------------------*/
   real(kind=4) function bdl_ffmc(ff_ffmc_in, hour_in, rh) result(adj_ffmc)
!/*----------------------------------------------------------------------------*
! This function first checks that the program inputs are within the
! appropriate ranges. If any of the input values are out of
! range an error message is printed to the screen and control is
! returned to the operating system.
! NOTE: Appropriate ranges for inputs are defined as follows:
!        hour    ( >= 1    and <= 2459  ) (integer)
!        RH      ( >= 0    and <= 100   ) (integer)
!    FF_FFMC ( >= 17.5 and <= 100.9 ) (one decimal place)
! Where: FF_FFMC = 17.5  corresponds to the lower limit of Van Wagner's
!            original diurnal adjustment graphs;
!    FF_FFMC = 100.9 corresponds to the theoretical upper limit of
!            the current FF_FFMC scale.
!
! The function then calculates the moisture content at
! 1600 (mc1600) based on the std. FF - FFMC value.
! This moisture content is then used in the appropriate equation
! to predict the moisture content at the desired time (adj_mc)
! between 1200 noon on the current day and 1159 the next day.
! Calculations for the various hours are based on the following:
!  1200 to 2000 : Equations for every hour (interpolated for minute resolution).
!  2001 to  559 : Interpolation between 2000 and 600 (using high RH equation)
!  600  to 1100 : Equations for every hour (interpolated for minute resolution)
!  1101 to 1159 : Extrapolation using 1100 and 1000 as end-points.
!
! The function then calculates the new FFMC from the new calculated
! moisture content and returns this time adjusted FFMC value.
!*-----------------------------------------------------------------------------*/
     implicit none

     integer(kind=4), intent (in) :: hour_in
     real(kind=4),    intent (in) :: ff_ffmc_in, rh

!/*------------- Declare Local Variables ------------------*/
     real(kind=4)    :: mc1600, adj_mc, ff_ffmc
     real(kind=4)    :: next_hr, current_hr, minutes, min_since_2000
     real(kind=4)    :: tran_crve, tran_dist, tran_end, tran_str
     integer(kind=4) :: hour, hrs_since_2000
     real(kind=4), parameter :: uemcb = 30.0    !/* Lower Equilibrium Moisture Content Bound */
     real(kind=4), parameter :: lemcb = 5.0     !/* Upper Equilibrium Moisture Content Bound */

     adj_mc   = 0.0
     adj_ffmc = -1.0
!/*------------------- Check validity of input data ---------------------------*/
      if (hour_in < 0 .or. hour_in > 4859 .or. rh < 0.0 .or. rh > 100.0 .or. &
          ff_ffmc_in < 0.0 .or. ff_ffmc_in >= 101.0) then
         return
      end if

      !/*--- Check Time
      hour = hour_in
      if (hour_in > 2459) hour = hour_in - 2400
      minutes = real(mod(hour, 100))

      if (minutes > 59.0) return

      !/*--- Check for low end of FFMC scale ---*/
      ff_ffmc = max(ff_ffmc_in, 17.5)

      !/*---------------- Calculate moisture content at 1600 ----------------*/
      mc1600 = 147.2 * (101.0 - ff_ffmc) / (59.5 + ff_ffmc)

      !/*--------------- Select the appropriate set of equations ------------*/
      if (hour >= 600 .and. hour < 700) then
         if (rh > 87.0) then
            current_hr = g4600hi(mc1600)
            next_hr    = g4700hi(mc1600)
            adj_mc     = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else if (rh >= 68.0 .and. rh <= 87.0) then
            if (rh > 77.0) then
               if (hour <= 630) then
                  adj_mc = g4600md(mc1600)
               else 
                  adj_mc = g4700hi(mc1600)
               end if
            else if (rh >= 58.0 .and. rh <= 77.0) then
               current_hr    = g4600md(mc1600)
               next_hr       = g4700md(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            end if
         else !if (rh < 68.0) then
            if (rh >= 58.0 .and. rh <= 77.0) then
               if (hour <= 630) then
                  adj_mc = g4600lo(mc1600)
               else
                  adj_mc = g4700md(mc1600)
               end if
            else
               current_hr = g4600lo(mc1600)
               next_hr    = g4700lo(mc1600)
               adj_mc     = current_hr + &
                            interp(current_hr, next_hr, minutes, 60.)
            end if
         end if
      else if (hour >= 700 .and. hour < 800) then
         if (rh > 77.0) then
            current_hr    = g4700hi(mc1600)
            next_hr       = g4800hi(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else if (rh >= 58.0 .and. rh <= 77.0) then
            if (rh > 67.0) then
               if (hour <= 730) then
                  adj_mc = g4700md(mc1600)
               else
                  adj_mc = g4800hi(mc1600)
               end if
            else !if (rh >= 48.0 && rh <= 67.0) then
              current_hr = g4700md(mc1600)
              next_hr    = g4800md(mc1600)
              adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
            end if
         else if (rh < 58.0) then
            if (rh >= 48.0) then !if (rh >= 48.0 .and. rh <= 67.0) then
               if (hour <= 730) then
                  adj_mc = g4700lo(mc1600)
               else
                 adj_mc = g4800md(mc1600)
               end if
            else if (rh < 48.0) then
               current_hr    = g4700lo(mc1600)
               next_hr       = g4800lo(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            end if
         end if

      else if (hour >= 800 .and. hour < 900) then
         if (rh > 67.0) then
            current_hr    = g4800hi(mc1600)
            next_hr       = g4900hi(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else if (rh >= 48.0) then !(rh >= 48.0 .and. rh <= 67.0) then
            if (rh > 62.0) then
               if (hour <= 830) then
                  adj_mc = g4800md(mc1600)
               else
                  adj_mc = g4900hi(mc1600)
               end if
            else if (rh >= 43.0) then !if (rh >= 43.0 .and. rh <= 62.0) then
               current_hr    = g4800md(mc1600)
               next_hr       = g4900md(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            end if
         else if (rh < 48.0) then
            if (rh >= 43.0) then
               if (hour <= 830) then
                  adj_mc = g4800lo(mc1600)
               else
                  adj_mc = g4900md(mc1600)
               end if
            else if (rh < 43.0) then
               current_hr    = g4800lo(mc1600)
               next_hr       = g4900lo(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            end if
         end if

      else if (hour >= 900 .and. hour < 1000) then
         if (rh > 62.0) then
            current_hr    = g4900hi(mc1600)
            next_hr       = g41000hi(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else if (rh >= 43.0) then ! && rh <= 62.0)
            if (rh > 57.0) then
               if (hour <= 930) then
                  adj_mc = g4900md(mc1600)
               else             
                  adj_mc = g41000hi(mc1600)
               end if
            else !if (rh >= 38.0 && rh <= 57.0)
               current_hr    = g4900md(mc1600)
               next_hr       = g41000md(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            end if
         else if (rh < 43.0) then
            if (rh >= 38.0) then ! && rh <= 57.0)
               if (hour <= 930) then
                  adj_mc = g4900lo(mc1600)
               else             
                  adj_mc = g41000md(mc1600)
               end if
            else !if (rh < 38.0) then
               current_hr    = g4900lo(mc1600)
               next_hr       = g41000lo(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            end if
         end if  

      else if (hour >= 1000 .and. hour < 1100) then
         if (mc1600 <= 5.0) then                     !/*  1000**() */
            if (rh > 57.0) then
               adj_mc = g41000hi(mc1600)
            else if (rh >= 38.0 .and. rh <= 57.0) then
               adj_mc = g41000md(mc1600)
            else !if (rh < 38.0) then
               adj_mc = g41000lo(mc1600)
            end if
         else if (mc1600 > 5.0 .and. mc1600 <= 30.0) then
            if (rh > 57.0) then                       !/* 1000hi() && 1100hi() */
               current_hr = g41000hi(mc1600)
               tran_str   = g41000hi(5.0)
               tran_end   = g31100hi(30.0)
               tran_dist  = ((mc1600 - lemcb))
               tran_crve  = tran_str + &
                            interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
               next_hr    = tran_crve
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
            else if (rh >= 38.0) then ! .and. rh <= 57.0) then
               if (rh > 54.5) then                     !/* 1000md() && 1100hi() */
                   current_hr = g41000md(mc1600)
                   tran_str   = g41000md(5.0)
                   tran_end   = g31100hi(30.0)
                   tran_dist  = ((mc1600 - lemcb))
                   tran_crve  = tran_str + &
                            interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                   next_hr    = tran_crve
                   adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
                else !if (rh >= 35.5 && rh <= 54.5)    !/* 1000md() && 1100md() */
                   current_hr = g41000md(mc1600)
                   tran_str   = g41000md(5.0)
                   tran_end   = g31100md(30.0)
                   tran_dist  = ((mc1600 - lemcb))
                   tran_crve  = tran_str + &
                             interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                   next_hr    = tran_crve
                   adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
               end if
            else if (rh < 38.0) then
               if (rh >= 35.5) then ! && rh <= 54.5)   !/* 1000lo() && 1100md() */
                  current_hr = g41000lo(mc1600)
                  tran_str   = g41000lo(5.0)
                  tran_end   = g31100md(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  next_hr    = tran_crve
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
               else !if (rh < 35.5)                    !/* 1000lo() && 1100lo() */
                  current_hr = g41000lo(mc1600)
                  tran_str   = g41000lo(5.0)
                  tran_end   = g31100lo(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  next_hr    = tran_crve
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
               end if    
            end if
         else if (mc1600 > 30.0) then
            if (rh > 57.0) then                         !/* 1000hi() && 1100hi() */
               current_hr    = g41000hi(mc1600)
               next_hr       = g31100hi(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
            else if (rh >= 38.0) then ! .and. rh <= 57) then
               if (rh > 54.5) then                     !/* 1000md() && 1100hi() */ 
                  if (hour <= 1030) then
                     adj_mc = g41000md(mc1600)
                  else
                     adj_mc = g31100hi(mc1600)
                  end if
               else ! if (rh >= 35.5 && rh <= 54.5)    !/* 1000md() && 1100md() */ 
                  current_hr    = g41000md(mc1600)
                  next_hr       = g31100md(mc1600)
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
               end if
            else if (rh < 38.0) then
               if (rh >= 35.5) then !rh <= 54.5)       !/* 1000()lo && 1000()md */ 
                  if (hour <= 1030) then
                     adj_mc = g41000lo(mc1600)
                  else
                     adj_mc = g31100md(mc1600)
                  end if
               else !if (rh < 35.5)                    !/* 1000lo() && 1100lo() */
                  current_hr    = g41000lo(mc1600)
                  next_hr       = g31100lo(mc1600)
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
               end if
            end if
         end if

      else if (hour >= 1100 .and. hour < 1200) then
         if (mc1600 <= 5.0) then                       !/* 1000**() */
            if (rh > 57.0) then
               adj_mc = g41000hi(mc1600)
            else if (rh >= 38.0) then ! && rh <=57.0) then
               adj_mc = g41000md(mc1600)
            else !if (rh < 38.0)
               adj_mc = g41000lo(mc1600)
            end if
         else if (mc1600 > 5.0 .and. mc1600 <= 30.0) then
            if (rh > 57.0) then                           !/* 1100hi() && 1159hi() */
               tran_str   = g41000hi(5.0)
               tran_end   = g31100hi(30.0)
               tran_dist  = ((mc1600 - lemcb))
               tran_crve  = tran_str + &
                        interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
               current_hr = tran_crve

               tran_end   = g41159hi(30.0)
               tran_dist  = ((mc1600 - lemcb))
               tran_crve  = tran_str + &
                        interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
               next_hr    = tran_crve
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
            else if (rh >= 38.0) then ! && rh <= 57.0)
               if (rh > 54.5) then                      !/* 1100md() && 1159hi() */
                  tran_str   = g41000md(5.0)
                  tran_end   = g31100hi(30.0)
                  tran_dist  = ((mc1600-lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  current_hr = tran_crve

                  tran_end   = g41159hi(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  next_hr    = tran_crve
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
               else !if (rh >= 35.5 && rh <= 54.5)     !/* 1100md() && 1159md() */
                  tran_str   = g41000md(5.0)
                  tran_end   = g31100md(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  current_hr = tran_crve

                  tran_end   = g41159md(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  next_hr    = tran_crve
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
               end if
            else !if (rh < 38.0)
               if (rh >= 35.5) then ! && rh <= 54.5)        !/* 1000lo() && 1100md() */
                  tran_str   = g41000lo(5.0)
                  tran_end   = g31100md(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  current_hr = tran_crve
                  
                  if (rh <= 52.0) then
                     tran_end   = g41159md(30.0)
                  else
                     tran_end   = g41159hi(30.0)
                  end if
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  next_hr    = tran_crve
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
               else  !if (rh < 35.5)                   !/* 1000lo() && 1100lo() */
                  tran_str   = g41000lo(5.0)
                  tran_end   = g31100lo(30.0)
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  current_hr = tran_crve

                  if (rh <  33.0) then
                     tran_end   = g41159lo(30.0)
                  else
                     tran_end   = g41159md(30.0)
                  end if
                  tran_dist  = ((mc1600 - lemcb))
                  tran_crve  = tran_str + &
                           interp(tran_str, tran_end, tran_dist, (uemcb - lemcb))
                  next_hr    = tran_crve
                  adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
               end if
            end if
         else if (mc1600 > 30.0) then
            if (rh < 35.5) then
               if (rh < 33.0) then
                  next_hr = g41159lo(mc1600)
               else  ! if (rh >= 33.0 && rh <= 52.0) then
                  next_hr = g41159md(mc1600)
               end if
               current_hr = g31100lo(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
               if (adj_mc > current_hr) adj_mc = current_hr
            else if (rh >=35.5 .and. rh <= 54.5) then
               if (rh > 52.0) then
                  next_hr = g41159hi(mc1600)
               else ! if (rh >=33.0 && rh <= 52.0)
                  next_hr = g41159md(mc1600)
               end if
               current_hr = g31100md(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
            else if (rh > 54.5) then
               next_hr    = g41159hi(mc1600)
               current_hr = g31100hi(mc1600)
               adj_mc = current_hr + interp(current_hr, next_hr, minutes, 59.0)
            end if
         end if
      else if (hour >= 1200 .and. hour < 1300) then
         if (mc1600 < 21.0) then
            next_hr       = g31300a(mc1600)
            current_hr    = g31200a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else if (mc1600 >= 21.0) then
            if (mc1600 < 22.0) then
               next_hr = g31300a(mc1600)
            else if (mc1600 >= 22.0) then
               next_hr = g31300b(mc1600)
            end if
            current_hr = g31200b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         end if
!
      else if (hour >= 1300 .and. hour < 1400) then
         if (mc1600 < 22.0) then
            next_hr       = g31400a(mc1600)
            current_hr    = g31300a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         else !if (mc1600 >=22)
            if (mc1600  < 23.0) then
               next_hr = g31400a(mc1600)
            else ! (mc1600 >= 23)
               next_hr = g31400b(mc1600)
            end if
            current_hr = g31300b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         end if
!
      else if (hour >= 1400 .and. hour < 1500) then
         if (mc1600 < 23.0) then
            next_hr       = g31500a(mc1600)
            current_hr    = g31400a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         else ! if (mc1600 >= 23.0) then
            next_hr       = g31500b(mc1600)
            current_hr    = g31400b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         end if
!
      else if (hour >= 1500 .and. hour < 1600) then
         if (mc1600 < 23.0) then
            next_hr       = mc1600
            current_hr    = g31500a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         else ! if (mc1600 >= 23.0) then
            next_hr       = mc1600
            current_hr    = g31500b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         end if
!
      else if (hour >= 1600 .and. hour < 1700) then
         if (mc1600 < 40.0) then
            next_hr     = g31700a(mc1600)
            current_hr    = mc1600
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         else !if (mc1600 >= 40.0)
            next_hr       = g31700b(mc1600)
            current_hr    = mc1600
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.0)
         end if
!
      else if (hour >= 1700 .and. hour < 1800) then
         if (mc1600 < 40.0) then
            next_hr       = g31800a(mc1600)
            current_hr    = g31700a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else !if (mc1600 >=40.0)
            next_hr       = g31800b(mc1600)
            current_hr    = g31700b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         end if
!
      else if (hour >= 1800 .and. hour < 1900) then
         if (mc1600 < 40.0) then
            next_hr       = g31900a(mc1600)
            current_hr    = g31800a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else ! if (mc1600 >=40.0)
            if (mc1600  < 42.0) then
               next_hr = g31900a(mc1600)
            else  !if (mc1600 >= 42.0)
               next_hr = g31900b(mc1600)
            end if
            current_hr = g31800b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         end if
!
      else if (hour >= 1900 .and. hour < 2000) then
         if (mc1600 < 42.0) then
            next_hr       = g32000a(mc1600)
            current_hr    = g31900a(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         else !if (mc1600 >=42.0)
            if (mc1600  < 49.0) then
                next_hr = g32000a(mc1600)
            else ! if (mc1600 >= 49.0)
                next_hr = g32000b(mc1600)
            end if
            current_hr  = g31900b(mc1600)
            adj_mc = current_hr + interp(current_hr, next_hr, minutes, 60.)
         end if
!
      else if (hour == 2000) then
         if (mc1600 < 49.0) then
            adj_mc = g32000a(mc1600)
         else ! if (mc1600 >=49.0)
            adj_mc = g32000b(mc1600)
         end if

      else if (hour > 2000 .or. hour < 600) then    !/*---- Extrapolate ------*/
          if (mc1600 < 49.0) then
             current_hr = g32000a(mc1600)
          else !if (mc1600 >=49.0)
             current_hr = g32000b(mc1600)
          end if

          next_hr = g4600hi(mc1600)      !/* Uses High RH Curve as Default */

          if (hour < 2000) hour = hour + 2400
          hrs_since_2000 = ((hour - 2000) - int(minutes)) / 100
          min_since_2000 = real(hrs_since_2000 * 60) + minutes
          adj_mc  = current_hr + &
                        interp(current_hr, next_hr, min_since_2000, 600.)

      else
         adj_ffmc = -1.0
         return
      end if

      adj_ffmc = 59.5 * (250 - adj_mc) / (147.2 + adj_mc)
      adj_ffmc = max(adj_ffmc, 0.0)
      adj_ffmc = min(adj_ffmc, 101.0)
 
      return
   end function bdl_ffmc

!/*--------------------------------------------------------------------------*/
!/*  EqFFMC.c
!This routine calculates an FFMC in equilibrium with the environment.
!It is likely more closer to reality than CVW's original hourly calculations

!Equation numbers are based on:

!Van Wagner, C.E., 1985: Equations and FORTRAN Program for the Canadian Forest Fire Weather Index System.
!    Can. For. Serv., Ottawa, Ont. For. Tech. Rep. 33. 18 pp.
!*/
   real(kind=4) function eq_ffmc(temp, rh, rf, oldffmc) result(adj_ffmc)
      implicit none

      real(kind=4), intent(in) :: temp, rh, rf, oldffmc

      real(kind=4) :: mo, m, mr, ed, ew, aa, bb

      if (oldffmc < 0.0 .or. oldffmc > 101.0 .or. rh < 0.0 .or. rh > 100.0 .or. &
          temp < -50. .or. temp > 50.0 .or. rf < 0.0 .or. rf > 200.0) then
         adj_ffmc = -99.0
!
      else
         mo = 147.2 * (101.0 - oldffmc) / (59.5 + oldffmc)          !  /* 1 */

         aa = 0.18 * (21.1 - temp) * (1.0 - exp(-0.115 * rh))
         bb = exp((rh - 100.) / 10.)
        
         ed = 0.942 * rh**0.679 + 11.0 * bb + aa                    ! /* 4 */
         ew = 0.618 * rh**0.753 + 10.0 * bb + aa                    ! /* 5 */

         m = mo

!/* Instantaneous recovery */
         if (mo > ed) m = ed                ! /* 8 with the decay term removed */
         if (mo < ew) m = ew                ! /* 9 with the decay term removed */

!/* wetting phase is now conducted on equilibrium mc following FFMC calculations*/
         if (rf > 0.0) then
            mo = m
            mr = mo + 42.5 * rf * exp(-100.0 / (251.0 - mo)) * &
                            (1.0 - exp(-6.93 / rf))               ! /* 3a */
            if (mo > 150.0) then 
               mr = mr + (0.0015*(mo-150.))**2 * sqrt(rf)         ! /* 3b */
            end if
            mr = min(mr, 250.0)

            m = mr
         end if

!/* new FFMC */
         adj_ffmc = 59.5 * (250.0 - m) / (147.2 + m)                  ! /* 10 */
      end if
      
      return
   end function eq_ffmc
!/*--------------------------------------------------------------------------*/

!/*--------------------------------------------------------------------------*/
!/*  hrffmc.c - Tue Mar 23 1993
!: mwotton@pnfi.forestry.ca (Mike Wotton)
!  kanderson@nofc.forestry.ca
!/*--------------------------------------------------------------------------*/
   real(kind=4) function hourly_ffmc(temp, rh, wind, rain, oldffmc) result(adj_ffmc)
      implicit none

      real(kind=4), intent(in)  :: temp, rh, wind, rain, oldffmc

      real(kind=4) :: mo, ed, ew, moew, moed, xm, a1, e, moe, xkd, aa, bb
      real(kind=4), parameter :: rf = 42.5, drf = 0.0579

      if ((oldffmc < 0.0 .or. oldffmc > 101.0) .or. &
          (temp < -50.0 .or. temp > 50.0) .or. (rh < 0.0 .or. rh > 100.0) .or. &
          (wind < 0.0 .or. wind > 200.0) .or. (rain < 0.0 .or. rain > 200.0)) then
         adj_ffmc = -999.0
!
      else
      
         mo = 147.2 * (101.0 - oldffmc) / (59.5 + oldffmc)

         if (rain /= 0.0) then
            mo = mo + rain * rf * exp(-100.0 / (251.0 - mo)) * &
                     (1.0 - exp(-6.93 / rain))
            mo = min(mo, 250.0)
         end if

         aa = 0.18 * (21.1 - temp) * (1.0 - 1.0 / exp(0.115 * rh))
         bb = exp((rh - 100.0) / 10.0)
         ed = 0.942 * (rh**0.679) + (11.0 * bb) + aa

         moed = mo - ed
           ew = 0.618 * (rh**0.753) + (10.0 * bb) + aa

         moew = mo - ew

         if (moed == 0.0 .or. (moew >= 0.0 .and. moed < 0.0)) then
            xm = mo
!            if (moed == 0.0) e = ed
!            if (moew >= 0.0) e = ew
         else
            if (moed > 0.0) then
               a1 = rh * 1.0e-2
               e = ed
               moe = moed
            else
                a1 = (100.0 - rh) * 1.0e-2
                 e = ew
               moe = moew
            end if

            xkd = 0.424 * (1.0 - (a1**1.7)) + (0.0694 * sqrt(wind) * (1.0 - (a1**8)))
            xkd = xkd * drf * exp(0.0365 * temp)
             xm = e + moe * exp(-2.303 * xkd)
!            mo = xm
         end if

         adj_ffmc = 59.5 * (250.0 - xm) / (147.2 + xm)
      end if
      
      return
   end function hourly_ffmc
!/*------------------------- Curve Equations ----------------------------------*/

!/*--------------------------------------------------------------------------*/
   real(kind=4) function interp(current, next, rmin, timestep) result(slope)
     implicit none
     real(kind=4), intent(in)  :: current, next, rmin, timestep
!/*--------------------------------------------------------------------------*/
      slope = 0.0
      slope = (next - current) / timestep
      slope = slope * rmin

     return
   end function interp

!/*--------------------------------------------------------------*/
   real(kind=4) function g31200a(x) result(y_sq)
!/*--------------------------------------------------------------*
!   Eqn# 7113  ûy=(a+cx+exý)/(1+bx+dxý)
! *--------------------------------------------------------------*/
     implicit none
     real(kind=4), intent(in)  :: x
     
     real(kind=4) :: y
      y = (1.460075955727318 + x * (0.2815668298889479 +    &
                               x * -0.01282068681081192)) / &
          (1.0 - x * (0.0003907916679793910 + x * 0.001539833829955819))
      y_sq = y * y
     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31200b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4740  y=a+bx+cx3+dûx+eexp-x
! *--------------------------------------------------------------*/
     implicit none
     real(kind=4), intent(in) :: x
     
     real(kind=4) :: x1, x2, x3, x4
      x1 = x
      x2 = x * x * x
      x3 = sqrt(x)
      x4 = exp(-x)
       y = -60.05817857808971 - 0.7922650736402593 * x1 + &
            1.049360630941214E-05 * x2 + 24.04228773333402 * x3 - &
            4790555393.429649 * x4
     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31300a(x) result(y_sq)
!/*--------------------------------------------------------------*
!   Eqn# 7114  ûy=(a+cx+exý)/(1+bx+dxý+fx3)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: y

     y = (1.255216373497914 + x * (0.3580951800872066 +    &
          x * -0.01642423056280615)) /                     &
         (1.0 + x * (0.02292170670740257 + x * (-0.003331105559008215 + &
          x * 3.056635174762784E-05)))
     y_sq = y * y
     
     return
   end function

!/*--------------------------------------------------------------*/
  real(kind=4) function g31300b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4560  y=a+bx+cxýlnx+d/ûx+elnx/x
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x * log(x)
     x3 = 1.0 / sqrt(x)
     x4 = log(x) / x
      y = 806.4657627391588 - 1.491623456269165 * x1 +          &
          0.0008873186770365859 * x2 - 11465.74581162260 * x3 + &
          12093.78039899465 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31400a(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4599  y=a+bx+cxýûx+dexpx+elnx
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x * sqrt(x)
     x3 = exp(x)
     x4 = log(x)
      y = 0.9082173869418283 + 0.9897247521440231 * x1 + &
          0.001041606201107530 * x2 + 4.633997000563253E-11 * x3 - &
          0.005581970517308027 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31400b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4874  y=a+bx+cûxlnx+dx/lnx+e/xý
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = sqrt(x) * log(x)
     x3 = x / log(x)
     x4 = 1.0 / (x * x)
      y = 6403.107752693009 + 352.7042531428901 * x1 +      &
          873.3642943916943 * x2 - 3766.492574730166 * x3 + &
          3580.933366117637 * x4

     return
   end function


!/*--------------------------------------------------------------*/
   real(kind=4) function g31500a(x) result(y_root)
!/*--------------------------------------------------------------*
!   Eqn# 6134  yý=a+bx+cxý+dx3+ex4+fx5
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x

     real(kind=4) :: y
     
     y = 0.2487113274993664 + x * (0.9002141389059909 +        &
         x * (0.9658994322020696 + x * (0.007692506399975233 + &
         x * (-0.0003031672886149768 + x * 1.121650772708818E-05))))
     y_root = sqrt(y)

     return
   end function


!/*--------------------------------------------------------------*/
   real(kind=4) function g31500b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4874  y=a+bx+cûxlnx+dx/lnx+e/xý
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = sqrt(x) * log(x)
     x3 = x / log(x)
     x4 = 1.0 / (x * x)
      y = 3201.553847063860 + 176.8521250328144 * x1 +      &
          436.6821438580977 * x2 - 1883.246271847646 * x3 + &
          1790.467302430289 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31700a(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4355  y=a+bx+cxý+dxýûx+eexp-x
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x
     x3 = x * x * sqrt(x)
     x4 = exp(-x)
      y = 0.3578377562139084 + 1.043214752557405 * x1 -  &
          0.001370303280372795 * x2 - 8.509201577454290E-05 * x3 + &
          0.1580591882752711 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31700b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4609  y=a+bx+cxýûx+dûxlnx+ex/lnx
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x * sqrt(x)
     x3 = sqrt(x) * log(x)
     x4 = x / log(x)
      y = 2776.473019059725 + 153.8288087656879 * x1 -          &
          0.0001010951827792066 * x2 + 371.9483315265095 * x3 - &
          1620.093039992771 * x4
      
     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31800a(x) result(y_root)
!/*--------------------------------------------------------------*
!   Eqn# 6132  yý=a+bx+cxý+dx3
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x
     
     real(kind=4) :: y

     y = 1.071980333417974 + x * (1.360477850187146 +    &
         x * (1.201854443540851 + x * -0.008273056193800057))
     y_root = sqrt(y)

     return
   end function


!/*--------------------------------------------------------------*/
   real(kind=4) function g31800b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4609  y=a+bx+cxýûx+dûxlnx+ex/lnx
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x * sqrt(x)
     x3 = sqrt(x) * log(x)
     x4 = x / log(x)
      y = 5552.947642901565 + 306.6577058442273 * x1 -          &
          0.0002021904025218839 * x2 + 743.8968800012295 * x3 - &
          3240.187020515272 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31900a(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4382  y=a+bx+cxý+dexpx+eexp-x
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x
     x3 = exp(x)
     x4 = exp(-x)
     y  = 1.948509314122392 + 1.124895722282613 * x1 -             &
          0.005100676601005492 * x2 + 8.905547777852358E-20 * x3 + &
          0.2620286583420982 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31900b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4173  y=a+bx+cxûx+dxý+exýûx
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * sqrt(x)
     x3 = x * x
     x4 = x * x2
     y = 28.76729089692713 - 1.511951568164489 * x1 +         &
         0.4217514053565216 * x2 - 0.02633183039679114 * x3 + &
         0.0005859070258838156 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g32000a(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4341  y=a+bx+cxý+dxýûx+ex3
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x
     x3 = x * x * sqrt(x)
     x4 = x * x * x
     y = 3.367449306255549 + 1.083974300132994 * x1 +            &
         0.007668482694768896 * x2 - 0.003614577217487989 * x3 + &
         0.0002675912148501817 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g32000b(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 4755  y=a+bx+cx3+d/lnx+eexp-x
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: x1, x2, x3, x4

     x1 = x
     x2 = x * x * x
     x3 = 1.0 / log(x)
     x4 = exp(-x)
     y = -111.6584390050989 + 1.238144218990064 * x1 -         &
         1.739987959934465E-06 * x2 + 379.1717488220770 * x3 - &
         5.512040174539097E+20 * x4

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4600lo(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
!  lognorm(x, 192.8242798917865, 1.748892433468639, &
!              6.966628145396413, 65.4192874088499)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 192.8242798917865) / 1.748892433468639
     y = 6.966628145396413 + 65.4192874088499 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4600md(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 1610.269344585912, 2.412647413894015, &
!            11.80584752343706, 145.1618675379673)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 1610.269344585912) / 2.412647413894015
     y = 11.80584752343706 + 145.1618675379673 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4600hi(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
!  lognorm(x, 2159.088827996704, 2.390534289077195, &
!             14.89281072738765, 194.5261398416724)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 2159.088827996704) / 2.390534289077195
     y = 14.89281072738765 + 194.5261398416724 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4700lo(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm (x, 216.2009556040568, 1.81202656162805, &
!             6.221403214623654, 61.83553855583578)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 216.2009556040568) / 1.81202656162805
     y = 6.221403214623654 + 61.83553855583578 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4700md(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 843.7712566576433, 2.143231970736562, &
!            10.62087345466205, 120.3071747968938)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 843.7712566576433) / 2.143231970736562
     y = 10.62087345466205 + 120.3071747968938 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4700hi(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
!  lognorm(x, 1308.435220755216, 2.269455130483546, &
!             12.5226863544613, 160.3933412357393)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 1308.435220755216) / 2.269455130483546
     y = 12.5226863544613 + 160.3933412357393 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4800lo(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 253.0830911423353, 1.896023727845998, &
!             5.45448266757961, 58.6461017619521)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 253.0830911423353) / 1.896023727845998
     y = 5.45448266757961 + 58.6461017619521 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4800md(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
!  lognorm(x, 547.1226761430743, 1.946001002988131, &
!             9.179219104642618, 105.6311973262788)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 547.1226761430743) / 1.946001002988131
     y = 9.179219104642618 + 105.6311973262788 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4800hi(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 848.3773713172653, 2.154869886473, &
!            10.21004190871538, 136.7485496998884)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 848.3773713172653) / 2.154869886473
     y = 10.21004190871538 + 136.7485496998884 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4900lo(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 206.2626505491258, 1.814962091988218, &
!            3.966946508633288, 47.66100216364151)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 206.2626505491258) / 1.814962091988218
     y = 3.966946508633288 + 47.66100216364151 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4900md(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 544.0978143522045, 2.000706807758964, &
!            6.381382418031537, 88.54320780795105)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 544.0978143522045) / 2.000706807758964
     y = 6.381382418031537 + 88.54320780795105 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g4900hi(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 1192.45753897097,  2.28873947108425, &
!            9.099751896903934, 127.6089430083535)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 1192.45753897097) / 2.28873947108425
     y = 9.099751896903934 + 127.6089430083535 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g41000lo(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
!  lognorm(x, 161.725408818154, 1.710574763630868, &
!              2.509991704520077, 37.42399134558683)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 161.725408818154) / 1.710574763630868
     y = 2.509991704520077 + 37.42399134558683 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g41000md(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 525.206855268004, 2.070941812038297, &
!            3.497497088121431, 71.24103374128234)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 525.206855268004) / 2.070941812038297
     y = 3.497497088121431 + 71.24103374128234 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g41000hi(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
! lognorm(x, 2357.682971305255, 2.538559054522592, &
!            7.891852884597797, 126.9570676813226)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x

     real(kind=4) :: n 

     n = log(x / 2357.682971305255) / 2.538559054522592
     y = 7.891852884597797 + 126.9570676813226 * exp(-0.5 * n * n)

     return
   end function
!/*--------------------------------------------------------------*/
   real(kind=4) function g31100lo(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 7203  y=(a+clnx+elnxý)/(1+blnx+dlnxý)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x
     real(kind=4) :: lx 
     
     lx = log(x)
     y = (1.291826915669442 + lx * (0.1581477301288331 +  &
          lx * 0.3560512548608463)) /                     &
         (1.0 + lx * (-0.3816865757039243 + lx * 0.05135364688196362))

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31100md(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 8005  y=a+bexp(-0.5(ln(x/c)/d)^2) [Log-Normal]
!  lognorm(x, 461.9583952258151, 2.149631748294819, &
!             0.5145364592257406, 53.63085253740166)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent (in) :: x

     real(kind=4) :: n 

     n = log(x / 461.9583952258151) / 2.149631748294819
     y = 0.5145364592257406 + 53.63085253740166 * exp(-0.5 * n * n)

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g31100hi(x) result(y)
!/*--------------------------------------------------------------*
!   Eqn# 7203  y=(a+clnx+elnxý)/(1+blnx+dlnxý)
! *--------------------------------------------------------------*/
     implicit none

     real(kind=4), intent(in) :: x
     real(kind=4) :: lx 

     lx = log(x)
     y = (7.934004974027670 + lx * (-0.2983586873533423 + &
          lx * 0.5901343665313362)) /                     &
         (1.0 + lx * (-0.2113457982134224 + lx * 0.01580693446540503))

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g41159lo(x) result(adjmc)
!/*--------------------------------------------------------------*
! Uses Extrapolation to estimate moisture content at 1159.
! *--------------------------------------------------------------*/
     implicit none
     real(kind=4), intent (in) :: x
     
     real(kind=4) :: diff

     diff = g41000lo(x) - g31100lo(x)
     if (diff < 0.0) diff = diff * -1.0
     adjmc = (g31100lo(x) - ((59.0 / 60.0) * diff))

     return
   end function

!/*--------------------------------------------------------------*/
   real(kind=4) function g41159md(x) result(adjmc)
!/*--------------------------------------------------------------*
! Uses Extrapolation to estimate moisture content at 1159.
! *--------------------------------------------------------------*/
     implicit none
     real(kind=4), intent (in) :: x
     
     real(kind=4) :: diff

     diff = g41000md(x) - g31100md(x)
     if (diff < 0.0) diff = diff * -1.0
     adjmc = (g31100md(x) - ((59.0 / 60.0) * diff))
 
     return
   end function


!/*--------------------------------------------------------------*/
  real(kind=4) function g41159hi(x) result(adjmc)
!/*--------------------------------------------------------------*
! Uses Extrapolation to estimate moisture content at 1159.
! *--------------------------------------------------------------*/
     implicit none
     real(kind=4), intent (in) :: x
     
     real(kind=4) :: diff

     diff = g41000hi(x) - g31100hi(x)
     if (diff < 0.0) diff = diff * -1.0
     adjmc = (g31100hi(x) - ((59.0 / 60.0) * diff))
 
     return
   end function


end module mach_cffeps_diurnal_mod
