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

!**s/p set_var - initialize list of variables to output
!
      integer function set_var (F_argc,F_argv_S,F_cmdtyp_S,F_v1,F_v2)

      use out_options
      use tdpack
      use glb_ld
      use lun
      use out3
      use levels
      use outp
      use outd
      use outc
      use outgrid
      use timestep
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_argc,F_v1,F_v2
      character(len=*), intent(in) :: F_argv_S(0:F_argc),F_cmdtyp_S

!object
!       initialization of the common blocks OUTD,OUTP. This function is
!       called when the keyword "sortie" is found in the first word
!       of the directives in the input file given in the statement
!       "process_f_callback". This feature is enabled by the
!       ARMNLIB "rpn_fortran_callback" routine (called in "srequet")
!       which allows a different way of passing user directives than
!       the conventional FORTRAN namelist. This function will process
!       the following example command read from the named input file.
!
! ie:   sortie([UU,VV,TT],levels,2,grid,3,steps,1)
!       sortie([PR,PC,RR],grid,3,steps,2,levels,1)
!
!       The "rpn_fortran_callback" routine will process the above
!       statement and return 5 arguments to this function. For more
!       information to how this is processed, see "SREQUET".
!
!
!arguments
!  Name        I/O                 Description
!----------------------------------------------------------------
! F_argc       I    - number of elements in F_argv_S
! F_argv_S     I    - array of elements received
!                     if F_argv_S(ii) contains "[", the value in this
!                     argument indicates number of elements following it
! F_cmdtyp_S   I    - character command type - not used
! F_v1         I    - integer parameter 1 - not used
! F_v2         I    - integer parameter 2 - not used
!----------------------------------------------------------------
!
!Notes:
!    ie:   sortie([UU,VV,TT],levels,2,grid,3,steps,1)
!          sortie([PR,PC,RR],grid,3,steps,2,levels,1)
!          sortie_p([AFSI,AFSV],grid, 1, levels, 1, steps, 1)
!          sortie_p([AFSI],grid, 1, levels, 1, steps, 3,average)
!          sortie_p([AFSV],grid, 1, levels, 1, steps, 3,accum)
!
! sortie([vr1,vr2,vr3,...],levels,[levelset],grid,[gridset],steps,[stepset])
!
!  vr1,vr2,vr3... - set of variable names to output (max of 60)
!  levelset - levelset number to use for this set of variables
!  gridset  - gridset number to use for this set of variables
!  stepset  - stepset number (timestep set) to use for this set of variables
!
!  For each "sortie" command, the levelset, gridset and stepset must be
!  specified or an error will occur.
!


      character(len=5) :: stuff_S
      character(len=4) :: string4
      character(len=16):: string16
      integer :: levset,stepset,gridset,varmax
      integer :: i, j, ii, jj
      logical :: accum_L,avg_L
!
!----------------------------------------------------------------
!
      if (Lun_out > 0) then
         write(Lun_out,*)
         write(Lun_out,*) F_argv_S
      end if
      set_var=0

      if (index(F_argv_S(1),'[') > 0) then
         stuff_S=F_argv_S(1)
         read(stuff_S(2:4),*) varmax
         if (varmax > MAXELEM) then
             if (Lun_out > 0) write(Lun_out,*) 'SET_VAR WARNING: too many variables (over 200) in this request set'
             set_var=1
             return
         endif
      else
         if (Lun_out > 0) then
           write(Lun_out,*) 'SET_VAR WARNING: syntax incorrect'
         end if
         set_var=1
         return
      end if
!
!     Check if chosen levels,grid and timestep sets are valid
!
      levset=-1
      gridset=-1
      stepset=-1
      accum_L=.false.
      avg_L=.false.
      do i=varmax+2, F_argc
         if (F_argv_S(i) == 'levels') then
            read(F_argv_S(i+1),*) levset
         else if (F_argv_S(i) == 'grid') then
            read(F_argv_S(i+1),*) gridset
         else if (F_argv_S(i) == 'steps') then
            read(F_argv_S(i+1),*) stepset
         else if (F_argv_S(i) == 'accum') then
            accum_L=.true.
         else if (F_argv_S(i) == 'average') then
            avg_L=.true.
         end if
      end do

      if (gridset < 0) then
         if (Lun_out > 0) then
            write(Lun_out,*) 'SET_VAR WARNING: no Grid chosen'
         end if
         set_var=1
         return
      else
         do i=1,OutGrid_sets
            if (gridset == OutGrid_id(i)) then
                gridset=i
                exit
            end if
         end do
         if (i > OutGrid_sets) then
            if (Lun_out > 0) then
                write(Lun_out,*) 'SET_VAR WARNING: invalid Grid set ID#'
            end if
            set_var=1
            return
         end if
      end if
      if (levset < 0) then
         if (Lun_out > 0) then
            write(Lun_out,*) 'SET_VAR WARNING: no Levels chosen'
         end if
         set_var=1
         return
      else
         do i=1,Level_sets
            if (levset == Level_id(i)) then
                levset=i
                exit
            end if
         end do
         if (i > Level_sets) then
            if (Lun_out > 0) then
               write(Lun_out,*) 'SET_VAR WARNING: invalid Level set ID#'
            end if
            set_var=1
            return
         end if
      end if
      if (stepset < 0) then
         if (Lun_out > 0) then
            write(Lun_out,*) 'SET_VAR WARNING: no Timesteps chosen'
         end if
         set_var=1
         return
      else
         do i=1,Timestep_sets
            if (stepset == Timestep_id(i)) then
                stepset=i
                exit
            end if
         end do
         if (i > Timestep_sets) then
            if (Lun_out > 0) then
               write(Lun_out,*) 'SET_VAR WARNING: invalid Timestep set ID#'
            end if
            set_var=1
            return
         end if
      end if
!
!     Store variables in variable sets
!
      if (F_argv_S(0) == 'sortie') then
         j = Outd_sets + 1
         if (j > MAXSET) then
            if (Lun_out > 0) then
               write(Lun_out,*) 'SET_VAR WARNING: too many OUTD sets'
            end if
            set_var=1
            return
         end if
!
         jj=0
         do ii=1,varmax
            jj = jj + 1
            call low2up  (F_argv_S(ii+1),string16)
            call low2up  (F_argv_S(ii+1),string4)
            Outd_varnm_S(jj,j)=string16
            Outd_var_S(jj,j)=string4
            if (Outd_var_S(jj,j)(1:3) == 'PW_') Outd_var_S(jj,j)= 'PW'//Outd_varnm_S(jj,j)(4:5)
            Outd_nbit(jj,j)    = Out3_nbitg
            Outd_filtpass(jj,j)= 0
            Outd_filtcoef(jj,j)= 0.0
            Outd_convmult(jj,j)= 1.0
            Outd_convadd (jj,j)= 0.0

            if (Outd_var_S(jj,j)(1:4) == 'LA  ') Outd_nbit(jj,j)= 32
            if (Outd_var_S(jj,j)(1:4) == 'LO  ') Outd_nbit(jj,j)= 32

            if (Outd_var_S(jj,j)(1:4) == 'GZ  ') Outd_convmult(jj,j)=0.1 / grav_8
            if (Outd_var_S(jj,j)(1:4) == 'ME  ') Outd_convmult(jj,j)=1.0 / grav_8
            if (Outd_var_S(jj,j)(1:4) == 'PX  '.or. &
                Outd_var_S(jj,j)(1:4) == 'P0  '.or. &
                Outd_var_S(jj,j)(1:4) == 'PT  '.or. &
                Outd_var_S(jj,j)(1:4) == 'PN  ') Outd_convmult(jj,j)= .01
            if (Outd_var_S(jj,j)(1:4) == 'VT  '.or. &
                Outd_var_S(jj,j)(1:4) == 'TT  '.or. &
                Outd_var_S(jj,j)(1:4) == 'TD  ') Outd_convadd (jj,j)= -tcdk_8
            if (Outd_var_S(jj,j)(1:4) == 'UU  '.or. &
                Outd_var_S(jj,j)(1:4) == 'VV  '.or. &
                Outd_var_S(jj,j)(1:4) == 'UV  ') Outd_convmult(jj,j)=1.0 / knams_8
         end do

         if (jj > 0) then
             Outd_sets       = j
             Outd_var_max(j) = jj
             Outd_grid(j)    = gridset
             Outd_lev(j)     = levset
             Outd_step(j)    = stepset
         else
             if (Lun_out > 0) write(Lun_out,1400)
         end if
      else if (F_argv_S(0) == 'sortie_p') then
         j = Outp_sets + 1
         if (j > MAXSET) then
            if (Lun_out > 0) then
               write(Lun_out,*) 'SET_VAR WARNING: too many OUTP sets'
            end if
            set_var=1
            return
         end if
!
         jj=0
         do ii=1,varmax
            jj = jj + 1
            call low2up  (F_argv_S(ii+1),string16)
            Outp_varnm_S(jj,j)=string16
            Outp_nbit(jj,j)  = Out3_nbitg
            Outp_filtpass(jj,j)= 0
            Outp_filtcoef(jj,j)= 0.0
            Outp_convmult(jj,j)= 1.0
            Outp_convadd (jj,j)= 0.0
            if (Outp_varnm_S(jj,j)(1:4) == 'LA') then
                Outp_nbit(jj,j)= 32
                Outp_convmult(jj,j)=180./pi_8
            end if
            if (Outp_varnm_S(jj,j)(1:4) == 'LO') then
                Outp_nbit(jj,j)= 32
                Outp_convmult(jj,j)=180./pi_8
            end if
            if (Outp_varnm_S(jj,j)(1:4) == 'SD') Outp_convmult(jj,j)=100.
            if (Outp_varnm_S(jj,j)(1:4) == 'SNDP') Outp_convmult(jj,j)=100.
            if (Outp_varnm_S(jj,j)(1:4) == 'SVDP') Outp_convmult(jj,j)=100.
         end do
         if (jj > 0) then
             Outp_sets         = j
             Outp_var_max(j)   = jj
             Outp_grid(j)      = gridset
             Outp_lev(j)       = levset
             Outp_step(j)      = stepset
             Outp_accum_L(j)   = accum_L
             Outp_avg_L  (j)   = avg_L
             if (Lun_out > 0) then
                write(Lun_out,*) '***PHY***Outp_sets=',Outp_sets
                write(Lun_out,*) 'Outp_var_max=',Outp_var_max(j)
                write(Lun_out,*) 'Outp_varnm_S=', &
                             (Outp_varnm_S(jj,j),jj=1,Outp_var_max(j))
                write(Lun_out,*) 'Outp_grid=',Outp_grid(j)
                write(Lun_out,*) 'Outp_lev=',Outp_lev(j)
                write(Lun_out,*) 'Outp_step=',Outp_step(j)
                if (Outp_accum_L(j))write(Lun_out,*)'Outp_accum_L=',Outp_accum_L(j)
                if (Outp_avg_L  (j))write(Lun_out,*)'Outp_avg_L='  ,Outp_avg_L  (j)
             end if
         else
             if (Lun_out > 0) write(Lun_out,1400)
         end if
      else if (F_argv_S(0) == 'sortie_c') then
         j = Outc_sets + 1
         if (j > MAXSET) then
            if (Lun_out > 0) then
               write(Lun_out,*) 'SET_VAR WARNING: too many OUTC sets'
            end if
            set_var=1
            return
         end if
!
         jj=0
         do ii=1,varmax
            jj = jj + 1
            Outc_varnm_S(jj,j)= F_argv_S(ii+1)
            Outc_nbit(jj,j)   = Out3_nbitg
            Outc_filtpass(jj,j)= 0
            Outc_filtcoef(jj,j)= 0.0
            Outc_convmult(jj,j)= 1.0
            Outc_convadd (jj,j)= 0.0
         end do
         if (jj > 0) then
             Outc_sets       = j
             Outc_var_max(j) = jj
             Outc_grid(j)    = gridset
             Outc_lev(j)     = levset
             Outc_step(j)    = stepset
             if (Lun_out > 0) then
                write(Lun_out,*) '***CHM***Outc_sets=',Outc_sets
                write(Lun_out,*) 'Outc_var_max=',Outc_var_max(j)
                write(Lun_out,*) 'Outc_varnm_S=', (Outc_varnm_S(jj,j),jj=1,Outc_var_max(j))
                write(Lun_out,*) 'Outc_grid=',Outc_grid(j)
                write(Lun_out,*) 'Outc_lev=',Outc_lev(j)
                write(Lun_out,*) 'Outc_step=',Outc_step(j)
             end if
         else
             if (Lun_out > 0) write(Lun_out,1400)
         end if
      end if
!
!----------------------------------------------------------------
!
 1400 format('SET_VAR - WARNING: NO VARIABLES DEFINED FOR THIS SET')
      return
      end
