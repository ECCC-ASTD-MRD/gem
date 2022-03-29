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

!**s/r nest_intt -- Linear interpolation in time of nesting data

      subroutine nest_intt()
      use cstv
      use gem_options
      use lam_options
      use mem_nest
      use glb_ld
      use step_options
      use omp_timing
      use ptopo
      use tr3d
      use, intrinsic :: iso_fortran_env
      implicit none

      include 'mpif.h'
      include 'rpn_comm.inc'

      integer,external ::  newdate
      
      character(len=16) :: datev, previous_S
      logical :: stag_L
      integer :: i,yy,mo,dd,hh,mm,ss,dum
      integer topo_diff,comm,err,gathV(Ptopo_numproc*Ptopo_ncolors)
      real(kind=REAL64) :: dayfrac,tx,dtf,a,b
      real(kind=REAL64), parameter :: one=1.0d0, &
                       sid=86400.0d0, rsid=one/sid
!
!     ---------------------------------------------------------------
!
      dayfrac = dble(Step_kount) * Cstv_dt_8 * rsid
      call incdatsd  (datev, Step_runstrt_S, dayfrac)
      call prsdate   (yy,mo,dd,hh,mm,ss,dum,datev)
      call pdfjdate2 (tx, yy,mo,dd,hh,mm,ss)

      dayfrac = Step_nesdt*rsid
      
      call gtmg_start (23, 'NEST_input', 20)

!$omp master
      if (tx < Lam_tdeb) then

         Lam_current_S  = Step_runstrt_S
         previous_S = Lam_current_S
         call prsdate   (yy,mo,dd,hh,mm,ss,dum,Step_runstrt_S)
         call pdfjdate2 (Lam_tdeb,yy,mo,dd,hh,mm,ss)
         do
            call incdatsd (datev, Lam_current_S,dayfrac)
            Lam_current_S = datev
            call prsdate   (yy,mo,dd,hh,mm,ss,dum,Lam_current_S)
            call pdfjdate2 (Lam_tfin, yy,mo,dd,hh,mm,ss)
            if (Lam_tfin >= tx) exit
            previous_S = Lam_current_S
            Lam_tdeb       = Lam_tfin
         end do
         Lam_current_S = previous_S
         Lam_tfin      = Lam_tdeb
         stag_L= .true.
         if ( previous_S == Step_runstrt_S) stag_L= .false.
         call nest_indata  (nest_u_fin, nest_v_fin , nest_w_fin, nest_t_fin ,&
                            nest_q_fin, nest_zd_fin, nest_s_fin, nest_tr_fin,&
                            nest_fullme_fin,stag_L,previous_S               ,&
                            l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr)
      end if

      dtf = 1.0d0
      if (tx > Lam_tfin) then
         dtf = (tx-Lam_tfin) * sid / Cstv_dt_8
         Lam_tdeb       = Lam_tfin

         call incdatsd (datev, Lam_current_S, dayfrac)
         Lam_current_S = datev
         call prsdate   (yy,mo,dd,hh,mm,ss,dum,Lam_current_S)
         call pdfjdate2 (Lam_tfin, yy,mo,dd,hh,mm,ss)

         nest_u_deb  = nest_u_fin
         nest_v_deb  = nest_v_fin
         nest_t_deb  = nest_t_fin
         nest_s_deb  = nest_s_fin
         nest_w_deb  = nest_w_fin
         nest_q_deb  = nest_q_fin
         nest_zd_deb = nest_zd_fin
         nest_fullme_deb = nest_fullme_fin
         nest_tr_deb = nest_tr_fin

         call nest_indata  (nest_u_fin, nest_v_fin , nest_w_fin, nest_t_fin ,&
                            nest_q_fin, nest_zd_fin, nest_s_fin, nest_tr_fin,&
                            nest_fullme_fin,.true.,Lam_current_S            ,&
                            l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr)

         topo_diff= 0
         if (maxval(abs(nest_fullme_fin-nest_fullme_deb))>0.) topo_diff= 1
         comm = RPN_COMM_comm ('MULTIGRID')
         call MPI_Allgather ( topo_diff,1,MPI_INTEGER,gathV,1,&
                              MPI_INTEGER,comm,err)
         Vtopo_mustadj_L= sum(gathV) > 0

      end if
!$omp end master
!$OMP BARRIER
      call gtmg_stop (23)
      call gtmg_start (24, 'NEST_tint', 20)

      b = (tx - Lam_tdeb) / (Lam_tfin - Lam_tdeb)
      a = one - b
      
!$omp do
      do i= lbound(nest_now,1), ubound(nest_now,1)
         nest_now(i) = a*nest_deb(i) + b*nest_fin(i)
      end do
!$omp enddo
      call gtmg_stop (24)
      
!      call nest_glbstat ((/'now','deb','fin'/),3)

!     
!     ---------------------------------------------------------------
!
      return
      end

