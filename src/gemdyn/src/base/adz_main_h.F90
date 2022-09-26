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

      subroutine adz_main_h ()
      use ISO_C_BINDING
      use dyn_fisl_options
      use cstv
      use gem_options
      use adz_mem
      use adz_interp_hlt_mod
      use dynkernel_options
      use mem_tstp
      implicit none

      integer :: n
      type(Adz_pntr_stack), dimension(3), target :: stack
!
!     ---------------------------------------------------------------
!
      call adz_traject (Cstv_dt_8)

      stack(1)%src =>  orhsu_ext
      stack(1)%dst =>  rhsu
      call adz_tricub_hlt ( stack,1,Adz_pmu,Adz_cpntr_q,Adz_num_u,&
                            Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0 )

      stack(1)%src => orhsv_ext
      stack(1)%dst =>  rhsv
      call adz_tricub_hlt ( stack,1,Adz_pmv,Adz_cpntr_q,Adz_num_v,&
                            Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0 )

      stack(1)%src => orhsc_ext
      stack(1)%dst =>  rhsc
      call adz_tricub_hlt ( stack,1,Adz_pm ,Adz_cpntr_q,Adz_num_q,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )

      stack(1)%src => orhst_ext
      stack(1)%dst =>  rhst
      stack(2)%src => orhsf_ext
      stack(2)%dst =>  rhsf
      n= 2
      if((.not.Dynamics_hydro_L) .or. Dynamics_hauteur_L) then
         n= 3
         stack(n)%src => orhsw_ext
         stack(n)%dst =>  rhsw
      endif
      call adz_tricub_hlt ( stack,n,Adz_pt ,Adz_cpntr_t,Adz_num_t,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0t )
!$OMP BARRIER

!$omp single
      call rpn_comm_xch_halo( rhsu, l_minx,l_maxx,l_miny,l_maxy, l_ni, l_nj,2*G_nk,&
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
!$omp end single
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_main_h
