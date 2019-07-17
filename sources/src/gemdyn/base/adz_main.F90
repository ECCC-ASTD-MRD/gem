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

      subroutine adz_main ( orhsu, rhsu, orhsv, rhsv, orhsc ,&
                            rhsc, orhst,  rhst, orhsf, rhsf ,&
                            orhsw, rhsw, Minx,Maxx,Miny,Maxy, Nk )
      use ISO_C_BINDING
      use adz_mem
      use adz_interp_rhs_mod
      use dynkernel_options
      use gem_timing
      use gmm_vt0
      use gmm_vt1
      use gmm_pw
      implicit none
#include <arch_specific.hf>

      integer,intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,NK),target,intent(in)  :: &
                             orhsu, orhsv, orhsc, orhst, orhsf, orhsw
      real, dimension(Minx:Maxx,Miny:Maxy,NK),target,intent(out) :: &
                             rhsu, rhsv, rhsc,  rhst, rhsf, rhsw
      integer :: n
!
!     ---------------------------------------------------------------
!
      call gemtime_start (30, 'ADZ_TRAJEC', 21)
      call adz_traject (pw_uu_moins, pw_vv_moins, zdt1, &
                           ut0     ,    vt0     , zdt0, &
                        l_minx,l_maxx,l_miny,l_maxy,l_nk)
      call gemtime_stop (30)

      call gemtime_start (31, 'ADZ_INTP_RH', 21)

      Adz_stack(1)%src => orhsu
      Adz_stack(1)%dst =>  rhsu
      call adz_tricub_rhs ( Adz_stack,1,Adz_pmu,Adz_cpntr_q,Adz_num_u,&
                            Adz_i0u,Adz_inu,Adz_j0,Adz_jn,Adz_k0 )

      Adz_stack(1)%src => orhsv
      Adz_stack(1)%dst =>  rhsv
      call adz_tricub_rhs ( Adz_stack,1,Adz_pmv,Adz_cpntr_q,Adz_num_v,&
                            Adz_i0,Adz_in,Adz_j0v,Adz_jnv,Adz_k0 )

      Adz_stack(1)%src => orhsc
      Adz_stack(1)%dst =>  rhsc
      call adz_tricub_rhs ( Adz_stack,1,Adz_pm ,Adz_cpntr_q,Adz_num_q,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0 )

      Adz_stack(1)%src => orhst
      Adz_stack(1)%dst =>  rhst
      Adz_stack(2)%src => orhsf
      Adz_stack(2)%dst =>  rhsf
      n= 2
      if(.not.Dynamics_hydro_L) then
         n= 3
         Adz_stack(n)%src => orhsw
         Adz_stack(n)%dst =>  rhsw
      endif
      call adz_tricub_rhs ( Adz_stack,n,Adz_pt ,Adz_cpntr_t,Adz_num_t,&
                            Adz_i0,Adz_in,Adz_j0,Adz_jn,Adz_k0t )

      call gemtime_stop (31)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_main
