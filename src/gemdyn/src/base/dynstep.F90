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

!**s/r dynstep -  Control of the dynamical timestep of the model

      subroutine dynstep()
      use dynkernel_options
      use dyn_fisl_options
      use gem_options
      use sol_options
      use gmm_vt1
      use gmm_vt2
      use gmm_pw
      use adz_mem
      use HORgrid_options
      use lun
      use step_options
      use omp_timing
      use theo_options
      implicit none

      logical, save :: done=.false.
      integer icn, keep_itcn, np
!
!     ---------------------------------------------------------------
!
      np= 3*G_nk
      if ( (.not.Dynamics_hydro_L) .or. (Dynamics_hauteur_L)) np= 4*G_nk+1
      call rpn_comm_xch_halo( ut1 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj ,np, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
      call rpn_comm_xch_halo( st1 , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,1   , &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
!      call rpn_comm_xch_halo( ut2 , l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj ,np, &
!                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
!      call rpn_comm_xch_halo( st2 , l_minx,l_maxx,l_miny,l_maxy,l_ni ,l_nj ,1   , &
!                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

      Adz_uu_ext(1:l_ni,1:l_nj,1:l_nk) = pw_uu_plus(1:l_ni,1:l_nj,1:l_nk)
      Adz_vv_ext(1:l_ni,1:l_nj,1:l_nk) = pw_vv_plus(1:l_ni,1:l_nj,1:l_nk)
      Adz_ww_ext(1:l_ni,1:l_nj,1:l_nk) =       zdt1(1:l_ni,1:l_nj,1:l_nk)
      call rpn_comm_xch_halo (Adz_uu_ext, Adz_lminx,Adz_lmaxx,Adz_lminy,Adz_lmaxy,&
                    l_ni,l_nj, 3*l_nk, Adz_halox,Adz_haloy, .false.,.false., l_ni,0)

      if (trim(Sol_type_S) == 'ITERATIVE_3D') then
         if (.not. done) then
!$omp parallel
            call matvec3d_init_hlt ()
!$omp end parallel
         endif
         done =.true.
      endif

      select case ( trim(Dynamics_Kernel_S) )
         case ('DYNAMICS_FISL_H')
            call fislh_dynstep()
            return
      end select

      if (Lun_debug_L) write(Lun_out,1000)
      call gtmg_start ( 10, 'DYNSTEP', 1 )

      keep_itcn = Schm_itcn

      if (Lun_debug_L) write(Lun_out,1005) Schm_itcn-1

      call psadj_init ( Step_kount )

      do icn = 1,Schm_itcn-1

         call tstpdyn (icn)

         call hzd_momentum()

      end do

      if (Lun_debug_L) write(Lun_out,1006)

      call tstpdyn (Schm_itcn)

      if (Ctrl_theoc_L .and. .not.Grd_yinyang_L) call theo_bndry ()

!$omp parallel
      call adz_tracers (.true.)

      call psadj_hlt ( Step_kount )

      call adz_tracers (.false.)
!$omp end parallel

      call t02t1()

      call HOR_bndry ()

      call canonical_cases ("VRD")

      call hzd_main ()

      if (Grd_yinyang_L) call yyg_blend()

      call pw_update_GW()
      call pw_update_UV()
      call pw_update_T()

      if ( Lctl_step-Vtopo_start == Vtopo_ndt) Vtopo_L = .false.

      Schm_itcn = keep_itcn

      call gtmg_stop ( 10 )

 1000 format( &
      /,'CONTROL OF DYNAMICAL STEP: (S/R DYNSTEP)', &
      /,'========================================'/)
 1005 format( &
      /3X,'##### Crank-Nicholson iterations: ===> PERFORMING',I3, &
          ' timestep(s) #####'/)
 1006 format( &
      /3X,'##### Crank-Nicholson iterations: ===> DONE... #####'/)
!
!     ---------------------------------------------------------------
!
      return
      end
