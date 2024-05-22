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

!**s/r indata - Read and process the input data at
!               beginning of integration

      subroutine indata()
      use dynkernel_options
      use gem_options
      use glb_ld
      use gmm_geof
      use gmm_pw
      use gmm_vt1
      use inp_mod
      use inp_options
      use lun
      use metric
      use step_options
      use theo_options
      use tr3d
      use outp
      use omp_timing
      use mem_tracers
      implicit none

      logical :: synthetic_data_L
      integer :: dimens
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,1000)

      synthetic_data_L = (Ctrl_theoc_L .or. Ctrl_testcases_L)

      if (synthetic_data_L) then
         call synthetic_data ()
      else
         call gtmg_start ( 71, 'INITIAL_input', 2)

         call get_topo ()

         call inp_data (pw_uu_plus,pw_vv_plus,wt1,pw_tt_plus,qt1       ,&
                        zdt1,st1,trt1,fis0,orols,.false.,Step_runstrt_S,&
                        l_minx,l_maxx,l_miny,l_maxy,G_nk,Tr3d_ntr )

         dimens=(l_maxx-l_minx+1)*(l_maxy-l_miny+1)*G_nk
         call bitflip ( pw_uu_plus, pw_vv_plus, pw_tt_plus, &
                        perturb_nbits, perturb_npts, dimens )
         call gtmg_stop  ( 71 )
      end if

      call gemtime ( Lun_out, 'AFTER INITIAL INPUT', .false. )

      if (Schm_sleve_L) then
         call update_sls (orols,sls,l_minx,l_maxx,l_miny,l_maxy)
      endif
     
      if (Grd_yinyang_L) then
         call yyg_int_xch_scal (fis0 , 1, .false., 'CUBIC', .true.)
         call yyg_int_xch_scal (orols, 1, .false., 'CUBIC', .true.)
         call yyg_int_xch_scal (sls  , 1, .false., 'CUBIC', .true. )
         call yyg_xchng_all() !????? plus tard????
      else
         call rpn_comm_xch_halo(fis0, l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo(orols, l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo(sls, l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      end if

      if (Dynamics_hauteur_L) then
         call vertical_metric (GVM, fis0, sls, l_minx,l_maxx,l_miny,l_maxy)
         me_full (1-G_halox+1:l_ni+G_halox-1,1-G_haloy+1:l_nj+G_haloy-1) = &
             fis0(1-G_halox+1:l_ni+G_halox-1,1-G_haloy+1:l_nj+G_haloy-1) / grav_8
         me_large(1-G_halox+1:l_ni+G_halox-1,1-G_haloy+1:l_nj+G_haloy-1) = &
              sls(1-G_halox+1:l_ni+G_halox-1,1-G_haloy+1:l_nj+G_haloy-1) / grav_8
      endif
      
      trt0 = trt1

      if (.not. synthetic_data_L) then
         call tt2virt2 (tt1, .true., l_minx,l_maxx,l_miny,l_maxy, G_nk,&
                         1-G_halox,l_ni+G_halox, 1-G_haloy,l_nj+G_haloy)
         call hwnd_stag2 ( ut1,vt1, pw_uu_plus,pw_vv_plus   ,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk   ,&
                         1-G_halox*west ,l_niu+G_halox*east ,&
                         1-G_haloy*south,l_njv+G_haloy*north, .true. )
         call derivate_data (zdt1,wt1, ut1,vt1,tt1,st1,qt1,fis0,sls,GVM,&
                             l_minx,l_maxx,l_miny,l_maxy, G_nk         ,&
                             .not.Inp_zd_L, .not.Inp_w_L, .not.Inp_qt_L )
      endif

      if (.not. Grd_yinyang_L) call nest_init()

      call pressure ( pw_pm_plus,pw_pt_plus,pw_p0_plus,pw_log_pm,&
                      pw_log_pt, pw_pm_plus_8,pw_p0_plus_8      ,&
                      l_minx,l_maxx,l_miny,l_maxy,l_nk,1 )

      call pw_update_GW()
      call pw_init     ()
      
      udiag(:,:) = pw_uu_plus(:,:,G_nk)
      vdiag(:,:) = pw_vv_plus(:,:,G_nk)
      tdiag(:,:) = pw_tt_plus(:,:,G_nk)
      qdiag(:,:) = tracers_P(Tr3d_hu)%pntr(:,:,G_nk)

      call out_outdir()

      call iau_apply (0)
      
      call nudge_read (0,Lctl_step)

      if ( Ctrl_phyms_L ) call itf_phy_step (0,Lctl_step)

      if ( Dynamics_FISL_L ) call firstguess()

!      call nest_glbstat ((/'now','deb','fin'/),3)
!
!     ---------------------------------------------------------------
!
 1000 format(/,'TREATING INITIAL CONDITIONS  (S/R INDATA)',/,41('='))

      return
      end subroutine indata

      subroutine bitflip (u,v,t,nbits,npts,n)
      implicit none

      integer, intent(in) :: n,nbits,npts
      real, dimension(n), intent(inout) :: u, v, t

      integer stride,i
      integer, dimension(n) :: u_bits, v_bits, t_bits
!
! ---------------------------------------------------------------------
!
      u_bits  = transfer(u, u_bits)
      v_bits  = transfer(v, v_bits)
      t_bits  = transfer(t, t_bits)

      if (nbits < 1) return

      stride = min(max(1,npts),n)

      do i=1,n,stride
         u_bits(i) = xor(u_bits(i), nbits)
         v_bits(i) = xor(v_bits(i), nbits)
         t_bits(i) = xor(t_bits(i), nbits)
      end do

      u = transfer(u_bits, u)
      v = transfer(v_bits, v)
      t = transfer(t_bits, t)
!
! ---------------------------------------------------------------------
!
      return
      end subroutine bitflip
