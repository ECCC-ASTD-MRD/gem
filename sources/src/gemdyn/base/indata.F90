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
      use gmm_itf_mod
      use gmm_pw
      use gmm_vt1
      use inp_mod
      use inp_options
      use lun
      use step_options
      use theo_options
      use tr3d
      use gem_timing
      implicit none
#include <arch_specific.hf>

      integer :: k, istat, dimens, err
      real, dimension(l_minx:l_maxx,l_miny:l_maxy) :: topo_large_scale
      real, dimension(:,:,:), pointer, contiguous :: plus, minus
!
!     ---------------------------------------------------------------
!
      if (Lun_out > 0) write (Lun_out,1000)

      istat = gmm_get (gmmk_pw_uu_plus_s, pw_uu_plus)
      istat = gmm_get (gmmk_pw_vv_plus_s, pw_vv_plus)
      istat = gmm_get (gmmk_pw_tt_plus_s, pw_tt_plus)
      istat = gmm_get (gmmk_ut1_s ,ut1 )
      istat = gmm_get (gmmk_vt1_s ,vt1 )
      istat = gmm_get (gmmk_wt1_s ,wt1 )
      istat = gmm_get (gmmk_tt1_s ,tt1 )
      istat = gmm_get (gmmk_zdt1_s,zdt1)
      istat = gmm_get (gmmk_st1_s ,st1 )
      istat = gmm_get (gmmk_sls_s ,sls )
      istat = gmm_get (gmmk_fis0_s,fis0)
      istat = gmm_get (gmmk_qt1_s ,qt1 )

      zdt1=0. ; wt1=0. ; qt1= 0.

      if ( Ctrl_theoc_L ) then
         call theo_data ( pw_uu_plus, pw_vv_plus, pw_tt_plus, &
                          st1, qt1, fis0, sls )
      else if ( Ctrl_canonical_williamson_L ) then
         call init_bar ( ut1,vt1,wt1,tt1,zdt1,st1,qt1,fis0,&
                          l_minx,l_maxx,l_miny,l_maxy,G_nk,&
                         .true.,'TR/',':P',Step_runstrt_S )
      else if ( Ctrl_canonical_dcmip_L ) then
         call dcmip_init( ut1,vt1,wt1,tt1,zdt1,st1,qt1,fis0,&
                           l_minx,l_maxx,l_miny,l_maxy,G_nk,&
                         .true.,'TR/',':P')
      else
         call gemtime_start ( 71, 'INITIAL_input', 2)
         istat= gmm_get (gmmk_topo_low_s , topo_low )
         istat= gmm_get (gmmk_topo_high_s, topo_high)

         call get_topo ( topo_high, topo_large_scale, &
                          l_minx,l_maxx,l_miny,l_maxy, 1,l_ni,1,l_nj )

         call get_s_large_scale ( topo_large_scale, &
                                  l_minx,l_maxx,l_miny,l_maxy )

         topo_low(1:l_ni,1:l_nj) = topo_high(1:l_ni,1:l_nj)
         dimens=(l_maxx-l_minx+1)*(l_maxy-l_miny+1)*G_nk

         call inp_data ( pw_uu_plus, pw_vv_plus, wt1, pw_tt_plus,&
                       zdt1,st1,fis0,l_minx,l_maxx,l_miny,l_maxy,&
                         G_nk,.false. ,'TR/',':P',Step_runstrt_S )

         call bitflip ( pw_uu_plus, pw_vv_plus, pw_tt_plus, &
                        perturb_nbits, perturb_npts, dimens )
         call gemtime_stop  ( 71 )
      end if

      call gemtime ( Lun_out, 'AFTER INITIAL INPUT', .false. )

      call set_dync ( .true., err )

      if (Grd_yinyang_L) then
         call yyg_xchng (fis0,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj,&
                         1, .false., 'CUBIC', .true.)
         call yyg_xchng (sls,l_minx,l_maxx,l_miny,l_maxy,l_ni,l_nj, &
                         1, .false., 'CUBIC', .true. )
         ! Exchange halos in pilot region

         call yyg_xchng_all()
      else
         call rpn_comm_xch_halo(fis0, l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
         call rpn_comm_xch_halo(sls, l_minx,l_maxx,l_miny,l_maxy,&
            l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      end if

      do k=1, Tr3d_ntr
         nullify (plus, minus)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(k))//':M',minus)
         istat = gmm_get('TR/'//trim(Tr3d_name_S(k))//':P',plus )
         minus = plus
      end do

      if (.not. Ctrl_testcases_L) then
         call tt2virt (tt1, .true., l_minx,l_maxx,l_miny,l_maxy, G_nk)
         call hwnd_stag ( ut1,vt1, pw_uu_plus,pw_vv_plus,&
                          l_minx,l_maxx,l_miny,l_maxy,G_nk,.true. )
      end if

      if( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' .or. &
           trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H' ) then
         call fislh_metric ()
         call fislh_diag_zd_w( zdt1,wt1, ut1,vt1,tt1,st1,qt1    ,&
                               l_minx,l_maxx,l_miny,l_maxy, G_nk,&
                               .not.Inp_zd_L, .not.Inp_w_L )
      end if

      if (trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P') then
         call diag_zd_w ( zdt1,wt1, ut1,vt1,tt1,st1   ,&
                          l_minx,l_maxx,l_miny,l_maxy, G_nk,&
                          .not.Inp_zd_L, .not.Inp_w_L )
      end if

      if (.not. Grd_yinyang_L) call nest_init()

      call pw_update_GPW()
      call pw_init()

      call out_outdir()

      call iau_apply (0)

      if ( Ctrl_phyms_L ) call itf_phy_step (0,Lctl_step)

      if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_P' .or. &
           trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' ) then
         call frstgss()
      end if

      call glbstat2 ( fis0,'ME',"indata",l_minx,l_maxx,l_miny,l_maxy, &
                      1,1, 1,G_ni,1,G_nj,1,1 )
!
!     ---------------------------------------------------------------
!
 1000 format(/,'TREATING INITIAL CONDITIONS  (S/R INDATA)',/,41('='))
 1002 format(/,' FILE ',A,'_gfilemap.txt IS NOT AVAILABLE --CONTINUE--',/,/)

      return
      end

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
      end
