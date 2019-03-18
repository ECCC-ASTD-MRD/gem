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
                            orhsw, rhsw, Minx,Maxx,Miny,Maxy,&
                            Nk, F_tracers_L )
      use adz_mem
      use dynkernel_options
      use gem_timing
      use gmm_vt0
      use gmm_vt1
      use gmm_pw
      use gmm_itf_mod
      use gmm_tracers
      use tr3d
      implicit none
#include <arch_specific.hf>

      logical,intent(in) :: F_tracers_L
      integer,intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(in)  :: &
                             orhsu, orhsv, orhsc, orhst, orhsf, orhsw
      real, dimension(Minx:Maxx,Miny:Maxy,NK),intent(out) :: &
                             rhsu, rhsv, rhsc,  rhst, rhsf, rhsw

      logical :: mono_L
      integer :: n, err
      real, pointer, contiguous, dimension (:,:,:) :: src,dst
!
!     ---------------------------------------------------------------
!
      call gemtime_start (30, 'ADZ_TRAJEC', 21)
      call adz_traject (pw_uu_moins, pw_vv_moins, zdt1, &
                           ut0     ,    vt0     , zdt0, &
                       l_minx,l_maxx,l_miny,l_maxy,l_nk)
      call gemtime_stop (30)

      call gemtime_start (31, 'ADZ_INTP_RH', 21)

      call adz_cubic ( rhsu, orhsu, Adz_pxyzmu                    ,&
                       l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                    Adz_i0u,Adz_inu, Adz_j0,Adz_jn, Adz_k0, 'm', .false.)

      call adz_cubic ( rhsv, orhsv, Adz_pxyzmv                    ,&
                       l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                    Adz_i0,Adz_in, Adz_j0v,Adz_jnv, Adz_k0, 'm', .false.)

      call adz_cubic ( rhsc, orhsc, Adz_pxyzm                    ,&
                       l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                    Adz_i0,Adz_in, Adz_j0,Adz_jn, Adz_k0, 'm', .false.)

      call adz_cubic ( rhst, orhst, Adz_pxyzt                     ,&
                       l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                    Adz_i0,Adz_in, Adz_j0,Adz_jn, Adz_k0t, 't', .false.)

      call adz_cubic ( rhsf, orhsf, Adz_pxyzt                     ,&
                       l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                    Adz_i0,Adz_in, Adz_j0,Adz_jn, Adz_k0t, 't', .false.)

      if(.not.Dynamics_hydro_L) then
         call adz_cubic ( rhsw, orhsw, Adz_pxyzt                     ,&
                          l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                    Adz_i0,Adz_in, Adz_j0,Adz_jn, Adz_k0t, 't', .false.)
      end if

      if (F_tracers_L) then
      do n=1, Tr3d_ntr
         if ( (Tr3d_mass(n) == 0) .and. (Tr3d_mono(n) < 2) ) then
            err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':P' ,src)
            err= gmm_get('TR/'//trim(Tr3d_name_S(n))//':M' ,dst)
            mono_L = (Tr3d_mono(n)==1)
            call adz_cubic (dst, src, Adz_pxyzt                   ,&
                       l_ni,l_nj,l_nk, l_minx,l_maxx,l_miny,l_maxy,&
                       Adz_i0,Adz_in, Adz_j0,Adz_jn, Adz_k0, 't', mono_L)
         end if
      end do
      end if

      call gemtime_stop (31)
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_main
