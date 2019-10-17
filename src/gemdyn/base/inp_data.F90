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

!**s/r inp_data  - Reads FST input files

      subroutine inp_data ( F_u, F_v, F_w, F_t, F_zd, F_s, F_topo   ,&
                            Mminx, Mmaxx, Mminy, Mmaxy, Nk, F_stag_L,&
                            F_trprefix_S, F_trsuffix_S, F_datev )
      use cstv
      use dyn_fisl_options
      use inp_options
      use gmm_itf_mod
      use inp_base, only: inp_get, inp_hwnd, inp_oro, inp_tv, &
                          inp_src_surface,inp_dst_surface   , &
                          inp_src_gz,inp_src_levels,inp_dst_levels
      use inp_mod
      use glb_ld
      use lun
      use tr3d
      use ver
      use vertical_interpolation, only: vertint2
      use vGrid_Descriptors
      implicit none
#include <arch_specific.hf>
#include <rmnlib_basics.hf>

      character(len=*), intent(in):: F_trprefix_S, F_trsuffix_S, F_datev
      logical, intent(in) :: F_stag_L
      integer, intent(in) :: Mminx, Mmaxx, Mminy, Mmaxy, Nk
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk),   intent(out) :: F_u, F_v, F_w, F_t, F_zd
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy),      intent(out) :: F_s, F_topo

      character(len=4) vname
      logical urt1_l, ut1_l
      integer n, nka, nka_tt, nka_hu, nka_gz, istat, err
      integer, dimension (:), pointer :: TT_ip1_list, HU_ip1_list,GZ_ip1_list,ip1_w
      real, dimension (:,:  ), pointer :: Sp0_q,Sp0_u,Sp0_v, &
                                          Dp0_q,Dp0_u,Dp0_v, &
                                          Dlsp0_q,Dlsp0_u,Dlsp0_v,meqr
      real, dimension (:,:  ), pointer :: dummy
      real, dimension (:,:,:), pointer :: dstlev,srclev
      real, dimension (:,:,:), pointer :: tv,ttr,hur,trp,gzr_q,gzr_u,gzr_v
!
!-----------------------------------------------------------------------
!
      if (Lun_out > 0) write(lun_out,9000) trim(F_datev)

      if (.not.Lun_debug_L) istat= fstopc ('MSGLVL','SYSTEM',RMN_OPT_SET)

      call inp_open ( F_datev )

      nullify ( meqr, tv, ttr, hur, TT_ip1_list, HU_ip1_list)

      call inp_oro ( F_topo, meqr, F_datev, l_minx,l_maxx,l_miny,l_maxy )

      call inp_tv (tv, ttr, hur, TT_ip1_list, HU_ip1_list, nka_tt,nka_hu)

      nullify ( Sp0_q,Sp0_u,Sp0_v,&
                Dp0_q,Dp0_u,Dp0_v,Dlsp0_q,Dlsp0_u,Dlsp0_v )
      nullify ( gzr_q,gzr_u,gzr_v,GZ_ip1_list )

      call inp_src_surface ( Sp0_q,Sp0_u,Sp0_v,gzr_q,gzr_u,gzr_v,&
                             F_topo,GZ_ip1_list,nka_gz,nka_tt )

      call inp_dst_surface ( F_s, Dp0_q, Dp0_u, Dp0_v,&
                             Dlsp0_q,Dlsp0_u,Dlsp0_v ,&
                             TT_ip1_list,Sp0_q,meqr,tv,F_topo        ,&
                             Schm_sleve_L,l_minx,l_maxx,l_miny,l_maxy,&
                             nka_tt,1,l_ni,1,l_nj )

      nullify (srclev,dstlev)
      allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,max(nka_tt,nka_hu)),&
                 dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk))

      nullify(dummy)
      call inp_src_levels ( srclev, nka, TT_ip1_list, Inp_vgd_src,Sp0_q,&
                    dummy,gzr_q,GZ_ip1_list, l_minx,l_maxx,l_miny,l_maxy)
      call inp_dst_levels (dstlev, Ver_vgdobj, Ver_ip1%t, Dp0_q, Dlsp0_q)

      if (F_stag_L) then
         call vertint2 ( F_t,dstlev,G_nk, tv,srclev,nka, &
                         l_minx,l_maxx,l_miny,l_maxy   , &
                         1,l_ni, 1,l_nj, varname='TT'  , &
                         levtype=Inp_levtype_S )
      else
         call vertint2 ( F_t,dstlev,G_nk,ttr,srclev,nka, &
                         l_minx,l_maxx,l_miny,l_maxy   , &
                         1,l_ni, 1,l_nj, varname='TT'  , &
                         levtype=Inp_levtype_S )
      end if

      nullify (trp)
      istat= gmm_get (trim(F_trprefix_S)//'HU'//trim(F_trsuffix_S),trp)

      call vertint2 ( trp,dstlev,G_nk, hur,srclev,nka            ,&
                      l_minx,l_maxx,l_miny,l_maxy, 1,l_ni, 1,l_nj,&
                inttype=Inp_vertintype_tracers_S, levtype=Inp_levtype_S )

      deallocate (tv,ttr,hur,srclev,dstlev)
      nullify (tv,ttr,hur,srclev,dstlev)
      if (associated(TT_ip1_list)) deallocate (TT_ip1_list)
      if (associated(HU_ip1_list)) deallocate (HU_ip1_list)

      NTR_Tr3d_ntr= 0
      do n=1,Tr3d_ntr
         nullify (trp)
         vname= trim(Tr3d_name_S(n))
         istat= gmm_get (&
               trim(F_trprefix_S)//trim(vname)//trim(F_trsuffix_S),trp)
         if (trim(vname) /= 'HU') then
            err = inp_get ( 'TR/'//trim(vname),'Q', Ver_ip1%t, Sp0_q,&
                            Dp0_q, Dlsp0_q,gzr_q,GZ_ip1_list,trp    ,&
                            l_minx,l_maxx,l_miny,l_maxy,G_nk        ,&
                            F_inttype_S=Inp_vertintype_tracers_S )
            if (err == 0) then
               NTR_Tr3d_ntr= NTR_Tr3d_ntr + 1
               NTR_Tr3d_name_S(NTR_Tr3d_ntr) = trim(vname)
            end if
         end if
         trp= max(min(trp,Tr3d_vmax(n)),Tr3d_vmin(n))
      end do

      allocate (ip1_w(1:G_nk))
      ip1_w(1:G_nk)= Ver_ip1%t(1:G_nk)
      err = inp_get ( 'WT1',  'Q', ip1_w,&
                   Sp0_q, Dp0_q, Dlsp0_q,gzr_q,GZ_ip1_list,F_w,&
                   l_minx,l_maxx,l_miny,l_maxy,G_nk,F_quiet_L=.true.)
      deallocate (ip1_w) ; nullify(ip1_w)
      Inp_w_L= ( err == 0 )

      err = inp_get ('ZDT1', 'Q', Ver_ip1%t,&
                   Sp0_q, Dp0_q, Dlsp0_q,gzr_q,GZ_ip1_list,F_zd,&
                   l_minx,l_maxx,l_miny,l_maxy,G_nk,F_quiet_L=.true.)
      Inp_zd_L= ( err == 0 )

      ut1_L= .false. ; urt1_L= .false. ; err = -1

      if (F_stag_L) then
         err = inp_get ( 'URT1', 'U', Ver_ip1%m               ,&
                   Sp0_u, Dp0_u, Dlsp0_u,gzr_u,GZ_ip1_list,F_u,&
                   l_minx,l_maxx,l_miny,l_maxy,G_nk,F_quiet_L=.true. )
         if ( err == 0 ) then
            err = inp_get ( 'VRT1', 'V', Ver_ip1%m            ,&
                   Sp0_v, Dp0_v, Dlsp0_v,gzr_v,GZ_ip1_list,F_v,&
                   l_minx,l_maxx,l_miny,l_maxy,G_nk,F_quiet_L=.true. )
         end if
         urt1_L= ( err == 0 )

         if (.not. urt1_L) then
            err = inp_get ( 'UT1', 'U', Ver_ip1%m              ,&
                   Sp0_u, Dp0_u, Dlsp0_u,gzr_u,GZ_ip1_list,F_u,&
                   l_minx,l_maxx,l_miny,l_maxy,G_nk,F_quiet_L=.true. )
            if ( err == 0 ) then
               err = inp_get ( 'VT1', 'V', Ver_ip1%m           ,&
                   Sp0_v, Dp0_v, Dlsp0_v,gzr_v,GZ_ip1_list,F_v,&
                   l_minx,l_maxx,l_miny,l_maxy,G_nk,F_quiet_L=.true. )
            end if
            ut1_L= ( err == 0 )
         ! Remove the .and. part of this test by 2021
            if ( ut1_L .and. Inp_ut1_is_urt1 == -1 ) then
              call image_to_real_winds ( F_u,F_v, l_minx,l_maxx,&
                                         l_miny,l_maxy, G_nk )
            end if
      end if
      end if

      if ((.not. urt1_L) .and. (.not. ut1_L)) then
         call inp_hwnd ( F_u,F_v, F_stag_L, Sp0_q,Sp0_u,Sp0_v     ,&
                         Dp0_q,Dp0_u,Dp0_v,Dlsp0_q,Dlsp0_u,Dlsp0_v,&
                         gzr_q,gzr_u,gzr_v,GZ_ip1_list            ,&
                         l_minx,l_maxx,l_miny,l_maxy,G_nk )
      end if

      deallocate (Sp0_q)
      if (associated(Sp0_u)) deallocate (Sp0_u)
      if (associated(Sp0_v)) deallocate (Sp0_v)
      nullify(Sp0_q,Sp0_u,Sp0_v)
      deallocate (Dp0_q,Dp0_u,Dp0_v) ; nullify(Dp0_q,Dp0_u,Dp0_v)
      if (Schm_sleve_L) then
         deallocate (Dlsp0_u,Dlsp0_v)
         nullify    (Dlsp0_q,Dlsp0_u,Dlsp0_v)
      endif

      call inp_close ()

      istat = fstopc ('MSGLVL','WARNIN',RMN_OPT_SET)

 9000 format(/,' TREATING INPUT DATA VALID AT: ',a,&
             /,' ===============================================')
!
!-----------------------------------------------------------------------
!
      return
      end
