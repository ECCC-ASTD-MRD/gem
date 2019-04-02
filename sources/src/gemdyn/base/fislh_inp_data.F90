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

!**s/r inp_data_H  - Reads FST input files
! Author: Andre Plante 2017 based on inp_data.F90 from Michel Desgagne.

    subroutine fislh_inp_data ( F_u, F_v, F_t, F_q, F_zd, F_w, &
         F_topo,F_sls,Mminx,Mmaxx,Mminy,Mmaxy, Nk, &
         F_stag_L,F_trprefix_S,F_trsuffix_S,F_datev )

      use cstv
      use dyn_fisl_options
      use inp_options
      use gem_options
      use lam_options
      use gmm_geof
      use gmm_itf_mod
      use gmm_pw
      use HORgrid_options
      use inp_base, only: inp_get, inp_read, inp_3dpres, inp_hwnd
      use fislh_inp_base
      use inp_mod
      use glb_ld
      use lun
      use nest_blending
      use step_options
      use tr3d
      use ver
      use vertical_interpolation
      use vgrid_wb
      use vGrid_Descriptors

      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_trprefix_S, F_trsuffix_S, F_datev
      logical :: sfcTT_L, F_stag_L

      integer, intent(in) :: Mminx, Mmaxx, Mminy, Mmaxy, Nk
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk), intent(out) :: &
                                                       F_u, F_v, F_w, F_t, F_zd
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy), intent(in) :: F_sls
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy), intent(out) :: F_topo
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk+1),intent(out) :: F_q

#include <rmnlib_basics.hf>

      character(len=4) vname
      logical urt1_l, ut1_l
      integer  nka_gz, nka_tt, nka_hu, istat, err, kt, kh, kind, i, n
      integer, dimension (:), pointer :: ip1_dum, GZ_ip1_list, HU_ip1_list, &
                                                     ip1_list, ip1_w
      real step_current, lev
      real, dimension (:,:), pointer :: topo_temp, topols, ssq0, S_q, dummy
      real, dimension (:,:,:), pointer :: dstlev,srclev
      real, dimension (:,:,:), pointer :: ssqr, meqr, ttr, tvr, tv, &
                                          hur, gzr, trp
      real*8 :: diffd
      type(vgrid_descriptor) :: vgd_src, vgd_dst
      character(len=128) :: message_S
      logical :: found_L
      integer, parameter :: INP_OK = 0, INP_ERROR = -1
!
!-----------------------------------------------------------------------
!

      if (Lun_out.gt.0) write(lun_out,9000) trim(F_datev)

      if (.not.Lun_debug_L) istat= fstopc ('MSGLVL','SYSTEM',.false.)

      nullify (ssqr, meqr, ttr, tvr, tv, hur, gzr, dstlev, srclev, ssq0, dummy, S_q)
      nullify (GZ_ip1_list, HU_ip1_list, ip1_list)

      allocate(topo_temp(l_minx:l_maxx,l_miny:l_maxy), &
           topols(l_minx:l_maxx,l_miny:l_maxy), &
           dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk+1), &
           S_q(l_minx:l_maxx,l_miny:l_maxy))

      call inp_open ( F_datev, vgd_src )

!=============================
! Get destination vgrid object
!-----------------------------
      nullify (ip1_dum)
      istat= vgrid_wb_get ('ref-m', vgd_dst, ip1_dum )
      deallocate (ip1_dum); nullify (ip1_dum)

!========================
! Source surface pressure
!------------------------
      err= inp_read ( 'SFCPRES'    , 'Q', ssqr,    ip1_list, i    )

!=================
! Source Orography
!-----------------
      err= inp_read ( 'OROGRAPHY', 'Q', meqr, ip1_list, i )
      if ( trim(F_datev) .eq. trim(Step_runstrt_S) ) then
         if ( associated(meqr) ) then
            istat= gmm_get(gmmk_topo_low_s , topo_low )
            topo_low(1:l_ni,1:l_nj)= meqr(1:l_ni,1:l_nj,1)
         else
            Vtopo_L= .false.
         end if
      end if
      call difdatsd (diffd,Step_runstrt_S,F_datev)
      step_current = diffd*86400.d0 / Step_dt + Step_initial
      call var_topo2 ( F_topo, step_current, &
           l_minx,l_maxx,l_miny,l_maxy )
      if ( associated(meqr) .and. .not. Grd_yinyang_L) then
         if ( Lam_blendoro_L ) then
            topo_temp(1:l_ni,1:l_nj)= meqr(1:l_ni,1:l_nj,1)
            call nest_blend ( F_topo, topo_temp, l_minx,l_maxx, &
                 l_miny,l_maxy, 'M', level=G_nk+1 )
         end if
      end if

!=============
! Geopotential
!-------------
      err= inp_read ( 'GEOPOTENTIAL', 'Q', gzr, GZ_ip1_list, nka_gz )
      if (nka_gz.lt.1) &
           call gem_error (-1,'inp_data_H','Missing field: GZ - geopotential')
      call convip (GZ_ip1_list(nka_gz), lev, kind, -1, message_S, .false.)
      if( abs(lev-1.) > epsilon(lev))&
           call gem_error (-1,'inp_data_H','Missing field: GZ - at hyb or eta = 1.0')

!============================================================
! Get temperature and humidity, calculate virtual temperature
!------------------------------------------------------------
      err= inp_read ( 'TEMPERATURE', 'Q', ttr , ip1_list, nka_tt )
      if (nka_tt.lt.1) &
           call gem_error (-1,'inp_data_H','Missing field: TT - temperature')
      allocate(tvr(l_minx:l_maxx,l_miny:l_maxy,nka_tt))
      tvr=ttr
      err= inp_read ('TR/HU', 'Q', hur, HU_ip1_list, nka_hu)
      if (Inp_kind /= 105) then
         LOOP_TT_IP1: do kt=1, nka_tt
            found_L = .false.
            LOOP_HU_IP1: do kh=1, nka_hu
               if (ip1_list(kt) == HU_ip1_list(kh)) then
                  found_L = .true.
                  call mfottv2 ( tvr(l_minx,l_miny,kt),tvr(l_minx,l_miny,kt),&
                       hur(l_minx,l_miny,kh),l_minx,l_maxx       ,&
                       l_miny,l_maxy,1, 1,l_ni,1,l_nj, .true. )
                  exit
               end if
            end do LOOP_HU_IP1
            if(.not. found_L)then
               write(message_S,*)'Missing field: HU for ip1 = ',&
                    ip1_list(kt)
               call gem_error (-1,'inp_data_H',message_S)
            end if
         end do LOOP_TT_IP1
      end if

!==========================================================
! Interpolate temperature, virtual temperature and humidity
!----------------------------------------------------------
! Topo in m
      topo_temp = F_topo / grav_8
      topols    = F_sls  / grav_8
      allocate (srclev(l_minx:l_maxx,l_miny:l_maxy,nka_tt))
      allocate (tv(l_minx:l_maxx,l_miny:l_maxy,G_nk))
      if ( inp_match_heights(srclev, gzr, ip1_list, GZ_ip1_list, nka_tt, nka_gz) &
           == INP_ERROR) call gem_error (-1,'inp_data_H','See error above.')
      call inp_3dhgts ( vgd_dst, Ver_ip1%t, topo_temp, topols, dstlev, 1, G_nk)
      call vertint2 ( tv,dstlev,G_nk,tvr,srclev,nka_tt, &
           l_minx,l_maxx,l_miny,l_maxy   , &
           1,l_ni, 1,l_nj, varname='TT' , levtype='H')
      if(F_stag_L)then
         F_t = tv
      else
         call vertint2 ( F_t,dstlev,G_nk,ttr,srclev,nka_tt, &
              l_minx,l_maxx,l_miny,l_maxy   , &
              1,l_ni, 1,l_nj, varname='TT' , levtype='H')
      end if
      nullify (trp)
      istat= gmm_get (trim(F_trprefix_S)//'HU'//trim(F_trsuffix_S),trp)
      if ( nka_hu > 1 ) then
         err= inp_read ( 'TR/HU', 'Q', hur, HU_ip1_list, nka_hu )
         call vertint2 ( trp,dstlev,G_nk, hur,srclev,nka_hu         ,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni, 1,l_nj,&
                         inttype=Inp_vertintype_tracers_S )
      else
         trp = 0.
      end if

!open(unit=62,file='source.txt',status='unknown')
!open(unit=63,file='desti.txt',status='unknown')
!do k=1,nka_tt
!   write(62,*)srclev(1,1,k),tvr(1,1,k)
!end do
!do k=1,G_nk
!   write(63,*)dstlev(1,1,k),F_t(1,1,k)
!end do
!close(62)
!close(63)
!stop

!============================
! Get source surface pressure
!----------------------------
      call convip (ip1_list(nka_tt), lev, i,-1, vname, .false.)
      sfcTT_L = (i == 4) .or. ( (i /= 2) .and. (abs(lev-1.) <= 1.e-5) )
      allocate (ssq0(l_minx:l_maxx,l_miny:l_maxy))
      if (.not.associated(ssqr)) then
         ! check for input on pressure vertical coordinate
         call gem_error (-1,'inp_date_H','TODO: input in presure levels.')
      else
         S_q(:,:) = ssqr(:,:,1)
         deallocate (ssqr)
         call inp_3dpres ( vgd_src, ip1_list, S_q, dummy, srclev, 1,nka_tt )
         if ( associated(meqr) .and. sfcTT_L ) then
            if (lun_out.gt.0) &
                 write(lun_out,'(" PERFORMING surface pressure adjustment")')
! TODO review the following note and code
! Note at nka_tt may not be at surface for diag level.
! So we replace it to make sure it is the surface pressure.

            srclev(1:l_ni,1:l_nj,nka_tt)= S_q(1:l_ni,1:l_nj)
            call adj_ss2topo ( ssq0, F_topo, srclev, meqr, tvr  , &
                 l_minx,l_maxx,l_miny,l_maxy, nka_tt, &
                 1,l_ni,1,l_nj )
            deallocate (meqr) ; nullify (meqr)
         else
            if (lun_out.gt.0) &
                 write(lun_out,'(" NO surface pressure adjustment")')
            ssq0(1:l_ni,1:l_nj)= S_q(1:l_ni,1:l_nj)
         end if
      end if

      if (.not.associated(S_q)) &
           call gem_error ( -1, 'inp_data', &
           'Missing input data: surface pressure')

!=====================================================
! Compute pressure variable:             q=RT*ln(p/p0)
! integrating the hydrostatic relation:  dq/dz=-gT*/Tv
! Compute perturbation pressure variable: q'=q-q*=q+gz
!-----------------------------------------------------
! TODO trap error

      if( inp_comp_pres_from_vt_p0(F_q, ssq0, tv , topo_temp, topols, vgd_dst, &
           l_minx,l_maxx,l_miny,l_maxy,G_nk) == INP_ERROR) &
           call gem_error (-1,'inp_data_H inp_comp_pres_from_vt_p0','See error above.')
      deallocate(S_q)

!=============
! Read tracers
!-------------
      NTR_Tr3d_ntr= 0
      do n=1,Tr3d_ntr
         nullify (trp)
         vname= trim(Tr3d_name_S(n))
         istat= gmm_get (&
              trim(F_trprefix_S)//trim(vname)//trim(F_trsuffix_S),trp)
         if (trim(vname) /= 'HU') then
            err= inp_get_h ( 'TR/'//trim(vname),'Q', Ver_ip1%t,&
                 GZ_ip1_list, gzr, vgd_dst, topo_temp, topols, trp,&
                 l_minx,l_maxx,l_miny,l_maxy,G_nk ,&
                 F_inttype_S=Inp_vertintype_tracers_S )
            if (err == 0) then
               NTR_Tr3d_ntr= NTR_Tr3d_ntr + 1
               NTR_Tr3d_name_S(NTR_Tr3d_ntr) = trim(vname)
            end if
         end if
         trp= max(trp,Tr3d_vmin(n))
      end do

!=================================================
! Try to read WT1 && interpolate, set Inp_w_L flag
!-------------------------------------------------
      allocate (ip1_w(1:G_nk))
      ip1_w(1:G_nk)= Ver_ip1%t(1:G_nk)
      err= inp_get_h ('WT1',  'Q', ip1_w,&
           GZ_ip1_list, gzr, vgd_dst, topo_temp, topols, F_w ,&
           l_minx,l_maxx,l_miny,l_maxy,G_nk)
      deallocate (ip1_w) ; nullify(ip1_w)
      Inp_w_L= ( err == 0 )

!===================================================
! Try to read ZDT1 && interpolate, set Inp_zd_L flag
!---------------------------------------------------
      err= inp_get_h ('ZDT1', 'Q', Ver_ip1%t        ,&
           GZ_ip1_list, gzr, vgd_dst, topo_temp, topols, F_zd,&
           l_minx,l_maxx,l_miny,l_maxy,G_nk)
      Inp_zd_L= ( err == 0 )

!===========================
! Read and interpolate winds
!---------------------------
! TODO ajout code comme dans inp_data.F90 pour URT1 et cie?
      urt1_L = .false. ! TODO a enlever lorsque ajout code pour URT1 et cie
      ut1_L = .false.  ! TODO a enlever lorsque ajout code pour URT1 et cie
      if ((.not. urt1_L) .and. (.not. ut1_L)) &
           call inp_hwnd_H(F_u, F_v, gzr, GZ_ip1_list, vgd_dst, F_stag_L, &
           topo_temp, topols, l_minx,l_maxx,l_miny,l_maxy,G_nk )

!================
! Clean up memory
!----------------
      deallocate (gzr, ttr, tvr, tv, srclev, dstlev, topo_temp, topols)
      deallocate (GZ_ip1_list, ip1_list)
      if(associated(HU_ip1_list))deallocate(HU_ip1_list)
      if(associated(hur))deallocate(hur)

 9000 format(/,' TREATING INPUT DATA VALID AT: ',a,&
             /,' ===============================================')
!
!-----------------------------------------------------------------------
!
      return

    end
