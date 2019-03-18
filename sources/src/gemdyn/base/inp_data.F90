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

      subroutine inp_data ( F_u, F_v, F_w, F_t, F_zd, F_s, F_q    , &
                            F_topo, Mminx, Mmaxx, Mminy, Mmaxy, Nk, &
                      F_stag_L, F_trprefix_S, F_trsuffix_S, F_datev )
      use cstv
      use dynkernel_options
      use dyn_fisl_options
      use inp_options
      use gem_options
      use lam_options
      use gmm_geof
      use gmm_itf_mod
      use gmm_pw
      use HORgrid_options
      use inp_base, only: inp_get, inp_read, inp_3dpres, inp_hwnd
      use inp_mod
      use glb_ld
      use lun
      use nest_blending, only: nest_blend
      use step_options
      use tr3d
      use ver
      use vertical_interpolation, only: vertint2
      use vgrid_wb, only: vgrid_wb_get
      use vGrid_Descriptors

      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in):: F_trprefix_S, F_trsuffix_S, F_datev
      logical, intent(in) :: F_stag_L
      integer, intent(in) :: Mminx, Mmaxx, Mminy, Mmaxy, Nk
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk),   intent(out) :: F_u, F_v, F_w, F_t, F_zd
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy),      intent(out) :: F_s, F_topo
      real, dimension(Mminx:Mmaxx,Mminy:Mmaxy,Nk+1), intent(out) :: F_q

#include <rmnlib_basics.hf>

      character(len=4) vname
      logical urt1_l, ut1_l, sfcTT_L
      integer i,j,k,n, nka, istat, err, kind,kt,kh
      integer nka_gz, nka_tt, nka_hu
      integer, dimension (:), pointer :: ip1_list, HU_ip1_list, &
                                         ip1_dum, ip1_w
      real topo_temp(l_minx:l_maxx,l_miny:l_maxy),step_current,lev
      real, dimension (:    ), allocatable  :: rna
      real, dimension (:,:  ), pointer :: S_q,S_u,S_v
      real, dimension (:,:  ), pointer :: ssq0,ssu0,ssv0
      real, dimension (:,:  ), pointer :: p0lsq0,p0lsu0,p0lsv0,dummy
      real, dimension (:,:,:), pointer :: dstlev,srclev
      real, dimension (:,:,:), pointer :: ssqr,ssur,ssvr,meqr,&
                                          tv,ttr,hur,gzr,trp
      logical :: using_qt1
      real*8 diffd
      type(vgrid_descriptor) :: vgd_src, vgd_dst
!
!-----------------------------------------------------------------------
!
      using_qt1 = ( .not. Dynamics_hydro_L ) .or. (trim(Dynamics_Kernel_S) == 'DYNAMICS_EXPO_H')

      if (Lun_out > 0) write(lun_out,9000) trim(F_datev)

      if (.not.Lun_debug_L) istat= fstopc ('MSGLVL','SYSTEM',.false.)

      call inp_open ( F_datev, vgd_src )

      nullify (meqr, ip1_list, HU_ip1_list)
      err = inp_read ( 'OROGRAPHY', 'Q', meqr, ip1_list, nka )
      if ( associated(ip1_list) ) then
         deallocate (ip1_list) ; nullify (ip1_list)
      end if

      if ( trim(F_datev) == trim(Step_runstrt_S) ) then
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

      nullify (ssqr,ssur,ssvr,ttr,hur,p0lsq0,p0lsu0,p0lsv0,dummy)

      err = inp_read ( 'SFCPRES'    , 'Q', ssqr,    ip1_list, nka    )

      if (associated(ssqr)) then
         err = inp_read ( 'SFCPRES'    , 'U', ssur,    ip1_list, nka    )
         if (nka < 1) then
            call gem_error (-1,'inp_data','Missing field: SFCPRES on hgrid U')
         end if

         err = inp_read ( 'SFCPRES'    , 'V', ssvr,    ip1_list, nka    )
         if (nka < 1) then
            call gem_error (-1,'inp_data','Missing field: SFCPRES on hgrid V')
         end if
      end if

      err = inp_read ( 'TR/HU'      , 'Q', hur , HU_ip1_list, nka_hu )
      if (nka_hu < 1) then
         call gem_error (-1,'inp_data','Missing field: TR/HU - humidity TR/HU')
      end if

      err = inp_read ( 'TEMPERATURE', 'Q', ttr ,    ip1_list, nka_tt )
      if (nka_tt < 1) then
         call gem_error (-1,'inp_data','Missing field: TT - temperature TT')
      end if

      nka= nka_tt
      allocate (tv(l_minx:l_maxx,l_miny:l_maxy,nka))
      tv=ttr

      if (Inp_kind /= 105) then
         do kt=1, nka_tt
inner:      do kh=1, nka_hu
               if (ip1_list(kt) == HU_ip1_list(kh)) then
                  call mfottv2 ( tv(l_minx,l_miny,kt),tv(l_minx,l_miny,kt),&
                                 hur(l_minx,l_miny,kh),l_minx,l_maxx       ,&
                                 l_miny,l_maxy,1, 1,l_ni,1,l_nj, .true. )
                  exit inner
               end if
            end do inner
         end do
      end if

      call convip (ip1_list(nka), lev, i,-1, vname, .false.)
      sfcTT_L = (i == 4) .or. ( (i /= 2) .and. (abs(lev-1.) <= 1.e-5) )

      allocate ( srclev(l_minx:l_maxx,l_miny:l_maxy,max(nka,nka_hu)) ,&
                 dstlev(l_minx:l_maxx,l_miny:l_maxy,G_nk),&
                 ssq0  (l_minx:l_maxx,l_miny:l_maxy)     ,&
                 ssu0  (l_minx:l_maxx,l_miny:l_maxy)     ,&
                 ssv0  (l_minx:l_maxx,l_miny:l_maxy) )
      nullify (S_q,S_u,S_v)

      if (.not.associated(ssqr)) then
         ! check for input on pressure vertical coordinate
         if (Inp_kind == 2) then
            nullify(gzr)
            err = inp_read ( 'GEOPOTENTIAL' , 'Q', gzr, ip1_list, nka_gz )
            if (nka_gz == nka) then
               allocate (rna(nka))
               do k=1,nka
                  call convip(ip1_list(k),rna(k),kind,-1,' ',.false.)
               end do
               call gz2p02 ( ssq0, gzr, F_topo, rna     ,&
                             l_minx,l_maxx,l_miny,l_maxy,&
                             nka,1,l_ni,1,l_nj )
               allocate (S_q(l_minx:l_maxx,l_miny:l_maxy))
               allocate (S_u(l_minx:l_maxx,l_miny:l_maxy))
               allocate (S_v(l_minx:l_maxx,l_miny:l_maxy))
               S_q(:,1:nka) = rna(nka)
               S_u(:,1:nka) = rna(nka)
               S_v(:,1:nka) = rna(nka)
               deallocate (rna,gzr) ; nullify (gzr)
            end if
         end if
      else
         allocate (S_q(l_minx:l_maxx,l_miny:l_maxy))
         allocate (S_u(l_minx:l_maxx,l_miny:l_maxy))
         allocate (S_v(l_minx:l_maxx,l_miny:l_maxy))
         S_q(:,:) = ssqr(:,:,1)
         S_u(:,:) = ssur(:,:,1)
         S_v(:,:) = ssvr(:,:,1)
         deallocate (ssqr,ssur,ssvr) ; nullify(ssqr,ssur,ssvr)

         call inp_3dpres ( vgd_src, ip1_list, S_q, dummy, srclev, 1,nka )

         if ( associated(meqr) .and. sfcTT_L ) then
            if (lun_out > 0) then
               write(lun_out,'(" PERFORMING surface pressure adjustment")')
            end if
            srclev(1:l_ni,1:l_nj,nka)= S_q(1:l_ni,1:l_nj)
            call adj_ss2topo ( ssq0, F_topo, srclev, meqr, tv  , &
                               l_minx, l_maxx, l_miny, l_maxy, nka, &
                               1, l_ni, 1, l_nj )
            deallocate (meqr) ; nullify (meqr)
         else
            if (lun_out > 0) then
               write(lun_out,'(" NO surface pressure adjustment")')
            end if
            ssq0(1:l_ni,1:l_nj)= S_q(1:l_ni,1:l_nj)
         end if
      end if

      if (.not.associated(S_q)) then
         call gem_error ( -1, 'inp_data', &
                          'Missing input data: surface pressure')
      end if

      call rpn_comm_xch_halo ( ssq0, l_minx,l_maxx,l_miny,l_maxy,&
           l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )

      do j=1,l_nj
      do i=1,l_ni-east
         ssu0(i,j)= (ssq0(i,j)+ssq0(i+1,j  ))*.5d0
      end do
      end do
      if (l_east) ssu0(l_ni,1:l_nj)= ssq0(l_ni,1:l_nj)

      do j=1,l_nj-north
      do i=1,l_ni
         ssv0(i,j)= (ssq0(i,j)+ssq0(i  ,j+1))*.5d0
      end do
      end do
      if (l_north) ssv0(1:l_ni,l_nj)= ssq0(1:l_ni,l_nj)

      F_s(1:l_ni,1:l_nj) = log(ssq0(1:l_ni,1:l_nj)/Cstv_pref_8)

      if ( Schm_sleve_L ) then
         allocate ( p0lsu0(l_minx:l_maxx,l_miny:l_maxy),&
                    p0lsv0(l_minx:l_maxx,l_miny:l_maxy) )
         istat = gmm_get (gmmk_pw_p0_ls_s, p0lsq0 )
         call rpn_comm_xch_halo ( p0lsq0, l_minx,l_maxx,l_miny,l_maxy,&
              l_ni,l_nj,1,G_halox,G_haloy,G_periodx,G_periody,l_ni,0 )
         do j=1,l_nj
         do i=1,l_ni-east
            p0lsu0(i,j)= (p0lsq0(i,j)+p0lsq0(i+1,j  ))*.5d0
         end do
         end do
         if (l_east) p0lsu0(l_ni,1:l_nj)= p0lsq0(l_ni,1:l_nj)

         do j=1,l_nj-north
         do i=1,l_ni
            p0lsv0(i,j)= (p0lsq0(i,j)+p0lsq0(i  ,j+1))*.5d0
         end do
         end do
         if (l_north) p0lsv0(1:l_ni,l_nj)= p0lsq0(1:l_ni,l_nj)
      end if

      nullify (ip1_dum)
      istat= vgrid_wb_get ('ref-m', vgd_dst, ip1_dum )
      deallocate (ip1_dum); nullify (ip1_dum)

      call inp_3dpres ( vgd_src,  ip1_list,  S_q,  dummy, srclev, &
                        1, nka, F_inlog_S='in_log')
      call inp_3dpres ( vgd_dst, Ver_ip1%t, ssq0, p0lsq0, dstlev, &
                        1,G_nk, F_inlog_S='in_log')

      if (F_stag_L) then
         call vertint2 ( F_t,dstlev,G_nk, tv,srclev,nka, &
                         l_minx,l_maxx,l_miny,l_maxy   , &
                         1,l_ni, 1,l_nj, varname='TT' )
      else
         call vertint2 ( F_t,dstlev,G_nk,ttr,srclev,nka, &
                         l_minx,l_maxx,l_miny,l_maxy   , &
                         1,l_ni, 1,l_nj, varname='TT' )
      end if

      nullify (trp)
      istat= gmm_get (trim(F_trprefix_S)//'HU'//trim(F_trsuffix_S),trp)
!      if ( nka_hu > 1 ) then
!         err = inp_read ( 'TR/HU'   , 'Q', hur , HU_ip1_list, nka_hu )
         call inp_3dpres ( vgd_src, HU_ip1_list, S_q, dummy, srclev, 1, &
                           nka_hu, F_inlog_S='in_log')
         call vertint2 ( trp,dstlev,G_nk, hur,srclev,nka_hu         ,&
                         l_minx,l_maxx,l_miny,l_maxy, 1,l_ni, 1,l_nj,&
                         inttype=Inp_vertintype_tracers_S )
!      else
!         trp = 0.
!      end if

      deallocate (tv,ttr,hur,srclev,dstlev) ; nullify (tv,ttr,hur,srclev,dstlev)
      if (associated(ip1_list   )) deallocate (ip1_list   )
      if (associated(HU_ip1_list)) deallocate (HU_ip1_list)

      NTR_Tr3d_ntr= 0
      do n=1,Tr3d_ntr
         nullify (trp)
         vname= trim(Tr3d_name_S(n))
         istat= gmm_get (&
               trim(F_trprefix_S)//trim(vname)//trim(F_trsuffix_S),trp)
         if (trim(vname) /= 'HU') then
            err = inp_get ( 'TR/'//trim(vname),'Q', Ver_ip1%t,&
                           vgd_src, vgd_dst, S_q, ssq0, p0lsq0, trp,&
                           l_minx,l_maxx,l_miny,l_maxy,G_nk ,&
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
      err = inp_get ('WT1',  'Q', ip1_w            ,&
                    vgd_src,vgd_dst,S_q,ssq0,p0lsq0,F_w ,&
                    l_minx,l_maxx,l_miny,l_maxy,G_nk)
      deallocate (ip1_w) ; nullify(ip1_w)
      Inp_w_L= ( err == 0 )

      err = inp_get ('ZDT1', 'Q', Ver_ip1%t        ,&
                    vgd_src,vgd_dst,S_q,ssq0,p0lsq0,F_zd,&
                    l_minx,l_maxx,l_miny,l_maxy,G_nk)
      Inp_zd_L= ( err == 0 )

      if (using_qt1) then
         allocate (ip1_w(1:G_nk+1))
         ip1_w(1:G_nk+1)=Ver_ip1%m(1:G_nk+1)
         err = inp_get ( 'QT1', 'Q', ip1_w,&
                        vgd_src,vgd_dst,S_q,ssq0,p0lsq0,F_q ,&
                        l_minx,l_maxx,l_miny,l_maxy,G_nk+1 )
         deallocate (ip1_w) ; nullify(ip1_w)
      end if

      ut1_L= .false. ; urt1_L= .false.

      if (F_stag_L) then
      err = inp_get ( 'URT1', 'U', Ver_ip1%m       ,&
                     vgd_src,vgd_dst,S_u,ssu0,p0lsu0,F_u,&
                     l_minx,l_maxx,l_miny,l_maxy,G_nk )
      if ( err == 0 ) then
         err = inp_get ( 'VRT1', 'V', Ver_ip1%m       ,&
                        vgd_src,vgd_dst,S_v,ssv0,p0lsv0,F_v,&
                        l_minx,l_maxx,l_miny,l_maxy,G_nk )
      end if
      urt1_L= ( err == 0 )

      if (.not. urt1_L) then
         err = inp_get ( 'UT1', 'U', Ver_ip1%m        ,&
                        vgd_src,vgd_dst,S_u,ssu0,p0lsu0,F_u,&
                        l_minx,l_maxx,l_miny,l_maxy,G_nk )
         if ( err == 0 ) then
            err = inp_get ( 'VT1', 'V', Ver_ip1%m        ,&
                           vgd_src,vgd_dst,S_v,ssv0,p0lsv0,F_v,&
                           l_minx,l_maxx,l_miny,l_maxy,G_nk )
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
         call inp_hwnd ( F_u,F_v, vgd_src,vgd_dst, F_stag_L,&
                          S_q,S_u,S_v, ssq0,ssu0,ssv0       ,&
                          p0lsq0,p0lsu0,p0lsv0              ,&
                          l_minx,l_maxx,l_miny,l_maxy,G_nk )
      end if

      deallocate (S_q,S_u,S_v,ssq0,ssu0,ssv0)
      nullify    (S_q,S_u,S_v,ssq0,ssu0,ssv0)
      if ( Schm_sleve_L ) then
         deallocate (p0lsu0,p0lsv0)
         nullify    (p0lsq0,p0lsu0,p0lsv0)
      end if

      err = vgd_free(vgd_dst)

      call inp_close ()

      istat = fstopc ('MSGLVL','INFORM',.false.)

 9000 format(/,' TREATING INPUT DATA VALID AT: ',a,&
             /,' ===============================================')
!
!-----------------------------------------------------------------------
!
      return
      end
