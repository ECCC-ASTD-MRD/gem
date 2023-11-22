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

!**s/r tracers - Establishes final list of tracers for GEMDM

      subroutine tracers ()
      use phy_itf, only: PHY_MAXNAMELENGTH,phymeta,phy_getmeta
      use adz_mem
      use geomh
      use gem_options
      use ctrl
      use lun
      use rstr
      use tr3d
      use clib_itf_mod
      use tracers_attributes_mod, only: tracers_attributes
      implicit none
#include <arch_specific.hf>

!        Tr3d_name_S acquires the list of tracers from itf_phy_inikey
!        In this subroutine, it will acquire what is introduced from
!        Tr3d_username_S (constructed from the entry routine)
!        ie: QC can be filled with either QC or QCT1 or QCT0 from the
!            given analysis (accessed by E_tr3dname_S of gement nml)
!            but this is converted via BMF/BCS/3DF as QC
!        Tracers requested under auto cascade requires only the first
!        2 letters

      character(len=512) :: varname,attributes
      character(len=PHY_MAXNAMELENGTH) :: varname_S,prefix_S, &
                                          basename_S,time_S,ext_S
      integer i,j,ind,n,wload,hzd,monot,massc,dejala,istat,nmeta
      real vmin,vmax
      character(len=12) intp_S
      type(phymeta)     , dimension(:), pointer :: pmeta
!
!     __________________________________________________________________
!
      nmeta= 0 ; nullify(pmeta)
      if ( Ctrl_phyms_L ) then
         nmeta = phy_getmeta(pmeta,' ',F_npath='V',F_bpath='D',F_quiet=.true.)
      end if

      do i=1,nmeta
         varname_S = pmeta(i)%vname
         istat = clib_toupper(varname_S)
         if (varname_S(1:3) /= 'TR/') cycle
         call gmmx_name_parts(varname_S,prefix_S,basename_S,time_S,ext_S)
         if (time_S/=':P') cycle
         dejala=0
         do j=1,Tr3d_ntr
            if (trim(Tr3d_name_S(j))==basename_S) dejala=j
         end do
         if (dejala==0) then
            Tr3d_ntr = Tr3d_ntr + 1
            dejala   = Tr3d_ntr
            Tr3d_name_S(dejala)= basename_S(1:4)
         end if
         Tr3d_hzd (dejala)= pmeta(i)%hzd   ; Tr3d_wload(dejala)= pmeta(i)%wload
         Tr3d_mono(dejala)= pmeta(i)%monot ; Tr3d_mass (dejala)= pmeta(i)%massc
         Tr3d_vmin(dejala)= pmeta(i)%vmin  ; Tr3d_vmax (dejala)= pmeta(i)%vmax
         Tr3d_intp(dejala)= 'TRICUB'
      end do

      do i=1,MAXTR3D
         if (Tr3d_list_s(i)=='') exit
         ind= index(Tr3d_list_s(i),",")
         if (ind == 0) then
            call low2up(Tr3d_list_s(i), varname)
            attributes = ''
         else
            call low2up(Tr3d_list_s(i)(1:ind-1),varname   )
            call low2up(Tr3d_list_s(i)(ind+1: ),attributes)
         end if
         istat = tracers_attributes(attributes,wload,hzd,monot,massc,vmin,vmax,intp_S)
         dejala=0
         do j=1,Tr3d_ntr
            if (trim(Tr3d_name_S(j))==trim(varname)) dejala=j
         end do
         if (dejala==0) then
            Tr3d_ntr = Tr3d_ntr + 1
            dejala   = Tr3d_ntr
            Tr3d_name_S(dejala)= trim(varname)
         end if
         Tr3d_hzd (dejala)= (hzd>0) ; Tr3d_wload(dejala)= (wload>0)
         Tr3d_mono(dejala)= monot   ; Tr3d_mass (dejala)= massc
         Tr3d_vmin(dejala)= vmin    ; Tr3d_vmax (dejala)= vmax
         Tr3d_intp(dejala)= intp_S
      end do

      dejala=0
      do j=1,Tr3d_ntr
         if (trim(Tr3d_name_S(j))=='HU') dejala=j
      end do
      if (dejala==0) then
         Tr3d_ntr = Tr3d_ntr + 1
         Tr3d_name_S(tr3d_ntr)(1:4) = 'HU  '
         Tr3d_hzd   (Tr3d_ntr)= .false. ; Tr3d_wload(Tr3d_ntr)= .false.
         Tr3d_mono  (Tr3d_ntr)= 0       ; Tr3d_mass (Tr3d_ntr)= 0
         Tr3d_vmin  (Tr3d_ntr)= 0.      ; Tr3d_vmax (Tr3d_ntr)= huge(1.)
         Tr3d_intp  (Tr3d_ntr)= 'TRICUB'
      else
         Tr3d_wload(dejala)= .false.
      end if

      Tr3d_ntrTRICUB_NT= 0 ; Tr3d_ntrTRICUB_WP= 0
      Tr3d_ntrBICHQV_NT= 0 ; Tr3d_ntrBICHQV_WP= 0
      allocate ( Tr_3CNT(Tr3d_ntr), Tr_3CWP(Tr3d_ntr), &
                 Tr_BQNT(Tr3d_ntr), Tr_BQWP(Tr3d_ntr) )

      if (Tr3d_ntr > MAXTR3D) &
      call gem_error (-1, 'TRACERS', 'Tr3d_ntr > MAXTR3D: WE STOP')

      do i=1, Tr3d_ntr
         if (Tr3d_intp(i) == 'TRICUB') then
            if ( (Tr3d_mass(i) < 1) .and. (Tr3d_mono(i) < 1)) then
               Tr3d_ntrTRICUB_NT= Tr3d_ntrTRICUB_NT+1
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%name = Tr3d_name_S(i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%intp = Tr3d_intp  (i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%wload= Tr3d_wload (i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%hzd  = Tr3d_hzd   (i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%mono = Tr3d_mono  (i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%mass = Tr3d_mass  (i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%vmin = Tr3d_vmin  (i)
               Tr_3CNT(Tr3d_ntrTRICUB_NT)%vmax = Tr3d_vmax  (i)
            else
               Tr3d_ntrTRICUB_WP= Tr3d_ntrTRICUB_WP+1
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%name = Tr3d_name_S(i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%intp = Tr3d_intp  (i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%wload= Tr3d_wload (i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%hzd  = Tr3d_hzd   (i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%mono = Tr3d_mono  (i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%mass = Tr3d_mass  (i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%vmin = Tr3d_vmin  (i)
               Tr_3CWP(Tr3d_ntrTRICUB_WP)%vmax = Tr3d_vmax  (i)
            endif
         endif
         if (Tr3d_intp(i) == 'BICUBH_QV') then
            if ( (Tr3d_mass(i) < 1) .and. (Tr3d_mono(i) < 1)) then
               Tr3d_ntrBICHQV_NT= Tr3d_ntrBICHQV_NT+1
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%name = Tr3d_name_S(i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%intp = Tr3d_intp  (i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%wload= Tr3d_wload (i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%hzd  = Tr3d_hzd   (i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%mono = Tr3d_mono  (i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%mass = Tr3d_mass  (i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%vmin = Tr3d_vmin  (i)
               Tr_BQNT(Tr3d_ntrBICHQV_NT)%vmax = Tr3d_vmax  (i)
            else
               Tr3d_ntrBICHQV_WP= Tr3d_ntrBICHQV_WP+1
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%name = Tr3d_name_S(i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%intp = Tr3d_intp  (i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%wload= Tr3d_wload (i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%hzd  = Tr3d_hzd   (i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%mono = Tr3d_mono  (i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%mass = Tr3d_mass  (i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%vmin = Tr3d_vmin  (i)
               Tr_BQWP(Tr3d_ntrBICHQV_WP)%vmax = Tr3d_vmax  (i)
            endif
         endif
      end do

      n=0
      Tr3d_debTRICUB_NT=n+1
      do i=1, Tr3d_ntrTRICUB_NT
         n=n+1
         Tr3d_name_S(n) = Tr_3CNT(i)%name
         Tr3d_intp  (n) = Tr_3CNT(i)%intp
         Tr3d_wload (n) = Tr_3CNT(i)%wload
         Tr3d_hzd   (n) = Tr_3CNT(i)%hzd
         Tr3d_mono  (n) = Tr_3CNT(i)%mono
         Tr3d_mass  (n) = Tr_3CNT(i)%mass
         Tr3d_vmin  (n) = Tr_3CNT(i)%vmin
         Tr3d_vmax  (n) = Tr_3CNT(i)%vmax
      end do
      Tr3d_debBICHQV_NT=n+1
      do i=1, Tr3d_ntrBICHQV_NT
         n=n+1
         Tr3d_name_S(n) = Tr_BQNT(i)%name
         Tr3d_intp  (n) = Tr_BQNT(i)%intp
         Tr3d_wload (n) = Tr_BQNT(i)%wload
         Tr3d_hzd   (n) = Tr_BQNT(i)%hzd
         Tr3d_mono  (n) = Tr_BQNT(i)%mono
         Tr3d_mass  (n) = Tr_BQNT(i)%mass
         Tr3d_vmin  (n) = Tr_BQNT(i)%vmin
         Tr3d_vmax  (n) = Tr_BQNT(i)%vmax
      end do
      Tr3d_debTRICUB_WP=n+1
      do i=1, Tr3d_ntrTRICUB_WP
         n=n+1
         Tr3d_name_S(n) = Tr_3CWP(i)%name
         Tr3d_intp  (n) = Tr_3CWP(i)%intp
         Tr3d_wload (n) = Tr_3CWP(i)%wload
         Tr3d_hzd   (n) = Tr_3CWP(i)%hzd
         Tr3d_mono  (n) = Tr_3CWP(i)%mono
         Tr3d_mass  (n) = Tr_3CWP(i)%mass
         Tr3d_vmin  (n) = Tr_3CWP(i)%vmin
         Tr3d_vmax  (n) = Tr_3CWP(i)%vmax
      end do
      Tr3d_debBICHQV_WP=n+1
      do i=1, Tr3d_ntrBICHQV_WP
         n=n+1
         Tr3d_name_S(n) = Tr_BQWP(i)%name
         Tr3d_intp  (n) = Tr_BQWP(i)%intp
         Tr3d_wload (n) = Tr_BQWP(i)%wload
         Tr3d_hzd   (n) = Tr_BQWP(i)%hzd
         Tr3d_mono  (n) = Tr_BQWP(i)%mono
         Tr3d_mass  (n) = Tr_BQWP(i)%mass
         Tr3d_vmin  (n) = Tr_BQWP(i)%vmin
         Tr3d_vmax  (n) = Tr_BQWP(i)%vmax
      end do
      do i=1,Tr3d_ntr
         if (trim(Tr3d_name_S(i))=='HU') Tr3d_hu=i
      end do

      Tr_3CWP(:)%BC_mass_deficit= 0.d0
      Tr_BQWP(:)%BC_mass_deficit= 0.d0
      if (Rstri_rstn_L) then
         Tr_3CWP(1:Tr3d_ntrTRICUB_WP)%BC_mass_deficit = &
         BCMD_3C(1:Tr3d_ntrTRICUB_WP)
         Tr_BQWP(1:Tr3d_ntrBICHQV_WP)%BC_mass_deficit = &
         BCMD_BQ(1:Tr3d_ntrBICHQV_WP)
      endif

      if (Lun_out > 0) then
         write (Lun_out,1001)
         do i=1,Tr3d_ntr
            write(Lun_out,1002) Tr3d_name_S(i),Tr3d_wload(i),&
                       Tr3d_hzd(i),Tr3d_mono(i),Tr3d_mass(i),&
                       Tr3d_vmin(i),Tr3d_vmax(i),Tr3d_intp(i)
         end do
      end if

      call ac_posi (geomh_longs(1),geomh_latgs(1),G_ni,G_nj,Lun_out > 0)

 1001 format (/' Final list of tracers:'/3x,' Name   Wload  Hzd   Mono  Mass    Min        Max     Intp')
 1002 format (4x,a4,2l6,2i6,3x,2(e10.3,1x),a12)
!
!     __________________________________________________________________
!
      return
      end
