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

module yyg_param

   use glb_ld
   use glb_pil
   use gem_options
   use geomh
   use tdpack
   use ptopo
   use, intrinsic :: iso_fortran_env
   implicit none
   public
   save

   type :: YYG_comm_param
      sequence
      integer :: send_all, recv_all, sendmaxproc, recvmaxproc, &
                 maxsend, maxrecv
      integer, dimension (:), pointer, contiguous :: &
            sendproc, recvproc, recv_len, send_len,  &
            recv_adr, send_adr, recv_i, recv_j    ,  &
            send_ixu, send_iyu, send_ixv, send_iyv
      real(kind=REAL64),  dimension (:  ), allocatable :: &
            send_xu, send_yu, send_xv, send_yv, send_s, &
            send_xxr, send_yyr
   end type YYG_comm_param

   type(YYG_comm_param) :: YYG_PILT_q2q, YYG_PILT_uv2u, YYG_PILT_uv2v
   type(YYG_comm_param) :: YYG_BLEN_q2q, YYG_BLEN_uv2u, YYG_BLEN_uv2v
   type(YYG_comm_param) :: YYG_HALO_q2q, YYG_NEAR_q2q

   integer :: YYG_low_CoL= 1, YYG_high_CoL= 3
   integer :: YYG_i0, YYG_j0, YYG_lni, YYG_lnj
   integer, dimension (:,:), allocatable :: YYG_pe_indx
   integer :: YYG_PILT_recv_uv, YYG_BLEN_recv_uv
   real  ,  dimension(:,:), allocatable :: YYG_usenduv,YYG_vsenduv,&
                                           YYG_urecvuv,YYG_vrecvuv
   real  ,  dimension(:,:), allocatable :: YYG_uvsend, YYG_uvrecv

   real(kind=REAL64),  dimension(:  ), allocatable :: YYG_xg_8 , YYG_yg_8, &
                                           YYG_xgu_8, YYG_ygv_8

contains

   subroutine yyg_indices ( F_x1,F_y1,F_x2,F_y2, F_xp,F_yp, F_s    ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            F_i0,F_j0, NIU, NJU, NIV, NJV )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_i0,F_j0,NIU, NJU, NIV, NJV
      integer, intent(out) :: F_x1,F_y1,F_x2,F_y2
      real(kind=REAL64) , intent(in) :: F_posx(1-G_ni:2*G_ni), F_posy(1-G_nj:2*G_nj) ! target positions
      real(kind=REAL64) , intent(in) :: F_pux (1-G_ni:2*G_ni), F_puy (1-G_nj:2*G_nj) ! u source positions
      real(kind=REAL64) , intent(in) :: F_pvx (1-G_ni:2*G_ni), F_pvy (1-G_nj:2*G_nj) ! v source positions
      real(kind=REAL64) , intent(out) :: F_xp,F_yp, F_s(4)
!
!----------------------------------------------------------------------
!
      ! find the global location (F_xp,F_yp) in this panel of the
      ! target interpolation point (F_posx,F_posy) of the other panel
      call smat ( F_s, F_xp,F_yp, F_posx(F_i0)-pi_8,F_posy(F_j0) )
      F_xp= F_xp+pi_8
      ! find the global lower indexes (imx1,imy1) for cubic
      ! interpolation of u-component to target location (x_a,y_a)
      call localise2 ( F_x1,F_y1, F_xp,F_yp,&
                      F_pux,F_puy, geomh_hx_8,geomh_hy_8,1,1,G_ni,G_nj )
      ! find the global lower indexes (imx2,imy2) for cubic
      ! interpolation of v-component to target location (x_a,y_a)
      call localise2 ( F_x2,F_y2, F_xp,F_yp,&
                      F_pvx,F_pvy, geomh_hx_8,geomh_hy_8,1,1,G_ni,G_nj )
      ! impose global limits to these indexes
      F_x1 = min(max(F_x1-YYG_low_CoL,glb_pil_w+1),&
                     NIU -glb_pil_e-YYG_high_CoL)
      F_y1 = min(max(F_y1-YYG_low_CoL,glb_pil_s+1),&
                     NJU -glb_pil_n-YYG_high_CoL)
      F_x2 = min(max(F_x2-YYG_low_CoL,glb_pil_w+1),&
                     NIV -glb_pil_e-YYG_high_CoL)
      F_y2 = min(max(F_y2-YYG_low_CoL,glb_pil_s+1),&
                     NJV -glb_pil_n-YYG_high_CoL)
!
!----------------------------------------------------------------------
!
   return
   end subroutine yyg_indices

   subroutine yyg_count ( F_recv,F_send, F_i0,F_j0, F_x1,F_y1,F_x2,F_y2 )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: F_i0,F_j0,F_x1,F_y1,F_x2,F_y2
      integer, intent(inout) :: F_recv(*), F_send(*)

      integer kk
!
!----------------------------------------------------------------------
!
! check to collect from who
      if (F_i0 >= YYG_i0 .and. F_i0 <= l_i0+YYG_lni-1 .and. &
          F_j0 >= YYG_j0 .and. F_j0 <= l_j0+YYG_lnj-1 ) then
         do kk= 1, Ptopo_numproc
            if (max(F_x1,F_x2) >= Ptopo_gindx(1,kk).and. &
                max(F_x1,F_x2) <= Ptopo_gindx(2,kk).and. &
                max(F_y1,F_y2) >= Ptopo_gindx(3,kk).and. &
                max(F_y1,F_y2) <= Ptopo_gindx(4,kk) ) then
               F_recv(kk)= F_recv(kk) + 1
            end if
         end do
      end if

! check to send to who
      if (max(F_x1,F_x2) >= l_i0.and.         &
          max(F_x1,F_x2) <= l_i0+l_ni-1 .and. &
          max(F_y1,F_y2) >= l_j0.and.         &
          max(F_y1,F_y2) <= l_j0+l_nj-1) then
         do kk=1,Ptopo_numproc
            if (F_i0 >= YYG_pe_indx(1,kk) .and. &
                F_i0 <= YYG_pe_indx(2,kk) .and. &
                F_j0 >= YYG_pe_indx(3,kk) .and. &
                F_j0 <= YYG_pe_indx(4,kk) )then
               F_send(kk)= F_send(kk) + 1
            end if
         end do
      end if
!
!----------------------------------------------------------------------
!
   return
   end subroutine yyg_count

   subroutine yyg_allocate ( F_comm, recv_len,send_len )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: recv_len(*), send_len(*)
      type(YYG_comm_param), intent(inout) :: F_comm

      integer kk,ksend,krecv
!
!----------------------------------------------------------------------
!
! Obtain sum of elements to send and receive for each processor
! and the total memory needed to store and receive for each processor

      F_comm%send_all=0
      F_comm%recv_all=0
      F_comm%sendmaxproc=0
      F_comm%recvmaxproc=0
      F_comm%maxsend=0
      F_comm%maxrecv=0

      do kk=1,Ptopo_numproc
         F_comm%send_all=send_len(kk)+F_comm%send_all
         F_comm%recv_all=recv_len(kk)+F_comm%recv_all

         if (send_len(kk) > 0) F_comm%sendmaxproc=F_comm%sendmaxproc+1
         if (recv_len(kk) > 0) F_comm%recvmaxproc=F_comm%recvmaxproc+1
         F_comm%maxsend= max(F_comm%maxsend,send_len(kk))
         F_comm%maxrecv= max(F_comm%maxrecv,recv_len(kk))
      end do

      allocate (F_comm%recvproc(F_comm%recvmaxproc))
      allocate (F_comm%recv_len(F_comm%recvmaxproc))
      allocate (F_comm%recv_adr(F_comm%recvmaxproc))

      allocate (F_comm%sendproc(F_comm%sendmaxproc))
      allocate (F_comm%send_len(F_comm%sendmaxproc))
      allocate (F_comm%send_adr(F_comm%sendmaxproc))
      F_comm%recv_len(:) = 0
      F_comm%send_len(:) = 0
      F_comm%recv_adr(:) = 0
      F_comm%send_adr(:) = 0

      ksend=0
      krecv=0
      F_comm%send_all=0
      F_comm%recv_all=0

! Fill the lengths and addresses for selected processors to communicate

      do kk=1,Ptopo_numproc
         if (send_len(kk) > 0) then
            ksend=ksend+1
            F_comm%sendproc(ksend)=kk
            F_comm%send_len(ksend)=send_len(kk)

            F_comm%send_adr(ksend)= F_comm%send_all
            F_comm%send_all= F_comm%send_all + F_comm%send_len(ksend)
         end if
         if (recv_len(kk) > 0) then
            krecv=krecv+1
            F_comm%recvproc(krecv)=kk
            F_comm%recv_len(krecv)=recv_len(kk)

            F_comm%recv_adr(krecv)= F_comm%recv_all
            F_comm%recv_all= F_comm%recv_all + F_comm%recv_len(krecv)
         end if

      end do

! Allocate the vectors needed for sending and receiving each processor

      if (F_comm%recv_all > 0) then
         allocate (F_comm%recv_i(F_comm%recv_all))
         allocate (F_comm%recv_j(F_comm%recv_all))
      end if

      if (F_comm%send_all > 0) then
         allocate (F_comm%send_ixu(F_comm%send_all))
         allocate (F_comm%send_iyu(F_comm%send_all))
         allocate (F_comm%send_ixv(F_comm%send_all))
         allocate (F_comm%send_iyv(F_comm%send_all))
         allocate (F_comm%send_xxr(F_comm%send_all))
         allocate (F_comm%send_yyr(F_comm%send_all))
         allocate (F_comm%send_xu(F_comm%send_all))
         allocate (F_comm%send_yu(F_comm%send_all))
         allocate (F_comm%send_xv(F_comm%send_all))
         allocate (F_comm%send_yv(F_comm%send_all))
         allocate (F_comm%send_s (F_comm%send_all*4))
      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine yyg_allocate

   subroutine yyg_setcomm ( F_comm, F_i0,F_j0, F_recv,F_send,&
                            F_x1,F_y1,F_x2,F_y2, F_xp,F_yp  ,&
                            F_s, F_pux,F_puy, F_pvx,F_pvy )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in)    :: F_i0,F_j0,F_x1,F_y1,F_x2,F_y2
      integer, intent(inout) :: F_recv(*), F_send(*)
      real(kind=REAL64) , intent(in) :: F_pux (1-G_ni:2*G_ni), F_puy (1-G_nj:2*G_nj) ! u source positions
      real(kind=REAL64) , intent(in) :: F_pvx (1-G_ni:2*G_ni), F_pvy (1-G_nj:2*G_nj) ! v source positions
      real(kind=REAL64) , intent(in) :: F_xp,F_yp, F_s(2,2)
      type(YYG_comm_param), intent(inout) :: F_comm

      integer kk,ki,iu,ju,iv,jv,next,nexts
!
!----------------------------------------------------------------------
!
! check to collect from who
      if (F_i0 >= YYG_i0 .and. F_i0 <= l_i0+YYG_lni-1 .and. &
          F_j0 >= YYG_j0 .and. F_j0 <= l_j0+YYG_lnj-1 ) then
         do kk= 1, F_comm%recvmaxproc
            ki=F_comm%recvproc(kk)
            if (max(F_x1,F_x2) >= Ptopo_gindx(1,ki).and. &
                max(F_x1,F_x2) <= Ptopo_gindx(2,ki).and. &
                max(F_y1,F_y2) >= Ptopo_gindx(3,ki).and. &
                max(F_y1,F_y2) <= Ptopo_gindx(4,ki) ) then
               F_recv(kk)= F_recv(kk) + 1
               F_comm%recv_i(F_comm%recv_adr(kk)+F_recv(kk))=F_i0-l_i0+1
               F_comm%recv_j(F_comm%recv_adr(kk)+F_recv(kk))=F_j0-l_j0+1
            end if
         end do
      end if

! check to send to who
      if (max(F_x1,F_x2) >= l_i0.and.         &
          max(F_x1,F_x2) <= l_i0+l_ni-1 .and. &
          max(F_y1,F_y2) >= l_j0.and.         &
          max(F_y1,F_y2) <= l_j0+l_nj-1) then
         do kk=1,F_comm%sendmaxproc
            ki=F_comm%sendproc(kk)
            if (F_i0 >= YYG_pe_indx(1,ki) .and. &
                F_i0 <= YYG_pe_indx(2,ki) .and. &
                F_j0 >= YYG_pe_indx(3,ki) .and. &
                F_j0 <= YYG_pe_indx(4,ki) )then
               F_send(kk)= F_send(kk) + 1
               iu= F_x1-l_i0+1 ; iv= F_x2-l_i0+1
               ju= F_y1-l_j0+1 ; jv= F_y2-l_j0+1
               next=F_comm%send_adr(kk)+F_send(kk)
               nexts=(next-1)*4+1
               F_comm%send_ixu(next)=iu
               F_comm%send_iyu(next)=ju
               F_comm%send_ixv(next)=iv
               F_comm%send_iyv(next)=jv
               F_comm%send_xxr(next)=F_xp
               F_comm%send_yyr(next)=F_yp
               F_comm%send_xu(next)= (F_xp-F_pux(F_x1))*geomh_inv_hx_8+F_x1
               F_comm%send_yu(next)= (F_yp-F_puy(F_y1))*geomh_inv_hy_8+F_y1
               F_comm%send_xv(next)= (F_xp-F_pvx(F_x2))*geomh_inv_hx_8+F_x2
               F_comm%send_yv(next)= (F_yp-F_pvy(F_y2))*geomh_inv_hy_8+F_y2
               F_comm%send_s (nexts  )=F_s(1,1)
               F_comm%send_s (nexts+1)=F_s(1,2)
               F_comm%send_s (nexts+2)=F_s(2,1)
               F_comm%send_s (nexts+3)=F_s(2,2)
            end if
         end do
      end if
!
!----------------------------------------------------------------------
!
   return
   end subroutine yyg_setcomm

      subroutine yyg_initcomm ( F_comm, F_posx, F_posy    , &
                                F_pux, F_puy, F_pvx, F_pvy, &
                                NI, NJ, F_halox, F_haloy  , &
                                NIU, NJU, NIV, NJV, F_inttype_S )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*), optional, intent(in) :: F_inttype_S
      integer, intent(in) :: NI, NJ, NIU, NJU, NIV, NJV, &
                             F_halox, F_haloy
      real(kind=REAL64) , intent(in) :: F_posx(1-G_ni:2*G_ni), F_posy(1-G_nj:2*G_nj) ! target positions
      real(kind=REAL64) , intent(in) :: F_pux (1-G_ni:2*G_ni), F_puy (1-G_nj:2*G_nj) ! u source positions
      real(kind=REAL64) , intent(in) :: F_pvx (1-G_ni:2*G_ni), F_pvy (1-G_nj:2*G_nj) ! v source positions
      type(YYG_comm_param) , intent(out) :: F_comm

      integer i,j,imx1,imx2,imy1,imy2
      integer, dimension (:), allocatable :: recv_len,send_len
      real(kind=REAL64)  s(2,2),x_a,y_a
!
!     ---------------------------------------------------------------
!
      if (present(F_inttype_S)) then
         if (trim(F_inttype_S) /= 'CUBIC') then
            YYG_low_CoL= 0 ; YYG_high_CoL= 1
         end if
      end if

      allocate (recv_len (Ptopo_numproc)) ; recv_len (:)=0
      allocate (send_len (Ptopo_numproc)) ; send_len (:)=0

! FIRST PASS is to find the number of processor to tag for
! communication and the number of items to send and receive for each
! processor before allocating the vectors

! WEST section
      do j=1-F_haloy, NJ+F_haloy
      do i=1-F_halox, glb_pil_w
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_count (recv_len,send_len,i,j,imx1,imy1,imx2,imy2)
      end do
      end do

! East section
      do j=1-F_haloy, NJ+F_haloy
      do i=NI-glb_pil_e+1,NI+F_halox
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_count (recv_len,send_len,i,j,imx1,imy1,imx2,imy2)
      end do
      end do

! South section
      do j=1-F_haloy,glb_pil_s
      do i=1+glb_pil_w,NI-glb_pil_e
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_count (recv_len,send_len,i,j,imx1,imy1,imx2,imy2)
      end do
      end do

! North section
      do j=NJ-glb_pil_n+1,NJ+F_haloy
      do i=1+glb_pil_w,NI-glb_pil_e
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_count (recv_len,send_len,i,j,imx1,imy1,imx2,imy2)
      end do
      end do

      call yyg_allocate ( F_comm, recv_len, send_len )

      recv_len(:)=0 ; send_len(:)=0

! SECOND PASS is to initialize the vectors with
! information for communication

! WEST section
      do j=1-F_haloy, NJ+F_haloy
      do i=1-F_halox, glb_pil_w
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_setcomm ( F_comm, i,j, recv_len,send_len,&
                            imx1,imy1,imx2,imy2, x_a,y_a,s,&
                            F_pux,F_puy, F_pvx,F_pvy )
      end do
      end do

! East section
      do j=1-F_haloy, NJ+F_haloy
      do i=NI-glb_pil_e+1,NI+F_halox
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_setcomm ( F_comm, i,j, recv_len,send_len,&
                            imx1,imy1,imx2,imy2, x_a,y_a,s,&
                            F_pux,F_puy, F_pvx,F_pvy )
      end do
      end do

! South section
      do j=1-F_haloy,glb_pil_s
      do i=1+glb_pil_w,NI-glb_pil_e
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_setcomm ( F_comm, i,j, recv_len,send_len,&
                            imx1,imy1,imx2,imy2, x_a,y_a,s,&
                            F_pux,F_puy, F_pvx,F_pvy )
      end do
      end do

! North section
      do j=NJ-glb_pil_n+1,NJ+F_haloy
      do i=1+glb_pil_w,NI-glb_pil_e
         call yyg_indices ( imx1,imy1,imx2,imy2, x_a,y_a,s         ,&
                            F_posx,F_posy, F_pux,F_puy, F_pvx,F_pvy,&
                            i,j, NIU, NJU, NIV, NJV )
         call yyg_setcomm ( F_comm, i,j, recv_len,send_len,&
                            imx1,imy1,imx2,imy2, x_a,y_a,s,&
                            F_pux,F_puy, F_pvx,F_pvy )
      end do
      end do

      deallocate (recv_len,send_len)
      YYG_low_CoL= 1 ; YYG_high_CoL= 3
!
!-------------------------------------------------------------------
!
      return
      end subroutine yyg_initcomm

      subroutine yyg_initblen ( F_comm, F_posx, F_posy    , &
                                F_pux, F_puy, F_pvx, F_pvy, &
                                NI, NJ, NIU, NJU, NIV, NJV )
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: NI, NJ, NIU, NJU, NIV, NJV
      real(kind=REAL64) , intent(in) :: F_posx(1-G_ni:2*G_ni), F_posy(1-G_nj:2*G_nj) ! target positions
      real(kind=REAL64) , intent(in) :: F_pux (1-G_ni:2*G_ni), F_puy (1-G_nj:2*G_nj) ! u source positions
      real(kind=REAL64) , intent(in) :: F_pvx (1-G_ni:2*G_ni), F_pvy (1-G_nj:2*G_nj) ! v source positions
      type(YYG_comm_param) , intent(out) :: F_comm

!author
!     Abdessamad Qaddouri/V.Lee - October 2009

      integer i,j,imx1,imx2,imy1,imy2
      integer, dimension (:), allocatable :: recv_len,send_len
      real(kind=REAL64)  s(2,2),x_a,y_a
!
!     ---------------------------------------------------------------
!
      allocate (recv_len (Ptopo_numproc)) ; recv_len (:)=0
      allocate (send_len (Ptopo_numproc)) ; send_len (:)=0

! FIRST PASS is to find the number of processor to tag for
! communication and the number of items to send and receive for each
! processor before allocating the vectors

      do j=1+glb_pil_s, NJ-glb_pil_n
      do i=1+glb_pil_w, NI-glb_pil_e

         call smat ( s, x_a,y_a, F_posx(i)-pi_8,F_posy(j) )
         x_a= x_a+pi_8
         call localise_blend(imx1,imy1,x_a,y_a,F_pux,F_puy,&
                             1-G_ni,2*NIU,1-G_nj,2*G_nj,&
                             geomh_hx_8,geomh_hy_8)
         call localise_blend(imx2,imy2,x_a,y_a,F_pvx,F_pvy,&
                             1-G_ni,2*G_ni,1-G_nj,2*NJV,&
                             geomh_hx_8,geomh_hy_8)
         if (imx1 > 1+glb_pil_w .and. imx1 < NIU-glb_pil_e .and. &
             imy1 > 1+glb_pil_s .and. imy1 < NJU-glb_pil_n .and. &
             imx2 > 1+glb_pil_w .and. imx2 < NIV-glb_pil_e .and. &
             imy2 > 1+glb_pil_s .and. imy2 < NJV-glb_pil_n) then

            imx1 = min(max(imx1-1,glb_pil_w+1),NIU-glb_pil_e-3)
            imy1 = min(max(imy1-1,glb_pil_s+1),NJU-glb_pil_n-3)
            imx2 = min(max(imx2-1,glb_pil_w+1),NIV-glb_pil_e-3)
            imy2 = min(max(imy2-1,glb_pil_s+1),NJV-glb_pil_n-3)

            call yyg_count (recv_len,send_len,i,j,imx1,imy1,imx2,imy2)
         end if

      end do
      end do

      call yyg_allocate ( F_comm, recv_len, send_len )

      recv_len(:)=0 ; send_len(:)=0

! SECOND PASS is to initialize the vectors with
! information for communication

      do j=1+glb_pil_s, NJ-glb_pil_n
      do i=1+glb_pil_w, NI-glb_pil_e

         call smat ( s, x_a,y_a, F_posx(i)-pi_8,F_posy(j) )
         x_a= x_a+pi_8
         call localise_blend(imx1,imy1,x_a,y_a,F_pux,F_puy,&
                             1-G_ni,2*NIU,1-G_nj,2*G_nj,&
                             geomh_hx_8,geomh_hy_8)
         call localise_blend(imx2,imy2,x_a,y_a,F_pvx,F_pvy,&
                             1-G_ni,2*G_ni,1-G_nj,2*NJV,&
                             geomh_hx_8,geomh_hy_8)
         if (imx1 > 1+glb_pil_w .and. imx1 < NIU-glb_pil_e .and. &
             imy1 > 1+glb_pil_s .and. imy1 < NJU-glb_pil_n .and. &
             imx2 > 1+glb_pil_w .and. imx2 < NIV-glb_pil_e .and. &
             imy2 > 1+glb_pil_s .and. imy2 < NJV-glb_pil_n) then

            imx1 = min(max(imx1-1,glb_pil_w+1),NIU-glb_pil_e-3)
            imy1 = min(max(imy1-1,glb_pil_s+1),NJU-glb_pil_n-3)
            imx2 = min(max(imx2-1,glb_pil_w+1),NIV-glb_pil_e-3)
            imy2 = min(max(imy2-1,glb_pil_s+1),NJV-glb_pil_n-3)

            call yyg_setcomm ( F_comm, i,j, recv_len,send_len,&
                               imx1,imy1,imx2,imy2, x_a,y_a,s,&
                               F_pux, F_puy, F_pvx, F_pvy )
         end if

      end do
      end do

      deallocate (recv_len,send_len)
!
!-------------------------------------------------------------------
!
      return
      end subroutine yyg_initblen

      subroutine yyg_SendRecv_s ( F_src, F_comm, Minx,Maxx,Miny,Maxy, &
                                  NK, F_interpo_S, mono_L)
      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in)  :: F_interpo_S
      logical, intent(in) :: mono_L
      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real   , intent(in) :: F_src (Minx:Maxx,Miny:Maxy,Nk)
      type(YYG_comm_param) , intent(inout) :: F_comm

      include "intrp_bicub_yx.inc"

      integer tag2, ireq, kk, m, mm, kk_proc, adr, ind1, ind2, ierr
      integer request(Ptopo_numproc*3)
!
!     ---------------------------------------------------------------
!
      call rpn_comm_xch_halo (F_src, Minx,Maxx,Miny,Maxy,l_ni,l_nj,Nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      tag2=14 ; ireq=0 ; request = 0

! Posting the receive requests
      do kk= 1, F_comm%recvmaxproc

         kk_proc= yyg_proc(F_comm%recvproc(kk))
         m= F_comm%recv_len(kk)
         if (m > 0) then
            ireq = ireq+1
            call RPN_COMM_IRecv (YYG_uvrecv(1,KK), m*NK          ,&
                                 'MPI_REAL', kk_proc,tag2+kk_proc,&
                                 'MULTIGRID',request(ireq),ierr)
         end if

      end do

! Interpolation and shipping
      do kk= 1, F_comm%sendmaxproc

         kk_proc= yyg_proc(F_comm%sendproc(kk))
         mm= F_comm%send_len(kk)
         if (mm > 0) then

             adr=F_comm%send_adr(kk)

             if (trim(F_interpo_S) == 'CUBIC') then
                if (mono_l) then
                   do m=1,mm
                      ind1= (m-1)*nk+1
                      ind2= adr+m
                      call intrp_bicub_yx_s_mono ( &
                                   F_src(1,1,1), YYG_uvsend(ind1,kk),&
                          (Maxx-Minx+1) ,(Maxx-Minx+1)*(Maxy-Miny+1),&
                          nk, F_comm%send_xu(ind2),F_comm%send_yu(ind2) )
                   end do
                else
                   do m=1,mm
                      ind1= (m-1)*nk+1
                      ind2= adr+m
                      call intrp_bicub_yx_s ( &
                                   F_src(1,1,1), YYG_uvsend(ind1,kk),&
                          (Maxx-Minx+1) ,(Maxx-Minx+1)*(Maxy-Miny+1),&
                          nk, F_comm%send_xu(ind2),F_comm%send_yu(ind2) )
                   end do
                end if
!!$             call int_cub_lag ( YYG_uvsend(1,kk), F_src,&
!!$                                F_comm%send_ixu(adr+1) ,&
!!$                                F_comm%send_iyu(adr+1) ,&
!!$                                geomh_x_8,geomh_y_8    ,&
!!$                                Minx,Maxx,Miny,Maxy,Nk ,&
!!$                                F_comm%send_xxr(adr+1) ,&
!!$                                F_comm%send_yyr(adr+1) ,&
!!$                                mm,mono_l )
             elseif (trim(F_interpo_S) == 'LINEAR') then
                call yyg_int_lin  ( YYG_uvsend(1,KK), F_src,&
                                    F_comm%send_ixu(adr+1:),&
                                    F_comm%send_iyu(adr+1:),&
                                    geomh_x_8,geomh_y_8    ,&
                                    Minx,Maxx,Miny,Maxy,Nk ,&
                                    F_comm%send_xxr(adr+1) ,&
                                    F_comm%send_yyr(adr+1), mm )
             elseif (trim(F_interpo_S) == 'NEAREST') then
                call yyg_int_near ( YYG_uvsend(1,KK), F_src,&
                                    F_comm%send_ixu(adr+1:),&
                                    F_comm%send_iyu(adr+1:),&
                                    geomh_x_8,geomh_y_8    ,&
                                    Minx,Maxx,Miny,Maxy,Nk ,&
                                    F_comm%send_xxr(adr+1) ,&
                                    F_comm%send_yyr(adr+1), mm )
             end if

             ireq = ireq+1
             call RPN_COMM_ISend ( YYG_uvsend(1,KK),mm*NK, 'MPI_REAL',&
                                   kk_proc, tag2+Ptopo_world_myproc  ,&
                                   'MULTIGRID', request(ireq), ierr )
         end if

      end do

      call RPN_COMM_waitall_nostat(ireq,request,ierr)
!
!-------------------------------------------------------------------
!
      return
      end subroutine yyg_SendRecv_s

      subroutine yyg_SendRecv_v ( F_u, F_v, F_commU, F_commV, &
                                  Minx,Maxx,Miny,Maxy,NK )
      use ISO_C_BINDING
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: Minx,Maxx,Miny,Maxy, Nk
      real, dimension(Minx:Maxx,Miny:Maxy,Nk), intent(inout) :: F_u, F_v
      type(YYG_comm_param) , intent(inout) :: F_commU, F_commV

      include "intrp_bicub_yx.inc"

      integer tag1, tag2, ireq, kk, m, mm, kk_proc, adr, ind1, ind2, ierr
      integer request(Ptopo_numproc*4)
      real dummy(NK)
!
!     ---------------------------------------------------------------
!
      call rpn_comm_xch_halo (F_u, Minx,Maxx,Miny,Maxy,l_ni,l_nj,Nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)
      call rpn_comm_xch_halo (F_v, Minx,Maxx,Miny,Maxy,l_ni,l_nj,Nk, &
                              G_halox,G_haloy,G_periodx,G_periody,l_ni,0)

      tag1= 12 ; tag2= 14 ; ireq= 0

! Posting the receive requests
      do kk= 1, F_commU%recvmaxproc
         kk_proc= yyg_proc(F_commU%recvproc(kk))
         m= F_commU%recv_len(kk)
         if (m > 0) then
            ireq = ireq+1
            call RPN_COMM_IRecv (YYG_urecvuv(1,KK), m*NK         ,&
                                 'MPI_REAL', kk_proc,tag1+kk_proc,&
                                 'MULTIGRID',request(ireq),ierr)
         end if
      end do

      do kk= 1, F_commV%recvmaxproc
         kk_proc= yyg_proc(F_commV%recvproc(kk))
         m= F_commV%recv_len(kk)
         if (m > 0) then
            ireq = ireq+1
            call RPN_COMM_IRecv (YYG_vrecvuv(1,KK), m*NK         ,&
                                 'MPI_REAL', kk_proc,tag2+kk_proc,&
                                 'MULTIGRID',request(ireq),ierr)
         end if
      end do

! Interpolation and shipping
      do kk= 1, F_commU%sendmaxproc

         kk_proc= yyg_proc(F_commU%sendproc(kk))
         mm= F_commU%send_len(kk)
         if (mm > 0) then

             adr= F_commU%send_adr(kk)
             do m=1, mm
                ind1= (m-1)*nk+1
                ind2= adr+m
                call intrp_bicub_yx_vec ( F_u(1,1,1), F_v(1,1,1) ,&
                                      YYG_usenduv(ind1,kk), dummy,&
                   (Maxx-Minx+1), (Maxx-Minx+1)*(Maxy-Miny+1), nk,&
                                      F_commU%send_xu(ind2),&
                                      F_commU%send_yu(ind2),&
                                      F_commU%send_xv(ind2),&
                                      F_commU%send_yv(ind2),&
                               F_commU%send_s((ind2-1)*4+1) )
             end do

             ireq = ireq+1
             call RPN_COMM_ISend ( YYG_usenduv(1,KK), mm*NK, 'MPI_REAL',&
                                   kk_proc, tag1+Ptopo_world_myproc,&
                                   'MULTIGRID', request(ireq), ierr )
         end if

      end do

      do kk= 1, F_commV%sendmaxproc

         kk_proc= yyg_proc(F_commV%sendproc(kk))
         mm= F_commV%send_len(kk)
         if (mm > 0) then

             adr=F_commV%send_adr(kk)
             do m=1, mm
                ind1= (m-1)*nk+1
                ind2= adr+m
                call intrp_bicub_yx_vec ( F_u(1,1,1), F_v(1,1,1) ,&
                                      dummy, YYG_vsenduv(ind1,kk),&
                   (Maxx-Minx+1), (Maxx-Minx+1)*(Maxy-Miny+1), nk,&
                                      F_commV%send_xu(ind2),&
                                      F_commV%send_yu(ind2),&
                                      F_commV%send_xv(ind2),&
                                      F_commV%send_yv(ind2),&
                               F_commV%send_s((ind2-1)*4+1) )
             end do

             ireq = ireq+1
             call RPN_COMM_ISend ( YYG_vsenduv(1,KK), mm*NK, 'MPI_REAL',&
                                   kk_proc, tag2+Ptopo_world_myproc,&
                                   'MULTIGRID', request(ireq), ierr )
         end if

      end do

      call RPN_COMM_waitall_nostat(ireq,request,ierr)
!
!-------------------------------------------------------------------
!
      return
      end subroutine yyg_SendRecv_v

      subroutine yyg_extend_grid ()
      use gem_options
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>

      integer i,j
      real(kind=REAL64), parameter :: TWO_8 = 2.0d0
!
!-------------------------------------------------------------------
!
      allocate (YYG_xg_8 (1-G_ni:2*G_ni),YYG_yg_8 (1-G_nj:2*G_nj),&
                YYG_xgu_8(1-G_ni:2*G_ni),YYG_ygv_8(1-G_nj:2*G_nj) )

      do i=1,G_ni
         YYG_xg_8(i) = G_xg_8(i)
      end do
      do j=1,G_nj
         YYG_yg_8(j) = G_yg_8(j)
      end do

      do i=-G_ni+1,0
         YYG_xg_8(i) = YYG_xg_8(1) + (i-1)*geomh_hx_8
      end do
      do i=G_ni+1,2*G_ni
         YYG_xg_8(i) = YYG_xg_8(G_ni) + (i-G_ni)*geomh_hx_8
      end do

      do j=-G_nj+1,0
         YYG_yg_8(j) = YYG_yg_8(1) + (j-1)*geomh_hy_8
      end do
      do j=G_nj+1,2*G_nj
         YYG_yg_8(j) = YYG_yg_8(G_nj) + (j-G_nj)*geomh_hy_8
      end do

      do i=1-G_ni,2*G_ni-1
         YYG_xgu_8(i)=0.5D0 *(YYG_xg_8(i+1)+YYG_xg_8(i))
      end do
      do j=1-G_nj,2*G_nj-1
         YYG_ygv_8(j)= 0.5D0*(YYG_yg_8(j+1)+YYG_yg_8(j))
      end do
!
!-------------------------------------------------------------------
!
      return
      end subroutine yyg_extend_grid

   integer function yyg_proc (F_proc)

      use, intrinsic :: iso_fortran_env
      implicit none
      integer, intent(in) :: F_proc
!
!----------------------------------------------------------------------
!
! MID: This logic will probably not work when we change
!      the logical processor layout
      yyg_proc = F_proc-1
      if (Ptopo_couleur == 0) then
         yyg_proc = F_proc+Ptopo_numproc-1
      end if
!
!----------------------------------------------------------------------
!
   return
   end function yyg_proc

end module yyg_param
