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
!*** s/r yyg_rhs_initscalbc - to initialize communication pattern for RHS data
!


      Subroutine yyg_rhs_initscalbc()
      use gem_options
      use glb_ld
      use glb_pil
      use ptopo
      use sol
      use tdpack
      use yyg_rhs
      use, intrinsic :: iso_fortran_env
      implicit none
#include <arch_specific.hf>
!
!author
!           Abdessamad Qaddouri/V.Lee - October 2009

      integer i,j,imx,imy,kk,ii,jj,stencil
      integer ki,ksend,krecv,maxsendlen_per_proc
      integer, dimension (:), pointer :: recv_len,send_len
      real(kind=REAL64)  xg_8(1-G_ni:2*G_ni),yg_8(1-G_nj:2*G_nj)
      real(kind=REAL64)  xx_8(G_ni,G_nj),yy_8(G_ni,G_nj)
      real(kind=REAL64)  s(2,2),h1,h2
      real(kind=REAL64)  x_d,y_d,x_a,y_a
      real(kind=REAL64) TWO_8
      parameter( TWO_8   = 2.0d0 )
!
!     Localise could get point way outside of the actual grid in search
!     So extend all global arrays: xg_8,yg_8
!

      do i=1,G_ni
         xg_8(i) = G_xg_8(i)
      end do
      do j=1,G_nj
         yg_8(j) = G_yg_8(j)
      end do

      do i=-G_ni+1,0
         xg_8(i) = xg_8(i+G_ni) - TWO_8*pi_8
      end do
      do i=G_ni+1,2*G_ni
         xg_8(i) = xg_8(i-G_ni) + TWO_8*pi_8
      end do

      yg_8( 0    ) = -(yg_8(1) + pi_8)
      yg_8(-1    ) = -TWO_8*pi_8 -  &
           (yg_8(0)+yg_8(1)+yg_8(2))
      yg_8(G_nj+1) =  pi_8 - yg_8(G_nj)
      yg_8(G_nj+2) =  TWO_8*pi_8 - &
           (yg_8(G_nj+1)+yg_8(G_nj)+yg_8(G_nj-1))
      do j=-2,-G_nj+1,-1
         yg_8(j) = 1.01*yg_8(j+1)
      end do
      do j=G_nj+3,2*G_nj
         yg_8(j) = 1.01*yg_8(j-1)
      end do

!
      do j=1,G_nj
      do i=1,G_ni
         xx_8(i,j)=xg_8(i)
      end do
      end do
      do j=1,G_nj
      do i=1,G_ni
         yy_8(i,j)=yg_8(j)
      end do
      end do

!Delta xg, yg is not identical between xg(i) and xg(i+1)
!h1, h2 used in this routine is ok as it is a close estimate for
!creating YY pattern exchange and it works on the global tile

      h1=xg_8(2)-xg_8(1)
      h2=yg_8(2)-yg_8(1)
!
! And allocate temp vectors needed for counting for each processor
!
      allocate (recv_len (Ptopo_numproc))
      allocate (send_len (Ptopo_numproc))
      recv_len (:)=0
      send_len (:)=0
!
!
! FIRST PASS is to find the number of processor to tag for
! communication and the number of items to send and receive for each
! processor before allocating the vectors
!
! WEST section

      do j=1+glb_pil_s,G_nj-glb_pil_n
         i=glb_pil_w
         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         call localise(imx,imy,x_a,y_a, &
                       xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)
! check to collect from who
         if (i+1 >= l_i0.and.i+1 <= l_i0+l_ni-1 .and. &
             j >= l_j0.and.j <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (imx >= Ptopo_gindx(1,kk).and.imx <= Ptopo_gindx(2,kk).and. &
                    imy >= Ptopo_gindx(3,kk).and.imy <= Ptopo_gindx(4,kk))then
                    recv_len(kk)=recv_len(kk)+1
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (i+1 >= Ptopo_gindx(1,kk).and.i+1 <= Ptopo_gindx(2,kk).and. &
                    j >= Ptopo_gindx(3,kk).and.j <= Ptopo_gindx(4,kk))then
                    send_len(kk)=send_len(kk)+1
                end if
             end do
         end if
      end do
!
!
! East section
      do j=1+glb_pil_s,G_nj-glb_pil_n
         i=G_ni-glb_pil_e+1
         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         call localise(imx,imy,x_a,y_a, &
                          xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)
! check to collect from who
         if (i-1 >= l_i0.and.i-1 <= l_i0+l_ni-1 .and. &
             j >= l_j0.and.j <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (imx >= Ptopo_gindx(1,kk).and.imx <= Ptopo_gindx(2,kk).and. &
                    imy >= Ptopo_gindx(3,kk).and.imy <= Ptopo_gindx(4,kk))then
                    recv_len(kk)=recv_len(kk)+1
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (i-1 >= Ptopo_gindx(1,kk).and.i-1 <= Ptopo_gindx(2,kk).and. &
                    j >= Ptopo_gindx(3,kk).and.j <= Ptopo_gindx(4,kk))then
                    send_len(kk)=send_len(kk)+1
                end if
             end do
         end if
      end do
!
! South section
      j=glb_pil_s
      do i=1+glb_pil_w,G_ni-glb_pil_e

         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         call localise(imx,imy,x_a,y_a, &
                          xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)
! check to collect from who
         if (i >= l_i0.and.i <= l_i0+l_ni-1 .and. &
             j+1 >= l_j0.and.j+1 <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (imx >= Ptopo_gindx(1,kk).and.imx <= Ptopo_gindx(2,kk).and. &
                    imy >= Ptopo_gindx(3,kk).and.imy <= Ptopo_gindx(4,kk))then
                    recv_len(kk)=recv_len(kk)+1
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (i >= Ptopo_gindx(1,kk).and.i <= Ptopo_gindx(2,kk).and. &
                    j+1 >= Ptopo_gindx(3,kk).and.j+1 <= Ptopo_gindx(4,kk))then
                    send_len(kk)=send_len(kk)+1
                end if
             end do
         end if
      end do
!
! North section
      j=G_nj-glb_pil_n+1
      do i=1+glb_pil_w,G_ni-glb_pil_e

         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         call localise(imx,imy,x_a,y_a, &
                          xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)

! check to collect from who
         if (i >= l_i0.and.i <= l_i0+l_ni-1 .and. &
             j-1 >= l_j0.and.j-1 <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (imx >= Ptopo_gindx(1,kk).and.imx <= Ptopo_gindx(2,kk).and. &
                    imy >= Ptopo_gindx(3,kk).and.imy <= Ptopo_gindx(4,kk))then
                    recv_len(kk)=recv_len(kk)+1
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Ptopo_numproc
                if (i >= Ptopo_gindx(1,kk).and.i <= Ptopo_gindx(2,kk).and. &
                    j-1 >= Ptopo_gindx(3,kk).and.j-1 <= Ptopo_gindx(4,kk))then
                    send_len(kk)=send_len(kk)+1
                end if
             end do
         end if
      end do
!
! Obtain sum of elements to send and receive for each processor
! and the total memory needed to store and receive for each processor
!
     Rhsx_send_all=0
     Rhsx_recv_all=0
     Rhsx_sendmaxproc=0
     Rhsx_recvmaxproc=0
     maxsendlen_per_proc=0

     do kk=1,Ptopo_numproc
        maxsendlen_per_proc = max(maxsendlen_per_proc,send_len(kk) )
        Rhsx_send_all=send_len(kk)+Rhsx_send_all
        Rhsx_recv_all=recv_len(kk)+Rhsx_recv_all

        if (send_len(kk) > 0) Rhsx_sendmaxproc=Rhsx_sendmaxproc+1
        if (recv_len(kk) > 0) Rhsx_recvmaxproc=Rhsx_recvmaxproc+1
     end do
!
!     print *,'Allocate common vectors'

      allocate (Sol_rhs(maxsendlen_per_proc*l_nk,1,Rhsx_sendmaxproc))
      Sol_rhs = 0.0
      allocate (Rhsx_recvproc(Rhsx_recvmaxproc))
      allocate (Rhsx_recv_len(Rhsx_recvmaxproc))
      allocate (Rhsx_recv_adr(Rhsx_recvmaxproc))

      allocate (Rhsx_sendproc(Rhsx_sendmaxproc))
      allocate (Rhsx_send_len(Rhsx_sendmaxproc))
      allocate (Rhsx_send_adr(Rhsx_sendmaxproc))
      Rhsx_recv_len(:) = 0
      Rhsx_send_len(:) = 0
      Rhsx_recv_adr(:) = 0
      Rhsx_send_adr(:) = 0


     ksend=0
     krecv=0
     Rhsx_send_all=0
     Rhsx_recv_all=0
!
! Fill the lengths and addresses for selected processors to communicate
!
     do kk=1,Ptopo_numproc
        if (send_len(kk) > 0) then
            ksend=ksend+1
            Rhsx_sendproc(ksend)=kk
            Rhsx_send_len(ksend)=send_len(kk)

            Rhsx_send_adr(ksend)= Rhsx_send_all
            Rhsx_send_all= Rhsx_send_all + Rhsx_send_len(ksend)
        end if
        if (recv_len(kk) > 0) then
            krecv=krecv+1
            Rhsx_recvproc(krecv)=kk
            Rhsx_recv_len(krecv)=recv_len(kk)

            Rhsx_recv_adr(krecv)= Rhsx_recv_all
            Rhsx_recv_all= Rhsx_recv_all + Rhsx_recv_len(krecv)
        end if

     end do
!    print *,'krecv=',krecv,'Rhsx_recvmaxproc=',Rhsx_recvmaxproc
!    print *,'ksend=',ksend,'Rhsx_sendmaxproc=',Rhsx_sendmaxproc

!     print *,'Summary of RHS comm procs'
!     do kk=1,Rhsx_recvmaxproc
!       print *,'From R proc:',Rhsx_recvproc(kk),'Rhsx_recv_len',Rhsx_recv_len(kk),'adr',Rhsx_recv_adr(kk)
!     end do
!     do kk=1,Rhsx_sendmaxproc
!       print *,'To R proc:',Rhsx_sendproc(kk),'Rhsx_send_len',Rhsx_send_len(kk),'adr',Rhsx_send_adr(kk)
!     end do

!
! Now allocate the vectors needed for sending and receiving each processor
!     print *,'yyg_rhs_initscalbc: Rhsx_recv_all=',Rhsx_recv_all
      if (Rhsx_recv_all > 0) then
          allocate (Rhsx_recv_i(Rhsx_recv_all))
          allocate (Rhsx_recv_j(Rhsx_recv_all))
          Rhsx_recv_i(:) = 0
          Rhsx_recv_j(:) = 0
      end if


!     print *,'yyg_rhs_initscalbc: Rhsx_send_all=',Rhsx_send_all
      if (Rhsx_send_all > 0) then
          allocate (Rhsx_send_imx(Rhsx_send_all))
          allocate (Rhsx_send_imy(Rhsx_send_all))
          allocate (Rhsx_send_xxr(Rhsx_send_all))
          allocate (Rhsx_send_yyr(Rhsx_send_all))
          allocate (Rhsx_send_sten(Rhsx_send_all))
          Rhsx_send_imx(:) = 0
          Rhsx_send_imy(:) = 0
          Rhsx_send_xxr(:) = 0.0
          Rhsx_send_yyr(:) = 0.0
          Rhsx_send_sten(:) = 0.0
      end if
!

      recv_len(:)=0
      send_len(:)=0
!
! SECOND PASS is to initialize the vectors with information for communication
!
! WEST section

      do j=1+glb_pil_s,G_nj-glb_pil_n
         stencil=j-glb_pil_s
         i=glb_pil_w
         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         call localise(imx,imy,x_a,y_a, &
                       xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)
! check to collect from who
         if (i+1 >= l_i0.and.i+1 <= l_i0+l_ni-1 .and. &
             j >= l_j0.and.j <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_recvmaxproc
                ki=Rhsx_recvproc(kk)
                if (imx >= Ptopo_gindx(1,ki).and.imx <= Ptopo_gindx(2,ki).and. &
                    imy >= Ptopo_gindx(3,ki).and.imy <= Ptopo_gindx(4,ki))then
                    recv_len(kk)=recv_len(kk)+1
                    ii=i-l_i0+2
                    jj=j-l_j0+1
                    Rhsx_recv_i(Rhsx_recv_adr(kk)+recv_len(kk))=ii
                    Rhsx_recv_j(Rhsx_recv_adr(kk)+recv_len(kk))=jj
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_sendmaxproc
                ki=Rhsx_sendproc(kk)
                if (i+1 >= Ptopo_gindx(1,ki).and.i+1 <= Ptopo_gindx(2,ki).and. &
                    j >= Ptopo_gindx(3,ki).and.j <= Ptopo_gindx(4,ki))then
                    send_len(kk)=send_len(kk)+1
                    Rhsx_send_imx(Rhsx_send_adr(kk)+send_len(kk))=imx-l_i0+1
                    Rhsx_send_imy(Rhsx_send_adr(kk)+send_len(kk))=imy-l_j0+1
                    Rhsx_send_xxr(Rhsx_send_adr(kk)+send_len(kk))=x_a
                    Rhsx_send_yyr(Rhsx_send_adr(kk)+send_len(kk))=y_a
                    Rhsx_send_sten(Rhsx_send_adr(kk)+send_len(kk))=Sol_stencil2_8(stencil)
                end if
             end do
         end if
      end do
!
!
! East section
      do j=1+glb_pil_s,G_nj-glb_pil_n
         i=G_ni-glb_pil_e+1
         stencil=j-glb_pil_s
         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         y_a= y_a
         call localise(imx,imy,x_a,y_a, &
                          xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)
! check to collect from who
         if (i-1 >= l_i0.and.i-1 <= l_i0+l_ni-1 .and. &
             j >= l_j0.and.j <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_recvmaxproc
                ki=Rhsx_recvproc(kk)
                if (imx >= Ptopo_gindx(1,ki).and.imx <= Ptopo_gindx(2,ki).and. &
                    imy >= Ptopo_gindx(3,ki).and.imy <= Ptopo_gindx(4,ki))then
                    recv_len(kk)=recv_len(kk)+1
                    ii=i-l_i0
                    jj=j-l_j0+1
                    Rhsx_recv_i(Rhsx_recv_adr(kk)+recv_len(kk))=ii
                    Rhsx_recv_j(Rhsx_recv_adr(kk)+recv_len(kk))=jj
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_sendmaxproc
                ki=Rhsx_sendproc(kk)
                if (i-1 >= Ptopo_gindx(1,ki).and.i-1 <= Ptopo_gindx(2,ki).and. &
                    j >= Ptopo_gindx(3,ki).and.j <= Ptopo_gindx(4,ki))then
                    send_len(kk)=send_len(kk)+1
                    Rhsx_send_imx(Rhsx_send_adr(kk)+send_len(kk))=imx-l_i0+1
                    Rhsx_send_imy(Rhsx_send_adr(kk)+send_len(kk))=imy-l_j0+1
                    Rhsx_send_xxr(Rhsx_send_adr(kk)+send_len(kk))=x_a
                    Rhsx_send_yyr(Rhsx_send_adr(kk)+send_len(kk))=y_a
                    Rhsx_send_sten(Rhsx_send_adr(kk)+send_len(kk))=Sol_stencil3_8(stencil)
                end if
             end do
         end if
      end do
!
! South section
      j=glb_pil_s
      do i=1+glb_pil_w,G_ni-glb_pil_e
         stencil=i-glb_pil_w

         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         y_a= y_a
         call localise(imx,imy,x_a,y_a, &
                          xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)
! check to collect from who
         if (i >= l_i0.and.i <= l_i0+l_ni-1 .and. &
             j+1 >= l_j0.and.j+1 <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_recvmaxproc
                ki=Rhsx_recvproc(kk)
                if (imx >= Ptopo_gindx(1,ki).and.imx <= Ptopo_gindx(2,ki).and. &
                    imy >= Ptopo_gindx(3,ki).and.imy <= Ptopo_gindx(4,ki))then
                    recv_len(kk)=recv_len(kk)+1
                    ii=i-l_i0+1
                    jj=j-l_j0+2
                    Rhsx_recv_i(Rhsx_recv_adr(kk)+recv_len(kk))=ii
                    Rhsx_recv_j(Rhsx_recv_adr(kk)+recv_len(kk))=jj
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_sendmaxproc
                ki=Rhsx_sendproc(kk)
                if (i >= Ptopo_gindx(1,ki).and.i <= Ptopo_gindx(2,ki).and. &
                    j+1 >= Ptopo_gindx(3,ki).and.j+1 <= Ptopo_gindx(4,ki))then
                    send_len(kk)=send_len(kk)+1
                    Rhsx_send_imx(Rhsx_send_adr(kk)+send_len(kk))=imx-l_i0+1
                    Rhsx_send_imy(Rhsx_send_adr(kk)+send_len(kk))=imy-l_j0+1
                    Rhsx_send_xxr(Rhsx_send_adr(kk)+send_len(kk))=x_a
                    Rhsx_send_yyr(Rhsx_send_adr(kk)+send_len(kk))=y_a
                    Rhsx_send_sten(Rhsx_send_adr(kk)+send_len(kk))=Sol_stencil4_8(stencil)
                end if
             end do
         end if
      end do
!
! North section
      j=G_nj-glb_pil_n+1
      do i=1+glb_pil_w,G_ni-glb_pil_e
         stencil = i-glb_pil_w

         x_d=xx_8(i,j)-pi_8
         y_d=yy_8(i,j)
         call smat(s,x_a,y_a,x_d,y_d)
         x_a=x_a+pi_8
         y_a= y_a
         call localise(imx,imy,x_a,y_a, &
                          xg_8(1),yg_8(1),h1,h2,1,1)
         imx = min(max(imx-1,1),G_ni-3)
         imy = min(max(imy-1,1),G_nj-3)

! check to collect from who
         if (i >= l_i0.and.i <= l_i0+l_ni-1 .and. &
             j-1 >= l_j0.and.j-1 <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_recvmaxproc
                ki=Rhsx_recvproc(kk)
                if (imx >= Ptopo_gindx(1,ki).and.imx <= Ptopo_gindx(2,ki).and. &
                    imy >= Ptopo_gindx(3,ki).and.imy <= Ptopo_gindx(4,ki))then
                    recv_len(kk)=recv_len(kk)+1
                    ii=i-l_i0+1
                    jj=j-l_j0
                    Rhsx_recv_i(Rhsx_recv_adr(kk)+recv_len(kk))=ii
                    Rhsx_recv_j(Rhsx_recv_adr(kk)+recv_len(kk))=jj
                end if
             end do
         end if

! check to send to who
         if (imx >= l_i0.and.imx <= l_i0+l_ni-1 .and. &
             imy >= l_j0.and.imy <= l_j0+l_nj-1      ) then
             do kk=1,Rhsx_sendmaxproc
                ki=Rhsx_sendproc(kk)
                if (i >= Ptopo_gindx(1,ki).and.i <= Ptopo_gindx(2,ki).and. &
                    j-1 >= Ptopo_gindx(3,ki).and.j-1 <= Ptopo_gindx(4,ki))then
                    send_len(kk)=send_len(kk)+1
                    Rhsx_send_imx(Rhsx_send_adr(kk)+send_len(kk))=imx-l_i0+1
                    Rhsx_send_imy(Rhsx_send_adr(kk)+send_len(kk))=imy-l_j0+1
                    Rhsx_send_xxr(Rhsx_send_adr(kk)+send_len(kk))=x_a
                    Rhsx_send_yyr(Rhsx_send_adr(kk)+send_len(kk))=y_a
                    Rhsx_send_sten(Rhsx_send_adr(kk)+send_len(kk))=Sol_stencil5_8(stencil)
                end if
             end do
         end if
      end do
!Check receive lengths from each processor
!     do ki=1,Rhsx_recvmaxproc
!        kk=Rhsx_recvproc(ki)
!        if (Ptopo_couleur == 0) then
!            kkproc = kk+Ptopo_numproc-1
!        else
!            kkproc = kk -1
!        end if
!    write(output_unit,1000) 'Rhsx_recv_len',kkproc,Rhsx_recv_len(kk),Rhsx_recv_adr(kk)
!   end do
!Check send lengths to each processor

!     do ki=1,Rhsx_sendmaxproc
!        kk=Rhsx_sendproc(ki)
!        if (Ptopo_couleur == 0) then
!            kkproc = kk+Ptopo_numproc-1
!        else
!            kkproc = kk -1
!        end if
! write(output_unit,1000) 'Rhsx_send_len',kkproc,Rhsx_send_len(kk),Rhsx_send_adr(kk)
!     end do
      deallocate (recv_len,send_len)


 1000 format(a15,i3,'=',i5,'bytes, addr=',i5)
 1001 format(a15,i3,'=',i4,'bytes   i:', i3,' j:',i3)
 1002 format(a15,i3,'=',i4,'bytes   i:', i3,' j:',i3,'sten=',i3)

      return
      end

