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
module out_collector

  use iso_c_binding
  implicit none
#include <arch_specific.hf>


      integer Bloc_npx   , Bloc_npy  , &
              Bloc_nblocs, Bloc_npes , &
              Bloc_mybloc, Bloc_me   , &
              Bloc_maxdim

      integer, dimension(:  ), pointer :: Bloc_peid, Bloc_pe0id
      integer, dimension(:,:), pointer :: Bloc_glbij

  private
  public :: block_collect_set, block_collect_fullp, Bloc_me

contains

      subroutine block_collect_set (F_npex, F_npey)
      use ptopo
      implicit none

      integer F_npex, F_npey

      include "rpn_comm.inc"

      integer, dimension(:,:) , allocatable :: block
      integer i,j,err,cnt,n1
      integer my_row, my_col, this_block
      integer mpx,irest,l_npex,i0,in,l_npey,j0,jn
      integer block_info(0:Ptopo_numproc-1,2),info(0:Ptopo_numproc-1,2)
!
!----------------------------------------------------------------------
!
      my_row= Ptopo_myproc         / Ptopo_npex
      my_col= Ptopo_myproc - my_row* Ptopo_npex

      Bloc_npx= F_npex ; Bloc_npy= F_npey
      Bloc_nblocs= Bloc_npx*Bloc_npy
      allocate (block(0:Bloc_nblocs-1,5))

      do i= 0, Bloc_npx-1

         mpx      = mod( i, Bloc_npx )
         l_npex   = Ptopo_npex / Bloc_npx
         irest    = Ptopo_npex - l_npex * Bloc_npx
         i0       = mpx * l_npex + 1
         if ( mpx < irest ) then
            l_npex = l_npex + 1
            i0     = i0 + mpx
         else
            i0 = i0 + irest
         end if
         i0= i0 - 1 ; in= i0 + l_npex - 1

         do j= 0, Bloc_npy-1
            mpx      = mod( j, Bloc_npy )
            l_npey   = Ptopo_npey / Bloc_npy
            irest    = Ptopo_npey - l_npey * Bloc_npy
            j0       = mpx * l_npey + 1
            if ( mpx < irest ) then
               l_npey = l_npey + 1
               j0     = j0 + mpx
            else
               j0 = j0 + irest
            end if
            j0= j0 - 1 ; jn= j0 + l_npey - 1

            this_block= j*Bloc_npx+i
            block(this_block,1) = l_npex*l_npey
            block(this_block,2) = Ptopo_gindx(1,j0*Ptopo_npex+i0+1)
            block(this_block,3) = Ptopo_gindx(2,jn*Ptopo_npex+in+1)
            block(this_block,4) = Ptopo_gindx(3,j0*Ptopo_npex+i0+1)
            block(this_block,5) = Ptopo_gindx(4,jn*Ptopo_npex+in+1)
            if ((my_col>=i0).and.(my_col<=(i0+l_npex-1)) .and. &
                (my_row>=j0).and.(my_row<=(j0+l_npey-1))) then
                Bloc_mybloc= this_block
                Bloc_me    = (my_row-j0)*l_npex + my_col - i0
            end if

         end do
      end do

      info=0
      info(Ptopo_myproc,1) = Bloc_mybloc
      info(Ptopo_myproc,2) = Bloc_me
      call rpn_comm_ALLREDUCE ( info, block_info, Ptopo_numproc*2, &
                                "MPI_INTEGER","MPI_SUM",'GRID',err )

      Bloc_npes= block(Bloc_mybloc,1)
      allocate (Bloc_peid(   0:Bloc_npes  -1), &
                Bloc_pe0id(  0:Bloc_nblocs-1), &
                Bloc_glbij(4,0:Bloc_nblocs-1) )

      cnt= -1
      do i=0,Ptopo_numproc-1
         if (block_info(i,1) == Bloc_mybloc) then
            cnt= cnt + 1
            Bloc_peid(cnt) = i
         end if
         if (block_info(i,2) == 0) then
            Bloc_pe0id(  block_info(i,1)) = i
            Bloc_glbij(1,block_info(i,1)) = block(block_info(i,1),2)
            Bloc_glbij(2,block_info(i,1)) = block(block_info(i,1),3)
            Bloc_glbij(3,block_info(i,1)) = block(block_info(i,1),4)
            Bloc_glbij(4,block_info(i,1)) = block(block_info(i,1),5)
         end if
      end do

      Bloc_maxdim= 0
      do i= 0, Bloc_nblocs-1
         n1= (Bloc_glbij(2,i)-Bloc_glbij(1,i)+1) * &
             (Bloc_glbij(4,i)-Bloc_glbij(3,i)+1)
         Bloc_maxdim= max(Bloc_maxdim,n1)
      end do
!
!----------------------------------------------------------------------
!
      return
      end subroutine block_collect_set

      subroutine block_collect_fullp ( src, lminx,lmaxx,lminy,lmaxy,Nk,&
                                       dst, nz, zlist )
      implicit none

      integer lminx,lmaxx,lminy,lmaxy,Nk,nz
      integer, dimension(:), pointer :: zlist
      real src(lminx:lmaxx,lminy:lmaxy,Nk)
      real, dimension(:,:,:), pointer :: dst

      real, dimension(:,:,:), pointer :: f2rc
!
!----------------------------------------------------------------------
!
      nullify(f2rc)
      call block_collect ( f2rc, src, lminx,lmaxx,lminy,lmaxy,Nk )
      call block_fplanes ( dst, zlist, nz, f2rc, Nk)
      if (associated(f2rc)) deallocate (f2rc)
!
!----------------------------------------------------------------------
!
      return
      end subroutine block_collect_fullp

      subroutine block_collect (f2rc,f2cc,lminx,lmaxx,lminy,lmaxy,Nk)
      use glb_ld
      use ptopo
      implicit none

      integer lminx,lmaxx,lminy,lmaxy,Nk
      real, dimension(:,:,:), pointer :: f2rc
      real f2cc(lminx:lmaxx,lminy:lmaxy,Nk)

      integer i, j, k, iproc, tag, err, status, ni, nj
      integer l_id, l_jd, len, tag_comm
      real, dimension(:,:,:), allocatable :: buf
!
!----------------------------------------------------------------------
!
      tag = 210

      if (Bloc_me == 0) then

! Copy local data (LD) segment
         allocate (&
            f2rc(Bloc_glbij(1,Bloc_mybloc):Bloc_glbij(2,Bloc_mybloc),&
                 Bloc_glbij(3,Bloc_mybloc):Bloc_glbij(4,Bloc_mybloc),Nk))
         l_id= Ptopo_gindx(1,Bloc_peid(0)+1)
         l_jd= Ptopo_gindx(3,Bloc_peid(0)+1)

         f2rc = 0.
         do k = 1, Nk
            do j = 1, l_nj
               do i = 1, l_ni
                  f2rc(i+l_id-1,j+l_jd-1,k) = f2cc(i,j,k)
               end do
            end do
         end do

! Receive local data (LD) segments from other processors of bloc

         do iproc = 1, Bloc_npes-1
            ni= Ptopo_gindx(2,Bloc_peid(iproc)+1) - &
                Ptopo_gindx(1,Bloc_peid(iproc)+1) + 1
            nj= Ptopo_gindx(4,Bloc_peid(iproc)+1) - &
                Ptopo_gindx(3,Bloc_peid(iproc)+1) + 1
            allocate (buf(ni,nj,Nk))
            len= ni*nj*Nk
            tag_comm= tag+iproc+Bloc_mybloc
            call RPN_COMM_recv ( buf, len, 'MPI_REAL', Bloc_peid(iproc),&
                                       tag_comm, 'GRID', status, err )
            l_id= Ptopo_gindx(1,Bloc_peid(iproc)+1)
            l_jd= Ptopo_gindx(3,Bloc_peid(iproc)+1)
            do k = 1, Nk
               do j = 1, nj
                  do i = 1, ni
                     f2rc(i+l_id-1,j+l_jd-1,k) = buf(i,j,k)
                  end do
               end do
            end do
            deallocate (buf)
         end do

      else

! Send local data (LD) segment to processor 0 of mybloc
         allocate (buf(l_ni,l_nj,Nk))
         buf(1:l_ni,1:l_nj,:) = f2cc(1:l_ni,1:l_nj,:)
         len= l_ni*l_nj*Nk
         tag_comm= tag+Bloc_me+Bloc_mybloc
         call RPN_COMM_send ( buf, len, 'MPI_REAL', &
             Bloc_pe0id(Bloc_mybloc), tag_comm, 'GRID',err )
         deallocate (buf)

      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine block_collect

      subroutine block_fplanes ( glb, zlist, nz, f2rc, Nk)
      use glb_ld
      use ptopo
      implicit none

      integer nz,Nk
      integer, dimension(:    ), pointer :: zlist
      real   , dimension(:,:,:), pointer :: glb, f2rc

      integer mpx,kstart,local_nk,irest,len
      real, dimension(:,:), allocatable :: buf

      integer i,j,k,b,ni,nj,nkk,err
      integer ireq, tag, tag_comm, request(Bloc_nblocs*2)
      integer blk_dist(2,0:Bloc_nblocs-1)
!
!----------------------------------------------------------------------
!
      nz= 0 ; ireq= 0 ; tag= 23

      if (Bloc_me == 0) then

! Distributing the work: Nk onto Bloc_nblocs
         do i= 0, Bloc_nblocs-1
            mpx      = mod( i, Bloc_nblocs )
            local_nk = Nk / Bloc_nblocs
            irest  = Nk  - local_nk * Bloc_nblocs
            kstart = mpx * local_nk + 1
            if ( mpx < irest ) then
               local_nk   = local_nk + 1
               kstart = kstart + mpx
            else
               kstart = kstart + irest
            end if
            blk_dist(1,i) = kstart
            blk_dist(2,i) = kstart+local_nk-1
         end do

         ni= Bloc_glbij(2,Bloc_mybloc) - Bloc_glbij(1,Bloc_mybloc) + 1
         nj= Bloc_glbij(4,Bloc_mybloc) - Bloc_glbij(3,Bloc_mybloc) + 1
! Sending
         do i= 0, Bloc_nblocs-1
            if (i /= Bloc_mybloc) then
            if (blk_dist(1,i) <= Nk) then
               nkk=(blk_dist(2,i)-blk_dist(1,i)+1)
               len= ni*nj*nkk
               ireq = ireq+1 ; tag_comm = tag+Ptopo_myproc
               call RPN_COMM_isend ( &
                       f2rc(Bloc_glbij(1,Bloc_mybloc)               ,&
                            Bloc_glbij(3,Bloc_mybloc),blk_dist(1,i)),&
                            len, 'MPI_REAL', Bloc_pe0id(i), tag_comm,&
                                          'GRID', request(ireq), err )
            end if
            end if
         end do

! Receiving
         if (blk_dist(1,Bloc_mybloc) <= Nk) then
            nkk= (blk_dist(2,Bloc_mybloc)-blk_dist(1,Bloc_mybloc)+1)
            nz = nkk
            allocate (zlist(nz), buf(Bloc_maxdim*nkk,0:Bloc_nblocs-1))
            buf = 0.
            do i=1,nz
               zlist(i) = blk_dist(1,Bloc_mybloc) + i - 1
            end do
            allocate (glb(G_ni,G_nj,nz))

            do k=blk_dist(1,Bloc_mybloc),blk_dist(2,Bloc_mybloc)
               do j=Bloc_glbij(3,Bloc_mybloc),Bloc_glbij(4,Bloc_mybloc)
                  do i=Bloc_glbij(1,Bloc_mybloc),Bloc_glbij(2,Bloc_mybloc)
                     glb(i,j,(k-blk_dist(1,Bloc_mybloc))+1)=f2rc(i,j,k)
                  end do
               end do
            end do

            do i= 0, Bloc_nblocs-1
               if (i /= Bloc_mybloc) then
                  ni= Bloc_glbij(2,i) - Bloc_glbij(1,i) + 1
                  nj= Bloc_glbij(4,i) - Bloc_glbij(3,i) + 1
                  len= ni*nj*nkk
                  ireq = ireq+1 ; tag_comm = tag + Bloc_pe0id(i)
                  call RPN_COMM_irecv ( &
                               buf(1,i), len, 'MPI_REAL', Bloc_pe0id(i),&
                               tag_comm, 'GRID', request(ireq), err )
               end if
            end do
         end if

      end if

      call RPN_COMM_waitall_nostat (ireq, request, err)

! Filling
      if ((Bloc_me == 0) .and. (blk_dist(1,Bloc_mybloc) <= Nk)) then
         do b= 0, Bloc_nblocs-1
            if (b /= Bloc_mybloc) then
               len= 0
               do k=1,nz
                  do j=Bloc_glbij(3,b),Bloc_glbij(4,b)
                     do i=Bloc_glbij(1,b),Bloc_glbij(2,b)
                        len= len+1
                        glb(i,j,k) = buf(len,b)
                     end do
                  end do
               end do
            end if
         end do
         deallocate (buf)
      end if
!
!----------------------------------------------------------------------
!
      return
      end subroutine block_fplanes

end module out_collector
