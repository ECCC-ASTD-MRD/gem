!/! RPN_COMM - Library of useful routines for C and FORTRAN programming
! ! Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! !                          Environnement Canada
! !
! ! This library is free software; you can redistribute it and/or
! ! modify it under the terms of the GNU Lesser General Public
! ! License as published by the Free Software Foundation,
! ! version 2.1 of the License.
! !
! ! This library is distributed in the hope that it will be useful,
! ! but WITHOUT ANY WARRANTY; without even the implied warranty of
! ! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! ! Lesser General Public License for more details.
! !
! ! You should have received a copy of the GNU Lesser General Public
! ! License along with this library; if not, write to the
! ! Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! ! Boston, MA 02111-1307, USA.
! !/

      subroutine RPN_COMM_globalsum(array,minx,maxx,miny,maxy,&
     &       nptsz,nil,njl,njlmax,gnj,sum)
      use rpn_comm
      implicit none

!      include 'rpn_comm.h'
!      include 'mpif.h'

      integer minx,maxx,miny,maxy,nptsz,nil,njl
      integer gnj,njlmax
      real array(minx:maxx,miny:maxy,nptsz)
      integer i,j,k,ierr,status(mpi_status_size)
      real sum(nptsz)
      real vsum(nptsz,njlmax)
      real vlsum(nptsz*max(0,pe_me-(pe_tot-2)),pe_ny*njlmax)

      do k=1,nptsz
          sum(k)=0.0
      enddo

      do j=1,njl
      do k=1,nptsz
         vsum(k,j)=0.0
      enddo
      enddo

      if(pe_tot.eq.1) then

         do k=1,nptsz
         do i=1,nil
         do j=1,njl
            vsum(k,j)=vsum(k,j)+array(i,j,k)
         enddo
         enddo
         enddo
         do j=1,njl
         do k=1,nptsz
            sum(k)=sum(k)+vsum(k,j)
         enddo
         enddo

         return
      endif

      if(pe_mex.eq.0)then

         do k=1,nptsz
         do i=1,nil
         do j=1,njl
            vsum(k,j)=vsum(k,j)+array(i,j,k)
         enddo
         enddo
         enddo

         if (pe_mex.eq.pe_nx-1) then

            call MPI_Gather(vsum,njlmax*nptsz,mpi_real,vlsum,&
     &         njlmax*nptsz,mpi_real,pe_ny-1,pe_mycol,ierr)
         
            if(pe_me.eq.pe_id(pe_nx-1,pe_ny-1)) then
               do j=1,gnj
               do k=1,nptsz
                  sum(k)=sum(k)+vlsum(k,j)
               enddo
               enddo
            endif
         else
            call mpi_send(vsum,njlmax*nptsz,mpi_real,1,pe_mex&
     &                   ,pe_myrow,ierr)
         endif
      else if(pe_mex.eq.pe_nx-1) then

         call mpi_recv(vsum,njlmax*nptsz,mpi_real,pe_mex-1,pe_mex-1&
     &          ,pe_myrow,status,ierr)
         do k=1,nptsz
         do i=1,nil
         do j=1,njl
            vsum(k,j)=vsum(k,j)+array(i,j,k)
         enddo
         enddo
         enddo
         if(pe_ny.eq.1) then
             do j=1,njl
             do k=1,nptsz
                sum(k)=sum(k)+vsum(k,j)
             enddo
             enddo
         else

            call MPI_Gather(vsum,njlmax*nptsz,mpi_real,vlsum&
     &         ,njlmax*nptsz,mpi_real,pe_ny-1,pe_mycol,ierr)
         
            if(pe_me.eq.pe_id(pe_nx-1,pe_ny-1)) then
               do j=1,gnj
               do k=1,nptsz
                  sum(k)=sum(k)+vlsum(k,j)
               enddo
               enddo
            endif
         endif
         
      else
      
         call mpi_recv(vsum,njlmax*nptsz,mpi_real,pe_mex-1,pe_mex-1&
     &      ,pe_myrow,status,ierr )
         
         do k=1,nptsz
         do i=1,nil
         do j=1,njl
            vsum(k,j)=vsum(k,j)+array(i,j,k)
         enddo
         enddo
         enddo

         call mpi_send(vsum,njlmax*nptsz,mpi_real,pe_mex+1,pe_mex,&
     &      pe_myrow,ierr)

      endif

      call mpi_bcast(sum,nptsz,mpi_real,pe_id(pe_nx-1,pe_ny-1),&
     &         pe_defcomm,ierr)

      return
      end
