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

!**s/r phyrdfile -- Reading file F_fichier_S for the physics package with callback
!                   routine F_read_cb

      subroutine phyrdfile (F_fichier_S, F_read_cb, F_messg_s, F_myproc)
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.inc>
!
      character(len=*) F_fichier_S, F_messg_s
      integer F_myproc
      external F_read_cb

!author
!     M. Desgagne  - January 2014

      integer, parameter :: max_ndim = 1000
      logical found_L
      integer iun,ilir,inbr,status,ierr
      integer dim(max_ndim)
      real, dimension(1) :: dummy
      real, dimension(:), allocatable :: rbuf
!
!-----------------------------------------------------------------
!
      status = 0

      if (F_myproc.eq.0) then

         inquire (FILE=trim(F_fichier_S),EXIST=found_L)

         if (found_L) then
            ilir = wkoffit(trim(F_fichier_S))
            if (  (ilir.eq.1) .or.(ilir.eq.2).or. &
                  (ilir.eq.33).or.(ilir.eq.34) ) then
               write (6,1001) trim(F_messg_s),trim(F_fichier_S)
            else
               print*, ' FILE ',trim(F_fichier_S)
               print*, ' NOT FST FILE FORMAT -- ABORT --'
               status = -1
            endif
         else
            print*
            print *,'********************************************'
            print *,'   CAN NOT FIND FILE: ',trim(F_fichier_S)
            print *,'********************************************'
            status = -1
         endif

      endif

      call handle_error(status,'itf_phy_rdfile','itf_phy_rdfile')

      status = 0
      if (F_myproc.eq.0) then
         iun  = 0
         ilir = fnom    (iun,trim(F_fichier_S),'STD+RND+OLD+R/O',0)
         ilir = fstouv  (iun,'RND')

         status = 200
         call F_read_cb (iun,dummy,dim,status)
         if (status.lt.0) goto 9977
         allocate (rbuf(dim(2)))
         rbuf = 0
         status = 250
         call F_read_cb (iun,rbuf,dim,status)
         inbr = fstfrm  (iun)
         inbr = fclos   (iun)
      endif

9977  call handle_error(status,'itf_phy_rdfile','itf_phy_rdfile')
      call RPN_COMM_bcast (dim,max_ndim,"MPI_INTEGER",0,"grid",ierr)
      if (F_myproc.gt.0) then
         allocate (rbuf(dim(2)))
         rbuf = 0
      endif
      call RPN_COMM_bcast (rbuf,dim(2),"MPI_REAL",0,"grid",ierr)
      status = 400
      call F_read_cb (iun,rbuf,dim,status)
      deallocate (rbuf) 

 9988 call handle_error(status,'itf_phy_rdfile','itf_phy_rdfile')

      inbr = fstopc ('MSGLVL','WARNIN',RMN_OPT_SET)

 1001 format (/'READING ',a,' FILE:'/a)
!
!-----------------------------------------------------------------
!
      return
      end
