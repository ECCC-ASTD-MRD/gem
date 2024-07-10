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

module phyrdfile
   implicit none
   private
   public :: phyrdfile1

   integer, parameter, public :: READRAD = 1
   integer, parameter, public :: READOZO = 2

contains
   
   subroutine phyrdfile1(F_fichier_S, F_rad_oz, F_messg_s, F_myproc)
      use rmn_fst24
      use rd_ozone, only: rd_ozone1
      use rd_radtab, only: rd_radtab1
      implicit none
!!!#include <arch_specific.hf>
#include <rmnlib_basics.inc>
!
      character(len=*), intent(in) :: F_fichier_S, F_messg_s
      integer, intent(in) :: F_rad_oz
      integer, intent(in) :: F_myproc

!author
!     M. Desgagne  - January 2014

      integer, parameter :: max_ndim = 1000
      logical found_L
      integer ilir,inbr,status,ierr
      integer dim(max_ndim)
      real, dimension(1) :: dummy
      real, dimension(:), allocatable :: rbuf

      type(fst_file) :: file
      logical        :: success
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

         success = file%open(trim(F_fichier_S),'STD+RND+OLD+R/O')

         status = 200
         if (F_rad_oz == READRAD) then
            call rd_radtab1(file,dummy,dim,status)
         else if (F_rad_oz == READOZO) then
            call rd_ozone1(file,dummy,dim,status)
         endif
         if (status.lt.0) goto 9977
         
         allocate (rbuf(dim(2)))
         rbuf = 0
         status = 250
         if (F_rad_oz == READRAD) then
            call rd_radtab1(file,rbuf,dim,status)
         else if (F_rad_oz == READOZO) then
            call rd_ozone1(file,rbuf,dim,status)
         endif
         success = file%close()
      endif

9977  call handle_error(status,'itf_phy_rdfile','itf_phy_rdfile')
      call RPN_COMM_bcast (dim,max_ndim,"MPI_INTEGER",0,"grid",ierr)
      if (F_myproc.gt.0) then
         allocate (rbuf(dim(2)))
         rbuf = 0
      endif
      call RPN_COMM_bcast (rbuf,dim(2),"MPI_REAL",0,"grid",ierr)
      
      status = 400
      if (F_rad_oz == READRAD) then
         call rd_radtab1(file,rbuf,dim,status)
      else if (F_rad_oz == READOZO) then
         call rd_ozone1(file,rbuf,dim,status)
      endif
      deallocate (rbuf) 

 9988 call handle_error(status,'itf_phy_rdfile','itf_phy_rdfile')

      inbr = fstopc ('MSGLVL','WARNIN',RMN_OPT_SET)

 1001 format (/'READING ',a,' FILE:'/a)
!
!-----------------------------------------------------------------
!
      return
   end subroutine phyrdfile1

end module phyrdfile
