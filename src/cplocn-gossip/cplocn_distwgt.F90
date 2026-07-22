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

!**s/r cplocn_distwgt

      subroutine cplocn_distwgt(F_weightfile_S,F_print_L,F_unout)
use iso_c_binding
use rmn_gmm
use, intrinsic :: iso_fortran_env

      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in) :: F_weightfile_S
      logical, intent(in)          :: F_print_L
      integer, intent(in)          :: F_unout

#include <rmn/msg.h>

      include "thermoconsts.inc"
      include "cpl.cdk"
      include "cplocn.cdk"
!
! v4_7  - Roy F. - initial version
!

! Purpose: Distributes interpolation weights for the coupling
!          with the ocean, and compress global arrays of
!          ocean points, wet points used for interpolation
!          are mapped into a 1-D array

      integer unout, errcode
      logical print_L

! First call global arrays ( de-allocated after initialization )
      integer,dimension (:,:,:),   pointer :: ocn_iwgt_g ! Ice/ocean model i weight indices
      integer,dimension (:,:,:),   pointer :: ocn_jwgt_g ! Ice/ocean model j weight indices
      integer,dimension (:,:),     pointer :: ocn_mwgt_g ! Mask coherent with weights
      real*8, dimension (:,:,:),   pointer :: ocn_wgt_g  ! Weights in the world of GEM
      real,   dimension (:,:),     pointer :: ocn_awgt_g ! Ice/ocean grid angle

      integer,dimension (:,:),     pointer :: atm_ijtoc  ! Compressed index sent to ocn

      real*8, dimension (cpl_minx:cpl_maxx,cpl_miny:cpl_maxy) :: wk8
      real,   dimension (cpl_minx:cpl_maxx,cpl_miny:cpl_maxy) :: wk4
      integer,dimension (cpl_minx:cpl_maxx,cpl_miny:cpl_maxy) :: wki

      integer ivar,ni,nj,nw,gni,gnj,i,j,k,k2,ic,isp,jsp,msk
      integer iproc, status, ier, gni0, gniE, gnj0, gnjE
      integer, external :: msg_getUnit
      logical found
      real*8, parameter :: epsw = 1.e-10
!     ________________________________________________________________

      errcode = -1

      ni=cpl_drv_lni
      nj=cpl_drv_lnj

      gni=cpl_drv_gni
      gnj=cpl_drv_gnj

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)

      allocate(ocn_ijtoc(cplocn_gni,cplocn_gnj))

      if (cplocn_myproc == 0) then

         nw=-1
         call cplao_cnte_w(trim(F_weightfile_S),gni,gnj,nw)
         if (nw <= 0) then
            write (F_unout,9900) ' PROBLEM WITH OCEAN COUPLING WEIGHTS COUNTING'
            call flush(F_unout)
            errcode=-1
            goto 998
         endif
         allocate (ocn_wgt_g (gni,gnj,nw) )
         allocate (ocn_iwgt_g(gni,gnj,nw) )
         allocate (ocn_jwgt_g(gni,gnj,nw) )
         allocate (ocn_mwgt_g(gni,gnj             ) )
         allocate (ocn_awgt_g(gni,gnj             ) )
         call cplao_read_w(trim(F_weightfile_S),                                   &
                           ocn_iwgt_g,ocn_jwgt_g,ocn_mwgt_g,ocn_awgt_g,ocn_wgt_g,  &
                           gni,gnj,nw,errcode)
         if (errcode < 0) then
            write (F_unout,9900) ' PROBLEM WITH OCEAN COUPLING WEIGHTS READING'
            call flush(F_unout)
            goto 998
         endif

! Compression of sent coupled global fields
         gni0=cplocn_dynphy_offset+1
         gnj0=cplocn_dynphy_offset+1
         gniE=gni-cplocn_dynphy_offset
         gnjE=gnj-cplocn_dynphy_offset
         cplocn_atm_nxchng=0
         do j=1,gnj
         do i=1,gni
            msk=0
            do jsp=max(j-cplocn_atm_nspread,1),min(j+cplocn_atm_nspread,gnj)
            do isp=max(i-cplocn_atm_nspread,1),min(i+cplocn_atm_nspread,gni)
               if (ocn_mwgt_g(isp,jsp) > 0) msk=1
            enddo
            enddo
            if ( msk == 1 .and.                  &
                 i >= gni0 .and. i <= gniE .and. &
                 j >= gnj0 .and. j <= gnjE ) cplocn_atm_nxchng = cplocn_atm_nxchng + 1
         enddo
         enddo
         allocate(atm_ijtoc(gni,gnj),atm_ctoi(cplocn_atm_nxchng), &
                                     atm_ctoj(cplocn_atm_nxchng))
         ic=0
         atm_ijtoc(:,:)=-9
         do j=1,gnj
         do i=1,gni
            msk=0
            do jsp=max(j-cplocn_atm_nspread,1),min(j+cplocn_atm_nspread,gnj)
            do isp=max(i-cplocn_atm_nspread,1),min(i+cplocn_atm_nspread,gni)
               if (ocn_mwgt_g(isp,jsp) > 0) msk=1
            enddo
            enddo
            if ( msk == 1 .and.                  &
                 i >= gni0 .and. i <= gniE .and. &
                 j >= gnj0 .and. j <= gnjE ) then
               ic=ic+1
               atm_ijtoc(i,j)=ic
               atm_ctoi(ic)=i
               atm_ctoj(ic)=j
            endif
         enddo
         enddo

         call cplao_xchng_cmap(atm_ijtoc,gni,gnj,                &
                                cplocn_atm_nxchng,               &
                                ocn_ijtoc,cplocn_gni,cplocn_gnj, &
                                cplocn_nxchng,errcode)

         if (errcode < 0) then
            write (F_unout,9900) ' PROBLEM WITH COMPRESSED MAP SENDING '
            call flush(F_unout)
            goto 998
         endif

         if (F_print_L) then
           write(F_unout,*) 'CPLOCN_DISTWGT: cmap topology gni,gnj,cplocn_atm_nxchng=', &
                                                           gni,gnj,cplocn_atm_nxchng
           write(F_unout,*) 'CPLOCN_DISTWGT: cmap topology cplocn_gni,cplocn_gnj,cplocn_nxchng=', &
                                                           cplocn_gni,cplocn_gnj,cplocn_nxchng
         endif

         deallocate(atm_ijtoc)
      
      endif
 
      call RPN_COMM_bcast (nw,            1, "MPI_INTEGER", 0,"GRID",ier)
      call RPN_COMM_bcast (cplocn_nxchng, 1, "MPI_INTEGER", 0,"GRID",ier)
      call RPN_COMM_bcast (ocn_ijtoc, cplocn_gni*cplocn_gnj, "MPI_INTEGER", 0,"GRID",ier)

      if (cplocn_myproc /= 0) then

         allocate (ocn_wgt_g (1,1,nw) )
         allocate (ocn_iwgt_g(1,1,nw) )
         allocate (ocn_jwgt_g(1,1,nw) )
         allocate (ocn_mwgt_g(1,1   ) )
         allocate (ocn_awgt_g(1,1   ) )

      endif

      allocate (ocn_wgt (ni,nj,nw) )
      allocate (ocn_mwgt(ni,nj   ) )
      allocate (ocn_awgt(ni,nj   ) )
      allocate (ocn_sint(ni,nj   ) )
      allocate (ocn_cost(ni,nj   ) )
      allocate (ocn_cwgt(ni,nj   ) )

      allocate (ocn_iwgt(ni,nj,nw) )
      allocate (ocn_jwgt(ni,nj,nw) )

      !garr Global array sending data  integer, real or real*8 I
      !gmini,gmaxi,gminj,gmaxj Size of garr integer I
      !nig,njg Domain size  integer I
      !nk Z axis size integer I
      !ghalox,ghaloy Halo size of garr for x and y axis integer I
      !size 1 for integer and real, 2 for real*8, etc. integer I
      !larr Local array that will get data integer, real ou real*8 O
      !mini,maxi,minj,maxj Size of larr integer I
      !halox,haloy Halo size of larr for x and y axis integer I
      !periodx,periody Global periodicity over x and y axis logical I
      !ierr ierr (0 if ok, non-0 if error) integer I

! Distribute weight interpolation variables

!drv_glb_ni (gni), drv_glb_nj (gnj), drv_lcl_ni (ni), drv_lcl_nj (ni)

      do k=1,nw

        call RPN_COMM_dist(                           &
                 ocn_wgt_g(:,:,k),                    &
                 1,gni,1,gnj,gni,gnj,1,               &
                 0,0,2,wk8,                           &
                 cpl_minx,cpl_maxx,cpl_miny,cpl_maxy, &
                 0,0,.false.,.false.,errcode)
        if (errcode < 0) goto 998
        ocn_wgt(1:ni,1:nj,k)=wk8(1:ni,1:nj)

        call RPN_COMM_dist(                           &
                 ocn_iwgt_g(:,:,k),                   &
                 1,gni,1,gnj,gni,gnj,1,               &
                 0,0,1,wki,                           &
                 cpl_minx,cpl_maxx,cpl_miny,cpl_maxy, &
                 0,0,.false.,.false.,errcode)
        if (errcode < 0) goto 998
        ocn_iwgt(1:ni,1:nj,k)=wki(1:ni,1:nj)

        call RPN_COMM_dist(                           &
                 ocn_jwgt_g(:,:,k),                   &
                 1,gni,1,gnj,gni,gnj,1,               &
                 0,0,1,wki,                           &
                 cpl_minx,cpl_maxx,cpl_miny,cpl_maxy, &
                 0,0,.false.,.false.,errcode)
        if (errcode < 0) goto 998
        ocn_jwgt(1:ni,1:nj,k)=wki(1:ni,1:nj)

      enddo

      if ( cplocn_myproc == 0 .and. cplocn_ocnf_nsprd > 0 ) then
         do j=1,gnj
         do i=1,gni
            if ( i <= cplocn_ocnf_nsprd .or. &
                 i >= (gni - cplocn_ocnf_nsprd + 1) .or. &
                 j <= cplocn_ocnf_nsprd .or. &
                 j >= (gnj - cplocn_ocnf_nsprd + 1) ) ocn_mwgt_g(i,j)=0
         enddo
         enddo
      endif

      call RPN_COMM_dist(                           &
               ocn_mwgt_g,                          &
               1,gni,1,gnj,gni,gnj,1,               &
               0,0,1,wki,                           &
               cpl_minx,cpl_maxx,cpl_miny,cpl_maxy, &
               0,0,.false.,.false.,errcode)
      if (errcode < 0) goto 998
      ocn_mwgt(1:ni,1:nj)=wki(1:ni,1:nj)

      call RPN_COMM_dist(                           &
               ocn_awgt_g,                          &
               1,gni,1,gnj,gni,gnj,1,               &
               0,0,1,wk4,                           &
               cpl_minx,cpl_maxx,cpl_miny,cpl_maxy, &
               0,0,.false.,.false.,errcode)
      if (errcode < 0) goto 998
      ocn_awgt(1:ni,1:nj)=wk4(1:ni,1:nj)

      deallocate (ocn_wgt_g,  ocn_iwgt_g, ocn_jwgt_g, &
                  ocn_mwgt_g, ocn_awgt_g)


! Adjust rotation angles
! (ocn_awgt originally represents dst grid relative to src grid, anticlockwise)

      do j=1,nj
      do i=1,ni
        ocn_awgt(i,j)=-ocn_awgt(i,j)*pi/180.
        ! (anticlockwise ==> clockwise, see rotation routine)
        ocn_sint(i,j)=sin(ocn_awgt(i,j))
        ocn_cost(i,j)=cos(ocn_awgt(i,j))
      enddo
      enddo
      deallocate(ocn_awgt)

! Compute number of valid weights for each tile

      cplocn_atm_nw_max=0
      do j=1,nj
      do i=1,ni
        ocn_cwgt(i,j) = 0
        if (ocn_mwgt(i,j).eq.1) then
          found = .false.
          do k=nw,1,-1
            if ( .not. found ) then
              if ( abs( ocn_wgt(i,j,k) ) > epsw ) then
                !Small negative weights are possible under limitroff point positions
                found = .true.
                ocn_cwgt(i,j) = k
              endif
            endif
          enddo
          cplocn_atm_nw_max = max(cplocn_atm_nw_max,ocn_cwgt(i,j))
        endif
      enddo
      enddo

! Quality control on compressed wet point indices
      if ( cplocn_atm_nw_max /= 0 ) then
         do j=1,nj
         do i=1,ni
           do k=1,ocn_cwgt(i,j)
             ic=ocn_ijtoc(ocn_iwgt(i,j,k),ocn_jwgt(i,j,k))
             if ( ic < 1 .or. ic > cplocn_nxchng ) then
                  write (F_unout,9900) 'PROBLEM WITH RECEIVED COMPRESSED INDICES ocn_ijtoc'
                  write (F_unout,*) 'ic=',ic
                  write (F_unout,*) 'cplocn_nxchng=',cplocn_nxchng
                  write (F_unout,*) 'i,j,k=',i,j,k
                  write (F_unout,*) 'ocn_iwgt(i,j,k)=',ocn_iwgt(i,j,k)
                  write (F_unout,*) 'ocn_jwgt(i,j,k)=',ocn_jwgt(i,j,k)
                  do k2=1,ocn_cwgt(i,j)
                    write (F_unout,*) 'i,j,k2=',i,j,k2
                    write (F_unout,*) 'ocn_iwgt(i,j,k2)=',ocn_iwgt(i,j,k2)
                    write (F_unout,*) 'ocn_jwgt(i,j,k2)=',ocn_jwgt(i,j,k2)
                  enddo
                  write (F_unout,*) 'TRY TO INCREASE cplocn_ocnf_nsprd IF REGIONAL CONFIGURATION'
                  write (F_unout,*) 'IN GEM_SETTINGS...'
                  call flush(F_unout)
                  errcode=-1
                  goto 998
             endif 
           enddo
         enddo
         enddo
      endif

      goto 999

 998  call handle_error(errcode,'cplocn_distwgt','Problems in weight distribution')
      if (F_print_L) write (F_unout,*) 'cplocn_distwgt: Finished'

 999  continue

9900 format (/,1x,a)

!     ________________________________________________________________
!
      return
      end
