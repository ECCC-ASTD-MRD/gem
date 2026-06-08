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

!**s/r cplocn_xchng

      subroutine cplocn_xchng(F_date_ou, F_date_in,  &
                              F_send, F_recv, F_stepdriver, F_err)
use rmn_gmm
use iso_c_binding
      implicit none
#include <arch_specific.hf>
!
!authors    Francois Roy - Fall 2014
! 
!revision
! v4_7  - Roy F.  - initial version

! Purpose: Performs exchange with ocean model
!          through gossip interface and apply 
!          distributed interpolation weights 
!          (including angles and vector rotation)

#include <rmn/msg.h>

      include "cpl.cdk"
      include "cplocn.cdk"

      character*16 F_date_ou, F_date_in
      integer F_send, F_recv, F_stepdriver, F_err

      real, pointer, dimension(:,:,:), save :: ocn_glb_busou
      real, pointer, dimension(:),     save :: ocn_cmp_busou
      real, pointer, dimension(:),     save :: ocn_cmp_busin

      real,   dimension (cpl_minx:cpl_maxx,cpl_miny:cpl_maxy) :: wk4

      integer unout
      logical print_L, debug_w

      integer, external :: msg_getUnit

      integer ivar,ni,nj,gni,gnj,gni_o,gnj_o,i,j,k
      integer iproc,ic,ier,status

!     ________________________________________________________________
!
      ier = gmm_get(gmmk_ocn_busin_s , ocn_busin)

      ni=cpl_drv_lni
      nj=cpl_drv_lnj

      gni=cpl_drv_gni
      gnj=cpl_drv_gnj

      gni_o=cplocn_gni   !Ocean model global grid
      gnj_o=cplocn_gnj

      unout   = msg_getUnit(MSG_INFO)
      print_L = (unout > 0)
      debug_w = cplocn_debug_L

      if ( .not. associated(ocn_glb_busou) ) then
        if ( cplocn_myproc == 0 ) then
          allocate(ocn_glb_busou(gni,gnj,1))
          allocate(ocn_cmp_busou(cplocn_atm_nxchng))
         else
          allocate(ocn_glb_busou(1,1,1))
        endif
      endif

! Send step
      do ivar = 1, cplocn_n_fldou

         !garr Global array receiving data  integer, real or real*8 O
         !gmini,gmaxi,gminj,gmaxj Size of garr integer I
         !nig,njg Domain size  integer I
         !nk Z axis size integer I
         !ghalox,ghaloy Halo size of garr for x and y axis integer I
         !size 1 for integer and real, 2 for real*8, etc. integer I
         !larr Local array that will be collected integer, real ou real*8 I
         !mini,maxi,minj,maxj Size of larr integer I
         !halox,haloy Halo size of larr for x and y axis integer I
         !ierr ierr (0 if ok, non-0 if error) integer I

         wk4(1:ni,1:nj)=ocn_busou(1:ni,1:nj,ivar)

         call RPN_COMM_coll (                             &
                  ocn_glb_busou,                          &
                  1,gni,1,gnj,gni,gnj,1,                  &
                  0,0,1,wk4,                              &
                  cpl_minx,cpl_maxx,cpl_miny,cpl_maxy,0,0,F_err)
         if (F_err /= 0)  then
           if (print_L) &
              write(unout,*) 'cplocn_xchng: RPN_COMM_coll error'
           F_err = -1
           goto 998
         endif
         if (cplocn_myproc.eq.0.and.debug_w) then
           call cplao_write (F_date_ou, ocn_glb_busou(:,:,1), &
                             cplocn_cvou_S(ivar), gni, gnj, F_stepdriver, 'put' )
         endif

         if (cplocn_myproc.eq.0) then
            do i=1,cplocn_atm_nxchng
              ocn_cmp_busou(i)=ocn_glb_busou(atm_ctoi(i),atm_ctoj(i),1)
            enddo
            if ( .not. debug_w ) then
              call cplao_xchng (F_date_ou, ocn_cmp_busou, cplocn_atm_nxchng, &
                                 'put', cplocn_cvou_S(ivar),                 &
                                 F_send, F_stepdriver, F_err)
            else
              call cplao_xchng_debug (F_date_ou, ocn_cmp_busou, cplocn_atm_nxchng, &
                                      'put', cplocn_cvou_S(ivar),                 &
                                      ocn_ijtoc, gni_o, gnj_o, F_send, debug_w,   &
                                      F_stepdriver, F_err)
              ! Here ocn_ijtoc, gni_o, gnj_o are not used
              ! (put mode)
            endif
         endif

         call RPN_COMM_bcast (F_err, 1, "MPI_INTEGER", 0,"GRID",ier)
         if (F_err /= 0)  then
           F_err = -1
           if (print_L) &
              write(unout,*) 'cplocn_xchng: cpl_ao_xchng put error'
           goto 998
         endif
 
      enddo

      if (.not.associated(ocn_cmp_busin)) &
         allocate(ocn_cmp_busin(cplocn_nxchng))

! Recieve step
      do ivar = 1, cplocn_n_fldin

         if ( cplocn_myproc == 0 ) then
           if ( .not. debug_w ) then
             call cplao_xchng (F_date_in,  ocn_cmp_busin, cplocn_nxchng, &
                               'get', cplocn_cvin_S(ivar),               &
                               F_recv, F_stepdriver, F_err)
           else
             call cplao_xchng_debug (F_date_in,  ocn_cmp_busin, cplocn_nxchng, &
                                     'get', cplocn_cvin_S(ivar),               &
                                     ocn_ijtoc, gni_o, gnj_o, F_recv, debug_w, &
                                     F_stepdriver, F_err)
             ! Here ocn_ijtoc, gni_o, gnj_o are used 
             ! only for debug mode, writing get mode fields
           endif
         endif

         call RPN_COMM_bcast (F_err, 1, "MPI_INTEGER", 0,"GRID",ier)
         if (F_err /= 0)  then
           F_err = -1
           if (print_L) &
              write(unout,*) 'cplocn_xchng: cpl_ao_xchng get error'
           goto 998
         endif
 
         call RPN_COMM_bcast (F_recv, 1, "MPI_INTEGER", 0,"GRID",ier)

         if (F_recv.eq.0) then
           if (ivar == 1 .and. print_L) &
              write(unout,*) 'cplocn_xchng: BUFFER UPDATE, F_date_in(1:15)=', F_date_in(1:15)

           ! distributes interpolation values

           call RPN_COMM_bcast (ocn_cmp_busin, cplocn_nxchng, "MPI_REAL", 0,"GRID",ier)

           ! distributed interpolation
           ocn_busin(:,:,ivar)=0.
           if ( cplocn_atm_nw_max /= 0 ) then
              do j=1,nj
              do i=1,ni
                do k=1,ocn_cwgt(i,j)
                    ic=ocn_ijtoc(ocn_iwgt(i,j,k),ocn_jwgt(i,j,k))
                    ocn_busin(i,j,ivar)=ocn_busin(i,j,ivar)+ &
     &                    ocn_wgt(i,j,k)*ocn_cmp_busin(ic)
                enddo
              enddo
              enddo

           endif

          else

           if (ivar == cplocn_n_fldin .and. print_L) &
              write(unout,*) 'cplocn_xchng: BUFFER NOT TOUCHED, F_date_in(1:15)=', F_date_in(1:15)

         endif

      enddo
  
      if (F_recv.eq.0) then 
         ! vector rotation
         if ( cplocn_atm_nw_max /= 0 ) then
           do ivar=1, cplocn_n_fldin
              if (cplocn_cvit_S(ivar) == 'U') then
                 if (ivar+1.gt.cplocn_n_fldin) then
                    if ( print_L ) write(unout,*) 'cplocn_xchng: Problems with vector ordering 1'
                    F_err = -1
                    goto 998
                 endif
                 if (cplocn_cvit_S(ivar+1) /= 'V') then
                    if ( print_L ) write(unout,*) 'cplocn_xchng: Problems with vector ordering 2'
                    F_err = -1
                    goto 998
                 endif
!                 call rotate_vector_ninj(ocn_busin(:,:,ivar),   &
!                                         ocn_busin(:,:,ivar+1), &
!                                         ocn_awgt,ni,nj)
                 call rotate_vector_ninj2(ocn_busin(:,:,ivar),   &
                                          ocn_busin(:,:,ivar+1), &
                                          ocn_sint,ocn_cost,ni,nj)
              endif
           enddo
         endif
      endif

      F_err = 0

      goto 999

 998  call handle_error(F_err,'cplocn_xchng','Problems')
 999  continue

!     ________________________________________________________________
!
      return
      end

      SUBROUTINE rotate_vector_ninj(u,v,theta,ni,nj)
      implicit none
      integer ni,nj,i,j
      real u(ni,nj),v(ni,nj),theta(ni,nj)
      real utmp(ni,nj),vtmp(ni,nj),ct,st

      do j=1,nj
      do i=1,ni
        st=sin(theta(i,j))
        ct=cos(theta(i,j))
        utmp(i,j) = ct*u(i,j) - st*v(i,j)
        vtmp(i,j) = st*u(i,j) + ct*v(i,j)
      enddo
      enddo

      u(:,:) = utmp(:,:)
      v(:,:) = vtmp(:,:)

      return
      END SUBROUTINE rotate_vector_ninj

      SUBROUTINE rotate_vector_ninj2(u,v,st,ct,ni,nj)
      implicit none
      integer ni,nj,i,j
      real u(ni,nj),v(ni,nj),st(ni,nj),ct(ni,nj)
      real utmp(ni,nj),vtmp(ni,nj)

      do j=1,nj
      do i=1,ni
        utmp(i,j) = ct(i,j)*u(i,j) - st(i,j)*v(i,j)
        vtmp(i,j) = st(i,j)*u(i,j) + ct(i,j)*v(i,j)
      enddo
      enddo

      u(:,:) = utmp(:,:)
      v(:,:) = vtmp(:,:)

      return
      END SUBROUTINE rotate_vector_ninj2
