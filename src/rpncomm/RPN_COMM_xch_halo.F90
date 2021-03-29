!/* RPN_COMM - Library of useful routines for C and FORTRAN programming
! * Copyright (C) 1975-2015  Division de Recherche en Prevision Numerique
! *                          Environnement Canada
! *
! * This library is free software; you can redistribute it and/or
! * modify it under the terms of the GNU Lesser General Public
! * License as published by the Free Software Foundation,
! * version 2.1 of the License.
! *
! * This library is distributed in the hope that it will be useful,
! * but WITHOUT ANY WARRANTY; without even the implied warranty of
! * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! * Lesser General Public License for more details.
! *
! * You should have received a copy of the GNU Lesser General Public
! * License along with this library; if not, write to the
! * Free Software Foundation, Inc., 59 Temple Place - Suite 330,
! * Boston, MA 02111-1307, USA.
! */
!=======================================================================
      integer function RPN_COMM_xch_halo_flip_test(nparams,params)
!=======================================================================
      use rpn_comm
      implicit none
      integer, intent(IN) :: nparams
      integer, intent(IN), dimension(nparams) :: params
      integer :: lni
      integer :: lnj
      integer :: nk
      integer :: gni, gnj
      logical :: periodx, periody
      integer, dimension(pe_nx) :: countx, offsetx
      integer, dimension(pe_ny) :: county, offsety
      integer :: halox, haloy, npol_row
      integer :: lminx, lmaxx, lminy, lmaxy
      integer :: minx, maxx, miny, maxy
      integer :: minx1, maxx1, miny1, maxy1
      integer, pointer, dimension(:,:,:) :: localarray,localarray2
      integer :: i, j, k, ierr
      integer, external :: RPN_COMM_limit
!
      RPN_COMM_xch_halo_flip_test=-1
      gni = params(1)
      gnj = params(2)
      nk = params(3)
      halox=params(4)
      haloy=params(5)
      periodx = .false.
      periody = .false.
      periodx = params(6).ne.0
      periody = params(7).ne.0

      npol_row = -999999
      ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      lni = countx(pe_mex+1)
      ierr = RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      lnj = county(pe_mey+1)
      minx = lminx-halox ; minx1 = minx - lminx + 1
      maxx = lmaxx+halox ; maxx1 = maxx - lminx + 1
      miny = lminy-haloy ; miny1 = miny - lminy + 1
      maxy = lmaxy+haloy ; maxy1 = maxy - lminy + 1

      allocate(localarray(minx:maxx,miny:maxy,nk))
      allocate(localarray2(minx:maxx,miny:maxy,nk))
      localarray = 1111
      localarray = 2222
      do k = 1,nk
      do j = lminy,lmaxy
      do i = lminx,lmaxx
        localarray(i,j,k) = (i - 1)*(360.0/gni)
        localarray2(i,j,k) = -90 + (j - .5) * (180.0/(gnj))
      enddo
      enddo
      enddo
      do j = maxy,miny,-1
        print 101,j,localarray(minx:maxx,j,1),-999,localarray2(minx:maxx,j,1)
      enddo
101   format(40I5)
      call rpn_comm_xch_halo(localarray,minx1,maxx1,miny1,maxy1,lni,lnj,nk,halox,haloy,periodx,periody,gni,npol_row)
      call rpn_comm_xch_halo(localarray2,minx1,maxx1,miny1,maxy1,lni,lnj,nk,halox,haloy,periodx,periody,gni,npol_row)
      do j = maxy,miny,-1
        print 101,j,localarray(minx:maxx,j,1),-999,localarray2(minx:maxx,j,1)
      enddo

      RPN_COMM_xch_halo_flip_test=0
      return
      end function RPN_COMM_xch_halo_flip_test
!
!=======================================================================
      function RPN_COMM_exchange_halo_test(nparams,params) result(status)
      use ISO_C_BINDING
      implicit none
      include 'RPN_COMM.inc'
      integer, intent(IN) :: nparams
      integer, intent(IN), dimension(nparams) :: params
      integer :: status
      status = 0
      end function RPN_COMM_exchange_halo_test
!=======================================================================
!=======================================================================
      integer function RPN_COMM_xch_halo_test(nparams,params)
!=======================================================================
      use rpn_comm
      implicit none
      include 'RPN_COMM_interfaces.inc'
      integer, intent(IN) :: nparams
      integer, intent(IN), dimension(nparams) :: params
!
      integer, pointer, dimension(:,:,:) :: localarray
      integer*8, pointer, dimension(:,:,:) :: localarray2
      integer :: lni
      integer :: lnj
      integer :: nk
      integer :: gni, gnj
      integer :: i, j, k, ierr
      integer :: lminx, lmaxx, lminy, lmaxy
      integer :: minx, maxx, miny, maxy
      integer :: minx1, maxx1, miny1, maxy1
      integer, dimension(pe_nx) :: countx, offsetx
      integer, dimension(pe_ny) :: county, offsety
      integer :: halox, haloy, npol_row, errors, value, ii, jj
      logical :: periodx, periody
!
      RPN_COMM_xch_halo_test=-1
      gni = params(1)
      gnj = params(2)
      nk = params(3)
      halox=params(4)
      haloy=params(5)
      periodx = .false.
      periody = .false.
      periodx = params(6).ne.0
      periody = params(7).ne.0
!
      ierr = RPN_COMM_limit(pe_mex,pe_nx,1,gni,lminx,lmaxx,countx,offsetx)
      lni = countx(pe_mex+1)
      ierr = RPN_COMM_limit(pe_mey,pe_ny,1,gnj,lminy,lmaxy,county,offsety)
      lnj = county(pe_mey+1)
      minx = lminx-halox ; minx1 = minx - lminx + 1
      maxx = lmaxx+halox ; maxx1 = maxx - lminx + 1
      miny = lminy-haloy ; miny1 = miny - lminy + 1
      maxy = lmaxy+haloy ; maxy1 = maxy - lminy + 1
      if(pe_me==pe_nx*pe_ny-1) write(rpn_u,100)  &
          'grid halo exchange test',  &
          pe_tot_grid,pe_nx,pe_ny,lminx,lmaxx,lminy,lmaxy,countx,county,  &
          minx1,maxx1,miny1,maxy1
100   format(A,25I5)
      allocate(localarray(minx:maxx,miny:maxy,nk))
      allocate(localarray2(minx:maxx,miny:maxy,nk))
!
      localarray = 99999
      localarray2 = 99999
      do k = 1,nk
      do j = lminy,lmaxy
      do i = lminx,lmaxx
        localarray(i,j,k) = k + 10*j + 1000*i
        localarray2(i,j,k) = k + 10*j + 1000*i
      enddo
      enddo
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
!
      npol_row = 0
!      return
      call RPN_COMM_xch_halo(localarray,minx1,maxx1,miny1,maxy1, &
                   lni,lnj,nk,halox,haloy,periodx,periody,  &
                  gni,npol_row)

      call mpi_barrier(MPI_COMM_WORLD,ierr)
      if(pe_mex==0 .and. pe_mey==0)then
        do j=lmaxy+haloy,lminy-haloy,-1
          write(rpn_u,90)j,localarray(:,j,1)
        enddo
90      format(X,I5.4,20(X,I5.5))
      endif
      errors=0
      do k=1,nk
      do j=lminy-haloy,lmaxy+haloy
!      do j=lminy,lmaxy+1
      do i=lminx-halox,lmaxx+halox
!      do i=lminx,lmaxx
        ii = i
        if(i>gni .and. .not. periodx) cycle
        if(i<1   .and. .not. periodx) cycle
        if(i>gni) ii=i-gni
        if(i<1)   ii=i+gni
        jj = j
        if(jj>gnj .and. .not. periody) cycle
        if(jj<1   .and. .not. periody) cycle
        if(jj>gnj) jj=j-gnj
        if(jj<1)   jj=j+gnj
        value = k+10*jj+1000*ii
        if(localarray(i,j,k) /= value)errors=errors+1
      enddo
      enddo
      enddo
      call mpi_barrier(MPI_COMM_WORLD,ierr)
      write(rpn_u,100)'errors=',errors,pe_mex,pe_mey
!
      RPN_COMM_xch_halo_test=0
      return
      end function RPN_COMM_xch_halo_test
!=======================================================================
!InTf!
      SUBROUTINE RPN_COMM_exchange_halo(pattern,array,periodx,periody,periodp,npol_row)  !InTf!
      use rpn_comm
!!    import ::  rpncomm_pattern, c_ptr                                                  !InTf!
      implicit none                                                                      !InTf!
      type(rpncomm_pattern), intent(IN) :: pattern                                       !InTf!
      type(c_ptr), intent(IN) :: array                                                   !InTf!
      logical, intent(IN), OPTIONAL  :: periodx,periody, periodp                         !InTf!
      integer, intent(IN), OPTIONAL :: npol_row                                          !InTf!
      integer, dimension(:), pointer :: g
      integer :: minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,gni
      type(rpncomm_field), pointer :: f

      call c_f_pointer(pattern%p,f)
      minx = f%x%lo
      maxx = f%x%hi
      miny = f%y%lo
      maxy = f%y%hi
      ni = f%x%lnp
      nj = f%y%lnp
      nk = f%z%lnp
      halox = f%hx
      haloy = f%hy
      gni = f%x%gnp
      call c_f_pointer(array,g,[(maxx-minx+1)*(maxy-miny+1)*nk])

      if(npol_row > 0) then  !  semi lag exchange, call old code
        call RPN_COMM_xch_halo(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row)
      else
        call RPN_COMM_exchange_ew(g,minx,maxx,miny,maxy,ni,nj,nk,halox,periodx,f%ew3d)
        call RPN_COMM_exchange_ns(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periody,f%ns3d)
        if(periodp) call RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,gni)
      endif

      end SUBROUTINE RPN_COMM_exchange_halo                                              !InTf!

      SUBROUTINE RPN_COMM_exchange_ew(g,minx,maxx,miny,maxy,ni,nj,nk,halox,periodx,ew3d)
      use rpn_comm
      implicit none
      integer, intent(IN) :: minx,maxx,miny,maxy,ni,nj,nk,halox,ew3d
      logical, intent(IN) :: periodx
      integer, dimension(minx:maxx,miny:maxy,nk), intent(INOUT) :: g
      logical :: east, west
      integer :: eastpe, westpe
      integer :: ierr
      integer, dimension(MPI_STATUS_SIZE) :: status

      east=(bnd_east) .and. (.not.periodx)
      eastpe=pe_id(pe_mex+1,pe_mey)
      west=(bnd_west) .and. (.not.periodx)
      westpe=pe_id(pe_mex-1,pe_mey)

      !   message tag is 1000 + pe_id of sender
      !   halo to east starts at g(ni-halox+1,1,1)
      !   halo from west starts at g(1-halox,1,1)
      !   halo to west starts at g(1,1,1)
      !   halo from east starts at g(ni+1,1,1)
      if(west) then                !   send to east_neighbor
        if(.not.east)then
          call MPI_SEND(g(ni-halox+1,1,1),1,ew3d,eastpe,1000+pe_medomm,PE_DEFCOMM,ierr)
        endif
      else if(east) then           !   receive from west_neighbor
        call MPI_RECV(g(1-halox,1,1),1,ew3d,westpe,1000+westpe,PE_DEFCOMM,status,ierr)
      else                         !   send to east_neighbor and receive from west_neighbor
        call MPI_SENDRECV(                                    &
               g(ni-halox+1,1,1),1,ew3d,eastpe,1000+pe_medomm,  &
               g(1-halox,1,1),1,ew3d,westpe,1000+westpe, &
               PE_DEFCOMM,status,ierr)
      endif
      if(east) then                !   send to west_neighbor
        if(.not.west) then
          call MPI_SEND(g(1,1,1),1,ew3d,westpe,1000+pe_medomm,PE_DEFCOMM,ierr)
        endif
      else if(west) then           !   receive from east_neighbor
        call MPI_RECV(g(ni+1,1,1),1,ew3d,eastpe,1000+eastpe,PE_DEFCOMM,status,ierr)
      else                         !   send to west_neighbor and receive from east_neighbor
        call MPI_SENDRECV(                                    &
               g(1,1,1),1,ew3d,westpe,1000+pe_medomm,  &
               g(ni+1,1,1),1,ew3d,eastpe,1000+eastpe, &
               PE_DEFCOMM,status,ierr)
      endif
      end SUBROUTINE RPN_COMM_exchange_ew

      SUBROUTINE RPN_COMM_exchange_ns(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periody,ns3d)
      use rpn_comm
      implicit none
      integer, intent(IN) :: minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,ns3d
      logical, intent(IN) :: periody
      integer, dimension(minx:maxx,miny:maxy,nk), intent(INOUT) :: g
      logical :: north, south
      integer :: northpe, southpe

      north=(bnd_north) .and. (.not.periody)
      northpe=pe_id(pe_mex,pe_mey+1)
      south=(bnd_south) .and. (.not.periody)
      southpe=pe_id(pe_mex,pe_mey-1)
      end SUBROUTINE RPN_COMM_exchange_ns
!=======================================================================
      SUBROUTINE RPN_COMM_xch_halo(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,gni,npol_row)
!=======================================================================
!
!  +=======================================================================================+___maxy
!  I       <halox>                                                           <halox>       I
!  I +-----+-----------------------------------------------------------------------+-----+ I
!  I :     :                                    |                                  :     : I
!  I :"cTL":               "iTLCR"            haloy                                :"cTR": I
!  I :     :                                    |                                  :     : I
!  I +-----+=====+===========================================================+=====+-----: I___nj
!  I :     I     :                              |                            :     I     : I
!  I :"iTL"I" TL":            "TC"            haloy                          :"TR "I"iTR": I
!  I :     I     :                              |                            :     I     : I
!  I :-----+-----+-----------------------------------------------------------+-----+-----: I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :"iCL"I" CL":                                                           :"CR "I"iCR": I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :     I     :                                                           :     I     : I
!  I :-----+-----+-----------------------------------------------------------+-----+-----: I
!  I :     I     :                              |                            :     I     : I
!  I :"iBL"I" BL":            "BC"            haloy                          :"BR "I"iBR": I
!  I :     I     :                              |                            :     I     : I
!  I +-----+=====+===========================================================+=====+-----: I___1
!  I :     :                                    |                                  :     : I
!  I :"cBL":               "iBLCR"            haloy                                :"cBR": I
!  I :     :                                    |                                  :     : I
!  I +-----+-----------------------------------------------------------------------+-----+ I
!  I <halox>                                                                       <halox> I
!  +=======================================================================================+___miny
!  |       |                                                                       |       |
!  |       |                                                                       |       |
!  |       1                                                                       ni      |
!  minx                                                                                 maxx
!
!    Description of the 8 internal arrays 
!
!  BR_CR_TR array     : 3 pieces  BR[nk*halox*haloy],  CR[nk*halox*(nj-2*haloy)],  TR[nk*halox*haloy]
!                       used for West to East send  (nk*halox*nj)
!  iBR_iCR_iTR array  : 3 pieces iBR[nk*halox*haloy], iCR[nk*halox*(nj-2*haloy)], iTR[nk*halox*haloy]
!                       used for East to West receive  (nk*halox*nj)
!                       (iBR, iTR then used for North/South send)
!  BL_CL_TL array     : 3 pieces  BL[nk*halox*haloy],  CL[nk*halox*(nj-2*haloy)],  TL[nk*halox*haloy]
!                       used for East to West send  (nk*halox*nj)
!  iBL_iCL_iTL array  : 3 pieces iBL[nk*halox*haloy], iCL[nk*halox*(nj-2*haloy)], iTL[nk*halox*haloy]
!                       used for West to East receive  (nk*halox*nj)
!                       (iBL, iTL then used for North/South send)
!  iBL_BLCR_iBR array : 3 pieces iBL[nk*halox*haloy],  BLCR[nk*haloy*ni],         iBR[nk*halox*haloy]
!                       used for North to South send (nk*haloy*(ni+2*halox))
!  iTL_TLCR_iTR array : 3 pieces iTL[nk*halox*haloy],  TLCR[nk*haloy*ni],         iTR[nk*halox*haloy]
!                       used for South to North send (nk*haloy*(ni+2*halox))
!  cBL_iBLCR_cBR array: 3 pieces cBL[nk*halox*haloy],  iBLCR[nk*haloy*ni],        cBR[nk*halox*haloy]
!                       used for South to North receive (nk*haloy*(ni+2*halox))
!  cTL_iTLCR_cTR array: 3 pieces cTL[nk*halox*haloy],  iTLCR[nk*haloy*ni],        cTR[nk*halox*haloy]
!                       used for North to South receive (nk*haloy*(ni+2*halox))
!
!    Communication tags (for iSEND/Irecv)
!
!  West  -> East  : sender_pe + 100000
!  West  <- East  : sender_pe
!  North -> South : sender_pe + 100000
!  North <- South : sender_pe
!
      use rpn_comm
      implicit none
!
!     exchange a halo with N/S/E/W neighbours
!
      integer minx,maxx,miny,maxy,ni,nj,nk,halox,haloy
      integer gni,npol_row
      logical periodx,periody
!     integer *8 mem_time, exch_time, ewtime
      integer g(minx:maxx,miny:maxy,nk)
!
        integer, dimension(halox*nj*nk) :: BR_CR_TR, BL_CL_TL  ! send to East/West
        integer, dimension(halox*nj*nk) :: iBR_iCR_iTR, iBL_iCL_iTL   ! recv from East/West
         ! send to South/North
         ! recv from South/North
        integer, dimension(haloy*(ni+2*halox)*nk) ::   &
                 iBL_BLCR_iBR, iTL_TLCR_iTR,  & 
                 cBL_iBLCR_cBR, cTL_iTLCR_cTR  
        integer :: bl,bc,br,cl,cr,tl,tc,tr
        integer :: cbl,iblcr,cbr,ibl,ibr,icl,icr
        integer :: itl,itr,ctl,itlcr,ctr
        integer :: to_east, to_west, to_north, to_south
        integer :: from_east, from_west, from_north, from_south
!     integer *8 time_base,temp_time
      integer i, j, k, m, m1, m2, m3
      integer nwds_ew, nwds_ns
      integer sendtag, gettag, ierr
      integer status(MPI_STATUS_SIZE)
      logical east,west,north,south
      integer eastpe,westpe,northpe,southpe
      integer :: east_m, east_m_n
      integer :: west_m, west_m_n
!
      integer globalni,polarrows, nilmax, jmin,jmax
        integer RPN_COMM_topo, mini,maxi,nil,ni0
!
      integer land_fill
      real r_land_fill
      equivalence(land_fill,r_land_fill)
      integer :: halo_status
!       integer, external :: rpn_comm_set_valid_halo
!
      globalni=abs(gni)
      polarrows=npol_row
      nwds_ew = size(BR_CR_TR)
      nwds_ns = size(iBL_BLCR_iBR)

1     continue
!       if(polarrows==0) then
!         halo_status = rpn_comm_set_valid_halo(g,halox,haloy,0)
!       else
!         halo_status = rpn_comm_set_valid_halo(g,-1,-1,-1)   ! if there are polar rows, invalidate halo table entry
!       endif
      
!     call RPN_COMM_tmg_in
      east=(bnd_east) .and. (.not.periodx)
      eastpe=pe_id(pe_mex+1,pe_mey)
      east_m = 0
      if(east) east_m=not(east_m)
      east_m_n = not(east_m)
!
      west=(bnd_west) .and. (.not.periodx)
      westpe=pe_id(pe_mex-1,pe_mey)
      west_m = 0
      if(west) west_m=not(west_m)
      west_m_n = not(west_m)
!
      north=(bnd_north) .and. (.not.periody)
      northpe=pe_id(pe_mex,pe_mey+1)
      south=(bnd_south) .and. (.not.periody)
      southpe=pe_id(pe_mex,pe_mey-1)

      jmin = 1
      jmax = nj 
      if(rpn_ew_ext_L) then            !  add haloy extra rows for EW halo exchange 
         if(north) jmax = nj+haloy     !  above for North tiles
         if(south) jmin = 1-haloy      !  below for South tiles
      endif
!
      if(pe_opcv(1) .ne. ' ') then !  fill halo option present
           r_land_fill=pe_oprv(1)
           if(pe_opcv(1) .eq. 'BAND') then
              if(iand(1,pe_opiv(1)) .ne. 0) then ! south band
                 do j=miny,0
                    do k=1,nk
                       do i=minx,maxx
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
              endif
            if(iand(2,pe_opiv(1)) .ne. 0) then ! east band
                 do i=ni+1,maxx
                    do j=miny,maxy
                       do k=1,nk
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
            endif
            if(iand(4,pe_opiv(1)) .ne. 0) then ! north band
                 do j=nj+1,maxy
                    do k=1,nk
                       do i=minx,maxx
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
            endif
            if(iand(8,pe_opiv(1)) .ne. 0) then ! west band
                 do i=minx,0
                    do j=miny,maxy
                       do k=1,nk
                          g(i,j,k)=land_fill
                       enddo
                    enddo
                 enddo
            endif
           endif
           if(pe_opcv(1) .eq. 'EDGE') then
              if(iand(1,pe_opiv(1)) .ne. 0) then ! south edge
                 do i=minx,0
                    do k=1,nk
                       g(i,1,k)=land_fill
                    enddo
                 enddo
                 do i=ni+1,maxx
                    do k=1,nk
                       g(i,1,k)=land_fill
                    enddo
                 enddo
            endif
            if(iand(2,pe_opiv(1)) .ne. 0) then ! east edge
                 do j=miny,0
                    do k=1,nk
                       g(ni,j,k)=land_fill
                    enddo
                 enddo
                 do j=nj+1,maxy
                    do k=1,nk
                       g(ni,j,k)=land_fill
                    enddo
                 enddo
            endif
            if(iand(4,pe_opiv(1)) .ne. 0) then ! north edge
                 do i=minx,0
                    do k=1,nk
                       g(i,nj,k)=land_fill
                    enddo
                 enddo
                 do i=ni+1,maxx
                    do k=1,nk
                       g(i,nj,k)=land_fill
                    enddo
                 enddo
            endif
            if(iand(8,pe_opiv(1)) .ne. 0) then ! west edge
                 do j=miny,0
                    do k=1,nk
                       g(1,j,k)=land_fill
                    enddo
                 enddo
                 do j=nj+1,maxy
                    do k=1,nk
                       g(1,j,k)=land_fill
                    enddo
                 enddo
            endif
           endif
        endif
!
!       use new fullly asynchronous code only if pe_nx and pe_ny both >1
!       and polarrows = 0
!
      if(pe_nx>1 .and. pe_ny>1              &
                 .and. polarrows<=0         &
                 .and. full_async_exch      &
                 .and. (.not. rpn_ew_ext_L) ) goto 2
!
!       if no halo along x, bypass
!     call tmg_start(90,'RPN_COMM_haloew')
      if (halox .gt. 0) then
           if (.not.(min(pe_mey+1,pe_ny-pe_mey).le.polarrows)) then
              call RPN_COMM_xch_haloew(g,minx,maxx,miny,maxy,ni,jmin,jmax,nk,halox,haloy,periodx,periody)
              endif
      endif
!     call tmg_stop(90)
!     call tmg_start(91,'RPN_COMM_halons')
      if (haloy .gt. 0) then
              call RPN_COMM_xch_halons(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody)
              if(periody .and. npol_row<-pe_ny .and. globalni>=ni) call RPN_COMM_haloflip(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,globalni)
      endif
!     call tmg_stop(91)
      if (min(pe_mey+1,pe_ny-pe_mey).le.polarrows) then
           ierr = RPN_COMM_topo(globalni,mini,maxi,nil,nilmax,halox,ni0,.TRUE.,.FALSE.)
!          call tmg_start(92,'RPN_COMM_xch_halosl')
           call RPN_COMM_xch_halosl(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody,globalni,npol_row,nilmax)
!          call tmg_stop(92)

      endif
!     call RPN_COMM_tmg_out
      return
!
!       version no 2 of this routine, fully asynchronous, overlapped communication /data movement
!       OOPS , rpn_ew_ext_L option broken, will have to fix it
!
!       step 1 post non blocking receives in the East/West direction
!
2       continue
        if(.not. west) then
          ! get from west neighbor unless i am west PE
          call MPI_IRECV(iBL_iCL_iTL,nwds_ew,MPI_INTEGER,westpe, &
               100000+westpe,PE_DEFCOMM,from_west,ierr)          ! sender was westpe, tag is westpe+100000
        endif
        if(.not. east) then
          ! get from east neighbor unless i am east PE
          call MPI_IRECV(iBR_iCR_iTR,nwds_ew,MPI_INTEGER,eastpe, &
               eastpe,PE_DEFCOMM,from_east,ierr)                 ! sender was eastpe, tag is eastpe
        endif
!
!       step 2 fill East/West send buffers from inner halo
!
        bl = 1                       ; br = bl ; ibl = bl ; ibr = br
        cl = 1 + nk*halox*haloy      ; cr = cl ; icl = cl ; icr = cr
        tl = 1 + nk*halox*(nj-haloy) ; tr = tl ; itl = tl ; itr = tr
        m = 0
        m1 = bl
        m2 = cl
        m3 = tl
        do k=1,nk 
         do j=1,haloy
          do i=1,halox
           BL_CL_TL(m1)=g(i,j,k)             ! BL/BR part
           BR_CR_TR(m1)=g(ni-halox+i,j,k)
           BL_CL_TL(m3)=g(i,nj-haloy+j,k)    ! TL/TR part
           BR_CR_TR(m3)=g(ni-halox+i,nj-haloy+j,k)
           m1=m1+1
           m3=m3+1
          enddo
         enddo
         m = m + 2*halox*haloy
!!        enddo
!!        do k=1,nk
         do j=1+haloy,nj-haloy
          do i=1,halox
           BL_CL_TL(m2)=g(i,j,k)             ! CL/CR part
           BR_CR_TR(m2)=g(ni-halox+i,j,k)
           m=m+1
           m2=m2+1
          enddo
         enddo
!!        enddo
!!        do k=1,nk  ! TL/TR part
!         do j=nj-haloy+1,nj
!         do i=1,halox
!           BL_CL_TL(m3)=g(i,j,k)
!           BR_CR_TR(m3)=g(ni-halox+i,j,k)
!           m=m+1
!           m3=m3+1
!         enddo
!         enddo
        enddo  ! k=1,nk
        if(m /= halox*nj*nk) print *,'OUCH EW send'
!      write(rpn_u,100),'DBG:',ni,nj,nk,halox,haloy,size(BL_CL_TL),
!     %      size(BR_CR_TR),m-1,halox*nj*nk
100   format(A,30I5)
!
!       step 3 non blocking send to East/West partners
!
        if(.not. east) then
          ! send to east neighbor unless i am east PE
          call MPI_ISEND(BR_CR_TR,nwds_ew,MPI_INTEGER,eastpe,   &
               100000+pe_medomm,PE_DEFCOMM,to_east,ierr)         ! tag is PE grid ordinal of sender+100000
        endif
        if(.not. west) then
          ! send to west neighbor unless i am west PE
          call MPI_ISEND(BL_CL_TL,nwds_ew,MPI_INTEGER,westpe,  &
               pe_medomm,PE_DEFCOMM,to_west,ierr)      ! tag is PE grid ordinal of sender
        endif
!
!       step 4 start filling the North/South send buffers
!
        cbl = 1                       ; ctl = 1
        iblcr = 1 + nk*halox*haloy    ; itlcr = iblcr
        cbr = 1 + nk*haloy*(ni+halox) ; ctr = cbr
        m = iblcr
        do k=1,nk  ! TL/TR part
         do j=1,haloy
          do i=1,ni
           iBL_BLCR_iBR(m)=g(i,j,k)
           iTL_TLCR_iTR(m)=g(i,nj-haloy+j,k)
           m=m+1
          enddo
         enddo
        enddo
        if(m /= iblcr+haloy*ni*nk) print *,'OUCH NS part buf'
!
!       step 5 get East/West inbound data, then finish filling the North/South send buffers
!
        if(.not. west) then
          call MPI_wait(from_west,status,ierr)  ! wait for inbound EAST <- West message to complete
          do i = 0,halox*haloy*nk-1             ! put iBL/iTL into N/S buffers
             iBL_BLCR_iBR(cbl+i) = iBL_iCL_iTL(bl+i)
             iTL_TLCR_iTR(ctl+i) = iBL_iCL_iTL(tl+i)
          enddo
        endif
        if(.not. east) then
          call MPI_wait(from_east,status,ierr)  ! wait for inbound East -> West message to complete
          do i = 0,halox*haloy*nk-1
             iBL_BLCR_iBR(cbr+i) = iBR_iCR_iTR(br+i)
             iTL_TLCR_iTR(ctr+i) = iBR_iCR_iTR(tr+i)
          enddo
        endif
!      write(rpn_u,100)'DBG',nwds_ns,minx,maxx,miny,maxy,
!     %                size(cTL_iTLCR_cTR),size(cBL_iBLCR_cBR),
!     %                size(iBL_BLCR_iBR),size(iTL_TLCR_iTR)
!
!       step 6 post North/South non blocking receives, post North/South non blocking sends
!
        if(.not. north) then
          ! recv from north neighbor unless I am north PE
          call MPI_IRECV(cTL_iTLCR_cTR,nwds_ns,MPI_INTEGER,northpe, &
               100000+northpe,PE_DEFCOMM,from_north,ierr)          ! sender was northpe therefore tag is northpe+100000
        endif
        if(.not. south) then
          ! recv from south neighbor unless I am south PE
          call MPI_IRECV(cBL_iBLCR_cBR,nwds_ns,MPI_INTEGER,southpe, &
               southpe,PE_DEFCOMM,from_south,ierr)                 ! sender was southpe therefore tag is southpe
        endif
        if(.not. south) then
          ! send to south neighbor unless I am south PE
          call MPI_ISEND(iBL_BLCR_iBR,nwds_ns,MPI_INTEGER,southpe, &
               100000+pe_medomm,PE_DEFCOMM,to_south,ierr)          ! tag is PE grid ordinal of sender+100000
        endif
        if(.not. north) then
          ! send to north neighbor unless I am north PE
          call MPI_ISEND(iTL_TLCR_iTR,nwds_ns,MPI_INTEGER,northpe, &
               pe_medomm,PE_DEFCOMM,to_north,ierr)                 ! tag is PE grid ordinal of sender
        endif
!
!       step 7  put East/West outer halos into array
!               (east/west boundary and periodx false, do not fill east/west outer halo)
!
        m = bl
        m1 = bl
        m2 = cl
        m3 = tl
        do k=1,nk   ! BL/BL part
         do j=1,haloy
          do i=1,halox
!           if(.not.east)g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iBR part
           g(ni+i,j,k)    = iand(east_m,  g(ni+i,j,k)) +  &
                            iand(east_m_n,iBR_iCR_iTR(m1))
!           if(.not.west)g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iBL part
           g(i-halox,j,k) = iand(west_m,  g(i-halox,j,k)) + &
                            iand(west_m_n,iBL_iCL_iTL(m1))
           m1=m1+1
!           if(.not.east)g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iTR part
           g(ni+i,nj-haloy+j,k)    =  &
                            iand(east_m,  g(ni+i,nj-haloy+j,k)) + &
                            iand(east_m_n,iBR_iCR_iTR(m3))
!           if(.not.west)g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iTL part
           g(i-halox,nj-haloy+j,k) = &
                            iand(west_m,  g(i-halox,nj-haloy+j,k)) + &
                            iand(west_m_n,iBL_iCL_iTL(m3))
           m3=m3+1
           m=m+2
          enddo
         enddo
        enddo
        do k=1,nk  ! CL/CR part
         do j=1+haloy,nj-haloy
          do i=1,halox
!           if(.not.east)g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iCR part
           g(ni+i,j,k)    = iand(east_m,  g(ni+i,j,k)) + &
                            iand(east_m_n,iBR_iCR_iTR(m2))
!           if(.not.west)g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iCL part
           g(i-halox,j,k) = iand(west_m,g(i-halox,j,k)) + &
                            iand(west_m_n,iBL_iCL_iTL(m2))
           m=m+1
           m2=m2+1
          enddo
         enddo
        enddo
!        do k=1,nk  ! TL/TR part
!         do j=nj-haloy+1,nj
!          do i=1,halox
!           if(.not.east)g(ni+i,j,k)   =iBR_iCR_iTR(m)      ! iTR part
!           g(ni+i,j,k)    = iand(east_m,  g(ni+i,j,k)) + 
!     %                      iand(east_m_n,iBR_iCR_iTR(m3))
!           if(.not.west)g(i-halox,j,k)=iBL_iCL_iTL(m)      ! iTL part
!           g(i-halox,j,k) = iand(west_m,  g(i-halox,j,k)) +
!     %                      iand(west_m_n,iBL_iCL_iTL(m3))
!           m=m+1
!           m3=m3+1
!          enddo
!         enddo
!        enddo
        if(m /= 1+halox*nj*nk) print *,'OUCH EW recv'
!
!       step 8 wait for inbound North -> South messages to complete 
!              and put into North outer halo
!
        if(.not. north) then   ! north boundary and periody is false, nothing to do
          call MPI_wait(from_north,status,ierr)
          m=1
          m1 = cbl
          m2 = iblcr
          m3 = cbr
!          if(west) then        ! west boundary and periodx false, do not fill west outer halo
!            m=m+nk*halox*haloy
!          else
            do k=1,nk
            do j=nj+1,nj+haloy
            do i=1,halox
!               g(i-halox,j,k) = cTL_iTLCR_cTR(m1)  ! cTL part
               g(i-halox,j,k) = iand(west_m,  g(i-halox,j,k)) + &
                                iand(west_m_n,cTL_iTLCR_cTR(m1))
               m1=m1+1
!               g(ni+i,j,k) = cTL_iTLCR_cTR(m1)     ! cTR part
               g(ni+i,j,k) = iand(east_m,  g(ni+i,j,k)) + &
                             iand(east_m_n,cTL_iTLCR_cTR(m3))
               m3=m3+1
               m=m+2
            enddo
            do i=1,ni
               g(i,j,k) = cTL_iTLCR_cTR(m2)  ! iTLCR part
               m2=m2+1
               m=m+1
            enddo
            enddo
            enddo     
!          endif
!          do k=1,nk
!          do j=nj+1,nj+haloy
!          do i=1,ni
!             g(i,j,k) = cTL_iTLCR_cTR(m2)  ! iTLCR part
!             m=m+1
!             m2=m2+1
!          enddo
!          enddo
!          enddo
!!          if(.not.east) then    ! east boundary and periodx false, do not fill east outer halo
!            do k=1,nk  ! cTR part
!            do j=nj+1,nj+haloy
!            do i=1,halox
!!               g(ni+i,j,k) = cTL_iTLCR_cTR(m)
!               g(ni+i,j,k) = iand(east_m,  g(ni+i,j,k)) +
!     %                       iand(east_m_n,cTL_iTLCR_cTR(m3))
!               m3=m3+1
!               m=m+1
!            enddo
!            enddo
!            enddo
!!          else
!!            m=m+nk*halox*haloy
!!          endif
          if(m /= 1+haloy*(ni+2*halox)*nk) print *,'OUCH N recv'
        endif
!
!       step 9 wait for inbound North <- South messages to complete 
!              and put into South outer halo
!

        if(.not. south) then   ! south boundary and periody is false, nothing to do
          call MPI_wait(from_south,status,ierr)
          m=1
          m1 = cbl
          m2 = iblcr
          m3 = cbr
!          if(west) then        ! west boundary and periodx is false, do not fill west outer halo
!            m=m+nk*halox*haloy
!          else
            do k=1,nk 
            do j=1-haloy,0
            do i=1,halox
!               g(i-halox,j,k) = cBL_iBLCR_cBR(m3) ! cBL part
               g(i-halox,j,k) = iand(west_m,  g(i-halox,j,k)) + &
                                iand(west_m_n,cBL_iBLCR_cBR(m1))
               m1=m1+1
!               g(ni+i,j,k) = cBL_iBLCR_cBR(m3)    ! cBR part
               g(ni+i,j,k) = iand(east_m,  g(ni+i,j,k)) + &
                             iand(east_m_n,cBL_iBLCR_cBR(m3))
               m3=m3+1
               m=m+3
            enddo
            do i=1,ni
               g(i,j,k) = cBL_iBLCR_cBR(m2)  ! iBLCR part
               m2=m2+1
               m=m+1
            enddo
            enddo
            enddo
!          endif
!          do k=1,nk
!          do j=1-haloy,0
!          do i=1,ni
!             g(i,j,k) = cBL_iBLCR_cBR(m2)  ! iBLCR part
!             m2=m2+1
!             m=m+1
!          enddo
!          enddo
!          enddo
!!          if(.not.east) then    ! east boundary and periodx false, do not fill east outer halo
!            do k=1,nk  ! cBR part
!            do j=1-haloy,0
!            do i=1,halox
!!               g(ni+i,j,k) = cBL_iBLCR_cBR(m3)
!               g(ni+i,j,k) = iand(east_m,  g(ni+i,j,k)) +
!     %                       iand(east_m_n,cBL_iBLCR_cBR(m3))
!               m3=m3+1
!               m=m+1
!            enddo
!            enddo
!            enddo
!!          else
!!            m=m+nk*halox*haloy
!!          endif  
          if(m /= 1+haloy*(ni+2*halox)*nk) print *,'OUCH S recv'
        endif
!
!       step 10 wait for all outbound messages to complete
!
        if(.not.east)call MPI_wait(to_east,status,ierr)
        if(.not.west)call MPI_wait(to_west,status,ierr)
        if(.not.south)call MPI_wait(to_south,status,ierr)
        if(.not.north)call MPI_wait(to_north,status,ierr)
!
        return
!
        entry xch_halo(g,minx,maxx,miny,maxy,ni,nj,nk,halox,haloy,periodx,periody)
        globalni=ni
        polarrows=0

        goto 1
        end
