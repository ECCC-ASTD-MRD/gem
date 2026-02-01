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
      subroutine OUTs_process_data (F_deb, F_fin, F_nplans)
      use iso_c_binding
      use vGrid_Descriptors
      use MiMd
      use IOs
      use OUTs
      use omp_timing
      implicit none

      integer, intent(IN) :: F_deb, F_fin, F_nplans

      character(len=1024) fn
      character(len=4) stag
      character(len=4) prefix,dumc4
      logical bool
      integer ORIGIN_COUNT, TARGET_COUNT, TARGET_RANK
      integer :: i,j,k,n,nko,nvar,skip,ii,jj,kk,offs,sk1,k01,v_stag,h_stag
      integer :: mpx,nklocal,irest,kstart,kend,k0,dim
      integer :: indo(10000),err,nfstecr,len,tag,cnt,ns
      integer :: sfc_lcl,sfc_i0,sfc_in
      integer(C_INTPTR_T) :: TARGET_DISP
      real :: gem_data(client_Hplane*G_nk*10)
      real, target, dimension(1) :: level_ground
      real, dimension (:,:), allocatable :: dataH
      type (sfc_var), dimension(:), allocatable :: vsfc
      logical :: bottom_L
!     
!--------------------------------------------------------------------
!
      if (F_nplans<1)then
         return
      endif
      if (OUTs_1o1_L) print*, 'Proceeding with marker:', F_deb, F_fin
      
      call gtmg_start ( 20, 'GET_data', 13)
      if (Lun_out>0) call clock ( Lun_out, 'Recieving data', .false. )
      dim= client_Hplane*F_nplans
      allocate (dataH(dim,client_pestart:client_peend))
      allocate (vsfc(8000))
      do n=client_pestart,client_peend
         tag= 1001
         call MPI_recv ( dataH(1,n), dim, MPI_REAL, n, &
                         tag, MPI_COMM_WORLD, MPI_STATUSES_IGNORE, err)
      end do
      call gtmg_stop ( 20 )
      if (Lun_out>0) call clock (Lun_out, 'Recieving ... DONE', .false.)

      skip=0 ; k0=0 ; ns=0
      call gtmg_start ( 21, 'Processing', 13)

 987  if (skip>=F_fin) goto 888
      ns= ns+1
      prefix= TRANSFER(metaG(skip+1), dumc4)
      call up2low (prefix, dumc4) ; prefix=dumc4
      call OUTs_whichFST (Out_unf, prefix)
      
      Out_i0  = metaG(skip+2)
      Out_in  = metaG(skip+3)
      Out_j0  = metaG(skip+4)
      Out_jn  = metaG(skip+5)
      Out_reduc_L = TRANSFER(metaG(skip+6), bool)
      nvar= metaG(skip+7)
      skip= skip+7 ; cnt=0

      LOOP_VAR: do i= 1, nvar
         bottom_L=.false.
         call gtmg_start ( 31, 'Assembling', 21)
         Out_nomvar= TRANSFER(metaG(skip+1), dumc4)
         stag      = TRANSFER(metaG(skip+2), dumc4)
         call low2up (stag ,dumc4)
         stag  = dumc4
         Out_kind= metaG(skip+3)
         Out_nbit= metaG(skip+4)
         nko     = metaG(skip+5)
         indo(1:nko) = metaG(skip+6:skip+6+nko-1)
         nullify (levels)
         if (stag(2:2) == 'S') levels => Ver_i
         if (stag(2:2) == 'M') levels => Ver_hybM
         if (stag(2:2) == 'T') levels => Ver_hybT
         if (stag(2:2) == 'G')then
            if(Out_kind == 1 .or. Out_kind == 5)then
               level_ground(1)=1.
            else if(Out_kind == 4 )then
               level_ground(1)=0.
            else
               if (OUTs_1o1_L) print*,'WARNING, ground level is not defined for kind ',&
                    Out_kind,', setting level to 0.0'
               level_ground(1)=0.
            endif
            levels => level_ground
         endif
         if (stag(2:2) == 'P') levels => Level_allpres
         if (stag(2:2) == 'H') levels => Level_allheights
         if (stag(3:3) == 'D') then
            if (stag(2:2) == 'M') levels => hybM_diag
            if (stag(2:2) == 'T') levels => hybT_diag
            nko= 1 ; indo(1) = G_nk+1
         endif
         if(OUTs_vgrid_usr_L(stag))then
            read(stag(2:2),*)v_stag
            ! Find matching user vgrid
            do k=1,size(vgd_usr)
               if(vgd_usr(k)%stag == v_stag)exit
            end do
            ! User level
            if(vgd_usr(k)%vcode .ne. 1002)then
               print*,'WARNING, user level ',v_stag,', of vcode ',vgd_usr(k)%vcode,' not supported, skipping'
               cycle
            endif
            levels => vgd_usr(k)%levels
            if(nko .ne. size(levels))then
               ! User asked for bottom nko levels
               bottom_L=.true.
               levels(1:nko) => vgd_usr(k)%levels(size(vgd_usr(k)%levels)-nko+1:size(vgd_usr(k)%levels))
            endif
         endif
         if (nko == 1 .and. (.not. bottom_L) ) then
            cnt=cnt+1
            vsfc(cnt)%nv   = Out_nomvar
            !if(Out_kind == 5 .or. Out_kind == 21)then
            !   ! For hyb pressure (kind=5) or hyb height (kind=21) the
            !   ! diag levels have kind 4, height with respect to ground level in m
            !   vsfc(cnt)%knd  = 4
            !else
               vsfc(cnt)%knd  = Out_kind
            !endif
            vsfc(cnt)%nbits= Out_nbit
            vsfc(cnt)%indx = i
            vsfc(cnt)%k0   = k0
            vsfc(cnt)%stag = stag
            vsfc(cnt)%skip = skip
            if ( stag(2:2) == 'S' ) then
               vsfc(cnt)%lvl= 0.
            else
               vsfc(cnt)%lvl= levels(indo(1))
            endif
            if(OUTs_hor_int_L(stag))then
               if(.not.UTL_get_usr_int_info(vsfc(cnt)%usr_int_info,vsfc(cnt)%nv,OUT_usr_int_info))then
                  print*,'Error in OUTs_process_data with, UTL_get_usr_int_info'
                  ! TODO handle error gracefully
               endif               
            endif
            skip= skip+6+nko-1
            k0  = k0+nko
            cycle
         endif
         skip  = skip+6+nko-1
         if(OUTs_hor_int_L(stag))then
            ! Find matching horizontal
            read(stag(4:4),*)h_stag
            do j=1,size(hgd_usr)
               if(hgd_usr(j)%usr_grid_index == h_stag)then
                  call OUTs_wrtref_usr (stag,hgd_usr(j))
                  exit
               end if
            enddo
         else   
            call OUTs_wrtref (stag)
         endif
         
         do n=client_pestart,client_peend
            offs=n-IOS_YIN*IOS_couleur-clients_npes(2,gem_id)+1
            call lcl2glb (dataH(1,n), model_gindx(1,offs), k0, nko)
         end do
         
         call MPI_barrier (MY_WORLD_COMM,err)
         call gtmg_stop ( 31 )
         call gtmg_start ( 32, 'FSTECR', 21)

         call splitW ( myproc_IOS, numproc_IOS,nko,1,&
                       sfc_lcl, kstart, kend )
         call OUTs_fstecr (levels,indo,nko,kstart,kend,k0,stag)

         call MPI_barrier (MY_WORLD_COMM,err)
         call gtmg_stop ( 32 )
         k0= k0+nko
      enddo LOOP_VAR
      
      if (cnt>=1) then
         call gtmg_start ( 33, 'OneLevel', 21)
         call splitW ( myproc_IOS, numproc_IOS,cnt,1,&
                       sfc_lcl, sfc_i0, sfc_in)
         do i=1,cnt
         do n=client_pestart,client_peend
            offs=n-IOS_YIN*IOS_couleur-clients_npes(2,gem_id)+1
            call lcl2glb_sfc (dataH(1,n), model_gindx(1,offs), vsfc(i)%k0, i)
         end do
         end do
         call MPI_barrier (MY_WORLD_COMM,err)
         call OUTs_fstecr_sfc (vsfc,sfc_i0, sfc_in,'MS'//stag(3:4))
         call MPI_barrier (MY_WORLD_COMM,err)
         call gtmg_stop ( 33 )
      endif
      goto 987
      
 888  deallocate (dataH,vsfc)
      call gtmg_stop ( 21 )

      if (OUTs_1o1_L) print*, 'Marker:', F_deb, F_fin, ' ...DONE'
      if (Lun_out>0) call clock ( Lun_out, 'Processing ...DONE', .false. )
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_process_data

      subroutine lcl2glb (src, lcl_indx, F_k0, F_nk)
      use IOs
      use OUTs
      implicit none
      
      integer, intent(IN) :: lcl_indx(4), F_k0, F_nk
      real, intent(IN) :: src(*)

      integer i,j,k,cnt
!
!--------------------------------------------------------------------
!
      if ( F_nk > ubound(IOs_glbdata,3) ) then
         if (Lun_out>0) &
         print*, 'Insufficient storage in array IOs_glbdata --ABORT',&
                  F_nk, ubound(IOs_glbdata,3)
         stop
      endif
      
      do k= 1, F_nk
         cnt= F_k0*client_Hplane + (k-1)*client_Hplane
         do j= lcl_indx(3), lcl_indx(4)
            do i= lcl_indx(1), lcl_indx(2)
               cnt=cnt+1
               IOs_glbdata(i,j,k,IOS_couleur+1)= src(cnt)
            end do
         end do
      end do
!
!--------------------------------------------------------------------
!
      return
      end subroutine lcl2glb
      
      subroutine lcl2glb_sfc (src, lcl_indx, F_k0, F_k)
      use IOs
      use OUTs
      implicit none
      
      integer, intent(IN) :: lcl_indx(4), F_k0, F_k
      real, intent(IN) :: src(*)

      integer i,j,cnt
!
!--------------------------------------------------------------------
!
      if ( F_k > ubound(IOs_glbdata,3) ) then
         if (Lun_out>0) &
         print*, 'Insufficient storage in array IOs_glbdata --ABORT',&
                  F_k, ubound(IOs_glbdata,3)
         stop
      endif
      
      cnt= F_k0*client_Hplane
      do j= lcl_indx(3), lcl_indx(4)
         do i= lcl_indx(1), lcl_indx(2)
            cnt=cnt+1
            IOs_glbdata(i,j,F_k,IOS_couleur+1)= src(cnt)
         end do
      end do
!
!--------------------------------------------------------------------
!
      return
      end subroutine lcl2glb_sfc
