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
module Openfst
      
   use, intrinsic :: iso_fortran_env
   use rmn_fst24
   use vGrid_Descriptors
   implicit none
   public
   save

contains
      subroutine open_fst (F_datev, F_file, F_vgd_src,&
                           F_TT_list, F_nka_tt, F_dir,F_updlst_L)
      use inp_base
      use path
      use ptopo
      implicit none

      character(len=*), intent(IN) :: F_datev, F_dir
      logical, intent(IN ) :: F_updlst_L
      integer, intent(OUT) :: F_nka_tt
      integer, dimension(:), pointer,intent(OUT) :: F_TT_list
      type(vgrid_descriptor), intent(OUT) :: F_vgd_src
      type(fst_file), intent(OUT) :: F_file

      character(len=2048) root,list_files,cmd,filename,fn
      logical :: success,found_L,done=.false.
      integer :: err,err_code,indx,cmcdate,unf,unit,kind,version
      integer, parameter :: n123_dim=3 , nlis = 1024
      integer n123(n123_dim), liste_sorted(nlis)
      real(kind=REAL64) :: pref_a_8
      real(kind=REAL64), pointer :: vtbl_8(:,:,:)
      type(fst_record) :: my_record
      type(fst_query)  :: my_query
      type(fst_record) :: recs(nlis) 
!     
!     ---------------------------------------------------------------
!
      found_L= .false. ; err_code= 0
      root= trim(Path_input_S)//'/'//trim(F_dir)
      indx= index(F_dir,"/")
      if (indx<1) indx=len(F_dir)
      list_files = F_dir(1:indx-1)
      list_files = trim(list_files)//'_files'
      call datp2f ( cmcdate, F_datev )

      if (.not.done .or. F_updlst_L) then
         if (Inp_iome >= 0) then
            cmd='ls -1 '//trim(root)//' > '//trim(list_files)
            call system(cmd)
         endif
         done= .true.
      endif
      
      if (Inp_iome >= 0) then
         unf= 0 ; err_code= -1
         if (fnom( unf,trim(list_files),'SEQ+FMT+OLD',0 ) /= 0) unf= 0
         if (unf == 0) goto 33
 55      read (unf,*,end=33) filename
         fn= trim(root)//'/'//trim(filename)

         if (F_file%open(trim(fn),'RND+OLD+R/O')) then
            unit=F_file%get_unit()
            my_query = F_file%new_query(nomvar = 'TT', datev=cmcdate)
            if (my_query%find_next(my_record)) then
               if (Lun_out > 0) print*, 'FOUND TT valid at: ', F_datev, ' in file ',trim(filename)
               found_L= .true.
               goto 33
            else
               success= F_file%close()
           endif
           call my_query%free()
           goto 55
         endif
 33      if (found_L) err_code= 0
      endif

      call gem_error ( err_code, 'open_fst', &
                       'Unable to find valid data in FST files' )

      nullify (vtbl_8) ; n123= -1
      if (Inp_iome >= 0) then
         err= vgd_new ( F_vgd_src, unit=F_file%get_unit(), &
                        format='fst', ip1=-1, ip2=-1 )
         if (err == 0) then
            err= vgd_get ( F_vgd_src, 'VTBL', vtbl_8, quiet=.true.)
            n123(1:3) = ubound(vtbl_8)
         end if
      endif

      call rpn_comm_bcast ( n123, n123_dim, "MPI_INTEGER", Inp_iobcast, &
                            "grid", err )

      kind= -1
      if (n123(1) > 0) then
         if (Inp_iome /= 0) allocate(vtbl_8(n123(1),n123(2),n123(3)))
         call rpn_comm_bcast ( vtbl_8,size(vtbl_8), &
             "MPI_DOUBLE_PRECISION", Inp_iobcast, "grid", err )
         if (Inp_iome /= 0) err= vgd_new ( F_vgd_src, vtbl_8 )
         deallocate (vtbl_8)
         err = vgd_get ( F_vgd_src, key='KIND',value=kind    )
         err = vgd_get ( F_vgd_src, key='VERS',value=version )
         if ( (kind == 5) .or. &
             ((kind == 1).and.(version == 3)) ) then
            err= vgd_get ( F_vgd_src, key='PREF',value=pref_a_8 )
         end if
         
         if (Lun_out > 0) write(lun_out,9000) kind, version
      else
         call gem_error ( -1, 'open_fst', &
                       'Unable to determine vertical structure')
      end if

      if (kind == 21) call gem_error ( -1, 'open_fst', &
                   'input on heights NOT yet available')

    ! using the Inp system
      Inp_file   = F_file
      Inp_cmcdate= cmcdate

      if (Inp_iome >= 0) then
         my_query = F_file%new_query(datev=cmcdate,&
                              nomvar='TT',ip1=-1)
         F_nka_tt= my_query%find_all(recs)
         call my_query%free()
         if (F_nka_tt > 1) then
            call record_sort_ip1 (recs,liste_sorted,F_nka_tt)
         endif
      endif
      call rpn_comm_bcast ( F_nka_tt, 1, "MPI_INTEGER", Inp_iobcast,&
                               "grid", err )
      call rpn_comm_bcast ( liste_sorted, nlis, "MPI_INTEGER",&
                               Inp_iobcast, "grid", err )
      allocate (F_TT_list(F_nka_tt))
      F_TT_list(1:F_nka_tt) = liste_sorted(1:F_nka_tt)

 9000 format(' Input vertical description obtained with vgd_get: '/&
            ' kind, version= ',i5,',',i5)
!                       
!     ---------------------------------------------------------------
!
      return
      end subroutine open_fst
      
end module Openfst
