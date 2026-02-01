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

!**s/r inp_Spn_openFST - Open+link all fst input files valid at F_datev
!                        and determine vertical structure (kind)
      subroutine INs_Spn_openFST ( F_datev, F_err )
      use, intrinsic :: iso_fortran_env
      use clib_itf_mod
      use INs_base
      use rmn_fst24
      implicit none

      character(len=*), intent(IN ) :: F_datev
      integer         , intent(OUT) :: F_err

      integer :: n1,n2,n3
      integer :: DISP_UNIT
      integer(KIND=MPI_ADDRESS_KIND) :: WINSIZE
      real :: dim
      type(C_PTR), save :: basepntr2
      
      logical :: success,found_L,done=.false.
      character(len=2048) :: fn,filename,root,spn_files,cmd
      integer :: i,j,err,unf,unit3,lislon
      integer, parameter :: nlis = 1024
      type(fst_record) :: my_record,recs3(nlis)
      type(fst_query)  :: my_query
!
!-----------------------------------------------------------------------
!
      if (INs_1o1_L) err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)
      fst(:)%fnt= -1 ; fst(:)%type= '0'

      F_err= -1 ; Ins_handle= -1 ; found_L= .false.

      Inp_datev= F_datev
      call datp2f ( Inp_cmcdate, F_datev )
      root= trim(Path_input_S)//'/NUDGE/'//INs_runstrt_S(1:8)//INs_runstrt_S(10:11)//'/'
      spn_files = trim(IOs_path_work_S)//'/iau_files'

      if (.not.done .or. INs_Spn_listfst_L) then
         if (INs_1o1_L) then
            cmd='cd '//trim(root)//&
            ' ; find ./ -maxdepth 1 -type f | sort > '//trim(spn_files)
            call system(trim(cmd))
         endif
         done= .true.
      endif

      call MPI_barrier (MY_WORLD_COMM,err)

      unf= 0
      if (fnom( unf,trim(spn_files),'SEQ+FMT+OLD',0 ) /= 0) unf= 0
      if (unf == 0) goto 33
 55   read (unf,'(a)',end=33) fn
      filename= trim(root)//trim(fn)
      if (SPN_file%open(trim(filename),'RND+OLD+R/O')) then
         Ins_handle=SPN_file%get_unit()
         my_query = SPN_file%new_query(nomvar = 'TT', datev=Inp_cmcdate)
         if (my_query%find_next(my_record)) then
            if (Lun_out > 0) print*, 'FOUND TT   valid at: ',&
                        F_datev, ' in file ',trim(fn)
            found_L= .true.
            goto 33
         else
            success= SPN_file%close()
         endif
         call my_query%free()
         goto 55
      endif
 33   if (found_L) F_err= 0
      err = fclos(unf)

      if (F_err < 0) return

      Inp_file = SPN_file
      unit3= SPN_file%get_unit()
      my_query = SPN_file%new_query(datev=Inp_cmcdate,nomvar= 'TT ')
      lislon   = my_query%find_all(recs3)
      call my_query%free()
         
      if (Lun_out>0) then
         write (Lun_out,'(" FILE: ",a," opened with unit: ",i6)')&
                                       trim(fn),Ins_handle
      endif
      if (lislon<=0) &
         stop 'Could NOT determine Input dimensions with var= TT'

      call INs_getvgd ()

      INs_nia= recs3(1)%ni
      INs_nja= recs3(1)%nj
      INs_nka= lislon
      if (INs_1o1_L) print*, 'Analysis at ',Inp_datev,&
              ' dimensions: ',INs_nia,INs_nja,INs_nka,&
              ' from TT in unit: ', unit3
      n1= INs_nid
      n2= INs_njd
      n3= (INs_nreq+2)*(INs_nka+2)
      INs_Spn_nkd= n3
      
      disp_unit = 4
      dim = 0.
      if (INs_hostmyproc == 0) dim= real(INs_hord)*real(n3*2)*real(disp_unit)

      WINSIZE = dim
      call MPI_Win_allocate_shared(WINSIZE, disp_unit, MPI_INFO_NULL,&
                                   INs_host, basepntr2, SPN_win, err)
      call MPI_Win_shared_query(SPN_win, MPI_PROC_NULL, WINSIZE  ,&
                                disp_unit, basepntr2, err)
      n2= n1*n2*n3
      call C_F_POINTER (basepntr2,SPN_DATA,[n2,2] )
      call MPI_barrier (MY_WORLD_COMM, err)
      
      allocate (DIP1(n3*2))
!
!-----------------------------------------------------------------------
!         
      return
      end subroutine INS_Spn_openFST
