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

!**s/r inp_Nes_openFST - Open+link all fst input files valid at F_datev
!                    and determine vertical structure (kind)
      subroutine INs_Nes_openFST ( F_datev, F_err )
      use, intrinsic :: iso_fortran_env
      use clib_itf_mod
      use INs_base
      use rmn_fst24
      implicit none

      character(len=*), intent(IN ) :: F_datev
      integer         , intent(OUT) :: F_err

      logical, external :: INs_shared_mem, INs_gz3d
      logical :: shared_mem_L
      logical :: success
      character(len=2048) :: fn,root,mesg
      integer :: i,j,err,err_code,unf
      integer :: lislon1,lislon2,lislon3,unit1,unit2,unit3
      integer, parameter :: nlis = 1024
      type(fst_file)   :: file(10)
      type(fst_query)  :: query
      type(fst_record) :: recs1(nlis),recs2(nlis), recs3(nlis) 
!
!-----------------------------------------------------------------------
!
      if (INs_1o1_L) err = fstopc('MSGLVL','INFORM',RMN_OPT_SET)
      fst(:)%fnt= -1 ; fst(:)%type= '0' ; SRC_GZ_L= .false.

      fn= trim(Path_input_S)//'/ANALYSIS'
      if (file(1)%open(trim(fn),'RND+OLD+R/O')) then
         Nes_fst(1) = file(1)
         unit1=file(1)%get_unit()
         fst(1)%fnt = unit1
         fst(1)%type= 'A'
         query = file(1)%new_query(datev=Inp_cmcdate,nomvar= 'TT ')
         lislon1 = query%find_all(recs1)
         if (INs_1o1_L) write (6,'("FILE: ",a," opened on unit:",i6)') trim(fn),unit1
      end if

      fn= trim(Path_input_S)//'/GEOPHY/Gem_geophy.fst'
      if (file(2)%open(trim(fn),'RND+OLD+R/O')) then 
            Nes_fst(2) = file(2)
            unit2=file(2)%get_unit()
            fst(2)%fnt = unit2
            fst(2)%type= 'G'
            if (INs_1o1_L) write (6,'("FILE: ",a," opened on unit:",i6)') trim(fn),unit2
            query = file(2)%new_query(datev=Inp_cmcdate,nomvar= 'TT ')
            lislon2 = query%find_all(recs2)
         endif
         
      unf= 0
      fn= trim(Path_input_S)//'/CLIMATO'
      ! Climato is a directory: don't know what to do
      unf= 0
      fn= trim(Path_input_S)//'/MODEL_INPUT'
      ! Don't know what to do with MODEL_INPUT

      F_err = 0
      Ins_handle= -1 ; Inp_nfiles= 0 ; i= 0

      if ( any(SRL(:)%src == 'R') ) then
         Inp_datev= F_datev
         call datp2f ( Inp_cmcdate, F_datev )

         root=trim(Path_input_S)//'/MODEL_INREP/VALID_'//trim(F_datev)
         err= clib_fileexist (trim(root)//'/content')

         if (err < 0) root=trim(Path_input_S)//&
                    '/MODEL_ANALYSIS/VALID_'//trim(F_datev)
         fn = trim(root)//'/content'
         unf= 0
         if (fnom( unf,trim(fn),'SEQ+FMT+OLD',0 ) /= 0) unf= 0
         if (unf == 0) goto 33
         read (unf,*,end=33) Inp_nfiles
         if (Inp_nfiles == 0) goto 33
         allocate (Inp_list_files(Inp_nfiles))

 55      read (unf,'(a)',end=33) fn
         i = i+1
         fn= trim(root)//'/'//trim(fn)

         if (.not. Inp_list_files(i)%open(trim(fn),'RND+OLD+R/O')) then
            F_err = -1
         end if
         if (F_err == 0) goto 55
  
 33      if ((Inp_nfiles == 0).or.(i /= Inp_nfiles)) F_err= -1
         if (unf > 0) err= fclos(unf)
      
         if (F_err == 0) then
            success = fst24_link(Inp_list_files(1:Inp_nfiles))
            Inp_file = Inp_list_files(1)
            Ins_handle=Inp_list_files(1)%get_unit()
            fst(3)%fnt = Ins_handle
            fst(3)%type= 'R'

            query = Inp_list_files(1)%new_query(datev=Inp_cmcdate,nomvar= 'TT ')
            lislon3 = query%find_all(recs3)
            unit3=fst(3)%fnt

            if (INs_1o1_L) write (6,'("FILE: ",a," opened on unit:",i6)') trim(root),Ins_handle
         endif
         if (F_err<0) return
      endif

      do i=1,INs_nreq
         do j=1,size(fst%type)
            if (fst(j)%type /= '0') then
               if (SRL(i)%src == fst(j)%type) then
                  SRL(i)%unf= fst(j)%fnt
                  cycle
               endif
            endif
         end do
      end do

! Establishing memory on first occurence of
! TT in an 'A' or 'R' input file
      shared_mem_L= .false.
      do i=1,INs_nreq
         if (shared_mem_L) exit
         if ( SRL(i)%src == 'A') then
            shared_mem_L= INs_shared_mem (unit1,recs1(1)%ni,recs1(1)%nj,lislon1, 'TT')
         elseif (SRL(i)%src == 'R') then
            shared_mem_L= INs_shared_mem (unit3,recs3(1)%ni,recs3(1)%nj,lislon3, 'TT')
         endif

         if (shared_mem_L) unf=SRL(i)%unf
      end do
      SRL(:)%src='' ! Not used in INs_read for nesting data
      
      SRC_GZ_L= .false.
      if (shared_mem_L) SRC_GZ_L= INs_gz3d () !Obtain GZ if possible
      if (.not.SRC_GZ_L) then
         GZcaract(:)= 0
         INs_nreq=0
         F_err= -1
      endif
      
      call query%free()
!
!-----------------------------------------------------------------------
!         
      return
      end subroutine INS_Nes_openFST
