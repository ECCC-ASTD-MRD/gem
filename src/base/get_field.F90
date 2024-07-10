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

!**s/r get_field - Read field F_nomvar_S from file F_filename_S and proceed
!                  with horizontal interpolation using scheme F_inttyp_S

      subroutine get_field ( F_dest, ni, nj, F_nomvar_S, F_filename_S, &
                             F_inttyp_S, F_xfi, F_yfi, F_status )
      use clib_itf_mod
      use glb_ld
      use glb_pil
      use hgc
      use lun
      use rmn_fst24

      implicit none
#include <arch_specific.hf>
!
      character(len=*), intent(in) :: F_nomvar_S, F_filename_S
      character(len=*), intent(out) :: F_inttyp_S
      integer, intent(in) :: ni,nj
      integer, intent(out) :: F_status
      real, dimension(ni,nj), intent(out) :: F_dest
      real, dimension(ni), intent(in) ::  F_xfi
      real, dimension(nj), intent(in) ::  F_yfi

      integer,external :: ezgdef_fmem,ezdefset,ezsetopt,ezsint,ezqkdef,ezgetopt,&
           ezsetival,ezgetival,ezgetval,samegrid_gid

      type(fst_file)   :: file
      type(fst_record) :: record, tic, tac
      type(fst_query)  :: query
      logical          :: success
      
      real(kind = real32), dimension(:  ), pointer :: ax,ay
      real(kind = real32), dimension(:,:), pointer :: source

      character(len=32):: onesubgrid_S, oldsubgrid_S
      integer unf,key,err,ierr,subid,oldsubid
      integer sgid,dgid,nic,njc

!-----------------------------------------------------------------------
!
      F_status = -1

      subid= -1 ; oldsubid= -1

      err = clib_isfile(trim(F_filename_S))
      if ( err < 0 ) then
         write(app_msg,9001) trim(F_filename_S)
         call lib_log(APP_LIBDYN,APP_ERROR,app_msg)
         goto 999
      end if

      success = file%open(trim(F_filename_S),'RND+OLD+R/O')
      if (.not. success) then
         write(app_msg,9002) trim(F_filename_S),'RND+OLD+R/O'
         call lib_log(APP_LIBDYN,APP_ERROR,app_msg)
         goto 999
      endif

   write(app_msg,1001) trim(F_filename_S)
      call lib_log(APP_LIBDYN,APP_INFO,app_msg)

      query = file%new_query(nomvar=F_nomvar_S)
      success = query%find_next(record)
      if (.not. success) then
         write(app_msg,9004) F_nomvar_S
         call lib_log(APP_LIBDYN,APP_ERROR,app_msg)
         goto 999
      end if

      success=record%read()

      dgid = ezgdef_fmem (ni, nj , 'Z', 'E', Hgc_ig1ro, &
                 Hgc_ig2ro, Hgc_ig3ro, Hgc_ig4ro, F_xfi , F_yfi )

      if (any(record%grtyp==(/'Z','z'/))) then

         query = file%new_query(ip1=record%ig1,ip2=record%ig2,nomvar=">>  ")
         success = query%find_next(tic)
         if (.not. success) then
            write(app_msg,9004) ">>"
            call lib_log(APP_LIBDYN,APP_ERROR,app_msg)
            goto 999
         end if
         success=tic%read()
         call tic%get_data_array(ax) 

         query = file%new_query(ip1=record%ig1,ip2=record%ig2,nomvar="^^  ")
         success = query%find_next(tac)
         if (.not. success) then
            write(app_msg,9004) "^^"
            call lib_log(APP_LIBDYN,APP_ERROR,app_msg)
            goto 999
         end if
         success=tac%read()
         call tac%get_data_array(ay)

         sgid = ezgdef_fmem (record%ni, record%nj, record%grtyp, tac%grtyp, &
                                   tac%ig1, tac%ig2, tac%ig3, tac%ig4, ax, ay)

         if (samegrid_gid (sgid, Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro,&
                           F_xfi, F_yfi, ni, nj) > -1) then
            F_inttyp_S = 'NEAREST'
         end if
      else

         sgid= ezqkdef(record%ni, record%nj, record%grtyp, tac%ig1, tac%ig2, tac%ig3, tac%ig4, file%get_unit())
         err = ezgetopt ('USE_1SUBGRID',oldsubgrid_S)
         err = ezgetival('SUBGRIDID',oldsubid)
         if (any(record%grtyp==(/'U','u'/))) then

! We only test samegrid on core domain
            nic = G_ni - Glb_pil_w - Glb_pil_e
            njc = G_nj - Glb_pil_s - Glb_pil_n
            if (record%ni >= nic .and. record%nj/2 >= njc) then
               subid= samegrid_gid ( sgid, Hgc_ig1ro,Hgc_ig2ro              ,&
                                     Hgc_ig3ro,Hgc_ig4ro, F_xfi(1+Glb_pil_w),&
                                     F_yfi(1+Glb_pil_s), nic, njc )
               if (subid >= 0) then
                  write(app_msg,*) 'U grid contains sub grid match to ',ni,nj
                  call lib_log(APP_LIBDYN,APP_INFO,app_msg)
               end if
            end if
            if (subid >= 0) then
               F_inttyp_S = 'NEAREST'
               onesubgrid_S = 'YES'
               err = ezsetival('SUBGRIDID',subid)
               err = ezsetopt ('USE_1SUBGRID',trim(onesubgrid_S))
            end if
         end if

      end if

      success = file%close()

      write (output_unit,1001) F_nomvar_S, trim(F_inttyp_S)

      call record%get_data_array(source)

      err = ezdefset ( dgid, sgid )
      err = ezsetopt ('INTERP_DEGREE', trim(F_inttyp_S))
      err = ezsint   (F_dest, source)
      if (err == 2) then
         write (output_unit,1002) F_nomvar_S, trim(F_inttyp_S)

      end if

      if (subid >= 0) then
!     Reset to original values
         ierr = ezsetival ('SUBGRIDID'   ,oldsubid)
         ierr = ezsetopt  ('USE_1SUBGRID',trim(oldsubgrid_S))
      end if

      if (err < 0) then
         write (output_unit,9005)
         goto 999
      end if

      F_status = 0

 999  call record%free()
      call tic%free()
      call tac%free()
      call query%free()

 1000 format (/'GET_FIELD: File ',a,' open')
 1001 format ( 'GET_FIELD: Interpolating ',a,' with scheme: ',a/)
 1002 format ( 'GET_FIELD: Extrapolating ',a,' with scheme: ',a/)
 9001 format ( 'GET_FIELD: File ',a,' not available'/)
 9002 format ( 'GET_FIELD: Problem opening file: ',a,' with attributes: ',a/)
 9004 format ( 'GET_FIELD: Variable ',a,' not available')
 9005 format ( 'GET_FIELD: Problem with ezscint'/)
!
!-----------------------------------------------------------------------
!
      return
      end
!

