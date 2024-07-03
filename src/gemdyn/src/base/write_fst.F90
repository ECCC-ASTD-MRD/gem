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

!**s/r write_fst

      SUBROUTINE write_fst (buf,nx,ny,nz,varname,maxrange, &
                              ig1,ig2,ig3,grdtype,filename)
      use step_options
      use out_options
      use geomh
      use glb_ld
      use out_mod
      use out3
      use hgc
      use path
      use ptopo
      use rmn_fst24

      implicit none
#include <arch_specific.hf>

      character(len=*) varname, filename, grdtype
      integer nx,ny,nz,ig1,ig2,ig3
      real buf(nx,ny,nz),maxrange

      type(fst_file)   :: file
      type(fst_record) :: tic,tac,rec
      type(fst_query)  :: query

      character(len=7) :: startindx
      character(len=1024) :: fn
      integer unf,pnip1,pnip3,i,j,k
      real, dimension(nx,ny), target :: wk3
!
!-------------------------------------------------------------------
!
      write (startindx,'((i3.3),a1,(i3.3))') Ptopo_mycol,'-',Ptopo_mycol
      fn = Path_input_S(1:len_trim(Path_input_S))// &
           filename(1:len_trim(filename))//'_'//startindx

      unf = 0
      success = file%open(fn,'RND+R/W')
 
      if (grdtype(1:1) == "Z") then

         query = file%new_query(ig1=ig1,ig2=ig2,ig3=ig3,nomvar=">>  ")
         success = query%find_next(tic)
         if ((.not. success) .or. (nx /= tic%ni)) then
           tic%data=c_loc(Geomh_longs(Ptopo_gindx(1,Ptopo_myproc+1)))
           tic%data_type = FST_TYPE_REAL_IEEE
           tic%data_bits=32
           tic%pack_bits=32
           tic%dateo=0
           tic%deet=0
           tic%npas=0
           tic%ni=nx
           tic%nj=1
           tic%nk=1
           tic%ip1=ig1
           tic%ip2=ig2
           tic%ip3=ig3
           tic%nomvar='>>'
           tic%typvar='X'
           tic%etiket=Out3_etik_S
           tic%grtyp=Hgc_gxtyp_s
           tic%ig1=Hgc_ig1ro
           tic%ig2=Hgc_ig2ro
           tic%ig4=Hgc_ig3ro
           tic%ig1=Hgc_ig4ro

           success=file%write(tic,rewrite=FST_SKIP);
         endif
         call query%free()

         query = file%new_query(ig1=ig1,ig2=ig2,ig3=ig3,nomvar="^^  ")
         success = query%find_next(tac)
         if ((.not. success) .or. (ny /= tac%nj)) then
           tac%data=c_loc(Geomh_latgs(Ptopo_gindx(3,Ptopo_myproc+1)))
           tic%data_type = FST_TYPE_REAL_IEEE
           tic%data_bits=32
           tic%pack_bits=32
           tac%dateo=0
           tac%deet=0
           tac%npas=0
           tac%ni=ny
           tac%nj=1
           tac%nk=1
           tac%ip1=ig1
           tac%ip2=ig2
           tac%ip3=ig3
           tac%nomvar='^^'
           tac%typvar='X'
           tac%etiket=Out3_etik_S
           tac%grtyp=Hgc_gxtyp_s
           tac%ig1=Hgc_ig1ro
           tac%ig2=Hgc_ig2ro
           tac%ig4=Hgc_ig3ro
           tac%ig1=Hgc_ig4ro

           success=file%write(tac,rewrite=FST_SKIP);
         endif
         call query%free()
      end if

      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               wk3(i,j)=buf(i,j,k)
               if (wk3(i,j) > maxrange) wk3(i,j)=-999999.
               if (wk3(i,j) < -maxrange) wk3(i,j)=-999999.
            end do
         end do

         pnip1 = k
         pnip3 = Out_ip3
         pnip3 = Lctl_step

         if (pnip3 < 0) pnip3 = Out_npas

           rec%data=c_loc(wk3)
           rec%pack_bits=32
           rec%dateo=Out_dateo
           rec%deet=int(Out_deet)
           rec%npas=Out_npas
           rec%ni=nx
           rec%nj=ny
           rec%nk=1
           rec%ip1=pnip1
           rec%ip2=Out_ip2
           rec%ip3=pnip3
           rec%nomvar=varname
           rec%typvar='P'
           rec%etiket=Out3_etik_S
           rec%grtyp=grdtype
           rec%ig1=ig1
           rec%ig2=ig2
           rec%ig4=ig3
           rec%ig1=0
           rec%data_type = FST_TYPE_REAL_IEEE

           success=file%write(rec,rewrite=FST_SKIP);
      end do

      success=file%close()
!
!-------------------------------------------------------------------
!
      return
      end
