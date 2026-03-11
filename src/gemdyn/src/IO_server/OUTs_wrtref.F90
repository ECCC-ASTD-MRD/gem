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

      subroutine OUTs_wrtref ( F_stag_S )
      use iso_c_binding
      use, intrinsic :: iso_fortran_env
      use mpi_f08
      use IOs
      use OUTs
      use vGrid_Descriptors
      use rmn_fst24
      implicit none

      character(len=4), intent(IN) :: F_stag_S

#include <rmnlib_basics.hf>
      
      character(len=1) :: familly_uencode_S
      integer :: i0,in,j0,jn,vesion_uencode,niyy,sindx
      integer :: k,err,dimx,dimy,nis,njs,c1,c2
      integer, dimension(:), allocatable :: ip1s
      real(kind=REAL64), dimension(:), allocatable :: zero
      real :: wk
      real, dimension(:), allocatable, target :: yy
      real, dimension(:), pointer :: posx,posy
      logical :: success
      type(fst_record) :: tictac_rec
!     
!--------------------------------------------------------------------
!
      i0 = max( 1   , Out_i0)
      in = min( G_ni, Out_in)
      j0 = max( 1   , Out_j0)
      jn = min( G_nj, Out_jn)
      nis = in - i0 + 1  ;  njs = jn - j0 + 1
      dimx= G_ni+2*G_halox ; dimy= G_nj+2*G_haloy
      
      if (F_stag_S(1:1) =='M') then
         posx => geomh_longs
         posy => geomh_latgs
         Out_ig3  = 1
      end if
      if (F_stag_S(1:1) =='U') then
         posx => geomh_longu
         posy => geomh_latgs
         Out_ig3  = 2
         in = min( G_ni-1, Out_in)
      end if
      if (F_stag_S(1:1) =='V') then
         posx => geomh_longs
         posy => geomh_latgv
         Out_ig3  = 3
         jn = min( G_nj-1, Out_jn)
      end if
      if (F_stag_S(1:1) =='F') then
         posx => geomh_longu
         posy => geomh_latgv
         Out_ig3  = 4
         in = min( G_ni-1, Out_in)
         jn = min( G_nj-1, Out_jn)
      end if
      Out_ig4 = 0

      call OUTs_igs ( Out_ig1, Out_ig2, posx(1), posy(1), G_ni,G_nj,&
                      Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4           ,&
                      i0,in,1, j0,jn,1 )
      
      tictac_rec%typvar='X'
      tictac_rec%etiket=Out_etik_S
      tictac_rec%dateo=Out_dateo
      tictac_rec%deet=0
      tictac_rec%npas=0
      tictac_rec%data_type = FST_TYPE_REAL_IEEE
      tictac_rec%data_bits=32!Out_nbit
      tictac_rec%pack_bits=32
      tictac_rec%nk=1
      tictac_rec%ip1=Out_ig1
      tictac_rec%ip2=Out_ig2
      tictac_rec%ip3=Out_ig3

      if (OUTs_1o1_L) then
         if ( (Grd_yinyang_L) .and. (.not.Out_reduc_l) ) then
            vesion_uencode    = 1
            familly_uencode_S = 'F'

            niyy=5+2*(10+nis+njs)
            allocate (yy(niyy))

            yy(1 ) = iachar(familly_uencode_S)
            yy(2 ) = vesion_uencode
            yy(3 ) = 2          ! 2 grids (Yin & Yang)
            yy(4 ) = 1          ! the 2 grids have same resolution
            yy(5 ) = 1          ! the 2 grids have same area extension
!YIN
            sindx  = 6
            yy(sindx  ) = nis
            yy(sindx+1) = njs
            yy(sindx+2) = posx(i0)
            yy(sindx+3) = posx(i0+nis-1)
            yy(sindx+4) = posy(j0)
            yy(sindx+5) = posy(j0+njs-1)
            yy(sindx+6) = Grd_xlat1
            yy(sindx+7) = Grd_xlon1
            yy(sindx+8) = Grd_xlat2
            yy(sindx+9) = Grd_xlon2
            yy(sindx+10    :sindx+9+nis    )= &
            posx(i0:i0+nis-1)
            yy(sindx+10+nis:sindx+9+nis+njs)= &
            posy(j0:j0+njs-1)
!YAN
            sindx  = sindx+10+nis+njs
            yy(sindx  ) = nis
            yy(sindx+1) = njs
            yy(sindx+2) = posx(i0)
            yy(sindx+3) = posx(i0+nis-1)
            yy(sindx+4) = posy(j0)
            yy(sindx+5) = posy(j0+njs-1)
            yy(sindx+6) = Grd_xlat1Y
            yy(sindx+7) = Grd_xlon1Y
            yy(sindx+8) = Grd_xlat2Y
            yy(sindx+9) = Grd_xlon2Y
            yy(sindx+10    :sindx+9+nis    )= &
            posx(i0:i0+nis-1)
            yy(sindx+10+nis:sindx+9+nis+njs)= &
            posy(j0:j0+njs-1)

            tictac_rec%nomvar='^>'
            tictac_rec%grtyp=familly_uencode_S
            tictac_rec%ni=niyy
            tictac_rec%nj=1
            tictac_rec%ig1=vesion_uencode
            tictac_rec%ig2=0
            tictac_rec%ig3=0
            tictac_rec%ig4=0

            tictac_rec%data=c_loc(yy)
            success = Out_file%write(tictac_rec,rewrite=FST_SKIP)

 		!field, work, npak, iun, dateo, deet, npas, ni,
                !         nj,nk, ip1, ip2, ip3, typvar, nomvar, etiket, grtyp,
                !         ig1, ig2, ig3, ig4, datyp, rewrit)
            !err= fstecr(yy,yy, -32, Out_unf,Out_dateo,0,0,niyy,1,1  ,&
            !            Out_ig1,Out_ig2,Out_ig3,'X','^>',Out_etik_S ,&
            !            familly_uencode_S,vesion_uencode,0,0,0      ,&
            !            5, .true.)
            deallocate (yy, STAT = err)
                         
         else
            tictac_rec%nomvar='>>'
            tictac_rec%grtyp='E'
            tictac_rec%ni=nis
            tictac_rec%nj=1
            tictac_rec%ig1=Rot_ig1
            tictac_rec%ig2=Rot_ig2
            tictac_rec%ig3=Rot_ig3
            tictac_rec%ig4=Rot_ig4
            tictac_rec%data=c_loc(posx(i0))
            success = Out_file%write(tictac_rec,rewrite=FST_SKIP)

         !err=fstecr(posx(i0),wk,-32,Out_unf,Out_dateo   ,&
         !           0,0, nis,1,1, Out_ig1,Out_ig2,Out_ig3,'X', '>>'    ,&
         !           Out_etik_S,'E',&
         !           Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, 5, .true.)

            tictac_rec%nomvar='^^'
            tictac_rec%ni=1
            tictac_rec%nj=njs
            tictac_rec%data=c_loc(posy(j0))
            success = Out_file%write(tictac_rec,rewrite=FST_SKIP)

         !err=fstecr(posy(j0),wk,-32,Out_unf,Out_dateo   ,&
         !           0,0, 1,njs,1,Out_ig1,Out_ig2,Out_ig3,'X', '^^'    ,&
         !           Out_etik_S,'E',&
         !           Rot_ig1, Rot_ig2, Rot_ig3, Rot_ig4, 5, .true.)
         endif
         
         if ( F_stag_S(2:2) == 'S' ) return
!if (( F_stag_S(2:2) == 'M' ) .or. ( F_stag_S(2:2) == 'T' )) then
         call Outs_wrtvref(F_stag_S,Out_ig1,Out_ig2)
      endif
!     
!--------------------------------------------------------------------
!
      return
      end subroutine OUTs_wrtref
