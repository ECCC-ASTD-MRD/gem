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

!*s/r out_href - write horizontal positional parameters

      subroutine out_href (F_arakawa_S ,F_x0, F_x1, F_stridex,&
                                        F_y0, F_y1, F_stridey )
      use geomh
      use glb_ld
      use hgc
      use HORgrid_options
      use lun
      use out_mod
      use out3
      use ptopo
      implicit none
#include <arch_specific.hf>

      character(len=*), intent(in ) :: F_arakawa_S
      integer,        intent(in ) :: F_x0,F_x1,F_stridex,&
                                     F_y0,F_y1,F_stridey
#include <rmnlib_basics.hf>

!!$      integer, external :: out_samegrd
      character(len=1) :: familly_uencode_S
      integer :: err,ni,nis,njs,niyy,ix1,ix2,ix3,ix4, &
               sindx,i0,in,j0,jn,vesion_uencode
      real :: wk
      real, dimension(:), pointer     :: posx,posy
      real, dimension(:), allocatable :: yy

      character(len=1) :: nomvar
      integer :: nbit, knd, lstep
      integer, dimension(1) :: ind_o
      real, dimension(1,1,1) :: fa
      real, dimension(1) :: rf

!
!----------------------------------------------------------------------
!
      call out_fstecr ( fa,1,1,1,1,rf,nomvar,wk,wk,knd,lstep,1,ind_o,1,nbit,.true. )

! to be completed if at all usefull
!      old_grid_L= out_samegrd ( F_arakawa_S ,F_x0, F_x1, F_stridex,&
!                                F_y0, F_y1, F_stridey, .false. )

      i0 = max( 1   , F_x0)
      in = min( G_ni, F_x1)
      j0 = max( 1   , F_y0)
      jn = min( G_nj, F_y1)

      if (F_arakawa_S =='Mass_point') then
         posx => geomh_longs
         posy => geomh_latgs
         Out_ig3  = 1
      end if
      if (F_arakawa_S =='U_point') then
         posx => geomh_longu
         posy => geomh_latgs
         Out_ig3  = 2
         in = min( G_ni-1, F_x1)
      end if
      if (F_arakawa_S =='V_point') then
         posx => geomh_longs
         posy => geomh_latgv
         Out_ig3  = 3
         jn = min( G_nj-1, F_y1)
      end if
      if (F_arakawa_S =='F_point') then
         posx => geomh_longu
         posy => geomh_latgv
         Out_ig3  = 4
         in = min( G_ni-1, F_x1)
         jn = min( G_nj-1, F_y1)
      end if

      Out_ig4 = 0

      call set_igs2 ( Out_ig1, Out_ig2, posx(1), posy(1), G_ni,G_nj,&
                      Hgc_ig1ro,Hgc_ig2ro,Hgc_ig3ro,Hgc_ig4ro,&
                      i0,in,F_stridex, j0,jn,F_stridey )

      Out_stride = 1 ! can only be one for now
      Out_gridi0 = i0
      Out_gridin = in
      Out_gridj0 = j0
      Out_gridjn = jn

      nis = in - i0 + 1  ;  njs = jn - j0 + 1
      if ( (nis <= 0) .or. (njs <= 0) ) then
         if (Lun_out > 0) write(Lun_out,9000)
         return
      end if

      nis = (in - i0) / Out_stride + 1
      njs = (jn - j0) / Out_stride + 1

      if ((Out3_iome >= 0) .and. (Ptopo_couleur == 0)) then

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
            yy(sindx+2) = posx(Out_gridi0)
            yy(sindx+3) = posx(Out_gridi0+nis-1)
            yy(sindx+4) = posy(Out_gridj0)
            yy(sindx+5) = posy(Out_gridj0+njs-1)
            yy(sindx+6) = Out_rot(1)
            yy(sindx+7) = Out_rot(2)
            yy(sindx+8) = Out_rot(3)
            yy(sindx+9) = Out_rot(4)
            yy(sindx+10    :sindx+9+nis    )= &
            posx(Out_gridi0:Out_gridi0+nis-1)
            yy(sindx+10+nis:sindx+9+nis+njs)= &
            posy(Out_gridj0:Out_gridj0+njs-1)

!YAN
            sindx  = sindx+10+nis+njs
            yy(sindx  ) = nis
            yy(sindx+1) = njs
            yy(sindx+2) = posx(Out_gridi0)
            yy(sindx+3) = posx(Out_gridi0+nis-1)
            yy(sindx+4) = posy(Out_gridj0)
            yy(sindx+5) = posy(Out_gridj0+njs-1)
            yy(sindx+6) = Out_rot(5)
            yy(sindx+7) = Out_rot(6)
            yy(sindx+8) = Out_rot(7)
            yy(sindx+9) = Out_rot(8)
            yy(sindx+10    :sindx+9+nis    )= &
            posx(Out_gridi0:Out_gridi0+nis-1)
            yy(sindx+10+nis:sindx+9+nis+njs)= &
            posy(Out_gridj0:Out_gridj0+njs-1)

            err= fstecr(yy,yy, -32, Out_unf,Out_dateo,0,0,niyy,1,1  ,&
                        Out_ig1,Out_ig2,Out_ig3,'X','^>',Out_etik_S ,&
                        familly_uencode_S,vesion_uencode,0,0,0      ,&
                        5, .true.)
            deallocate (yy, STAT = err)

         else

            if ( Out_stride <= 1 ) then
               ni = nis
               if ( (Grd_typ_S(1:2) == 'GU') .and. (nis == G_ni) ) &
               ni = G_ni+1
               ix1= Ptopo_couleur*4+1
               ix2= Ptopo_couleur*4+2
               ix3= Ptopo_couleur*4+3
               ix4= Ptopo_couleur*4+4
               err=fstecr(posx(Out_gridi0),wk,-32,Out_unf,Out_dateo   ,&
                    0,0, ni,1,1, Out_ig1,Out_ig2,Out_ig3,'X', '>>'    ,&
                    Out_etik_S,Out_gridtyp_S,Out_ixg(ix1),Out_ixg(ix2),&
                    Out_ixg(ix3), Out_ixg(ix4), 5, .true.)
               err=fstecr(posy(Out_gridj0),wk,-32,Out_unf,Out_dateo   ,&
                    0,0, 1,njs,1,Out_ig1,Out_ig2,Out_ig3,'X', '^^'    ,&
                    Out_etik_S,Out_gridtyp_S,Out_ixg(ix1),Out_ixg(ix2),&
                    Out_ixg(ix3), Out_ixg(ix4), 5, .true.)
            end if

         end if
      end if

 9000 format(/,'OUT_HREF - no grid to output'/)
!
!----------------------------------------------------------------------
!
      return
      end
