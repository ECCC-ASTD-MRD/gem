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
      implicit none
#include <arch_specific.hf>

      character(len=*) varname, filename, grdtype
      integer nx,ny,nz,ig1,ig2,ig3
      real buf(nx,ny,nz),maxrange

      integer  fnom,fstinl,fstecr,fclos,longueur
      external fnom,fstinl,fstecr,fclos,longueur

      character(len=7) :: startindx
      character(len=1024) :: fn
      integer unf,err,n1,n2,n3,nlis,list_x,list_y, &
              nipos,njpos,pnip1,pnip3,i,j,k
      parameter (nlis = 1024)
      integer liste (nlis)
      real wk3(nx,ny)
!
!-------------------------------------------------------------------
!
      write (startindx,'((i3.3),a1,(i3.3))') Ptopo_mycol,'-',Ptopo_mycol
      fn = Path_input_S(1:longueur(Path_input_S))// &
           filename(1:longueur(filename))//'_'//startindx
!
      unf = 0
      err = fnom  (unf, fn, 'rnd', 0)
      call fstouv (unf,'rnd')
!
      if (grdtype(1:1) == "Z") then
!
      err = fstinl (unf, nipos,n2,n3, -1, ' ', ig1, ig2,ig3,' ', &
                                      '>>', liste, list_x, nlis)
      err = fstinl (unf, n1,njpos,n3, -1, ' ', ig1, ig2,ig3,' ', &
                                      '^^', liste, list_y, nlis)
!
      if ((list_x < 1).or.(list_y < 1).or. &
          (nipos /= nx).or.(njpos /= ny)) then
         err= fstecr(Geomh_longs(Ptopo_gindx(1,Ptopo_myproc+1)),wk3, &
                     -32,unf,0,0,0,nx,1,1,ig1,ig2,ig3,'X', '>>', &
                     Out3_etik_S,Hgc_gxtyp_s,Hgc_ig1ro,Hgc_ig2ro, &
                     Hgc_ig3ro,Hgc_ig4ro, 5, .true.)
         err= fstecr(Geomh_latgs(Ptopo_gindx(3,Ptopo_myproc+1)),wk3, &
                     -32,unf,0,0,0,1,ny,1,ig1,ig2,ig3,'X', '^^', &
                     Out3_etik_S,Hgc_gxtyp_s,Hgc_ig1ro,Hgc_ig2ro, &
                     Hgc_ig3ro,Hgc_ig4ro, 5, .true.)
      end if
!
      end if
!
      do k = 1, nz
!
      do j = 1, ny
      do i = 1, nx
        wk3(i,j)=buf(i,j,k)
        if (wk3(i,j) > maxrange) wk3(i,j)=-999999.
        if (wk3(i,j) < -maxrange) wk3(i,j)=-999999.
      end do
      end do
!
      pnip1 = k
      pnip3 = Out_ip3
      pnip3 = Lctl_step
      if (pnip3 < 0) pnip3 = Out_npas

      err = fstecr (wk3,wk3,-32,unf,Out_dateo,int(Out_deet),Out_npas, &
                    nx,ny,1,pnip1,Out_ip2,pnip3,'P',varname,Out3_etik_S, &
                    grdtype,ig1,ig2,ig3,0,5,.false.)
!
      end do
!
      call  fstfrm (unf)
      err = fclos  (unf)
!
!-------------------------------------------------------------------
!
      return
      end
