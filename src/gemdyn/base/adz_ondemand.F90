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

      subroutine adz_ondemand ( F_export, F_i0,F_in,F_j0,F_jn,F_k0,F_kn,&
                                F_x, F_y, F_z, Ni, Nj, Nk )
      use adz_mem
      use HORgrid_options
      use glb_pil
      use geomh
      use tdpack
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer, external :: adz_os_send
      integer, intent(IN) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn,Ni, Nj, Nk
      real(kind=REAL64) , dimension(-1:Ni+2,-1:Nj+2,Nk), intent(INOUT) :: F_x,F_y,F_z
      type(ADZ_SLOD), intent(INOUT) :: F_export

      integer i,j,k,cnt,dim,ijk(l_ni*l_nj*l_nk,8)
      integer err, col,ecol, row,erow, pe, np, n3, ind
      integer, dimension (l_ni*l_nj*l_nk) :: pe_target, pe_sort, ijk_target, exp_ijk
      real   , dimension (3*l_ni*l_nj*l_nk) :: exp_target, exp_pos
!
!     ---------------------------------------------------------------
!
return
      dim= l_ni*l_nj*l_nk
      cnt=0
      do k= F_k0, F_kn
         do j= F_j0, F_jn
            do i= F_i0, F_in
               if ( (F_x(i,j,k)<Adz_iminposx) .or. (F_x(i,j,k)>Adz_imaxposx) .or.&
                    (F_y(i,j,k)<Adz_iminposy) .or. (F_y(i,j,k)>Adz_imaxposy) ) then
                    ecol= Ptopo_npex
                    do col=0, Ptopo_npex
                       if (F_x(i,j,k)<Adz_Xlim(1,col)) then
                          ecol= col-1
                          exit
                       endif
                    end do
                    erow= Ptopo_npey
                    do row=0, Ptopo_npey
                       if (F_y(i,j,k)<Adz_Ylim(1,row)) then
                          erow= row-1
                          exit
                       endif
                    end do
                    if (Ptopo_colrow(Ptopo_couleur,ecol,erow)/=Ptopo_myproc) then
                       cnt=cnt+1
                       n3 =(cnt-1)*3+1
                       pe_target (cnt )= Ptopo_colrow(Ptopo_couleur,ecol,erow)
                       ijk_target(cnt )= (k-1)*Adz_2dnh+(j-1)*l_ni+i
                       exp_target(n3  )= F_x(i,j,k)
                       exp_target(n3+1)= F_y(i,j,k)
                       exp_target(n3+2)= F_z(i,j,k)
                    endif
               endif
            end do
         end do
      end do

      call ipsorti(pe_sort,pe_target,cnt)

      F_export%send=0
      pe= -1 ; np= 0 ; err= 0
      do i=1,cnt
         if (pe_target(pe_sort(i)) /= pe) then
            if (np>0) then
               if (adz_os_send(F_export,ijk,exp_pos,exp_ijk,dim,np,pe)<0) err=-1
            end if
            np= 0 ; pe= pe_target(pe_sort(i))
         endif
         np= np+1
         n3 =(np-1)*3+1
         ind = (pe_sort(i)-1)*3+1
         exp_ijk(np)   = ijk_target(pe_sort(i)  )
         exp_pos(n3  ) = exp_target(ind)
         exp_pos(n3+1) = exp_target(ind+1)
         exp_pos(n3+2) = exp_target(ind+2)
      end do

      if (np>0) then
         if (adz_os_send (F_export,ijk,exp_pos,exp_ijk,dim,np,pe,err)<0) err=-1
      endif

      call gem_error(err,'adz_ondemand','PROBLEM with adz_os_send')

      cnt = 1
      do k= 1, 8
         if (F_export%send(k)>0) then
            F_export%ijk(cnt:cnt+F_export%send(k)-1) = ijk(1:F_export%send(k),k)
         endif
         cnt = cnt+F_export%send(k)
      end do

      call adz_export (F_export)

!!$      if (Grd_yinyang_L) then
!!$      cnt= 0
!!$      do k= F_kn/2, F_kn/2
!!$         do j= F_j0, F_jn
!!$            do i= F_i0, F_in
!!$               if ( ((Ptopo_myrow==0).and.(F_y(i,j,k)<Adz_yyminposy)&
!!$                                     .and.(F_x(i,j,k)<Adz_yymaxposx)) .or.&
!!$                    ((Ptopo_mycol==0).and.(F_y(i,j,k)>Adz_yyminposy)&
!!$                                     .and.(F_x(i,j,k)<Adz_yyminposx)) .or.&
!!$                    ((Ptopo_myrow==Ptopo_npey-1).and.(F_y(i,j,k)>Adz_yymaxposy)&
!!$                                                .and.(F_x(i,j,k)>Adz_yyminposx)) .or.&
!!$                    ((Ptopo_mycol==Ptopo_npex-1).and.(F_y(i,j,k)<Adz_yymaxposy)&
!!$                                                .and.(F_x(i,j,k)>Adz_yymaxposx)) ) then
!!$                  ii= F_x(i,j,k)
!!$                  lon= G_xg_8(ii) + ( F_x(i,j,k)-dble(ii) ) * geomh_hx_8
!!$                  jj= F_y(i,j,k)
!!$                  lat= G_yg_8(jj) + ( F_y(i,j,k)-dble(jj) ) * geomh_hy_8
!!$                  call smat ( rot, x, y, lon-pi_8, lat)
!!$                  x= x+pi_8
!!$
!!$                  cnt=cnt+1
!!$               endif
!!$            end do
!!$         end do
!!$      end do
!!$      endif
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_ondemand

      integer function adz_os_send (F_export,F_ijk,F_exp_pos,F_exp_ijk,F_n,F_np,F_pe)
      use adz_mem
      implicit none
      integer, intent(IN) :: F_n,F_np,F_pe
      integer, dimension (F_np), intent(IN)    :: F_exp_ijk
      integer, dimension (F_n,8), intent(OUT)  :: F_ijk
      real   , dimension (3*F_np), intent (IN) :: F_exp_pos
      type(ADZ_SLOD), intent(INOUT) :: F_export

      integer, external :: which_comm
      integer comm_id
!
!     ---------------------------------------------------------------
!
      adz_os_send= -1
      comm_id= which_comm(F_pe)
      if (comm_id>0) then
         F_export%send(comm_id) = F_np
         F_ijk(1:F_np,comm_id) = F_exp_ijk(1:F_np)
         F_export%gpos(1:3*F_np,comm_id) = F_exp_pos
         adz_os_send= 0
      endif
!
!     ---------------------------------------------------------------
!
      return
      end function adz_os_send

integer function which_comm(pe)
use ptopo
implicit none
integer pe
      which_comm=-1
      if ((Ptopo_mycol<Ptopo_npex-1).and.(Ptopo_myrow>0)) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow-1)) which_comm=1
      endif
      if (Ptopo_myrow>0) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol  ,Ptopo_myrow-1)) which_comm=2
      end if
      if ((Ptopo_mycol>0).and.(Ptopo_myrow>0)) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow-1)) which_comm=3
      endif
      if (Ptopo_mycol>0) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow  )) which_comm=4
      endif
      if ((Ptopo_mycol>0).and.(Ptopo_myrow<Ptopo_npey-1)) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol-1,Ptopo_myrow+1)) which_comm=5
      endif
      if (Ptopo_myrow<Ptopo_npey-1) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol  ,Ptopo_myrow+1)) which_comm=6
      end if
      if ((Ptopo_mycol<Ptopo_npex-1).and.(Ptopo_myrow<Ptopo_npey-1)) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow+1)) which_comm=7
      endif
      if (Ptopo_mycol<Ptopo_npex-1) then
         if (pe==Ptopo_colrow(Ptopo_couleur,Ptopo_mycol+1,Ptopo_myrow  )) which_comm=8
      end if

return
end function which_comm
