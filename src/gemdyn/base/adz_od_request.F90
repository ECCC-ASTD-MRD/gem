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

      subroutine adz_od_request ( F_export, F_i0,F_in,F_j0,F_jn,&
                                  F_k0,F_kn, F_x,F_y,F_z, Ni,Nj,Nk )
      use adz_mem
      use HORgrid_options
      use glb_pil
      use geomh
      use tdpack
      use ptopo
      implicit none
#include <arch_specific.hf>

      integer, external :: adz_stkr, find_col,find_row, findpeyy
      integer, intent(IN) :: F_i0,F_in,F_j0,F_jn,F_k0,F_kn,Ni, Nj, Nk
      real(kind=REAL64) , dimension(-1:Ni+2,-1:Nj+2,Nk), intent(INOUT) :: F_x,F_y,F_z
      type(ADZ_SLOD), intent(INOUT) :: F_export

      integer i,j,k,cnt
      integer err, errF, ecol, erow, pe, np, n3, ind
      integer, dimension (l_ni*l_nj*l_nk) :: pe_target, pe_sort, ijk_target, exp_ijk
      real   , dimension (3*l_ni*l_nj*l_nk) :: exp_target, exp_pos
      real(kind=REAL64) :: x_a,y_a
!
!     ---------------------------------------------------------------
!
      if (.not.ADZ_OD_L) return
      
! Determine points that are outside my coverage and store them in exp_target
      cnt=0
      do k= F_k0, F_kn
         do j= F_j0, F_jn
            do i= F_i0, F_in
               
! given that the upstream positions are fractional global indices
! this next section should be replaced by a rpn_comm function
! that would return (ecol,erow) of a particular target F_x(i,j,k),F_y(i,j,k)
! If the target position is outside the LAM domain, that function
! will return (-1,-1). In a regular LAM this will mean clipping of
! the trajectories. In a GY, the (-1,-1) would trigger a rotation
! of the F_x(i,j,k),F_y(i,j,k) positions to the rotation of the other colour
! and a subsequent search of (ecol,erow) in that other colour.

               ecol= Ptopo_mycol ; erow= Ptopo_myrow
               if ((F_x(i,j,k)<Adz_iminposx).or.(F_x(i,j,k)>Adz_imaxposx)) ecol= find_col(F_x(i,j,k))
               if ((F_y(i,j,k)<Adz_iminposy).or.(F_y(i,j,k)>Adz_imaxposy)) erow= find_row(F_y(i,j,k))
               if ((Grd_yinyang_L).and.((ecol<0).or.(erow<0).or.&
                   (ecol>Ptopo_npex-1).or.(erow>Ptopo_npey-1))) then
                  pe = findpeyy(F_x(i,j,k),F_y(i,j,k),x_a,y_a)
                  cnt= cnt+1
                  n3 = (cnt-1)*3+1
                  pe_target (cnt )= pe
                  ijk_target(cnt )= (k-1)*Adz_2dnh+(j-1)*l_ni+i
                  exp_target(n3  )= x_a
                  exp_target(n3+1)= y_a
                  exp_target(n3+2)= F_z(i,j,k)
               else
                  pe= Ptopo_colrow(Ptopo_couleur,ecol,erow)
                  if (pe/=Ptopo_myproc) then
                     cnt= cnt+1
                     n3 = (cnt-1)*3+1
                     pe_target (cnt )= pe + Ptopo_couleur*Ptopo_numproc
                     ijk_target(cnt )= (k-1)*Adz_2dnh+(j-1)*l_ni+i
                     exp_target(n3  )= F_x(i,j,k)
                     exp_target(n3+1)= F_y(i,j,k)
                     exp_target(n3+2)= F_z(i,j,k)
                  endif
               endif
            end do
         end do
      end do

      call ipsorti (pe_sort, pe_target, cnt)
      
      pe= -1 ; np= 0 ; F_export%npe= 0 ; err= 0
      do i=1,cnt
         if (pe_target(pe_sort(i)) /= pe) then
            if (np>0) then
               err= min(err, adz_stkr (F_export,exp_pos,exp_ijk,np,pe))
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
         err= min(err, adz_stkr (F_export,exp_pos,exp_ijk,np,pe))
      endif
      
      F_export%from= -999
      ORIGIN_COUNT = 2
      TARGET_COUNT = ORIGIN_COUNT
      ORIGIN_DATATYPE = INTEGER_DATATYPE
      TARGET_DATATYPE = INTEGER_DATATYPE
      
      call MPI_Win_fence(0, F_export%winreqs, errF)
      call gem_error(min(err,errF),'adz_od_request','NOT ENNOUGH MEMORY - Adz_MAX_MPI_OS_SIZE')
 
! Deposit requests to proper target PEs (one-sided)   
      do i=1,F_export%npe
         TARGET_RANK = F_export%stk(1,i)
         TARGET_DISP = (Ptopo_world_myproc*2)
         call MPI_Put( F_export%stk(2,i), ORIGIN_COUNT, ORIGIN_DATATYPE, TARGET_RANK,&
                       TARGET_DISP, TARGET_COUNT, TARGET_DATATYPE     ,&
                       F_export%winreqs, err )
      end do
!
!     ---------------------------------------------------------------
!
      return
      end subroutine adz_od_request

      integer function adz_stkr ( F_export,F_exp_pos,F_exp_ijk,F_np,F_pe )
      use adz_mem
      implicit none
      integer, intent(IN) :: F_np,F_pe
      integer, dimension (F_np)  , intent(IN)  :: F_exp_ijk
      real   , dimension (3*F_np), intent(IN)  :: F_exp_pos
      type(ADZ_SLOD), intent(INOUT) :: F_export

      integer offset
!
!     ---------------------------------------------------------------
!
      adz_stkr= 0
      F_export%npe= F_export%npe+1
      if (F_export%npe==1) then
         offset= 0
      else
         offset= F_export%stk(4,F_export%npe-1)
      endif
      F_export%stk(1,F_export%npe)= F_pe
      F_export%stk(2,F_export%npe)= F_np
      F_export%stk(3,F_export%npe)= offset
      F_export%stk(4,F_export%npe)= offset+F_np*3 !3 words fx,fy,fz
      if ( (F_np>Adz_MAX_MPI_OS_SIZE) .or. &
           (offset+3*F_np>size(F_export%requests))) then
         adz_stkr= -1
      else
         F_export%dest(1:F_np,F_pe)  = F_exp_ijk(1:  F_np)
! F_export%requests is exposed to all PEs in F_export%wintraj window
         F_export%requests(offset+1:offset+3*F_np) = F_exp_pos(1:3*F_np)
      endif
!
!     ---------------------------------------------------------------
!
      return
      end function adz_stkr
