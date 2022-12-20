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

!**s/r itf_phy_glbstat - Global statistics on physical (PW) quantities
!
      subroutine itf_phy_glbstat (name_S)
      use step_options
      use lun
      use HORgrid_options
      use glb_ld
      use tr3d
      use phy_itf, only: phy_get,phymeta,phy_getmeta
      use rmn_gmm
      implicit none

      character(len=*) name_S
#include <arch_specific.hf>

      type(phymeta),pointer :: pmeta(:)

      integer :: err,nvars,ivar,ii
      integer gphy_i0, gphy_in, gphy_j0, gphy_jn,offset,cnt
      integer p_li0,p_li1,p_lj0,p_lj1
      character(len=1) :: bus_S(4)
      character(len=4) :: busname_S(4)
      logical flag
      real, dimension(:,:,:), pointer :: ptr3d
      real, dimension(:,:,:), allocatable, target:: data3d
!     ________________________________________________________________
!
      flag = .false.
      if (Step_gstat > 0) flag = (mod(Lctl_step,Step_gstat) == 0)
      if (.not. flag) return
      offset = 2
      gphy_i0 = 1 + offset
      gphy_in = G_ni - offset
      gphy_j0 = 1 + offset
      gphy_jn = G_nj - offset
      p_li0= Grd_lphy_i0 ; p_li1=Grd_lphy_in
      p_lj0= Grd_lphy_j0 ; p_lj1=Grd_lphy_jn
      !  F_npath='O',F_bpath='PVED', F_quiet=.true.)
      bus_S(1)='E'; busname_S(1)='ENTR'
      bus_S(2)='D'; busname_S(2)='DYNA'
      bus_S(3)='P'; busname_S(3)='PERM'
      bus_S(4)='V'; busname_S(4)='VOLA'

    do ii = 1,4
      nullify(ptr3d,pmeta)
      nvars =phy_getmeta (pmeta, ' ', &
                          F_npath='O',F_bpath=bus_S(ii), F_quiet=.true.)
      cnt=0
      do ivar = 1,nvars
         cnt= max(cnt,pmeta(ivar)%fmul * pmeta(ivar)%nk)
      end do

      allocate ( data3d(l_ni,l_nj,cnt))

      if (Lun_out > 0) write(Lun_out,*)'PHYBLOC STAT on BUS', busname_S(ii)
       do ivar = 1,nvars
          data3d=0.
          cnt= pmeta(ivar)%fmul * pmeta(ivar)%nk
          ptr3d => data3d(p_li0:p_li1,p_lj0:p_lj1,1:cnt)
          err =  phy_get ( ptr3d, trim(pmeta(ivar)%oname), &
                                  F_npath='O', F_bpath='PVED')
         if (err == 0) then
          call pglbstat ( data3d,trim(pmeta(ivar)%vname),trim(name_S),&
                                  1,l_ni,1,l_nj, &
                       1,cnt, gphy_i0,gphy_in,gphy_j0,gphy_jn,1,cnt )
         else
          if (Lun_out > 0) write (Lun_out,*) &
                    'error in getting ',trim(pmeta(ivar)%vname)
         end if
       end do
      deallocate(data3d)
    end do

!For indepth glbstats on particular variables, and layer by layer
!      do ivar = 1,nvars
!       if (Lun_out > 0) print *,'ivar=',ivar,'var=',pmeta(ivar)%vname
!       if ((trim(pmeta(ivar)%vname)=="tr/hu:p").or. &
!           (trim(pmeta(ivar)%vname)=="tr/qc:p").or. &
!           (trim(pmeta(ivar)%vname)=="tr/tt:p") ) then
!       if ((trim(pmeta(ivar)%vname)=="tr/hu:p") ) then
!            data3d=0.
!            cnt= pmeta(ivar)%fmul * pmeta(ivar)%nk
!            ptr3d => data3d(p_li0:p_li1,p_lj0:p_lj1,1:cnt)
!            err =  phy_get ( ptr3d, trim(pmeta(ivar)%oname), &
!                                 F_npath='O', F_bpath='PVED')
!            if (err == 0) then
!               do i=1,cnt
!                call pglbstat ( data3d,trim(pmeta(ivar)%vname),trim(name_S),&
!                                1,l_ni,1,l_nj,1,cnt,&
!             !        gphy_i0,gphy_in,gphy_j0,gphy_jn,1,cnt )
!                      gphy_i0,gphy_in,gphy_j0,gphy_jn,i,i )
!               end do
!            else
!        if (Lun_out > 0) write (Lun_out,*) 'error in getting ',trim(pmeta(ivar)%vname)
!            end if
!       end if
!      end do


!     ________________________________________________________________
!
      return
      end
