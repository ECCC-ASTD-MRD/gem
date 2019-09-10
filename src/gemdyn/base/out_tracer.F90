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

!**s/r out_tracer - output tracer

      subroutine out_tracer (levset, set)
      use out_options
      use gmm_itf_mod
      use gmm_pw
      use glb_ld
      use HORgrid_options
      use levels
      use outd
      use outp
      use out3
      use phy_itf, only: phy_get
      use tr3d
      use vertical_interpolation, only: vertint2
      use vGrid_Descriptors, only: vgrid_descriptor, vgd_get, vgd_free, VGD_OK, VGD_ERROR
      use vgrid_wb, only: vgrid_wb_get

      implicit none
#include <arch_specific.hf>

      integer, intent(in) :: levset, set

      type(vgrid_descriptor) :: vcoord
      character(len=512) :: fullname
      integer i,j,k,ii,n,nko,model_nk,knd,istat,indxtr
      integer, dimension(:), allocatable::indo
      integer, dimension(:), pointer :: ip1t
      real ,dimension(:), allocatable::prprlvl,rf
      real, dimension(:), pointer :: hybt
      save hybt
      real,dimension(l_minx:l_maxx,l_miny:l_maxy,G_nk+1), target :: w4
      real,dimension(:,:,:), allocatable :: cible
      real, dimension(:,:  ), pointer :: qdiag
      real, dimension(:,:,:), pointer :: tr1,tr5,ptr3d
      logical :: write_diag_lev,near_sfc_L,outvar_L
      real hybt_gnk2(1)
      integer ind0(1)
      integer, dimension(3), save :: lijk = [ -1,-1,-1 ]
      integer, dimension(3), save :: uijk = [ -1,-1,-1 ]
!
!----------------------------------------------------------------------
!
      ptr3d => w4(Grd_lphy_i0:Grd_lphy_in,Grd_lphy_j0:Grd_lphy_jn,1:G_nk+1)

      if (Level_typ_S(levset) == 'M') then  ! output tracers on model levels

         knd=Level_kind_ip1
!        Setup the indexing for output
         allocate (indo( min(Level_max(levset),Level_thermo) ))
         call out_slev2(Level(1,levset), Level_max(levset), &
                       Level_thermo,indo,nko,near_sfc_L)
         write_diag_lev= near_sfc_L .and. out3_sfcdiag_L

!        Retreieve vertical coordinate description
         if ( .not. associated(hybt) ) then
            nullify(ip1t,hybt)
            istat = vgrid_wb_get('ref-t',vcoord,ip1t)
            deallocate(ip1t); nullify(ip1t)
            if (vgd_get(vcoord,'VCDT - vertical coordinate (t)',hybt) /= VGD_OK) istat = VGD_ERROR
         end if
         istat = vgd_free(vcoord)
         hybt_gnk2(1)=hybt(G_nk+2)
         ind0(1)=1

         do ii=1,Outd_var_max(set)
            outvar_L=.false.
            do n=1,Tr3d_ntr
               if (Outd_var_S(ii,set) == trim(Tr3d_name_S(n))) then
                  nullify (tr1)
                  fullname= 'TR/'//trim(Tr3d_name_S(n))//':P'
                  indxtr=n
                  istat = gmm_get(fullname,tr1)
                  if (.not. GMM_IS_ERROR(istat)) then
                     outvar_L=.true.
                  end if
                  exit
               end if
            end do

            if (outvar_L) then
               w4(:,:,1:G_nk) = tr1(:,:,1:G_nk)
               model_nk = G_nk
               if (write_diag_lev) then
                  if (trim(Tr3d_name_S(indxtr))=='HU') then
                     istat = gmm_get(gmmk_diag_hu_s,qdiag)
                     if (istat == 0) w4(:,:,G_nk+1) = qdiag(:,:)
                  else
                     w4(:,:,G_nk+1) = tr1(:,:,G_nk)
                     istat = phy_get (ptr3d, trim(fullname), F_npath='VO', F_bpath='D',&
                                      F_start=lijk, F_end=uijk, F_quiet=.true.)
                  end if
               else
                  istat=-1
               end if
               if (istat == 0) model_nk = G_nk + 1

               if (Out3_cliph_L) then
                  do k=1,model_nk
                  do j=1,l_nj
                  do i=1,l_ni
                     w4(i,j,k) = max ( w4(i,j,k), 0. )
                  end do
                  end do
                  end do
               end if
               call out_fstecr3(w4,l_minx,l_maxx,l_miny,l_maxy,hybt,&
                       Outd_var_S(ii,set),Outd_convmult(ii,set)    ,&
                       Outd_convadd(ii,set),knd,-1,G_nk,indo,nko   ,&
                       Outd_nbit(ii,set),.false. )
               if (model_nk > G_nk) &
                  call out_fstecr3 ( w4(l_minx,l_miny,G_nk+1)       ,&
                            l_minx,l_maxx,l_miny,l_maxy,hybt_gnk2   ,&
                            Outd_var_S(ii,set),Outd_convmult(ii,set),&
                            Outd_convadd(ii,set),Level_kind_diag    ,&
                            -1,1,ind0,1,Outd_nbit(ii,set),.false. )

            end if
         end do
         deallocate(indo)

      else  ! output tracers on pressure levels

         knd=2

!        Setup the indexing for output
         nko=Level_max(levset)
         allocate ( indo(nko), rf(nko) , prprlvl(nko) , &
                  tr5(l_minx:l_maxx,l_miny:l_maxy,nko), &
               cible(l_minx:l_maxx,l_miny:l_maxy,nko))

         istat= gmm_get(gmmk_pw_log_pt_s  , pw_log_pt)

         do i = 1, nko
            indo(i)=i
            rf(i)= Level(i,levset)
            prprlvl(i) = rf(i) * 100.0
            cible(:,:,i) = log(prprlvl(i))
         end do

         do ii=1,Outd_var_max(set)

            outvar_L=.false.
            do n=1,Tr3d_ntr
               if (Outd_var_S(ii,set) == trim(Tr3d_name_S(n))) then
                  nullify (tr1)
                  fullname= 'TR/'//trim(Tr3d_name_S(n))//':P'
                  indxtr=n
                  istat = gmm_get(fullname,tr1)
                  if (.not.GMM_IS_ERROR(istat)) then
                     outvar_L=.true.
                  end if
                  exit
               end if
            end do

            if (outvar_L) then

               w4(:,:,1:G_nk) = tr1(:,:,1:G_nk)
               if (out3_sfcdiag_L) then
                  if (trim(Tr3d_name_S(indxtr))=='HU') then
                     istat = gmm_get(gmmk_diag_hu_s,qdiag)
                     if (istat == 0) w4(:,:,G_nk+1) = qdiag(:,:)
                  else
                     w4(:,:,G_nk+1) = tr1(:,:,G_nk)
                     istat = phy_get (ptr3d, trim(fullname), F_npath='VO', F_bpath='D',&
                                      F_start=lijk, F_end=uijk, F_quiet=.true.)
                  end if
               else
                  istat=-1
               end if
               if (istat == 0) then
                  call vertint2 ( tr5,cible,nko, w4 ,pw_log_pt,G_nk+1,&
                                   l_minx,l_maxx,l_miny,l_maxy       ,&
                           1,l_ni,1,l_nj, inttype=Out3_vinterp_type_S )
               else
                  call vertint2 ( tr5,cible,nko, w4,pw_log_pt,G_nk,&
                                  l_minx,l_maxx,l_miny,l_maxy     ,&
                           1,l_ni,1,l_nj, inttype=Out3_vinterp_type_S )
               end if

               if (Outd_filtpass(ii,set) > 0) &
                    call filter2( tr5,Outd_filtpass(ii,set), &
                                      Outd_filtcoef(ii,set), &
                            l_minx,l_maxx,l_miny,l_maxy, nko)
               if (Out3_cliph_L) then
                  do k=1,nko
                  do j=1,l_nj
                  do i=1,l_ni
                     tr5(i,j,k) = max(tr5(i,j,k), 0. )
                  end do
                  end do
                  end do
               end if

               call out_fstecr3 ( tr5,l_minx,l_maxx,l_miny,l_maxy,rf      , &
                                  Outd_var_S(ii,set),Outd_convmult(ii,set), &
                                  Outd_convadd(ii,set),knd,-1 , &
                                  nko, indo, nko               , &
                                  Outd_nbit(ii,set) , .false. )

            end if

         end do

         deallocate(indo,rf,prprlvl,tr5,cible)

      end if
!
!----------------------------------------------------------------------
!
      return
      end
