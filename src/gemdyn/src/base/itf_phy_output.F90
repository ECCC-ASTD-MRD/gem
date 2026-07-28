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

!**s/r itf_phy_output

      subroutine itf_phy_output (stepno)
      use vertical_interpolation
      use vGrid_Descriptors
      use vgrid_wb
      use out_vref
      use phy_itf
      use step_options
      use gmm_pw
      use HORgrid_options
      use gem_options
      use out_options
      use glb_ld
      use svro_mod
      use lun
      use levels
      use out3
      use outgrid
      use outp
      use outusrdir
      use out_listes
      use out_mod
      use out_meta
      use rmn_gmm
      use set_level_mod, only: set_level_usr_val,set_level_ERROR
      use, intrinsic :: iso_fortran_env     
      implicit none

      integer, intent(IN) :: stepno

#include <rmnlib_basics.hf>

      type(phymeta) :: pmeta
      type(vgrid_descriptor) :: vcoord
      character(len=15) prefix, model_var_stag_S
      integer i,ii,jj,kk,levset,usrdirset,nko,nko_pres,cnt,istat,&
              gridset,mult, knd , kind,&
              p_li0,p_li1,p_lj0,p_lj1,last_timestep
      integer grille_x0,grille_x1,grille_y0,grille_y1
      integer, dimension(:), allocatable :: indo,irff
      integer, dimension(:), pointer     :: ip1m,indo_pres
      logical flag_clos, write_diag_lev, accum_L
      real(kind=REAL64) avgfact
      real, dimension (l_ni,l_nj,G_nk+1), target :: wlnpi_m,wlnpi_t
      real, dimension(:), pointer    :: hybm,hybt,rf
      real, dimension(:), allocatable:: rff
      real, dimension(:,:,:), pointer :: lnpres,ptr3d,cible,usr_src,cible_dyn,usr_src_dyn
      real, dimension(:,:,:), allocatable         :: buso_pres
      real, dimension(:,:,:), allocatable, target :: data3d, zero
      real hybt_gnk2(1),hybm_gnk2(1)
      integer ind0(1)
!
!----------------------------------------------------------------------
!
      nullify(cible,usr_src,cible_dyn,usr_src_dyn,indo_pres,rf)
      if (outp_sorties(0,stepno) <= 0) then
         return
      else
         if (Lun_out > 0) then
            write(Lun_out,7001) stepno,trim(Out_laststep_S)
         end if
      end if
      call gtmg_start ( 48, 'PHY_output', 40 )

      istat = fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)
      out_type_S   = 'REGPHY'

      p_li0= Grd_lphy_i0 ; p_li1=Grd_lphy_in
      p_lj0= Grd_lphy_j0 ; p_lj1=Grd_lphy_jn

!     setup domain extent to retrieve physics data
      allocate ( data3d(l_ni,l_nj,G_nk+1), zero(p_li0:p_li1,p_lj0:p_lj1,G_nk+1), &
                 rff(Outp_multxmosaic), irff(Outp_multxmosaic))
      data3d= 0. ; zero= 0.

      istat= gmm_get(gmmk_pw_log_pm_s, pw_log_pm)
      istat= gmm_get(gmmk_pw_log_pt_s, pw_log_pt)
      wlnpi_m(1:l_ni,1:l_nj,:)= pw_log_pm(1:l_ni,1:l_nj,:)
      wlnpi_t(1:l_ni,1:l_nj,:)= pw_log_pt(1:l_ni,1:l_nj,:)

!     Retrieeve vertical coordinate description
      nullify(ip1m,hybm,hybt)
      istat = vgrid_wb_get('ref-m',vcoord,ip1m)
      deallocate(ip1m); nullify(ip1m)
      if (vgd_get(vcoord,'VCDM - vertical coordinate (m)',hybm, quiet=.true.) /= VGD_OK) istat = VGD_ERROR
      if (vgd_get(vcoord,'VCDT - vertical coordinate (t)',hybt, quiet=.true.) /= VGD_OK) istat = VGD_ERROR
      istat = vgd_free(vcoord)
      hybt_gnk2(1)= hybt(G_nk+2)
      hybm_gnk2(1)= hybm(G_nk+2)
      ind0(1) = 1

      do jj=1, outp_sorties(0,stepno)
         
         Out_nfstecr= 0
         kk       = outp_sorties(jj,stepno)
         gridset  = Outp_grid(kk)
         levset   = Outp_lev(kk)
         usrdirset= Outp_usrdir(kk)
         accum_L  = (Outp_avg_L(kk).or.Outp_accum_L(kk))
         avgfact  = 1.d0/dble(max(1,lctl_step-Outp_lasstep(kk,stepno)))
         last_timestep= -1
         if (accum_L) then
            last_timestep= max(0,Outp_lasstep(kk,stepno))
         end if
         allocate ( indo( min(Level_max(levset),Level_momentum ) ) )

         call out_slev (Level(1,levset), Level_max(levset), &
                         Level_momentum,indo,nko,write_diag_lev)

         if (Level_typ_S(levset) == 'P') then
            nko_pres = Level_max(levset)
         end if
         
         if(OutGrid_hgrid_usr(gridset)%usr_grid_L .or. Level_vgrid_usr(levset)%usr_grid_L )then
            Out_prefix_S(1:2) = 'u'//OutUsrdir_name_S(usrdirset)(4:4)
            Out_prefix_S(3:4) = '  '
         else
            Out_prefix_S(1:1) = 'p'
            Out_prefix_S(2:2) = Level_typ_S(levset)
            Out_prefix_S(3:3) = ' '
            Out_prefix_S(4:4) = OutGrid_hgrid_usr(gridset)%usr_grid_index_S
         endif
         call up2low (Out_prefix_S(1:2),prefix)
         Out_reduc_l       = OutGrid_reduc(gridset)

         grille_x0 = max( 1   +Grd_bsc_ext1, OutGrid_x0(gridset) )
         grille_x1 = min( G_ni-Grd_bsc_ext1, OutGrid_x1(gridset) )
         grille_y0 = max( 1   +Grd_bsc_ext1, OutGrid_y0(gridset) )
         grille_y1 = min( G_nj-Grd_bsc_ext1, OutGrid_y1(gridset) )

         Out_stride = 1         ! can only be one for now
         Out_gridi0 = max( 1   , grille_x0)
         Out_gridin = min( G_ni, grille_x1)
         Out_gridj0 = max( 1   , grille_y0)
         Out_gridjn = min( G_nj, grille_y1)

         if ( .not. OUTs_server_L) then

            call out_open_file (trim(prefix))

            call out_href ( 'Mass_point',grille_x0,grille_x1,1,&
                                       grille_y0,grille_y1,1 )

            if (Level_typ_S(levset) == 'M') then
               call out_vref_itf (etiket=Out_etik_S)
            elseif (Level_typ_S(levset) == 'P') then
               call out_vref_itf (Level_allpres(1:Level_npres),&
                               etiket=Out_etik_S)
            elseif (Level_typ_S(levset) == 'H') then
               call out_vref_itf (Level_allheights(1:Level_nheights),&
                               etiket=Out_etik_S,agl_L=.true.)
            end if
         endif
         
         call OUTs_metaS ()
         
         PHYSICS_VARS: do ii=1, Outp_var_max(kk)

            WRITE_FIELD: if (phy_getmeta (pmeta, Outp_var_S(ii,kk), &
                             F_npath='O',F_bpath='PVED', F_quiet=.true.)&
                             > 0 ) then
               FIELD_SHAPE: if (pmeta%nk == 1) then ! 2D field
                  ! For surface fields, encode the ip1 style ('N' or 'O') in Out_stag_S(3:3)
                  Out_stag_S= 'MSN'//OutGrid_hgrid_usr(gridset)%usr_grid_index_S
                  if(Level_vgrid_usr(levset)%vcode .eq. 1002 )Out_stag_S(3:3)='O'
                  rff(1)= 0. ; irff(1)= 1 ; knd= 2
                  if ( pmeta%fmul > 1 ) then
                     do mult=1,pmeta%fmul
                        rff(mult)= mult
                        irff(mult)= mult
                     end do
                     knd= 3
                  end if
                  cnt= pmeta%fmul

                  ptr3d => data3d(p_li0:p_li1,p_lj0:p_lj1,1:cnt)
                  istat = phy_get ( ptr3d, Outp_var_S(ii,kk), &
                                    F_npath='O', F_bpath='PVED')
                  if (Outp_avg_L(kk)) data3d = data3d*avgfact
                  call out_fstecr ( data3d, 1,l_ni, 1,l_nj, rff  ,&
                          Outp_var_S(ii,kk),Outp_convmult(ii,kk)  ,&
                          Outp_convadd(ii,kk), knd, last_timestep,&
                          cnt,irff,cnt,Outp_nbit(ii,kk),.false. )

                  if (accum_L) then
                      ptr3d => zero(:,:,1:cnt)
                      istat = phy_put ( ptr3d, Outp_var_S(ii,kk),&
                                        F_npath='O', F_bpath='PV')
                  end if

               else ! 3D field

                  ptr3d => data3d(p_li0:p_li1,p_lj0:p_lj1,:)
                  istat = phy_get (ptr3d,Outp_var_S(ii,kk),F_npath='O', &
                                   F_bpath='PVD')
                  if (Outp_avg_L(kk)) data3d = data3d*avgfact

                  if (Level_typ_S(levset) == 'M') then
                     if (pmeta%stag > 0) then ! thermo
                        Out_stag_S= 'MT '//OutGrid_hgrid_usr(gridset)%usr_grid_index_S
                        call out_fstecr (data3d                       ,&
                                 1,l_ni, 1,l_nj, hybt                  ,&
                                 Outp_var_S(ii,kk),Outp_convmult(ii,kk),&
                                 Outp_convadd(ii,kk),Level_kind_ip1,last_timestep,&
                                 G_nk,indo,nko,Outp_nbit(ii,kk),.false. )
                        if (write_diag_lev) then
                           Out_stag_S(3:3)= 'D '//OutGrid_hgrid_usr(gridset)%usr_grid_index_S
                           call out_fstecr (data3d(1,1,G_nk+1)        ,&
                                 1,l_ni, 1,l_nj, hybt_gnk2             ,&
                                 Outp_var_S(ii,kk),Outp_convmult(ii,kk),&
                                 Outp_convadd(ii,kk),Level_kind_diag,last_timestep,&
                                 1,ind0,1,Outp_nbit(ii,kk),.false. )
                        end if
                     else  ! momentum
                        Out_stag_S= 'MM '//OutGrid_hgrid_usr(gridset)%usr_grid_index_S
                        call out_fstecr (data3d                       ,&
                                 1,l_ni, 1,l_nj, hybm                  ,&
                                 Outp_var_S(ii,kk),Outp_convmult(ii,kk),&
                                 Outp_convadd(ii,kk),Level_kind_ip1,last_timestep,&
                                 G_nk,indo,nko,Outp_nbit(ii,kk),.false. )
                        if (write_diag_lev) then
                           Out_stag_S(3:3)= 'D '//OutGrid_hgrid_usr(gridset)%usr_grid_index_S
                           call out_fstecr (data3d(1,1,G_nk+1)        ,&
                                 1,l_ni, 1,l_nj, hybm_gnk2             ,&
                                 Outp_var_S(ii,kk),Outp_convmult(ii,kk),&
                                 Outp_convadd(ii,kk),Level_kind_diag,last_timestep,&
                                 1,ind0,1,Outp_nbit(ii,kk),.false. )
                        end if
                     end if

                  else
                     ! Output on pressure, heights AGL or user levels
                     
                     ! Note: The pointers cible_dyn, usr_src_dyn, indo_pres, and rf
                     ! will point to memory allocated within set_level_usr_val.
                     ! This "internal" memory allocation will persist 
                     ! throughout the model integration and will be expanded if necessary. 
                     ! Therefore, cible_dyn, usr_src_dyn, indo_pres, and rf must not
                     ! be deallocated, but they can be nullified if needed.
                     model_var_stag_S='MOMENTUM'
                     if ( pmeta%stag > 0 )model_var_stag_S='THERMO'
                     Out_stag_S='M???'

                     if( set_level_usr_val(cible_dyn, indo_pres, rf, kind, nko_pres, usr_src_dyn, &
                          levset,Level,Level_max, model_var_stag_S,&
                          l_minx,l_maxx,l_miny,l_maxy, G_nk, &
                          Out_stag_S, Level_typ_S(levset), Outp_grid(kk)) &
                          == set_level_ERROR )then
                        print*,'TODO itf_phy_output handle error gracefully 1'
                        stop
                        return
                     end if

                     ! 3D Arrays returned by function set_level_usr_val
                     ! have the following scope, l_minx,l_maxx,l_miny,l_maxy
                     ! We take on the 1:l_ni,1:l_nj part.
                     cible   =>   cible_dyn(1:l_ni,1:l_nj,1:nko_pres)
                     usr_src => usr_src_dyn(1:l_ni,1:l_nj,1:G_nk)
                     
                     allocate(buso_pres(l_ni,l_nj,nko_pres))

                     call vertint2 ( buso_pres, cible, nko_pres, data3d,&
                                     usr_src, G_nk, 1,l_ni, 1,l_nj      ,&
                                     1,l_ni, 1,l_nj, inttype=Out3_vinterp_type_S,&
                                     levtype=Level_vgrid_usr(levset)%class_S)
                     
                     call out_fstecr ( buso_pres, 1,l_ni, 1,l_nj      ,&
                           rf,Outp_var_S(ii,kk)           ,&
                           Outp_convmult(ii,kk),Outp_convadd(ii,kk),kind,last_timestep,&
                           nko_pres,indo_pres,nko_pres,Outp_nbit(ii,kk),&
                           .false. )

                     nullify(cible_dyn, cible, usr_src_dyn, usr_src, indo_pres, rf)
                     deallocate(buso_pres)
                     
                  end if
                  if (accum_L) then
                      ptr3d => zero
                      istat = phy_put (ptr3d,Outp_var_S(ii,kk),F_npath='O', &
                                       F_bpath='PV')
                  end if
               end if FIELD_SHAPE
            end if WRITE_FIELD
         end do PHYSICS_VARS

         deallocate (indo)

         if ( .not. OUTs_server_L) then
            flag_clos= .true.
            if (jj < outp_sorties(0,stepno)) then
               flag_clos= .not.( (gridset == Outp_grid(outp_sorties(jj+1,stepno))).and. &
                 (Level_typ_S(levset) == Level_typ_S(Outp_lev(outp_sorties(jj+1,stepno)))))
            end if
            if (flag_clos) call out_cfile ()
         endif

         call OUTs_metaF (Out_nfstecr, OUTs_nvar_indx)

      end do

      deallocate(rff,irff,data3d,zero)
      deallocate(hybm,hybt); nullify(hybm,hybt)
      call gtmg_stop  ( 48 )

      istat = fstopc('MSGLVL','WARNIN',RMN_OPT_SET)

 7001 format(/,' OUT_PHY- WRITING PHYSICS OUTPUT FOR STEP (',I8,') in directory: ',a)
!
!----------------------------------------------------------------------
!
      return
      end


subroutine itf_phy_output_list(stepno)
   use wb_itf_mod, only: wb_get, wb_put, WB_REWRITE_MANY
   use clib_itf_mod, only: clib_tolower
   use dimout, only: MAXELEM, MAXSET
   use out_listes, only: outp_sorties
   use outp, only: Outp_var_max, Outp_var_S
   implicit none
#include <arch_specific.hf>

   integer, intent(in) :: stepno

#include <rmnlib_basics.hf>
#include <rmn/msg.h>

   character(len=32) :: varlist_S(MAXELEM*MAXSET)
   integer :: nn, jj, kk, ii, istat

   integer, save :: n0 = -1
   !----------------------------------------------------------------------
   if (n0 == -1) then
      istat = wb_get('itf_phy/PHYOUT', varlist_S, n0)
      n0 = max(1, n0)
   endif

   ! Collect the list of requested output physics vars
   nn = 0
   varlist_S(1:n0) = ' '
   do jj=1, outp_sorties(0,stepno)
      kk = outp_sorties(jj,stepno)
      do ii=1, Outp_var_max(kk)
         nn = nn + 1
         varlist_S(nn) = Outp_var_S(ii,kk)
         istat = clib_tolower(varlist_S(nn))
         if (nn > 1) then
            if (any(varlist_S(nn) == varlist_S(1:nn-1))) nn = nn - 1
         endif
      enddo
   enddo

   istat = wb_put('itf_phy/PHYSTEPOUT_N', nn, WB_REWRITE_MANY)
   istat = wb_put('itf_phy/PHYSTEPOUT_V', varlist_S(1:n0), WB_REWRITE_MANY)

   !----------------------------------------------------------------------
   return
end subroutine itf_phy_output_list
