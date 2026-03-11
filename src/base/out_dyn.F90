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

!**s/r out_dyn - perform dynamic output

      subroutine out_dyn ( F_reg_out, F_casc_L )
      use dynkernel_options
      use gem_options
      use init_options
      use grdc_options
      use glb_ld
      use svro_mod
      use levels
      use lun
      use outd
      use outgrid
      use out_listes
      use out_mod
      use out_vref
      use outusrdir
      use step_options
      use omp_timing
      implicit none

      logical F_reg_out, F_casc_L

#include <rmnlib_basics.hf>

      character(len=15) prefix
      logical ontimec,flag_clos,near_sfc_L
      integer kk,jj,levset,gridset,usrdirset,istat
      integer :: nko,indo(10000)
!
!----------------------------------------------------------------------
!
      if (.not.Lun_debug_L) istat= fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)

      Out_type_S   = 'REGDYN'
!
!########## REGULAR OUTPUT #######################################
!
      if (F_reg_out) then

         if (outd_sorties(0,Lctl_step) < 1) goto 887

         if (Lun_out > 0) then
            write(Lun_out,7001) Lctl_step,trim(Out_laststep_S)
         end if

         call gtmg_start ( 80, 'OUT_DYN', 1)

         call canonical_cases ("OUT")
         
         do jj=1, outd_sorties(0,Lctl_step)
            Out_nfstecr= 0
            kk       = outd_sorties(jj,Lctl_step)
            gridset  = Outd_grid(kk)
            levset   = Outd_lev(kk)
            usrdirset= Outd_usrdir(kk)
            if(OutGrid_hgrid_usr(gridset)%usr_grid_L .or. Level_vgrid_usr(levset)%usr_grid_L )then
               Out_prefix_S(1:2) = 'u'//OutUsrdir_name_S(usrdirset)(4:4)
            else
               Out_prefix_S(1:1) = 'd'
               Out_prefix_S(2:2) = Level_typ_S(levset)
            endif
            Out_prefix_S(3:4) = '  '
            call up2low (Out_prefix_S ,prefix)
            Out_prefix_S = prefix
            Out_reduc_l       = OutGrid_reduc(gridset)

            Out_stride = 1      ! can only be one for now
            Out_gridi0 = max( 1   , OutGrid_x0 (gridset))
            Out_gridin = min( G_ni, OutGrid_x1 (gridset))
            Out_gridj0 = max( 1   , OutGrid_y0 (gridset))
            Out_gridjn = min( G_nj, OutGrid_y1 (gridset))

            if ( .not. OUTs_server_L) then
               call out_open_file (trim(prefix))

               call out_href ( 'Mass_point'                , &
                  OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                  OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )

            
               if (Level_typ_S(levset) == 'M') then
                  call out_vref_itf (etiket=Out_etik_S)
                  call out_slev (Level(1,levset), Level_max(levset), &
                              G_nk,indo,nko,near_sfc_L)
               else if (Level_typ_S(levset) == 'P') then
                  call out_vref_itf (Level_allpres(1:Level_npres),etiket=Out_etik_S)
               else if (Level_typ_S(levset) == 'H') then
                  call out_vref_itf (Level_allheights(1:Level_nheights),etiket=Out_etik_S,agl_L=.true.)
               end if
            endif

            call OUTs_metaS ()
            
            call out_tracer (levset, kk)

            if ( trim(Dynamics_Kernel_S) == 'DYNAMICS_FISL_H' ) then
               call out_thm_hlt (levset, kk)
            else
               call out_thm    (levset, kk)
            endif

            call out_uv     (levset, kk)

            call out_dq     (levset, kk)

            call out_gmm    (levset, kk)

            if ( .not. OUTs_server_L) then
               flag_clos= .true.
               if (jj < outd_sorties(0,Lctl_step)) then
                  flag_clos= .not.( (gridset == Outd_grid(outd_sorties(jj+1,Lctl_step))).and. &
                  (Level_typ_S(levset) == Level_typ_S(Outd_lev(outd_sorties(jj+1,Lctl_step)))))
               end if
               if (flag_clos) call out_cfile ()
            endif

            call OUTs_metaF (Out_nfstecr, OUTs_nvar_indx)

            !if (Lun_out > 0) write (6,'(a,3i5,(" ",a4))')  'TREATING: ',jj,Out_nfstecr,Out_nplans,Outd_varnm_S(1:Outd_var_max(kk),kk),Out_prefix_S
         end do
         call gtmg_stop ( 80 )
      end if

!#################################################################
!
 887  continue
!
!########## SPECIAL OUTPUT FOR CASCADE ###########################

      if ((F_casc_L) .and. (Grdc_ndt > 0)) then

         ontimec = .false.
         if ( Lctl_step >= Grdc_start.and.Lctl_step <= Grdc_end ) then
            ontimec = (mod(Lctl_step+Grdc_start,Grdc_ndt) == 0)
         end if

         if ( Init_mode_L .and. (Step_kount >= Init_halfspan) ) then
            ontimec = .false.
         end if

         if ( ontimec ) then

            call gtmg_start ( 80, 'OUT_DYN', 1)
            if ( .not. OUTs_server_L) call out_open_file ('casc')
            call out_dyn_casc()
            if ( .not. OUTs_server_L) call out_cfile()
            call gtmg_stop ( 80 )
            
         end if

         ontimec = .false.

      end if

      istat = fstopc('MSGLVL','WARNIN',RMN_OPT_SET)

 7001 format(/,' OUT_DYN- WRITING DYNAMIC OUTPUT FOR STEP (',I8,') in directory: ',a)
!
!--------------------------------------------------------------------
!
      return
      end

