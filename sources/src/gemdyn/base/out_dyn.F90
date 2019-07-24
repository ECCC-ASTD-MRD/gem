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
      use gem_options
      use init_options
      use grdc_options
      use levels
      use lun
      use outd
      use outgrid
      use out_listes
      use out_mod
      use out_vref, only: out_vref_itf
      use step_options
      use gem_timing
      implicit none
#include <arch_specific.hf>

      logical F_reg_out, F_casc_L

#include <rmnlib_basics.hf>

      character(len=15) prefix
      logical ontimec,flag_clos
      integer kk,jj,levset,gridset,istat
!
!----------------------------------------------------------------------
!
      call gemtime_start ( 80, 'OUT_DYN', 1)
      if (.not.Lun_debug_L) istat= fstopc('MSGLVL','SYSTEM',RMN_OPT_SET)

      Out_type_S   = 'REGDYN'
!
!########## REGULAR OUTPUT #######################################
!
      if (F_reg_out) then

         if (outd_sorties(0,Lctl_step) < 1) then
            if (Lun_out > 0) write(Lun_out,7002) Lctl_step
            goto 887
         end if

         if (Lun_out > 0) then
            write(Lun_out,7001) Lctl_step,trim(Out_laststep_S)
         end if

         call canonical_cases ("OUT")

         ! Precompute diagnostic level values
         call itf_phy_diag ()

         do jj=1, outd_sorties(0,Lctl_step)

            kk       = outd_sorties(jj,Lctl_step)
            gridset  = Outd_grid(kk)
            levset   = Outd_lev(kk)

            Out_prefix_S(1:1) = 'd'
            Out_prefix_S(2:2) = Level_typ_S(levset)
            call up2low (Out_prefix_S ,prefix)
            Out_reduc_l       = OutGrid_reduc(gridset)

            call out_open_file (trim(prefix))

            call out_href3 ( 'Mass_point'                , &
                  OutGrid_x0 (gridset), OutGrid_x1 (gridset), 1, &
                  OutGrid_y0 (gridset), OutGrid_y1 (gridset), 1 )

            if (Level_typ_S(levset) == 'M') then
               call out_vref_itf (etiket=Out_etik_S)
            else if (Level_typ_S(levset) == 'P') then
               call out_vref_itf (Level_allpres(1:Level_npres),etiket=Out_etik_S)
            end if

            call out_tracer (levset, kk)

            call out_thm    (levset, kk)

            call out_uv     (levset, kk)

            call out_dq     (levset, kk)

            call out_gmm    (levset, kk)

            flag_clos= .true.
            if (jj < outd_sorties(0,Lctl_step)) then
              flag_clos= .not.( (gridset == Outd_grid(outd_sorties(jj+1,Lctl_step))).and. &
              (Level_typ_S(levset) == Level_typ_S(Outd_lev(outd_sorties(jj+1,Lctl_step)))))
            end if

            if (flag_clos) call out_cfile

         end do

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

            call out_open_file ('casc')
            call out_dyn_casc
            call out_cfile

         end if

         ontimec = .false.

      end if

      istat = fstopc('MSGLVL','WARNIN',RMN_OPT_SET)
      call gemtime_stop ( 80 )

 7001 format(/,' OUT_DYN- WRITING DYNAMIC OUTPUT FOR STEP (',I8,') in directory: ',a)
 7002 format(/,' OUT_DYN- NO DYNAMIC OUTPUT FOR STEP (',I8,')')
!
!--------------------------------------------------------------------
!
      return
      end

