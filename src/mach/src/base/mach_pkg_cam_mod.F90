!---------------------------------- LICENCE BEGIN -------------------------------
! GEM-MACH - Atmospheric chemistry library for the GEM numerical atmospheric model
! Copyright (C) 2007-2013 - Air Quality Research Division &
!                           National Prediction Operations division
!                           Environnement Canada
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
! Lesser General Public License for more details.
!
! You should have received a copy of the GNU Lesser General Public
! License along with this library; if not, write to the Free Software
! Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
!---------------------------------- LICENCE END ---------------------------------


!============================================================================!
!         Environnement Canada         |        Environment Canada           !
!                                      |                                     !
! - Service meteorologique du Canada   | - Meteorological Service of Canada  !
! - Direction generale des sciences    | - Science and Technology Branch     !
!   et de la technologie               |                                     !
!============================================================================!
!                            http://www.ec.gc.ca                             !
!============================================================================!
!
! Projet/Project : GEM-MACH
! Fichier/File   : mach_pkg_cam_mod.ftn90
! Description    : Defines the PM species starting index (sp_AERO), and
!                  metadata for all the fields that are part of the CAM 
!                  aerosol package:
!                  "SULPHATE", "SEA-SALT", "OMCARBON", "NITRATES",
!                  "AMMONIUM", "SOILDUST", "BLCARBON", "PMCARBON",
!                  with the aliases:
!                  SU, SS, OC, NI, AM, CM, EC, PC
!
!============================================================================

module mach_pkg_cam_mod
   use chm_nml_mod,          only: chm_get_ae_emis_l, chm_get_mj_emis_l,  &
                                   chm_seaflux_s, chm_mono_s, chm_mass_s, &
                                   chm_get_fd_emis_l, chm_pkg_gas_s
   use chm_utils_mod,        only: pre_increment
   use mach_cam_utils_mod,   only: isize, icom, aeroname, aero_sname, mwt_aero
   use chm_species_idx_mod,  only: sp_AERO, sp_H2O2, sp_HNO3, sp_NH3, sp_O3, &
                                   sp_ROOH, sp_SO2, sp_SO4
   use chm_species_info_mod, only: species_master

   contains
!============================================================================
! Name           : pkg_cam_idxinit
!
! Description    : Initialize the PM species starting index; using the
!                  the index passed in as a dummy argument.
!
! Arguments:  IN/OUT
!                idx  -> index to start from
!============================================================================
      integer function pkg_cam_idxinit(idx)
         implicit none

         integer(kind=4), intent(inout) :: idx ! the index from where to start

         integer(kind=4) :: nb_fields, ngas

         ! Gas species required by CAM, if no gas-phase mechanism was selected
         ngas = 0
         if (chm_pkg_gas_s == 'NIL') then
            sp_H2O2 = pre_increment(idx)
            sp_HNO3 = pre_increment(idx)
            sp_NH3  = pre_increment(idx)
            sp_O3   = pre_increment(idx)
            sp_ROOH = pre_increment(idx)
            sp_SO2  = pre_increment(idx)
            sp_SO4  = pre_increment(idx)
            ngas = 7
         end if

         sp_AERO = idx + 1
         nb_fields = icom * isize
         
         idx = nb_fields + sp_AERO - 1

         pkg_cam_idxinit = nb_fields + ngas

      end function pkg_cam_idxinit


!============================================================================
! Name           : pkg_cam_metainit
!
! Description    : Initialize the meta information for each species
!
! Arguments:  None
!
! Extra info: Naming convention for bins:
!             - The aerosol specie short name are used to identify the species
!             - The other 2 are be used for the bin number
!
!============================================================================
      subroutine pkg_cam_metainit()
         implicit none
         integer           :: nbin, cam_emiss_nb_bins
         character(len=1)  :: zbin
         character(len=2)  :: cbin
         character(len=3)  :: sname
         character(len=10) :: lname
!
!        local variable
         integer(kind=4)   :: nsp, iae
         character(len=22) :: tr_prop_s

         tr_prop_s = ';min=0'
         select case (chm_mono_s)
           case ('CLIP')
             tr_prop_s = trim(tr_prop_s)//';monot=1'
           case ('ILMC')
             tr_prop_s = trim(tr_prop_s)//';monot=2'
           case ('NIL')
             continue
         end select
         select case (chm_mass_s)
           case ('BC')
             tr_prop_s = trim(tr_prop_s)//';massc=1'
           case ('NIL')
             continue
         end select

         cam_emiss_nb_bins = max(2, (isize - 2))

         do nbin = 1, isize
            write(zbin, '(Z1)') nbin
            write(cbin, '(i2.2)') nbin
            do iae = 1, icom
               sname = aero_sname(iae) // zbin
               lname = 'TR/T'//sname//':P'
               nsp = (iae - 1) * isize + nbin + sp_AERO - 1
               species_master(nsp) % dyn_name   = 'T' // sname
               species_master(nsp) % dyn_string = 'VN='//trim(lname)//trim(tr_prop_s)//     &
                                                & '; ON=T'//sname//'; VD='//aeroname(iae)// &
                                                & ' Bin '//cbin//' ug/kg; VS=T;VB=D1'
               species_master(nsp) % mol_wt     = mwt_aero(iae)
            end do
         end do

         if (chm_seaflux_s == 'GONG_MONAHAN_F') then
           ! iae = maxloc(merge(1, 0, aeroname == "SEA-SALT"), dim = 1)
            iae = findloc(aeroname, "SEA-SALT", dim = 1)
            do nbin = 1, isize
               write(zbin, '(Z1)') nbin
               write(cbin, '(i2.2)') nbin
               nsp = (iae - 1) * isize + nbin + sp_AERO - 1
               species_master(nsp) % ae_name   = 'ESS'//zbin
               species_master(nsp) % ae_string = 'VN=sp_SS_ESS'//zbin//'; ON=ESS'//zbin//   &
                                               & '; VD=Sea-salt emissions for Bin '//cbin// &
                                               & '(g/s); VS=A;VB=p0'
            end do
         end if

         if (chm_get_ae_emis_l .or. chm_get_fd_emis_l .or. chm_get_mj_emis_l) then
            do iae = 1, icom
               ! No emissions for secondary organic carbon, and
               ! sea-salt emissions are parameterized in the model
               if (any(aeroname(iae) == (/"SEA-SALT", "OMCARBON"/))) cycle

               do nbin = 1, cam_emiss_nb_bins  !no emissions for bins B and C
                  write(zbin, '(Z1)') nbin
                  write(cbin, '(i2.2)') nbin
                  sname = aero_sname(iae) // zbin
                  if (aero_sname(iae) == 'NI') sname = 'NT' // zbin
                  nsp = (iae - 1) * isize + nbin + sp_AERO - 1

                  ! Area emissions
                  if (chm_get_ae_emis_l) then
                     lname = 'sp_'// aero_sname(iae) // '_E' // sname
                     species_master(nsp) % ae_name   = 'E' // sname
                     species_master(nsp) % ae_string = 'VN='//lname//'; ON=E'//sname//'; VD='// &
                                                  & aeroname(iae)//' area emissions for Bin '//  &
                                                  & cbin//'(g/s); VS=A;VB=p1'
                  end if

                  ! Fugitive aerosol emissions
                  if (chm_get_fd_emis_l) then
                     lname = 'sp_'// aero_sname(iae) // '_F' // sname
                     species_master(nsp) % fae_name   = 'F' // sname
                     species_master(nsp) % fae_string = 'VN='//lname//'; ON=F'//sname//'; VD='// &
                                                  & aeroname(iae)//' fugitive emissions for Bin '//  &
                                                  & cbin//'(g/s); VS=A;VB=p1'
                  end if

                  ! Major point emissions
                  if (chm_get_mj_emis_l) then
                     species_master(nsp) % me_name = 'E' // sname
                  end if
               end do
            end do
         end if

         return
      end subroutine pkg_cam_metainit

end module mach_pkg_cam_mod
