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
! Fichier/File   : chm_load_emissions.ftn90
! Creation       : A. Kallaur (MSC/ARQI), Octobre 2006
! Description    : Read from file (if need be) for major point sources.
!
! Modification1   : December 2010 (A. Kallaur)
!                   Adaptation to GEM4.4.0+ new interface.
!                   Chemical buses no longer exist, so all emissions data
!                   is written to the permanent Physics bus.
!
!                   Oct. 2016 Jack Chen
!                   if (chm_get_wf_emis_l=.true.) open separate mjr-point
!                   file specific for fire emissions from
!                   CFFEPS output including new hourly plumeirse parameters
!
!                   November 2017 (Verica Savic-Jovcic)
!                   To remove hard-coded name of the file with the linoz
!                   coefficients, introduced here a call to linoz_litcoeff
!                   that was previously called from chm_nml.
!                   To do: Move this call to new subroutine that will be called
!                          from physics. Get GEM group involved in that change.
!
!                   November 2019 (Deji Akingunola)
!                   Add area emission(s) variables to the (Physics) input list
!
!==============================================================================
!!if_on
subroutine chm_load_emissions2(F_basedir_S, gem_tstep_num, inputobj, nbvar_input)
   use inputio_mod,          only: INPUTIO_T
!!if_off
   use chm_utils_mod,        only: global_debug, chm_lun_out
   use chm_nml_mod,          only: chm_master, chm_emis_master, chm_get_mj_emis_l, &
                                   chm_get_wf_emis_l, em_nesdt, chm_do_mjpts_l,    &
                                   chm_step_factor, chm_cffeps_online_l
   use chm_datime_mod,       only: ijul_day, imonth, chm_dttim_s
   use chm_mjrpts_sortinfo_mod, only:  nb_me_species, lnb_sources, ndata,      &
                                   i_gilc, i_gjlc, i_hgt, i_dia, i_tem, i_vel, &
                                   i_typ, Lstack_Info, Lstack_emis
   use chm_headers_mod,      only: chm_mjrpts_get_stkinfo, chm_mjrpts_get_size, &
                                   chm_mjrpts_get_emissions, chm_cffeps_init,   &
                                   chm_fst_openfile, chm_fst_closefile
   use chm_species_info_mod, only: species_master, nb_dyn_tracers
   use ezgrid_mod,           only: ezgrid_find_ij0
   use phygridmap,           only: phy_lcl_gid, phy_glb_gid,       &
                                   mapmod2phy, phy_lcl_ni, phydim_ni
   use mu_jdate_mod,         only: jdate_day_of_year, jdate_to_cmc, jdate_to_print
   use phy_options,          only: delt, jdateo, MU_JDATE_HALFDAY
   use inputio_mod,          only: inputio_add, inputio_nbvar
   use rpn_comm_itf_mod

   implicit none
!!if_on
   character(len=*)               :: F_basedir_S  !- base path for input data file
   integer(kind=4), intent(in)    :: gem_tstep_num
   type(INPUTIO_T), intent(inout) :: inputobj
   integer(kind=4), intent(inout) :: nbvar_input
!!if_off
!
! Local variables
!
   character(len=512)                           :: Major_point_file, Fire_point_file
   integer(kind=4)                              :: datev, ier, err, int_x, int_y, istatus
   real(kind=4),    save                        :: nsec = 0.
   real(kind=4)                                 :: nsec_tmp
   logical(kind=4)                              :: local_dbg, did_stack_info, init_done
   logical(kind=4), save                        :: print_once = .true.
   logical(kind=4), save                        :: chm_incfg_init_l = .false.
   integer(kind=4)                              :: mjp_unf !Major_point_file unit number
   integer(kind=4)                              :: fir_unf !Fire_point_file unit number
   integer(kind=4)                              :: i, il
   integer(kind=4)                              :: i0, j0, lni, lnj, ij(2)
   integer(kind=4)                              :: lx_lbound, lx_ubound, ly_lbound, ly_ubound
   integer(kind=4), save                        :: my_pe
   integer(kind=4)                              :: numproc
   integer(kind=4)                              :: lnb_sources_x
   integer(kind=4), save                        :: na_sources    ! Anthropogenic sources
   integer(kind=4), save                        :: nf_sources    ! Fire sources
   integer(kind=4), save                        :: nb_sources
   integer(kind=4), dimension(:),   allocatable :: pe_info, loc_src
   integer(kind=4), dimension(:,:), allocatable :: pi0_j0, rpi0_j0
   real(kind=4),    dimension(:),   allocatable :: lat, lon, hgt, dia, tem, vel, typ
   real(kind=4),    dimension(:),   allocatable :: lat2x, lon2y
! Time dependent fire plume rise parameters from CFFEPS
   real(kind=4),    dimension(:),   allocatable :: zplm, rsmk
   real(kind=4),    dimension(:,:), allocatable :: Stack_Emis_a, Stack_Emis_f
   real(kind=4),    dimension(:,:), pointer     :: Stack_Info, Stack_Emis, Vstack_info
   type(rpncomm_context), save                  :: mjr_context=NULL_rpncomm_context
   integer(kind=4)                              :: ihour
   integer(kind=8)                              :: jdatev, last_read, nsec_i8
   character(len=15)                            :: emiss_read_time
!
   character(len=256)                           :: incfg_S
   character(len=4)                             :: inname_S, em_nesdt_c
!  Declaration of external functions and subroutines
   external physeterror, rpn_comm_reduce
   integer(kind=4), external                    :: gdxyfll
!
!
!  Detect Master switch. If false, NORMAL EXIT WITH MESSAGE
!
   if (.not. chm_master) then
      if (print_once) then
         if (chm_lun_out > 0) write(chm_lun_out, *) &
           'CHM_LOAD_EMISSIONS -> DETECTED CHEMICAL MASTER KILL: NO EMISSIONS, NORMAL EXIT'
         print_once = .false.
      end if
      return
   end if
!
!  Set the dates for which the emissions are to be read
!
   nsec_i8  = gem_tstep_num * int(delt)
   jdatev   = jdateo + nsec_i8
   datev    = jdate_to_cmc(jdatev)
   ijul_day = jdate_day_of_year(jdatev + MU_JDATE_HALFDAY)
   chm_dttim_s = adjustl(jdate_to_print(jdatev))
   read(chm_dttim_s(5:6), '(I2)', iostat=ier) imonth
   if (ier < 0 ) then
      call physeterror('chm_load_emissions', 'problem with setting date')
      return
   end if
   read(chm_dttim_s(10:11), '(I2)', iostat=ier) ihour

   if ((mod(gem_tstep_num, chm_step_factor) == 0) .and. (chm_lun_out > 0)) then
      write (chm_lun_out, 1001) gem_tstep_num
   end if

!  Detect EMISSIONS Master switch. If false, NORMAL EXIT WITH MESSAGE -> NO EMISSIONS READ IN
!
   if (.not. chm_emis_master) then
      if (print_once) then
         if (chm_lun_out > 0) write(chm_lun_out, *) &
           'CHM_LOAD_EMISSIONS -> DETECTED EMISSIONS MASTER KILL: NO EMISSIONS, NORMAL EXIT'
         print_once = .false.
      end if
      return
   end if
!
   if (.not. chm_incfg_init_l) then
      ! Add area emission species to the physics input list
      incfg_S = '(chm_load_emissions) Adding area emission species to the input list'
      if (chm_lun_out > 0) write(chm_lun_out, *) trim(incfg_S)
      write(em_nesdt_c, '(i4.4)') nint(em_nesdt / 60.0)
      incfg_S = 'freq=MINUTE,0,'//em_nesdt_c//'; search=ANAL,MINPUT; hinterp=near;'
      istatus = 0
      do i = 1, nb_dyn_tracers
! Skip if specie is sea-salt
         if (species_master(i) % ae_name(1:3) == 'ESS') cycle

         if (species_master(i) % ae_offset > 0) then
            inname_S = species_master(i) % ae_name
            istatus = min(inputio_add(inputobj%cfg, 'in='//inname_S//'; '// &
                                  & trim(incfg_S)), istatus)
         end if
         if (species_master(i) % fae_offset > 0) then
            inname_S = species_master(i) % fae_name
            istatus = min(inputio_add(inputobj%cfg, 'in='//inname_S//'; '// &
                                  & trim(incfg_S)), istatus)
         end if
         if (species_master(i) % mae_offset > 0) then
            inname_S = species_master(i) % mae_name
            istatus = min(inputio_add(inputobj%cfg, 'in='//inname_S//'; '// &
                                  & trim(incfg_S)), istatus)
         end if
!
      end do
!
      nbvar_input = inputio_nbvar(inputobj)
      if (istatus < 0) then
         call physeterror('chm_load_emissions', 'error adding area emissions inputs')
         return
      end if

      chm_incfg_init_l = .true.
   end if
!
!  Detect Major Points and Wildfire Emissions switch
!
   if (.not. (chm_get_mj_emis_l .or. chm_get_wf_emis_l)) then
      if (print_once) then
         if (chm_lun_out > 0) write(chm_lun_out, *) &
           'CHM_LOAD_EMISSIONS -> DETECTED NO MAJOR POINTS EMISSIONS, NORMAL EXIT'
         print_once = .false.
      end if
      return
   else
!   Detect if Anthropogenic emissions are enabled
      if (chm_get_mj_emis_l) then
         Major_point_file = trim(F_basedir_S)//'/MODEL_INPUT/Mjrpts_emissions.fst'
         mjp_unf = 11
         if (print_once) then
            if (chm_lun_out > 0) write(chm_lun_out, *) &
              'CHM_LOAD_EMISSIONS -> MAJOR POINTS EMISSIONS FILE is ', trim(Major_point_file)
         end if
      end if
!
!   Detect if get fire emissions as major point source
      if (chm_get_wf_emis_l) then
         Fire_point_file = trim(F_basedir_S)//'/MODEL_INPUT/Mjrpts_fire_emissions.fst'
         fir_unf = 12
         if (print_once) then
            if (chm_lun_out > 0) write(chm_lun_out, *) &
         'CHM_LOAD_EMISSIONS -> MAJOR POINTS _FIRE_ EMISSIONS ', trim(Fire_point_file)
         end if
      end if
      if (print_once) print_once = .false.
   end if

! Set debug switch

   local_dbg = ((.false. .or. global_debug) .and. (chm_lun_out > 0))
!
   if (local_dbg) then
      write(*,*) 'Inside chm_load_emissions with chm_lun_out= ',chm_lun_out
   end if
!
! Get stack info
!
   did_stack_info = .false.
   IF_FIRST: if (RPN_COMM_is_null(mjr_context)) then

      did_stack_info = .true.

      nullify(Lstack_info)

      call rpn_comm_rank (RPN_COMM_GRID, my_pe, err)
      call rpn_comm_size (RPN_COMM_GRID, numproc, err)
      allocate (pi0_j0(5, numproc), rpi0_j0(5, numproc))
!
      pi0_j0 = 0
      rpi0_j0 = 0

      err = ezgrid_find_ij0(phy_lcl_gid, phy_glb_gid, i0, j0, lni, lnj)

      pi0_j0(1, my_pe+1) = i0
      pi0_j0(2, my_pe+1) = j0
      pi0_j0(3, my_pe+1) = lni
      pi0_j0(4, my_pe+1) = lnj
      pi0_j0(5, my_pe+1) = my_pe
      call rpn_comm_reduce(pi0_j0, rpi0_j0, 5*numproc, "MPI_INTEGER", &
                           "MPI_SUM", 0, RPN_COMM_GRID, err)

      IF_PE0: if (my_pe == 0) then

         if (local_dbg) then
            write(chm_lun_out, *) ' '
            write(chm_lun_out, *) '-------BEGINNING of TIME=0 setup of chm_load_emissions ------- '
            write(chm_lun_out, *) '==============================================================='
         end if

         if (chm_get_mj_emis_l) then
            call chm_fst_openfile(mjp_unf, Major_point_file, 'RND+OLD+R/O', 'RND', istatus)
            if (istatus < 0) then
               call physeterror('chm_load_emissions', 'problem with Major_point_file')
               return
            end if

            call chm_mjrpts_get_size(mjp_unf, na_sources, istatus)
            if (istatus < 0 ) then
               call physeterror('chm_load_emissions', 'problem with chm_mjrpts_get_size')
               return
            end if
         else
            na_sources = 0
         end if
!
         if (chm_get_wf_emis_l) then
            call chm_fst_openfile(fir_unf, Fire_point_file, 'RND+OLD+R/O', 'RND', istatus)
            if (istatus < 0) then
               call physeterror('chm_load_emissions', 'problem with Fire_point_file')
               return
            end if

            call chm_mjrpts_get_size(fir_unf, nf_sources, istatus)
            if (istatus < 0) then
               call physeterror('chm_load_emissions', 'error getting fire emission sources')
               return
            end if
         else
            nf_sources = 0
         end if

         nb_sources = na_sources + nf_sources

         if (local_dbg) then
            write(chm_lun_out, *) 'In CHM_LOAD_EMISSIONS '
            write(chm_lun_out, *) 'number of anth. major point sources: ', na_sources
            write(chm_lun_out, *) 'number of FIRE point sources: ', nf_sources
            write(chm_lun_out, *) 'Total number of sources: ', nb_sources
         end if

         allocate(lat(nb_sources), lon(nb_sources), hgt(nb_sources), &
                  dia(nb_sources), tem(nb_sources), vel(nb_sources), typ(nb_sources))
         allocate(lat2x(nb_sources), lon2y(nb_sources))

         if (chm_get_mj_emis_l) then
            call chm_mjrpts_get_stkinfo(mjp_unf, na_sources, istatus,         &
                                        lat(1:na_sources), lon(1:na_sources), &
                                        tem(1:na_sources), hgt(1:na_sources), &
                                        dia(1:na_sources), vel(1:na_sources), &
                                        typ(1:na_sources))
            if (istatus < 0) then
               call physeterror('chm_load_emissions', 'error getting stack info')
               return
            end if
!  To prevent double counting of fire emissions, do not allow for fire emissions to be present in both files
            if (chm_get_wf_emis_l .and. (any(typ(1:na_sources) > 490.0))) then
               write(chm_lun_out, *) '*** ERROR ***'
               write(chm_lun_out, *) 'No fire emissions allowed in Majro_point_file when Fire_point_file provided.'
               write(chm_lun_out, *) '*** ABORT ***'
               call physeterror('chm_load_emissions', 'possible double counting of fire emissions')
               return
            end if
         end if

         if (chm_get_wf_emis_l) then
            call chm_mjrpts_get_stkinfo(fir_unf, nf_sources, istatus, &
                                        lat(na_sources+1:nb_sources), &
                                        lon(na_sources+1:nb_sources))
            if (istatus < 0) then
               call physeterror('chm_load_emissions', 'error getting fire locations')
               return
            end if
      ! Set stack height, diameter, temperature and exit velocity to -1.0
      ! for fire plumes. Time dependent plume height and smoke to be set height
      ! and exit velocity (respectively) later below.
            do i = na_sources+1, nb_sources
               dia(i) = -1.0
               tem(i) = -1.0
               hgt(i) = -1.0
               vel(i) = -1.0
               typ(i) = 500.0
            end do
         end if

         if (local_dbg) then
            write(chm_lun_out, *) ' '
            write(chm_lun_out, *) 'stack info read        '
         end if

         where(lon < 0.0) lon = lon + 360.0
         ier = gdxyfll(phy_glb_gid, lat2x, lon2y, lat, lon, nb_sources)

         allocate(loc_src(numproc))
         allocate(pe_info(nb_sources))
         allocate(Stack_Info(nb_sources, ndata))
         loc_src = 0
         pe_info = -1
         Stack_Info = -1.0
         do i = 1, nb_sources
!          Round to the nearest i-j point
            int_x = nint(lat2x(i))
            int_y = nint(lon2y(i))
!
            do il = 1, numproc
               lx_lbound = rpi0_j0(1,il)                        ! i0
               lx_ubound = rpi0_j0(1,il) + rpi0_j0(3,il) - 1    ! in
               ly_lbound = rpi0_j0(2,il)                        ! j0
               ly_ubound = rpi0_j0(2,il) + rpi0_j0(4,il) - 1    ! jn

               if (int_x .ge. lx_lbound .and. int_x .le. lx_ubound  .and. &
                   int_y .ge. ly_lbound .and. int_y .le. ly_ubound) then
                  pe_info(i)  = rpi0_j0(5, il) ! pe
                  loc_src(il) = loc_src(il) + 1
                  exit
               end if
            end do

            Stack_Info(i, i_gilc) = real(int_x)   ! i-index on full physics grid
            Stack_Info(i, i_gjlc) = real(int_y)   ! j-index on full physics grid
            Stack_Info(i, i_hgt)  = hgt(i)        ! height of stack
            Stack_Info(i, i_dia)  = dia(i)        ! diameter of stack
            Stack_Info(i, i_tem)  = tem(i)        ! exit temperature of stack plume
            Stack_Info(i, i_vel)  = vel(i)        ! exit velocity of stack plume
            Stack_Info(i, i_typ)  = typ(i)        ! Emission type
         end do
         if (local_dbg) then
            write(chm_lun_out,*)' -- for each of the ', numproc,    &
                                ' processors, respective local sources are ', loc_src
            do il = 1, numproc
               lx_lbound = rpi0_j0(1, il)                        ! i0
               lx_ubound = rpi0_j0(1, il) + rpi0_j0(3, il) - 1    ! in
               ly_lbound = rpi0_j0(2, il)                        ! j0
               ly_ubound = rpi0_j0(2, il) + rpi0_j0(4, il) - 1    ! jn
               write(chm_lun_out, *)' -- for processor ', il - 1 ,  &
                                    ' g_i0, g_in, g_j0, g_jn are ', &
                                    lx_lbound, lx_ubound, ly_lbound, ly_ubound
            end do
         end if
         if (chm_get_mj_emis_l) call chm_fst_closefile(mjp_unf)
         if (chm_get_wf_emis_l) call chm_fst_closefile(fir_unf)

         deallocate(lat, lon, hgt, dia, tem, vel, typ, lat2x, lon2y)
         deallocate(loc_src)
      else ! IF_PE0:
         allocate(pe_info(1))
         allocate(Stack_Info(1, 1))
      end if IF_PE0

!
      istatus = 0
      err = RPN_COMM_spread_context(mjr_context, RPN_COMM_GRID, 0, pe_info, &
                                    nb_sources)
      if (err /= 0) then
         call physeterror('chm_load_emissions', 'problem in spread_context')
         return
      end if

      if (local_dbg) then
         if (.not. RPN_COMM_is_null(mjr_context)) then
            write(chm_lun_out,*)'pointer mjr_context is now associated '
         end if
      end if

!    From global array Stack_Info of dimension (nb_sources,ndata) obtain
!         local POINTER Lstack_Info for dimension (ndata,lnb_sources)
!         since some sources may not be part of the domain Sum(lnb_sources) /= nb_sources

      lnb_sources = RPN_COMM_spread(mjr_context, Stack_Info, nb_sources, &
                                    ndata, Lstack_info)
      if (local_dbg) then
         if (associated(Lstack_info)) then
            write(chm_lun_out, *) 'Lstack_info is now associated with dimension ', shape(Lstack_info)
         end if
      end if

      if (lnb_sources < 0) then
         call physeterror('chm_load_emissions', 'problem in rpn_comm_spread')
         return
      end if
      chm_do_mjpts_l = (lnb_sources > 0)

      if (local_dbg) then
         write(chm_lun_out, *)' -- for processor ', my_pe, ' chm_do_mjpts_l is ', chm_do_mjpts_l
         write(chm_lun_out, *)' -- for processor ', my_pe, ' RPN_COMM_spread returned ', lnb_sources , ' sources'
      end if

      deallocate(pe_info)
      deallocate(Stack_Info)

!    remap global i/j indices to local indices
      do i = 1, lnb_sources
         int_x = int(lstack_info(i_gilc, i)) - i0 + 1
         int_y = int(lstack_info(i_gjlc, i)) - j0 + 1
         ij = mapmod2phy(int_x, int_y, phy_lcl_ni, phydim_ni)
         lstack_info(i_gilc, i) = real(ij(1))
         lstack_info(i_gjlc, i) = real(ij(2))
      end do
!    Ensure stack height is at least 0.8m
      lstack_info(i_hgt, :) = max(lstack_info(i_hgt, :), 0.8)
!
     ! Initialize online CFFEPS
      if (chm_cffeps_online_l .and. (gem_tstep_num == 0)) then
         istatus = chm_cffeps_init(F_basedir_S, my_pe, numproc, rpi0_j0)
         if (istatus < 0) then
            call physeterror('chm_load_emissions', 'problem in chm_cffeps init')
            return
         end if
      end if

      deallocate (pi0_j0, rpi0_j0)
      nullify(Lstack_Emis)
      if (local_dbg) then
         write(chm_lun_out, *) ' '
         write(chm_lun_out, *) '-------END of TIME=0 setup of chm_load_emissions ------- '
         write(chm_lun_out, *) '============================================================'
      end if
!
   end if IF_FIRST

   nsec_tmp = real(nsec_i8)
   init_done = nsec_tmp < nsec
   nsec = nsec_tmp

! read emissions at em_nesdt interval. Except for restart when we may have to read previous emissions.

   IF_READ : if (mod(nsec, em_nesdt) == 0 .or. did_stack_info .or. init_done) then

      if (local_dbg) then
         write(chm_lun_out, *) ' '
         write(chm_lun_out, *) '-------BEGINNING of emission read in chm_load_emissions ------- '
         write(chm_lun_out, *) '-------for step ', gem_tstep_num
         write(chm_lun_out, *) '==============================================================='
         if (associated(Lstack_emis)) then
            write(chm_lun_out, *) 'Lstack_emis is already associated with dimension ', shape(Lstack_emis)
         end if
      end if

      err = 0

      if (my_pe == 0) then
!
         emiss_read_time = chm_dttim_s
         if (mod(nsec, em_nesdt) /= 0.0) then
            write(chm_lun_out, *) '-------possible restart or end of digital filter......  '
            last_read = floor(nsec / em_nesdt)
            jdatev = jdateo + int(last_read * em_nesdt)
            datev  = jdate_to_cmc(jdatev)
            emiss_read_time = jdate_to_print(jdatev)
         end if
!
         allocate(Stack_Emis(nb_sources, nb_me_species))
!        if (local_dbg) then
         write(chm_lun_out, *) ' '
         write(chm_lun_out, *) '-------reading major point emissions for date.time = ', emiss_read_time
!        end if

         if (chm_get_mj_emis_l) then
            call chm_fst_openfile(mjp_unf, Major_point_file, 'RND+OLD+R/O', 'RND', ier)

            allocate(Stack_Emis_a(na_sources, nb_me_species))
            call chm_mjrpts_get_emissions(mjp_unf, Stack_Emis_a, na_sources, datev, err)
            Stack_Emis(1:na_sources, 1:nb_me_species) = Stack_Emis_a(1:na_sources, 1:nb_me_species)
            deallocate(Stack_Emis_a)
            if (err < 0) then
               call physeterror('chm_load_emissions', 'problem with reading Anthropogenic emissions')
               return
            end if

            call chm_fst_closefile(mjp_unf)
         end if

         if (chm_get_wf_emis_l) then
            call chm_fst_openfile(fir_unf, Fire_point_file, 'RND+OLD+R/O', 'RND', ier)

            allocate(Stack_Emis_f(nf_sources, nb_me_species))
            ! for time dependent fire plume parameter
            allocate(zplm(nf_sources), rsmk(nf_sources))
            allocate(Stack_Info(nb_sources, 2))
!
!  Read emissions from Fire_point_file along with the plume height and smoke parameters
            call chm_mjrpts_get_emissions(fir_unf, Stack_Emis_f, nf_sources, datev, err, &
                                          rsmk, zplm, 'RSMK', 'ZPLM')
            if (err < 0) then
               call physeterror('chm_load_emissions', 'problem with reading fire emissions')
               return
            end if
            Stack_Info = -1.0
            do i = 1, nf_sources
               Stack_Emis(na_sources + i, 1:nb_me_species) = Stack_Emis_f(i, 1:nb_me_species)
               Stack_Info(na_sources + i, 1) = zplm(i)
               Stack_Info(na_sources + i, 2) = rsmk(i)
            end do
            deallocate(Stack_Emis_f, zplm, rsmk)

            call chm_fst_closefile(fir_unf)
         end if

      else
         if (chm_get_wf_emis_l) allocate(Stack_Info(1, 1))
         allocate(Stack_Emis(1, 1))
      end if

!
!    From global array Stack_Emis of dimension (nb_sources,nb_me_species) obtain
!         local POINTER Lstack_emis for dimension (nb_me_species,lnb_sources_x)
!         lnb_sources should be equal to lnb_sources_x
      lnb_sources_x = RPN_COMM_spread(mjr_context, Stack_Emis, nb_sources, nb_me_species, lstack_emis)
      if (local_dbg) then
         if (associated(Lstack_emis)) then
            write(chm_lun_out, *) 'Lstack_emis is associated with dimension ', shape(Lstack_emis)
         end if
      end if
      if (lnb_sources_x /= lnb_sources) then
         call physeterror('chm_load_emissions', 'inconsistent number of local sources')
         return
      end if

      deallocate(Stack_Emis)
!
      if (chm_get_wf_emis_l) then
!====Time varying fire plume info
!
         nullify(Vstack_Info)
         lnb_sources_x = RPN_COMM_spread(mjr_context, Stack_Info, nb_sources, &
                                         2, Vstack_Info)
         if (local_dbg) then
            if (associated(Vstack_info)) then
               write(chm_lun_out, *)'PE:', my_pe, 'Vstack_info is now associated with dimension ', shape(Vstack_info)
            end if
         end if
         if (lnb_sources_x /= lnb_sources) then
            call physeterror('chm_load_emissions', 'inconsistent number of local sources')
            return
         end if
         do i = 1, lnb_sources
            if (Lstack_info(i_dia, i) < 0.0) then
               Lstack_info(i_hgt, i) = Vstack_info(1, i)
               Lstack_info(i_vel, i) = Vstack_info(2, i)
            end if
         end do

         deallocate(Stack_Info)
         deallocate(Vstack_Info)
      end if

   end if IF_READ

 1001 format(/,'CHEMISTRY : PERFORMING TIMESTEP #',I9, &
             /,'========================================')
   return

end subroutine chm_load_emissions2
