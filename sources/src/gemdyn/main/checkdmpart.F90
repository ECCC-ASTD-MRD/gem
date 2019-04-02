      subroutine checkdmpart()
      use clib_itf_mod
      use gem_options
      use dyn_fisl_options
      use step_options
      use iso_c_binding
      use glb_ld
      use HORgrid_options
      use VERgrid_options
      use dynkernel_options
      use geomh
      use lun
      use path
      use ptopo
      use version
      implicit none

#include <rmnlib_basics.hf>
      include "rpn_comm.inc"

      external dummy_checkdm
      integer, external :: domain_decomp, sol_transpose, &
                           set_fft, gemdm_config

      character(len=16) ndomains_S,npex_S,npey_S
      character(len=2048) cdm_eigen_S,fn
      logical cdm_grid_L
      integer cdm_npex(2), cdm_npey(2), unf, cnt
      integer pe_xcoord(1000), pe_ycoord(1000)
      integer err,ierr(4),npex,npey,i,max_io_pes

      namelist /cdm_cfgs/ cdm_npex,cdm_npey,cdm_grid_L,cdm_eigen_S
!
!-------------------------------------------------------------------
!
      call init_component

      if (Ptopo_couleur == 0) then
         call open_status_file3 (trim(Path_input_S)//'/../checkdmpart_status.dot')
         call write_status_file3 ('checkdmpart_status=ABORT')
      endif

      cdm_npex   = 0
      cdm_npey   = 0
      cdm_grid_L = .false.
      cdm_eigen_S= 'NONE@#$%'

      Grd_yinyang_L = .false.
      Grd_yinyang_S = ''
      if (clib_getenv ('GEM_YINYANG',ndomains_S) >= 0) &
      Grd_yinyang_L = .true.
      Lun_out=6

      unf = 0
      if (fnom (unf,'./checkdm.nml', 'SEQ+OLD', 0) /= 0) goto 9110
      rewind(unf)
      read (unf, nml=cdm_cfgs, end = 9120, err = 9120)
      err= fclos (unf) ; goto 8888
 9110 if (Lun_out > 0) then
         write (Lun_out, 9050) 'checkdm.nml'
         write (Lun_out, 8000)
      endif
      goto 9999

 9120 err= fclos (unf)
      if (Lun_out >= 0) then
         write (Lun_out, 9150) 'cdm_cfgs','checkdm.nml'
         write (Lun_out, 8000)
      endif
      goto 9999
 8888 write (Lun_out,nml=cdm_cfgs)

      fn  = trim(Path_input_S)//'/model_settings.nml'
      unf = 0
      if (fnom (unf,fn, 'SEQ+OLD', 0) == 0) then
         if (HORgrid_nml (unf) < 0) goto 9999
         if (HORgrid_config (1) < 0) goto 9999
         err = HORgrid_nml (-1)

         if (VERgrid_nml (unf) < 0) goto 9999
         if (VERgrid_config () < 0) goto 9999
         err = VERgrid_nml (-1)

         if (step_nml  (unf) < 0) goto 9999
         err = step_nml  (-1)
         Step_runstrt_S='2011020300'

         if (dynKernel_nml (unf) < 0) goto 9999
         err = dynKernel_nml (-1)

         if ( Dynamics_Kernel_S(1:13) == 'DYNAMICS_FISL' ) then
            if (dyn_fisl_nml (unf) < 0) goto 9999
            err = dyn_fisl_nml (-1)
         endif

         if (gemdm_config () < 0) goto 9999
         err= fclos(unf)
      else
         print *, ' Namelist FILE: ',trim(fn),' NOT AVAILABLE'
         goto 9999
      endif

      cnt=0
      if ((cdm_npex(1) > 0) .and. (cdm_npey(1) > 0)) then
         if (cdm_npex(2)==0) cdm_npex(2)=cdm_npex(1)
         if (cdm_npey(2)==0) cdm_npey(2)=cdm_npey(1)
         do npey=cdm_npey(1), cdm_npey(2)
         do npex=cdm_npex(1), cdm_npex(2)
            cnt=cnt+1
            ierr=0
            ierr(1) = domain_decomp ( npex, npey, .true. )
            ierr(2) = sol_transpose ( npex, npey, .true. )
            if (minval(ierr) >= 0 ) then
               write(npex_S,'(i6)') npex
               write(npey_S,'(i6)') npey
               if (Ptopo_couleur == 0) &
               call write_status_file3 ('topo_allowed="'//trim(npex_S)//'x'//trim(npey_S)//'"')
            endif
         end do
         end do
         call gemtim4 ( Lun_out, 'AFTER domain_decomp', .false. )
      endif

      if (cnt == 1) then
         npex= cdm_npex(1) ; npey=cdm_npey(1)
         max_io_pes=npex*npey
         do i= 1, npex*npey
            err= RPN_COMM_io_pe_valid_set (pe_xcoord,pe_ycoord,i,&
            npex,npey,.false.,0)
            if (err /= 0) then
               max_io_pes= i-1
               exit
            endif
         end do
         write(npex_S,'(i6.6)') max_io_pes
         if (Ptopo_couleur == 0) call write_status_file3 ('MAX_PES_IO='//trim(npex_S))
         call gemtim4 ( Lun_out, 'AFTER io_pe_valid', .false. )
      endif

      if ( (cdm_grid_L) .or. (cdm_eigen_S /= 'NONE@#$%') ) then
         err = domain_decomp ( 1, 1, .false. )
         err = sol_transpose ( 1, 1, .false. )
         call glbpos
         call set_geomh
         call canonical_cases ("SET_GEOM")
         call write_status_file3 ('Fft_fast_L=OK')
         if (Ptopo_couleur == 0) then
            call write_status_file3 ('Fft_fast_L=OK')! TODO : no longer needed
         endif

         if (cdm_eigen_S /= 'NONE@#$%') then
            call set_opr (cdm_eigen_S)
            call gemtim4 ( Lun_out, 'AFTER set_opr', .false. )
         endif

         if (cdm_grid_L) then
            call set_params
            call set_dync ( .false., err )
            call gemtim4 ( Lun_out, 'AFTER set_dync', .false. )
            if (Ptopo_couleur == 0) then
               if (err == 0) then
                  call write_status_file3 ('SOLVER=OK')
               else
                  call write_status_file3 ('SOLVER=ABORT')
               endif
            endif
         endif

      endif

      if (Ptopo_couleur == 0) then
         call write_status_file3 ('checkdmpart_status=OK')
         call close_status_file3 ()
      endif

      call gemtim4 ( Lun_out, 'END OF CHECKDMPART', .true. )
      call memusage (Lun_out)

      if (Lun_out > 0) &
      err = exfin (trim(Version_title_S),trim(Version_number_S), 'OK')

 9999 call rpn_comm_FINALIZE(err)

 8000 format (/,'========= ABORT ============='/)
 9050 format (/,' FILE: ',A,' NOT AVAILABLE'/)
 9150 format (/,' NAMELIST ',A,' INVALID IN FILE: ',A/)
!
!-------------------------------------------------------------------
!
      return
      end

subroutine dummy_checkdm
return
end
