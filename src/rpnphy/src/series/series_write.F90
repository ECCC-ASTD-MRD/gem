!-------------------------------------- LICENCE BEGIN -------------------------
!Environment Canada - Atmospheric Science and Technology License/Disclaimer,
!                     version 3; Last Modified: May 7, 2008.
!This is free but copyrighted software; you can use/redistribute/modify it under the terms
!of the Environment Canada - Atmospheric Science and Technology License/Disclaimer
!version 3 or (at your option) any later version that should be found at:
!http://collaboration.cmc.ec.gc.ca/science/rpn.comm/license.html
!
!This software is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
!without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!See the above mentioned License/Disclaimer for more details.
!You should have received a copy of the License/Disclaimer along with this software;
!if not, you can write to: EC-RPN COMM Group, 2121 TransCanada, suite 500, Dorval (Quebec),
!CANADA, H9P 1J3; or send e-mail to service.rpn@ec.gc.ca
!-------------------------------------- LICENCE END ---------------------------

module series_write_mod
   use, intrinsic :: iso_fortran_env, only: REAL64
   use iso_c_binding
   use rpn_comm_itf_mod
   use wb_itf_mod, only: wb_get
   use tdpack_const, only: GRAV, RGASD
   use series_options
   use phygridmap, only: phydim_nk
   use vGrid_Descriptors, only: vgrid_descriptor,vgd_get,vgd_free,VGD_OK,VGD_ERROR
   use vgrid_wb, only: vgrid_wb_get
   use ptopo_utils, only: ptopo_init_var, ptopo_grid_ipe, ptopo_isblocmaster_L
   use mu_jdate_mod, only: jdate_to_print
   implicit none
   private

   public :: series_write_geo, series_write

!!!#include <arch_specific.hf>
#include <rmnlib_basics.hf>
#include <rmn/msg.h>

contains

   !@author Ron McTaggart-Cowan, 2009-04
   !@revision 2017-01, S.Chamberland - Major revision

   !#TODO: subroutine series_write_fst(F_force_L)

   !/@*
   subroutine series_write_geo()
      implicit none
      !@objective write geo fields to binary file
      !*@/
      integer :: istat, nn, k, l, m, nk, istep, series_ngeop
      real(REAL64) :: heure_8
      character(len=1024) :: msg_S
      !---------------------------------------------------------------
      if (series_paused_L .or. .not.(series_initok_L .and. series_on_L)) return
      if (series_nstnb <= 0) return

      nk = phydim_nk - 1
      series_ngeop = 0
      heure_8 = 0.

      call priv_open()

      if (ptopo_isblocmaster_L) series_gdata_out = 0.
      call rpn_comm_reduce( &
           series_gdata, series_gdata_out, size(series_gdata),  &
           RPN_COMM_INTEGER, RPN_COMM_BOR, RPN_COMM_MASTER, &
           RPN_COMM_BLOC_COMM, istat)

      IF_MASTER: if (ptopo_isblocmaster_L) then

         call msg(MSG_INFO, PKGNAME_S//'Write Geop header')
         write(msg_S, *) series_ngeo, ' : ', (trim(P_serg_srgeo_s(nn))//', ', nn=1,series_ngeo)
         call msg(MSG_INFO, PKGNAME_S//'Geop Fields: '//trim(msg_S))

         call priv_write_hdr(P_serg_srgeo_s, P_serg_srgeo_s, series_ngeo, series_ngeop, nk)

         write(msg_S,'(i6,a,f9.2,a,i4,a,i4,a,i4,a,i4,a)') series_kount, ",", &
              0., "h (n=", series_nstnb, '; s=', series_ngeo, '; p=', 0, &
              '; nk=',nk,')'
         call msg(MSG_INFO, PKGNAME_S//'Write geophy data at:'//trim(msg_S))

         !# Testing code
!!$            write(series_fileid) heure_8
!!$            do m=1,series_ngeo
!!$               print *,PKGNAME_S//'('//P_serg_srgeo_s(m)//'',minval(series_gdata_out(m,:)),maxval(series_gdata_out(m,:))
!!$               write(msg_S,'(i6.6,a,a4,a)') series_kount, ' (', P_serg_srgeo_s(m), ') '
!!$               write(series_fileid) '----surf-0-'//msg_S(1:13)//'----'
!!$               write(series_fileid) (series_gdata_out(m,l), l=1,series_nstnb)
!!$               write(series_fileid) '----surf-1-'//msg_S(1:13)//'----'
!!$            enddo
!!$            do m=1,series_ngeop
!!$               write(msg_S,'(i6.6,a,a4,a)') series_kount, ' (', prof_S(m), ') '
!!$               write(series_fileid) '----prof-0-'//msg_S(1:13)//'----'
!!$               write(series_fileid) ((series_pdata_out(m,l,k,istep), k=1,nk), l=1,series_nstnb)
!!$               write(series_fileid) '----prof-1-'//msg_S(1:13)//'----'
!!$            enddo

         !# Old Series replication code
         !# Write Date
         write(series_fileid) heure_8, &
              ((series_gdata_out(m,l), l=1,series_nstnb), m=1,series_ngeo), &
              (((series_pdata_out(m,l,k,istep), k=1,nk), l=1,series_nstnb), m=1,series_ngeop)


         !# Future code
!!$            write(series_fileid) -1, VERSION_S
!!$            write(series_fileid) &
!!$                 series_stng, series_nstnb, &
!!$                 series_ngeo, series_nsurf, series_nprof, nk, &
!!$                 size(vtbl_t,1), size(vtbl_t,2), size(vtbl_t,3), &
!!$                 size(ip1_t), size(ip1_m), &
!!$                 series_moyhr, series_delt, series_jdate, series_satuco_L
!!$            !#TODO: F_etiket,(F_ig(k),k=1,4), F_dgrw,F_rgas,F_grav,F_satues_L,F_satuco_L,tsmoyhr,srwri
!!$            write(series_fileid) &
!!$                 (series_stnb(nn)%name, &
!!$                 series_stnb(nn)%lat, series_stnb(nn)%lon, &
!!$                 series_stnb(nn)%ig, series_stnb(nn)%jg, series_stnb(nn)%gidx, &
!!$                 nn=1,series_nstnb)
!!$            write(series_fileid) &
!!$                 (P_serg_srgeo_s(nn), nn=1,series_ngeo), &
!!$                 (P_serg_srsrf_s(nn), nn=1,series_nsurf), &
!!$                 (P_serg_srprf_s(nn), nn=1,series_nprof)
!!$            write(series_fileid) &
!!$                 vtbl_t, ip1_t, ip1_m

      endif IF_MASTER

      !# init buffers
      series_sdone = 0
      series_pdone = 0
      series_sdone2 = 0
      series_pdone2 = 0
      series_sdata = 0.
      series_pdata = 0.
      !---------------------------------------------------------------
      return
   end subroutine series_write_geo


   !/@*
   subroutine series_write(F_force_L)
      implicit none
      !@objective Generate a time series binary dump file.
      !@arguments
      logical, intent(in) :: F_force_L
      !*@/

      integer :: istat, nn, istep, nsteps, kount, k, l, m, nk
      real(REAL64) :: heure_8
      character(len=1024) :: msg_S
      character(len=SER_STRLEN_VAR) :: surf_S(NVARMAX), prof_S(NVARMAX)
      !---------------------------------------------------------------
      !#TODO: protect from multi thread

      if (series_paused_L .or. .not.(series_initok_L .and. series_on_L)) return
      if (series_nstnb <= 0) return
      !#TODO: adjust for not writing at every extraction
!!$      if (mod(series_kount, series_out_nsteps) /= 0 .and. .not.F_force_L) return
      if (mod(series_kount, P_serg_srwri) /= 0 .and. .not.F_force_L) return

      nk = phydim_nk - 1

      call priv_open()

      !# Write header for surf+prof data
      !#TODO: adjust for not writing at every extraction
!!$      IF_KOUNT_1: if (series_kount == series_out_nsteps) then
      IF_KOUNT_1: if (series_kount == P_serg_srwri) then
         if (ptopo_isblocmaster_L) series_sdone_out = 0.
         call rpn_comm_reduce( &
              series_sdone, series_sdone_out, size(series_sdone),  &
              RPN_COMM_INTEGER, RPN_COMM_BOR, RPN_COMM_MASTER, &
              RPN_COMM_BLOC_COMM, istat)

         if (ptopo_isblocmaster_L) series_pdone_out = 0.
         call rpn_comm_reduce( &
              series_pdone, series_pdone_out, size(series_pdone),  &
              RPN_COMM_INTEGER, RPN_COMM_BOR, RPN_COMM_MASTER, &
              RPN_COMM_BLOC_COMM, istat)

         IF_MASTER1: if (ptopo_isblocmaster_L) then

            surf_S(:) = ' '
            where(series_sdone_out > 0)
               surf_S = P_serg_srsrf_s
            end where
            prof_S(:) = ' '
            where(series_pdone_out > 0)
               prof_S = P_serg_srprf_s
            end where

            call msg(MSG_INFO, PKGNAME_S//'Write Surf+Prof header')
            write(msg_S, *) series_nsurf, ' : ', (trim(surf_S(nn))//', ', nn=1,series_nsurf)
            call msg(MSG_INFO, PKGNAME_S//'Surface Fields: '//trim(msg_S))
            write(msg_S, *) series_nprof, ' : ', (trim(prof_S(nn))//', ', nn=1,series_nprof)
            call msg(MSG_INFO, PKGNAME_S//'Profile Fields: '//trim(msg_S))

            call priv_write_hdr(surf_S, prof_S, series_nsurf, series_nprof, nk)
         endif IF_MASTER1

      endif IF_KOUNT_1


      !# Write time series data

      if (ptopo_isblocmaster_L) series_sdata_out = 0.
      call rpn_comm_reduce( &
           series_sdata, series_sdata_out, size(series_sdata),  &
           RPN_COMM_INTEGER, RPN_COMM_BOR, RPN_COMM_MASTER, &
           RPN_COMM_BLOC_COMM, istat)

      if (ptopo_isblocmaster_L) series_pdata_out = 0.
      call rpn_comm_reduce( &
           series_pdata, series_pdata_out, size(series_pdata),  &
           RPN_COMM_INTEGER, RPN_COMM_BOR, RPN_COMM_MASTER, &
           RPN_COMM_BLOC_COMM, istat)

      IF_MASTER2: if (ptopo_isblocmaster_L) then
         nsteps = series_out_nsteps
         !#TODO: if not mod(series_kount, series_out_nsteps) == 0, nstep < series_out_nsteps
         do istep = 1, nsteps
            kount = series_kount - nsteps + istep  !#TODO: check this
            heure_8 = dble(kount) * dble(series_delt)/3600.d0

            write(msg_S,'(i6,a,f9.2,a,i4,a,i4,a,i4,a,i4,a)') kount, ",", &
                 heure_8, "h (n=", series_nstnb, '; s=', series_nsurf, &
                 '; p=', series_nprof, '; nk=', nk,')'
            call msg(MSG_INFO, PKGNAME_S//'Write surf+prof data at: '//trim(msg_S))

            !# Testing code
!!$            write(series_fileid) heure_8
!!$            do m=1,series_nsurf
!!$               write(msg_S,'(i6.6,a,a4,a)') series_kount, ' (', surf_S(m), ') '
!!$               write(series_fileid) '----surf-0-'//msg_S(1:13)//'----'
!!$               write(series_fileid) (series_sdata_out(m,l,istep), l=1,series_nstnb)
!!$               write(series_fileid) '----surf-1-'//msg_S(1:13)//'----'
!!$            enddo
!!$            do m=1,series_nprof
!!$               write(msg_S,'(i6.6,a,a4,a)') series_kount, ' (', prof_S(m), ') '
!!$               write(series_fileid) '----prof-0-'//msg_S(1:13)//'----'
!!$               write(series_fileid) ((series_pdata_out(m,l,k,istep), k=1,nk), l=1,series_nstnb)
!!$               write(series_fileid) '----prof-1-'//msg_S(1:13)//'----'
!!$            enddo

            !# Old Series replication code
            write(series_fileid) heure_8, &
                 ((series_sdata_out(m,l,istep), l=1,series_nstnb), m=1,series_nsurf), &
                 (((series_pdata_out(m,l,k,istep), k=1,nk), l=1,series_nstnb), m=1,series_nprof)

            !# Future code
!!$            write(series_fileid) heure_8, kount, &
!!$                 series_sdata(:, :, istep), series_pdata(:, :, 1:nk, istep)
!!$            !#TODO: write reduced series_sdone, series_pdone
         enddo
      endif IF_MASTER2

      !# Reset buffer after writing
      series_sdone = 0
      series_pdone = 0
      series_sdone2 = 0
      series_pdone2 = 0
      series_sdata = 0.
      series_pdata = 0.
      !---------------------------------------------------------------
      return
   end subroutine series_write


   !#---- Private functions ------------------------------------------


   !/@*
   subroutine priv_open()
      implicit none
      !@objective Open series files
      !*@/
      integer :: istat
      character(len=128) :: filename_S
      !---------------------------------------------------------------
      call ptopo_init_var()

      if (series_fileid == -1 .and. ptopo_isblocmaster_L) then
         
         write(filename_S,'(a,i6.6)') '../time_series_'// &
              trim(jdate_to_print(series_jdate))//'.bin_', ptopo_grid_ipe
         series_fileid = 0
         if (series_kount == 0) then
            call msg(MSG_INFO, PKGNAME_S//'Open: '//trim(filename_S))
            istat = fnom(series_fileid, trim(filename_S), 'SEQ+FTN+UNF', 0)
         else
            call msg(MSG_INFO, PKGNAME_S//'Re-Open: '//trim(filename_S))
            istat = fnom(series_fileid, trim(filename_S), 'SEQ+FTN+UNF+OLD+APPEND', 0)
         endif
         !#TODO: handle error
      endif
      !---------------------------------------------------------------
      return
   end subroutine priv_open


   !/@*
   subroutine priv_write_hdr(surf_S, prof_S, nsurf, nprof, nk)
      implicit none
      !@objective 
      !@arguments
      character(len=*), intent(in)  :: surf_S(:), prof_S(:)
      integer, intent(in) :: nsurf, nprof, nk
      !*@/
      integer, dimension(:), pointer :: ip1_t, ip1_m
      real(REAL64), dimension(:,:,:), pointer :: vtbl_t
      integer :: istat, nn, hgc(4)
      type(vgrid_descriptor) :: vcoord
      real :: dgrw
      character(len=12) :: etk_S
      !---------------------------------------------------------------
      nullify(ip1_m, ip1_t, vtbl_t)
      istat = vgrid_wb_get('phy-t', vcoord, ip1_t)
      istat = vgd_get(vcoord, 'VTBL - vgrid_descriptor table', vtbl_t)
      istat = vgd_free(vcoord)
      if (associated(ip1_t)) deallocate(ip1_t,stat=istat)
      !#TODO: why next 3 lines
!!$      istat = vgrid_wb_get('phy-t', vcoord, ip1_m)
!!$      istat = vgd_free(vcoord)
!!$      if (associated(ip1_m)) deallocate(ip1_m,stat=istat)

      dgrw = 0.
      etk_S = ' '
      hgc(:) = 0
      istat = wb_get('model/Output/etik', etk_S)
      istat = wb_get('model/Hgrid/hgcrot', hgc, nn)
!!$      istat = wb_get('model/Vgrid/size-hybm', series_nkm)
!!$      istat = wb_get('model/Vgrid/size-hybt', series_nkt)

      ! Put the character 'G' in front of etiket to allow correct rotation
      ! of wind related variables by feseri program
      etk_S(1:1) = 'G'
 
      write(series_fileid) series_nstnb, nsurf, nprof, nk
      write(series_fileid) size(vtbl_t,1),size(vtbl_t,2),size(vtbl_t,3)
      write(series_fileid) vtbl_t
      write(series_fileid) series_nstnb, &
           (series_stnb(nn)%name, series_stnb(nn)%ig, &
           series_stnb(nn)%jg, nn=1, series_nstnb), &
           nsurf, (surf_S(nn), nn=1,nsurf), &
           nprof, (prof_S(nn), nn=1,nprof), &
           (series_cmcdateo14(nn), nn=1,14), etk_S, &
           (hgc(nn), nn=1,4), dgrw, RGASD, GRAV, &
           series_satues_L, series_satuco_L, series_moyhr, &
           series_interval_sec
      if (associated(vtbl_t)) deallocate(vtbl_t,stat=istat)
      !---------------------------------------------------------------
      return
   end subroutine priv_write_hdr

end module series_write_mod
