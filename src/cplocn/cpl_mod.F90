module cpl_mod

      logical cpl_ocn_L,          & ! Activate coupling with ocean
              cpl_wav_L,          & ! Activate coupling with waves
              cpl_dgflt_H,        & ! Digital filter second half if true
              cpl_rstn_L            ! Restart logic

      integer cpl_drv_gni, cpl_drv_gnj, cpl_drv_lni, cpl_drv_lnj,&
              cpl_drv_i0,  cpl_drv_j0,  cpl_drv_in,  cpl_drv_jn ,&
              cpl_drv_gnk

      integer cpl_minx, cpl_maxx, cpl_miny, cpl_maxy


      real    cpl_drv_delt

      character(len=17) cpl_rstn_S

      integer cplocn_ocnf_nsprd
      logical cplocn_debug_L,     &
              cplocn_iweight_L

end module cpl_mod

