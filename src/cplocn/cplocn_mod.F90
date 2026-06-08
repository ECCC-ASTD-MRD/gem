module cplocn_mod

      use rmn_gmm, only: GMM_MAXNAMELENGTH

      logical cplocn_1st_L,      & ! Only one exchange at start
              cplocn_off_L,      & ! Simulate offline run
              cplocn_bzone_L       ! Activate buffer zone in coupling mask

      integer cplao_dt
      integer cplocn_n_fldou, cplocn_n_fldin
      integer cplocn_rap_dt,  cplocn_it
      integer cplocn_myproc,  cplocn_mycol, cplocn_myrow,  &
              cplocn_nprocs,                               &
              cplocn_atm_nw_max,                           &
              cplocn_gni, cplocn_gnj,                      &
              cplocn_nxchng, cplocn_atm_nxchng,            &
              cplocn_atm_nspread,                          &
              cplocn_dynphy_offset
                         ! Offset between GEM global model grid and physics global grid
                         ! in LAM mode
      integer cplao_xchg_mode          ! 0=gossip like, i.e., senc/rev at the same time, including step=0
                                       ! 1=send done at the end of the physical time step, except at step=0
                                       ! 2=same as 1 with no exchange at start of step=0, just a send at the end

      real,    parameter :: cplocn_missval=-9999.

      character(len=16) cplocn_runstrt_S

      parameter (cplocn_n_fldou=10)
      character(len=3), dimension(cplocn_n_fldou) :: cplocn_cvou_S
      character(len=1), dimension(cplocn_n_fldou) :: cplocn_cvot_S
      parameter (cplocn_n_fldin=38)
      character(len=3), dimension(cplocn_n_fldin) :: cplocn_cvin_S
      character(len=1), dimension(cplocn_n_fldin) :: cplocn_cvit_S

! Atmospheric model sends
      data cplocn_cvou_S /'FBA','FIA','RTA','TTA',            &
                           'UUA','VVA','QQA','PMA','PTA','P0A'/
      data cplocn_cvot_S /'S'  ,'S'  ,'S'  ,'S'  ,            &
                           'U'  ,'V'  ,'S'  ,'S'  ,'S'  ,'S'  /
      character(len=39) :: cplocn_iris_provides
      data cplocn_iris_provides /'FBA,FIA,RTA,TTA,UUA/VVA,QQA,PMA,PTA,P0A'/

! Ocean model received
      data cplocn_cvin_S /'MCP','ALO','ALI','T4O','T4I',  &
                          'SHO','SHI','LHO','LHI',        &
                          'TXO','TYO','TXI','TYI',        &
                          'ZTO','ZQO','ZUO','ZVO',        &
                          'ZTI','ZQI','ZUI','ZVI',        &
                          'GLI','I8I','SDI','TMO',        &
                          'UUO','VVO','I7I','UUI','VVI',  &
                          'QSO','QSI','ILO','ILI',        &
                          'ZMO','ZMI','ZHO','ZHI'      /
      data cplocn_cvit_S /'S'  ,'S'  ,'S'  ,'S',  'S'  ,  &
                          'S'  ,'S'  ,'S'  ,'S'  ,        &
                          'U'  ,'V'  ,'U'  ,'V'  ,        &
                          'S'  ,'S'  ,'U'  ,'V'  ,        &
                          'S'  ,'S'  ,'U'  ,'V'  ,        &
                          'S'  ,'S'  ,'S'  ,'S'  ,        &
                          'U'  ,'V'  ,'S'  ,'U'  ,'V',    &
                          'S'  ,'S'  ,'S'  ,'S'  ,        &
                          'S'  ,'S'  ,'S'  ,'S'        /

      character(len=151) :: cplocn_iris_consumes
      data cplocn_iris_consumes /'MCP,ALO,ALI,T40,T4I,SHO,SHI,LHO,LHI,TXO/TYO,TXI/TYI,ZTO,ZQO,ZUO,ZVO,ZTI,ZQI,ZUI/ZVI,GLI,I8I,SDI,TMO,UUO/VVO,I7I,UUI/VVI,QSO,QSI,ILO,ILI,ZMO,ZMI,ZHO,ZHI'/

      real, dimension (:,:,:),   pointer :: ocn_busou

! GMM access section (busin)
      real, dimension (:,:,:),   pointer :: ocn_busin  => null()
      real, dimension (:,:  ),   pointer :: gli_0      => null()
      real, dimension (:,:  ),   pointer :: i8i_0      => null()
      real, dimension (:,:  ),   pointer :: sdi_0      => null()
      real, dimension (:,:  ),   pointer :: tmo_0      => null()

      character(len=GMM_MAXNAMELENGTH) :: gmmk_ocn_busin_s, &
                              gmmk_gli_0_s, gmmk_i8i_0_s, &
                              gmmk_sdi_0_s, gmmk_tmo_0_s

      integer :: icvin_MCP, icvin_ALO, icvin_ALI,            &
                 icvin_T4O, icvin_T4I, icvin_SHO, icvin_SHI, &
                 icvin_LHO, icvin_LHI, icvin_TXO, icvin_TYO, &
                 icvin_TXI, icvin_TYI, icvin_ZTO, icvin_ZQO, &
                 icvin_ZUO, icvin_ZVO, icvin_ZTI, icvin_ZQI, &
                 icvin_ZUI, icvin_ZVI, icvin_GLI, icvin_I8I, &
                 icvin_SDI, icvin_TMO, icvin_UUO, icvin_VVO, &
                 icvin_I7I, icvin_UUI, icvin_VVI,            &
                 icvin_QSO, icvin_QSI, icvin_ILO, icvin_ILI, &
                 icvin_ZMO, icvin_ZMI, icvin_ZHO, icvin_ZHI

! GMM access (busou)
      character(len=1), dimension(cplocn_n_fldou) :: cplocn_cvou_G
      character(len=GMM_MAXNAMELENGTH), &
                        dimension(cplocn_n_fldou) :: cplocn_cvou_N
      integer, dimension(2,cplocn_n_fldou)        :: cplocn_cvou_K

end module cplocn_mod
