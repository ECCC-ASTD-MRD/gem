module testphy
  ! Implementation of RPN physics test stub
  use, intrinsic :: iso_fortran_env, only: INT64, REAL64
  use vgrid_descriptors, only: vgrid_descriptor
  use clib_itf_mod, only: clib_toupper
  use wb_itf_mod
  use mu_jdate_mod
  use rmn_gmm
  implicit none
  private

  ! External definitions
#include <rmnlib_basics.hf>

  ! External parameters
  character(len=*), dimension(2), parameter, public :: PT_TEST_TYPES = (/ &
       'full', &
       'itf '  &
       /)

  ! External subprograms
  public :: pt_run

  ! Internal parameters
  integer, parameter :: NI=5,NJ=5,NK=10                                         !domain extent (i,j,k)
  real, parameter :: GRID_DX=1.                                                 !grid spacing (deg)
  real, parameter :: STEP_DT=900.                                               !timestep (sec)
  real, parameter :: STATE_TT=270.                                              !prescribed state temperature (K)
  real, dimension(2), parameter :: GRID_CENTRE=(/45.,270./)                     !grid centre coordiantes (deg lat/lon)
  character(len=GMM_MAXNAMELENGTH), parameter, private :: &
       gmmk_pw_tt_plus_s='PW_TT:P', &                                           !temperature (K, time plus)
       gmmk_pw_uu_plus_s='PW_UU:P', &                                           !west wind component (m/s, time plus)
       gmmk_pw_vv_plus_s='PW_VV:P', &                                           !south wind component (m/s, time plus)
       gmmk_pw_wz_plus_s='PW_WZ:P', &                                           !vertical wind (m/s, time plus)
       gmmk_pw_gz_moins_s='PW_GZ:M', &                                          !geopotential (m^2/s^2 ASL, time minus)
       gmmk_pw_pm_moins_s='PW_PM:M', &                                          !momentum level pressure (Pa, time minus)
       gmmk_pw_pt_moins_s='PW_PT:M', &                                          !thermodynamic level pressure (Pa, time minus)
       gmmk_pw_tt_moins_s='PW_TT:M', &                                          !temperature (K, time minus)
       gmmk_pw_uu_moins_s='PW_UU:M', &                                          !west wind component (m/s, time minus)
       gmmk_pw_vv_moins_s='PW_VV:M'                                             !south wind component (m/s, time minus)
  character(len=*), parameter :: GRID_NAME='grid/phy_local'                     !name of physics grid

  ! Internal module variables
  character(len=2048) :: input_path
  type(vgrid_descriptor), save :: vcoord

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function pt_run(F_type) result(F_istat)
    ! Main driver for the test stub

    ! Input arguments
    character(len=*), intent(in) :: F_type              !test type

    ! Output variables
    integer :: F_istat                                  !return status

    ! Internal variables
    integer :: istat

    ! Set error return status
    F_istat = RMN_ERR

    ! Set up basic environment
    istat = set_env()
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by set_env:'//trim(F_type))

    ! Set up grid and coordinate
    istat = set_grid()
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by set_grid')

    ! Initialize the physics
    istat = itf_phy_init()
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by itf_phy_init')

    ! Set up GMM space
    istat = set_vt()
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by set_vt')

    ! Fill GMM state variables
    istat = set_state()
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by set_state')

    ! Take a kount-0 (initialization) physics step
    istat = itf_phy_step(0)
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by kount=0 itf_phy_step')

    ! Take a full physics step
    istat = itf_phy_step(1)
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by itf_phy_step')

    ! Generate output
    istat = itf_phy_output()
    call handle_error_l(istat==RMN_OK,'pt_run','Error returned by itf_phy_output')

    ! Successful completion
    F_istat = RMN_OK
    return
  end function pt_run

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function set_env() result(F_istat)
    ! Set up GMM memory space

    ! Output variables
    integer :: F_istat                                  !return status

    ! External declarations

    ! Internal variables
    integer :: istat

    ! Set error return status
    F_istat = RMN_ERR

    ! Retrieve required environment
    call get_environment_variable('TASK_INPUT',value=input_path,status=istat)
    call handle_error_l(istat==0,'set_env','Environment variable PHY_DFILES must be defined')

    ! Establish writing master for stdout of physics
    istat = WB_OK
    istat = wb_put('model/outout/pe_master',0)
    call handle_error_l(WB_IS_OK(istat),'set_env','Cannot specify model master for stdout')

    ! Successful completion
    F_istat = RMN_OK
    return
  end function set_env

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function set_vt() result(F_istat)
    use phy_itf, only: phymeta,phy_getmeta
    ! Set up GMM memory space

    ! Output variables
    integer :: F_istat                                  !return status

    ! External declarations

    ! Internal parameters
    integer, parameter :: MAX_TRACERS=1000

    ! Internal variables
    integer :: i,ntr,istat
!!$    real, dimension(:,:), pointer :: ptr2d
    real, dimension(:,:,:), pointer :: ptr3d
    character(len=GMM_MAXNAMELENGTH) :: vname,prefix,basename,time,ext
    character(len=GMM_MAXNAMELENGTH), dimension(MAX_TRACERS) :: tr_list
    type(gmm_metadata) :: meta3d_nk,meta3d_nk1
    type(phymeta), dimension(:), pointer :: pmeta

    ! Set error return status
    F_istat = RMN_ERR

    ! Extablish metadata structures
    call gmm_build_meta3D(meta3d_nk,1,NI,0,0,NI,1,NJ,0,0,NJ,1,NK,0,0,NK,0,GMM_NULL_FLAGS)
    call gmm_build_meta3D(meta3d_nk1,1,NI,0,0,NI,1,NJ,0,0,NJ,1,NK+1,0,0,NK+1,0,GMM_NULL_FLAGS)

    ! Create state GMM space
    nullify(ptr3d)
    istat = GMM_OK
    istat = min(istat,gmm_create(gmmk_pw_tt_plus_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_uu_plus_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_vv_plus_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_wz_plus_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_gz_moins_s,ptr3d,meta3d_nk1),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_tt_moins_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_pm_moins_s,ptr3d,meta3d_nk1),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_pt_moins_s,ptr3d,meta3d_nk1),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_uu_moins_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    istat = min(istat,gmm_create(gmmk_pw_vv_moins_s,ptr3d,meta3d_nk),istat); nullify(ptr3d)
    call handle_error_l(GMM_IS_OK(istat),'set_vt','Cannot create state GMM space')

    ! Create tracer GMM space and initialize to zero
    nullify(pmeta)
    istat = phy_getmeta(pmeta,' ',F_npath='V',F_bpath='D',F_quiet=.true.)
    call handle_error_l(RMN_IS_OK(istat),'set_vt','Error retrieving dynamics bus list')
    tr_list = ''
    ntr = 1
    tr_list(ntr) = 'HU'
    do i=1,size(pmeta)
       vname = pmeta(i)%vname
       istat = clib_toupper(vname)
       if (vname(1:3) /= 'TR/') cycle
       call gmmx_name_parts(vname,prefix,basename,time,ext)
       if (any(tr_list(1:ntr) == basename)) cycle
       ntr = ntr + 1
       tr_list(ntr) = basename
    enddo
    istat = GMM_OK; nullify(ptr3d)
    do i=1,ntr
       istat = min(istat,gmm_create('TR/'//trim(tr_list(i))//':M',ptr3d,meta3d_nk,GMM_FLAG_IZER)); nullify(ptr3d)
       istat = min(istat,gmm_create('TR/'//trim(tr_list(i))//':P',ptr3d,meta3d_nk,GMM_FLAG_IZER)); nullify(ptr3d)
    enddo
    call handle_error_l(GMM_IS_OK(istat),'set_vt','Cannot create tracer GMM space')

    ! Successful completion
    F_istat = RMN_OK
    return
  end function set_vt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function set_grid() result(F_istat)
    ! Set up horizontal and vertical grid
    use hgrid_wb, only: hgrid_wb_put
    use vgrid_descriptors, only: VGD_OK,vgd_new,vgd_get
    use vgrid_wb, only: vgrid_wb_put

    ! Output variables
    integer :: F_istat                                  !return status
 
    ! Local variables
    integer :: i,j,k,gid,ig1,ig2,ig3,ig4,istat
    integer, dimension(:), pointer :: ip1m,ip1t
    real, dimension(NI) :: ax
    real, dimension(NJ) :: ay
    real, dimension(NK) :: hyb
    real(kind=8) :: ptop
    character(len=1) :: gtype
 
    ! External subprograms
    integer, external :: ezgdef_fmem

    ! Set error return status
    F_istat = RMN_ERR

    ! Fill positional vectors
    do i=1,size(ax)
       ax(i) = GRID_CENTRE(2) + ((i-1) - NI/2)
    enddo
    do j=1,size(ay)
       ay(j) = GRID_CENTRE(1) + ((j-1) - NJ/2)
    enddo

    ! Create horizontal grid
    gtype = 'L'
    call cxgaig(gtype,ig1,ig2,ig3,ig4,ay(1),ax(1),GRID_DX,GRID_DX)
    gid = ezgdef_fmem(NI,NJ,'Z',gtype,ig1,ig2,ig3,ig4,ax,ay)

    ! Store horizontal grid in whiteboard
    istat = hgrid_wb_put(GRID_NAME,gid,F_i0=1,F_j0=1,F_lni=NI,F_lnj=NJ)
    call handle_error_l(RMN_IS_OK(istat),'set_grid','Error saving horizontal grid')

    ! Construct vertical grid
    hyb(1) = 0.0001; hyb(NK) = 0.995
    do k=2,NK-1
       hyb(k) = hyb(1) + (k-1)*(hyb(NK)-hyb(1))/real(NK-1)
    enddo
    istat = vgd_new(vcoord,kind=5,version=4,hyb=hyb,rcoef1=1.,rcoef2=1., &
         ptop_8=7.5d0,pref_8=1d5,ptop_out_8=ptop)
    call handle_error_l(istat==VGD_OK,'set_grid','Error constructing vertical coordinate')
 
    ! Save vertical grid for physics input
    nullify(ip1m,ip1t)
    istat = vgd_get(vcoord,'VIPM - IP1 MOMENTUM',ip1m)
    call handle_error_l(istat==VGD_OK,'set_grid','Retrieving IP1(M)')
    istat = vgrid_wb_put('ref-m',vcoord,ip1m)
    call handle_error(istat,'set_grid','Creating vgrid_wb entry for IP1(M)')
    istat = vgd_get(vcoord,'VIPT - IP1 THERMO',ip1t)
    call handle_error_l(istat==VGD_OK,'set_grid','Retrieving IP1(T)')
    istat = vgrid_wb_put('ref-t',vcoord,ip1t)
    call handle_error(istat,'set_grid','Creating vgrid_wb entry for IP1(T)')
    deallocate(ip1m,ip1t); nullify(ip1m,ip1t)

    ! Successful completion
    F_istat = RMN_OK
    return
  end function set_grid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function set_state() result(F_istat)
    ! Create an atmospheric state
    use vgrid_descriptors, only: VGD_OK,vgd_get,vgd_levels

    ! Output variables
    integer :: F_istat                                  !return status
    
    ! Internal parameters
    real, parameter :: GRAV=9.81,R=287.

    ! Internal variables
    integer :: k,istat
    integer, dimension(:), pointer :: ip1m,ip1t
    real :: p0
    real, dimension(NK) :: tt
    real, dimension(NK+1) :: gz
    real, dimension(:), pointer :: presm,prest
    real, dimension(:,:,:), pointer :: ptr3d

    ! Set error return status
    F_istat = RMN_ERR

    ! Set up pressure profile for momentum and thermodynamic levels (hPa)
    p0 = 1e5
    nullify(ip1m,ip1t,presm,prest)
    istat = vgd_get(vcoord,'VIPM - IP1 MOMENTUM',ip1m)
    call handle_error_l(istat==VGD_OK,'set_state','Retrieving VIPM')
    istat = vgd_levels(vcoord,ip1m,presm,sfc_field=p0)
    call handle_error_l(istat==VGD_OK,'set_state','Computing momentum-level pressures')
    istat = vgd_get(vcoord,'VIPT - IP1 THERMO',ip1t)
    call handle_error_l(istat==VGD_OK,'prof_update_pres','Retrieving VIPT')
    istat = vgd_levels(vcoord,ip1t,prest,sfc_field=p0)
    call handle_error_l(istat==VGD_OK,'prof_update_pres','Computing thermo-level pressures')

    ! Set up temperature profile (K)
    tt = STATE_TT

    ! Set up hydrostatic geopotential profile on momentum levels (m^2/s^2 ASL)
    gz(NK+1) = 100.
    do k=NK,1,-1
       gz(k) = gz(k+1) + R*tt(k)/GRAV * log(presm(k+1)/presm(k))
    enddo
    gz = gz*GRAV

    ! Create a highly idealized initial state at time-minus
    istat = GMM_OK
    nullify(ptr3d)
    istat = min(istat,gmm_get(gmmk_pw_uu_moins_s,ptr3d),istat)
    ptr3d = 10.
    istat = min(istat,gmm_get(gmmk_pw_vv_moins_s,ptr3d),istat)
    ptr3d = 0.
    istat = min(istat,gmm_get(gmmk_pw_tt_moins_s,ptr3d),istat)
    do k=1,size(tt)
       ptr3d(:,:,k) = tt(k)
    enddo
    istat = min(istat,gmm_get(gmmk_pw_pm_moins_s,ptr3d),istat)
    do k=1,size(presm)
       ptr3d(:,:,k) = presm(k)
    enddo
    istat = min(istat,gmm_get(gmmk_pw_pt_moins_s,ptr3d),istat)
    do k=1,size(prest)
       ptr3d(:,:,k) = prest(k)
    enddo
    istat = min(istat,gmm_get(gmmk_pw_gz_moins_s,ptr3d),istat)
    do k=1,size(gz)
       ptr3d(:,:,k) = gz(k)
    enddo

    ! Create a highly idealized initial state at time-plus
    istat = GMM_OK
    nullify(ptr3d)
    istat = min(istat,gmm_get(gmmk_pw_uu_plus_s,ptr3d),istat)
    ptr3d = 10.
    istat = min(istat,gmm_get(gmmk_pw_vv_plus_s,ptr3d),istat)
    ptr3d = 0.
    istat = min(istat,gmm_get(gmmk_pw_wz_plus_s,ptr3d),istat)
    ptr3d = 0.01
    istat = min(istat,gmm_get(gmmk_pw_tt_plus_s,ptr3d),istat)
    do k=1,size(tt)
       ptr3d(:,:,k) = tt(k)
    enddo

    ! Successful completion
    F_istat = RMN_OK
    return
  end function set_state

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function itf_phy_init() result(F_istat)
    use phy_itf, only: PHY_COMPATIBILITY_LVL,phy_nml,phy_init
    use vgrid_descriptors, only: VGD_OK,vgd_get,vgd_levels
    ! Initialize the physics

    ! Output variables
    integer :: F_istat                                  !return status

    ! External declarations
 
    ! Internal parameters
    !#TODO: this modules needs to be updated, we are now at API 17
    integer, parameter :: COMPATIBILITY_LVL=20 !#TODO: was 6... may need to review

    ! Internal variables
    integer(INT64) :: jdateo
    integer :: istat,dateo
    integer, dimension(:), pointer :: ip1m
    real, dimension(:), pointer :: std_p_prof
    character(len=MU_JDATE_PDF_LEN) :: dateo_S

    ! Set error return status
    F_istat = RMN_ERR

    ! Check for physics interface compatibility
    if (PHY_COMPATIBILITY_LVL /= COMPATIBILITY_LVL) then
       write(RMN_STDERR,*) '(itf_phy_init) Incompatible physics API'
       write(RMN_STDERR,*) '  Current: ',COMPATIBILITY_LVL
       write(RMN_STDERR,*) '  Required: ',PHY_COMPATIBILITY_LVL
       return
    endif

    ! Start physics initialization with mandatory parameters
    istat = WB_OK
    istat= min(wb_put('itf_phy/TLIFT',0),istat)
    call handle_error_l(WB_IS_OK(istat),'itf_phy_init','Cannot fill required whiteboard entries')

    ! Initialize physics configuration
    istat = phy_nml(trim(input_path)//'/configs/phy_settings.nml')
    call handle_error_l(RMN_IS_OK(istat),'itf_phy_init','Cannot process physics namelist')

    ! Create standard pressure profile for the physics on momentum
    istat = vgd_get(vcoord,'VIPM - IP1 MOMENTUM',ip1m)
    call handle_error_l(istat==VGD_OK,'itf_phy_init','Retrieving VIPM')
    istat = vgd_levels(vcoord,ip1m,std_p_prof,sfc_field=1e5)
    call handle_error_l(istat==VGD_OK,'itf_phy_init','Cannot create a standard pressure profile')

    ! Obtain required values for final initialization
    dateo_S = '20090427.000000'
    jdateo = jdate_from_print(dateo_S)
    dateo  = jdate_to_cmc(jdateo)

    ! Complete physics initialization
    istat = phy_init(trim(input_path)//'/MODEL_INPUT/',dateo,STEP_DT,GRID_NAME,GRID_NAME,NK+1,std_p_prof)
    call handle_error_l(RMN_IS_OK(istat),'itf_phy_init','Cannot complete physics initialization')

    ! Successful completion
    F_istat = RMN_OK
    return
  end function itf_phy_init

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function itf_phy_step(F_kount) result(F_istat)
    ! Take a physics step
    use phy_itf, only: phy_input1, phy_step

    ! Input arguments
    integer, intent(in) :: F_kount                      !step number

    ! Output variables
    integer :: F_istat                                  !return status
    
    ! Internal variables
    integer :: istat

    ! Set error return status
    F_istat = RMN_ERR

    ! Create grid-specific fields
    if (F_kount == 0) then
       istat = itf_phy_geom()
       call handle_error_l(RMN_IS_OK(istat),'itf_phy_init','Cannot create grid-specific entries')
    endif
       
    ! Process physics inputs
    istat = phy_input1(prefold_opr,F_kount,trim(input_path)//'/configs/physics_input_table', &
         trim(input_path)//'/MODEL_INPUT/','geophy/Phy_geophy.fst')
    call handle_error_l(RMN_IS_OK(istat),'itf_phy_init','Cannot complete physics input')
    
    ! Take a physics step
    write(RMN_STDOUT,1002) F_kount
    istat = phy_step(F_kount,F_kount)
    call handle_error_l(RMN_IS_OK(istat),'itf_phy_init','Cannot complete physics step')

    ! Info messages
1002 format(/'PERFORM A PHYSICS STEP: stepno= ',i6,' (S/R itf_phy_step)'/58('='))

    ! Successful completion
    F_istat = RMN_OK
    return
  end function itf_phy_step

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function itf_phy_geom() result(F_istat)
    ! Specify grid-specific fields for the physics
    use hgrid_wb, only: hgrid_wb_get

    ! Output variables
    integer :: F_istat                                  !return status

    ! Internal parameters
    real, parameter :: RAYT=6371e3
    
    ! Internal variables
    integer :: gid,istat
    real :: deg2rad
    real, dimension(NI,NJ) :: lat,lon
    real, dimension(:,:), pointer :: ptr2d
    type(gmm_metadata) :: meta2d

    ! External subprograms
    integer, external :: gdll

    ! Set error return status
    F_istat = RMN_ERR

    ! Basic setup
    deg2rad=acos(-1d0)/180d0

    ! Establish GMM metadata structure for 2D fields
    call gmm_build_meta2D(meta2d,1,NI,0,0,NI,1,NJ,0,0,NJ,0,GMM_NULL_FLAGS)

    ! Retrieve grid information
    istat = hgrid_wb_get(GRID_NAME,gid)
    istat = min(gdll(gid,lat,lon),istat)
    call handle_error_l(RMN_IS_OK(istat),'itf_phy_geom','Unable to retrieve grid info for '//trim(GRID_NAME))

    ! Create GMM entry for latitude
    nullify(ptr2d)
    istat = gmm_create('DLAT',ptr2d,meta2d)
    call handle_error_l(associated(ptr2d),'itf_phy_geom','Cannot create latitude entry (DLAT)')
    ptr2d = deg2rad * lat

    ! Create GMM entry for longitude
    nullify(ptr2d)
    istat = gmm_create('DLON',ptr2d,meta2d)
    call handle_error_l(associated(ptr2d),'itf_phy_geom','Cannot create longitude entry (DLON)')
    ptr2d = deg2rad * lon

    ! Create GMM entry for grid area
    nullify(ptr2d)
    istat = gmm_create('DXDY',ptr2d,meta2d)
    call handle_error_l(associated(ptr2d),'itf_phy_geom','Cannot create cell area entry (DXDY)')
    ptr2d = (GRID_DX*deg2rad*RAYT) * (GRID_DX*deg2rad*RAYT*cos(deg2rad*lat))
 
    ! Create GMM entry for tendency mask (0. to 1.)
    nullify(ptr2d)
    istat = gmm_create('TDMASK',ptr2d,meta2d)
    call handle_error_l(associated(ptr2d),'itf_phy_geom','Cannot create tendency mask entry (TDMASK)')
    ptr2d = STEP_DT

    ! Successful completion
    F_istat = RMN_OK
    return
  end function itf_phy_geom

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function prefold_opr(F_data,F_name_S,F_horiz_interp_S, &
       F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn) result(F_istat)
    ! Perform necessary steps before physics folding.  For this
    ! test driver, this simply means returning the data in-place.
 
    ! Input arguments
    character(len=*),intent(in) :: F_name_S,F_horiz_interp_S
    integer,intent(in) :: F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn

    ! Input/output arguments
    real,intent(inout) :: F_data(F_minx:F_maxx,F_miny:F_maxy,F_k0:F_kn)

    ! Output arguments
    integer :: F_istat                                  !return status

    ! No folding required for test stub
    F_istat = RMN_OK

    if (F_kn<0) print *,'testphy, prefold_opr called with:',F_name_S,F_horiz_interp_S, &
       F_minx,F_maxx,F_miny,F_maxy,F_k0,F_kn,F_data
 
  end function prefold_opr

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  function itf_phy_output() result(F_istat)
    ! Take a physics step
    use phy_itf, only: phymeta,phy_get

    ! Output arguments
    integer :: F_istat                                  !return status

    ! External declarations

    ! Internal variables
    integer :: istat
    real, dimension(:,:,:), pointer :: ptr3d
    logical :: test_phy
    type(phymeta) :: pmeta

    ! Set error return status
    F_istat = RMN_OK

    ! Retrieve information about simplified physics tests
    call handle_error_l(RMN_IS_OK(wb_get('phy/test_phy',test_phy)),'itf_phy_output','Cannot retrieve test information')

    ! Validate test physics output
    PHY_TEST_MODE: if (test_phy) then

       ! Retrieve computed temperature and compare
       nullify(ptr3d)
       istat = phy_get(ptr3d,gmmk_pw_tt_plus_s)
       call handle_error_l(RMN_IS_OK(istat),'itf_phy_output','Cannot retrieve temperature')
       if (any(abs(ptr3d(:,:,1:NK)-STATE_TT-1.) > epsilon(ptr3d))) then
          write(RMN_STDOUT,'(a,1x,f6.3)') '*FAIL: Test double returned an unexpected temperature difference of', &
               maxval(abs(ptr3d(:,:,1:NK)-STATE_TT-1.))
          F_istat = RMN_ERR
       endif
       deallocate(ptr3d); nullify(ptr3d)

       ! Retrieve prescribed PBL height and compare
       nullify(ptr3d)
       istat = phy_get(ptr3d,'H',F_npath='O',F_bpath='P',F_meta=pmeta)
       call handle_error_l(RMN_IS_OK(istat),'itf_phy_output','Cannot retrieve temperature')
       if (pmeta%nk /= 1) then
          write(RMN_STDOUT,'(a,1x,i3)') '*FAIL: Test double returned an unexpected PBL height dimension of',pmeta%nk
          F_istat = RMN_ERR
       endif
       if (any(abs(ptr3d-100.) > epsilon(ptr3d))) then
          write(RMN_STDOUT,'(a,1x,f6.3)') '*FAIL: Test double returned an unexpected PBL height error of',maxval(abs(ptr3d-100.))
          F_istat = RMN_ERR
       endif
       deallocate(ptr3d); nullify(ptr3d)

       ! Emit message confirming that the test passed
       if (F_istat == RMN_OK) write(RMN_STDOUT,*) '*PASS'

    endif PHY_TEST_MODE

    ! Successful completion
    return
  end function itf_phy_output

end module testphy

