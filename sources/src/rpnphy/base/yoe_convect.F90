MODULE YOE_CONVECT 
IMPLICIT NONE
 
SAVE
! maybe specifiy values later in su0phy.F90
character(len=8):: CCUMFSCHEM 

integer :: NBDIA            ! starting level (from bottom=1) for convection computat.
integer :: NENSM            ! number of additional ensemble members for deep
integer :: NCH1  ! number of chemical tracers
integer :: NICE             ! use ice computations or only liquid (1=with ice 0 without)

LOGICAL   :: LDEEP ! LMFPEN   ! switch for deep convection
LOGICAL   :: LSHAL ! LSHCV    ! switch for shallow convection
LOGICAL   :: LDOWN ! LMFDD    ! take or not convective downdrafts into account 
LOGICAL   :: LREFRESH         ! refresh convective tendencies at every time step
LOGICAL   :: LSETTADJ         ! logical to set convective adjustment time by user
LOGICAL   :: LUVTRANS         ! flag to compute convective transport for hor. wind
LOGICAL   :: LCHTRANS         ! flag to compute convective transport for chemical tracers
real    :: RTADJD             ! user defined deep  adjustment time (s) (if LSETTADJ)
real    :: RTADJS             ! user defined shallow adjustment time (s) (if LSETTADJ)
LOGICAL   :: LDIAGCONV2D      ! logical for convection 2D diagnostics
LOGICAL   :: LDIAGCONV3D      ! logical for convection 3D diagnostics

END MODULE YOE_CONVECT
