module outp
   use dimout, only: MAXELEM,MAXSET
   use rmn_gmm
   implicit none
   public
   save
!
!**comdeck Outp.cdk
!
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH OUTPUT FOR PHYSICS VARIABLES              |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Outp_sets          | The total number of sets of output variables    |
!                    | requested from the physics                      |
! Outp_var_S         | A table containing "Outp_sets" sets of output   |
!                    | variables requested from the physics            |
! Out_varnm_S        | A table containing "Outp_sets" sets of output   |
!                    | variables (longname) requested from the physics |
! Outp_nbit          | A table containing "Outp_sets" sets of compact  |
!                    | bit corresponding to Outp_var_S                 |
! Outp_filtpass      | A table containing "Outp_sets" sets of filter   |
!                    | pass corresponding to Outp_var_S                |
! Outp_filtcoef      | A table containing "Outp_sets" sets of filter   |
!                    | coef corresponding to Outp_var_S                |
! Outp_var_max       | Outp_var_max(j) is the total number of requested|
!                    | output variables in each set Outp_var(*,j)      |
! Outp_lev           | Outp_lev(j) contains the index that indicates   |
!                    | the Level set to use for Outp_var(*,j)          |
! Outp_grid          | Outp_grid(j) contains the index that indicates  |
!                    | the Grid set to use for Outp_var(*,j)           |
! Outp_step          | Outp_step(j) contains the index that indicates  |
!                    | the Step set to use for Outp_var(*,j)           |
! Outp_numstep       | Outp_numstep(j) contains the number of steps    |
!                    | elapsed since the last avg/acc for Outp_var(*,j)|
! Outp_avg_L         | Outp_avg_L(j)    contains TRUE or FALSE         |
!                    | for Outp_var(*,j)                               |
! Outp_accum_L       | Outp_accum_L(j)  contains TRUE or FALSE         |
!                    | for Outp_var(*,j)                               |
! gmmk_diag_tt_s     | GMM name for diagnostic level temperature       |
! gmmk_diag_hu_s     | GMM name for diagnostic level specific humidity |
! gmmk_diag_uu_s     | GMM name for diagnostic level u-wind            |
! gmmk_diag_vv_s     | GMM name for diagnostic level v-wind            |
!----------------------------------------------------------------------
!
   character(len=4) Outp_var_S(MAXELEM,MAXSET)
   character(len=16) Outp_varnm_S(MAXELEM,MAXSET)
   character(len=GMM_MAXNAMELENGTH) :: gmmk_diag_tt_s, gmmk_diag_hu_s
   character(len=GMM_MAXNAMELENGTH) :: gmmk_diag_uu_s, gmmk_diag_vv_s
   logical Outp_avg_L(MAXSET),Outp_accum_L(MAXSET)
   integer Outp_nbit(MAXELEM,MAXSET)
   integer Outp_filtpass(MAXELEM,MAXSET)
   integer, dimension(:,:), pointer :: Outp_lasstep
   real    Outp_filtcoef(MAXELEM,MAXSET)
   real    Outp_convmult(MAXELEM,MAXSET)
   real    Outp_convadd (MAXELEM,MAXSET)
   integer, dimension(MAXSET), target :: Outp_lev,Outp_grid,Outp_step
   integer Outp_var_max(MAXSET)
   integer Outp_sets, Outp_multxmosaic
   real, dimension(:,:  ), pointer :: udiag,vdiag,tdiag,qdiag
   real, dimension(:,:,:), allocatable :: diag_dgf
      
end module outp
