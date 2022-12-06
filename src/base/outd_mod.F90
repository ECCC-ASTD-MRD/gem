module outd
   use dimout, only: MAXELEM,MAXSET
   implicit none
   public
   save
!
!**comdeck outd.cdk
!
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH OUTPUT FOR DYNAMIC VARIABLES              |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Outd_sets          | The total number of sets of output variables    |
!                    | requested from the dynamics                     |
! Outd_var_S         | A table containing "Outd_sets" sets of output   |
!                    | variables requested from the dynamics           |
! Outd_varnm_S       | A table containing "Outd_sets" sets of output   |
!                    | variables (longname) requested from the dynamics|
! Outd_nbit          | A table containing "Outd_sets" sets of compact  |
!                    | bit corresponding to Outd_var_S                 |
! Outd_filtpass      | A table containing "Outd_sets" sets of filter   |
!                    | pass corresponding to Outd_var_S                |
! Outd_filtcoef      | A table containing "Outd_sets" sets of filter   |
!                    | coef corresponding to Outd_var_S                |
! Outd_var_max       | Outd_var_max(j) is the total number of requested|
!                    | output variables in each set Outd_var(*,j)      |
! Outd_lev           | Outd_lev(j) contains the index that indicates   |
!                    | the Level set to use for Outd_var(*,j)          |
! Outd_grid          | Outd_grid(j) contains the index that indicates  |
!                    | the Grid set to use for Outd_var(*,j)           |
! Outd_step          | Outd_step(i) contains the index that indicates  |
!                    | the Step set to use for Outd_var(*,j)           |
!----------------------------------------------------------------------
!
!
   character(len=4) Outd_var_S(MAXELEM,MAXSET)
   character(len=16) Outd_varnm_S(MAXELEM,MAXSET)
   integer Outd_nbit(MAXELEM,MAXSET)
   integer Outd_filtpass(MAXELEM,MAXSET)
   real    Outd_filtcoef(MAXELEM,MAXSET)
   real    Outd_convmult(MAXELEM,MAXSET)
   real    Outd_convadd (MAXELEM,MAXSET)
   integer, dimension(MAXSET), target ::  Outd_lev,Outd_grid,Outd_step
   integer Outd_var_max(MAXSET)
   integer Outd_sets
end module outd
