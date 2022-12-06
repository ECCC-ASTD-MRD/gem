module outc
   use dimout, only: MAXELEM,MAXSET
   implicit none
   public
   save
!
!**comdeck outc.cdk
!
!______________________________________________________________________
!                                                                      |
!  VARIABLES ASSOCIATED WITH OUTPUT FOR CHEMISTRY VARIABLES            |
!______________________________________________________________________|
!                    |                                                 |
! NAME               | DESCRIPTION                                     |
!--------------------|-------------------------------------------------|
! Outc_sets          | The total number of sets of output variables    |
!                    | requested from the chemistry                    |
! Outc_var_S         | A table containing "Outc_sets" sets of output   |
!                    | variables requested from the chemistry          |
! Outc_varnm_S       | A table containing "Outc_sets" sets of output   |
!                    | variables (longname) requested from the chem    |
! Outc_nbit          | A table containing "Outc_sets" sets of compact  |
!                    | bit corresponding to Outc_var_S                 |
! Outc_filtpass      | A table containing "Outc_sets" sets of filter   |
!                    | pass corresponding to Outc_var_S                |
! Outc_filtcoef      | A table containing "Outc_sets" sets of filter   |
!                    | coef corresponding to Outc_var_S                |
! Outc_var_max       | Outc_var_max(j) is the total number of requested|
!                    | output variables in each set Outc_var(*,j)      |
! Outc_lev           | Outc_lev(j) contains the index that indicates   |
!                    | the Level set to use for Outc_var(*,j)          |
! Outc_grid          | Outc_grid(j) contains the index that indicates  |
!                    | the Grid set to use for Outc_var(*,j)           |
! Outc_step          | Outc_step(i) contains the index that indicates  |
!                    | the Step set to use for Outc_var(*,j)           |
!----------------------------------------------------------------------
!

   character(len=4) Outc_var_S(MAXELEM,MAXSET)
   character(len=16) Outc_varnm_S(MAXELEM,MAXSET)
   logical Outc_dostep_L(MAXSET)
   integer Outc_nbit(MAXELEM,MAXSET)
   integer Outc_filtpass(MAXELEM,MAXSET)
   real    Outc_filtcoef(MAXELEM,MAXSET)
   real    Outc_convmult(MAXELEM,MAXSET)
   real    Outc_convadd (MAXELEM,MAXSET)
   integer, dimension(MAXSET), target :: Outc_lev,Outc_grid,Outc_step
   integer Outc_var_max(MAXSET)
   integer Outc_sets

end module outc
